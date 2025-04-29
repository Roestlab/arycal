use anyhow::{anyhow, Result as AnyHowResult};
use flate2::read::ZlibDecoder;
use ordered_float::OrderedFloat;
use r2d2::Pool;
use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::Error as RusqliteError;
use rusqlite::Result;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::io::Read;
use rayon::prelude::*;

use arycal_common::chromatogram::Chromatogram;
use crate::msnumpress::{decode_linear, decode_slof};
use crate::osw::PrecursorIdData;
use crate::util::extract_basename;

/// Create a type alias for the connection pool used for SQLite connections.
pub type DbPool = Pool<SqliteConnectionManager>;

/// Represents a precursor entry with its associated data.
#[derive(Debug)]
pub struct Precursor {
    pub chromatogram_id: i32,
    pub native_id: String,
    pub peptide_sequence: String,
    pub charge: i32,
}


/// Represents a group of related chromatograms.
#[derive(Debug, Clone)]
pub struct TransitionGroup {
    /// A unique identifier for the transition group.
    pub group_id: String,
    /// A map of chromatograms, keyed by their native_id.
    pub chromatograms: HashMap<String, Chromatogram>,
    /// Additional metadata for the transition group (if needed). Current keys include "file" and "basename".
    pub metadata: HashMap<String, String>,
}

impl TransitionGroup {
    /// Creates a new TransitionGroup with a given group_id.
    pub fn new(group_id: String) -> Self {
        TransitionGroup {
            group_id,
            chromatograms: HashMap::new(),
            metadata: HashMap::new(),
        }
    }

    /// Adds a chromatogram to the group.
    pub fn add_chromatogram(&mut self, chromatogram: Chromatogram) {
        self.chromatograms
            .insert(chromatogram.native_id.clone(), chromatogram);
    }

    /// Retrieves a chromatogram by its native_id.
    pub fn get_chromatogram(&self, native_id: &str) -> Option<&Chromatogram> {
        self.chromatograms.get(native_id)
    }

    /// Adds metadata to the group.
    pub fn add_metadata(&mut self, key: String, value: String) {
        self.metadata.insert(key, value);
    }

    /// Retrieves metadata by key.
    pub fn get_metadata(&self, key: &str) -> Option<&String> {
        self.metadata.get(key)
    }

    /// Returns the number of chromatograms in the group.
    pub fn chromatogram_count(&self) -> usize {
        self.chromatograms.len()
    }

    /// Calculates the Total Ion Current (TIC) chromatogram by summing intensities across all chromatograms.
    ///
    /// # Returns
    /// A new `Chromatogram` instance representing the TIC.
    pub fn calculate_tic(&self) -> Chromatogram {
        let mut rt_intensity_map: HashMap<OrderedFloat<f64>, f64> = HashMap::new();

        // Iterate through all chromatograms in the group
        for chromatogram in self.chromatograms.values() {
            for (&rt, &intensity) in chromatogram
                .retention_times
                .iter()
                .zip(chromatogram.intensities.iter())
            {
                *rt_intensity_map.entry(OrderedFloat(rt)).or_insert(0.0) += intensity;
            }
        }

        // Convert the HashMap to sorted vectors of retention times and intensities
        let mut sorted_data: Vec<(f64, f64)> = rt_intensity_map
            .into_iter()
            .map(|(rt, intensity)| (rt.into_inner(), intensity))
            .collect();
        sorted_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let retention_times: Vec<f64> = sorted_data.iter().map(|&(rt, _)| rt).collect();
        let intensities: Vec<f64> = sorted_data
            .iter()
            .map(|&(_, intensity)| intensity)
            .collect();

        // Create a new Chromatogram for the TIC
        let mut metadata = HashMap::new();
        metadata.insert(
            "calculation_method".to_string(),
            "Total Ion Current".to_string(),
        );
        metadata.insert(
            "source_chromatograms".to_string(),
            format!("{}", self.chromatograms.len()),
        );

        // Add file name from TransitionGroup metadata
        if let Some(file) = self.get_metadata("file") {
            metadata.insert("file".to_string(), file.clone());
        }
        if let Some(basename) = self.get_metadata("basename") {
            metadata.insert("basename".to_string(), basename.clone());
        }

        Chromatogram {
            id: 0, // TODO: assign a special ID for the TIC?
            native_id: format!("{}_TIC", self.group_id),
            retention_times,
            intensities,
            metadata,
        }
    }
}

/// Define a custom error type for SQ Mass access
#[derive(Debug)]
pub enum SqMassSqliteError {
    DatabaseError(String),
    GeneralError(String),
    RusqliteError(RusqliteError),
}

impl fmt::Display for SqMassSqliteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SqMassSqliteError::DatabaseError(msg) => write!(f, "[SqMassSqliteError] Database Error: {}", msg),
            SqMassSqliteError::GeneralError(msg) => write!(f, "[SqMassSqliteError] Error: {}", msg),
            SqMassSqliteError::RusqliteError(err) => write!(f, "[SqMassSqliteError] Rusqlite Error: {}", err),
        }
    }
}

/// Implement From for SqMassSqliteError to convert rusqlite errors
impl From<RusqliteError> for SqMassSqliteError {
    fn from(err: RusqliteError) -> SqMassSqliteError {
        SqMassSqliteError::RusqliteError(err)
    }
}

/// Implement std::error::Error for SqMassSqliteError
impl Error for SqMassSqliteError {}

/// Define the SqMassAccess struct
pub struct SqMassAccess {
    pub file: String,
    pool: Pool<SqliteConnectionManager>,
}

impl SqMassAccess {
    /// Constructor to create a new SqMassAccess instance with a connection pool
    pub fn new(db_path: &str) -> Result<Self, SqMassSqliteError> {
        log::trace!("Creating SqMassAccess instance with db_path: {}", db_path);
        let manager = SqliteConnectionManager::file(db_path);
        let pool =
            Pool::new(manager).map_err(|e| SqMassSqliteError::DatabaseError(e.to_string()))?;

        // Create the database indices if they don't exist
        log::trace!("Creating indices for the database");
        Self::create_indices(&pool)?;

        Ok(SqMassAccess {
            file: db_path.to_string(),
            pool,
        })
    }

    pub fn create_indices(pool: &Pool<SqliteConnectionManager>) -> Result<(), SqMassSqliteError> {
        // Get a connection from the pool
        let conn = pool
            .get()
            .map_err(|e| SqMassSqliteError::DatabaseError(e.to_string()))?;

        // Create indices for the CHROMATOGRAM and DATA tables
        conn.execute_batch(
            r#"
            CREATE INDEX IF NOT EXISTS idx_chromatogram_native_id ON CHROMATOGRAM(NATIVE_ID);
            CREATE INDEX IF NOT EXISTS idx_data_chromatogram_id ON DATA(CHROMATOGRAM_ID);
            "#,
        )?;

        Ok(())
    }

    /// Reads chromatograms from the database based on specified filter criteria.
    ///
    /// Note: it is a assumed that you would extract chromtaograms for a single transition group at a time. (i.e. a precursor chromatorgram and its associated product ion chromatograms)
    ///
    /// # Parameters
    /// - `filter_type`: A string that specifies the filter type; either "NATIVE_ID" or "CHROMATOGRAM_ID".
    /// - `filter_values`: A vector of string slices that represent the values to filter by. These values should correspond to the filter_type.
    /// - `group_id`: A string representing the group ID for the transition group. Could be peptide charge state, etc.
    ///
    /// # Returns
    /// - A `Result` containing a vector of `Chromatogram` instances on success,
    ///   or an error of type `rusqlite::Error` on failure.
    pub fn read_chromatograms(
        &self,
        filter_type: &str,
        filter_values: Vec<&str>,
        group_id: String,
    ) -> Result<TransitionGroup, SqMassSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| SqMassSqliteError::DatabaseError(e.to_string()))?;

        let placeholders = filter_values
            .iter()
            .map(|id| format!("'{}'", id))
            .collect::<Vec<_>>()
            .join(",");

        let filter_column = match filter_type {
            "NATIVE_ID" => "CHROMATOGRAM.NATIVE_ID",
            "CHROMATOGRAM_ID" => "DATA.CHROMATOGRAM_ID",
            _ => {
                return Err(SqMassSqliteError::GeneralError(
                    "Invalid filter type".to_string(),
                ))
            }
        };

        let query = format!(
            "SELECT DATA.CHROMATOGRAM_ID, CHROMATOGRAM.NATIVE_ID, DATA.COMPRESSION, DATA.DATA_TYPE, DATA.DATA 
             FROM DATA 
             INNER JOIN CHROMATOGRAM ON CHROMATOGRAM.ID = DATA.CHROMATOGRAM_ID 
             WHERE {} IN ({})", 
             filter_column, placeholders
        );

        let mut stmt = conn.prepare(&query).map_err(SqMassSqliteError::from)?;

        let mut transition_group = TransitionGroup::new(group_id);

        // Add file name to metadata
        transition_group.add_metadata("file".to_string(), self.file.to_string());
        let basename = extract_basename(&self.file);
        transition_group.add_metadata("basename".to_string(), basename);

        let data_iter = stmt
            .query_map([], |row| {
                let id: i32 = row.get(0)?;
                let native_id: String = row.get(1)?;
                let compression: i32 = row.get(2)?;
                let data_type: i32 = row.get(3)?;
                let encoded_data: Vec<u8> = row.get(4)?;

                let decoded_data = decompress_data(&encoded_data, compression)?;

                match data_type {
                    2 => Ok((id, native_id, decoded_data, Vec::new())),
                    1 => Ok((id, native_id, Vec::new(), decoded_data)),
                    _ => Err(rusqlite::Error::InvalidColumnType(
                        0,
                        String::from("Expected RT or intensity data"),
                        rusqlite::types::Type::Blob,
                    )),
                }
            })
            .map_err(SqMassSqliteError::from)?;

        for result in data_iter {
            let (id, native_id, retention_times, intensities) = result?;
            let chromatogram = transition_group
                .chromatograms
                .entry(native_id.clone())
                .or_insert(Chromatogram {
                    id,
                    native_id,
                    retention_times: Vec::new(),
                    intensities: Vec::new(),
                    metadata: HashMap::new(),
                });

            if !retention_times.is_empty() {
                chromatogram.retention_times = retention_times;
            } else if !intensities.is_empty() {
                chromatogram.intensities = intensities;
            }
        }

        Ok(transition_group)
    }

    pub fn read_chromatograms_for_precursors(
        &self,
        precursors: &[PrecursorIdData],
        include_precursor: bool,
        num_isotopes: usize,
    ) -> Result<HashMap<i32, TransitionGroup>, SqMassSqliteError> {
        // Collect all native IDs for all precursors in parallel
        let (all_native_ids, precursor_to_native_ids): (Vec<String>, HashMap<i32, Vec<String>>) = precursors
        .par_iter()
        .map(|precursor| {
            let native_ids = precursor.extract_native_ids_for_sqmass(include_precursor, num_isotopes);
            (precursor.precursor_id, native_ids)
        })
        .fold(
            || (Vec::new(), HashMap::new()),
            |(mut all_ids, mut prec_map), (prec_id, native_ids)| {
                all_ids.extend(native_ids.clone());
                prec_map.insert(prec_id, native_ids);
                (all_ids, prec_map)
            }
        )
        .reduce(
            || (Vec::new(), HashMap::new()),
            |(mut all_ids1, mut prec_map1), (all_ids2, prec_map2)| {
                all_ids1.extend(all_ids2);
                prec_map1.extend(prec_map2);
                (all_ids1, prec_map1)
            }
        );
    
        // Get all chromatograms in one query
        let placeholders = all_native_ids
            .iter()
            .map(|id| format!("'{}'", id))
            .collect::<Vec<_>>()
            .join(",");
    
        let query = format!(
            "SELECT DATA.CHROMATOGRAM_ID, CHROMATOGRAM.NATIVE_ID, DATA.COMPRESSION, DATA.DATA_TYPE, DATA.DATA 
             FROM DATA 
             INNER JOIN CHROMATOGRAM ON CHROMATOGRAM.ID = DATA.CHROMATOGRAM_ID 
             WHERE CHROMATOGRAM.NATIVE_ID IN ({})", 
             placeholders
        );

        // log::trace!("Query: {}", query);
    
        let conn = self.pool.get().map_err(|e| {
            SqMassSqliteError::DatabaseError(e.to_string())
        })?;
        let mut stmt = conn.prepare(&query)?;
    
        // Process all results and group by precursor
        let mut precursor_groups = HashMap::new();
        
        let rows = stmt.query_map([], |row| {
            let id: i32 = row.get(0)?;
            let native_id: String = row.get(1)?;
            let compression: i32 = row.get(2)?;
            let data_type: i32 = row.get(3)?;
            let encoded_data: Vec<u8> = row.get(4)?;
            
            let decoded_data = decompress_data(&encoded_data, compression)?;
            
            Ok((id, native_id, decoded_data, data_type))
        })?;
    
        for row in rows {
            let (id, native_id, decoded_data, data_type) = row?;
            
            // Find which precursor this chromatogram belongs to
            for (precursor_id, precursor_native_ids) in &precursor_to_native_ids {
                if precursor_native_ids.contains(&native_id) {
                    let group_id = format!("{}_{}", 
                        precursors.iter().find(|p| p.precursor_id == *precursor_id)
                            .map(|p| p.modified_sequence.clone())
                            .unwrap_or_default(),
                        precursors.iter().find(|p| p.precursor_id == *precursor_id)
                            .map(|p| p.precursor_charge.to_string())
                            .unwrap_or_default()
                    );
                    
                    let group = precursor_groups.entry(*precursor_id)
                        .or_insert_with(|| TransitionGroup {
                            group_id: group_id.clone(),
                            chromatograms: HashMap::new(),
                            metadata: HashMap::new(),
                        });
                    
                    let chrom = group.chromatograms.entry(native_id.clone())
                        .or_insert_with(|| Chromatogram {
                            id,
                            native_id: native_id.clone(),
                            retention_times: Vec::new(),
                            intensities: Vec::new(),
                            metadata: HashMap::new(),
                        });
                    
                    match data_type {
                        2 => chrom.retention_times = decoded_data,
                        1 => chrom.intensities = decoded_data,
                        _ => continue,
                    }
                    
                    break;
                }
            }
        }
    
        // Add metadata to each group
        for group in precursor_groups.values_mut() {
            group.metadata.insert("file".to_string(), self.file.to_string());
            group.metadata.insert("basename".to_string(), extract_basename(&self.file));
        }
    
        Ok(precursor_groups)
    }

    /// Extracts chromatogram IDs along with their corresponding peptides and charge states.
    ///
    /// # Returns
    /// - A `Result` containing a HashMap where keys are tuples of `(peptide_sequence, charge)`
    ///   and values are vectors of `chromatogram_id`s on success,
    ///   or an error of type `rusqlite::Error` on failure.
    pub fn extract_precursors(
        &self,
    ) -> Result<HashMap<(String, i32), Vec<i32>>, SqMassSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| SqMassSqliteError::DatabaseError(e.to_string()))?;

        // SQL query to fetch data from the PRECURSOR table
        let query = "
            SELECT CHROMATOGRAM_ID, NATIVE_ID, PEPTIDE_SEQUENCE, CHARGE
            FROM PRECURSOR
            INNER JOIN CHROMATOGRAM ON CHROMATOGRAM.ID = PRECURSOR.CHROMATOGRAM_ID
            ORDER BY PEPTIDE_SEQUENCE, CHARGE";

        // Prepare the SQL statement
        let mut stmt = conn.prepare(query)?;

        // Create a HashMap to group chromatogram IDs by (peptide_sequence, charge)
        let mut precursor_map: HashMap<(String, i32), Vec<i32>> = HashMap::new();

        // Query the database and process each row
        let data_iter = stmt.query_map([], |row| {
            let chromatogram_id: i32 = row.get(0)?;
            let native_id: String = row.get(1)?;
            let peptide_sequence: String = row.get(2)?;
            let charge: i32 = row.get(3)?;

            Ok(Precursor {
                chromatogram_id,
                native_id,
                peptide_sequence,
                charge,
            })
        })?;

        // Populate the HashMap with grouped data
        for precursor in data_iter {
            let precursor = precursor?;
            let key = (precursor.peptide_sequence.clone(), precursor.charge);

            precursor_map
                .entry(key)
                .or_insert_with(Vec::new)
                .push(precursor.chromatogram_id);
        }

        Ok(precursor_map)
    }

    /// Extracts precursor data for a specific peptide sequence and charge state.
    /// This function returns a vector of `Precursor` instances containing the chromatogram IDs and native IDs.
    ///
    /// # Parameters
    /// - `peptide_sequence`: A string representing the peptide sequence.
    /// - `charge`: An integer representing the charge state.
    ///
    /// # Returns
    /// - A `Result` containing a vector of `Precursor` instances on success,
    pub fn extract_precursor_data(
        &self,
        peptide_sequence: &str,
        charge: i32,
    ) -> Result<Vec<Precursor>, SqMassSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| SqMassSqliteError::DatabaseError(e.to_string()))?;

        // SQL query to fetch data from the PRECURSOR table for the specific peptide and charge
        let query = "
            SELECT CHROMATOGRAM_ID, NATIVE_ID, PEPTIDE_SEQUENCE, CHARGE
            FROM PRECURSOR
            INNER JOIN CHROMATOGRAM ON CHROMATOGRAM.ID = PRECURSOR.CHROMATOGRAM_ID
            WHERE PEPTIDE_SEQUENCE = ?1 AND CHARGE = ?2
            ORDER BY CHROMATOGRAM_ID";

        // Prepare the SQL statement
        let mut stmt = conn.prepare(query)?;

        // Query the database and process each row
        let precursors = stmt
            .query_map([peptide_sequence, &charge.to_string()], |row| {
                Ok(Precursor {
                    chromatogram_id: row.get(0)?,
                    native_id: row.get(1)?,
                    peptide_sequence: row.get(2)?,
                    charge: row.get(3)?,
                })
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(precursors)
    }
}

/// Decompresses encoded data based on the specified compression type.
///
/// # Parameters
/// - `encoded_data`: A byte slice representing compressed data.
/// - `compression`: An integer indicating the compression type. Compression types in sqMass are: 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
///              For this function, currently only zlib compression (type 1) is supported.
///
/// # Returns
/// - A `Result` containing a vector of f64 values on success,
///   or an error of type `rusqlite::Error` on failure.
fn decompress_data(encoded_data: &[u8], compression: i32) -> Result<Vec<f64>, rusqlite::Error> {
    match compression {
        1 => {
            // zlib decompression
            let mut d = ZlibDecoder::new(encoded_data);
            let mut decoded = Vec::new();

            // Use map_err to convert std::io::Error into rusqlite::Error
            d.read_to_end(&mut decoded).map_err(|_e| {
                rusqlite::Error::InvalidColumnType(
                    0,
                    "Decompression Error".to_string(),
                    rusqlite::types::Type::Blob,
                )
            })?;

            // Convert byte data to f64 vector
            let floats: Vec<f64> = decoded
                .chunks_exact(8)
                .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
                .collect();
            Ok(floats)
        }
        5 => {
            // np-linear + zlib
            // zlib decompression
            let mut d = ZlibDecoder::new(encoded_data);
            let mut decoded = Vec::new();

            // Use map_err to convert std::io::Error into rusqlite::Error
            d.read_to_end(&mut decoded).map_err(|_e| {
                rusqlite::Error::InvalidColumnType(
                    0,
                    "Decompression Error".to_string(),
                    rusqlite::types::Type::Blob,
                )
            })?;

            // Decode linear data
            let decoded_res = decode_linear(&decoded).unwrap();

            Ok(decoded_res)
        }
        6 => {
            // np-slof + zlib
            // zlib decompression
            let mut d = ZlibDecoder::new(encoded_data);
            let mut decoded = Vec::new();

            // Use map_err to convert std::io::Error into rusqlite::Error
            d.read_to_end(&mut decoded).map_err(|_e| {
                rusqlite::Error::InvalidColumnType(
                    0,
                    "Decompression Error".to_string(),
                    rusqlite::types::Type::Blob,
                )
            })?;

            // Decode slof data
            let decoded_res = decode_slof(&decoded).unwrap();

            Ok(decoded_res)
        }
        _ => Err(rusqlite::Error::InvalidColumnType(
            0,
            "Unsupported compression type".to_string(),
            rusqlite::types::Type::Blob,
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_extract_precursors() {
        // Use the actual path to your test sqMass database file
        let db_path = "/home/singjc/Documents/github/hackbio2024/data/xics/hroest_K120808_Strep0%PlasmaBiolRepl1_R02_SW.chrom.sqMass";

        // Ensure the database file exists
        assert!(fs::metadata(db_path).is_ok(), "Test database not found");

        // Create an instance of SqMassAccess
        let sq_mass_access = SqMassAccess::new(db_path).expect("Failed to create SqMassAccess");

        // Call the extract_precursors function
        match sq_mass_access.extract_precursors() {
            Ok(precursor_map) => {
                // Assert that the precursor map is not empty
                assert!(
                    !precursor_map.is_empty(),
                    "No precursors were found in the database"
                );
                println!("Length of precursor map: {}", precursor_map.len());

                // Optionally, you can print or verify specific entries in the precursor_map
                for ((peptide_sequence, charge), chromatogram_ids) in &precursor_map {
                    // println!("Peptide: {}, Charge: {}, Chromatogram IDs: {:?}", peptide_sequence, charge, chromatogram_ids);

                    // Assert that each entry has at least one chromatogram ID
                    assert!(
                        !chromatogram_ids.is_empty(),
                        "Chromatogram IDs should not be empty for peptide: {}, charge: {}",
                        peptide_sequence,
                        charge
                    );
                }
            }
            Err(e) => {
                panic!("Error extracting precursors: {:?}", e);
            }
        }
    }

    #[test]
    fn test_extract_precursor_data() {
        // Use the actual path to your test sqMass database file
        let db_path = "/home/singjc/Documents/github/hackbio2024/data/xics/hroest_K120808_Strep0%PlasmaBiolRepl1_R02_SW.chrom.sqMass";

        // Ensure the database file exists
        assert!(fs::metadata(db_path).is_ok(), "Test database not found");

        // Create an instance of SqMassAccess
        let sq_mass_access = SqMassAccess::new(db_path).expect("Failed to create SqMassAccess");

        // Test parameters
        let peptide_sequence = "LTPEAIR";
        let charge = 2;

        // Execute the function
        let result = sq_mass_access.extract_precursor_data(peptide_sequence, charge);

        println!("result: {:?}", result);

        // Assertions
        match result {
            Ok(precursors) => {
                assert!(
                    !precursors.is_empty(),
                    "No precursors found for the given peptide and charge"
                );

                for precursor in precursors {
                    assert_eq!(precursor.peptide_sequence, peptide_sequence);
                    assert_eq!(precursor.charge, charge);
                    assert!(
                        !precursor.native_id.is_empty(),
                        "Native ID should not be empty"
                    );
                    assert!(
                        precursor.chromatogram_id > 0,
                        "Chromatogram ID should be positive"
                    );
                }

                // You might want to add more specific assertions if you know exactly how many
                // precursors or which specific chromatogram IDs you expect
            }
            Err(e) => panic!("Failed to extract precursor data: {}", e),
        }
    }

    #[test]
    fn test_read_chromatograms() {
        // Use the actual path to your test sqMass database file
        let db_path = "/home/singjc/Documents/github/hackbio2024/data/xics/hroest_K120808_Strep0%PlasmaBiolRepl1_R02_SW.chrom.sqMass";

        // Ensure the database file exists
        assert!(fs::metadata(db_path).is_ok(), "Test database not found");

        // Create an instance of SqMassAccess
        let sq_mass_access = SqMassAccess::new(db_path).expect("Failed to create SqMassAccess");

        // --- Test using NATIVE_ID ---
        let native_ids = vec![
            "442_Precursor_i0",
            "23876",
            "23877",
            "23878",
            "23879",
            "23880",
            "23881",
        ];
        let group_id_native = "TestGroupNative".to_string();

        match sq_mass_access.read_chromatograms("NATIVE_ID", native_ids, group_id_native) {
            Ok(transition_group) => {
                // Assert that chromatograms are read correctly
                assert!(
                    !transition_group.chromatograms.is_empty(),
                    "No chromatograms were found for the specified NATIVE_IDs"
                );

                // Print chromatogram details for verification
                for (_native_id, chrom) in transition_group.chromatograms.iter() {
                    // Assert that retention times and intensities are populated
                    assert!(
                        !chrom.retention_times.is_empty(),
                        "Retention times should not be empty"
                    );
                    assert!(
                        !chrom.intensities.is_empty(),
                        "Intensities should not be empty"
                    );
                }

                // Assert that the group ID is correct
                assert_eq!(transition_group.group_id, "TestGroupNative");
            }
            Err(e) => {
                panic!("Error reading chromatograms by NATIVE_ID: {:?}", e);
            }
        }

        // --- Test using CHROMATOGRAM_ID ---
        let chrom_ids = vec![
            "47471", "88888", "90126", "90652", "92822", "93412", "94037",
        ];
        let group_id_chrom = "TestGroupChromID".to_string();

        match sq_mass_access.read_chromatograms("CHROMATOGRAM_ID", chrom_ids, group_id_chrom) {
            Ok(transition_group) => {
                // Assert that chromatograms are read correctly
                assert!(
                    !transition_group.chromatograms.is_empty(),
                    "No chromatograms were found for the specified CHROMATOGRAM_IDs"
                );

                // Print chromatogram details for verification
                for (_native_id, chrom) in transition_group.chromatograms.iter() {
                    // Assert that retention times and intensities are populated
                    assert!(
                        !chrom.retention_times.is_empty(),
                        "Retention times should not be empty"
                    );
                    assert!(
                        !chrom.intensities.is_empty(),
                        "Intensities should not be empty"
                    );
                }

                // Assert that the group ID is correct
                assert_eq!(transition_group.group_id, "TestGroupChromID");
            }
            Err(e) => {
                panic!("Error reading chromatograms by CHROMATOGRAM_ID: {:?}", e);
            }
        }
    }

    #[test]
    fn test_decompress_zlib_data() {
        // Example zlib compressed data
        let encoded_data = vec![
            120, 156, 99, 96, 192, 0, 11, 142, 44, 249, 103, 15, 164, 15, 68, 117, 74, 58, 32, 137,
            55, 228, 5, 8, 58, 96, 42, 199, 9, 28, 78, 36, 61, 183, 71, 226, 31, 184, 171, 242, 8,
            196, 95, 96, 188, 248, 46, 178, 248, 130, 247, 221, 255, 236, 25, 8, 131, 5, 83, 111,
            16, 165, 14, 174, 254, 18, 15, 155, 3, 50, 63, 110, 129, 52, 50, 223, 33, 236, 137, 56,
            50, 159, 225, 136, 146, 8, 50, 63, 33, 248, 204, 15, 112, 56, 28, 212, 255, 142, 205,
            222, 132, 170, 55, 220, 32, 245, 14, 125, 31, 101, 64, 116, 66, 220, 62, 30, 16, 189,
            160, 201, 66, 12, 108, 78, 145, 37, 170, 121, 17, 51, 153, 64, 252, 3, 85, 161, 119,
            73, 241, 7, 220, 189, 108, 23, 81, 204, 195, 0, 188, 182, 224, 240, 78, 56, 240, 249,
            17, 57, 230, 147, 10, 26, 254, 236, 69, 137, 223, 5, 107, 52, 159, 80, 98, 239, 3, 25,
            78, 148, 244, 213, 192, 225, 203, 15, 230, 71, 181, 113, 129, 232, 132, 144, 30, 85,
            16, 237, 224, 38, 225, 15, 162, 31, 124, 16, 79, 5, 209, 10, 206, 218, 185, 32, 250,
            128, 254, 196, 44, 48, 237, 160, 156, 1, 86, 215, 118, 47, 1, 68, 47, 88, 230, 18, 13,
            86, 255, 206, 62, 16, 172, 126, 241, 77, 43, 48, 205, 243, 90, 7, 108, 126, 6, 183, 29,
            178, 189, 164, 2, 7, 137, 194, 119, 224, 116, 173, 232, 103, 142, 108, 142, 130, 240,
            14, 117, 48, 191, 215, 22, 146, 143, 246, 237, 83, 2, 211, 213, 188, 138, 96, 247, 172,
            61, 5, 142, 207, 7, 94, 102, 95, 64, 250, 21, 228, 22, 49, 144, 227, 142, 134, 84, 45,
            101, 124, 250, 22, 188, 169, 122, 134, 146, 15, 239, 120, 11, 129, 212, 43, 176, 133,
            128, 221, 117, 32, 176, 5, 76, 43, 60, 220, 43, 15, 162, 23, 28, 250, 14, 118, 87, 194,
            161, 62, 1, 48, 45, 113, 5, 226, 143, 192, 108, 21, 48, 253, 219, 215, 16, 44, 190,
            222, 216, 19, 236, 254, 166, 213, 238, 96, 253, 75, 166, 59, 129, 249, 238, 133, 102,
            32, 218, 225, 143, 38, 88, 95, 66, 104, 183, 44, 152, 190, 218, 35, 12, 214, 191, 105,
            6, 56, 158, 27, 122, 83, 89, 145, 221, 237, 48, 125, 19, 7, 88, 29, 255, 69, 118, 176,
            248, 148, 133, 96, 119, 46, 144, 53, 224, 65, 86, 183, 32, 228, 35, 214, 114, 232, 192,
            20, 7, 78, 176, 59, 50, 246, 115, 163, 168, 215, 40, 103, 4, 187, 75, 255, 195, 95, 80,
            56, 60, 88, 86, 2, 46, 23, 96, 229, 81, 67, 67, 235, 111, 228, 240, 105, 56, 27, 13,
            54, 63, 193, 170, 6, 172, 254, 64, 142, 243, 23, 100, 249, 7, 15, 23, 126, 64, 230, 43,
            156, 179, 197, 90, 46, 53, 176, 59, 126, 70, 9, 247, 182, 207, 76, 40, 254, 181, 93,
            248, 199, 30, 0, 253, 128, 181, 173,
        ];

        let expected_decoded_data = vec![
            0.0,
            0.0,
            1.9152265787124634,
            6.3841352462768555,
            0.0,
            4.328546524047852,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.7308083772659302,
            0.0,
            0.5670002698898315,
            0.4630860388278961,
            0.0,
            1.9091640710830688,
            0.0,
            0.0,
            0.0,
            0.0,
            1.9278770685195923,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2.756261110305786,
            0.0,
            6.90661096572876,
            0.0,
            5.97298526763916,
            0.0,
            5.0339508056640625,
            0.0,
            1.5498842000961304,
            1.4491593837738037,
            0.0,
            0.0,
            3.4904677867889404,
            7.235894203186035,
            3.5929534435272217,
            5.555185794830322,
            5.056098937988281,
            0.0,
            2.324875593185425,
            0.4583422541618347,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            5.20412540435791,
            0.0,
            0.0,
            0.0,
            0.7262024879455566,
            0.5922548174858093,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.7419416904449463,
            0.0,
            0.6300871968269348,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            4.258899211883545,
            0.0,
            3.912613868713379,
            3.315601348876953,
            10.774081230163574,
            62.18964385986328,
            168.74815368652344,
            233.35194396972656,
            212.53707885742188,
            193.10165405273438,
            134.95388793945313,
            109.07266998291016,
            68.99895477294922,
            26.85014533996582,
            14.459076881408691,
            30.0445556640625,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.9513055086135864,
            23.305200576782227,
            0.0,
            11.85952091217041,
            6.310108184814453,
            9.372543334960938,
            8.526329040527344,
            5.1979289054870605,
            1.263255000114441,
            2.0791590213775635,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            9.582805633544922,
            0.0,
            0.0,
            0.0,
            0.7025054097175598,
            0.0,
            4.574084281921387,
            6.332054615020752,
            6.379218101501465,
            7.935429096221924,
            5.241953372955322,
            4.139413356781006,
            11.9142484664917,
            10.209602355957031,
            17.304611206054688,
            50.40378952026367,
            47.33993148803711,
            37.184696197509766,
            22.442502975463867,
            11.582002639770508,
            7.386067867279053,
            4.887532711029053,
            4.399116516113281,
            2.674586296081543,
            0.0,
            3.087202548980713,
            2.9770801067352295,
            4.6577911376953125,
            3.523494005203247,
            0.0,
            4.4856743812561035,
            0.0,
            0.0,
            3.156533718109131,
            3.4684603214263916,
            0.0,
            2.183182954788208,
            1.871139407157898,
            2.8069589138031006,
            0.5670002698898315,
            1.7200932502746582,
            0.0,
            4.339651107788086,
            1.8428291082382202,
            1.2664611339569092,
            0.0,
            1.0395220518112183,
            0.0,
            1.8900891542434692,
            0.0,
            0.0,
            1.203376293182373,
            0.0,
            2.3689093589782715,
            0.0,
            1.7893650531768799,
        ];

        match decompress_data(&encoded_data, 1) {
            Ok(decoded_data) => {
                // Assert that the data is decompressed correctly
                assert!(
                    !decoded_data.is_empty(),
                    "Decompressed data should not be empty"
                );
                assert_eq!(
                    decoded_data.len() % 8,
                    0,
                    "Decompressed data should be valid f64 data"
                );
                // println!("Decoded data: {:?}", decoded_data);
                assert_eq!(
                    decoded_data, expected_decoded_data,
                    "Decompressed data does not match expected data"
                );
            }
            Err(e) => {
                panic!("Error decompressing data: {:?}", e);
            }
        }
    }
}
