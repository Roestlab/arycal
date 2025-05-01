use duckdb::{params, Connection, Result};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use crate::osw::PrecursorIdData;
use crate::sqmass::{decompress_data, TransitionGroup};
use crate::util::extract_basename;
use crate::ChromatogramReader;
use arycal_common::chromatogram::Chromatogram;

use duckdb::{self, Error as DuckDBError};
use std::error::Error;
use std::fmt;

/// Define a custom error type for DuckDB Parquet access
#[derive(Debug)]
pub enum DuckDBParquetError {
    DatabaseError(String),
    GeneralError(String),
    DuckDBError(DuckDBError),
    DecompressionError(String),
}

impl fmt::Display for DuckDBParquetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DuckDBParquetError::DatabaseError(msg) => {
                write!(f, "[DuckDBParquetError] Database Error: {}", msg)
            }
            DuckDBParquetError::GeneralError(msg) => {
                write!(f, "[DuckDBParquetError] Error: {}", msg)
            }
            DuckDBParquetError::DuckDBError(err) => {
                write!(f, "[DuckDBParquetError] DuckDB Error: {}", err)
            }
            DuckDBParquetError::DecompressionError(msg) => {
                write!(f, "[DuckDBParquetError] Decompression Error: {}", msg)
            }
        }
    }
}

/// Implement From for DuckDBParquetError to convert DuckDB errors
impl From<DuckDBError> for DuckDBParquetError {
    fn from(err: DuckDBError) -> DuckDBParquetError {
        DuckDBParquetError::DuckDBError(err)
    }
}

/// Implement std::error::Error for DuckDBParquetError
impl Error for DuckDBParquetError {}
unsafe impl Send for DuckDBParquetError {}
unsafe impl Sync for DuckDBParquetError {}

/// Reader for chromatograms stored in Parquet format using DuckDB
pub struct DuckDBParquetChromatogramReader {
    file: String,
    conn: Connection,
}

unsafe impl Send for DuckDBParquetChromatogramReader {}
unsafe impl Sync for DuckDBParquetChromatogramReader {}

impl ChromatogramReader for DuckDBParquetChromatogramReader {
    /// Create a new reader for the given Parquet file
    fn new(db_path: &str) -> anyhow::Result<Self> {
        log::trace!(
            "Creating DuckDBParquetChromatogramReader instance with path: {}",
            db_path
        );

        let conn = Connection::open_in_memory().map_err(|e| {
            DuckDBParquetError::DatabaseError(format!("Failed to create DuckDB connection: {}", e))
        })?;

        // // Enable Parquet extension
        // conn.execute("INSTALL 'parquet';", []).map_err(|e| {
        //     DuckDBParquetError::DatabaseError(format!("Failed to install parquet extension: {}", e))
        // })?;
        // conn.execute("LOAD 'parquet';", []).map_err(|e| {
        //     DuckDBParquetError::DatabaseError(format!("Failed to load parquet extension: {}", e))
        // })?;

        // let parquet_supported: bool = conn
        //     .query_row(
        //         "SELECT count(*) > 0 FROM duckdb_extensions() WHERE extension_name = 'parquet'", 
        //         [], 
        //         |row| row.get(0)
        //     )
        //     .unwrap_or(false);
        
        // if !parquet_supported {
        //     // Try loading explicitly if not detected
        //     if let Err(e) = conn.execute("LOAD 'parquet'", []) {
        //         return Err(anyhow::anyhow!(
        //             "Parquet support not available and failed to load: {}", e
        //         ));
        //     }
        // }

        Ok(DuckDBParquetChromatogramReader {
            file: db_path.to_string(),
            conn,
        })
    }

    /// Reads chromatograms from the Parquet file based on specified filter criteria.
    ///
    /// Note: Only supports "NATIVE_ID" filter type (unlike SQLite version which also supports "CHROMATOGRAM_ID")
    ///
    /// # Parameters
    /// - `filter_type`: Must be "NATIVE_ID" (other values will return error)
    /// - `filter_values`: A vector of string slices representing the native IDs to filter by
    /// - `group_id`: A string representing the group ID for the transition group
    ///
    /// # Returns
    /// - A `Result` containing a `TransitionGroup` on success,
    ///   or an error of type `DuckDBParquetError` on failure
    fn read_chromatograms(
        &self,
        filter_type: &str,
        filter_values: Vec<&str>,
        group_id: String,
    ) -> anyhow::Result<TransitionGroup> {
        // Validate filter type (only NATIVE_ID supported for Parquet)
        if filter_type != "NATIVE_ID" {
            return Err(DuckDBParquetError::GeneralError(
                "Only NATIVE_ID filter type is supported for Parquet files".to_string(),
            ).into());
        }

        // Create placeholders for the query
        let placeholders = filter_values
            .iter()
            .map(|id| format!("'{}'", id))
            .collect::<Vec<_>>()
            .join(",");

        // Query the Parquet file directly using DuckDB
        let query = format!(
            "SELECT 
                NATIVE_ID,
                RT_DATA, 
                INTENSITY_DATA, 
                RT_COMPRESSION, 
                INTENSITY_COMPRESSION 
             FROM read_parquet('{}') 
             WHERE NATIVE_ID IN ({})",
            self.file, placeholders
        );

        let mut stmt = self.conn.prepare(&query)?;

        let mut transition_group = TransitionGroup::new(group_id);

        // Add file metadata
        transition_group.add_metadata("file".to_string(), self.file.to_string());
        let basename = extract_basename(&self.file);
        transition_group.add_metadata("basename".to_string(), basename);

        let rows = stmt.query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,  // native_id
                row.get::<_, Vec<u8>>(1)?, // rt_data
                row.get::<_, Vec<u8>>(2)?, // intensity_data
                row.get::<_, i32>(3)?,     // rt_compression
                row.get::<_, i32>(4)?,     // intensity_compression
            ))
        })?;

        for row in rows {
            let (native_id, rt_data, intensity_data, rt_compression, intensity_compression) = row?;

            // Decode retention times
            let retention_times = decompress_data(&rt_data, rt_compression)
                .map_err(|e| DuckDBParquetError::DecompressionError(e.to_string()))?;

            // Decode intensities
            let intensities = decompress_data(&intensity_data, intensity_compression)
                .map_err(|e| DuckDBParquetError::DecompressionError(e.to_string()))?;

            // Create or update chromatogram
            let chromatogram = transition_group
                .chromatograms
                .entry(native_id.clone())
                .or_insert(Chromatogram {
                    id: 0, // No chromatogram ID in Parquet format
                    native_id,
                    retention_times: Vec::new(),
                    intensities: Vec::new(),
                    metadata: HashMap::new(),
                });

            chromatogram.retention_times = retention_times;
            chromatogram.intensities = intensities;
        }

        Ok(transition_group)
    }

    /// Read chromatograms for the given precursors
    fn read_chromatograms_for_precursors(
        &self,
        precursors: &[PrecursorIdData],
        include_precursor: bool,
        num_isotopes: usize,
    ) -> anyhow::Result<HashMap<i32, TransitionGroup>> {
        // Build native ID to precursor mapping
        let precursor_to_native_ids: HashMap<i32, Vec<String>> = precursors
            .iter()
            .map(|precursor| {
                let native_ids =
                    precursor.extract_native_ids_for_sqmass(include_precursor, num_isotopes);
                (precursor.precursor_id, native_ids)
            })
            .collect();

        // Collect all native IDs we need to load
        let all_native_ids: Vec<String> = precursor_to_native_ids
            .values()
            .flatten()
            .cloned()
            .collect();

        if all_native_ids.is_empty() {
            return Ok(HashMap::new());
        }

        // Create placeholders for the query
        let placeholders = all_native_ids
            .iter()
            .map(|id| format!("'{}'", id))
            .collect::<Vec<_>>()
            .join(",");

        // Query the Parquet file directly using DuckDB
        let query = format!(
            "SELECT 
                PRECURSOR_ID, 
                NATIVE_ID, 
                RT_DATA, 
                INTENSITY_DATA, 
                RT_COMPRESSION, 
                INTENSITY_COMPRESSION 
             FROM read_parquet('{}') 
             WHERE NATIVE_ID IN ({})",
            self.file, placeholders
        );

        let mut stmt = self.conn.prepare(&query)?;

        let rows = stmt.query_map([], |row| {
            Ok((
                row.get::<_, i32>(0)?,     // precursor_id
                row.get::<_, String>(1)?,  // native_id
                row.get::<_, Vec<u8>>(2)?, // rt_data
                row.get::<_, Vec<u8>>(3)?, // intensity_data
                row.get::<_, i32>(4)?,     // rt_compression
                row.get::<_, i32>(5)?,     // intensity_compression
            ))
        })?;

        let mut precursor_groups = HashMap::new();

        for row in rows {
            let (
                precursor_id,
                native_id,
                rt_data,
                intensity_data,
                rt_compression,
                intensity_compression,
            ) = row?;

            // Find which precursor this chromatogram belongs to
            if let Some(precursor_native_ids) = precursor_to_native_ids.get(&precursor_id) {
                if precursor_native_ids.contains(&native_id) {
                    let group_id = format!(
                        "{}_{}",
                        precursors
                            .iter()
                            .find(|p| p.precursor_id == precursor_id)
                            .map(|p| p.modified_sequence.clone())
                            .unwrap_or_default(),
                        precursors
                            .iter()
                            .find(|p| p.precursor_id == precursor_id)
                            .map(|p| p.precursor_charge.to_string())
                            .unwrap_or_default()
                    );

                    let group =
                        precursor_groups
                            .entry(precursor_id)
                            .or_insert_with(|| TransitionGroup {
                                group_id: group_id.clone(),
                                chromatograms: HashMap::new(),
                                metadata: HashMap::new(),
                            });

                    let chrom = group
                        .chromatograms
                        .entry(native_id.clone())
                        .or_insert_with(|| Chromatogram {
                            id: precursor_id, // Using precursor_id as chromatogram ID
                            native_id: native_id.clone(),
                            retention_times: Vec::new(),
                            intensities: Vec::new(),
                            metadata: HashMap::new(),
                        });

                    // Decode retention times
                    chrom.retention_times =
                        decompress_data(&rt_data, rt_compression).map_err(|e| {
                            DuckDBParquetError::DecompressionError(format!(
                                "RT decompression failed: {}",
                                e
                            ))
                        })?;

                    // Decode intensities
                    chrom.intensities = decompress_data(&intensity_data, intensity_compression)
                        .map_err(|e| {
                            DuckDBParquetError::DecompressionError(format!(
                                "Intensity decompression failed: {}",
                                e
                            ))
                        })?;
                }
            }
        }

        // Add metadata to each group
        for group in precursor_groups.values_mut() {
            group
                .metadata
                .insert("file".to_string(), self.file.to_string());
            group
                .metadata
                .insert("basename".to_string(), extract_basename(&self.file));
        }

        Ok(precursor_groups)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_read_chromatograms() {
        let db_path = "/home/singjc/Documents/github/arycal/tests/synth_phospho/xics/chludwig_K150309_013_SW_0.chrom.parquet";

        let precursors = vec![PrecursorIdData {
            precursor_id: 29,
            unmodified_sequence: "ANSSPTTNIDHLK".to_string(),
            modified_sequence: "ANS(UniMod:21)SPTTNIDHLK(UniMod:259)".to_string(),
            precursor_charge: 2,
            transition_ids: vec![174, 175, 176, 177, 178, 179],
            identifying_transition_ids: Vec::new(),
            decoy: false,
        }];

        let reader =
            DuckDBParquetChromatogramReader::new(db_path).expect("Failed to create reader");

        let chromatograms = reader.read_chromatograms_for_precursors(&precursors, true, 1);
        println!("Chromatograms: {:?}", chromatograms);
    }
}
