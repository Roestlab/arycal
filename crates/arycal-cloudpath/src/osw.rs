use r2d2::Pool;
use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::{params, Connection, Error as RusqliteError, Result};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, time::Duration};
use std::collections::HashSet;
use std::error::Error;
use std::fmt;
use std::path::Path;
use deepsize::DeepSizeOf;

use arycal_common::{AlignedTransitionScores, FullTraceAlignmentScores, PeakMapping};

/// Define a custom error type
#[derive(Debug)]
pub enum OpenSwathSqliteError {
    DatabaseError(String),
    GeneralError(String),
    RusqliteError(RusqliteError),
    NotFoundError(String),
}

impl fmt::Display for OpenSwathSqliteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OpenSwathSqliteError::DatabaseError(msg) => write!(f, "[OpenSwathSqliteError] Database Error: {}", msg),
            OpenSwathSqliteError::GeneralError(msg) => write!(f, "[OpenSwathSqliteError] Error: {}", msg),
            OpenSwathSqliteError::RusqliteError(err) => write!(f, "[OpenSwathSqliteError] Rusqlite Error: {}", err),
            OpenSwathSqliteError::NotFoundError(msg) => write!(f, "[OpenSwathSqliteError] Not Found Error: {}", msg),
        }
    }
}

/// Implement From for MyError to convert rusqlite errors
impl From<RusqliteError> for OpenSwathSqliteError {
    fn from(err: RusqliteError) -> OpenSwathSqliteError {
        OpenSwathSqliteError::RusqliteError(err)
    }
}

/// Implement std::error::Error for OpenSwathSqliteError
impl Error for OpenSwathSqliteError {}

use std::ops::Deref;

/// Define the ValueEntryType enum to store single or multiple values
#[derive(Debug, Clone, Serialize, Deserialize, DeepSizeOf )]
pub enum ValueEntryType<T> {
    Single(T),
    Multiple(Vec<T>),
}

impl<T: Default> Default for ValueEntryType<T> {
    fn default() -> Self {
        ValueEntryType::Single(T::default())
    }
}

impl<T> Deref for ValueEntryType<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        match self {
            ValueEntryType::Single(ref value) => value,
            ValueEntryType::Multiple(ref vec) => &vec[0],
        }
    }
}

impl<T: Default> ValueEntryType<T> {
    // Method to push a value into the Multiple variant
    pub fn push(&mut self, value: T) {
        match self {
            ValueEntryType::Single(single_value) => {
                // Convert Single to Multiple
                let old_value = std::mem::take(self); // Take the current value
                *self = ValueEntryType::Multiple(vec![old_value.into_single().unwrap(), value]);
            }
            ValueEntryType::Multiple(vec) => {
                vec.push(value);
            }
        }
    }

    // Method to get a reference to the inner value if it's a Single
    pub fn as_single(&self) -> Option<&T> {
        if let ValueEntryType::Single(ref value) = self {
            Some(value)
        } else {
            None
        }
    }

    // Method to get a reference to the inner vector if it's Multiple
    pub fn as_multiple(&self) -> Option<&Vec<T>> {
        if let ValueEntryType::Multiple(ref vec) = self {
            Some(vec)
        } else {
            None
        }
    }

    // Method to convert from Single to Multiple
    pub fn into_single(self) -> Option<T> {
        if let ValueEntryType::Single(value) = self {
            Some(value)
        } else {
            None
        }
    }
}

/// Struct to store feature data for a precursor in a single run i.e. identified peaks, and peak boundaries.
#[derive(Debug, Clone, Serialize, Deserialize, DeepSizeOf )]
pub struct FeatureData {
    pub filename: String,
    pub basename: String,
    pub run_id: i64,
    pub precursor_id: i32,
    pub feature_id: Option<ValueEntryType<i64>>,
    pub exp_rt: ValueEntryType<f64>,
    pub left_width: Option<ValueEntryType<f64>>,
    pub right_width: Option<ValueEntryType<f64>>,
    pub intensity: Option<ValueEntryType<f64>>,
    pub rank: Option<ValueEntryType<i32>>,
    pub qvalue: Option<ValueEntryType<f64>>,

    pub normalized_summed_intensity: Option<ValueEntryType<f64>>,
}

impl FeatureData {
    /// Creates a new FeatureData instance and extracts the basename from the filename.
    pub fn new(
        filename: String,
        run_id: i64,
        precursor_id: i32,
        feature_id: Option<ValueEntryType<i64>>,
        exp_rt: ValueEntryType<f64>,
        left_width: Option<ValueEntryType<f64>>,
        right_width: Option<ValueEntryType<f64>>,
        intensity: Option<ValueEntryType<f64>>,
        rank: Option<ValueEntryType<i32>>,
        qvalue: Option<ValueEntryType<f64>>,
        normalized_summed_intensity: Option<ValueEntryType<f64>>,
    ) -> Self {
        let basename = Self::extract_basename(&filename);

        FeatureData {
            filename,
            basename,
            run_id,
            precursor_id,
            feature_id,
            exp_rt,
            left_width,
            right_width,
            intensity,
            rank,
            qvalue,
            normalized_summed_intensity,
        }
    }

    /// Extracts the basename without any extensions.
    fn extract_basename(filename: &str) -> String {
        // Use Path to get the file stem and remove all extensions
        let path = Path::new(filename);
        let mut stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or_default()
            .to_string();

        // Remove any additional extensions
        while let Some(pos) = stem.rfind('.') {
            stem.truncate(pos); // Remove the last extension
        }

        stem
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, DeepSizeOf )]
pub struct PrecursorIdData {
    pub precursor_id: i32,
    pub unmodified_sequence: String,
    pub modified_sequence: String,
    pub precursor_charge: i32,
    pub transition_ids: Vec<i32>,
    pub identifying_transition_ids: Vec<i32>,
    pub decoy: bool,
}

// Ensure all fields are Send + Sync
unsafe impl Send for PrecursorIdData {}
unsafe impl Sync for PrecursorIdData {}

impl PrecursorIdData {
    pub fn new(
        precursor_id: i32,
        unmodified_sequence: String,
        modified_sequence: String,
        precursor_charge: i32,
        decoy: bool,
    ) -> Self {
        PrecursorIdData {
            precursor_id,
            unmodified_sequence,
            modified_sequence,
            precursor_charge,
            transition_ids: Vec::new(),
            identifying_transition_ids: Vec::new(),
            decoy,
        }
    }

    // Get number of transition ids
    pub fn n_transitions(&self) -> usize {
        self.transition_ids.len()
    }

    // Get number of identifying transition ids
    pub fn n_identifying_transitions(&self) -> usize {
        self.identifying_transition_ids.len()
    }

    /// Method to add a transition ID
    pub fn add_transition(&mut self, transition_id: i32) {
        self.transition_ids.push(transition_id);
    }

    /// Method to add identifying transition IDs
    pub fn add_identifying_transitions(&mut self, transition_ids: i32) {
        self.identifying_transition_ids.push(transition_ids);
    }

    /// Method to extract native IDs for SqMass
    ///
    /// # Parameters
    /// - `include_precursor`: A boolean flag to include precursor IDs
    /// - `n_isotopes`: The number of isotopes to consider for the precursor
    ///
    /// # Returns
    /// A vector of strings containing the native IDs for SqMass
    pub fn extract_native_ids_for_sqmass(
        &self,
        include_precursor: bool,
        n_isotopes: usize,
    ) -> Vec<String> {
        let mut result = Vec::new();

        // Step 1: Add the precursor ID with "_Precursor_i{N}"
        if include_precursor {
            for i in 0..n_isotopes {
                let precursor_string = format!("{}_Precursor_i{}", self.precursor_id, i);
                result.push(precursor_string);
            }
        }

        // Step 2: Add the transition IDs as strings
        for &transition_id in &self.transition_ids {
            result.push(transition_id.to_string());
        }

        result
    }

    /// Method to extract identifying native IDs for SqMass
    ///
    /// # Returns
    /// A vector of strings containing the native IDs for SqMass
    pub fn extract_identifying_native_ids_for_sqmass(&self) -> Vec<String> {
        self.identifying_transition_ids
            .iter()
            .map(|id| id.to_string())
            .collect()
    }
}

/// Define the OSW access structure
#[derive(Clone)]
pub struct OswAccess {
    pool: Pool<SqliteConnectionManager>,
    run_table: HashMap<i64, String>,  // RUN_ID -> basename
    filename_to_id: HashMap<String, i64>,  // basename -> RUN_ID
}

impl OswAccess {
    /// Constructor to create a new OswAccess instance with a connection pool
    /// 
    /// # Parameters
    /// - `db_path`: Path to the OSW database file
    /// - `init_run_table`: Whether to initialize and cache the RUN table
    pub fn new(db_path: &str, init_run_table: bool) -> Result<Self, OpenSwathSqliteError> {
        let manager = SqliteConnectionManager::file(db_path);
        
        let pool = Pool::new(manager)
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Ensure indexes exist
        Self::ensure_indexes(&pool)?;

        // Initialize run tables if requested
        let (run_table, filename_to_id) = if init_run_table {
            Self::load_run_table(&pool).unwrap_or_else(|e| {
                log::warn!("Failed to load RUN table: {}. Using empty tables.", e);
                (HashMap::new(), HashMap::new())
            })
        } else {
            (HashMap::new(), HashMap::new())
        };

        Ok(OswAccess {
            pool,
            run_table,
            filename_to_id,
        })
    }

    /// Verify or create necessary indexes
    fn ensure_indexes(pool: &Pool<SqliteConnectionManager>) -> Result<(), OpenSwathSqliteError> {
        let conn = pool.get().map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        
        // Check which tables exist
        let existing_tables: Vec<String> = conn
        .prepare("SELECT name FROM sqlite_master WHERE type='table'")?
        .query_map([], |row| row.get(0))?
        .collect::<Result<_, _>>()?;

        // Create indexes only for existing tables
        let mut index_queries = Vec::new();

        if existing_tables.iter().any(|t| t == "FEATURE") {
            index_queries.push("CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE(PRECURSOR_ID)");
            index_queries.push("CREATE INDEX IF NOT EXISTS idx_feature_run_id ON FEATURE(RUN_ID)");
        }

        if existing_tables.iter().any(|t| t == "FEATURE_MS2") {
            index_queries.push("CREATE INDEX IF NOT EXISTS idx_feature_ms2_feature_id ON FEATURE_MS2(FEATURE_ID)");
        }

        if existing_tables.iter().any(|t| t == "SCORE_MS2") {
            index_queries.push("CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2(FEATURE_ID)");
        }

        // Execute all index creation queries in a transaction
        let tx = conn.unchecked_transaction()?;
        for query in index_queries {
            if let Err(e) = tx.execute(query, []) {
                log::warn!("Failed to create index: {} - {}", query, e);
            }
        }
        tx.commit()?;
        
        Ok(())
    }

    /// Loads the RUN table into memory
    fn load_run_table(pool: &Pool<SqliteConnectionManager>) -> Result<(HashMap<i64, String>, HashMap<String, i64>), OpenSwathSqliteError> {
        let conn = pool.get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        let query = "SELECT ID, FILENAME FROM RUN";
        let mut stmt = conn.prepare(query)?;

        let mut run_table = HashMap::new();
        let mut filename_to_id = HashMap::new();

        let rows = stmt.query_map([], |row| {
            Ok((row.get::<_, i64>(0)?, row.get::<_, String>(1)?))
        })?;

        for row in rows {
            let (id, filename) = row?;
            let basename = Path::new(&filename)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or(&filename)
                .to_string();
            
            // Handle cases where file_stem() might still have extension (like .gz)
            let basename = if basename.ends_with(".mzML") {
                basename.trim_end_matches(".mzML").to_string()
            } else {
                basename
            };

            run_table.insert(id, basename.clone());
            filename_to_id.insert(basename, id);
        }

        log::trace!("Loaded RUN table with {} entries", run_table.len());

        Ok((run_table, filename_to_id))
    }

    /// Get RUN_IDs for given basenames
    fn get_run_ids(&self, basenames: &[String]) -> Vec<i64> {
        basenames.iter()
            .filter_map(|name| self.filename_to_id.get(name))
            .copied()
            .collect()
    }

    /// Method to fetch precursor id and detecting transition id data from the OSW database
    ///
    /// Parameters
    /// - `filter_decoys`: A boolean flag to filter out decoy precursors.
    /// - `include_identifying_transitions`: A boolean flag to include identifying transitions.
    /// - `precursor_ids`: An optional vector of precursor IDs to filter by.
    ///
    /// Returns
    /// A vector of `PrecursorIdData` instances containing the precursor and transition IDs.
    pub fn fetch_transition_ids(
        &self,
        filter_decoys: bool,
        include_identifying_transitions: bool,
        precursor_ids: Option<Vec<u32>>,
    ) -> Result<Vec<PrecursorIdData>, OpenSwathSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Prepare the base SQL query
        let base_query = r#"
            SELECT 
                PRECURSOR.ID AS PRECURSOR_ID,
                TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID,
                UNMODIFIED_SEQUENCE,
                MODIFIED_SEQUENCE,
                PRECURSOR.CHARGE AS PRECURSOR_CHARGE,
                PRECURSOR.DECOY
            FROM PRECURSOR
            INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
            INNER JOIN PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID 
            INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
            INNER JOIN TRANSITION ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID
            WHERE TRANSITION.DETECTING=1
        "#;

        // Append condition for filtering decoys if applicable
        let query = if filter_decoys {
            format!("{} AND PRECURSOR.DECOY=0", base_query)
        } else {
            base_query.to_string()
        };

        // Append condition for filtering by precursor IDs if applicable
        let query = if let Some(ids) = precursor_ids {
            log::trace!("Filtering for {} precursor IDs", ids.len());
            let id_list = ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join(",");
            format!("{} AND PRECURSOR.ID IN ({})", query, id_list)
        } else {
            query
        };

        // Execute the query and map results to tuples of precursor and transition data
        let mut stmt = conn
            .prepare(&query)
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        let transition_data_iter = stmt
            .query_map([], |row| {
                let precursor_id: i32 = row.get(0)?;
                let transition_id: i32 = row.get(1)?;
                let unmodified_sequence: String = row.get(2)?;
                let modified_sequence: String = row.get(3)?;
                let precursor_charge: i32 = row.get(4)?;
                let decoy: i32 = row.get(5)?;

                // Return a tuple of all the fetched data
                Ok((
                    precursor_id,
                    transition_id,
                    unmodified_sequence,
                    modified_sequence,
                    precursor_charge,
                    decoy == 1,
                ))
            })
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Create a HashMap to collect precursor data
        let mut precursor_map: HashMap<i32, PrecursorIdData> = HashMap::new();

        // Iterate over each row (tuple) and update the HashMap
        for result in transition_data_iter {
            let (
                precursor_id,
                transition_id,
                unmodified_sequence,
                modified_sequence,
                precursor_charge,
                decoy,
            ) = result?;

            // log::trace!(
            //     "Precursor ID: {}, Transition ID: {}, Unmodified Sequence: {}, Modified Sequence: {}, Charge: {}, Decoy: {}",
            //     precursor_id,
            //     transition_id,
            //     unmodified_sequence,
            //     modified_sequence,
            //     precursor_charge,
            //     decoy
            // );

            // Insert into the map or update existing entry
            precursor_map
                .entry(precursor_id)
                .or_insert_with(|| {
                    PrecursorIdData::new(
                        precursor_id,
                        unmodified_sequence.clone(),
                        modified_sequence.clone(),
                        precursor_charge,
                        decoy,
                    )
                })
                .add_transition(transition_id);
        }

        // If identifying transitions are requested, fetch them and add to the PrecursorIdData
        if include_identifying_transitions {
            // Prepare the SQL query
            let identifying_query = r#"
                SELECT 
                    PRECURSOR.ID AS PRECURSOR_ID,
                    TRANSITION.ID AS TRANSITION_ID
                FROM PRECURSOR
                INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
                INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
                INNER JOIN TRANSITION ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID
                WHERE TRANSITION.DETECTING=0
                -- TRANSITION.DECOY=0
                ORDER BY TRANSITION_ID
            "#;

            // Prepare the statement
            let mut stmt = conn
                .prepare(identifying_query)
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

            // Fetch rows based on the query
            let identifying_data_iter = stmt
                .query_map([], |row| {
                    let precursor_id: i32 = row.get(0)?;
                    let transition_id: i32 = row.get(1)?;

                    Ok((precursor_id, transition_id))
                })
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

            // Iterate over each row (tuple) and update the HashMap
            for result in identifying_data_iter {
                let (precursor_id, transition_id) = result?;

                // log::trace!(
                //     "Precursor ID: {}, Identifying Transition ID: {}",
                //     precursor_id,
                //     transition_id
                // );

                // Add the transition ID to the existing PrecursorIdData
                if let Some(ref mut data) = precursor_map.get_mut(&precursor_id) {
                    data.add_identifying_transitions(transition_id);
                }
            }
        }

        // Collect results into a vector from the HashMap
        let precursor_data_vec: Vec<PrecursorIdData> = precursor_map.into_values().collect();

        Ok(precursor_data_vec)
    }

    /// Method to fetch precursor id and transition id data from the OSW database for a specific MODIFIED_SEQUENCE and PRECURSOR_CHARGE
    ///
    /// # Parameters
    /// - `modified_sequence`: A string representing the modified sequence.
    /// - `precursor_charge`: An integer representing the precursor charge.
    ///
    /// # Returns
    /// A single `PrecursorIdData` instance containing the precursor and transition IDs.
    pub fn fetch_detecting_transition_ids_for_sequence(
        &self,
        modified_sequence: &str,
        precursor_charge: i32,
    ) -> Result<PrecursorIdData, OpenSwathSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Prepare the SQL query
        let query = r#"
            SELECT 
                PRECURSOR.ID AS PRECURSOR_ID,
                TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID,
                UNMODIFIED_SEQUENCE,
                MODIFIED_SEQUENCE,
                PRECURSOR.CHARGE AS PRECURSOR_CHARGE,
                PRECURSOR.DECOY
            FROM PRECURSOR
            INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
            INNER JOIN PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID 
            INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
            INNER JOIN TRANSITION ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID
            WHERE PRECURSOR.DECOY=0
            AND TRANSITION.DETECTING=1
            AND MODIFIED_SEQUENCE = ?1
            AND PRECURSOR.CHARGE = ?2
        "#;

        // Prepare the statement
        let mut stmt = conn
            .prepare(query)
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Fetch rows based on the query
        let mut precursor_data: Option<PrecursorIdData> = None;

        let transition_data_iter = stmt
            .query_map([modified_sequence, &precursor_charge.to_string()], |row| {
                let precursor_id: i32 = row.get(0)?;
                let transition_id: i32 = row.get(1)?;
                let unmodified_sequence: String = row.get(2)?;
                let modified_sequence: String = row.get(3)?;
                let precursor_charge: i32 = row.get(4)?;
                let decoy: i32 = row.get(5)?;

                // If this is the first row, initialize the PrecursorIdData struct
                if precursor_data.is_none() {
                    precursor_data = Some(PrecursorIdData::new(
                        precursor_id,
                        unmodified_sequence,
                        modified_sequence,
                        precursor_charge,
                        decoy == 1,
                    ));
                }

                // Add the transition ID to the existing PrecursorIdData
                if let Some(ref mut data) = precursor_data {
                    data.add_transition(transition_id);
                }

                Ok(())
            })
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Ensure the iterator is fully consumed, as query_map is lazy
        for result in transition_data_iter {
            result?; // Check for errors during row processing
        }

        // Return the fetched data or an error if no data was found
        match precursor_data {
            Some(data) => Ok(data),
            None => Err(OpenSwathSqliteError::NotFoundError(format!(
                "No precursor found for modified sequence '{}' and charge '{}'",
                modified_sequence, precursor_charge
            ))),
        }
    }

    /// Method to fetch RT values from OSW database
    ///
    /// This method fetches the filename, precursor ID, and expected retention time (RT) values
    /// from the OSW database. It also allows filtering by decoys, MS2 rank, and maximum q-value.
    ///  
    /// # Parameters
    /// - `filter_decoys`: A boolean flag to filter out decoy precursors.
    /// - `score_ms2_rank`: An integer to filter up to MS2 rank.
    /// - `max_qvalue`: A floating-point value to filter by maximum q-value.
    /// - `runs`: A vector of strings containing the filenames of the runs to filter by.
    ///
    /// # Returns
    /// A vector of `FeatureData` instances containing the filename, precursor ID, and RT values.
    pub fn fetch_feature_data_for_runs(
        &self,
        filter_decoys: bool,
        score_ms2_rank: i32,
        max_qvalue: f64,
        runs: Vec<String>,
    ) -> Result<Vec<FeatureData>, OpenSwathSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Check if SCORE_MS2 table exists
        let table_exists: bool = conn
            .query_row(
                "SELECT COUNT(*) > 0 FROM sqlite_master WHERE type='table' AND name='SCORE_MS2'",
                params![],
                |row| row.get(0),
            )
            .unwrap_or(false);

        if !table_exists {
            return Err(OpenSwathSqliteError::GeneralError(
                "MS2 scoring needs to be performed. The SCORE_MS2 table does not exist."
                    .to_string(),
            ));
        }

        // Start building the SQL query
        let mut sql_query = r#"
            SELECT 
                FILENAME,
                RUN.ID AS RUN_ID,
                FEATURE.PRECURSOR_ID,
                EXP_RT
            FROM FEATURE
            INNER JOIN RUN ON RUN.ID=FEATURE.RUN_ID
            INNER JOIN PRECURSOR ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
            INNER JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
            WHERE 1=1
        "#
        .to_string(); // Convert to String for modification

        // Add filter for decoys
        if filter_decoys {
            sql_query.push_str(" AND PRECURSOR.DECOY = 0");
        }

        // Add filter for SCORE_MS2 rank if specified
        sql_query.push_str(&format!(" AND SCORE_MS2.RANK <= {}", score_ms2_rank));

        // Add filter for SCORE_MS2 QVALUE
        sql_query.push_str(&format!(" AND SCORE_MS2.QVALUE <= {}", max_qvalue));

        // Add filter for specific runs
        let mut run_filter = String::new();
        for run in runs {
            run_filter.push_str(&format!("FILENAME LIKE \"%{}%\" OR ", run));
        }
        run_filter.pop(); // Remove the last space
        run_filter.pop(); // Remove the last O
        run_filter.pop(); // Remove the last R
        run_filter.pop(); // Remove the last space
        sql_query.push_str(&format!(" AND ({})", run_filter));

        // println!("SQL Query: {}", sql_query);

        // Prepare and execute the SQL query
        let mut stmt = conn
            .prepare(&sql_query)
            .map_err(OpenSwathSqliteError::from)?;

        // Execute the query and collect results into a vector of FeatureData
        let feature_data_iter = stmt
            .query_map(params![], |row| {
                let filename: String = row.get(0)?;
                let run_id: i64 = row.get(1)?;
                let precursor_id: i32 = row.get(2)?;
                let exp_rt: f64 = row.get(3)?;

                Ok(FeatureData::new(
                    filename,
                    run_id,
                    precursor_id,
                    None,
                    ValueEntryType::Single(exp_rt),
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                ))
            })
            .map_err(OpenSwathSqliteError::from)?;

        // Collect results into a Vec<FeatureData>
        let feature_data: Vec<FeatureData> = feature_data_iter
            .collect::<Result<Vec<_>, _>>()
            .map_err(OpenSwathSqliteError::from)?;

        Ok(feature_data)
    }

    /// Method to fetch RT values from OSW database
    ///
    /// This method fetches the filename, precursor ID, and expected retention time (RT) values
    /// from the OSW database. It also allows filtering by decoys, MS2 rank, and maximum q-value.
    ///  
    /// # Parameters
    /// - `filter_decoys`: A boolean flag to filter out decoy precursors.
    /// - `score_ms2_rank`: An integer to filter up to MS2 rank.
    /// - `max_qvalue`: A floating-point value to filter by maximum q-value.
    ///
    /// # Returns
    /// A vector of `FeatureData` instances containing the filename, precursor ID, and RT values.
    pub fn fetch_feature_data(
        &self,
        filter_decoys: bool,
        score_ms2_rank: i32,
        max_qvalue: f64,
    ) -> Result<Vec<FeatureData>, OpenSwathSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Check if SCORE_MS2 table exists
        let table_exists: bool = conn
            .query_row(
                "SELECT COUNT(*) > 0 FROM sqlite_master WHERE type='table' AND name='SCORE_MS2'",
                params![],
                |row| row.get(0),
            )
            .unwrap_or(false);

        if !table_exists {
            return Err(OpenSwathSqliteError::GeneralError(
                "MS2 scoring needs to be performed. The SCORE_MS2 table does not exist."
                    .to_string(),
            ));
        }

        // Start building the SQL query
        let mut sql_query = r#"
            SELECT 
                FILENAME,
                RUN.ID AS RUN_ID,
                FEATURE.PRECURSOR_ID,
                EXP_RT
            FROM FEATURE
            INNER JOIN RUN ON RUN.ID=FEATURE.RUN_ID
            INNER JOIN PRECURSOR ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
            INNER JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
            WHERE 1=1
        "#
        .to_string(); // Convert to String for modification

        // Add filter for decoys
        if filter_decoys {
            sql_query.push_str(" AND PRECURSOR.DECOY = 0");
        }

        // Add filter for SCORE_MS2 rank if specified
        sql_query.push_str(&format!(" AND SCORE_MS2.RANK <= {}", score_ms2_rank));

        // Add filter for SCORE_MS2 QVALUE
        sql_query.push_str(&format!(" AND SCORE_MS2.QVALUE <= {}", max_qvalue));

        // Prepare and execute the SQL query
        let mut stmt = conn
            .prepare(&sql_query)
            .map_err(OpenSwathSqliteError::from)?;

        // Execute the query and collect results into a vector of FeatureData
        let feature_data_iter = stmt
            .query_map(params![], |row| {
                let filename: String = row.get(0)?;
                let run_id: i64 = row.get(1)?;
                let precursor_id: i32 = row.get(2)?;
                let exp_rt: f64 = row.get(3)?;

                Ok(FeatureData::new(
                    filename,
                    run_id,
                    precursor_id,
                    None,
                    ValueEntryType::Single(exp_rt),
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                ))
            })
            .map_err(OpenSwathSqliteError::from)?;

        // Collect results into a Vec<FeatureData>
        let feature_data: Vec<FeatureData> = feature_data_iter
            .collect::<Result<Vec<_>, _>>()
            .map_err(OpenSwathSqliteError::from)?;

        Ok(feature_data)
    }

    pub fn fetch_full_precursor_feature_data_for_runs(
        &self,
        precursor_id: i32,
        runs: Vec<String>,
    ) -> Result<Vec<FeatureData>, OpenSwathSqliteError> {
        // Get a connection from the pool
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Convert basenames to RUN_IDs
        let run_ids = self.get_run_ids(&runs);
        log::trace!("RUN_IDs: {:?}", run_ids);

        // Start building the SQL query
        let mut sql_query = r#"
            SELECT 
                FILENAME,
                RUN.ID AS RUN_ID,
                FEATURE.PRECURSOR_ID,
                FEATURE.ID AS FEATURE_ID,
                EXP_RT,
                LEFT_WIDTH,
                RIGHT_WIDTH,
                FEATURE_MS2.AREA_INTENSITY AS INTENSITY
        "#
        .to_string(); // Convert to String for modification

        // Check if SCORE_MS2 table exists
        let score_ms2_exists: bool = conn
            .query_row(
                "SELECT COUNT(*) > 0 FROM sqlite_master WHERE type='table' AND name='SCORE_MS2'",
                params![],
                |row| row.get(0),
            )
            .unwrap_or(false);

        if score_ms2_exists {
            // Add SCORE_MS2 columns to the query
            sql_query.push_str(", SCORE_MS2.RANK, SCORE_MS2.QVALUE");
        }

        sql_query.push_str(
            r#"
            FROM FEATURE
            INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
            INNER JOIN PRECURSOR ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
            INNER JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
        "#,
        );

        if score_ms2_exists {
            sql_query.push_str(" INNER JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID");
        }

        sql_query.push_str(
            r#"
            WHERE 1=1
            AND FEATURE.PRECURSOR_ID = ?1
        "#,
        );

        // Add filter for specific runs using RUN_IDs
        if !run_ids.is_empty() {
            let placeholders = run_ids.iter()
                .map(|_| "?")
                .collect::<Vec<_>>()
                .join(",");
            sql_query.push_str(&format!(" AND RUN.ID IN ({})", placeholders));
        }

        // Order by EXP_RT
        sql_query.push_str(" ORDER BY FILENAME, EXP_RT");

        log::trace!("SQL Query: {}", sql_query);

        // Prepare and execute the SQL query
        let mut stmt = conn
            .prepare(&sql_query)
            .map_err(OpenSwathSqliteError::from)?;

        // Build parameters - precursor_id first, then run_ids
        let mut params = vec![rusqlite::types::Value::from(precursor_id)];
        params.extend(run_ids.iter().map(|&id| rusqlite::types::Value::from(id)));

        // log::debug!("SQL Query: {}", sql_query);

        // Execute the query and collect results into a vector of FeatureData
        let feature_data_list: Vec<(
            String,
            i64,
            i32,
            i64,
            f64,
            f64,
            f64,
            f64,
            Option<i32>,
            Option<f64>,
        )> = stmt
            .query_map(rusqlite::params_from_iter(params.iter()), |row| {  // Use params_from_iter
                let filename: String = row.get(0)?;
                let run_id: i64 = row.get(1)?;
                let precursor_id: i32 = row.get(2)?;
                let feature_id: i64 = row.get(3)?;
                let exp_rt: f64 = row.get(4)?;
                let left_width: f64 = row.get(5)?;
                let right_width: f64 = row.get(6)?;
                let intensity: f64 = row.get(7)?;

                // Optional fields (if SCORE_MS2 exists)
                let rank_option = if score_ms2_exists {
                    Some(row.get::<_, i32>(8)?)
                } else {
                    None
                };

                let qvalue_option = if score_ms2_exists {
                    Some(row.get::<_, f64>(9)?)
                } else {
                    None
                };
                
                Ok((
                    filename,
                    run_id,
                    precursor_id,
                    feature_id,
                    exp_rt,
                    left_width,
                    right_width,
                    intensity,
                    rank_option,
                    qvalue_option,
                ))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        // Collect results into FeatureData structs grouped by filename (run)
        let mut feature_data_map: HashMap<String, FeatureData> = HashMap::new();

        for (
            filename,
            run_id,
            precursor_id,
            feature_id,
            exp_rt,
            left_width,
            right_width,
            intensity,
            rank_option,
            qvalue_option,
        ) in feature_data_list
        {
            let entry = feature_data_map.entry(filename.clone()).or_insert_with(|| {
                FeatureData::new(
                    filename,
                    run_id,
                    precursor_id,
                    Some(ValueEntryType::Multiple(vec![])),
                    ValueEntryType::Multiple(vec![]),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                )
            });

            // Push values into their respective vectors
            if let Some(ValueEntryType::Multiple(ref mut ids)) = entry.feature_id {
                ids.push(feature_id);
            }

            if let ValueEntryType::Multiple(ref mut exps) = entry.exp_rt {
                exps.push(exp_rt);
            }

            if let Some(ValueEntryType::Multiple(ref mut widths)) = entry.left_width {
                widths.push(left_width);
            }

            if let Some(ValueEntryType::Multiple(ref mut widths)) = entry.right_width {
                widths.push(right_width);
            }

            if let Some(ValueEntryType::Multiple(ref mut intensities)) = entry.intensity {
                intensities.push(intensity);
            }

            if let Some(ValueEntryType::Multiple(ref mut ranks)) = entry.rank {
                if let Some(rank) = rank_option {
                    ranks.push(rank);
                }
            }

            if let Some(ValueEntryType::Multiple(ref mut qvalues)) = entry.qvalue {
                if let Some(qvalue) = qvalue_option {
                    qvalues.push(qvalue);
                }
            }
        }

        Ok(feature_data_map.into_values().collect())
    }

    pub fn fetch_feature_data_for_precursor_batch(
        &self,
        precursor_run_sets: &[(i32, Vec<String>)],
    ) -> Result<HashMap<i32, Vec<FeatureData>>, OpenSwathSqliteError> {
        let conn = self.pool.get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        
        // Get all unique precursor IDs
        let precursor_ids: Vec<i32> = precursor_run_sets.iter()
            .map(|(id, _)| *id)
            .collect();
    
        // Get all unique run IDs
        let all_runs: HashSet<String> = precursor_run_sets.iter()
            .flat_map(|(_, runs)| runs.iter().cloned())
            .collect();
        let run_ids = self.get_run_ids(&all_runs.into_iter().collect::<Vec<_>>());
    
        // Build the query
        let mut sql_query = r#"
            SELECT 
                FILENAME,
                RUN.ID AS RUN_ID,
                FEATURE.PRECURSOR_ID,
                FEATURE.ID AS FEATURE_ID,
                EXP_RT,
                LEFT_WIDTH,
                RIGHT_WIDTH,
                FEATURE_MS2.AREA_INTENSITY AS INTENSITY
        "#.to_string();
    
        // Check if SCORE_MS2 table exists
        let score_ms2_exists = conn.query_row(
            "SELECT COUNT(*) > 0 FROM sqlite_master WHERE type='table' AND name='SCORE_MS2'",
            params![],
            |row| row.get(0),
        ).unwrap_or(false);
    
        if score_ms2_exists {
            sql_query.push_str(", SCORE_MS2.RANK, SCORE_MS2.QVALUE");
        }
    
        sql_query.push_str(r#"
            FROM FEATURE
            INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
            INNER JOIN PRECURSOR ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
            INNER JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
        "#);
    
        if score_ms2_exists {
            sql_query.push_str(" INNER JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID");
        }
    
        sql_query.push_str(" WHERE FEATURE.PRECURSOR_ID IN (");
        sql_query.push_str(&precursor_ids.iter().map(|_| "?").collect::<Vec<_>>().join(","));
        sql_query.push_str(")");
    
        if !run_ids.is_empty() {
            sql_query.push_str(" AND RUN.ID IN (");
            sql_query.push_str(&run_ids.iter().map(|_| "?").collect::<Vec<_>>().join(","));
            sql_query.push_str(")");
        }
    
        sql_query.push_str(" ORDER BY FEATURE.PRECURSOR_ID, FILENAME, EXP_RT");
    
        // Prepare and execute
        let mut stmt = conn.prepare(&sql_query)?;
        
        // Build parameters - precursor_ids first, then run_ids
        let mut params: Vec<rusqlite::types::Value> = precursor_ids.iter()
            .map(|&id| rusqlite::types::Value::from(id))
            .collect();
        params.extend(run_ids.iter().map(|&id| rusqlite::types::Value::from(id)));
    
        // Process results
        let mut feature_data_map: HashMap<i32, HashMap<String, FeatureData>> = HashMap::new();
        let rows = stmt.query_map(rusqlite::params_from_iter(params.iter()), |row| {
            let filename: String = row.get(0)?;
            let run_id: i64 = row.get(1)?;
            let precursor_id: i32 = row.get(2)?;
            let feature_id: i64 = row.get(3)?;
            let exp_rt: f64 = row.get(4)?;
            let left_width: f64 = row.get(5)?;
            let right_width: f64 = row.get(6)?;
            let intensity: f64 = row.get(7)?;
    
            // Optional fields (if SCORE_MS2 exists)
            let rank_option = if score_ms2_exists {
                Some(row.get::<_, i32>(8)?)
            } else {
                None
            };
    
            let qvalue_option = if score_ms2_exists {
                Some(row.get::<_, f64>(9)?)
            } else {
                None
            };
            
            Ok((
                precursor_id,
                filename,
                run_id,
                feature_id,
                exp_rt,
                left_width,
                right_width,
                intensity,
                rank_option,
                qvalue_option,
            ))
        })?;
    
        for row in rows {
            let (
                precursor_id,
                filename,
                run_id,
                feature_id,
                exp_rt,
                left_width,
                right_width,
                intensity,
                rank_option,
                qvalue_option,
            ) = row?;
            
            // Get or create the precursor's feature data map
            let precursor_data = feature_data_map.entry(precursor_id)
                .or_default();
            
            // Get or create the FeatureData for this filename
            let feature_data = precursor_data.entry(filename.clone())
                .or_insert_with(|| FeatureData::new(
                    filename,
                    run_id,
                    precursor_id,
                    Some(ValueEntryType::Multiple(vec![])),
                    ValueEntryType::Multiple(vec![]),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                    Some(ValueEntryType::Multiple(vec![])),
                ));
    
            // Push values into their respective vectors
            if let Some(ValueEntryType::Multiple(ref mut ids)) = feature_data.feature_id {
                ids.push(feature_id);
            }
    
            if let ValueEntryType::Multiple(ref mut exps) = feature_data.exp_rt {
                exps.push(exp_rt);
            }
    
            if let Some(ValueEntryType::Multiple(ref mut widths)) = feature_data.left_width {
                widths.push(left_width);
            }
    
            if let Some(ValueEntryType::Multiple(ref mut widths)) = feature_data.right_width {
                widths.push(right_width);
            }
    
            if let Some(ValueEntryType::Multiple(ref mut intensities)) = feature_data.intensity {
                intensities.push(intensity);
            }
    
            if let Some(ValueEntryType::Multiple(ref mut ranks)) = feature_data.rank {
                if let Some(rank) = rank_option {
                    ranks.push(rank);
                }
            }
    
            if let Some(ValueEntryType::Multiple(ref mut qvalues)) = feature_data.qvalue {
                if let Some(qvalue) = qvalue_option {
                    qvalues.push(qvalue);
                }
            }
        }
    
        // Convert to the final HashMap<i32, Vec<FeatureData>> structure
        let result = feature_data_map.into_iter()
            .map(|(precursor_id, data_map)| {
                (precursor_id, data_map.into_values().collect())
            })
            .collect();
    
        Ok(result)
    }

    /// Create the FEATURE_ALIGNMENT table if it doesn't exist
    pub fn create_feature_alignment_table(&self) -> Result<(), OpenSwathSqliteError> {
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Check if the table exists
        let table_exists: bool = conn.query_row(
            "SELECT EXISTS (SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'FEATURE_ALIGNMENT');",
            [],
            |row| row.get(0),
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // If the table exists, drop it and log a warning
        if table_exists {
            log::warn!("Table FEATURE_ALIGNMENT seems to already exist. Dropping it to create a new table for incoming data.");
            conn.execute(
                "DROP TABLE FEATURE_ALIGNMENT;",
                [],
            )
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        }
        
        conn.execute(
            r#"
                CREATE TABLE IF NOT EXISTS FEATURE_ALIGNMENT (
                    REFERENCE_FILENAME TEXT,
                    ALIGNED_FILENAME TEXT,
                    XCORR_COELUTION_TO_REFERENCE REAL,
                    XCORR_SHAPE_TO_REFERENCE REAL,
                    MI_TO_REFERENCE REAL,
                    XCORR_COELUTION_TO_ALL REAL,
                    XCORR_SHAPE_TO_ALL REAL,
                    MI_TO_ALL REAL
                );
                "#,
            [],
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // // Create indices for faster lookups
        // conn.execute(
        //     r#"
        //     CREATE INDEX IF NOT EXISTS idx_alignment_id ON FEATURE_ALIGNMENT (ALIGNMENT_ID);
        //     CREATE INDEX IF NOT EXISTS idx_reference_feature_id ON FEATURE_ALIGNMENT (REFERENCE_FEATURE_ID);
        //     CREATE INDEX IF NOT EXISTS idx_aligned_feature_id ON FEATURE_ALIGNMENT (ALIGNED_FEATURE_ID);
        //     "#,
        //     [],
        // )
        // .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        Ok(())
    }

    /// Insert batch of feature alignment data into the FEATURE_ALIGNMENT table
    pub fn insert_feature_alignment_batch(
        &self,
        scores: &Vec<&FullTraceAlignmentScores>,
    ) -> Result<(), OpenSwathSqliteError> {
        {
            // Get a connection from the pool
            let mut conn = self
                .pool
                .get()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
            // Begin a transaction
            let tx = conn
                .transaction()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
            {
                let mut stmt = tx
                    .prepare(
                        r#"
                        INSERT INTO FEATURE_ALIGNMENT (
                            reference_filename, aligned_filename,
                            xcorr_coelution_to_reference, xcorr_shape_to_reference, mi_to_reference,
                            xcorr_coelution_to_all, xcorr_shape_to_all, mi_to_all
                        ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)
                        "#,
                    )
                    .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
                for peak_mapping in scores {
                    stmt.execute(params![
                        peak_mapping.reference_filename,
                        peak_mapping.aligned_filename,
                        peak_mapping.xcorr_coelution_to_ref,
                        peak_mapping.xcorr_shape_to_ref,
                        peak_mapping.mi_to_ref,
                        peak_mapping.xcorr_coelution_to_all,
                        peak_mapping.xcorr_shape_to_all,
                        peak_mapping.mi_to_all
                    ])
                    .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
                }
            } // `stmt` is dropped here
    
            tx.commit()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        } // `tx` and `conn` are dropped here, returning the connection to the pool
    
        Ok(())
    }
    

    /// Create the FEATURE_MS2_ALIGNMENT table if it doesn't exist
    pub fn create_feature_ms2_alignment_table(&self) -> Result<(), OpenSwathSqliteError> {
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Check if the table exists
        let table_exists: bool = conn.query_row(
            "SELECT EXISTS (SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'FEATURE_MS2_ALIGNMENT');",
            [],
            |row| row.get(0),
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // If the table exists, drop it and log a warning
        if table_exists {
            log::warn!("Table FEATURE_MS2_ALIGNMENT seems to already exist. Dropping it to create a new table for incoming data.");
            conn.execute(
                "DROP TABLE FEATURE_MS2_ALIGNMENT;",
                [],
            )
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        }

        conn.execute(
            r#"
            CREATE TABLE IF NOT EXISTS FEATURE_MS2_ALIGNMENT (
                ALIGNMENT_ID INTEGER,
                PRECURSOR_ID INTEGER,
                RUN_ID INTEGER,
                REFERENCE_FEATURE_ID INTEGER,
                ALIGNED_FEATURE_ID INTEGER,
                REFERENCE_RT REAL,
                ALIGNED_RT REAL,
                REFERENCE_LEFT_WIDTH REAL,
                REFERENCE_RIGHT_WIDTH REAL,
                ALIGNED_LEFT_WIDTH REAL,
                ALIGNED_RIGHT_WIDTH REAL,
                REFERENCE_FILENAME TEXT,
                ALIGNED_FILENAME TEXT,
                XCORR_COELUTION_TO_REFERENCE REAL,
                XCORR_SHAPE_TO_REFERENCE REAL,
                MI_TO_REFERENCE REAL,
                XCORR_COELUTION_TO_ALL REAL,
                XCORR_SHAPE_TO_ALL REAL,
                MI_TO_ALL REAL,
                RETENTION_TIME_DEVIATION REAL,
                PEAK_INTENSITY_RATIO REAL,
                LABEL INTEGER
            );
            "#,
            [],
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Create indices for faster lookups
        conn.execute(
            r#"
            CREATE INDEX IF NOT EXISTS idx_alignment_id ON FEATURE_MS2_ALIGNMENT (ALIGNMENT_ID);
            CREATE INDEX IF NOT EXISTS idx_reference_feature_id ON FEATURE_MS2_ALIGNMENT (REFERENCE_FEATURE_ID);
            CREATE INDEX IF NOT EXISTS idx_aligned_feature_id ON FEATURE_MS2_ALIGNMENT (ALIGNED_FEATURE_ID);
            CREATE INDEX IF NOT EXISTS idx_precursor_id ON FEATURE_MS2_ALIGNMENT (PRECURSOR_ID);
            CREATE INDEX IF NOT EXISTS idx_run_id ON FEATURE_MS2_ALIGNMENT (RUN_ID);
            "#,
            [],
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        Ok(())
    }

    /// Insert batch of feature MS2 alignment data into the FEATURE_MS2_ALIGNMENT table
    pub fn insert_feature_ms2_alignment_batch(
        &self,
        peak_mappings: &[PeakMapping], // Accepts a slice of PeakMapping for batch insertion
    ) -> Result<(), OpenSwathSqliteError> {
        {
            // Get a connection from the pool
            let mut conn = self
                .pool
                .get()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
            // Begin a transaction
            let tx = conn
                .transaction()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
            {
                let mut stmt = tx
                    .prepare(
                        r#"
                        INSERT INTO FEATURE_MS2_ALIGNMENT (
                            alignment_id, precursor_id, run_id, reference_feature_id, aligned_feature_id,
                            reference_rt, aligned_rt, reference_left_width, reference_right_width,
                            aligned_left_width, aligned_right_width, reference_filename, aligned_filename,
                            xcorr_coelution_to_reference, xcorr_shape_to_reference, mi_to_reference,
                            xcorr_coelution_to_all, xcorr_shape_to_all, mi_to_all,
                            retention_time_deviation, peak_intensity_ratio, label
                        ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16, ?17, ?18, ?19, ?20, ?21, ?22)
                        "#,
                    )
                    .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
                for peak_mapping in peak_mappings {
                    stmt.execute(params![
                        peak_mapping.alignment_id,
                        peak_mapping.precursor_id,
                        peak_mapping.run_id,
                        peak_mapping.reference_feature_id,
                        peak_mapping.aligned_feature_id,
                        peak_mapping.reference_rt,
                        peak_mapping.aligned_rt,
                        peak_mapping.reference_left_width,
                        peak_mapping.reference_right_width,
                        peak_mapping.aligned_left_width,
                        peak_mapping.aligned_right_width,
                        peak_mapping.reference_filename,
                        peak_mapping.aligned_filename,
                        peak_mapping.xcorr_coelution_to_ref,
                        peak_mapping.xcorr_shape_to_ref,
                        peak_mapping.mi_to_ref,
                        peak_mapping.xcorr_coelution_to_all,
                        peak_mapping.xcorr_shape_to_all,
                        peak_mapping.mi_to_all,
                        peak_mapping.rt_deviation,
                        peak_mapping.intensity_ratio,
                        peak_mapping.label,
                    ])
                    .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
                }
            } // `stmt` goes out of scope here and is dropped
    
            tx.commit()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        } // `tx` and `conn` go out of scope here, ensuring they are returned to the pool
    
        Ok(())
    }
    

    /// Create the FEATURE_TRANSITION_ALIGNMENT table if it doesn't exist
    pub fn create_feature_transition_alignment_table(&self) -> Result<(), OpenSwathSqliteError> {
        let conn = self
            .pool
            .get()
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Check if the table exists
        let table_exists: bool = conn.query_row(
            "SELECT EXISTS (SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'FEATURE_TRANSITION_ALIGNMENT');",
            [],
            |row| row.get(0),
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // If the table exists, drop it and log a warning
        if table_exists {
            log::warn!("Table FEATURE_TRANSITION_ALIGNMENT seems to already exist. Dropping it to create a new table for incoming data.");
            conn.execute(
                "DROP TABLE FEATURE_TRANSITION_ALIGNMENT;",
                [],
            )
            .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        }

        conn.execute(
            r#"
            CREATE TABLE IF NOT EXISTS FEATURE_TRANSITION_ALIGNMENT (
                FEATURE_ID INTEGER,
                TRANSITION_ID INTEGER,
                LABEL INTEGER,
                XCORR_COELUTION_TO_REFERENCE REAL,
                XCORR_SHAPE_TO_REFERENCE REAL,
                MI_TO_REFERENCE REAL,
                XCORR_COELUTION_TO_ALL REAL,
                XCORR_SHAPE_TO_ALL REAL,
                MI_TO_ALL REAL,
                RETENTION_TIME_DEVIATION REAL,
                PEAK_INTENSITY_RATIO REAL
            );
            "#,
            [],
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        // Create indices for faster lookups
        conn.execute(
            r#"
            CREATE INDEX IF NOT EXISTS idx_alignment_feature_id ON FEATURE_TRANSITION_ALIGNMENT (FEATURE_ID);
            CREATE INDEX IF NOT EXISTS idx_alignment_transition ON FEATURE_TRANSITION_ALIGNMENT (TRANSITION_ID);
            "#,
            [],
        )
        .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;

        Ok(())
    }

    /// Insert batch of feature transition alignment data into the FEATURE_TRANSITION_ALIGNMENT table
    pub fn insert_feature_transition_alignment_batch(
        &self,
        peak_mappings: &[AlignedTransitionScores],
    ) -> Result<(), OpenSwathSqliteError> {
        {
            // Get a connection from the pool
            let mut conn = self
                .pool
                .get()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
            // Begin a transaction
            let tx = conn
                .transaction()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
            {
                let mut stmt = tx
                    .prepare(
                        r#"
                        INSERT INTO FEATURE_TRANSITION_ALIGNMENT (
                            feature_id, transition_id, label,
                            xcorr_coelution_to_reference, xcorr_shape_to_reference, mi_to_reference,
                            xcorr_coelution_to_all, xcorr_shape_to_all, mi_to_all,
                            retention_time_deviation, peak_intensity_ratio
                        ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)
                        "#,
                    )
                    .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
    
                for peak_mapping in peak_mappings {
                    stmt.execute(params![
                        peak_mapping.feature_id,
                        peak_mapping.transition_id,
                        peak_mapping.label,
                        peak_mapping.xcorr_coelution_to_ref,
                        peak_mapping.xcorr_shape_to_ref,
                        peak_mapping.mi_to_ref,
                        peak_mapping.xcorr_coelution_to_all,
                        peak_mapping.xcorr_shape_to_all,
                        peak_mapping.mi_to_all,
                        peak_mapping.rt_deviation,
                        peak_mapping.intensity_ratio,
                    ])
                    .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
                }
            } // `stmt` is dropped here
    
            tx.commit()
                .map_err(|e| OpenSwathSqliteError::DatabaseError(e.to_string()))?;
        } // `tx` and `conn` are dropped here, returning the connection to the pool
    
        Ok(())
    }
    
}
