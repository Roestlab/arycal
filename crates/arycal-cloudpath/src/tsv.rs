use std::fs::File;
use anyhow::Error;
use csv::ReaderBuilder;

/// Load precursor IDs from a TSV file
/// 
/// TSV file should have a single column with precursor IDs
/// 
/// # Example
/// 
/// precursor_id
/// 1
/// 2
/// 3
/// 
/// # Arguments
/// 
/// * `file_path` - A string slice that holds the path to the TSV file
/// 
/// # Returns
/// 
/// A vector of u32 values that represent the precursor IDs
pub fn load_precursor_ids_from_tsv(file_path: &str) -> Result<Vec<u32>, Error> {
    let file = File::open(file_path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t') 
        .has_headers(true) 
        .from_reader(file);

    let mut precursor_ids = Vec::new();

    for result in rdr.records() {
        let record = result?;
        let precursor_id: u32 = record[0].parse()?; // Parse the first column as u32
        precursor_ids.push(precursor_id);
    }

    Ok(precursor_ids)
}