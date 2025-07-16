use arycal_cloudpath::osw::{OswAccess, OpenSwathSqliteError};



fn main() -> Result<(), OpenSwathSqliteError> {
    // Specify the path to your SQLite database file
    let db_path = "/home/singjc/Documents/github/ptms_align/data/fragpipe_oswbench_20220512/data/osw/merged.osw";

    // Create an instance of OswAccess
    let osw_access = OswAccess::new(db_path)?;

    // Fetch feature data from the database
    match osw_access.fetch_feature_data(true, 1, 0.01) {
        Ok(data) => {
            for feature in data.iter().take(5) {
                println!("{:?}", feature);
            }
        }
        Err(e) => eprintln!("Error fetching data: {}", e),
    }

    let prec_tr_ids = osw_access.fetch_transition_ids(true, false)?;
    println!("Precursor ID data length: {}", prec_tr_ids.len());
    println!("Precursor ID data first 5: {:?}", prec_tr_ids.iter().take(5).collect::<Vec<_>>());

    let peptide_sequence = "ANS(UniMod:21)SPTTNIDHLK(UniMod:259)";
    let charge = 3;
    let peptide_charge = format!("{}_{}", peptide_sequence, charge);

    let tgt_prec = osw_access.fetch_detecting_transition_ids_for_sequence(peptide_sequence, charge)?;

    println!("Target precursor IDs: {:?}", tgt_prec);

    Ok(())
}