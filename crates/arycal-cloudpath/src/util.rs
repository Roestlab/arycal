use std::path::Path;

/// Extracts the basename without any extensions.
pub fn extract_basename(filename: &str) -> String {
    // Use Path to get the file stem and remove all extensions
    let path = Path::new(filename);
    let mut stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or_default().to_string();

    // Remove any additional extensions
    while let Some(pos) = stem.rfind('.') {
        stem.truncate(pos); // Remove the last extension
    }

    stem
}