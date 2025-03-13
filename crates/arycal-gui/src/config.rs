// arycal-gui/src/config.rs
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct XicConfig {
    #[serde(rename = "include-precursor")]
    pub include_precursor: bool,
    #[serde(rename = "num-isotopes")]
    pub num_isotopes: usize,
    pub r#type: String,
    #[serde(rename = "file-paths")]
    pub file_paths: Vec<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FeaturesConfig {
    pub r#type: String,
    #[serde(rename = "file-paths")]
    pub file_paths: Vec<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FiltersConfig {
    pub decoy: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SmoothingConfig {
    #[serde(rename = "sgolay_window")]
    pub sgolay_window: usize,
    #[serde(rename = "sgolay_order")]
    pub sgolay_order: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AlignmentConfig {
    pub method: String,
    pub reference_type: String,
    #[serde(rename = "use_tic")]
    pub use_tic: bool,
    pub smoothing: SmoothingConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AppConfig {
    pub xic: XicConfig,
    pub features: FeaturesConfig,
    pub filters: FiltersConfig,
    pub alignment: AlignmentConfig,
}