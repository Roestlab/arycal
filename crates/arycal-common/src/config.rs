use anyhow::Result;
use serde::{Deserialize, Deserializer, Serialize};
use std::path::PathBuf;

#[derive(Debug, Clone, Serialize, PartialEq)]
pub enum XicFileType {
    SqMass,
    Unknown,
}

impl Default for XicFileType {
    fn default() -> Self {
        XicFileType::SqMass
    }
}

impl<'de> Deserialize<'de> for XicFileType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        match s.to_lowercase().as_str() {
            "sqmass" => Ok(XicFileType::SqMass),
            _ => Ok(XicFileType::Unknown),
        }
    }
}

impl XicFileType {
    pub fn as_str(&self) -> &str {
        match self {
            XicFileType::SqMass => "sqMass",
            XicFileType::Unknown => "Unknown",
        }
    }
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub enum FeaturesFileType {
    OSW,
    Unknown,
}

impl Default for FeaturesFileType {
    fn default() -> Self {
        FeaturesFileType::OSW
    }
}

impl<'de> Deserialize<'de> for FeaturesFileType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        match s.to_lowercase().as_str() {
            "osw" => Ok(FeaturesFileType::OSW),
            _ => Ok(FeaturesFileType::Unknown),
        }
    }
}

impl FeaturesFileType {
    pub fn as_str(&self) -> &str {
        match self {
            FeaturesFileType::OSW => "osw",
            FeaturesFileType::Unknown => "Unknown",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct XicConfig {
    #[serde(rename = "include-precursor")]
    pub include_precursor: bool,
    #[serde(rename = "num-isotopes")]
    pub num_isotopes: usize,
    #[serde(rename = "file-type")]
    pub file_type: Option<XicFileType>,
    #[serde(rename = "file-paths")]
    pub file_paths: Vec<PathBuf>,
}

impl XicConfig {
    /// Get the number of file paths.
    pub fn len(&self) -> usize {
        self.file_paths.len()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FeaturesConfig {
    #[serde(rename = "file-type")]
    pub file_type: Option<FeaturesFileType>,
    #[serde(rename = "file-paths")]
    pub file_paths: Vec<PathBuf>,
}

impl FeaturesConfig {
    /// Get the number of file paths.
    pub fn len(&self) -> usize {
        self.file_paths.len()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FiltersConfig {
    pub decoy: bool,
    pub include_identifying_transitions: Option<bool>,
    /// TSV file containing the list of precursors to filter for.
    pub precursor_ids: Option<String>
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SmoothingConfig {
    #[serde(rename = "sgolay_window")]
    pub sgolay_window: usize,
    #[serde(rename = "sgolay_order")]
    pub sgolay_order: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentConfig {
    /// Number of precursors to process in a batch per thread.
    pub batch_size: Option<usize>,
    /// Method to use for alignment. Current options are "FFT", "DTW", "FFTDTW"
    pub method: String,
    /// Type of reference to use for alignment. Current options are "star", "mst", "progressive"
    pub reference_type: String,
    /// Name of the run to use as the reference. If not provided, a run will be used. Only used when the reference type is "star".
    pub reference_run: Option<String>,
    /// Whether to use TIC for alignment. Currently not used. We always use TIC.
    #[serde(rename = "use_tic")]
    pub use_tic: bool,
    /// Smoothing configuration for the chromatograms.
    pub smoothing: SmoothingConfig,
    /// Retention time mapping tolerance in seconds for mapping aligned query peak to reference peak.
    pub rt_mapping_tolerance: Option<f64>,
    /// Method to use for mapping decoy peaks. Current options are "shuffle" and "random_region".
    #[serde(rename = "decoy_peak_mapping_method")]
    pub decoy_peak_mapping_method: String,
    /// Size of the window to use for the decoy peak mapping. Only used when the method is "random_region".
    pub decoy_window_size: Option<usize>,
    /// Optionally compute alignment scores for the full trace alignment and peak mapping. Default is true.
    pub compute_scores: Option<bool>,
    /// Optionally output the scores to a separate OSW file. Otherwise, the scores are added to the input OSW file.
    pub scores_output_file: Option<String>,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        AlignmentConfig {
            batch_size: Some(1000),
            method: "fft_dtw".to_string(),
            reference_type: "star".to_string(),
            reference_run: None,
            use_tic: true,
            smoothing: SmoothingConfig {
                sgolay_window: 11,
                sgolay_order: 3,
            },
            rt_mapping_tolerance: Some(20.0),
            decoy_peak_mapping_method: "shuffle".to_string(),
            decoy_window_size: Some(30),
            compute_scores: Some(true),
            scores_output_file: None,
        }
    }
}

impl std::fmt::Display for AlignmentConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "\n---- Alignment Config ----\n\
            batch_size: {}\n\
            method: {}\n\
            reference_type: {}\n\
            reference_run: {:?}\n\
            use_tic: {}\n\
            sgolay_window: {}\n\
            sgolay_order: {}\n\
            rt_mapping_tolerance: {}\n\
            decoy_peak_mapping_method: {}\n\
            decoy_window_size: {:?}\n\
            compute_scores: {:?}\n\
            scores_output_file: {:?}\n\
            -------------------------",
            self.batch_size.unwrap_or_default(),
            self.method,
            self.reference_type,
            self.reference_run,
            self.use_tic,
            self.smoothing.sgolay_window,
            self.smoothing.sgolay_order,
            self.rt_mapping_tolerance.unwrap_or_default(),
            self.decoy_peak_mapping_method,
            self.decoy_window_size.unwrap_or_default(),
            self.compute_scores.unwrap_or_default(),
            self.scores_output_file
        )
    }
}