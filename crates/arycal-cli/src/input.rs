use anyhow::{Context, Result};
use clap::ArgMatches;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use arycal_common::config::{XicConfig, FeaturesConfig, FiltersConfig, AlignmentConfig, XicFileType, FeaturesFileType, VisualizationConfig};


#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Input {
    pub xic: XicConfig,
    pub features: FeaturesConfig,
    pub filters: FiltersConfig,
    pub alignment: AlignmentConfig,

    pub visualization: Option<VisualizationConfig>, 
}

impl Input {
    /// Load parameters from a JSON file and validate them.
    pub fn from_arguments(matches: &ArgMatches) -> Result<Self> {
        let path = matches
            .get_one::<String>("parameters")
            .expect("required parameters");

        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

        // Infer types if not provided
        input.infer_types()?;

        // Validate the parameters
        input.validate()?;

        // Handle additional command-line arguments for overrides
        if let Some(xic_paths) = matches.get_many::<String>("xic_paths") {
            input.xic.file_paths.extend(xic_paths.map(|p| p.into()));
        }

        log::info!("Loaded parameters from: {}", path);
        log::info!("Features files: {}", input.features.len());
        log::info!("XIC files: {}", input.xic.len());
        if input.xic.len() > 1 && input.features.len() == 1 {
            log::warn!("Multiple XIC files passed and only one feature file passed. Assuming the feature file contains features for all XIC files.");
        }

        // Check if optional filters.include_identifying_transitions is set to true, if it is, ensure alignment.retain_alignment_path is set to true as well if it's not, then set it to true
        if let Some(include_identifying_transitions) = input.filters.include_identifying_transitions {
            if include_identifying_transitions && !input.alignment.retain_alignment_path {
                log::warn!("`filters.include-identifying-transitions` is set to true, but `alignment.retain-alignment-path` is not set. Setting `alignment.retain-alignment-path` to true.");
                input.alignment.retain_alignment_path = true;
            }
        }

        Ok(input)
    }

    /// Load parameters from a JSON file.
    pub fn load(file_path: &str) -> Result<Self> {
        let mut file = File::open(file_path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;

        // Deserialize JSON into Input struct
        let params: Input = serde_json::from_str(&contents).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("JSON parse error: {}", e),
            )
        })?;

        Ok(params)
    }

    /// Infer types from the first file path if not provided.
    fn infer_types(&mut self) -> Result<()> {
        if self.xic.file_type.is_none() && !self.xic.file_paths.is_empty() {
            let first_file = &self.xic.file_paths[0];
            if first_file.ends_with(".sqMass") {
                self.xic.file_type = Some(XicFileType::SqMass);
            } else {
                self.xic.file_type = Some(XicFileType::Unknown);
            }
        }

        if self.features.file_type.is_none() && !self.features.file_paths.is_empty() {
            let first_file = &self.features.file_paths[0];
            if first_file.ends_with(".osw") {
                self.features.file_type = Some(FeaturesFileType::OSW);
            } else {
                self.features.file_type = Some(FeaturesFileType::Unknown);
            }
        }

        Ok(())
    }

    /// Validate the parameters.
    fn validate(&self) -> Result<()> {
        // Validate xic type
        if self.xic.file_type != Some(XicFileType::SqMass) && self.xic.file_type != Some(XicFileType::parquet) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid xic type; expected 'sqMass'",
            )
            .into());
        }

        // Validate features type
        if self.features.file_type != Some(FeaturesFileType::OSW) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Invalid features type {:?}; expected 'osw'",
                    self.features.file_type
                ),
            )
            .into());
        }

        // Validate file paths for xic
        for path in &self.xic.file_paths {
            if !Path::new(path).exists() {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("File not found: {:?}", path),
                )
                .into());
            }
        }

        // Validate file paths for features
        for path in &self.features.file_paths {
            if !Path::new(path).exists() {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("File not found: {:?}", path),
                )
                .into());
            }
        }

        Ok(())
    }
}
