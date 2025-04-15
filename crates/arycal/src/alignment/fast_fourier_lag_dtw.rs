use anyhow::Error as AnyHowError;
use ndarray::Array1;
use rand::prelude::IndexedRandom;
use std::collections::HashMap;
use std::f64;
use dtw_rs::{Algorithm, DynamicTimeWarping};

use crate::alignment::fast_fourier_lag::{find_lag_with_max_correlation, shift_chromatogram};
use arycal_common::chromatogram::{Chromatogram, AlignedChromatogram};
use arycal_common::config::AlignmentConfig;
use arycal_cloudpath::util::extract_basename;

/// Creates a mapping between the original retention times (RT) of two chromatograms based on the lag and DTW alignment.
///
/// The mapping is stored as a list of dictionaries, where each dictionary contains the following keys:
///
/// # Parameters
/// - `lag`: The lag between the two chromatograms.
/// - `chrom1`: The reference chromatogram.
/// - `chrom2`: The chromatogram to align.
///
/// # Returns
/// A list of dictionaries containing the following
pub fn create_fft_dtw_rt_mapping(
    lag: isize,
    chrom1: &Chromatogram,
    chrom2: &Chromatogram,
) -> Vec<HashMap<String, String>> {
    let mut mapping = Vec::new();

    // Iterate over retention times in chrom1
    for (i, &rt1) in chrom1.retention_times.iter().enumerate() {
        // Calculate the corresponding index in chrom2 using the lag and DTW alignment
        let j = (i as isize + lag) as usize;

        // Check if the index is valid within chrom2
        if j < chrom2.retention_times.len() {
            let rt2 = chrom2.retention_times[j];

            // Create a mapping entry
            let mut entry = HashMap::new();
            entry.insert("rt1".to_string(), rt1.to_string());
            entry.insert("rt2".to_string(), rt2.to_string());
            entry.insert("alignment".to_string(), format!("({}, {})", rt1, rt2));
            entry.insert(
                "run1".to_string(),
                chrom1
                    .metadata
                    .get("basename")
                    .unwrap_or(&chrom1.native_id)
                    .clone(),
            );
            entry.insert(
                "run2".to_string(),
                chrom2
                    .metadata
                    .get("basename")
                    .unwrap_or(&chrom2.native_id)
                    .clone(),
            );

            mapping.push(entry);
        }
    }

    mapping
}

/// Aligns a series of chromatograms using FFT-based cross-correlation with local refinement via DTW.
///
/// This function aligns a series of chromatograms to a randomly selected reference chromatogram.
/// The alignment is performed in two steps:
///
/// 1. FFT-based cross-correlation to find the lag between the chromatograms.
/// 2. Local refinement using DTW to align the chromatograms based on the computed lag.
///
/// # Parameters
/// - `smoothed_tics`: A vector of chromatograms to align.
/// - `params`: Extra alignment parameters.
///
/// # Returns
/// A vector of aligned chromatograms.
pub fn star_align_tics_fft_with_local_refinement(
    smoothed_tics: Vec<Chromatogram>, params: &AlignmentConfig
) -> Result<Vec<AlignedChromatogram>, AnyHowError> {
    // Ensure we have at least two chromatograms to align
    if smoothed_tics.len() < 2 {
        return Err(AnyHowError::msg(
            "At least two chromatograms are required for alignment",
        ));
    }

    // Randomly pick a reference chromatogram if not specified in params
    let reference_chrom = if let Some(ref_chrom) = &params.reference_run {
        smoothed_tics.iter().find(|x| x.metadata.get("basename").unwrap_or(&x.native_id) == &extract_basename(ref_chrom)).unwrap()
    } else {
        let mut rng = rand::rng();
        let binding = (0..smoothed_tics.len()).collect::<Vec<_>>();
        let reference_idx = binding.choose(&mut rng).unwrap();
        &smoothed_tics[*reference_idx]
    };

    // println!(
    //     "Selected reference run: {:?}",
    //     reference_chrom.metadata.get("basename").unwrap_or(&reference_chrom.native_id)
    // );

    // Initialize storage for aligned chromatograms
    let mut aligned_chromatograms = vec![];

    // Align each chromatogram to the reference
    for (_idx, chrom) in smoothed_tics.iter().enumerate() {
        // if idx == *reference_idx {
        //     // Skip aligning the reference chromatogram to itself
        //     continue;
        // }

        log::trace!(
            "Aligning run: {:?} to reference run",
            chrom.metadata.get("basename").unwrap_or(&chrom.native_id)
        );

        // Step 1: Perform FFT-based cross-correlation
        let cross_corr = fftconvolve::fftcorrelate(
            &Array1::from(reference_chrom.intensities.clone()),
            &Array1::from(chrom.intensities.clone()),
            fftconvolve::Mode::Full,
        )
        .unwrap()
        .to_vec();

        // Step 2: Find the lag with the maximum correlation
        let lag = find_lag_with_max_correlation(&cross_corr);

        // Step 3: Shift chromatogram retention times by the computed lag
        let mut aligned_chrom = shift_chromatogram(chrom, lag);

        // Step 4: Local refinement using DTW
        let (_ref_rt, ref_intensities) = (
            reference_chrom.retention_times.clone(),
            reference_chrom.intensities.clone(),
        );
        let (query_rt, query_intensities) = (
            aligned_chrom.retention_times.clone(),
            aligned_chrom.intensities.clone(),
        );

        // Compute the DTW alignment
        let dtw = DynamicTimeWarping::between(&ref_intensities, &query_intensities);

        // Apply the DTW alignment to the query chromatogram
        let refined_rt: Vec<f64> = dtw.path().iter().map(|&(_, j)| query_rt[j]).collect();

        let refined_intensities: Vec<f64> = dtw
            .path()
            .iter()
            .map(|&(_, j)| query_intensities[j])
            .collect();

        // Update the aligned chromatogram with the refined retention times and intensities
        aligned_chrom.retention_times = refined_rt;
        aligned_chrom.intensities = refined_intensities;

        // Store the aligned chromatogram
        aligned_chromatograms.push(AlignedChromatogram {
            chromatogram: aligned_chrom.clone(),
            alignment_path: dtw.path().to_vec(), // Store the DTW alignment path
            lag: Some(lag),
            rt_mapping: create_fft_dtw_rt_mapping(lag, reference_chrom, &aligned_chrom),
        });
    }

    Ok(aligned_chromatograms)
}
