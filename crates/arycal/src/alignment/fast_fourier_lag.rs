use anyhow::Error as AnyHowError;
use ndarray::Array1;
use rand::prelude::IndexedRandom;
use std::collections::HashMap;
use std::collections::HashSet;
use std::f64;

use super::alignment::calculate_distance;
use super::alignment::construct_mst;
use super::alignment::AlignedChromatogram;
use arycal_cloudpath::sqmass::{Chromatogram, TransitionGroup};
use arycal_common::config::AlignmentConfig;
use arycal_cloudpath::util::extract_basename;

/// Finds the lag with the maximum correlation in a cross-correlation vector.
///
/// # Parameters
/// - `cross_corr`: A vector of cross-correlation values
///
/// # Returns
/// - The lag with the maximum correlation
pub fn find_lag_with_max_correlation(cross_corr: &[f64]) -> isize {
    let (max_idx, _) = cross_corr
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();
    let center = cross_corr.len() / 2;
    max_idx as isize - center as isize
}

/// Shifts the retention times of a chromatogram by a given lag.
///
/// # Parameters
/// - `chrom` A chromatogram to shift
/// - `lag` The lag to apply to the retention times
///
/// # Returns
/// - A new chromatogram with shifted retention times
pub fn shift_chromatogram(chrom: &Chromatogram, lag: isize) -> Chromatogram {
    let mut shifted_chrom = chrom.clone();
    shifted_chrom.retention_times = chrom
        .retention_times
        .iter()
        .map(|&rt| rt + lag as f64) // Apply the lag
        .collect();
    shifted_chrom
}

/// Creates a mapping between the original retention times (RT) of two chromatograms based on the lag.
///
/// # Parameters
/// - `lag`: The lag between the two chromatograms
/// - `chrom1`: The first chromatogram
/// - `chrom2`: The second chromatogram
///
/// # Returns
/// - A vector of hashmaps containing the RT mapping
pub fn create_fft_rt_mapping(
    lag: isize,
    chrom1: &Chromatogram,
    chrom2: &Chromatogram,
) -> Vec<HashMap<String, String>> {
    // Shift the retention times of chrom2 by the lag
    let shifted_chrom2 = shift_chromatogram(chrom2, lag);

    let mut mapping = Vec::new();

    // Iterate over retention times in chrom1
    for (i, &rt1) in shifted_chrom2.retention_times.iter().enumerate() {
        // Calculate the corresponding index in chrom2 using the lag
        // let j = (i as isize + lag) as usize;

        // Check if the index is valid within chrom2
        // if j < chrom2.retention_times.len() {
        let rt2 = chrom2.retention_times[i];

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
        // }
    }

    mapping
}

/// Aligns a series of chromatograms by picking a random run as the reference and aligning all others to it,
/// using FFT-based cross-correlation.
///
/// # Parameters
/// - `smoothed_tics`: A vector of smoothed chromatograms
/// - `params`: Extra alignment configuration parameters
///
/// # Returns
/// - A vector of aligned chromatograms with their alignment offsets
pub fn star_align_tics_fft(
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
        log::trace!(
            "Aligning run: {:?} to reference run",
            chrom.metadata.get("basename").unwrap_or(&chrom.native_id)
        );

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
        let aligned_chrom = shift_chromatogram(chrom, lag.clone());

        let rt_mapping = create_fft_rt_mapping(lag, reference_chrom, chrom);

        // Store the aligned chromatogram
        aligned_chromatograms.push(AlignedChromatogram {
            chromatogram: aligned_chrom.clone(),
            alignment_path: vec![], // No path for FFT-based alignment
            lag: Some(lag), 
            rt_mapping: rt_mapping, 
        });
    }

    Ok(aligned_chromatograms)
}

/// Progressively aligns a series of chromatograms using FFT-based cross-correlation.
/// Instead of aligning all chromatograms to a single reference, it updates the reference progressively.
///
/// # Parameters
/// - smoothed_tics: A vector of smoothed chromatograms
///
/// # Returns
/// - A vector of aligned chromatograms with their alignment offsets
pub fn progressive_align_tics_fft(
    smoothed_tics: Vec<Chromatogram>,
) -> Result<Vec<AlignedChromatogram>, AnyHowError> {
    // Ensure we have at least two chromatograms to align
    if smoothed_tics.len() < 2 {
        return Err(AnyHowError::msg(
            "At least two chromatograms are required for alignment",
        ));
    }

    // Initialize the first chromatogram as the reference
    let mut aligned_sum = smoothed_tics[0].clone();

    // println!(
    //     "Using run: {:?} as the initial reference",
    //     aligned_sum.metadata.get("basename").unwrap_or(&aligned_sum.native_id)
    // );

    let mut aligned_chromatograms = vec![];

    // Progressively align chromatograms
    for (_idx, current_tic) in smoothed_tics.iter().enumerate() {
        // println!(
        //     "Aligning run: {:?} to progressive reference",
        //     current_tic.metadata.get("basename").unwrap_or(&current_tic.native_id)
        // );

        // Step 1: Perform FFT-based cross-correlation
        let cross_corr = fftconvolve::fftcorrelate(
            &Array1::from(aligned_sum.intensities.clone()),
            &Array1::from(current_tic.intensities.clone()),
            fftconvolve::Mode::Full,
        )
        .unwrap()
        .to_vec();

        // Step 2: Find the lag with the maximum correlation
        let lag = find_lag_with_max_correlation(&cross_corr);

        // Step 3: Shift chromatogram retention times by the computed lag
        let aligned_chrom = shift_chromatogram(current_tic, lag);

        // Step 4: Update the progressive reference by averaging the aligned chromatograms
        for (j, intensity) in aligned_chrom.intensities.iter().enumerate() {
            aligned_sum.intensities[j] = (aligned_sum.intensities[j] + intensity) / 2.0;
        }

        let rt_mapping = create_fft_rt_mapping(lag, &aligned_sum, current_tic);

        // Store the aligned chromatogram
        aligned_chromatograms.push(AlignedChromatogram {
            chromatogram: aligned_chrom.clone(),
            alignment_path: vec![], // No path for FFT-based alignment
            lag: Some(lag),
            rt_mapping: rt_mapping, 
        });
    }

    Ok(aligned_chromatograms)
}

/// Aligns a series of chromatograms using a Minimum Spanning Tree (MST) approach
/// with FFT-based cross-correlation.
///
/// # Parameters
/// - `smoothed_tics`: A vector of smoothed chromatograms
///
/// # Returns
/// - A vector of aligned chromatograms with their alignment paths
pub fn mst_align_tics_fft(
    smoothed_tics: Vec<Chromatogram>,
) -> Result<Vec<AlignedChromatogram>, AnyHowError> {
    // Ensure we have at least two chromatograms to align
    if smoothed_tics.len() < 2 {
        return Err(AnyHowError::msg(
            "At least two chromatograms are required for alignment",
        ));
    }

    // Calculate the pairwise distance matrix between chromatograms
    let num_chromatograms = smoothed_tics.len();
    let mut distances = Vec::new();

    for i in 0..num_chromatograms {
        for j in (i + 1)..num_chromatograms {
            // Use a distance function (e.g., DTW or Euclidean distance) to compute dissimilarity
            let distance = calculate_distance(&smoothed_tics[i], &smoothed_tics[j]);
            distances.push((i, j, distance));
        }
    }

    // Generate an MST using Kruskal's algorithm based on the distance matrix
    let mst_edges = construct_mst(&distances, num_chromatograms);

    // Initialize result storage
    let mut aligned_chromatograms = vec![];
    let mut trgrp_aligned = TransitionGroup::new("Aligned MST (FFT)".to_string());

    // HashSet to track chromatograms that have been added to aligned_chromatograms
    let mut aligned_chromatogram_indices = HashSet::new();

    // Align chromatograms based on the edges in the MST
    for (chrom1_idx, chrom2_idx, _) in mst_edges {
        let chrom1 = &smoothed_tics[chrom1_idx];
        let chrom2 = &smoothed_tics[chrom2_idx];

        // println!("Aligning runs: {:?} and {:?}", chrom1.metadata.get("basename"), chrom2.metadata.get("basename"));

        // Step 1: Perform FFT-based cross-correlation
        let cross_corr = fftconvolve::fftcorrelate(
            &Array1::from(chrom1.intensities.clone()),
            &Array1::from(chrom2.intensities.clone()),
            fftconvolve::Mode::Full,
        )
        .unwrap()
        .to_vec();

        let lag = find_lag_with_max_correlation(&cross_corr);

        // Step 2: Shift chromatogram retention times by the computed lag
        let aligned_chrom2 = shift_chromatogram(chrom2, lag);
        let rt_mapping = create_fft_rt_mapping(lag, chrom1, chrom2);

        // Store aligned chromatograms of the first chromatogram if not already added
        if !aligned_chromatogram_indices.contains(&chrom1_idx) {
            aligned_chromatograms.push(AlignedChromatogram {
                chromatogram: chrom1.clone(),
                alignment_path: vec![],
                lag: Some(lag),
                rt_mapping: rt_mapping.clone(),
            });
            aligned_chromatogram_indices.insert(chrom1_idx);
        }

        // Store aligned chromatograms of the second chromatogram
        aligned_chromatograms.push(AlignedChromatogram {
            chromatogram: aligned_chrom2.clone(),
            alignment_path: vec![],
            lag: Some(lag),
            rt_mapping,
        });
        aligned_chromatogram_indices.insert(chrom2_idx);

        // Add to aligned transition group for visualization
        let mut tmp_chrom = aligned_chrom2.clone();
        tmp_chrom.native_id = tmp_chrom
            .metadata
            .get("basename")
            .unwrap_or(&tmp_chrom.native_id)
            .clone();
        trgrp_aligned.add_chromatogram(tmp_chrom);
    }

    Ok(aligned_chromatograms)
}
