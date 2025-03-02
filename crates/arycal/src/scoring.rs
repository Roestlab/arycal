use rayon::prelude::*;
use ndarray::Array1;
use std::collections::HashMap;
use rand::{seq::SliceRandom, Rng};

use arycal_cloudpath::sqmass::{Chromatogram, TransitionGroup};
use arycal_common::{AlignedTransitionScores, FullTraceAlignmentScores, PeakMapping};
use crate::{alignment::alignment::AlignedChromatogram, stats::{calc_mi_score, calc_mi_to_many_score, calc_xcorr_coelution_score, calc_xcorr_shape_score, calc_xcorr_shape_to_many_score, calc_xcorr_to_many_score}};


pub fn compute_alignment_scores(
    aligned_chromatograms: Vec<AlignedChromatogram>
) -> HashMap<String, FullTraceAlignmentScores> {
    let mut alignment_scores = HashMap::new();

    // Iterate over each aligned chromatogram
    for aligned_chrom in &aligned_chromatograms {
        let aligned_filename = aligned_chrom.chromatogram.metadata.get("basename")
            .unwrap_or(&aligned_chrom.chromatogram.native_id)
            .to_string();

        // For current aligned_chrom, get the corresponding reference chromatogram
        let reference_filename = aligned_chrom.rt_mapping[0].get("run1").unwrap();
        let reference_chrom = aligned_chromatograms.iter()
            .find(|chrom| chrom.chromatogram.metadata.get("basename").unwrap() == reference_filename)
            .unwrap();

        let xcorr_coelution_to_ref = calc_xcorr_coelution_score(&Array1::from(reference_chrom.chromatogram.intensities.clone()), &Array1::from(aligned_chrom.chromatogram.intensities.clone()));

        let xcorr_shape_to_ref = calc_xcorr_shape_score(&Array1::from(reference_chrom.chromatogram.intensities.clone()), &Array1::from(aligned_chrom.chromatogram.intensities.clone()));

        let mi_to_ref = calc_mi_score(&Array1::from(reference_chrom.chromatogram.intensities.clone()), &Array1::from(aligned_chrom.chromatogram.intensities.clone()));

        let xcorr_coelution_to_all = calc_xcorr_to_many_score(&Array1::from(aligned_chrom.chromatogram.intensities.clone()), &aligned_chromatograms.iter().map(|chrom| Array1::from(chrom.chromatogram.intensities.clone())).collect::<Vec<Array1<f64>>>());

        let xcorr_shape_to_all = calc_xcorr_shape_to_many_score(&Array1::from(aligned_chrom.chromatogram.intensities.clone()), &aligned_chromatograms.iter().map(|chrom| Array1::from(chrom.chromatogram.intensities.clone())).collect::<Vec<Array1<f64>>>());

        let mi_to_all = calc_mi_to_many_score(&Array1::from(aligned_chrom.chromatogram.intensities.clone()), &aligned_chromatograms.iter().map(|chrom| Array1::from(chrom.chromatogram.intensities.clone())).collect::<Vec<Array1<f64>>>());

        let alignment_score = FullTraceAlignmentScores {
            reference_filename: reference_filename.to_string(),
            aligned_filename: aligned_filename.to_string(),
            xcorr_coelution_to_ref,
            xcorr_shape_to_ref,
            mi_to_ref,
            xcorr_coelution_to_all,
            xcorr_shape_to_all,
            mi_to_all,
        };

        alignment_scores.insert(aligned_filename, alignment_score);        
    }
    alignment_scores
}


/// Computes scores for the PeakMapping struct based on aligned chromatograms and peak mappings.
/// 
/// # Parameters
/// - `aligned_chromatograms`: A list of aligned chromatograms.
/// - `peak_mappings`: A dictionary of peak mappings for each run.
/// 
/// # Returns
/// A dictionary of peak mappings with updated scores.
pub fn compute_peak_mapping_scores(
    aligned_chromatograms: Vec<AlignedChromatogram>,
    peak_mappings: HashMap<String, Vec<PeakMapping>>,
) -> HashMap<String, Vec<PeakMapping>> {
    let mut scored_peak_mappings = peak_mappings.clone();

    // Iterate over each aligned chromatogram
    for aligned_chrom in &aligned_chromatograms {
        let aligned_filename = aligned_chrom.chromatogram.metadata.get("basename")
            .unwrap_or(&aligned_chrom.chromatogram.native_id)
            .to_string();

        // For current aligned_chrom, get the corresponding reference chromatogram
        let reference_filename = aligned_chrom.rt_mapping[0].get("run1").unwrap();
        let reference_chrom = aligned_chromatograms.iter()
            .find(|chrom| chrom.chromatogram.metadata.get("basename").unwrap() == reference_filename)
            .unwrap();

        // Get the peak mappings for this aligned chromatogram
        if let Some(peak_mappings_for_run) = scored_peak_mappings.get_mut(&aligned_filename) {
            // Iterate over each peak mapping
            for peak_mapping in peak_mappings_for_run.iter_mut() {
                log::trace!("Scoring peak mapping for peak_id: {} for run: {}", peak_mapping.alignment_id, aligned_filename);

                let reference_intensities = get_peak_intensities(&reference_chrom.chromatogram, peak_mapping.reference_left_width, peak_mapping.reference_right_width);
                let aligned_intensities = get_peak_intensities(&aligned_chrom.chromatogram, peak_mapping.aligned_left_width, peak_mapping.aligned_right_width);

                // Compute cross-correlation to reference
                let xcorr_coelution_to_ref = calc_xcorr_coelution_score(&Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));
                peak_mapping.xcorr_coelution_to_ref = Some(xcorr_coelution_to_ref);

                let xcorr_shape_to_ref = calc_xcorr_shape_score(&Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));
                peak_mapping.xcorr_shape_to_ref = Some(xcorr_shape_to_ref);

                let mi_to_ref = calc_mi_score(&Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));
                peak_mapping.mi_to_ref = Some(mi_to_ref);

                // Subset peak_mappings for peaks with the same alignment_id across runs
                let peak_mappings_for_alignment_id: Vec<PeakMapping> = peak_mappings.clone().into_iter()
                .filter(|(_, mappings)| mappings.iter().any(|m| m.alignment_id == peak_mapping.alignment_id))
                .flat_map(|(_, mappings)| mappings)
                .collect();

                let xcorr_coelution_to_all = calc_xcorr_to_many_score(&Array1::from(aligned_intensities.clone()), &get_array_peak_intensities(aligned_chromatograms.clone(), peak_mappings_for_alignment_id.clone()));
                peak_mapping.xcorr_coelution_to_all = Some(xcorr_coelution_to_all);

                let xcorr_shape_to_all = calc_xcorr_shape_to_many_score(&Array1::from(aligned_intensities.clone()), &get_array_peak_intensities(aligned_chromatograms.clone(), peak_mappings_for_alignment_id.clone()));
                peak_mapping.xcorr_shape_to_all = Some(xcorr_shape_to_all);

                let mi_to_all = calc_mi_to_many_score(&Array1::from(aligned_intensities.clone()), &get_array_peak_intensities(aligned_chromatograms.clone(), peak_mappings_for_alignment_id.clone()));
                peak_mapping.mi_to_all = Some(mi_to_all);

                // Compute retention time deviation
                let rt_deviation = (peak_mapping.aligned_rt - peak_mapping.reference_rt).abs();
                peak_mapping.rt_deviation = Some(rt_deviation);

                // Compute peak intensity ratio
                let intensity_ratio = compute_peak_intensity_ratio(
                    &Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));

                peak_mapping.intensity_ratio = Some(intensity_ratio);
            }
        }
    }

    scored_peak_mappings
}

pub fn compute_peak_mapping_transitions_scores(
    aligned_identifying_trgrps: Vec<TransitionGroup>,
    aligned_chromatograms: Vec<AlignedChromatogram>,
    peak_mappings: HashMap<String, Vec<PeakMapping>>
) -> HashMap<String, Vec<AlignedTransitionScores>> {
    let mut scored_peak_mapped_transitions = HashMap::new();

    // Iterate over each transition group 
    for identifying_trgrp in &aligned_identifying_trgrps {
        let current_filename = identifying_trgrp.metadata.get("basename").unwrap()
            .to_string();

        let reference_filename = aligned_chromatograms.iter()
            .find(|chrom| chrom.chromatogram.metadata.get("basename").unwrap() == &current_filename)
            .unwrap().rt_mapping[0].get("run1").unwrap();
        let reference_chrom = aligned_chromatograms.iter()
            .find(|chrom| chrom.chromatogram.metadata.get("basename").unwrap() == reference_filename)
            .unwrap();

        // Create a temporary HashMap to store the results
        let temp_results: HashMap<String, Vec<AlignedTransitionScores>> = identifying_trgrp.chromatograms.par_iter().map(|(transition_id, transition_chrom)| {
            let mut identifying_peak_mapped_scores = Vec::new();
            if let Some(peak_mappings_for_run) = peak_mappings.get(&current_filename) {
                for peak_mapping in peak_mappings_for_run.iter() {
                    let reference_intensities = get_peak_intensities(&reference_chrom.chromatogram, peak_mapping.reference_left_width, peak_mapping.reference_right_width);
                    let aligned_intensities = get_peak_intensities(&transition_chrom, peak_mapping.aligned_left_width, peak_mapping.aligned_right_width);

                    // Compute cross-correlation to reference
                    let xcorr_coelution_to_ref = calc_xcorr_coelution_score(&Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));

                    let xcorr_shape_to_ref = calc_xcorr_shape_score(&Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));

                    let mi_to_ref = calc_mi_score(&Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));

                    // Subset peak_mappings for peaks with the same alignment_id across runs
                    let peak_mappings_for_alignment_id: Vec<PeakMapping> = peak_mappings.clone().into_iter()
                        .filter(|(_, mappings)| mappings.iter().any(|m| m.alignment_id == peak_mapping.alignment_id))
                        .flat_map(|(_, mappings)| mappings)
                        .collect();

                    // Compute cross-correlation to all (average across all aligned chromatograms)
                    let xcorr_coelution_to_all = calc_xcorr_to_many_score(&Array1::from(aligned_intensities.clone()), &get_array_peak_intensities(aligned_chromatograms.clone(), peak_mappings_for_alignment_id.clone()));

                    let xcorr_shape_to_all = calc_xcorr_shape_to_many_score(&Array1::from(aligned_intensities.clone()), &get_array_peak_intensities(aligned_chromatograms.clone(), peak_mappings_for_alignment_id.clone()));

                    let mi_to_all = calc_mi_to_many_score(&Array1::from(aligned_intensities.clone()), &get_array_peak_intensities(aligned_chromatograms.clone(), peak_mappings_for_alignment_id.clone()));

                    // Compute retention time deviation
                    let rt_deviation = (peak_mapping.aligned_rt - peak_mapping.reference_rt).abs();

                    // Compute peak intensity ratio
                    let intensity_ratio = compute_peak_intensity_ratio(
                        &Array1::from(reference_intensities.clone()), &Array1::from(aligned_intensities.clone()));

                    // Store the scores for the peak mapping
                    let aligned_identifying_transition_scores = AlignedTransitionScores {
                        feature_id: peak_mapping.alignment_id,
                        transition_id: transition_id.clone().parse::<i64>().unwrap(),
                        label: peak_mapping.label,
                        xcorr_coelution_to_ref: Some(xcorr_coelution_to_ref),
                        xcorr_shape_to_ref: Some(xcorr_shape_to_ref),
                        mi_to_ref: Some(mi_to_ref),
                        xcorr_coelution_to_all: Some(xcorr_coelution_to_all),
                        xcorr_shape_to_all: Some(xcorr_shape_to_all),
                        mi_to_all: Some(mi_to_all),
                        rt_deviation: Some(rt_deviation),
                        intensity_ratio: Some(intensity_ratio),
                    };
                    identifying_peak_mapped_scores.push(aligned_identifying_transition_scores);
                }
            }
            (current_filename.clone(), identifying_peak_mapped_scores)
        }).collect();

        // Merge the temporary results into the main HashMap
        for (key, value) in temp_results {
            scored_peak_mapped_transitions.entry(key).or_insert_with(Vec::new).extend(value);
        }
    }
    scored_peak_mapped_transitions
}


/// Computes the peak intensity ratio between two intensity arrays around a given retention time.
/// 
/// # Parameters
/// - `intensities1`: Intensity array of the first chromatogram.
/// - `intensities2`: Intensity array of the second chromatogram.
/// 
/// # Returns
/// The peak intensity ratio between the two intensity arrays.
fn compute_peak_intensity_ratio(
    intensities1: &Array1<f64>,
    intensities2: &Array1<f64>
) -> f64 {

    // Compute peak intensity ratio
    let intensity1: f64 = intensities1.iter().sum();
    let intensity2: f64 = intensities2.iter().sum();

    // Account for NaN values
    if intensity1 == 0.0 {
        0.0
    } else {
        intensity2 / intensity1
    }
}

/// Get peak intensities given the left and right width.
/// 
/// # Parameters
/// - `aligned_chrom`: An aligned chromatogram.
/// - `left_width`: The left width.
/// - `right_width`: The right width.
/// 
/// # Returns
/// A vector of peak intensities.
fn get_peak_intensities(
    chrom: &Chromatogram,
    left_width: f64,
    right_width: f64,
) -> Vec<f64> {
    let left_boundary_idx = find_closest_index(&chrom.retention_times, left_width).unwrap_or(0);
    let right_boundary_idx = find_closest_index(&chrom.retention_times, right_width).unwrap_or(chrom.retention_times.len() - 1);
    let peak_intensities = trim_vector(&chrom.intensities, left_boundary_idx, right_boundary_idx);
    peak_intensities
}


/// Get peak intensities for all aligned chromatograms given the aligned left and right width mappings
/// 
/// # Parameters
/// - `aligned_chromatograms`: A list of aligned chromatograms.
/// - `peak_mappings`: A dictionary of peak mappings for each run.
/// 
/// # Returns
/// A list of peak intensities for each aligned peak.
fn get_array_peak_intensities(
    aligned_chromatograms: Vec<AlignedChromatogram>,
    peak_mappings: Vec<PeakMapping>,
) -> Vec<Array1<f64>> {

    let mut peak_intensities_array = Vec::new();

    for chrom in aligned_chromatograms {
        let filename = chrom.chromatogram.metadata.get("basename").unwrap_or(&chrom.chromatogram.native_id).to_string();
        if let Some(peak_mapping_for_run) = peak_mappings.iter().find(|m| m.aligned_filename == filename) {
            let peak_intensities = get_peak_intensities(&chrom.chromatogram, peak_mapping_for_run.aligned_left_width, peak_mapping_for_run.aligned_right_width);
            peak_intensities_array.push(Array1::from(peak_intensities));
        }
    }

    peak_intensities_array
}

/// Trim vector given the start and end indices.
/// 
/// # Parameters
/// - `x`: A vector of values.
/// - `start_idx`: The start index.
/// - `end_idx`: The end index.
///
/// # Returns
/// A trimmed vector.
fn trim_vector(x: &Vec<f64>, start_idx: usize, end_idx: usize) -> Vec<f64> {
    // Ensure the indices are valid
    if start_idx >= x.len() || end_idx >= x.len() {
        panic!("Invalid indices, start_idx: {}, end_idx: {}, vector length: {}", start_idx, end_idx, x.len());
    }

    // Check if the start index is greater than the end index
    if start_idx > end_idx {
        panic!("Start index is greater than end index");
    }

    x[start_idx..end_idx].to_vec()
}

/// Finds the index of the closest value to the target value in an array.
/// 
/// # Parameters
/// - `x`: A vector of values.
/// - `target`: The target value.
/// 
/// # Returns
/// The index of the closest value to the target value.
fn find_closest_index(x: &Vec<f64>, target: f64) -> Option<usize> {
    if x.is_empty() {
        panic!("Empty array");
    }

    let mut min_diff = f64::MAX;
    let mut closest_index = 0;

    for (index, time) in x.iter().enumerate() {
        let diff = (time - target).abs();
        if diff < min_diff {
            min_diff = diff;
            closest_index = index;
        }
    }

    Some(closest_index)
}


/// Creates decoy peaks by shuffling the aligned peaks within each run.
/// 
/// # Parameters
/// - `peak_mappings`: A dictionary of peak mappings for each run.
/// 
/// # Returns
/// A dictionary of peak mappings with shuffled aligned peaks.
pub fn create_decoy_peaks_by_shuffling(
    peak_mappings: &HashMap<String, Vec<PeakMapping>>,
) -> HashMap<String, Vec<PeakMapping>> {
    let mut decoy_peak_mappings = peak_mappings.clone();
    let mut rng = rand::rng();

    // Iterate over each run (filename) and shuffle aligned peaks within the run
    for (_, peaks) in decoy_peak_mappings.iter_mut() {
        // Collect the aligned peak attributes to shuffle
        let mut aligned_rts: Vec<f64> = peaks.iter().map(|p| p.aligned_rt).collect();
        let mut aligned_left_widths: Vec<f64> = peaks.iter().map(|p| p.aligned_left_width).collect();
        let mut aligned_right_widths: Vec<f64> = peaks.iter().map(|p| p.aligned_right_width).collect();

        // Shuffle the aligned peak attributes
        aligned_rts.shuffle(&mut rng);
        aligned_left_widths.shuffle(&mut rng);
        aligned_right_widths.shuffle(&mut rng);

        // Assign the shuffled aligned peak attributes back to the peaks
        for (i, peak) in peaks.iter_mut().enumerate() {
            peak.aligned_rt = aligned_rts[i];

            // Ensure the left width is less than the right width
            if aligned_left_widths[i] > aligned_right_widths[i] {
                let temp = aligned_left_widths[i];
                aligned_left_widths[i] = aligned_right_widths[i];
                aligned_right_widths[i] = temp;
            }

            peak.aligned_left_width = aligned_left_widths[i];
            peak.aligned_right_width = aligned_right_widths[i];
            peak.aligned_feature_id = -1; // Mark as decoy
            peak.label = -1; // Mark as decoy
        }
    }

    decoy_peak_mappings
}

/// Creates decoy peaks by randomly selecting regions in the aligned chromatogram for each peak.
/// 
/// # Parameters
/// - `aligned_chromatograms`: A list of aligned chromatograms.
/// - `peak_mappings`: A dictionary of peak mappings for each run.
/// - `window_size`: Size of the decoy peak in retention time points.
/// 
/// # Returns
/// A dictionary of peak mappings with randomly selected regions for the aligned peaks.
pub fn create_decoy_peaks_by_random_regions(
    aligned_chromatograms: &[AlignedChromatogram],
    peak_mappings: &HashMap<String, Vec<PeakMapping>>,
    window_size: usize, 
) -> HashMap<String, Vec<PeakMapping>> {
    let mut decoy_peak_mappings = peak_mappings.clone();
    let mut rng = rand::rng();

    // Iterate over each run (filename) and create decoy peaks
    for (run_id, peaks) in decoy_peak_mappings.iter_mut() {
        // Find the corresponding chromatogram for this run
        let chromatogram = aligned_chromatograms
            .iter()
            .find(|chrom| chrom.chromatogram.metadata.get("basename").unwrap() == run_id)
            .expect("Chromatogram not found for run");

        let retention_times = &chromatogram.chromatogram.retention_times;

        // Ensure the window size is valid
        if window_size >= retention_times.len() {
            panic!("Window size is larger than the retention time array");
        }

        // Create decoy peaks by randomly selecting regions for the aligned peak
        for peak in peaks.iter_mut() {
            // Randomly select a start index for the decoy peak in the aligned chromatogram
            let start_idx = rng.random_range(0..retention_times.len() - window_size);
            let end_idx = start_idx + window_size;

            // Update only the aligned peak information
            peak.aligned_rt = retention_times[start_idx];
            peak.aligned_left_width = retention_times[start_idx];
            peak.aligned_right_width = retention_times[end_idx];
            peak.aligned_feature_id = -1; // Mark as decoy
            peak.label = -1; // Mark as decoy
        }
    }

    decoy_peak_mappings
}

