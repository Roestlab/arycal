use rayon::prelude::*;
use ndarray::{Array1, Array2};
use std::{collections::HashMap, sync::Arc};
use rand::{seq::SliceRandom, Rng};

use arycal_cloudpath::sqmass::TransitionGroup;
use arycal_common::{chromatogram::{Chromatogram, AlignedChromatogram}, AlignedTransitionScores, FullTraceAlignmentScores, PeakMapping};
use crate::{alignment::alignment::validate_widths, stats::{calc_mi_score, calc_mi_to_many_score, calc_xcorr_coelution_score, calc_xcorr_shape_score, calc_xcorr_shape_to_many_score, calc_xcorr_to_many_score}};


pub fn compute_alignment_scores(
    aligned_chromatograms: Vec<AlignedChromatogram>
) -> HashMap<String, FullTraceAlignmentScores> {
    // Pre-convert all intensities to Array1 once
    let chrom_arrays: Vec<Array1<f64>> = aligned_chromatograms.iter()
        .map(|chrom| Array1::from(chrom.chromatogram.intensities.clone()))
        .collect();

    // Create lookup maps
    let chrom_lookup: HashMap<_, _> = aligned_chromatograms.iter()
        .map(|chrom| (
            chrom.chromatogram.metadata.get("basename")
                .unwrap_or(&chrom.chromatogram.native_id)
                .as_str(),
            chrom
        ))
        .collect();

    // Process in parallel
    aligned_chromatograms.par_iter()
        .map(|aligned_chrom| {
            let aligned_filename = aligned_chrom.chromatogram.metadata.get("basename")
                .unwrap_or(&aligned_chrom.chromatogram.native_id)
                .to_string();

            let reference_filename = match aligned_chrom.rt_mapping[0].get("run1") {
                Some(name) => name,
                None => return (aligned_filename, None),
            };

            let reference_chrom = match chrom_lookup.get(reference_filename.as_str()) {
                Some(c) => c,
                None => return (aligned_filename, None),
            };

            let aligned_array = Array1::from(aligned_chrom.chromatogram.intensities.clone());
            let reference_array = Array1::from(reference_chrom.chromatogram.intensities.clone());

            // Compute all scores
            let xcorr_coelution_to_ref = calc_xcorr_coelution_score(
                &reference_array,
                &aligned_array
            );

            let xcorr_shape_to_ref = calc_xcorr_shape_score(
                &reference_array,
                &aligned_array
            );

            let mi_to_ref = calc_mi_score(
                &reference_array,
                &aligned_array
            );

            // Convert chrom_arrays to slice for "to_all" functions
            let all_arrays_slice = &chrom_arrays;

            let xcorr_coelution_to_all = calc_xcorr_to_many_score(
                &aligned_array,
                all_arrays_slice
            );

            let xcorr_shape_to_all = calc_xcorr_shape_to_many_score(
                &aligned_array,
                all_arrays_slice
            );

            let mi_to_all = calc_mi_to_many_score(
                &aligned_array,
                all_arrays_slice
            );

            let alignment_score = FullTraceAlignmentScores {
                reference_filename: reference_filename.to_string(),
                aligned_filename: aligned_filename.clone(),
                xcorr_coelution_to_ref,
                xcorr_shape_to_ref,
                mi_to_ref,
                xcorr_coelution_to_all,
                xcorr_shape_to_all,
                mi_to_all,
            };

            (aligned_filename, Some(alignment_score))
        })
        .filter_map(|(k, v)| v.map(|score| (k, score)))
        .collect()
}


pub fn compute_peak_mapping_scores(
    aligned_chromatograms: Vec<AlignedChromatogram>,
    peak_mappings: HashMap<String, Vec<PeakMapping>>,
) -> HashMap<String, Vec<PeakMapping>> {
    // Create Arc wrappers for shared access
    let chroms_arc = Arc::new(aligned_chromatograms);
    let peak_maps_arc = Arc::new(peak_mappings);

    // Create lookup structure (using Arc reference)
    let chrom_lookup: HashMap<_, _> = chroms_arc.iter()
        .map(|c| (
            c.chromatogram.metadata.get("basename")
                .unwrap_or(&c.chromatogram.native_id)
                .as_str(),
            c
        ))
        .collect();

    // Process in parallel using Arc references
    peak_maps_arc.par_iter()
        .map(|(aligned_filename, peak_mappings_for_run)| {
            // Clone the Vec to mutate it
            let mut peak_mappings_for_run = peak_mappings_for_run.clone();

            let aligned_chrom = match chrom_lookup.get(aligned_filename.as_str()) {
                Some(c) => c,
                None => return (aligned_filename.clone(), peak_mappings_for_run),
            };

            let reference_filename = match aligned_chrom.rt_mapping[0].get("run1") {
                Some(name) => name,
                None => return (aligned_filename.clone(), peak_mappings_for_run),
            };

            let reference_chrom = match chrom_lookup.get(reference_filename.as_str()) {
                Some(c) => c,
                None => return (aligned_filename.clone(), peak_mappings_for_run),
            };

            for peak_mapping in &mut peak_mappings_for_run {
                // Get intensities once and reuse
                let reference_intensities = Array1::from(get_peak_intensities(
                    &reference_chrom.chromatogram,
                    peak_mapping.reference_left_width,
                    peak_mapping.reference_right_width
                ));
                
                let aligned_intensities = Array1::from(get_peak_intensities(
                    &aligned_chrom.chromatogram,
                    peak_mapping.aligned_left_width,
                    peak_mapping.aligned_right_width
                ));

                // Compute reference scores
                peak_mapping.xcorr_coelution_to_ref = Some(calc_xcorr_coelution_score(
                    &reference_intensities, 
                    &aligned_intensities
                ));
                
                peak_mapping.xcorr_shape_to_ref = Some(calc_xcorr_shape_score(
                    &reference_intensities, 
                    &aligned_intensities
                ));
                
                peak_mapping.mi_to_ref = Some(calc_mi_score(
                    &reference_intensities, 
                    &aligned_intensities
                ));

                // Compute "to_all" scores using Arc references
                if let Some(all_intensities) = get_all_intensities_for_alignment(
                    &chroms_arc,
                    &peak_maps_arc,
                    peak_mapping.alignment_id.try_into().unwrap()
                ) {
                    peak_mapping.xcorr_coelution_to_all = Some(calc_xcorr_to_many_score(
                        &aligned_intensities, 
                        &all_intensities
                    ));
                    
                    peak_mapping.xcorr_shape_to_all = Some(calc_xcorr_shape_to_many_score(
                        &aligned_intensities, 
                        &all_intensities
                    ));
                    
                    peak_mapping.mi_to_all = Some(calc_mi_to_many_score(
                        &aligned_intensities, 
                        &all_intensities
                    ));
                }

                peak_mapping.rt_deviation = Some((peak_mapping.aligned_rt - peak_mapping.reference_rt).abs());
                peak_mapping.intensity_ratio = Some(compute_peak_intensity_ratio(
                    &reference_intensities, 
                    &aligned_intensities
                ));
            }

            (aligned_filename.clone(), peak_mappings_for_run)
        })
        .collect()
}

pub fn compute_peak_mapping_transitions_scores(
    aligned_identifying_trgrps: Vec<TransitionGroup>,
    aligned_chromatograms: Vec<AlignedChromatogram>,
    peak_mappings: HashMap<String, Vec<PeakMapping>>,
) -> HashMap<String, Vec<AlignedTransitionScores>> {
    // Wrap data in Arc for shared access
    let chroms_arc = Arc::new(aligned_chromatograms);
    let peak_maps_arc = Arc::new(peak_mappings);
    
    // Create lookup maps
    let chrom_lookup: HashMap<_, _> = chroms_arc.iter()
        .map(|c| (
            c.chromatogram.metadata.get("basename")
                .unwrap_or(&c.chromatogram.native_id)
                .as_str(),
            c
        ))
        .collect();

    // Process transition groups in parallel
    aligned_identifying_trgrps.par_iter()
        .flat_map(|identifying_trgrp| {
            let current_filename = identifying_trgrp.metadata.get("basename").unwrap();
            
            // Get reference chromatogram once per group
            let (reference_chrom, peak_mappings_for_run) = {
                let current_chrom = match chrom_lookup.get(current_filename.as_str()) {
                    Some(c) => c,
                    None => return Vec::new(),
                };
                
                let reference_filename = match current_chrom.rt_mapping[0].get("run1") {
                    Some(name) => name,
                    None => return Vec::new(),
                };
                
                let reference_chrom = match chrom_lookup.get(reference_filename.as_str()) {
                    Some(c) => c,
                    None => return Vec::new(),
                };
                
                let peak_mappings_for_run = match peak_maps_arc.get(current_filename) {
                    Some(mappings) => mappings,
                    None => return Vec::new(),
                };
                
                (reference_chrom, peak_mappings_for_run)
            };

            // Process transitions in parallel
            identifying_trgrp.chromatograms.par_iter()
                .flat_map(|(transition_id, transition_chrom)| {
                    peak_mappings_for_run.par_iter()
                        .map(|peak_mapping| {
                            // Get intensities once
                            let reference_intensities = Array1::from(get_peak_intensities(
                                &reference_chrom.chromatogram,
                                peak_mapping.reference_left_width,
                                peak_mapping.reference_right_width
                            ));
                            
                            let aligned_intensities = Array1::from(get_peak_intensities(
                                transition_chrom,
                                peak_mapping.aligned_left_width,
                                peak_mapping.aligned_right_width
                            ));

                            // Compute reference scores
                            let xcorr_coelution_to_ref = calc_xcorr_coelution_score(
                                &reference_intensities, 
                                &aligned_intensities
                            );
                            
                            let xcorr_shape_to_ref = calc_xcorr_shape_score(
                                &reference_intensities, 
                                &aligned_intensities
                            );
                            
                            let mi_to_ref = calc_mi_score(
                                &reference_intensities, 
                                &aligned_intensities
                            );

                            // Compute "to_all" scores
                            let all_intensities = get_all_intensities_for_alignment(
                                &chroms_arc,
                                &peak_maps_arc,
                                peak_mapping.alignment_id.clone()
                            );
                            
                            let xcorr_coelution_to_all = all_intensities.as_ref()
                                .map(|ints| calc_xcorr_to_many_score(&aligned_intensities, ints))
                                .unwrap_or(0.0);
                            
                            let xcorr_shape_to_all = all_intensities.as_ref()
                                .map(|ints| calc_xcorr_shape_to_many_score(&aligned_intensities, ints))
                                .unwrap_or(0.0);
                            
                            let mi_to_all = all_intensities.as_ref()
                                .map(|ints| calc_mi_to_many_score(&aligned_intensities, ints))
                                .unwrap_or(0.0);

                            // Create scores struct
                            AlignedTransitionScores {
                                feature_id: peak_mapping.aligned_feature_id,
                                transition_id: transition_id.parse().unwrap_or(0),
                                run_id: peak_mapping.run_id,
                                aligned_filename: current_filename.clone(),
                                label: peak_mapping.label,
                                xcorr_coelution_to_ref: Some(xcorr_coelution_to_ref),
                                xcorr_shape_to_ref: Some(xcorr_shape_to_ref),
                                mi_to_ref: Some(mi_to_ref),
                                xcorr_coelution_to_all: Some(xcorr_coelution_to_all),
                                xcorr_shape_to_all: Some(xcorr_shape_to_all),
                                mi_to_all: Some(mi_to_all),
                                rt_deviation: Some((peak_mapping.aligned_rt - peak_mapping.reference_rt).abs()),
                                intensity_ratio: Some(compute_peak_intensity_ratio(
                                    &reference_intensities, 
                                    &aligned_intensities
                                )),
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .fold(
            || HashMap::new(),
            |mut acc, score| {
                let filename = score.aligned_filename.clone(); // Adjust based on your actual data
                acc.entry(filename).or_insert_with(Vec::new).push(score);
                acc
            }
        )
        .reduce(
            || HashMap::new(),
            |mut a, b| {
                for (k, v) in b {
                    a.entry(k).or_insert_with(Vec::new).extend(v);
                }
                a
            }
        )
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


// /// Get peak intensities for all aligned chromatograms given the aligned left and right width mappings
// /// 
// /// # Parameters
// /// - `aligned_chromatograms`: A list of aligned chromatograms.
// /// - `peak_mappings`: A dictionary of peak mappings for each run.
// /// 
// /// # Returns
// /// A list of peak intensities for each aligned peak.
// fn get_array_peak_intensities(
//     aligned_chromatograms: Vec<AlignedChromatogram>,
//     peak_mappings: Vec<PeakMapping>,
// ) -> Vec<Array1<f64>> {

//     let mut peak_intensities_array = Vec::new();

//     for chrom in aligned_chromatograms {
//         let filename = chrom.chromatogram.metadata.get("basename").unwrap_or(&chrom.chromatogram.native_id).to_string();
//         if let Some(peak_mapping_for_run) = peak_mappings.iter().find(|m| m.aligned_filename == filename) {
//             let peak_intensities = get_peak_intensities(&chrom.chromatogram, peak_mapping_for_run.aligned_left_width, peak_mapping_for_run.aligned_right_width);
//             peak_intensities_array.push(Array1::from(peak_intensities));
//         }
//     }

//     peak_intensities_array
// }

/// Optimized version of get_array_peak_intensities
fn get_array_peak_intensities(
    aligned_chromatograms: &[AlignedChromatogram],
    peak_mappings: &[&PeakMapping],  // Changed to accept references
) -> Vec<Array1<f64>> {
    // Create lookup map for faster filename searching
    let peak_map_lookup: HashMap<_, _> = peak_mappings.iter()
        .map(|m| (&m.aligned_filename, *m))  // Dereference here
        .collect();

    aligned_chromatograms.iter()
        .filter_map(|chrom| {
            let filename = chrom.chromatogram.metadata.get("basename")
                .unwrap_or(&chrom.chromatogram.native_id);
            
            peak_map_lookup.get(filename).map(|mapping| {
                Array1::from(get_peak_intensities(
                    &chrom.chromatogram,
                    mapping.aligned_left_width,
                    mapping.aligned_right_width
                ))
            })
        })
        .collect()
}

/// Helper function to get all intensities for an alignment ID
fn get_all_intensities_for_alignment(
    chromatograms: &Arc<Vec<AlignedChromatogram>>,
    peak_mappings: &Arc<HashMap<String, Vec<PeakMapping>>>,
    alignment_id: i64,
) -> Option<Vec<Array1<f64>>> {
    let relevant_mappings: Vec<_> = peak_mappings.values()
        .flat_map(|mappings| mappings.iter())
        .filter(|m| m.alignment_id == alignment_id)
        .collect();

    if relevant_mappings.is_empty() {
        return None;
    }

    Some(get_array_peak_intensities(
        chromatograms.as_ref(),
        &relevant_mappings
    ))
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
    // // Ensure the indices are valid
    // if start_idx >= x.len() || end_idx >= x.len() {
    //     panic!("Invalid indices, start_idx: {}, end_idx: {}, vector length: {}", start_idx, end_idx, x.len());
    // }

    // // Check if the start index is greater than the end index
    // if start_idx > end_idx {
    //     panic!("Start index is greater than end index");
    // }

    let (valid_start_idx, valid_end_index) = validate_widths(start_idx as f64, end_idx as f64);
    let start_idx = valid_start_idx as usize;
    let end_idx = valid_end_index as usize;

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
        let mut aligned_feature_ids: Vec<i64> = peaks.iter().map(|p| p.aligned_feature_id).collect();

        // Shuffle the aligned peak attributes
        aligned_rts.shuffle(&mut rng);
        aligned_left_widths.shuffle(&mut rng);
        aligned_right_widths.shuffle(&mut rng);
        aligned_feature_ids.shuffle(&mut rng);

        // Assign the shuffled aligned peak attributes back to the peaks
        for (i, peak) in peaks.iter_mut().enumerate() {
            peak.aligned_rt = aligned_rts[i];

            // Ensure the left width is less than the right width
            if aligned_left_widths[i] > aligned_right_widths[i] {
                let temp = aligned_left_widths[i];
                aligned_left_widths[i] = aligned_right_widths[i];
                aligned_right_widths[i] = temp;
            }
            // log::debug!("Decoy feature_id: {}", aligned_feature_ids[i]);
            peak.aligned_left_width = aligned_left_widths[i];
            peak.aligned_right_width = aligned_right_widths[i];
            peak.aligned_feature_id = aligned_feature_ids[i];
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
            peak.aligned_rt = (retention_times[start_idx] + retention_times[end_idx]) / 2.0;
            peak.aligned_left_width = retention_times[start_idx];
            peak.aligned_right_width = retention_times[end_idx];
            peak.aligned_feature_id = -1; // Mark as decoy
            peak.label = -1; // Mark as decoy
        }
    }

    decoy_peak_mappings
}


/// Creates a feature matrix and labels for all the peak mappings.
pub fn create_feature_matrix(
    peak_mappings: &HashMap<String, Vec<PeakMapping>>
) -> (Array2<f64>, Array1<i32>) {

    // feature matrix should be of shape len(peak_mappings) + inner len of each peak_mapping by 7 features
    let nrows: usize = peak_mappings.iter().map(|(_, mappings)| mappings.len()).sum::<usize>() + peak_mappings.len();
    let ncols = 8;
    let mut feature_matrix = Array2::zeros((nrows, ncols));
    let mut labels = Array1::zeros(nrows);

    let mut row_idx = 0;
    for (_filename, mappings) in peak_mappings.iter() {
        for mapping in mappings {
            feature_matrix.row_mut(row_idx).assign(&Array1::from(vec![
                mapping.xcorr_coelution_to_ref.unwrap_or(0.0),
                mapping.xcorr_shape_to_ref.unwrap_or(0.0),
                mapping.mi_to_ref.unwrap_or(0.0),
                mapping.xcorr_coelution_to_all.unwrap_or(0.0),
                mapping.xcorr_shape_to_all.unwrap_or(0.0),
                mapping.mi_to_all.unwrap_or(0.0),
                mapping.rt_deviation.unwrap_or(-1.0), // TODO: Should this be -1.0 or some large value?
                mapping.intensity_ratio.unwrap_or(0.0),
            ]));
            labels[row_idx] = mapping.label;
            row_idx += 1;
        }
    }

    (feature_matrix, labels)
}