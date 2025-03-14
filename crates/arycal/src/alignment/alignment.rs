use std::cmp::Ordering;
use std::collections::HashMap;
use std::f64;
use arycal_common::config::{AlignmentConfig, SmoothingConfig};
use rand::seq::IndexedRandom;
use union_find::{QuickFindUf, UnionBySize, UnionFind};

use arycal_cloudpath::osw::{FeatureData, ValueEntryType};
use arycal_cloudpath::sqmass::{apply_common_rt_space_single, Chromatogram, TransitionGroup};
use arycal_common::PeakMapping;

use super::fast_fourier_lag::shift_chromatogram;

#[derive(Debug)]
pub enum ValueType {
    Int32(i32),
    Int64(i64),
    Float(f64),
    Text(String),
}

impl ValueType {
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            ValueType::Int32(v) => Some(*v as i64),
            ValueType::Int64(v) => Some(*v),
            ValueType::Float(v) => Some(*v as i64),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            ValueType::Int32(v) => Some(*v as f64),
            ValueType::Int64(v) => Some(*v as f64),
            ValueType::Float(v) => Some(*v),
            _ => None,
        }
    }

    pub fn as_string(&self) -> Option<String> {
        match self {
            ValueType::Text(s) => Some(s.clone()),
            _ => None,
        }
    }
}

impl From<ValueType> for i64 {
    fn from(value: ValueType) -> Self {
        match value {
            ValueType::Int32(v) => v as i64,
            ValueType::Int64(v) => v,
            ValueType::Float(v) => v as i64,
            ValueType::Text(_) => 0,
        }
    }
}

impl From<ValueType> for f64 {
    fn from(value: ValueType) -> Self {
        match value {
            ValueType::Int32(v) => v as f64,
            ValueType::Int64(v) => v as f64,
            ValueType::Float(v) => v,
            ValueType::Text(_) => 0.0,
        }
    }
}

impl From<ValueType> for String {
    fn from(value: ValueType) -> Self {
        match value {
            ValueType::Text(s) => s,
            _ => String::new(),
        }
    }
}


/// Enum for the alignment method. Either DTW, FFT or FFT-DTW.
/// FFT-DTW is a hybrid method that uses FFT for cross-correlation and DTW for local refinement.
#[derive(Debug, Clone)]
pub enum AlignmentMethod {
    DTW,
    FFT,
    FFTDTW,
}

impl Default for AlignmentMethod {
    fn default() -> Self {
        AlignmentMethod::FFT
    }
    
}

impl AlignmentMethod {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentMethod::DTW => "dtw",
            AlignmentMethod::FFT => "ffw",
            AlignmentMethod::FFTDTW => "fftdtw",
        }
    }
    
}

/// Enum for reference method. Either STAR, MST, or PROGRESSIVE.
/// - STAR uses a single rrandomly selected reference chromatogram.
/// - MST constructs a minimum spanning tree from pairwise distances and selects the centroid.
/// - PROGRESSIVE aligns chromatograms in a progressive manner. Selects the first chromatogram as the reference, and then aligns the next chromatogram to the reference. The aligned chromatogram is averaged with the reference, to be used as the new reference for the next chromatogram. This process is repeated until all chromatograms are aligned.
#[derive(Debug, Clone)]
pub enum ReferenceMethod {
    STAR,
    MST,
    PROGRESSIVE,
}

impl Default for ReferenceMethod {
    fn default() -> Self {
        ReferenceMethod::STAR
    }
    
}

impl ReferenceMethod {
    pub fn as_str(&self) -> &str {
        match self {
            ReferenceMethod::STAR => "star",
            ReferenceMethod::MST => "mst",
            ReferenceMethod::PROGRESSIVE => "progressive",
        }
    }
    
}

/// Represents the mapping of peaks across chromatograms.
#[derive(Debug, Clone)]
pub struct AlignedChromatogram {
    /// Aligned chromatogram
    pub chromatogram: Chromatogram,
    /// Optimal alignment path between the reference and query chromatograms (Only for DTW and FFT-DTW)
    pub alignment_path: Vec<(usize, usize)>,
    /// Lag between the reference and query chromatograms (Only for FFT and FFT-DTW)
    pub lag: Option<isize>,
    /// Mapping of retention times between the original and aligned chromatograms
    pub rt_mapping: Vec<HashMap<String, String>>,
}

/// Calculates the Euclidean distance between two chromatograms.
///
/// # Parameters
/// - `chrom1`: First chromatogram
/// - `chrom2`: Second chromatogram
///
/// # Returns
/// - The Euclidean distance between the two chromatograms
pub fn calculate_distance(chrom1: &Chromatogram, chrom2: &Chromatogram) -> f64 {
    // Ensure the chromatograms have the same length
    if chrom1.intensities.len() != chrom2.intensities.len() {
        panic!("Chromatograms must have the same number of intensities");
    }

    // Calculate the Euclidean distance between the intensities
    let mut sum = 0.0;
    for (i1, i2) in chrom1.intensities.iter().zip(&chrom2.intensities) {
        sum += (i1 - i2).powi(2);
    }

    sum.sqrt() // Return the Euclidean distance (square root of sum of squares)
}

/// Constructs a minimum spanning tree (MST) from a list of distances between chromatograms.
///
/// # Parameters
/// - `distances`: A list of distances between chromatograms (chrom1, chrom2, distance)
/// - `num_chromatograms`: The total number of chromatograms
///
/// # Returns
/// - A list of edges in the MST (chrom1, chrom2, distance)
pub fn construct_mst(
    distances: &[(usize, usize, f64)],
    num_chromatograms: usize,
) -> Vec<(usize, usize, f64)> {
    // Use QuickFindUf with UnionBySize as the union strategy
    let mut uf: QuickFindUf<UnionBySize> = QuickFindUf::new(num_chromatograms);
    let mut mst_edges = Vec::new();

    // Sort edges by distance
    let mut sorted_distances = distances.to_vec();
    sorted_distances.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

    // Kruskal's algorithm
    for (i, j, dist) in sorted_distances {
        if uf.find(i) != uf.find(j) {
            uf.union(i, j);
            mst_edges.push((i, j, dist)); // Now return the distance as well
        }
    }

    mst_edges
}

/// Validates the widths of two peaks and switches the order if necessary.
pub fn validate_widths(left_width: f64, right_width: f64) -> (f64, f64) {
    if left_width >= right_width {
        log::trace!(
            "Invalid widths: left_width ({}) is not smaller than right_width ({}). Switching the order.",
            left_width,
            right_width
        );
        (right_width, left_width)
    } else {
        (left_width, right_width)
    }
}

/// Maps peaks across aligned chromatograms using the alignment information.
///
/// # Parameters
///
/// - `aligned_chrom`: The aligned chromatogram.
/// - `reference_features`: The features in the reference chromatogram.
/// - `aligned_features`: The features in the query chromatogram that are to be mapped to the aligned chromatogram.
/// - `rt_tolerance`: The retention time tolerance for mapping peaks.
pub fn map_peaks_across_runs(
    aligned_chrom: &AlignedChromatogram,
    reference_features: Vec<FeatureData>,
    aligned_features: Vec<FeatureData>,
    rt_tolerance: f64,
    alignment_config: &AlignmentConfig,
) -> Vec<PeakMapping> {
    let mut peak_mappings = Vec::new();

    log::trace!("There are {} reference features and {} query features to map", reference_features[0].feature_id.clone().unwrap().as_multiple().unwrap().len(), aligned_features[0].feature_id.clone().unwrap().as_multiple().unwrap().len() );

    // Step 1: Map peaks from reference to query chromatogram
    for ref_feature in &reference_features { // TODO: This is only a Vec of one element which contains inner Vecs of features
        let mut alignment_id = 0;

        for (i, &rt) in ref_feature.exp_rt.as_multiple().unwrap().iter().enumerate() {
            // Get the aligned query rt that maps to the reference rt
            let target_rt = map_retention_time(rt, &aligned_chrom.rt_mapping);

            // Reverse the rt_mapping to get the aligned query rt in it's original space
            // let original_target_rt = reverse_rt_mapping(target_rt, &aligned_chrom, alignment_config).unwrap();


            log::trace!("Mapping closest aligned query RT: {:?} to reference RT feature: {:?}", target_rt, rt);

            // Generate alignment_id based on reference_rt (or another unique identifier)
            // let alignment_id = rt.to_bits() as i64; // Use the bits of the reference_rt as alignment_id

            if let Some((aligned_feature_id, aligned_rt, aligned_left_width, aligned_right_width)) =
                find_closest_feature(target_rt, &aligned_features, rt_tolerance)
            {
                log::trace!("Found query feature (id: {}) mapping to reference feature (id: {}): {} -> {}", aligned_feature_id, ref_feature.feature_id.clone().unwrap().as_multiple().unwrap()[i], aligned_rt, rt);

                // TODO: Really shouldn't need to have to validate widths, as these are derived from OpenSwath 
                let (validated_left_width_ref, validated_right_width_ref) = validate_widths(
                    ref_feature.left_width.as_ref().unwrap().as_multiple().unwrap()[i],
                    ref_feature.right_width.as_ref().unwrap().as_multiple().unwrap()[i],
                );

                let (validated_left_width_aligned, validated_right_width_aligned) = validate_widths(
                    aligned_left_width,
                    aligned_right_width,
                );
                
                peak_mappings.push(PeakMapping {
                    alignment_id, 
                    precursor_id: ref_feature.precursor_id.clone(),
                    run_id: aligned_features[0].run_id.clone(),
                    reference_feature_id: ref_feature.feature_id.clone().unwrap().as_multiple().unwrap()[i],
                    aligned_feature_id,
                    reference_rt: rt,
                    aligned_rt,
                    reference_left_width: validated_left_width_ref,
                    reference_right_width: validated_right_width_ref,
                    aligned_left_width: validated_left_width_aligned,
                    aligned_right_width: validated_right_width_aligned,
                    reference_filename: ref_feature.basename.clone(),
                    aligned_filename: aligned_features[0].basename.clone(),
                    label: 1,
                    xcorr_coelution_to_ref: None,
                    xcorr_shape_to_ref: None,
                    mi_to_ref: None,
                    xcorr_coelution_to_all: None,
                    xcorr_shape_to_all: None,
                    mi_to_all: None,
                    rt_deviation: None,
                    intensity_ratio: None,
                });
            } else {
                log::trace!("Couldn't find a matching feature for reference RT: {:?} with id: {}", rt, ref_feature.feature_id.clone().unwrap().as_multiple().unwrap()[i]);
                // // Recover missing peak in the query chromatogram
                // log::trace!("Recovering missing peak in query chromatogram for reference RT: {:?}", rt);
                // let (validated_left_width_ref, validated_right_width_ref) = validate_widths(
                //     ref_feature.left_width.as_ref().unwrap().as_multiple().unwrap()[i],
                //     ref_feature.right_width.as_ref().unwrap().as_multiple().unwrap()[i],
                // );
                // peak_mappings.push(PeakMapping {
                //     alignment_id, // Use the same alignment_id for the same peak across runs
                //     reference_feature_id: *ref_feature.feature_id.clone().unwrap(),
                //     aligned_feature_id: -1, // Use -1 to indicate a missing peak
                //     reference_rt: rt,
                //     aligned_rt: target_rt,
                //     reference_left_width: validated_left_width_ref,
                //     reference_right_width: validated_right_width_ref,
                //     aligned_left_width: validated_left_width_ref, // Use reference values as placeholder
                //     aligned_right_width: validated_right_width_ref, // Use reference values as placeholder
                //     reference_filename: ref_feature.basename.clone(),
                //     aligned_filename: aligned_features[0].basename.clone(),
                //     label: 1,
                //     xcorr_coelution_to_ref: None,
                //     xcorr_shape_to_ref: None,
                //     mi_to_ref: None,
                //     xcorr_coelution_to_all: None,
                //     xcorr_shape_to_all: None,
                //     mi_to_all: None,
                //     rt_deviation: None,
                //     intensity_ratio: None,
                // });
            }
            alignment_id += 1;
        }
    }

    // // Step 2: Map peaks from query to reference chromatogram (to recover missing peaks in the reference)
    // for aligned_feature in &aligned_features {
    //     for (i, &rt) in aligned_feature.exp_rt.as_multiple().unwrap().iter().enumerate() {
    //         let reference_rt = map_retention_time(rt, &aligned_chrom.rt_mapping);

    //         // println!("Mapping reference RT: {:?} to query aligned RT: {:?}", reference_rt, rt);

    //         // Generate alignment_id based on reference_rt (or another unique identifier)
    //         let alignment_id = reference_rt.to_bits() as i64; // Use the bits of the reference_rt as alignment_id

    //         if let Some((reference_feature_id, reference_rt, reference_left_width, reference_right_width)) =
    //             find_closest_feature(reference_rt, &reference_features, rt_tolerance)
    //         {
    //             // Skip if the peak is already mapped in Step 1
    //             if !peak_mappings.iter().any(|m| m.aligned_feature_id == *aligned_feature.feature_id.clone().unwrap()) {
    //                 // println!("Found reference feature mapping to query feature: {:?}", reference_feature_id);
    //                 peak_mappings.push(PeakMapping {
    //                     alignment_id, // Use the same alignment_id for the same peak across runs
    //                     reference_feature_id,
    //                     aligned_feature_id: *aligned_feature.feature_id.clone().unwrap(),
    //                     reference_rt,
    //                     aligned_rt: rt,
    //                     reference_left_width,
    //                     reference_right_width,
    //                     aligned_left_width: aligned_feature.left_width.as_ref().unwrap().as_multiple().unwrap()[i],
    //                     aligned_right_width: aligned_feature.right_width.as_ref().unwrap().as_multiple().unwrap()[i],
    //                     reference_filename: reference_features[0].basename.clone(),
    //                     aligned_filename: aligned_feature.basename.clone(),
    //                 });
    //             }
    //         } else {
    //             // Recover missing peak in the reference chromatogram
    //             // println!("Recovering missing peak in reference chromatogram for query RT: {:?}", rt);
    //             peak_mappings.push(PeakMapping {
    //                 alignment_id, // Use the same alignment_id for the same peak across runs
    //                 reference_feature_id: -1, // Use -1 to indicate a missing peak
    //                 aligned_feature_id: *aligned_feature.feature_id.clone().unwrap(),
    //                 reference_rt,
    //                 aligned_rt: rt,
    //                 reference_left_width: aligned_feature.left_width.as_ref().unwrap().as_multiple().unwrap()[i], // Use query values as placeholder
    //                 reference_right_width: aligned_feature.right_width.as_ref().unwrap().as_multiple().unwrap()[i], // Use query values as placeholder
    //                 aligned_left_width: aligned_feature.left_width.as_ref().unwrap().as_multiple().unwrap()[i],
    //                 aligned_right_width: aligned_feature.right_width.as_ref().unwrap().as_multiple().unwrap()[i],
    //                 reference_filename: reference_features[0].basename.clone(),
    //                 aligned_filename: aligned_feature.basename.clone(),
    //             });
    //         }
    //     }
    // }

    // Step 3: Remove overlapping peaks
    // let filtered_peaks = remove_overlapping_peaks(peak_mappings);

    peak_mappings
}

/// Removes overlapping peaks within the same run by comparing peak boundaries.
fn remove_overlapping_peaks(peak_mappings: Vec<PeakMapping>) -> Vec<PeakMapping> {
    let mut filtered_peaks = Vec::new();

    // Group peaks by run (filename or chromatogram ID)
    let mut peaks_grouped: HashMap<String, Vec<PeakMapping>> = HashMap::new();
    for peak in peak_mappings {
        let run_id = peak.aligned_filename.clone(); // Assuming `aligned_filename` is part of PeakMapping
        peaks_grouped
            .entry(run_id)
            .or_insert_with(Vec::new)
            .push(peak);
    }

    // Process each group of peaks separately
    for (_, peaks) in peaks_grouped {
        let mut non_overlapping_peaks = Vec::new();

        // Sort peaks by retention time (reference_rt or aligned_rt)
        let mut sorted_peaks = peaks.clone();
        sorted_peaks.sort_by(|a, b| a.reference_rt.partial_cmp(&b.reference_rt).unwrap());

        // Iterate through sorted peaks and remove overlaps
        let mut prev_peak: Option<PeakMapping> = None;
        for peak in sorted_peaks {
            if let Some(prev) = &prev_peak {
                log::trace!("Checking for overlapping peaks: Peak 1 (RT: {}), Peak 2 (RT: {})", prev.reference_rt, peak.reference_rt);
                // Check if the current peak overlaps with the previous peak
                let current_left = peak.reference_left_width;
                let current_right = peak.reference_right_width;
                let prev_left = prev.reference_left_width;
                let prev_right = prev.reference_right_width;

                if current_left <= prev_right && current_right >= prev_left {
                    // Overlapping peaks detected
                    // println!(
                    //     "Overlapping peaks detected: Peak 1 (RT: {}-{}), Peak 2 (RT: {}-{})",
                    //     prev_left, prev_right, current_left, current_right
                    // );

                    // Resolve overlapping peaks based on feature IDs
                    if peak.aligned_feature_id != -1 || peak.reference_feature_id != -1 {
                        // Prefer the peak with a valid feature ID
                        if prev.aligned_feature_id == -1 && prev.reference_feature_id == -1 {
                            log::trace!("Removing overlapping peaks with missing feature IDs");
                            // Replace the previous peak with the current one
                            non_overlapping_peaks.pop();
                            non_overlapping_peaks.push(peak.clone());
                        }
                    } else {
                        log::trace!("Both peaks have missing feature IDs: computing a consensus peak");
                        // Both peaks have missing IDs: compute a consensus peak
                        let consensus_peak = PeakMapping {
                            alignment_id: peak.alignment_id,
                            precursor_id: peak.precursor_id,
                            run_id: peak.run_id,
                            reference_feature_id: -1,
                            aligned_feature_id: -1,
                            reference_rt: (prev.reference_rt + peak.reference_rt) / 2.0,
                            aligned_rt: (prev.aligned_rt + peak.aligned_rt) / 2.0,
                            reference_left_width: (prev.reference_left_width
                                + peak.reference_left_width)
                                / 2.0,
                            reference_right_width: (prev.reference_right_width
                                + peak.reference_right_width)
                                / 2.0,
                            aligned_left_width: (prev.aligned_left_width + peak.aligned_left_width)
                                / 2.0,
                            aligned_right_width: (prev.aligned_right_width
                                + peak.aligned_right_width)
                                / 2.0,
                            reference_filename: peak.reference_filename.clone(),
                            aligned_filename: peak.aligned_filename.clone(),
                            label: 1,
                            xcorr_coelution_to_ref: None,
                            xcorr_shape_to_ref: None,
                            mi_to_ref: None,
                            xcorr_coelution_to_all: None,
                            xcorr_shape_to_all: None,
                            mi_to_all: None,
                            rt_deviation: None,
                            intensity_ratio: None,
                        };
                        non_overlapping_peaks.pop();
                        non_overlapping_peaks.push(consensus_peak);
                    }
                } else {
                    // No overlap: add the current peak
                    non_overlapping_peaks.push(peak.clone());
                }
            } else {
                // First peak: add it
                non_overlapping_peaks.push(peak.clone());
            }

            // Update the previous peak
            prev_peak = Some(peak);
        }

        // Add non-overlapping peaks to the final list
        filtered_peaks.extend(non_overlapping_peaks);
    }

    filtered_peaks
}

/// Maps a retention time from the reference chromatogram to the aligned chromatogram.
fn map_retention_time(rt: f64, rt_mapping: &[HashMap<String, String>]) -> f64 {
    // If the mapping is empty, return the original RT (no mapping)
    // TODO: Should we panic here instead?
    if rt_mapping.is_empty() {
        return rt;
    }

    // Find the nearest points in the mapping for interpolation
    let mut lower_idx = 0;
    let mut upper_idx = rt_mapping.len() - 1;

    // Binary search to find the nearest lower and upper points
    for (idx, map) in rt_mapping.iter().enumerate() {
        let rt1 = map
            .get("rt1")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        if rt1 <= rt {
            lower_idx = idx;
        } else {
            upper_idx = idx;
            break;
        }
    }

    // Get the lower and upper RT pairs
    let lower_map = &rt_mapping[lower_idx];
    let upper_map = &rt_mapping[upper_idx];

    let rt1_lower = lower_map
        .get("rt1")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);
    let rt2_lower = lower_map
        .get("rt2")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);

    let rt1_upper = upper_map
        .get("rt1")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);
    let rt2_upper = upper_map
        .get("rt2")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);

    // If the RT matches exactly, return the mapped RT
    if rt == rt1_lower {
        return rt2_lower;
    }
    if rt == rt1_upper {
        return rt2_upper;
    }

    // Interpolate between the nearest points
    let slope = (rt2_upper - rt2_lower) / (rt1_upper - rt1_lower);
    let mapped_rt = rt2_lower + slope * (rt - rt1_lower);

    mapped_rt
}

/// Finds the closest retention time in the aligned feature's exp_rt vector.
fn find_closest_feature(
    target_rt: f64,
    features: &[FeatureData],
    tolerance: f64,
) -> Option<(i64, f64, f64, f64)> {
    // Returns (feature_id, rt, left_width, right_width)
    let mut closest_match = None;
    let mut min_diff = f64::MAX;

    // Iterate over all features
    for feature in features {
        // Iterate over each retention time in the feature's exp_rt
        for (i, &rt) in feature.exp_rt.as_multiple().unwrap().iter().enumerate() {
            let diff = (rt - target_rt).abs();

            // Check if the retention time is within the tolerance window
            if diff <= tolerance && diff < min_diff {
                min_diff = diff;
                closest_match = Some((
                    feature.feature_id.clone().unwrap().as_multiple().unwrap()[i], // feature_id
                    rt,                                   // retention time
                    feature.left_width.as_ref().unwrap().as_multiple().unwrap()[i], // left_width
                    feature.right_width.as_ref().unwrap().as_multiple().unwrap()[i], // right_width
                ));
                log::trace!("Found closest potential matching feature to {} target rt: {:?}", target_rt, closest_match);
            }
        }
    }
    log::trace!("Returning  final closest match: {:?}", closest_match);
    closest_match
}


/// Applies alignment to a list of transition groups using an existing alignment result.
/// 
/// # Parameters
/// - `transition_groups`: The transition groups to be aligned
/// - `aligned_chromatograms`: The aligned chromatograms
/// - `common_rt_space`: The common retention time space
/// - `alignment_config`: The alignment configuration
pub fn apply_post_alignment_to_trgrp(
    transition_groups: Vec<TransitionGroup>,
    aligned_chromatograms: Vec<AlignedChromatogram>,
    common_rt_space: Vec<f64>,
    alignment_config: &AlignmentConfig
) -> Vec<TransitionGroup> {

    let mut aligned_transition_groups = Vec::new();

    for trgrp in transition_groups {
        let mut aligned_trgrp = trgrp.clone();

        let current_filename = trgrp.metadata.get("basename").unwrap();

        let mut aligned_identifying_chromatograms= HashMap::new();
        // Apply post-alignment to each transition in the group
        for (transition_id, transition_xic) in &mut aligned_trgrp.chromatograms {

            let common_rt_chrom = apply_common_rt_space_single(transition_xic.clone(), &common_rt_space);

            let smooth_chrom = common_rt_chrom.smooth_sgolay(alignment_config.smoothing.sgolay_window, alignment_config.smoothing.sgolay_order).unwrap().normalize().unwrap();

            let query_aligned_chrom = aligned_chromatograms.iter().find(|chrom| chrom.chromatogram.metadata.get("basename").unwrap() == current_filename).unwrap();

            let aligned_chrom = apply_post_alignment_to_chromatogram(smooth_chrom, query_aligned_chrom.clone(), &alignment_config);

            aligned_identifying_chromatograms.insert(transition_id.clone(), aligned_chrom.chromatogram);
        }

        aligned_trgrp.chromatograms = aligned_identifying_chromatograms;
        aligned_transition_groups.push(aligned_trgrp);
    }

    aligned_transition_groups
}


/// Applies alignment to a chromatogram using an existing alignment result
/// 
/// # Parameters
/// - `chromatogram`: The chromatogram to be aligned
/// - `aligned_chromatogram`: The aligned chromatogram
/// - `alignment_config`: The alignment configuration
/// 
/// # Returns
/// - The aligned chromatogram
pub fn apply_post_alignment_to_chromatogram(
    chromatogram: Chromatogram,
    aligned_chromatogram: AlignedChromatogram,
    alignment_config: &AlignmentConfig,
) -> AlignedChromatogram {
    let aligned_chromatogram = match alignment_config.method.to_lowercase().as_str() {
        "dtw" => {
            let alignment_path = aligned_chromatogram.alignment_path.clone();

            let mut aligned_chrom = AlignedChromatogram { chromatogram: chromatogram.clone(), alignment_path: alignment_path.clone(), lag: None, rt_mapping: aligned_chromatogram.rt_mapping.clone() };

            // Apply the DTW alignment to the query chromatogram
            let (query_rt, query_intensities) = (
                chromatogram.retention_times.clone(),
                chromatogram.intensities.clone(),
            );

            let refined_rt: Vec<f64> = alignment_path.iter().map(|&(_, j)| query_rt[j]).collect();

            let refined_intensities: Vec<f64> = alignment_path
                .iter()
                .map(|&(_, j)| query_intensities[j])
                .collect();
            
            aligned_chrom.chromatogram.retention_times = refined_rt;
            aligned_chrom.chromatogram.intensities = refined_intensities;

            aligned_chrom
        },
        "fft" => {
            let lag = aligned_chromatogram.lag.unwrap();
            let aligned_chrom = AlignedChromatogram { chromatogram: shift_chromatogram(&chromatogram.clone(), lag), alignment_path: aligned_chromatogram.alignment_path.clone(), lag: Some(lag), rt_mapping: aligned_chromatogram.rt_mapping.clone() };
            
            aligned_chrom
        },
        "fftdtw" => {
            let alignment_path = aligned_chromatogram.alignment_path.clone();
            let lag = aligned_chromatogram.lag.unwrap();
            let mut aligned_chrom = shift_chromatogram(&chromatogram, lag);
            
            let (query_rt, query_intensities) = (
                aligned_chrom.retention_times.clone(),
                aligned_chrom.intensities.clone(),
            );

            // Apply the DTW alignment to the query chromatogram
            let refined_rt: Vec<f64> = alignment_path.iter().map(|&(_, j)| query_rt[j]).collect();

            let refined_intensities: Vec<f64> = alignment_path
                .iter()
                .map(|&(_, j)| query_intensities[j])
                .collect();

            aligned_chrom.retention_times = refined_rt;
            aligned_chrom.intensities = refined_intensities;

            AlignedChromatogram { chromatogram: aligned_chrom, alignment_path: alignment_path, lag: Some(lag), rt_mapping: aligned_chromatogram.rt_mapping.clone() }
        },
        _ => {
            let alignment_path = aligned_chromatogram.alignment_path.clone();

            let mut aligned_chrom = AlignedChromatogram { chromatogram: chromatogram.clone(), alignment_path: alignment_path.clone(), lag: None, rt_mapping: aligned_chromatogram.rt_mapping.clone() };

            // Apply the DTW alignment to the query chromatogram
            let (query_rt, query_intensities) = (
                chromatogram.retention_times.clone(),
                chromatogram.intensities.clone(),
            );

            let refined_rt: Vec<f64> = alignment_path.iter().map(|&(_, j)| query_rt[j]).collect();

            let refined_intensities: Vec<f64> = alignment_path
                .iter()
                .map(|&(_, j)| query_intensities[j])
                .collect();
            
            aligned_chrom.chromatogram.retention_times = refined_rt;
            aligned_chrom.chromatogram.intensities = refined_intensities;

            aligned_chrom
        },
    };
    aligned_chromatogram
}

/// Reverses the RT mapping to convert an aligned RT back to the original RT space.
///
/// # Parameters
/// - `aligned_rt`: The RT value in the aligned space.
/// - `aligned_chromatogram`: The aligned chromatogram containing the alignment path, lag, and RT mapping.
/// - `alignment_config`: The alignment configuration to determine the alignment method and parameters.
///
/// # Returns
/// - The original RT value corresponding to the aligned RT.
pub fn reverse_rt_mapping(
    aligned_rt: f64,
    aligned_chromatogram: &AlignedChromatogram,
    alignment_config: &AlignmentConfig,
) -> Option<f64> {
    
    match alignment_config.method.to_lowercase().as_str() {
        "dtw" => {
            log::debug!("Getting original RT for aligned RT: {} using DTW alignment", aligned_rt); 

            // Use rt_mapping to get index where target_rt is closest to 'rt1'
            let ref_rts = aligned_chromatogram.rt_mapping.iter().map(|m| m.get("rt1").unwrap().parse::<f64>().unwrap()).collect::<Vec<f64>>();

            let closest_index = find_closest_index(&ref_rts, aligned_rt)?;

            // Map back to the original RT using the alignment path
            let query_rts = aligned_chromatogram.rt_mapping.iter().map(|m| m.get("rt2").unwrap().parse::<f64>().unwrap()).collect::<Vec<f64>>();

            Some(query_rts[closest_index])
        }
        "fft" => {
            log::debug!("Getting original RT for aligned RT: {} using FFT alignment", aligned_rt);
            // For FFT, use the lag to reverse the mapping
            let lag = aligned_chromatogram.lag? as f64;
            Some(aligned_rt + lag)
        }
        "fftdtw" => {
            log::debug!("Getting original RT for aligned RT: {} using FFT-DTW alignment", aligned_rt);
            // For FFT-DTW, first reverse the FFT shift, then reverse the DTW mapping
            let lag = aligned_chromatogram.lag? as f64;
            let shifted_rt = aligned_rt + lag;

            // Use the alignment path to reverse the DTW mapping
            let ref_rts = aligned_chromatogram.rt_mapping.iter().map(|m| m.get("rt1").unwrap().parse::<f64>().unwrap()).collect::<Vec<f64>>();

            let closest_index = find_closest_index(&ref_rts, shifted_rt)?;

            // Map back to the original RT using the alignment path
            let query_rts = aligned_chromatogram.rt_mapping.iter().map(|m| m.get("rt2").unwrap().parse::<f64>().unwrap()).collect::<Vec<f64>>();

            Some(query_rts[closest_index])
        }
        _ => {
            // Default to DTW behavior if the method is unknown
            log::debug!("Getting original RT for aligned RT: {} using default DTW alignment", aligned_rt);

            // Use rt_mapping to get index where target_rt is closest to 'rt1'
            let ref_rts = aligned_chromatogram.rt_mapping.iter().map(|m| m.get("rt1").unwrap().parse::<f64>().unwrap()).collect::<Vec<f64>>();

            let closest_index = find_closest_index(&ref_rts, aligned_rt)?;

            // Map back to the original RT using the alignment path
            let query_rts = aligned_chromatogram.rt_mapping.iter().map(|m| m.get("rt2").unwrap().parse::<f64>().unwrap()).collect::<Vec<f64>>();

            Some(query_rts[closest_index])
        }
    }
}

/// Helper function to find the index of the closest value in a vector.
fn find_closest_index(values: &[f64], target: f64) -> Option<usize> {
    values
        .iter()
        .enumerate()
        .min_by(|(_, &a), (_, &b)| {
            let diff_a = (a - target).abs();
            let diff_b = (b - target).abs();
            diff_a.partial_cmp(&diff_b).unwrap_or(Ordering::Equal)
        })
        .map(|(index, _)| index)
}