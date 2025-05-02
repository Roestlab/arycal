use std::collections::{HashMap, HashSet};
use arycal_common::chromatogram::AlignedRTPointPair;
use rand::prelude::IndexedRandom;
use anyhow::Error as AnyHowError;
use dtw_rs::{Algorithm, DynamicTimeWarping};
use rayon::prelude::*;

use arycal_cloudpath::sqmass::TransitionGroup;

use super::alignment::calculate_distance;
use super::alignment::construct_mst;
use arycal_common::chromatogram::{Chromatogram, AlignedChromatogram, pad_chromatograms};
use arycal_common::config::AlignmentConfig;
use arycal_cloudpath::util::extract_basename;


/// Aligns two chromatograms based on the optimal path.
/// 
/// # Parameters
/// - `chrom1`: The first chromatogram
/// - `chrom2`: The second chromatogram
/// - `optimal_path`: The optimal path calculated by `dtw.path`
/// 
/// # Returns
/// - A tuple of two aligned chromatograms
pub fn align_chromatograms(
    chrom1: &Chromatogram,
    chrom2: &Chromatogram,
    optimal_path: &[(usize, usize)],
) -> (Chromatogram, Chromatogram) {
    let mut aligned_rt = Vec::with_capacity(optimal_path.len());
    let mut aligned_data1 = Vec::with_capacity(optimal_path.len());
    let mut aligned_data2 = Vec::with_capacity(optimal_path.len());

    for &(i, j) in optimal_path {
        match (i.checked_sub(1), j.checked_sub(1)) {
            (Some(i), Some(j)) => {
                aligned_rt.push((chrom1.retention_times[i] + chrom2.retention_times[j]) / 2.0);
                aligned_data1.push(chrom1.intensities[i]);
                aligned_data2.push(chrom2.intensities[j]);
            }
            (Some(i), None) => {
                aligned_rt.push(chrom1.retention_times[i]);
                aligned_data1.push(chrom1.intensities[i]);
                aligned_data2.push(0.0);
            }
            (None, Some(j)) => {
                aligned_rt.push(chrom2.retention_times[j]);
                aligned_data1.push(0.0);
                aligned_data2.push(chrom2.intensities[j]);
            }
            _ => continue,
        }
    }

    // Reuse metadata with additional alignment info
    let mut meta1 = chrom1.metadata.clone();
    meta1.insert("is_aligned".to_string(), "true".to_string());
    meta1.insert("reference_run".to_string(), 
        chrom1.metadata.get("basename").unwrap_or(&chrom1.native_id).clone());

    let mut meta2 = chrom2.metadata.clone();
    meta2.insert("is_aligned".to_string(), "true".to_string());
    meta2.insert("reference_run".to_string(),
        chrom1.metadata.get("basename").unwrap_or(&chrom1.native_id).clone());

    (
        Chromatogram {
            retention_times: aligned_rt.clone(),
            intensities: aligned_data1,
            metadata: meta1,
            ..chrom1.clone()
        },
        Chromatogram {
            retention_times: aligned_rt,
            intensities: aligned_data2,
            metadata: meta2,
            ..chrom2.clone()
        }
    )
}


/// Creates a mapping between the original retention times (RT) of two chromatograms based on the optimal path.
// fn create_rt_mapping(
//     optimal_path: &[(usize, usize)],
//     chrom1: &Chromatogram,
//     chrom2: &Chromatogram,
// ) -> Vec<HashMap<String, f64>> {
//     // let run1_name = chrom1.metadata.get("basename").unwrap_or(&chrom1.native_id);
//     // let run2_name = chrom2.metadata.get("basename").unwrap_or(&chrom2.native_id);

//     optimal_path.iter()
//         .filter(|&&(i, j)| i > 0 && j > 0)
//         .map(|&(i, j)| {
//             let rt1 = chrom1.retention_times[i - 1];
//             let rt2 = chrom2.retention_times[j - 1];
            
//             let mut entry = HashMap::with_capacity(0);
//             entry.insert("rt1".to_string(), rt1);
//             entry.insert("rt2".to_string(), rt2);
//             // entry.insert("alignment".to_string(), format!("({}, {})", rt1, rt2));
//             // entry.insert("run1".to_string(), run1_name.clone());
//             // entry.insert("run2".to_string(), run2_name.clone());
//             entry
//         })
//         .collect()
// }

fn create_rt_mapping(
    optimal_path: &[(usize, usize)],
    chrom1: &Chromatogram,
    chrom2: &Chromatogram,
) -> Vec<AlignedRTPointPair> {
    optimal_path.iter()
        .filter(|&&(i, j)| i > 0 && j > 0)
        .map(|&(i, j)| AlignedRTPointPair {
            rt1: chrom1.retention_times[i - 1] as f32,
            rt2: chrom2.retention_times[j - 1] as f32,
        })
        .collect()
}


/// Aligns a series of chromatograms by picking a random run as the reference and aligning all others to it.
/// 
/// # Parameters
/// - `smoothed_tics`: A vector of smoothed chromatograms
/// - `params`: Extra alignment configuration parameters
/// 
/// # Returns
/// - A vector of aligned chromatograms with their alignment paths
pub fn star_align_tics(
    smoothed_tics: &Vec<Chromatogram>,
    params: &AlignmentConfig,
) -> Result<Vec<AlignedChromatogram>, AnyHowError> {
    // Early return if insufficient chromatograms
    if smoothed_tics.len() < 2 {
        // return Err(AnyHowError::msg("At least two chromatograms are required for alignment"));
        log::warn!("At least two chromatograms are required for alignment - returning empty result");
        return Ok(Vec::new());
    }

    // Select reference chromatogram
    let reference_chrom = if let Some(ref_chrom) = &params.reference_run {
        smoothed_tics.iter()
            .find(|x| x.metadata.get("basename").unwrap_or(&x.native_id) == &extract_basename(ref_chrom))
            .ok_or_else(|| AnyHowError::msg("Reference chromatogram not found"))?
    } else {
        let mut rng = rand::rng();
        let binding = (0..smoothed_tics.len()).collect::<Vec<_>>();
        let reference_idx = binding.choose(&mut rng).unwrap();
        &smoothed_tics[*reference_idx]
    };

    let ref_basename = reference_chrom.metadata.get("basename").unwrap_or(&reference_chrom.native_id);
    // log::trace!("Selected reference run: {}", ref_basename);

    // Pre-compute reference intensities once
    let ref_intensities = &reference_chrom.intensities;

    // Process chromatograms in parallel
    let aligned_chromatograms: Vec<_> = smoothed_tics.par_iter()
        .enumerate()
        .filter_map(|(_idx, chrom)| {
            let chrom_basename = chrom.metadata.get("basename").unwrap_or(&chrom.native_id);
            // if chrom_basename == ref_basename {
            //     return None; // Skip reference chromatogram
            // }

            // log::trace!("Aligning run: {} to reference", chrom_basename);

            // Compute DTW alignment
            let dtw = DynamicTimeWarping::between(ref_intensities, &chrom.intensities);
            let path = dtw.path();

            // Align chromatograms
            let (_, aligned_chrom) = align_chromatograms(reference_chrom, chrom, &path);
            let rt_mapping = create_rt_mapping(&path, reference_chrom, chrom);

            // Check if we want to retain the alignment path
            let path_out = if params.retain_alignment_path {
                path.to_vec()
            } else {
                Vec::new()
            };

            Some(AlignedChromatogram {
                chromatogram: aligned_chrom,
                alignment_path: path_out,
                lag: None,
                rt_mapping,
                reference_basename: ref_basename.clone(),
            })
        })
        .collect();

    Ok(aligned_chromatograms)
}


/// Aligns a series of chromatograms progressively.
/// 
/// # Parameters
/// - `smoothed_tics`: A vector of smoothed chromatograms
/// 
/// # Returns
/// - A vector of aligned chromatograms with their alignment paths
pub fn progressive_align_tics(smoothed_tics: &Vec<Chromatogram>) -> Result<Vec<AlignedChromatogram>, AnyHowError> {
    // Ensure we have at least two tics to align
    if smoothed_tics.len() < 2 {
        return Err(AnyHowError::msg("At least two chromatograms are required for alignment"));
    }

    // Initialize a vector to hold results
    let mut results: Vec<AlignedChromatogram> = Vec::new();

    // Start with the first tic
    let mut aligned_sum = smoothed_tics[0].clone();

    // println!("Using run: {:?} as the first reference", smoothed_tics[0].metadata.get("basename"));

    let mut trgrp_aligned = TransitionGroup::new("Aligned".to_string());
    let mut tmp_chrom = aligned_sum.clone();
    tmp_chrom.native_id = format!("First ref: {}", tmp_chrom.metadata.get("basename").unwrap_or(&tmp_chrom.native_id).clone());
    trgrp_aligned.add_chromatogram(tmp_chrom);
    

    // Iterate over each subsequent tic
    for i in 0..smoothed_tics.len() {
        let current_tic = &smoothed_tics[i];
        // println!("Aligning run: {:?}", current_tic.metadata.get("basename"));

        // Create a common RT space
        let common_aligned_rt_space = pad_chromatograms(vec![aligned_sum.clone(), current_tic.clone()]);
        aligned_sum = common_aligned_rt_space[0].clone();
        let current_tic = &common_aligned_rt_space[1].clone();

        let (_ref_rt, ref_intensities) = (
            aligned_sum.retention_times.clone(),
            aligned_sum.intensities.clone(),
        );
        let (_query_rt, query_intensities) = (
            current_tic.retention_times.clone(),
            current_tic.intensities.clone(),
        );

        // Compute the DTW alignment
        let dtw = DynamicTimeWarping::between(&ref_intensities, &query_intensities);

        // Align the two chromatograms
        let (mut aligned_sum_temp, aligned_current_tic) = align_chromatograms(&aligned_sum, current_tic, &dtw.path());

        // Create_rt_mapping
        let rt_mapping = create_rt_mapping(&dtw.path(), &aligned_sum, &current_tic);

        // Add alignment path for the first TIC used as reference
        if i == 1 {
            results[0].alignment_path = dtw.path().clone(); // Store the alignment path of the first TIC
            results[0].rt_mapping = rt_mapping.clone(); // Store the RT mapping of the first TIC
        }

        // Store the aligned current TIC and its alignment path
        results.push(AlignedChromatogram {
            chromatogram: aligned_current_tic.clone(),
            alignment_path: dtw.path(),
            lag: None,
            rt_mapping: rt_mapping,
            reference_basename: aligned_sum.metadata.get("basename").unwrap().to_string(),
        });

        // average aligned chromatograms
        for (j, intensity) in aligned_current_tic.intensities.iter().enumerate() {
            aligned_sum_temp.intensities[j] = (aligned_sum_temp.intensities[j] + intensity) / 2.0;
        }

        // Update aligned_sum for the next iteration
        aligned_sum = aligned_sum_temp.clone();

        
    }

    Ok(results)
}


/// Aligns a series of chromatograms using a Minimum Spanning Tree (MST) approach.
/// 
/// # Parameters
/// - `smoothed_tics`: A vector of smoothed chromatograms
/// 
/// # Returns
/// - A vector of aligned chromatograms with their alignment paths
pub fn mst_align_tics(smoothed_tics: &Vec<Chromatogram>) -> Result<Vec<AlignedChromatogram>, AnyHowError> {
    // Ensure we have at least two chromatograms to align
    if smoothed_tics.len() < 2 {
        return Err(AnyHowError::msg("At least two chromatograms are required for alignment"));
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

    // HashSet to track chromatograms that have been added to aligned_chromatograms
    let mut aligned_chromatogram_indices = HashSet::new();

    // Align chromatograms based on the edges in the MST
    for (chrom1_idx, chrom2_idx, _) in mst_edges {
        let chrom1 = &smoothed_tics[chrom1_idx].clone();
        let chrom2 = &smoothed_tics[chrom2_idx].clone();

        // println!("Aligning runs: {:?} and {:?}", chrom1.metadata.get("basename"), chrom2.metadata.get("basename"));

        let (_ref_rt, ref_intensities) = (
            chrom1.retention_times.clone(),
            chrom1.intensities.clone(),
        );
        let (_query_rt, query_intensities) = (
            chrom2.retention_times.clone(),
            chrom2.intensities.clone(),
        );

        // Compute the DTW alignment
        let dtw = DynamicTimeWarping::between(&ref_intensities, &query_intensities);

        // Align the chromatograms
        let (aligned_chrom1, aligned_chrom2) = align_chromatograms(chrom1, chrom2, &dtw.path());
        let rt_mapping = create_rt_mapping(&dtw.path(), chrom1, chrom2);

        // Store aligned chromatograms of 1st chromatogram only if not already added
        if !aligned_chromatogram_indices.contains(&chrom1_idx) {
            aligned_chromatograms.push(AlignedChromatogram {
                chromatogram: aligned_chrom1.clone(),
                alignment_path: dtw.path().clone(),
                lag: None,
                rt_mapping: rt_mapping.clone(),
                reference_basename: chrom1.metadata.get("basename").unwrap().to_string(),
            });
            aligned_chromatogram_indices.insert(chrom1_idx);
        }

        // Store aligned chromatograms of 2nd chromatogram only if not already added
        if !aligned_chromatogram_indices.contains(&chrom2_idx) {
            aligned_chromatograms.push(AlignedChromatogram {
                chromatogram: aligned_chrom2.clone(),
                alignment_path: dtw.path().clone(),
                lag: None,
                rt_mapping: rt_mapping.clone(),
                reference_basename: chrom1.metadata.get("basename").unwrap().to_string(),
            });
            aligned_chromatogram_indices.insert(chrom2_idx);
        }
    }

    Ok(aligned_chromatograms)
}