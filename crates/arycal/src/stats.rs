use std::sync::atomic::{AtomicPtr, Ordering};
use ndarray::{Array1, Array2, s};
use rayon::prelude::*;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};


pub fn initialize_xcorr_matrix(intensities: &[Array1<f64>]) -> (Array2<f64>, Array2<f64>, Array2<f64>) {
    let n = intensities.len();
    let mut xcorr_matrix = Array2::<f64>::zeros((n, n));
    let mut xcorr_matrix_max_peak = Array2::<f64>::zeros((n, n));
    let mut xcorr_matrix_max_peak_sec = Array2::<f64>::zeros((n, n));

    // Precompute stats
    let stats: Vec<_> = intensities.par_iter()
        .map(|arr| {
            let mean = arr.mean().unwrap_or(0.0);
            let std = arr.std(0.0);
            (mean, std)
        })
        .collect();

    // Convert to atomic pointers
    let xcorr_ptr = AtomicPtr::new(xcorr_matrix.as_mut_ptr());
    let max_peak_ptr = AtomicPtr::new(xcorr_matrix_max_peak.as_mut_ptr());
    let max_peak_sec_ptr = AtomicPtr::new(xcorr_matrix_max_peak_sec.as_mut_ptr());

    (0..n).into_par_iter().for_each(|i| {
        unsafe {
            let xcorr_ptr = xcorr_ptr.load(Ordering::Relaxed);
            let max_peak_ptr = max_peak_ptr.load(Ordering::Relaxed);
            let max_peak_sec_ptr = max_peak_sec_ptr.load(Ordering::Relaxed);

            for j in i..n {
                let xcorr = normalized_cross_correlation_optimized(
                    &intensities[i],
                    &intensities[j],
                    stats[i].0, stats[i].1,
                    stats[j].0, stats[j].1
                );
                
                let (max_peak, max_peak_sec) = xcorr_array_get_max_peak_optimized(&xcorr);
                
                let idx_ij = i * n + j;
                let idx_ji = j * n + i;
                
                *xcorr_ptr.add(idx_ij) = max_peak;
                *max_peak_ptr.add(idx_ij) = max_peak;
                *max_peak_sec_ptr.add(idx_ij) = max_peak_sec;
                
                if i != j {
                    *xcorr_ptr.add(idx_ji) = max_peak;
                    *max_peak_ptr.add(idx_ji) = max_peak;
                    *max_peak_sec_ptr.add(idx_ji) = max_peak_sec;
                }
            }
        }
    });
    
    (xcorr_matrix, xcorr_matrix_max_peak, xcorr_matrix_max_peak_sec)
}


/// Get the max peak from a cross-correlation array.
pub fn xcorr_array_get_max_peak(xcorr: &Array1<f64>) -> (f64, f64) {
    let mut max_index = 0;
    let mut max_value = f64::NEG_INFINITY;

    for (i, &value) in xcorr.iter().enumerate() {
        if value.is_nan() {
            continue; // Skip NaN values
        }
        if value > max_value {
            max_value = value;
            max_index = i;
        }
    }

    (max_index as f64, max_value)
}


/// Optimized max peak finder
pub fn xcorr_array_get_max_peak_optimized(xcorr: &Array1<f64>) -> (f64, f64) {
    let (max_i, max_v) = xcorr.iter()
        .enumerate()
        .filter(|(_, &v)| !v.is_nan())
        .fold((0, f64::NEG_INFINITY), |(max_i, max_v), (i, &v)| {
            if v > max_v { (i, v) } else { (max_i, max_v) }
        });
    
    (max_i as f64, max_v)
}

/// Optimized normalized cross-correlation
pub fn normalized_cross_correlation_optimized(
    intensities1: &Array1<f64>,
    intensities2: &Array1<f64>,
    mean1: f64,
    std1: f64,
    mean2: f64,
    std2: f64,
) -> Array1<f64> {
    let min_len = intensities1.len().min(intensities2.len());
    let intensities1 = intensities1.slice(s![..min_len]);
    let intensities2 = intensities2.slice(s![..min_len]);
    
    // Pre-normalized slices
    let norm1 = (&intensities1 - mean1) / std1;
    let norm2 = (&intensities2 - mean2) / std2;

    let len = intensities1.len();
    let mut xcorr = Array1::<f64>::zeros(len);

    for i in 0..len {
        let valid_len = len - i;
        if valid_len == 0 {
            break;
        }
        
        let slice1 = norm1.slice(s![..valid_len]);
        let slice2 = norm2.slice(s![i..]);

        // Manual dot product optimized
        let dot_product = slice1.iter()
            .zip(slice2.iter())
            .fold(0.0, |acc, (&a, &b)| acc + a * b);
        
        xcorr[i] = dot_product / (valid_len as f64);  // Normalize by length
    }
    xcorr
}

// /// Compute the normalized cross-correlation between two intensity arrays.
// pub fn normalized_cross_correlation(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> Array1<f64> {
//     let mean1 = intensities1.mean().unwrap_or(0.0);
//     let mean2 = intensities2.mean().unwrap_or(0.0);
//     let std1 = intensities1.std(0.0);
//     let std2 = intensities2.std(0.0);
    
//     let norm1 = (intensities1 - mean1) / std1;
//     let norm2 = (intensities2 - mean2) / std2;
    
//     let mut xcorr = Array1::<f64>::zeros(intensities1.len());
//     for i in 0..intensities1.len() {
//         xcorr[i] = norm1.slice(s![..(intensities1.len() - i)]).dot(&norm2.slice(s![i..]));
//     }
//     xcorr
// }

/// Compute the normalized cross-correlation between two intensity arrays.
pub fn normalized_cross_correlation(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> Array1<f64> {
    let min_len = intensities1.len().min(intensities2.len());
    let intensities1 = &intensities1.slice(s![..min_len]);
    let intensities2 = &intensities2.slice(s![..min_len]);
    let mean1 = intensities1.mean().unwrap_or(0.0);
    let mean2 = intensities2.mean().unwrap_or(0.0);
    let std1 = intensities1.std(0.0);
    let std2 = intensities2.std(0.0);

    let norm1 = (intensities1 - mean1) / std1;
    let norm2 = (intensities2 - mean2) / std2;

    let len = intensities1.len();
    let mut xcorr = Array1::<f64>::zeros(len);

    for i in 0..len {
        let valid_len = len - i;
        if valid_len == 0 {
            break;
        }
        
        let slice1 = &norm1.slice(s![..valid_len]);
        let slice2 = &norm2.slice(s![i..(i + valid_len)]);

        let dot_product: f64 = slice1.iter().zip(slice2.iter()).map(|(a, b)| a * b).sum();
        
        xcorr[i] = dot_product;
    }
    xcorr
}


// /// Compute cross-correlation coelution score.
// pub fn calc_xcorr_coelution_score(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> f64 {
//     let (xcorr_matrix, _, _) = initialize_xcorr_matrix(&vec![intensities1.clone(), intensities2.clone()]);
//     let mut deltas = Vec::new();
//     for i in 0..xcorr_matrix.nrows() {
//         for j in i..xcorr_matrix.ncols() {
//             deltas.push(xcorr_matrix[(i, j)]);
//         }
//     }
//     let mean = deltas.iter().copied().sum::<f64>() / deltas.len() as f64;
//     let stddev = (deltas.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / deltas.len() as f64).sqrt();
//     mean + stddev
// }

/// Optimized cross-correlation coelution score calculation
pub fn calc_xcorr_coelution_score(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> f64 {
    // Return 0 if either input is empty
    if intensities1.is_empty() || intensities2.is_empty() {
        return 0.0;
    }

    // Direct calculation without matrix overhead
    let mean1 = intensities1.mean().unwrap_or(0.0);
    let mean2 = intensities2.mean().unwrap_or(0.0);
    let std1 = intensities1.std(0.0);
    let std2 = intensities2.std(0.0);
    
    let xcorr = normalized_cross_correlation_optimized(
        intensities1,
        intensities2,
        mean1, std1,
        mean2, std2
    );
    
    let (max_peak, _) = xcorr_array_get_max_peak_optimized(&xcorr);
    
    // Check for NaN or infinite values and return 0.0 in those cases
    if max_peak.is_nan() || max_peak.is_infinite() {
        0.0
    } else {
        max_peak
    }
}

/// Compute cross-correlation to many
// pub fn calc_xcorr_to_many_score(intensities1: &Array1<f64>, intensities2: &Vec<Array1<f64>>) -> f64 {
//     let mut total_xcorr = 0.0;
//     let mut count = 0;
//     for intensities in intensities2 {
//         let xcorr = calc_xcorr_coelution_score(&intensities1, intensities);
//         total_xcorr += xcorr;
//         count += 1;
//     }
//     total_xcorr / count as f64
// }

/// Compute cross-correlation to many
pub fn calc_xcorr_to_many_score(
    intensities1: &Array1<f64>,
    intensities2: &[Array1<f64>],
) -> f64 {
    // Return 0 if either input is empty
    if intensities1.is_empty() || intensities2.is_empty() {
        return 0.0;
    }

    // Precompute stats for intensities1 (only once!)
    let mean1 = intensities1.mean().unwrap_or(0.0);
    let std1 = intensities1.std(0.0);

    // Process intensities2 in parallel
    let (total_xcorr, count) = intensities2.par_iter()
        .map(|intensities2| {
            // Skip empty arrays
            if intensities2.is_empty() {
                return 0.0;
            }

            let mean2 = intensities2.mean().unwrap_or(0.0);
            let std2 = intensities2.std(0.0);
            
            let xcorr = normalized_cross_correlation_optimized(
                intensities1,
                intensities2,
                mean1, std1,
                mean2, std2,
            );
            
            let (_, max_peak) = xcorr_array_get_max_peak_optimized(&xcorr);
            
            // Return 0 for NaN or infinite values
            if max_peak.is_nan() || max_peak.is_infinite() {
                0.0
            } else {
                max_peak
            }
        })
        .fold(
            || (0.0, 0),
            |(sum, count), xcorr| (sum + xcorr, count + 1),
        )
        .reduce(
            || (0.0, 0),
            |(sum1, count1), (sum2, count2)| (sum1 + sum2, count1 + count2),
        );

    // Avoid division by zero and return 0 if no valid scores were computed
    if count == 0 {
        0.0
    } else {
        let avg = total_xcorr / count as f64;
        // Final check in case the average itself is NaN/inf
        if avg.is_nan() || avg.is_infinite() {
            0.0
        } else {
            avg
        }
    }
}

/// Compute cross-correlation shape score.
// pub fn calc_xcorr_shape_score(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> f64 {
//     let (xcorr_matrix, _, _) = initialize_xcorr_matrix(&vec![intensities1.clone(), intensities2.clone()]);
//     let mut sum_intensities = 0.0;
//     let mut count = 0;
//     for i in 0..xcorr_matrix.nrows() {
//         for j in i..xcorr_matrix.ncols() {
//             sum_intensities += xcorr_matrix[(i, j)];
//             count += 1;
//         }
//     }
//     sum_intensities / count as f64
// }

pub fn calc_xcorr_shape_score(
    intensities1: &Array1<f64>,
    intensities2: &Array1<f64>,
) -> f64 {
    // Early return for empty inputs
    if intensities1.is_empty() || intensities2.is_empty() {
        return 0.0;
    }

    // Precompute stats with fallbacks
    let mean1 = intensities1.mean().unwrap_or(0.0);
    let mean2 = intensities2.mean().unwrap_or(0.0);
    let std1 = intensities1.std(0.0);
    let std2 = intensities2.std(0.0);

    // Compute cross-correlation
    let xcorr = normalized_cross_correlation_optimized(
        intensities1,
        intensities2,
        mean1, std1,
        mean2, std2,
    );

    // Compute mean of valid values with comprehensive checks
    let (sum, count) = xcorr.iter()
        .filter(|&&v| !v.is_nan() && !v.is_infinite())  // Filter both NaN and infinite
        .fold((0.0, 0), |(sum, count), &v| (sum + v, count + 1));

    // Final result with multiple safeguards
    if count == 0 {
        0.0
    } else {
        let result = sum / count as f64;
        if result.is_nan() || result.is_infinite() {
            0.0  // Final sanity check
        } else {
            result
        }
    }
}

// Compute cross-correlation shape score to many
// pub fn calc_xcorr_shape_to_many_score(intensities1: &Array1<f64>, intensities2: &Vec<Array1<f64>>) -> f64 {
//     let mut total_xcorr = 0.0;
//     let mut count = 0;
//     for intensities in intensities2 {
//         let xcorr = calc_xcorr_shape_score(&intensities1, intensities);
//         total_xcorr += xcorr;
//         count += 1;
//     }
//     total_xcorr / count as f64
// }

// Compute cross-correlation shape score to many
pub fn calc_xcorr_shape_to_many_score(
    intensities1: &Array1<f64>,
    intensities2: &[Array1<f64>],
) -> f64 {
    // Early return for empty inputs
    if intensities1.is_empty() || intensities2.is_empty() {
        return 0.0;
    }

    // Precompute stats for the reference array with fallbacks
    let mean1 = intensities1.mean().unwrap_or(0.0);
    let std1 = intensities1.std(0.0);

    // Process all target arrays in parallel
    let (total_xcorr, count) = intensities2.par_iter()
        .map(|intensities2| {
            // Skip empty arrays
            if intensities2.is_empty() {
                return (0.0, 0);
            }

            // Precompute stats with fallbacks
            let mean2 = intensities2.mean().unwrap_or(0.0);
            let std2 = intensities2.std(0.0);

            // Compute cross-correlation
            let xcorr = normalized_cross_correlation_optimized(
                intensities1,
                intensities2,
                mean1, std1,
                mean2, std2
            );

            // Calculate mean of valid correlation values (excluding NaN and infinite)
            let (sum, count) = xcorr.iter()
                .filter(|&&v| v.is_finite())  // Filters out both NaN and infinite
                .fold((0.0, 0), |(sum, count), &v| (sum + v, count + 1));

            // Return partial results with validation
            if count == 0 {
                (0.0, 0)
            } else {
                let partial_avg = sum / count as f64;
                if partial_avg.is_finite() {
                    (partial_avg * count as f64, count)  // Return as sum to avoid precision loss
                } else {
                    (0.0, 0)
                }
            }
        })
        .reduce(
            || (0.0, 0),  // Identity element for reduction
            |(sum1, count1), (sum2, count2)| (sum1 + sum2, count1 + count2)
        );

    // Compute final average with multiple safeguards
    if count == 0 {
        0.0
    } else {
        let result = total_xcorr / count as f64;
        if result.is_finite() {
            result
        } else {
            0.0
        }
    }
}




// pub fn initialize_mi_matrix(intensities: &Vec<Array1<f64>>) -> Array2<f64> {
//     let n = intensities.len();
//     let mut mi_matrix = Array2::<f64>::zeros((n, n));

//     let rank_vec = compute_rank_vector(intensities);
//     let max_rank_vec = compute_max_rank_vector(&rank_vec);

//     for i in 0..n {
//         for j in i..n {
//             mi_matrix[(i, j)] = ranked_mutual_information(&rank_vec[i], &rank_vec[j], max_rank_vec[i], max_rank_vec[j]);
//         }
//     }
//     mi_matrix
// }





// /// Compute rank vectors for mutual information calculation.
// pub fn compute_rank_vector(intensities: &Vec<Array1<f64>>) -> Vec<Vec<u32>> {
//     let mut rank_vec = Vec::new();
//     for intensity in intensities {
//         let mut sorted_indices: Vec<usize> = (0..intensity.len()).collect();
//         sorted_indices.sort_by(|&a, &b| intensity[a].partial_cmp(&intensity[b]).unwrap());
//         let ranks: Vec<u32> = sorted_indices.iter().enumerate().map(|(i, _)| i as u32 + 1).collect();
//         rank_vec.push(ranks);
//     }
//     rank_vec
// }

// /// Compute the max rank vector for normalizing mutual information.
// pub fn compute_max_rank_vector(rank_vec: &Vec<Vec<u32>>) -> Vec<u32> {
//     rank_vec.iter().map(|ranks| *ranks.iter().max().unwrap_or(&1)).collect()
// }

// pub fn ranked_mutual_information(ranks1: &Vec<u32>, ranks2: &Vec<u32>, max_rank1: u32, max_rank2: u32) -> f64 {
//     let mut mi_score = 0.0;
//     let min_len = ranks1.len().min(ranks2.len()); 
//     for i in 0..min_len {
//         if ranks1[i] == ranks2[i] {
//             mi_score += 1.0;
//         }
//     }
//     mi_score / ((max_rank1 + max_rank2) as f64 / 2.0)
// }


// /// Compute mutual information score.
// pub fn calc_mi_score(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> f64 {
//     let mi_matrix = initialize_mi_matrix(&vec![intensities1.clone(), intensities2.clone()]);
//     let mi_sum = mi_matrix.sum();
//     let count = mi_matrix.nrows() * mi_matrix.ncols() / 2 + (mi_matrix.nrows() + 1) / 2;
//     mi_sum / count as f64
// }

// /// Compute mutual information to many
// pub fn calc_mi_to_many_score(intensities1: &Array1<f64>, intensities2: &Vec<Array1<f64>>) -> f64 {
//     let mut total_mi = 0.0;
//     let mut count = 0;
//     for intensities in intensities2 {
//         let mi = calc_mi_score(&intensities1, intensities);
//         total_mi += mi;
//         count += 1;
//     }
//     total_mi / count as f64
// }

/// Optimized mutual information matrix initialization
pub fn initialize_mi_matrix(intensities: &[Array1<f64>]) -> Array2<f64> {
    let n = intensities.len();
    let mut mi_matrix = Array2::<f64>::zeros((n, n));

    // Parallel rank computation
    let (rank_vec, max_rank_vec): (Vec<_>, Vec<_>) = intensities.par_iter()
        .map(|intensity| {
            let mut sorted_indices: Vec<usize> = (0..intensity.len()).collect();
            sorted_indices.sort_by(|&a, &b| intensity[a].partial_cmp(&intensity[b]).unwrap());
            let max_rank = sorted_indices.len() as u32;
            let ranks: Vec<u32> = sorted_indices.iter().enumerate()
                .map(|(i, _)| i as u32 + 1)
                .collect();
            (ranks, max_rank)
        })
        .unzip();

    // Parallel matrix filling (upper triangle)
    // Using ndarray's indexed parallel iterator instead
    mi_matrix.indexed_iter_mut().par_bridge().for_each(|((i, j), elem)| {
        if j >= i {
            *elem = ranked_mutual_information_optimized(
                &rank_vec[i],
                &rank_vec[j],
                max_rank_vec[i],
                max_rank_vec[j]
            );
        }
    });

    // Fill lower triangle (symmetric)
    for i in 0..n {
        for j in 0..i {
            mi_matrix[(i, j)] = mi_matrix[(j, i)];
        }
    }

    mi_matrix
}

/// Optimized ranked mutual information calculation
pub fn ranked_mutual_information_optimized(
    ranks1: &[u32],
    ranks2: &[u32],
    max_rank1: u32,
    max_rank2: u32,
) -> f64 {
    let min_len = ranks1.len().min(ranks2.len());
    let matches = ranks1.iter().zip(ranks2.iter())
        .take(min_len)
        .filter(|(&r1, &r2)| r1 == r2)
        .count() as f64;
    
    matches / ((max_rank1 + max_rank2) as f64 / 2.0)
}

/// Optimized mutual information score
pub fn calc_mi_score(
    intensities1: &Array1<f64>,
    intensities2: &Array1<f64>,
) -> f64 {
    // Early return for empty inputs
    if intensities1.is_empty() || intensities2.is_empty() {
        return 0.0;
    }

    // Compute ranks with validation
    let (rank1, max_rank1) = match compute_rank_and_max(intensities1) {
        (r, m) if r.is_empty() => return 0.0,
        result => result,
    };
    let (rank2, max_rank2) = match compute_rank_and_max(intensities2) {
        (r, m) if r.is_empty() => return 0.0,
        result => result,
    };

    // Calculate mutual information with validation
    let mi = ranked_mutual_information_optimized(&rank1, &rank2, max_rank1, max_rank2);
    
    // Validate result
    if mi.is_nan() || mi.is_infinite() || mi < 0.0 {
        0.0
    } else {
        mi
    }
}

/// Optimized mutual information to many
pub fn calc_mi_to_many_score(
    intensities1: &Array1<f64>,
    intensities2: &[Array1<f64>],
) -> f64 {
    // Early return for empty inputs
    if intensities1.is_empty() || intensities2.is_empty() {
        return 0.0;
    }

    // Precompute reference ranks with validation
    let (rank1, max_rank1) = match compute_rank_and_max(intensities1) {
        (r, m) if r.is_empty() => return 0.0,
        result => result,
    };

    // Process in parallel with full validation
    let (total_mi, count) = intensities2.par_iter()
        .filter_map(|intensity| {
            // Skip empty arrays
            if intensity.is_empty() {
                return None;
            }

            // Compute ranks with validation
            let (rank2, max_rank2) = match compute_rank_and_max(intensity) {
                (r, m) if r.is_empty() => return None,
                result => result,
            };

            // Calculate MI with validation
            let mi = ranked_mutual_information_optimized(&rank1, &rank2, max_rank1, max_rank2);
            
            // Only include valid MI scores
            if mi.is_finite() && mi >= 0.0 {
                Some(mi)
            } else {
                None
            }
        })
        .fold(
            || (0.0, 0),
            |(sum, count), mi| (sum + mi, count + 1),
        )
        .reduce(
            || (0.0, 0),
            |(sum1, count1), (sum2, count2)| (sum1 + sum2, count1 + count2),
        );

    // Final validation
    if count == 0 {
        0.0
    } else {
        let avg = total_mi / count as f64;
        if avg.is_finite() && avg >= 0.0 {
            avg
        } else {
            0.0
        }
    }
}

/// Helper function to compute rank and max rank together
fn compute_rank_and_max(intensity: &Array1<f64>) -> (Vec<u32>, u32) {
    let mut sorted_indices: Vec<usize> = (0..intensity.len()).collect();
    sorted_indices.sort_by(|&a, &b| intensity[a].partial_cmp(&intensity[b]).unwrap());
    let max_rank = sorted_indices.len() as u32;
    let ranks = sorted_indices.iter().enumerate()
        .map(|(i, _)| i as u32 + 1)
        .collect();
    (ranks, max_rank)
}