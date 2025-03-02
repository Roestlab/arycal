use ndarray::{Array1, Array2, s};


pub fn initialize_xcorr_matrix(intensities: &Vec<Array1<f64>>) -> (Array2<f64>, Array2<f64>, Array2<f64>) {
    let n = intensities.len();
    let mut xcorr_matrix = Array2::<f64>::zeros((n, n));
    let mut xcorr_matrix_max_peak = Array2::<f64>::zeros((n, n));
    let mut xcorr_matrix_max_peak_sec = Array2::<f64>::zeros((n, n));

    for i in 0..n {
        for j in i..n {
            let xcorr = normalized_cross_correlation(&intensities[i], &intensities[j]);
            let (max_peak, max_peak_sec) = xcorr_array_get_max_peak(&xcorr);
            
            // Store only the max value of the cross-correlation in the matrix
            xcorr_matrix[(i, j)] = max_peak;  
            xcorr_matrix_max_peak[(i, j)] = max_peak;
            xcorr_matrix_max_peak_sec[(i, j)] = max_peak_sec;
        }
    }
    (xcorr_matrix, xcorr_matrix_max_peak, xcorr_matrix_max_peak_sec)
}


pub fn initialize_mi_matrix(intensities: &Vec<Array1<f64>>) -> Array2<f64> {
    let n = intensities.len();
    let mut mi_matrix = Array2::<f64>::zeros((n, n));

    let rank_vec = compute_rank_vector(intensities);
    let max_rank_vec = compute_max_rank_vector(&rank_vec);

    for i in 0..n {
        for j in i..n {
            mi_matrix[(i, j)] = ranked_mutual_information(&rank_vec[i], &rank_vec[j], max_rank_vec[i], max_rank_vec[j]);
        }
    }
    mi_matrix
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
    let mean1 = intensities1.mean().unwrap_or(0.0);
    let mean2 = intensities2.mean().unwrap_or(0.0);
    let std1 = intensities1.std(0.0);
    let std2 = intensities2.std(0.0);

    let norm1 = (intensities1 - mean1) / std1;
    let norm2 = (intensities2 - mean2) / std2;

    let mut xcorr = Array1::<f64>::zeros(intensities1.len());
    for i in 0..intensities1.len() {
        let slice1 = norm1.slice(s![..(intensities1.len() - i)]);
        let slice2 = norm2.slice(s![i..]);

        // Compute dot product manually
        let mut dot_product = 0.0;
        for (a, b) in slice1.iter().zip(slice2.iter()) {
            dot_product += a * b;
        }

        xcorr[i] = dot_product;
    }
    xcorr
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

/// Compute rank vectors for mutual information calculation.
pub fn compute_rank_vector(intensities: &Vec<Array1<f64>>) -> Vec<Vec<u32>> {
    let mut rank_vec = Vec::new();
    for intensity in intensities {
        let mut sorted_indices: Vec<usize> = (0..intensity.len()).collect();
        sorted_indices.sort_by(|&a, &b| intensity[a].partial_cmp(&intensity[b]).unwrap());
        let ranks: Vec<u32> = sorted_indices.iter().enumerate().map(|(i, _)| i as u32 + 1).collect();
        rank_vec.push(ranks);
    }
    rank_vec
}

/// Compute the max rank vector for normalizing mutual information.
pub fn compute_max_rank_vector(rank_vec: &Vec<Vec<u32>>) -> Vec<u32> {
    rank_vec.iter().map(|ranks| *ranks.iter().max().unwrap_or(&1)).collect()
}

/// Compute ranked mutual information between two ranked intensity vectors.
pub fn ranked_mutual_information(ranks1: &Vec<u32>, ranks2: &Vec<u32>, max_rank1: u32, max_rank2: u32) -> f64 {
    let mut mi_score = 0.0;
    for i in 0..ranks1.len() {
        if ranks1[i] == ranks2[i] {
            mi_score += 1.0;
        }
    }
    mi_score / ((max_rank1 + max_rank2) as f64 / 2.0)
}

/// Compute cross-correlation coelution score.
pub fn calc_xcorr_coelution_score(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> f64 {
    let (xcorr_matrix, _, _) = initialize_xcorr_matrix(&vec![intensities1.clone(), intensities2.clone()]);
    let mut deltas = Vec::new();
    for i in 0..xcorr_matrix.nrows() {
        for j in i..xcorr_matrix.ncols() {
            deltas.push(xcorr_matrix[(i, j)]);
        }
    }
    let mean = deltas.iter().copied().sum::<f64>() / deltas.len() as f64;
    let stddev = (deltas.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / deltas.len() as f64).sqrt();
    mean + stddev
}

/// Compute cross-correlation to many
pub fn calc_xcorr_to_many_score(intensities1: &Array1<f64>, intensities2: &Vec<Array1<f64>>) -> f64 {
    let mut total_xcorr = 0.0;
    let mut count = 0;
    for intensities in intensities2 {
        let xcorr = calc_xcorr_coelution_score(&intensities1, intensities);
        total_xcorr += xcorr;
        count += 1;
    }
    total_xcorr / count as f64
}

/// Compute cross-correlation shape score.
pub fn calc_xcorr_shape_score(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> f64 {
    let (xcorr_matrix, _, _) = initialize_xcorr_matrix(&vec![intensities1.clone(), intensities2.clone()]);
    let mut sum_intensities = 0.0;
    let mut count = 0;
    for i in 0..xcorr_matrix.nrows() {
        for j in i..xcorr_matrix.ncols() {
            sum_intensities += xcorr_matrix[(i, j)];
            count += 1;
        }
    }
    sum_intensities / count as f64
}

// Compute cross-correlation shape score to many
pub fn calc_xcorr_shape_to_many_score(intensities1: &Array1<f64>, intensities2: &Vec<Array1<f64>>) -> f64 {
    let mut total_xcorr = 0.0;
    let mut count = 0;
    for intensities in intensities2 {
        let xcorr = calc_xcorr_shape_score(&intensities1, intensities);
        total_xcorr += xcorr;
        count += 1;
    }
    total_xcorr / count as f64
}

/// Compute mutual information score.
pub fn calc_mi_score(intensities1: &Array1<f64>, intensities2: &Array1<f64>) -> f64 {
    let mi_matrix = initialize_mi_matrix(&vec![intensities1.clone(), intensities2.clone()]);
    let mi_sum = mi_matrix.sum();
    let count = mi_matrix.nrows() * mi_matrix.ncols() / 2 + (mi_matrix.nrows() + 1) / 2;
    mi_sum / count as f64
}

/// Compute mutual information to many
pub fn calc_mi_to_many_score(intensities1: &Array1<f64>, intensities2: &Vec<Array1<f64>>) -> f64 {
    let mut total_mi = 0.0;
    let mut count = 0;
    for intensities in intensities2 {
        let mi = calc_mi_score(&intensities1, intensities);
        total_mi += mi;
        count += 1;
    }
    total_mi / count as f64
}