pub mod input;
pub mod output;

#[cfg(feature = "mpi")]
use mpi::traits::*;

use anyhow::Result;
use log::info;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::time::Instant;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};

use arycal_cloudpath::{
    tsv::load_precursor_ids_from_tsv,
    osw::{OswAccess, PrecursorIdData},
    sqmass::{create_common_rt_space, SqMassAccess},
};
use arycal_common::{logging::Progress, AlignedTransitionScores, ArycalError, FullTraceAlignmentScores, PeakMapping, PrecursorAlignmentResult};
use arycal_core::{alignment::{self, alignment::{apply_post_alignment_to_trgrp, AlignedChromatogram}, dynamic_time_warping::align_chromatograms, fast_fourier_lag::shift_chromatogram}, scoring::{compute_alignment_scores, compute_peak_mapping_scores, compute_peak_mapping_transitions_scores}};
use arycal_core::{
    alignment::alignment::map_peaks_across_runs,
    alignment::dynamic_time_warping::{star_align_tics, mst_align_tics, progressive_align_tics},
    alignment::fast_fourier_lag::{star_align_tics_fft, mst_align_tics_fft, progressive_align_tics_fft},
    alignment::fast_fourier_lag_dtw::star_align_tics_fft_with_local_refinement,
    scoring::{create_decoy_peaks_by_random_regions, create_decoy_peaks_by_shuffling},
};
use input::Input; 

pub struct Runner {
    precursor_map: Vec<PrecursorIdData>,
    parameters: input::Input,
    feature_access: Vec<OswAccess>,
    xic_access: Vec<SqMassAccess>,
    start: Instant,
    progress_num: Option<Arc<Mutex<f32>>>, 
}

impl Runner {
    pub fn new(parameters: Input, progress_num: Option<Arc<Mutex<f32>>>) -> anyhow::Result<Self> {
        let start = Instant::now();

        // TODO: Currently only supports a single OSW file
        let start_io = Instant::now();
        let osw_access = OswAccess::new(&parameters.features.file_paths[0].to_str().unwrap())?;

        // Check if precursor_ids tsv file is provided
        let mut precursor_ids: Option<Vec<u32>> = None;
        if let Some(precursor_ids_file) = &parameters.filters.precursor_ids {
            precursor_ids = Some(load_precursor_ids_from_tsv(precursor_ids_file)?);
        }

        let precursor_map: Vec<PrecursorIdData> = osw_access.fetch_transition_ids(parameters.filters.decoy, parameters.filters.include_identifying_transitions.unwrap_or_default(), precursor_ids)?;
        let run_time = (Instant::now() - start_io).as_millis();

        info!(
            "Loaded {} target precursors and {} decoy precursors identifiers - took {}ms",
            precursor_map.iter().filter(|v| !v.decoy).count(),
            precursor_map.iter().filter(|v| v.decoy).count(),
            run_time
        );

        let xic_accessors: Result<Vec<SqMassAccess>, anyhow::Error> = parameters
            .xic
            .file_paths
            .iter()
            .map(|path| SqMassAccess::new(path.to_str().unwrap()).map_err(anyhow::Error::from))
            .collect();

        Ok(Self {
            precursor_map: precursor_map.clone(),
            parameters,
            feature_access: vec![osw_access],
            xic_access: xic_accessors?,
            start,
            progress_num: Some(progress_num.expect("Progress number is not set")),
        })
    }

    #[cfg(not(feature = "mpi"))]
    pub fn run(&mut self) -> anyhow::Result<()> {

        // Print log info of alignment params from self.parameters.alignment
        log::debug!("{}", self.parameters.alignment);

        let precursor_map = &self.precursor_map; 

        // let results: Vec<
        //     Result<HashMap<i32, PrecursorAlignmentResult>, ArycalError>,
        // > = precursor_map
        //     .par_chunks(self.parameters.alignment.batch_size.unwrap_or_default()) 
        //     .map(|batch| {
        //         let mut batch_results = Vec::new();
        //         let mut progress = None;
        //         if !self.parameters.disable_progress_bar {
        //             progress = Some(Progress::new(
        //                 batch.len(),
        //                 format!(
        //                     "[arycal] Thread {:?} - Aligning precursors",
        //                     rayon::current_thread_index().unwrap_or(0)
        //                 )
        //                 .as_str(),
        //             ));
        //         }
        //         let start_time = Instant::now();
        //         for precursor in batch {
        //             let result = self.process_precursor(precursor).map_err(|e| ArycalError::Custom(e.to_string()));
        //             batch_results.push(result);
        //             if !self.parameters.disable_progress_bar {
        //                 progress.as_ref().expect("The Progess tqdm logger is not enabled").inc();
        //             }
        //             // Update the shared progress if progress_num is Some
        //             if let Some(progress_num) = &self.progress_num {
        //                 let mut progress_num = progress_num.lock().unwrap();
        //                 *progress_num += 1.0 / precursor_map.len() as f32;
        //             }
        //         }
        //         let end_time = Instant::now();
        //         if self.parameters.disable_progress_bar {
        //             log::info!("Batch of {} precursors aligned in {:?}", batch.len(), end_time.duration_since(start_time));
        //         }
        //         Ok(batch_results)
        //     })
        //     .collect::<Result<Vec<_>, ArycalError>>()? 
        //     .into_iter()
        //     .flatten() 
        //     .collect();

        // Create a custom thread pool with limited parallelism
        let pool = ThreadPoolBuilder::new()
        .num_threads(self.parameters.alignment.batch_size.unwrap_or(rayon::max_num_threads() - 4)) 
        .build()
        .unwrap();

        let processed_count = Arc::new(AtomicUsize::new(0));
        let batch_size = self.parameters.alignment.batch_size.unwrap_or(rayon::max_num_threads() - 4);
        let total_count = precursor_map.len();
        let global_start = Instant::now();
        let batch_start = Arc::new(std::sync::Mutex::new(Instant::now())); // For batch timing

        log::info!("Will align {} precursors in parallel at a time", batch_size);
        log::info!("Starting alignment for {} precursors", total_count);

        let results = pool.install(|| {
            precursor_map
                .par_iter()
                .map(|precursor| {
                    let result = self.process_precursor(precursor)
                        .map_err(|e| ArycalError::Custom(e.to_string()));

                    // Update progress
                    if let Some(progress_num) = &self.progress_num {
                        let mut progress_num = progress_num.lock().unwrap();
                        *progress_num += 1.0 / total_count as f32;
                    }

                    // Get and increment count
                    let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;

                    // Log timing for every 500 precursors
                    if count % 500 == 1 {
                        // First item in new batch - reset timer
                        *batch_start.lock().unwrap() = Instant::now();
                    }
                    else if count % 500 == 0 || count == total_count {
                        // Last item in batch - log batch duration
                        let batch_time = batch_start.lock().unwrap().elapsed();
                        log::info!(
                            "Batch of {} precursors aligned and scored in {:?} ({:.4}/min) | Total: {}/{}",
                            count.min(500), // Handle final partial batch
                            batch_time,
                            500 as f64 / (batch_time.as_secs_f64() / 60.0),
                            count,
                            total_count
                        );
                    }

                    result
                })
                .collect()
        });

        // Final stats
        let total_elapsed = global_start.elapsed();
        log::info!(
            "Aligned and scored {} precursors in {:?} ({:.2}/sec)",
            total_count,
            total_elapsed,
            total_count as f64 / total_elapsed.as_secs_f64()
        );

        // // Filter precursor map for modified sequence = LTPEAIR and charge = 2
        // let filtered_precursors: Vec<PrecursorIdData> = self.precursor_map
        //     .iter()
        //     .filter(|p| p.modified_sequence == "LTPEAIR" && p.precursor_charge == 2)
        //     .cloned()
        //     .collect();

        // // Process the filtered precursors
        // let results: Vec<Result<HashMap<i32, HashMap<std::string::String, Vec<PeakMapping>>>, anyhow::Error>> = filtered_precursors
        //     .par_iter()
        //     .map(|precursor| {
        //         self.process_precursor(&precursor)
        //     })
        //     .collect();

        // Write out to separate file if scores_output_file is provided
        let feature_access = if let Some(scores_output_file) = &self.parameters.alignment.scores_output_file {
            let scores_output_file = scores_output_file.clone();
            let osw_access = OswAccess::new(&scores_output_file)?;
            vec![osw_access]
        } else {
            self.feature_access.clone()
        };
        

        if self.parameters.alignment.compute_scores.unwrap_or_default() {
            // Write FEATURE_ALIGNMENT results to the database
            self.write_aligned_score_results_to_db(&feature_access, &results)?;
        }
        
        // Write FEATURE_MS2_ALIGNMENT results to the database
        self.write_ms2_alignment_results_to_db(&feature_access, &results)?;

        // Write FEATURE_TRANSITION_ALIGNMENT results to the database
        if self.parameters.filters.include_identifying_transitions.unwrap_or_default() && self.parameters.alignment.compute_scores.unwrap_or_default() {
            self.write_transition_alignment_results_to_db(&feature_access, &results)?;   
        }

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);
        Ok(())
    }

    #[cfg(feature = "mpi")]
    pub fn run(&self) -> anyhow::Result<()> {
        // Initialize MPI
        let universe = mpi::initialize().unwrap();
        let world = universe.world();
        let rank = world.rank();  // Current process rank
        let size = world.size();  // Total number of processes

        let hostname = std::env::var("HOSTNAME").unwrap_or_else(|_| "unknown".to_string());
        log::info!("MPI initialized: rank = {}, size = {}, hostname = {}", rank, size, hostname);

        // Determine chunking for distributing precursors
        let total_precursors = self.precursor_map.len();
        let chunk_size = (total_precursors + size as usize - 1) / size as usize; // Ensure even distribution
        let start = rank as usize * chunk_size;
        let end = std::cmp::min(start + chunk_size, total_precursors);

        let local_chunk = &self.precursor_map[start..end];
        log::info!(
            "Process {}: Processing chunk {} to {} ({} precursors)",
            rank, start, end, local_chunk.len()
        );

        // Parallel processing using Rayon within each node
        let local_results: Vec<Result<HashMap<i32, PrecursorAlignmentResult>, ArycalError>> = local_chunk
            .par_chunks(self.parameters.alignment.batch_size.unwrap_or_default())
            .map(|batch| {
                let mut batch_results = Vec::new();
                let mut progress = None;
                if !self.parameters.disable_progress_bar {
                    progress = Some(Progress::new(
                        batch.len(),
                        format!(
                            "[arycal] Node {} - Rank {} - Thread {:?} - Aligning precursors",
                            hostname,
                            rank,
                            rayon::current_thread_index().unwrap()
                        )
                        .as_str(),
                    ));
                }
                let start_time = Instant::now();
                for precursor in batch {
                    let result = self.process_precursor(precursor).map_err(|e| ArycalError::Custom(e.to_string()));
                    batch_results.push(result);
                    if !self.parameters.disable_progress_bar {
                        progress.as_ref().expect("The Progess tqdm logger is not enabled").inc();
                    }
                }
                let end_time = Instant::now();
                if self.parameters.disable_progress_bar {
                    log::info!("Node {} - Rank {} - Thread {:?} - Batch of {} precursors aligned in {:?}", hostname, rank, rayon::current_thread_index().unwrap_or(0), batch.len(), end_time.duration_since(start_time));
                }
                Ok(batch_results)
            })
            .collect::<Result<Vec<_>, ArycalError>>()?
            .into_iter()
            .flatten()
            .collect();

        // Serialize results
        let serialized_results = bincode::serialize(&local_results)?;

        // Gather results at root (rank 0)
        let gathered_results = if rank == 0 {
            let mut all_results = Vec::new();
            all_results.extend(local_results);

            for process_rank in 1..size {
                log::info!("Root process receiving results from process {}", process_rank);
                let (received_bytes, _status) = world.process_at_rank(process_rank).receive_vec::<u8>();
                let received_results: Vec<Result<HashMap<i32, PrecursorAlignmentResult>, ArycalError>> =
                    bincode::deserialize(&received_bytes)?;
                all_results.extend(received_results);
            }

            all_results
        } else {
            log::info!("Process {} sending results to root process", rank);
            world.process_at_rank(0).send(&serialized_results);
            Vec::new() // Non-root processes do not keep results
        };

        // Only root process writes results to the database
        if rank == 0 {

            // Write out to separate file if scores_output_file is provided
            let feature_access = if let Some(scores_output_file) = &self.parameters.alignment.scores_output_file {
                let scores_output_file = scores_output_file.clone();
                let osw_access = OswAccess::new(&scores_output_file)?;
                vec![osw_access]
            } else {
                self.feature_access.clone()
            };

            if self.parameters.alignment.compute_scores.unwrap_or_default() {
                // Write FEATURE_ALIGNMENT results to the database
                self.write_aligned_score_results_to_db(&feature_access, &gathered_results)?;
            }
            self.write_ms2_alignment_results_to_db(&feature_access, &gathered_results)?;

            if self.parameters.filters.include_identifying_transitions.unwrap_or_default() && self.parameters.alignment.compute_scores.unwrap_or_default() {
                self.write_transition_alignment_results_to_db(&feature_access, &gathered_results)?;   
            }

            let run_time = (Instant::now() - self.start).as_secs();
            info!("Alignment completed in {}s", run_time);
        }

        Ok(())
    }


    pub fn process_precursor(
        &self,
        precursor: &PrecursorIdData,
    ) -> Result<HashMap<i32, PrecursorAlignmentResult>, anyhow::Error> {
        let native_ids: Vec<String> = precursor.clone().extract_native_ids_for_sqmass(
            self.parameters.xic.include_precursor,
            self.parameters.xic.num_isotopes,
        );
        let native_ids_str: Vec<&str> = native_ids.iter().map(|s| s.as_str()).collect();

        log::trace!("modified_sequence: {:?}, precursor_charge: {:?}, detecting transitions: {:?}, identifying transitions: {:?}", precursor.modified_sequence, precursor.precursor_charge, precursor.n_transitions(), precursor.n_identifying_transitions());

        // log::trace!("native_ids: {:?}", native_ids);

        let group_id =
            precursor.modified_sequence.clone() + "_" + &precursor.precursor_charge.to_string();
    
        /* ------------------------------------------------------------------ */
        /* Step 1. Extract and transform XICs                                 */
        /* ------------------------------------------------------------------ */

        // Extract chromatograms from the XIC files
        let chromatograms: Vec<_> = self
            .xic_access
            .iter()
            .map(|access| {
                access.read_chromatograms("NATIVE_ID", native_ids_str.clone(), group_id.clone())
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Check length of the first chromatogram, should be at least more than 10 points
        if chromatograms[0]
            .chromatograms
            .iter()
            .map(|chromatogram| chromatogram.1.intensities.len())
            .sum::<usize>()
            < 10
        {
            log::trace!("The first chromatogram has less than 10 points, skipping precursor: {:?}", precursor.precursor_id);
            return Ok(HashMap::new());
        }

        // Check that there are no NaN values in the chromatograms
        for chrom in chromatograms.iter() {
            for (_, chrom_data) in chrom.chromatograms.iter() {
                if chrom_data.intensities.iter().any(|&x| x.is_nan()) || chrom_data.retention_times.iter().any(|&x| x.is_nan()) {
                    log::trace!("NaN values detected in chromatograms, skipping precursor: {:?}", precursor.precursor_id);
                    return Ok(HashMap::new());
                }
            }
        }

        // Compute TICs
        let tics: Vec<_> = chromatograms
            .iter()
            .map(|chromatogram| chromatogram.calculate_tic())
            .collect();

        // Create common retention time space
        let common_rt_space = create_common_rt_space(tics);
        // let common_rt_space = tics.clone();

        // Smooth and normalize TICs
        let smoothed_tics: Vec<_> = common_rt_space
            .iter()
            .map(|tic| {
                tic.smooth_sgolay(
                    self.parameters.alignment.smoothing.sgolay_window,
                    self.parameters.alignment.smoothing.sgolay_order,
                )?
                .normalize()
            })
            .collect::<Result<Vec<_>, _>>()?;

        /* ------------------------------------------------------------------ */
        /* Step 2. Pair-wise Alignment of TICs                                */
        /* ------------------------------------------------------------------ */

        log::debug!("Aligning TICs using {:?} using reference type: {:?}", self.parameters.alignment.method.as_str(), self.parameters.alignment.reference_type);
        let start_time = Instant::now();
        let aligned_chromatograms = match self.parameters.alignment.method.to_lowercase().as_str() {
            "dtw" => {
                match self.parameters.alignment.reference_type.as_str() {
                    "star" => star_align_tics(smoothed_tics.clone(), &self.parameters.alignment)?,
                    "mst" => mst_align_tics(smoothed_tics.clone())?,
                    "progressive" => progressive_align_tics(smoothed_tics.clone())?,
                    _ => star_align_tics(smoothed_tics.clone(), &self.parameters.alignment)?,
                }
            },
            "fft" => {
                match self.parameters.alignment.reference_type.as_str() {
                    "star" => star_align_tics_fft(smoothed_tics.clone(), &self.parameters.alignment)?,
                    "mst" => mst_align_tics_fft(smoothed_tics.clone())?,
                    "progressive" => progressive_align_tics_fft(smoothed_tics.clone())?,
                    _ => star_align_tics_fft(smoothed_tics.clone(), &self.parameters.alignment)?,
                }
            },
            "fftdtw" => star_align_tics_fft_with_local_refinement(smoothed_tics.clone(), &self.parameters.alignment)?,
            _ => star_align_tics(smoothed_tics.clone(), &self.parameters.alignment)?,
        };
        log::debug!("Alignment took: {:?}", start_time.elapsed());

        /* ------------------------------------------------------------------ */
        /* Step 3. Score Algined TICs                                         */
        /* ------------------------------------------------------------------ */
        let alignment_scores = if self.parameters.alignment.compute_scores.unwrap_or_default() {
            log::debug!("Computing full trace alignment scores");
            let start_time = Instant::now();
            let alignment_scores = compute_alignment_scores(aligned_chromatograms.clone());
            log::debug!("Scoring took: {:?}", start_time.elapsed());
            alignment_scores
        } else {
            HashMap::new()
        };

        // let output_path = "aligned_chromatograms.parquet";
        // println!("Writing aligned chromatograms to: {:?}", output_path);
        // output::write_aligned_chromatograms_to_parquet(&aligned_chromatograms.clone(), output_path)?;

        /* ------------------------------------------------------------------ */
        /* Step 4. Aligned Peak Mapping                                       */
        /* ------------------------------------------------------------------ */
        let start_time = Instant::now();
        // fetch feature data from the database
        // TODO: Currently only supports a single merged OSW file
        let prec_feat_data = self.feature_access[0]
            .fetch_full_precursor_feature_data_for_runs(
                precursor.precursor_id,
                common_rt_space
                    .clone()
                    .iter()
                    .map(|tic| tic.metadata.get("basename").unwrap().to_string())
                    .collect(),
            )?;

        // First collect all the mapping results in parallel
        let peak_mapping_results: Vec<(String, Vec<arycal_common::PeakMapping>)> = aligned_chromatograms
        .par_iter()
        .filter_map(|chrom| {
            let current_run = chrom.chromatogram.metadata.get("basename").unwrap();
            log::trace!(
                "Mapping peaks from reference run: {} to current run: {}",
                chrom.rt_mapping[0].get("run1").unwrap(),
                current_run
            );

            // Filter prec_feat_data for current run
            let current_run_feat_data: Vec<_> = prec_feat_data
                .iter()
                .filter(|f| &f.basename == current_run)
                .cloned()
                .collect();

            // Get reference run feature data
            let ref_run_feat_data: Vec<_> = prec_feat_data
                .iter()
                .filter(|f| &f.basename == chrom.rt_mapping[0].get("run1").unwrap())
                .cloned()
                .collect();

            if current_run_feat_data.is_empty() || ref_run_feat_data.is_empty() {
                log::trace!(
                    "Current run feature data or reference run feature data is empty, skipping peak mapping for run: {}",
                    current_run
                );
                return None;
            }

            let mapped_peaks = map_peaks_across_runs(
                chrom,
                ref_run_feat_data,
                current_run_feat_data,
                self.parameters.alignment.rt_mapping_tolerance.unwrap_or_default(),
                &self.parameters.alignment,
            );

            Some((
                chrom.chromatogram.metadata.get("basename").unwrap().to_string(),
                mapped_peaks,
            ))
        })
        .collect();

        // Then insert into HashMap serially
        let mut mapped_prec_peaks: HashMap<String, Vec<arycal_common::PeakMapping>> = HashMap::new();
        for (key, value) in peak_mapping_results {
            mapped_prec_peaks.insert(key, value);
        }

        log::debug!("Peak mapping took: {:?}", start_time.elapsed());

        /* ------------------------------------------------------------------ */
        /* Step 5. Score Aligned Peaks                                        */
        /* ------------------------------------------------------------------ */
        let (scored_peak_mappings, all_peak_mappings) = if self.parameters.alignment.compute_scores.unwrap_or_default() {
            log::debug!("Computing peak mapping scores");
            let start_time = Instant::now();
            let scored_peak_mappings =
                compute_peak_mapping_scores(aligned_chromatograms.clone(), mapped_prec_peaks.clone());

            // Create decoy aligned peaks based on the method specified in the parameters
            let mut decoy_peak_mappings: HashMap<String, Vec<PeakMapping>> = HashMap::new();
            if self.parameters.alignment.decoy_peak_mapping_method == "shuffle" {
                log::debug!("Creating decoy peaks by shuffling query peaks");
                decoy_peak_mappings = create_decoy_peaks_by_shuffling(&mapped_prec_peaks.clone());
            } else if self.parameters.alignment.decoy_peak_mapping_method == "random_regions" {
                log::debug!("Creating decoy peaks by picking random regions in the query XIC");
                decoy_peak_mappings = create_decoy_peaks_by_random_regions(&aligned_chromatograms.clone(), &mapped_prec_peaks.clone(), self.parameters.alignment.decoy_window_size.unwrap_or_default());
            }
            log::debug!("Computing peak mapping scores for decoy peaks");
            let scored_decoy_peak_mappings =
                compute_peak_mapping_scores(aligned_chromatograms.clone(), decoy_peak_mappings.clone());

            // Combine true and decoy peaks for analysis into HashMap<String, Vec<PeakMapping>>
            let all_peak_mappings: HashMap<String, Vec<PeakMapping>> = {
                let mut all_peak_mappings = HashMap::new();
                for (key, value) in scored_peak_mappings
                    .iter()
                    .chain(scored_decoy_peak_mappings.iter())
                {
                    all_peak_mappings
                        .entry(key.clone())
                        .or_insert_with(Vec::new)
                        .extend(value.clone());
                }
                all_peak_mappings
            };
            log::debug!("Peak mapping scoring took: {:?}", start_time.elapsed());
            (scored_peak_mappings, all_peak_mappings)
        } else {
            (HashMap::new(), mapped_prec_peaks)
        };

        /* ------------------------------------------------------------------ */
        /* Step 6. Optional Step: Align and Score Identifying Transitions     */
        /* ------------------------------------------------------------------ */
        let identifying_peak_mapping_scores: HashMap<String, Vec<AlignedTransitionScores>> = if self.parameters.filters.include_identifying_transitions.unwrap_or_default() && self.parameters.alignment.compute_scores.unwrap_or_default() {
            let start_time = Instant::now();
            log::debug!("Processing identifying transitions - aligning and scoring");
            let id_peak_scores = self.process_identifying_transitions(group_id.clone(), precursor, aligned_chromatograms.clone(), scored_peak_mappings.clone(), smoothed_tics[0].retention_times.clone());
            log::debug!("Identifying peak mapping scoring took: {:?}", start_time.elapsed());
            id_peak_scores
        } else {
            HashMap::new()
        };
        
        // output::write_mapped_peaks_to_parquet(all_peak_mappings, "mapped_peaks.parquet")?;

        let mut result = HashMap::new();
        result.insert(precursor.precursor_id.clone(), PrecursorAlignmentResult{
            alignment_scores,
            detecting_peak_mappings: all_peak_mappings,
            identifying_peak_mapping_scores,
        });

        Ok(result)
    }

    fn process_identifying_transitions(
        &self,
        group_id: String,
        precursor: &PrecursorIdData,
        aligned_chromatograms: Vec<AlignedChromatogram>,
        peak_mappings: HashMap<String, Vec<PeakMapping>>,
        common_rt_space: Vec<f64>
    ) -> HashMap<String, Vec<AlignedTransitionScores>> {
        // Extract identifying transition ids
        let identifying_transitions_ids: Vec<String> = precursor.clone().extract_identifying_native_ids_for_sqmass();
        let identifying_transitions_ids_str: Vec<&str> = identifying_transitions_ids.iter().map(|s| s.as_str()).collect();
        log::trace!("identifying_transitions_ids: {:?}", identifying_transitions_ids);

        // Extract chromatograms for identifying transitions
        let identifying_chromatograms: Vec<_> = self
            .xic_access
            .iter()
            .map(|access| {
                access.read_chromatograms("NATIVE_ID", identifying_transitions_ids_str.clone(), group_id.clone())
            })
            .collect::<Result<Vec<_>, _>>().unwrap_or(Vec::new());

        // Check if identifying_chromatograms is empty
        if identifying_chromatograms.is_empty() {
            log::trace!("Identifying chromatograms are empty");
            return HashMap::new();
        }
        
        // Score Identifying transitions
        let aligned_identifying_trgrps = apply_post_alignment_to_trgrp(identifying_chromatograms, aligned_chromatograms.clone(), common_rt_space, &self.parameters.alignment);

        let scored_aligned_identifying_transitions = compute_peak_mapping_transitions_scores(aligned_identifying_trgrps, aligned_chromatograms, peak_mappings);

        scored_aligned_identifying_transitions
    }

    fn align_precursor(
        &self,
        precursor: &PrecursorIdData,
    ) -> Result<Vec<AlignedChromatogram>, anyhow::Error> 
    {
        let native_ids: Vec<String> = precursor.clone().extract_native_ids_for_sqmass(
            self.parameters.xic.include_precursor,
            self.parameters.xic.num_isotopes,
        );
        let native_ids_str: Vec<&str> = native_ids.iter().map(|s| s.as_str()).collect();

        log::debug!("modified_sequence: {:?}, precursor_charge: {:?}, detecting transitions: {:?}, identifying transitions: {:?}", precursor.modified_sequence, precursor.precursor_charge, precursor.n_transitions(), precursor.n_identifying_transitions());

        log::trace!("native_ids: {:?}", native_ids);

        let group_id =
            precursor.modified_sequence.clone() + "_" + &precursor.precursor_charge.to_string();
    
        /* ------------------------------------------------------------------ */
        /* Step 1. Extract and transform XICs                                 */
        /* ------------------------------------------------------------------ */

        // Extract chromatograms from the XIC files
        let chromatograms: Vec<_> = self
            .xic_access
            .iter()
            .map(|access| {
                access.read_chromatograms("NATIVE_ID", native_ids_str.clone(), group_id.clone())
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Check length of the first chromatogram, should be at least more than 10 points
        if chromatograms[0]
            .chromatograms
            .iter()
            .map(|chromatogram| chromatogram.1.intensities.len())
            .sum::<usize>()
            < 10
        {
            return Ok(Vec::new());
        }

        // Check that there are no NaN values in the chromatograms
        for chrom in chromatograms.iter() {
            for (_, chrom_data) in chrom.chromatograms.iter() {
                if chrom_data.intensities.iter().any(|&x| x.is_nan()) || chrom_data.retention_times.iter().any(|&x| x.is_nan()) {
                    return Ok(Vec::new());
                }
            }
        }

        // Compute TICs
        let tics: Vec<_> = chromatograms
            .iter()
            .map(|chromatogram| chromatogram.calculate_tic())
            .collect();

        // Create common retention time space
        let common_rt_space = create_common_rt_space(tics);
        // let common_rt_space = tics.clone();

        // Smooth and normalize TICs
        let smoothed_tics: Vec<_> = common_rt_space
            .iter()
            .map(|tic| {
                tic.smooth_sgolay(
                    self.parameters.alignment.smoothing.sgolay_window,
                    self.parameters.alignment.smoothing.sgolay_order,
                )?
                .normalize()
            })
            .collect::<Result<Vec<_>, _>>()?;

        /* ------------------------------------------------------------------ */
        /* Step 2. Pair-wise Alignment of TICs                                */
        /* ------------------------------------------------------------------ */

        log::debug!("Aligning TICs using {:?} using reference type: {:?}", self.parameters.alignment.method.as_str(), self.parameters.alignment.reference_type);
        let aligned_chromatograms = match self.parameters.alignment.method.as_str() {
            "dtw" => {
                match self.parameters.alignment.reference_type.as_str() {
                    "star" => star_align_tics(smoothed_tics.clone(), &self.parameters.alignment)?,
                    "mst" => mst_align_tics(smoothed_tics.clone())?,
                    "progressive" => progressive_align_tics(smoothed_tics.clone())?,
                    _ => star_align_tics(smoothed_tics.clone(), &self.parameters.alignment)?,
                }
            },
            "fft" => {
                match self.parameters.alignment.reference_type.as_str() {
                    "star" => star_align_tics_fft(smoothed_tics.clone(), &self.parameters.alignment)?,
                    "mst" => mst_align_tics_fft(smoothed_tics.clone())?,
                    "progressive" => progressive_align_tics_fft(smoothed_tics.clone())?,
                    _ => star_align_tics_fft(smoothed_tics.clone(), &self.parameters.alignment)?,
                }
            },
            "fft_dtw" => star_align_tics_fft_with_local_refinement(smoothed_tics.clone(), &self.parameters.alignment)?,
            _ => star_align_tics(smoothed_tics.clone(), &self.parameters.alignment)?,
        };
        Ok(aligned_chromatograms)
    }

    fn write_aligned_score_results_to_db(
        &self,
        feature_access: &Vec<OswAccess>,
        results: &Vec<Result<HashMap<i32, PrecursorAlignmentResult>, ArycalError>>, 
    ) -> Result<()> {
        // Ensure the FEATURE_ALIGNMENT table exists
        for osw_access in feature_access {
            osw_access.create_feature_alignment_table()?;
        }
    
        let mut batch = Vec::new();
        let batch_size = 1000; 
    
        let mut progress = None;
        if !self.parameters.disable_progress_bar {
            progress = Some(Progress::new(
                results.len(),
                "[arycal] Writing FEATURE_ALIGNMENT table to the database",
            ));
        }
    
        for result in results.iter() { // Iterate by reference
            match result {
                Ok(precursor_alignments) => {
                    for (_, run_alignments) in precursor_alignments {
                        for (_, peak_mappings) in &run_alignments.alignment_scores {
                            batch.push(peak_mappings); // Clone values to avoid moving ownership
                        }
                    }
    
                    // Insert in batches for better performance
                    if batch.len() >= batch_size {
                        for osw_access in feature_access {
                            log::debug!("Inserting batch of {} of {} full trace aligned features and scores", batch.len(), results.len());
                            osw_access.insert_feature_alignment_batch(&batch)?;
                        }
                        batch.clear();
                    }
                }
                Err(err) => {
                    log::warn!("[arycal::write_aligned_score_results_to_db] Skipping result due to error: {:?}", err);
                }
            }
            if !self.parameters.disable_progress_bar {
                progress.as_ref().expect("The Progess tqdm logger is not enabled").inc();
            }
        }
    
        // Insert any remaining records
        if !batch.is_empty() {
            for osw_access in feature_access {
                log::debug!("Inserting remaining batch of {} full trace aligned features and scores", batch.len());
                osw_access.insert_feature_alignment_batch(&batch)?;
            }
        }
    
        Ok(())
    }

    fn write_ms2_alignment_results_to_db(
        &self,
        feature_access: &Vec<OswAccess>,
        results: &Vec<Result<HashMap<i32, PrecursorAlignmentResult>, ArycalError>>, 
    ) -> Result<()> {
        // Ensure the FEATURE_MS2_ALIGNMENT table exists
        for osw_access in feature_access {
            osw_access.create_feature_ms2_alignment_table()?;
        }
    
        let mut batch = Vec::new();
        let batch_size = 1000; 
    
        let mut progress = None;
        if !self.parameters.disable_progress_bar {
            progress = Some(Progress::new(
                results.len(),
                "[arycal] Writing FEATURE_MS2_ALIGNMENT table to the database",
            ));
        }
    
        for result in results.iter() { // Iterate by reference
            match result {
                Ok(precursor_alignments) => {
                    for (_, run_alignments) in precursor_alignments {
                        for (_, peak_mappings) in &run_alignments.detecting_peak_mappings {
                            batch.extend(peak_mappings.iter().cloned()); // Clone values to avoid moving ownership
                        }
                    }
    
                    // Insert in batches for better performance
                    if batch.len() >= batch_size {
                        for osw_access in feature_access {
                            log::debug!("Inserting batch of {} of {} ms2 aligned features and scores", batch.len(), results.len());
                            osw_access.insert_feature_ms2_alignment_batch(&batch)?;
                        }
                        batch.clear();
                    }
                }
                Err(err) => {
                    log::warn!("[arycal::write_ms2_alignment_results_to_db] Skipping result due to error: {:?}", err);
                }
            }
            if !self.parameters.disable_progress_bar {
                progress.as_ref().expect("The Progess tqdm logger is not enabled").inc();
            }
        }
    
        // Insert any remaining records
        if !batch.is_empty() {
            for osw_access in feature_access {
                log::debug!("Inserting remaining batch of {} ms2 aligned features and scores", batch.len());
                osw_access.insert_feature_ms2_alignment_batch(&batch)?;
            }
        }
    
        Ok(())
    }
    
    

    fn write_transition_alignment_results_to_db(
        &self,
        feature_access: &Vec<OswAccess>,
        results: &Vec<Result<HashMap<i32, PrecursorAlignmentResult>, ArycalError>>, 
    ) -> Result<()> {
        // Ensure the FEATURE_TRANSITION_ALIGNMENT table exists
        for osw_access in feature_access {
            osw_access.create_feature_transition_alignment_table()?;
        }
    
        let mut batch = Vec::new();
        let batch_size = 1000; 
        
        let mut progress: Option<Progress> = None;
        if !self.parameters.disable_progress_bar {
            progress = Some(Progress::new(
                results.len(),
                "[arycal] Writing FEATURE_TRANSITION_ALIGNMENT table to the database",
            ));
        }   
    
        for result in results.iter() { // Iterate by reference
            match result {
                Ok(precursor_alignments) => {
                    for (_, run_alignments) in precursor_alignments {
                        for (_, peak_mappings) in &run_alignments.identifying_peak_mapping_scores {
                            batch.extend(peak_mappings.iter().cloned()); // Clone values to avoid moving ownership
                        }
                    }
    
                    // Insert in batches for better performance
                    if batch.len() >= batch_size {
                        for osw_access in feature_access {
                            log::debug!("Inserting batch of {} of {} transition aligned features and scores", batch.len(), results.len());
                            osw_access.insert_feature_transition_alignment_batch(&batch)?;
                        }
                        batch.clear();
                    }
                }
                Err(err) => {
                    log::warn!(
                        "[arycal::write_transition_alignment_results_to_db] Skipping transition alignment result due to error: {:?}",
                        err
                    );
                }
            }
            if !self.parameters.disable_progress_bar {
                progress.as_ref().expect("The Progess tqdm logger is not enabled").inc();
            }
        }
    
        // Insert any remaining records
        if !batch.is_empty() {
            for osw_access in feature_access {
                log::debug!("Inserting remaining batch of {} transition aligned features and scores", batch.len());
                osw_access.insert_feature_transition_alignment_batch(&batch)?;
            }
        }
    
        Ok(())
    }
    

}
