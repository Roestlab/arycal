use arycal_cli::input::Input;
use arycal_cloudpath::util::{add_suffix_and_ext, extract_basename, extract_directory};
use arycal_common::config::PQPConfig;
use chrono::{format, Local};
use eframe::egui::Ui;
use egui::{Button, Color32, ScrollArea, TextBuffer};
use rfd::FileDialog;
use std::{fs::File, path::PathBuf, sync::{
    atomic::{AtomicBool, Ordering}, mpsc::{channel, Receiver, Sender}, Arc, Mutex
}, time::{Duration, Instant}};
use std::thread::{self, JoinHandle};
use std::process::Child;
use std::io::{BufRead, BufReader, Write};

use super::open_swath_tab::{spawn_batch_runner, spawn_process_runner, Progress};

// ------------------ Arg list builders ------------------


fn tfc_build_arg_lists_from_config(cfg: &PQPConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();
    let n_inputs = cfg.file_paths.len().max(1);
    for idx in 0..n_inputs {
        // pick input path
        let input = if !cfg.file_paths.is_empty() {
            &cfg.file_paths[idx.min(cfg.file_paths.len() - 1)]
        } else {
            continue; // nothing to run
        };
        let mut args = Vec::new();
        // input file
        args.push("-in".into());
        args.push(input.to_string_lossy().into_owned());


        // output
        let output_file = if n_inputs == 1 {
            cfg.output_path.clone().to_string_lossy().into_owned()
        } else {
            let run_basename = extract_basename(input.to_string_lossy().as_ref());
            let base_directory = cfg.output_path.clone().to_string_lossy().into_owned();
            // Add suffix and extension to the output file
            add_suffix_and_ext(&base_directory, &run_basename, "", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
        };
        
        args.push("-out".into());
        args.push(output_file);

        // flags
        if cfg.tfc_legacy_traml_id {
            args.push("-tfc_legacy_traml_id".into());
        }
        if cfg.tfc_legacy_traml_id {
            args.push("true".into());
        }

        // threads 
        args.push("-threads".into());
        args.push(cfg.tfc_threads.to_string());

        // advanced params
        if !cfg.tfc_advanced.is_empty() {
            args.extend(
                cfg.tfc_advanced
                    .split_whitespace()
                    .map(|s| s.to_string()),
            );
        }
        all_args.push(args);
    }
    all_args
}

fn osg_build_arg_lists_from_config(cfg: &PQPConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();
    let n_inputs = cfg.file_paths.len().max(1);
    for idx in 0..n_inputs {
        // pick input path
        let input = if !cfg.file_paths.is_empty() {
            &cfg.file_paths[idx.min(cfg.file_paths.len() - 1)]
        } else {
            continue; // nothing to run
        };
        let mut args = Vec::new();
        // input file
        args.push("-in".into());
        args.push(input.to_string_lossy().into_owned());

        // output
        let output_file = if n_inputs == 1 {
            if cfg.odg_enabled & cfg.main_mode {
                let run_basename = extract_basename(cfg.output_path.clone().to_string_lossy().as_str());
                let base_directory = extract_directory(cfg.output_path.clone().to_string_lossy().as_str())
                    .unwrap_or_else(|| PathBuf::from("."));
                let targets_path = add_suffix_and_ext(&base_directory, &run_basename, "_targets", cfg.pqp_out_type.as_str());
                // cfg.odg_input = targets_path.clone();
                targets_path.to_string_lossy().into_owned()
            } else {
                cfg.output_path.clone().to_string_lossy().into_owned()
            }
        } else {
            let run_basename = extract_basename(input.to_string_lossy().as_ref());
            // let base_directory = extract_directory(input.to_string_lossy().as_ref())
                // .unwrap_or_else(|| PathBuf::from("."));
            let base_directory = cfg.output_path.clone().to_string_lossy().into_owned();
            let output_path = if cfg.main_mode {
                add_suffix_and_ext(&base_directory, &run_basename, "_targets", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
            } else {
                add_suffix_and_ext(&base_directory, &run_basename, "_nonlinear_irt", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
            };
            output_path
        };
        
        args.push("-out".into());
        args.push(output_file);

        // flags
        args.push("-min_transitions".into());
        args.push(cfg.osg_min_transitions.to_string());
        args.push("-max_transitions".into());
        args.push(cfg.osg_max_transitions.to_string());
        args.push("-allowed_fragment_types".into());
        args.push(cfg.osg_allowed_fragment_types.clone());
        args.push("-allowed_fragment_charges".into());
        args.push(cfg.osg_allowed_fragment_charges.clone());

        if cfg.main_mode {
            if cfg.osg_enable_detection_specific_losses {
                args.push("-enable_detection_specific_losses".into());
                args.push("true".into());
            }
    
            if cfg.osg_enable_detection_unspecific_losses {
                args.push("-enable_detection_unspecific_losses".into());
                args.push("true".into());
            }
    
            if cfg.osg_enable_ipf {
                args.push("-enable_ipf".into());
                args.push("true".into());
            }
    
            // Check if osg_unimod_file is not an empty PathBuf
            if !cfg.osg_unimod_file.as_os_str().is_empty() {
                args.push("-unimod_file".into());
                args.push(cfg.osg_unimod_file.to_string_lossy().into_owned());
            }
        }
        
        // advanced params
        if !cfg.osg_advanced.is_empty() {
            args.extend(
                cfg.osg_advanced
                    .split_whitespace()
                    .map(|s| s.to_string()),
            );
        }
        all_args.push(args);
    }
    all_args
}

fn odg_build_arg_lists_from_config(cfg: &PQPConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();
    let n_inputs = cfg.file_paths.len().max(1);
    for idx in 0..n_inputs {
        // pick input path
        let input = if !cfg.file_paths.is_empty() {
            &cfg.file_paths[idx.min(cfg.file_paths.len() - 1)]
        } else {
            continue; // nothing to run
        };
        let mut args = Vec::new();
        // input file

        // If OSG is enabled, use the output from OSG as input for ODG
        let input_file: String = if cfg.osg_enabled && cfg.main_mode {
            let inner_input = if n_inputs == 1 {
                let run_basename = extract_basename(cfg.output_path.clone().to_string_lossy().into_owned());
                // Extract directory of the output path
                let base_directory = extract_directory(cfg.output_path.clone().to_string_lossy().as_str())
                    .unwrap_or_else(|| PathBuf::from("."));
                add_suffix_and_ext(&base_directory, &run_basename, "_targets", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
            } else {
                let run_basename = extract_basename(input.to_string_lossy().as_ref());
                // let base_directory = extract_directory(input.to_string_lossy().as_ref())
                    // .unwrap_or_else(|| PathBuf::from("."));
                let base_directory = cfg.output_path.clone().to_string_lossy().into_owned();
                add_suffix_and_ext(&base_directory, &run_basename, "_targets", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
            };
            inner_input
        } else {
            input.to_string_lossy().into_owned()
        };
            
        args.push("-in".into());
        args.push(input_file);

        // output
        let output_file = if n_inputs == 1 {
            cfg.output_path.clone().to_string_lossy().into_owned()
        } else {
            let run_basename = extract_basename(input.to_string_lossy().as_ref());
            let base_directory = cfg.output_path.clone().to_string_lossy().into_owned();
            // Add suffix and extension to the output file
            add_suffix_and_ext(&base_directory, &run_basename, "_decoys", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
        };
        
        args.push("-out".into());
        args.push(output_file);

        // flags
        args.push("-method".into());
        args.push(cfg.odg_method.clone());
        args.push("-decoy_tag".into());
        args.push(cfg.odg_decoy_tag.clone());
        args.push("-min_decoy_fraction".into());
        args.push(cfg.odg_min_decoy_fraction.to_string());
        args.push("-aim_decoy_fraction".into());
        args.push(cfg.odg_aim_decoy_fraction.to_string());
        args.push("-shuffle_max_attempts".into());
        args.push(cfg.odg_shuffle_max_attempts.to_string());
        args.push("-shuffle_sequence_identity_threshold".into());
        args.push(cfg.odg_shuffle_sequence_identity_threshold.to_string());
        args.push("-shift_precursor_mz_shift".into());
        args.push(cfg.odg_shift_precursor_mz_shift.to_string());
        args.push("-shift_product_mz_shift".into());
        args.push(cfg.odg_shift_product_mz_shift.to_string());
        
        // advanced params
        if !cfg.odg_advanced.is_empty() {
            args.extend(
                cfg.odg_advanced
                    .split_whitespace()
                    .map(|s| s.to_string()),
            );
        }
        
        all_args.push(args);
    }
    all_args
}

fn easypqp_reduce_build_arg_lists_from_config(cfg: &PQPConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();
    let n_inputs = cfg.file_paths.len().max(1);
    for idx in 0..n_inputs {
        // pick input path
        let input = if !cfg.file_paths.is_empty() {
            &cfg.file_paths[idx.min(cfg.file_paths.len() - 1)]
        } else {
            continue; // nothing to run
        };

        // If OSG is enabled, use the output from OSG as input for EasyPQP reduce
        let input_file = if cfg.osg_enabled && !cfg.main_mode {
            let inner_input = if n_inputs == 1 {
                let run_basename = extract_basename(cfg.output_path.clone().to_string_lossy().into_owned());
                // Extract directory of the output path
                let base_directory = extract_directory(cfg.output_path.clone().to_string_lossy().as_str())
                    .unwrap_or_else(|| PathBuf::from("."));
                add_suffix_and_ext(&base_directory, &run_basename, "_nonlinear_irt", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
            } else {
                let run_basename = extract_basename(input.to_string_lossy().as_ref());
                // let base_directory = extract_directory(input.to_string_lossy().as_ref())
                //     .unwrap_or_else(|| PathBuf::from("."));
                let base_directory = cfg.output_path.clone().to_string_lossy().into_owned();
                add_suffix_and_ext(&base_directory, &run_basename, "_nonlinear_irt", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
            };
            inner_input
        } else {
            input.to_string_lossy().into_owned()
        };


        let mut args = Vec::new();

        // Add reduce subcommand as the first argument
        args.push("reduce".into());

        // input file
        args.push("--in".into());
        args.push(input_file);

        // output
        let output_file = if n_inputs == 1 {
            cfg.output_path.clone().to_string_lossy().into_owned()
        } else {
            let run_basename = extract_basename(input.to_string_lossy().as_ref());
            // let base_directory = extract_directory(input.to_string_lossy().as_ref())
                // .unwrap_or_else(|| PathBuf::from("."));
            let base_directory = cfg.output_path.clone().to_string_lossy().into_owned();
            // Add suffix and extension to the output file
            add_suffix_and_ext(&base_directory, &run_basename, "_linear_irt", cfg.pqp_out_type.as_str()).to_string_lossy().into_owned()
        };
        
        args.push("--out".into());
        args.push(output_file);

        // flags
        args.push("--bins".into());
        args.push(cfg.irt_bins.to_string());
        args.push("--peptides".into());
        args.push(cfg.irt_num_peptides.to_string());
        
        all_args.push(args);
    }
    all_args
}

// ------------------ Pipeline definitions ------------------

#[derive(Clone, Debug)]
enum PQPTool {
    TFC,
    OSG,
    ODG,
    EasyPQPReduce,
}

impl PQPTool {
    fn as_str(&self) -> &'static str {
        match self {
            PQPTool::TFC => "TargetedFileConverter",
            PQPTool::OSG => "OpenSwathAssayGenerator",
            PQPTool::ODG => "OpenSwathDecoyGenerator",
            PQPTool::EasyPQPReduce => "EasyPQP reduce",
        }
    }
}

#[derive(Clone)]
pub struct Stage {
    tool: PQPTool,
    bin: PathBuf,
    arg_lists: Vec<Vec<String>>,
    pub display_name: &'static str,
}

pub struct StageRunner {
    handle: JoinHandle<()>,
    progress_rx: Receiver<Progress>,
    pub stage: Stage,
}

pub struct PQPState {
    pub output_lines: Arc<Mutex<Vec<String>>>,
    pub cancel_flag: Arc<AtomicBool>,
    pub n_concurrent_processes: usize, // Number of concurrent processes to run
    pub current_children: Arc<Mutex<Vec<Arc<Mutex<Child>>>>>,

    // Current stage runner (replaces previous runner/progress_rx)
    pub current_stage_runner: Option<StageRunner>,

    // Sequential pipeline
    pub pipeline: Vec<Stage>,
    pub current_stage_index: usize,

    // Progress bookkeeping (per stage + global)
    pub total_number_of_runs: usize,          // runs in current stage
    pub processing_run_n: usize,              // current run index inside stage
    pub successful_runs: usize,               // global successes
    pub failed_runs: usize,                   // global failures
    pub processed_fraction: f32,              // per-stage fraction
    pub processing_time: Vec<Duration>,       // per-stage durations
    pub avg_processing_time: Duration,
    pub initial_start_time: Instant,
    pub total_processing_time: Duration,

    // Per-stage counts
    pub stage_successful_runs: usize,
    pub stage_failed_runs: usize,

    // Global total (sum of all stage arg counts)
    pub overall_total_runs: usize,
}

impl PQPState {
    pub fn default() -> Self {
        Self {
            output_lines: Arc::new(Mutex::new(Vec::new())),
            cancel_flag: Arc::new(AtomicBool::new(false)),
            n_concurrent_processes: 1, // Default to 1 concurrent process
            current_children: Arc::new(Mutex::new(Vec::new())),
            current_stage_runner: None,
            pipeline: Vec::new(),
            current_stage_index: 0,
            total_number_of_runs: 0,
            processing_run_n: 0,
            successful_runs: 0,
            failed_runs: 0,
            processed_fraction: 0.0,
            processing_time: Vec::new(),
            avg_processing_time: Duration::ZERO,
            initial_start_time: Instant::now(),
            total_processing_time: Duration::ZERO,
            stage_successful_runs: 0,
            stage_failed_runs: 0,
            overall_total_runs: 0,
        }
    }

    pub fn ui(&mut self, ui: &mut Ui, config: &Input) {
        // update n_concurrent_processes from config
        self.n_concurrent_processes = config.n_concurrent_processes;
        if let Some(cfg) = config.pqp.as_ref() {
            let mut missing = Vec::new();
            if cfg.tfc_enabled && (cfg.tfc_binary_path.as_os_str().is_empty() || !cfg.tfc_binary_path.exists()) {
                missing.push("TargetedFileConverter");
            }
            if cfg.osg_enabled && (cfg.osg_binary_path.as_os_str().is_empty() || !cfg.osg_binary_path.exists()) {
                missing.push("OpenSwathAssayGenerator");
            }
            if cfg.odg_enabled && (cfg.odg_binary_path.as_os_str().is_empty() || !cfg.odg_binary_path.exists()) {
                missing.push("OpenSwathDecoyGenerator");
            }
            if cfg.irt_reduce_enabled && (cfg.easypqp_binary_path.as_os_str().is_empty() || !cfg.easypqp_binary_path.exists()) {
                missing.push("EasyPQP");
            }

            if !missing.is_empty() {
                let names = missing.join(" / ");
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 0.0;
                    ui.colored_label(
                        Color32::RED,
                        format!(
                            "Could not detect enabled tool(s): {}. Please set their binary paths in the config panel. See ",
                            names
                        )
                    );
                    ui.hyperlink_to(
                        "OpenMS installation instructions",
                        "https://openms.readthedocs.io/en/latest/about/installation.html",
                    );
                });
                ui.colored_label(
                    Color32::RED,
                    "Or install via conda: \
                     `conda create -n openms python=3.10 && \
                     conda config --add channels defaults && \
                     conda config --add channels bioconda && \
                     conda config --add channels conda-forge && \
                     conda install openms`",
                );
                if missing.contains(&"EasyPQP") {
                    ui.colored_label(
                        Color32::RED,
                        "EasyPQP missing: install with `pip install easypqp`"
                    );
                }
            } else {
                self.console_ui(ui, config);
            }
        }
    }

    pub fn is_running(&self) -> bool {
        self.current_stage_runner.is_some()
    }

    fn build_pipeline(cfg: &PQPConfig) -> Vec<Stage> {
        let mut stages = Vec::new();

        if cfg.tfc_enabled {
            stages.push(Stage {
                tool: PQPTool::TFC,
                bin: cfg.tfc_binary_path.clone(),
                arg_lists: tfc_build_arg_lists_from_config(cfg),
                display_name: "TargetedFileConverter",
            });
        }
        if cfg.osg_enabled {
            stages.push(Stage {
                tool: PQPTool::OSG,
                bin: cfg.osg_binary_path.clone(),
                arg_lists: osg_build_arg_lists_from_config(cfg),
                display_name: "OpenSwathAssayGenerator",
            });
        }
        if cfg.main_mode {
            if cfg.odg_enabled {
                stages.push(Stage {
                    tool: PQPTool::ODG,
                    bin: cfg.odg_binary_path.clone(),
                    arg_lists: odg_build_arg_lists_from_config(cfg),
                    display_name: "OpenSwathDecoyGenerator",
                });
            }
        } else {
            if cfg.irt_reduce_enabled {
                // Append 'reduce' subcommand at spawn time by *embedding* in bin path? (Simplify: modify bin PathBuf)
                // We keep bin as original path; we will handle subcommand in spawn if needed OR you can
                // store a modified PathBuf here:
                // let reduce_bin = PathBuf::from(format!("{} reduce", cfg.easypqp_binary_path.display()));
                stages.push(Stage {
                    tool: PQPTool::EasyPQPReduce,
                    bin: cfg.easypqp_binary_path.clone(),
                    arg_lists: easypqp_reduce_build_arg_lists_from_config(cfg),
                    display_name: "EasyPQP reduce",
                });
            }
        }

        stages
    }

    fn start_next_stage(&mut self) {
        if self.current_stage_index >= self.pipeline.len() {
            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
            self.output_lines.lock().unwrap().push(
                format!("[{}] Workflow complete ({} stages).", ts, self.pipeline.len())
            );
            self.current_stage_runner = None;
            return;
        }

        let stage = self.pipeline[self.current_stage_index].clone();

        // Per-stage reset
        self.total_number_of_runs = stage.arg_lists.len();
        self.processing_run_n = 0;
        self.processing_time.clear();
        self.avg_processing_time = Duration::ZERO;
        self.total_processing_time = Duration::ZERO;
        self.stage_successful_runs = 0;
        self.stage_failed_runs = 0;
        self.cancel_flag.store(false, Ordering::Relaxed);
        self.processed_fraction = 0.0;

        let (tx, rx) = channel();
        let handle = spawn_batch_runner(
            stage.bin.clone(),
            stage.arg_lists.clone(),
            Arc::clone(&self.output_lines),
            tx,
            Arc::clone(&self.cancel_flag),
            Arc::clone(&self.current_children),
            self.n_concurrent_processes,
            None  
        );

        let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
        self.output_lines.lock().unwrap().push(
            format!(
                "[{}] Starting stage {} ({}/{}) with {} run(s)",
                ts,
                stage.display_name,
                self.current_stage_index + 1,
                self.pipeline.len(),
                self.total_number_of_runs
            )
        );

        self.current_stage_runner = Some(StageRunner {
            handle,
            progress_rx: rx,
            stage,
        });
    }

    fn abort_pipeline(&mut self, user: bool) {
        self.cancel_flag.store(true, Ordering::Relaxed);
        // kill every running child
        for child in self.current_children.lock().unwrap().drain(..) {
            let _ = child.lock().unwrap().kill();
        }
        let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
        self.output_lines.lock().unwrap().push(
            if user {
                format!("[{}] Pipeline cancelled by user.", ts)
            } else {
                format!("[{}] Pipeline aborted.", ts)
            }
        );
        self.current_stage_runner = None;
        self.pipeline.clear();
    }

    pub fn overall_progress_fraction(&self) -> f32 {
        if self.overall_total_runs == 0 { return 0.0; }
        let done = (self.successful_runs + self.failed_runs) as f32;
        done / (self.overall_total_runs as f32)
    }

    fn console_ui(&mut self, ui: &mut Ui, config: &Input) {
        let is_running = self.is_running();

        if let Some(pqp_cfg) = config.pqp.as_ref() {
            ui.horizontal(|ui| {

                // Run Workflow
                if ui.add_enabled(!is_running, Button::new("▶︎ Run Workflow")).clicked() {
                    self.pipeline = Self::build_pipeline(pqp_cfg);
                    self.current_stage_index = 0;
                    self.successful_runs = 0;
                    self.failed_runs = 0;
                    self.overall_total_runs = self.pipeline.iter().map(|s| s.arg_lists.len()).sum();
                    self.initial_start_time = Instant::now();
                    if self.pipeline.is_empty() {
                        self.output_lines.lock().unwrap()
                            .push("No enabled tools to run.".into());
                    } else {
                        self.start_next_stage();
                    }
                }

                // Show Help (spawns help for each enabled — unchanged basic logic, sequential output)
                if ui.add_enabled(!is_running, Button::new("❔ Show Help")).clicked() {
                    if pqp_cfg.tfc_enabled {
                        let (_h,_c)=spawn_process_runner(pqp_cfg.tfc_binary_path.clone(), vec!["--helphelp".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }
                    if pqp_cfg.osg_enabled {
                        let (_h,_c)=spawn_process_runner(pqp_cfg.osg_binary_path.clone(), vec!["--helphelp".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }
                    if pqp_cfg.odg_enabled {
                        let (_h,_c)=spawn_process_runner(pqp_cfg.odg_binary_path.clone(), vec!["--helphelp".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }
                    if pqp_cfg.irt_reduce_enabled {
                        let (_h,_c)=spawn_process_runner(pqp_cfg.easypqp_binary_path.clone(), vec!["reduce".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }
                }

                // Stop
                if ui.add_enabled(is_running, Button::new("⏹ Stop")).clicked() {
                    self.abort_pipeline(true);
                }

                // Clear Console
                if ui.add_enabled(!is_running, Button::new("🗑 Clear")).clicked() {
                    self.output_lines.lock().unwrap().clear();
                }

                // Save Console Log
                if ui.add_enabled(!is_running, Button::new("💾 Save")).clicked() {
                    let ts = Local::now().format("%Y%m%d_%H%M%S");
                    let default_name = format!("peptide_query_parameter_generation_log_{}.txt", ts);
                    if let Some(path) = FileDialog::new()
                        .set_title("Save Console Log")
                        .set_file_name(&default_name)
                        .save_file()
                    {
                        if let Ok(mut f) = File::create(&path) {
                            for line in self.output_lines.lock().unwrap().iter() {
                                let _ = writeln!(f, "{}", line);
                            }
                        }
                    }
                }
            });

            // Status / progress bars
            if !self.pipeline.is_empty() {
                if let Some(runner) = &self.current_stage_runner {
                    ui.label(format!(
                        "Stage {}/{}: {}",
                        self.current_stage_index + 1,
                        self.pipeline.len(),
                        runner.stage.display_name
                    ));
                    ui.label(format!(
                        "Stage progress: {}/{} (success {}, failed {})",
                        self.stage_successful_runs + self.stage_failed_runs,
                        self.total_number_of_runs,
                        self.stage_successful_runs,
                        self.stage_failed_runs
                    ));
                }
                // Overall
                ui.label(format!(
                    "Overall progress: {}/{} (success {}, failed {})",
                    self.successful_runs + self.failed_runs,
                    self.overall_total_runs,
                    self.successful_runs,
                    self.failed_runs
                ));
            }

            ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
                for line in self.output_lines.lock().unwrap().iter() {
                    ui.label(line);
                }
            });

            // Process progress messages for current stage
            let mut need_repaint = false;
            if let Some(stage_runner) = self.current_stage_runner.as_ref() {
                let rx = &stage_runner.progress_rx;
                for msg in rx.try_iter() {
                    match msg {
                        Progress::Started { index, total } => {
                            self.processing_run_n = index;
                            self.total_number_of_runs = total;
                        }
                        Progress::Finished { duration, .. } => {
                            self.successful_runs += 1;
                            self.stage_successful_runs += 1;
                            self.processing_time.push(duration);
                            let sum: Duration = self.processing_time.iter().sum();
                            // self.total_processing_time = sum;
                            self.total_processing_time = self.initial_start_time.elapsed();
                            self.avg_processing_time = sum / (self.processing_time.len() as u32);
                        }
                        Progress::Failed { index, code } => {
                            self.failed_runs += 1;
                            self.stage_failed_runs += 1;
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(
                                format!(
                                    "[{}] Stage {} run {}/{} failed (exit {})",
                                    ts,
                                    stage_runner.stage.display_name,
                                    index,
                                    self.total_number_of_runs,
                                    code
                                )
                            );
                        }
                        Progress::Cancelled { index } => {
                            self.failed_runs += 1;
                            self.stage_failed_runs += 1;
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(
                                format!(
                                    "[{}] Stage {} run {}/{} cancelled",
                                    ts,
                                    stage_runner.stage.display_name,
                                    index,
                                    self.total_number_of_runs
                                )
                            );
                            // Cancel whole pipeline
                            self.abort_pipeline(false);
                            break;
                        }
                        Progress::AllDone => {
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(
                                format!(
                                    "[{}] Stage {} complete: {} success, {} failed.",
                                    ts,
                                    stage_runner.stage.display_name,
                                    self.stage_successful_runs,
                                    self.stage_failed_runs
                                )
                            );
                            // Advance to next
                            self.current_stage_index += 1;
                            self.current_stage_runner = None;
                            self.start_next_stage();
                            break;
                        }
                    }
                    // per-stage fraction update
                    let stage_done = (self.stage_successful_runs + self.stage_failed_runs) as f32;
                    self.processed_fraction = if self.total_number_of_runs > 0 {
                        stage_done / (self.total_number_of_runs as f32)
                    } else { 0.0 };
                }
                need_repaint = self.is_running();
            }

            if need_repaint {
                ui.ctx().request_repaint();
            }
        }
    }
}