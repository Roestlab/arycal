use arycal_cli::input::Input;
use arycal_cloudpath::util::extract_basename;
use arycal_common::config::{OpenSwathConfig, RawFileType};
use chrono::Local;
use eframe::egui::{ScrollArea, Ui};
use egui::{Button, Color32};
use rfd::FileDialog;
use std::ffi::OsStr;
use std::fs::File;
use std::sync::atomic::AtomicUsize;
use std::time::{Duration, Instant};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::{Arc, Mutex, atomic::{AtomicBool, Ordering},
mpsc::{channel, Sender, Receiver}};
use std::thread::{self, JoinHandle};

/// Progress messages for batch runs
pub enum Progress {
    Started { index: usize, total: usize },
    Finished { index: usize, duration: Duration },
    Cancelled { index: usize },
    Failed { index: usize, code: i32 },
    AllDone,
}

/// Spawn a subprocess, capture its stdout/stderr in *parallel*, and return
/// both the I/O‐draining thread handle and a handle to the child.
pub fn spawn_process_runner(
    bin: PathBuf,
    args: Vec<String>,
    output_lines: Arc<Mutex<Vec<String>>>,
    envs: Option<Vec<(String, String)>>,
) -> (JoinHandle<()>, Arc<Mutex<Child>>) {
    let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
    let cmd_str = format!("{} {}", bin.display(), args.join(" "));
    output_lines.lock().unwrap()
        .push(format!("[{}] Running: {}", ts, cmd_str));

    let mut cmd = Command::new(&bin);

    // inject any requested env vars
    if let Some(envs) = envs {
        for (k, v) in envs {
            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
            output_lines.lock().unwrap()
                .push(format!("[{}] Setting env: {}={}", ts, k, v));
            cmd.env(k, v);
        }
    }

    // capture both stdout and stderr
    let mut child = cmd
        .args(&args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect(&format!("Failed to start process for: {}", cmd_str));

    let child_arc = Arc::new(Mutex::new(child));

    // take the pipes out
    let stdout = child_arc.lock().unwrap().stdout.take();
    let stderr = child_arc.lock().unwrap().stderr.take();

    let buf_out = Arc::clone(&output_lines);
    let buf_err = Arc::clone(&output_lines);

    // thread #1: drain stdout
    let t1 = thread::spawn(move || {
        if let Some(out) = stdout {
            for line in BufReader::new(out).lines().flatten() {
                buf_out.lock().unwrap().push(line);
            }
        }
    });

    // thread #2: drain stderr
    let t2 = thread::spawn(move || {
        if let Some(err) = stderr {
            for line in BufReader::new(err).lines().flatten() {
                buf_err.lock().unwrap().push(line);
            }
        }
    });

    // thread #3: join on both drains
    let handle = thread::spawn(move || {
        let _ = t1.join();
        let _ = t2.join();
    });

    (handle, child_arc)
}

/// Spawn up to `concurrency` worker threads to run `all_args` (each a Vec<String>)
/// in parallel, sending Progress back on `progress_tx`. Cancelling or failure
/// stops all.
///
/// - `current_children` collects all live Child handles so Stop can kill them.
/// - `cancel_flag` when set to true causes all workers to bail out as soon as possible.
/// - `envs` is any iterator over (key, value) pairs to set in each child's environment.
pub fn spawn_batch_runner(
    bin: PathBuf,
    all_args: Vec<Vec<String>>,
    output_lines: Arc<Mutex<Vec<String>>>,
    progress_tx: Sender<Progress>,
    cancel_flag: Arc<AtomicBool>,
    current_children: Arc<Mutex<Vec<Arc<Mutex<Child>>>>>,
    concurrency: usize,
    envs: Option<Vec<(String, String)>>,
) -> thread::JoinHandle<()> {
    // Make sure we have at least one worker.
    let concurrency = concurrency.max(1);
    let total = all_args.len();
    let args_arc = Arc::new(all_args);
    let next_index = Arc::new(AtomicUsize::new(0));

    // Wrap the Option<Vec<…>> in an Arc so we can cheaply clone it into each closure.
    let envs_arc = Arc::new(envs);

    // Spawn the worker threads.
    let mut workers = Vec::with_capacity(concurrency);
    for _ in 0..concurrency {
        let bin         = bin.clone();
        let args_arc    = Arc::clone(&args_arc);
        let cancel_flag = Arc::clone(&cancel_flag);
        let output_lines= Arc::clone(&output_lines);
        let current_children = Arc::clone(&current_children);
        let progress_tx = progress_tx.clone();
        let next_index  = Arc::clone(&next_index);
        // **clone the Arc** for this particular worker:
        let envs_for_worker = Arc::clone(&envs_arc);

        workers.push(thread::spawn(move || {
            while !cancel_flag.load(Ordering::Relaxed) {
                let idx = next_index.fetch_add(1, Ordering::Relaxed);
                if idx >= total {
                    break;
                }
                let one_based = idx + 1;
                let _ = progress_tx.send(Progress::Started { index: one_based, total });

                let start = Instant::now();
                let args = args_arc[idx].clone();

                // Here we _move_ out the Arc<Option<…>> by cloning the inner Option<Vec<…>>:
                let worker_envs: Option<Vec<(String, String)>> = (*envs_for_worker).clone();

                let (io_handle, child_arc) =
                    spawn_process_runner(bin.clone(), args, Arc::clone(&output_lines), worker_envs);

                // Register the child so Stop can kill it:
                {
                    let mut kids = current_children.lock().unwrap();
                    kids.push(child_arc.clone());
                }

                // Wait for I/O thread to finish
                let _ = io_handle.join();
                // Then wait on the child process:
                let exit_code = {
                    let mut guard = child_arc.lock().unwrap();
                    guard.wait().map(|s| s.code().unwrap_or(-1)).unwrap_or(-1)
                };

                // Unregister this child
                {
                    let mut kids = current_children.lock().unwrap();
                    if let Some(pos) = kids.iter().position(|c| Arc::ptr_eq(c, &child_arc)) {
                        kids.remove(pos);
                    }
                }

                // Handle cancellation / failure / success
                if cancel_flag.load(Ordering::Relaxed) {
                    let _ = progress_tx.send(Progress::Cancelled { index: one_based });
                    break;
                } else if exit_code != 0 {
                    let _ = progress_tx.send(Progress::Failed { index: one_based, code: exit_code });
                    // Stop all other workers:
                    cancel_flag.store(true, Ordering::Relaxed);
                    break;
                } else {
                    let dur = start.elapsed();
                    let _ = progress_tx.send(Progress::Finished { index: one_based, duration: dur });
                }
            }
        }));
    }

    // Finally, join all workers and fire AllDone.
    thread::spawn(move || {
        for w in workers {
            let _ = w.join();
        }
        let _ = progress_tx.send(Progress::AllDone);
    })
}


/// Given an OpenSwathConfig, build one CLI-arg vector per input file
fn build_arg_lists_from_config(cfg: &OpenSwathConfig) -> Vec<Vec<String>> {
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

        // spectral library (one-to-one or first)
        if !cfg.spectral_library_paths.is_empty() {
            let lib = &cfg.spectral_library_paths[idx.min(cfg.spectral_library_paths.len() - 1)];
            args.push("-tr".into());
            args.push(lib.to_string_lossy().into_owned());
        }
        // linear iRT
        if !cfg.linear_irt_library_paths.is_empty() {
            let lir =
                &cfg.linear_irt_library_paths[idx.min(cfg.linear_irt_library_paths.len() - 1)];
            args.push("-tr_irt".into());
            args.push(lir.to_string_lossy().into_owned());
        }
        // non-linear iRT
        if cfg.include_non_linear_irt && !cfg.nonlinear_irt_library_paths.is_empty() {
            let nlir = &cfg.nonlinear_irt_library_paths
                [idx.min(cfg.nonlinear_irt_library_paths.len() - 1)];
            args.push("-tr_irt_nonlinear".into());
            args.push(nlir.to_string_lossy().into_owned());
        }
        // output
        let run_basename = extract_basename(input.to_string_lossy().as_ref());
        let output_directory = cfg.output_path.clone();
        let output_filename = format!("{}.{}", run_basename, cfg.output_file_type.as_str());
        let output_file = output_directory.join(&output_filename).to_string_lossy().into_owned();
        args.push("-out_features".into());
        args.push(output_file);
        args.push("-out_chrom".into());
        args.push(output_directory.join(format!("{}.chrom.sqMass", run_basename))
            .to_string_lossy()
            .into_owned());

        if cfg.output_debug_files {
            args.push("-Debugging:irt_mzml".into());
            args.push(output_directory.join(format!("{}.irt.mzML", run_basename))
                .to_string_lossy()
                .into_owned());
            args.push("-Debugging:irt_trafo".into());
            args.push(output_directory.join(format!("{}.irt.trafoXML", run_basename))
                .to_string_lossy()
                .into_owned());
            args.push("-Calibration:debug_mz_file".into());
            args.push(output_directory.join(format!("{}.irt.mz.txt", run_basename))
                .to_string_lossy()
                .into_owned());
            if cfg.is_pasef {
                args.push("-Calibration:debug_im_file".into());
                args.push(output_directory.join(format!("{}.irt.im.txt", run_basename))
                    .to_string_lossy()
                    .into_owned());
            }
        }

        // flags
        if cfg.enable_ms1 {
            args.push("-enable_ms1".into());
        }
        if cfg.enable_ms1 {
            args.push("true".into());
        }
        if cfg.enable_ipf {
            args.push("-enable_ipf".into());
        }
        if cfg.enable_ipf {
            args.push("true".into());
        }
        if cfg.is_pasef {
            args.push("-pasef".into());
        }
        // windows
        args.push("-rt_extraction_window".into());
        args.push(cfg.rt_extraction_window.to_string());
        if cfg.is_pasef {
            args.push("-ion_mobility_window".into());
            args.push(cfg.ion_mobility_window.to_string());
        }
        args.push("-mz_extraction_window".into());
        args.push(cfg.ms2_mz_extraction_window.to_string());
        args.push("-mz_extraction_window_unit".into());
        args.push(cfg.ms2_mz_extraction_window_unit.clone());
        args.push("-mz_extraction_window_ms1".into());
        args.push(cfg.ms1_mz_extraction_window.to_string());
        args.push("-mz_extraction_window_ms1_unit".into());
        args.push(cfg.ms1_mz_extraction_window_unit.clone());
        if cfg.is_pasef {
            args.push("-im_extraction_window_ms1".into());
            args.push(cfg.ms1_ion_mobility_extraction_window.to_string());
        }
        args.push("-irt_mz_extraction_window".into());
        args.push(cfg.irt_mz_extraction_window.to_string());
        args.push("-irt_mz_extraction_window_unit".into());
        args.push(cfg.irt_mz_extraction_window_unit.clone());
        if cfg.is_pasef {
            args.push("-irt_im_extraction_window".into());
            args.push(cfg.irt_ion_mobility_extraction_window.to_string());
        }
        // mass correction
        args.push("-mz_correction_function".into());
        args.push(cfg.mz_correction_function.clone());
        // read options
        args.push("-readOptions".into());
        args.push(cfg.read_options.clone());
        if let Some(base_tmp_dir) = &cfg.temp_directory {
            args.push("-tempDirectory".into());
            // Clone the base temp dir and push the run‐specific subdirectory.
            let mut tmp = base_tmp_dir.clone();
            tmp.push(&run_basename);

            log::info!("Using temporary directory: {}", tmp.display());

            // Make sure it exists
            std::fs::create_dir_all(&tmp)
                .expect("Failed to create temp directory");

            // Finally tell PyProphet to use that
            args.push(tmp.to_string_lossy().into_owned());
        }
        // threads & batch
        args.push("-batchSize".into());
        args.push(cfg.batch_size.to_string());
        args.push("-threads".into());
        args.push(cfg.threads.to_string());
        args.push("-outer_loop_threads".into());
        args.push(cfg.outer_loop_threads.to_string());
        // RT normalization
        args.push("-RTNormalization:alignmentMethod".into());
        args.push(cfg.rt_normalization_alignment_method.clone());
        args.push("-RTNormalization:outlierMethod".into());
        args.push(cfg.rt_normalization_outlier_method.clone());
        if cfg.rt_normalization_estimate_best_peptides {
            args.push("-RTNormalization:estimateBestPeptides".into());
        }
        args.push("-RTNormalization:lowess:span".into());
        args.push(cfg.rt_normalization_lowess_span.to_string());
        // advanced params
        if !cfg.advanced_params.is_empty() {
            args.extend(
                cfg.advanced_params
                    .split_whitespace()
                    .map(|s| s.to_string()),
            );
        }
        all_args.push(args);
    }
    all_args
}

#[derive(PartialEq)]
enum OpenSwathTab {
    Console,
    QCPlots,
}

pub struct OpenSwathState {
    pub active_tab: OpenSwathTab,
    pub output_lines: Arc<Mutex<Vec<String>>>,
    pub n_concurrent_processes: usize, // Number of concurrent processes to run
    pub runner: Option<JoinHandle<()>>,
    pub progress_rx: Option<Receiver<Progress>>,
    pub cancel_flag: Arc<AtomicBool>,
    pub current_children: Arc<Mutex<Vec<Arc<Mutex<Child>>>>>,
    // progress bookkeeping
    pub total_number_of_runs: usize,
    /// Index of the run currently in progress (1‐based), or last started
    pub processing_run_n: usize,
    /// how many completed successfully
    pub successful_runs: usize,
    /// how many failed
    pub failed_runs: usize,
    /// processed fraction = (successful_runs + failed_runs) / total_number_of_runs
    pub processed_fraction: f32,
    pub processing_time: Vec<std::time::Duration>,
    pub avg_processing_time: std::time::Duration,
    pub initial_start_time: Instant,
    pub total_processing_time: std::time::Duration,
}

impl OpenSwathState {
    pub fn default() -> Self {
        OpenSwathState {
            active_tab: OpenSwathTab::Console,
            output_lines: Arc::new(Mutex::new(Vec::new())),
            n_concurrent_processes: 1, // Default to 1 concurrent process
            runner: None,
            progress_rx: None,
            cancel_flag: Arc::new(AtomicBool::new(false)),
            current_children: Arc::new(Mutex::new(Vec::new())),
            total_number_of_runs: 0,
            processing_run_n: 0,
            successful_runs: 0,
            failed_runs: 0,
            processed_fraction: 0.0,
            processing_time: Vec::new(),
            avg_processing_time: std::time::Duration::ZERO,
            initial_start_time: Instant::now(),
            total_processing_time: std::time::Duration::ZERO,
        }
    }

    pub fn ui(&mut self, ui: &mut Ui, config: &Input) {
        // update n_concurrent_processes from config
        self.n_concurrent_processes = config.n_concurrent_processes;
        // Warn if binary not found
        if let Some(cfg) = config.openswath.as_ref() {
            if cfg.binary_path.to_string_lossy().is_empty() || !cfg.binary_path.exists() {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 0.0;
                    ui.colored_label(Color32::RED, "Could not find OpenSwathWorkflow executable on system path. Please set the path in the config panel. If you do not have OpenSwathWorkflow installed, see ");
                    ui.hyperlink_to("OpenMS Read the Docs for installation instructions.", "https://openms.readthedocs.io/en/latest/about/installation.html");
                });
                ui.colored_label(Color32::RED, "If you have a conda environment, you can install OpenMS tools including OpenSwathWorkflow with:");
                ui.code("conda create -n openms python=3.10\nconda config --add channels defaults\nconda config --add channels bioconda\nconda config --add channels conda-forge\nconda install openms");
            } else {
                ui.horizontal(|ui| {
                    if ui
                        .selectable_label(self.active_tab == OpenSwathTab::Console, "Console Log")
                        .clicked()
                    {
                        self.active_tab = OpenSwathTab::Console;
                    }
                    if cfg.output_debug_files {
                        if ui
                            .selectable_label(self.active_tab == OpenSwathTab::QCPlots, "Quality Control Plots")
                            .clicked()
                        {
                            self.active_tab = OpenSwathTab::QCPlots;
                        }
                    }
                });
                ui.separator();
        
                match self.active_tab {
                    OpenSwathTab::Console => self.console_ui(ui, config),
                    OpenSwathTab::QCPlots => self.qc_plots_ui(ui, config),
                }
            }
        }

    }

    fn console_ui(&mut self, ui: &mut Ui, config: &Input) {
        let is_running = self.runner.is_some();
        if let Some(cfg) = &config.openswath {
            // Buttons
            ui.horizontal(|ui| {
                if ui.add_enabled(!is_running, Button::new("▶︎ Run Workflow")).clicked() {
                    let all_args = build_arg_lists_from_config(cfg);
                    self.total_number_of_runs = all_args.len();
                    self.processing_run_n = 0;
                    self.processing_time.clear();
                    self.avg_processing_time = Duration::ZERO;
                    self.total_processing_time = Duration::ZERO;
                    self.successful_runs = 0;
                    self.failed_runs = 0;
                    self.processed_fraction = 0.0;
                    self.cancel_flag.store(false, Ordering::Relaxed);

                    let (tx, rx) = channel();
                    self.initial_start_time = Instant::now();
                    let runner = spawn_batch_runner(
                        cfg.binary_path.clone(),
                        all_args,
                        Arc::clone(&self.output_lines),
                        tx,
                        Arc::clone(&self.cancel_flag),
                        Arc::clone(&self.current_children),
                        self.n_concurrent_processes,  
                        None
                    );
                    self.progress_rx = Some(rx);
                    self.runner = Some(runner);
                }

                if ui.add_enabled(!is_running, Button::new("❔ Show Help")).clicked() {
                    let bin = cfg.binary_path.clone();
                    let args = vec!["--helphelp".into()];
                    let out = Arc::clone(&self.output_lines);
                    let (_h, _c) = spawn_process_runner(bin, args, out, None);
                }

                if ui.add_enabled(is_running, Button::new("⏹ Stop")).clicked() {
                    self.cancel_flag.store(true, Ordering::Relaxed);
                    for child in self.current_children.lock().unwrap().drain(..) {
                        let _ = child.lock().unwrap().kill();
                    }
                }

                if ui.add_enabled(!is_running, Button::new("🗑 Clear")).clicked() {
                    self.output_lines.lock().unwrap().clear();
                }

                if ui.add_enabled(!is_running, Button::new("💾 Save")).clicked() {
                    let ts = Local::now().format("%Y%m%d_%H%M%S");
                    let name = format!("open_swath_log_{}.txt", ts);
                    if let Some(path) = FileDialog::new().set_title("Save Console Log")
                        .set_file_name(&name).save_file() {
                        if let Ok(mut f) = File::create(&path) {
                            for line in self.output_lines.lock().unwrap().iter() {
                                let _ = writeln!(f, "{}", line);
                            }
                        }
                    }
                }
            });

            ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
                for line in self.output_lines.lock().unwrap().iter() {
                    ui.label(line);
                }
            });

            // Process progress messages
            if let Some(rx) = self.progress_rx.take() {
                let mut keep = false;
                for msg in rx.try_iter() {
                    match msg {
                        Progress::Started { index, total } => {
                            self.processing_run_n = index;
                            self.total_number_of_runs = total;
                        }
                        Progress::Finished { duration, .. } => {
                            self.successful_runs += 1;
                            self.processing_time.push(duration);
                            let sum: std::time::Duration = self.processing_time.iter().sum();
                            // self.total_processing_time = sum;
                            self.total_processing_time = self.initial_start_time.elapsed();
                            self.avg_processing_time = sum / (self.processing_time.len() as u32);
                        }
                        Progress::Failed { index, code } => {
                            self.failed_runs += 1;
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(
                                format!("[{}] Run {}/{} failed with exit code {}", ts, index, self.total_number_of_runs, code)
                            );
                        }
                        Progress::Cancelled { index } => {
                            self.failed_runs += 1;
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(
                                format!("[{}] Run {}/{} cancelled by user", ts, index, self.total_number_of_runs)
                            );
                            self.runner = None;
                        }
                        Progress::AllDone => {
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(
                                format!("[{}] All runs complete", ts)
                            );
                            self.runner = None;
                        }
                    }
                    // update fraction after each event
                    let done = (self.successful_runs + self.failed_runs) as f32;
                    self.processed_fraction = if self.total_number_of_runs > 0 {
                        done / (self.total_number_of_runs as f32)
                    } else {
                        0.0
                    };
                }
                if self.runner.is_some() {
                    keep = true;
                }
                if keep {
                    self.progress_rx = Some(rx);
                }
            }

            if self.runner.is_some() {
                ui.ctx().request_repaint();
            }
        }
    }

    fn qc_plots_ui(&mut self, ui: &mut Ui, _config: &Input) {
        ui.label("Quality control plots visualization placeholder");
        // TODO: render any loaded results (e.g. plots, tables)
    }
}
