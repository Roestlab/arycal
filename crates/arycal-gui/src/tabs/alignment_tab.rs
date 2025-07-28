use arycal_cli::input::Input;
use arycal_cloudpath::util::extract_basename;
use chrono::Local;
use eframe::egui::{ScrollArea, Ui};
use egui::{Button, Color32};
use rfd::FileDialog;
use std::env;
use std::fs::File;
use std::sync::atomic::AtomicUsize;
use std::time::{Duration, Instant};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::{Arc, Mutex, atomic::{AtomicBool, Ordering},
mpsc::{channel, Sender, Receiver}};
use std::thread::{self, JoinHandle};

use super::open_swath_tab::{spawn_batch_runner, spawn_process_runner, Progress};


pub struct AlignmentState {
    /// Captured stdout/stderr lines
    pub output_lines: Arc<Mutex<Vec<String>>>,

    /// The handle returned by `spawn_batch_runner`
    pub io_handle: Option<JoinHandle<()>>,

    /// Receiver for progress events
    pub progress_rx: Option<Receiver<Progress>>,

    /// Cancellation flag
    pub cancel_flag: Arc<AtomicBool>,

    /// The one-and-only active child process (for kill)
    pub current_child: Arc<Mutex<Vec<Arc<Mutex<Child>>>>>,

    /// When we actually launched alignment
    pub start_time: Option<Instant>,

    /// How many files to process at once (for future use)
    pub batch_size: usize,
    /// How many input files in total
    pub n_runs: usize,

    /// How long the run took (filled after finish)
    pub total_processing_time: Duration,
}

impl AlignmentState {
    pub fn default() -> Self {
        AlignmentState {
            output_lines: Arc::new(Mutex::new(Vec::new())),
            io_handle: None,
            progress_rx: None,
            cancel_flag: Arc::new(AtomicBool::new(false)),
            current_child: Arc::new(Mutex::new(Vec::new())),
            start_time: None,
            batch_size: 1,
            n_runs: 0,
            total_processing_time: Duration::ZERO,
        }
    }

    pub fn ui(&mut self, ui: &mut Ui, config: &Input) {
        // refresh GUI-only fields
        self.batch_size = config.alignment.batch_size.unwrap_or(1);
        self.n_runs = config.xic.file_paths.len();

        // find the arycal binary
        let bin = if let Some(p) = &config.arycal_binary_path {
            p.clone()
        } else {
            ui.colored_label(Color32::RED, "No `arycal` binary configured.");
            return;
        };
        if !bin.exists() {
            ui.horizontal(|ui| {
                ui.colored_label(
                    Color32::RED,
                    "Could not find `arycal`; set it in Config → Alignment.",
                );
                ui.hyperlink_to("Releases →", "https://github.com/singjc/arycal/releases");
            });
            return;
        }

        let is_running = self.io_handle.is_some();

        ui.horizontal(|ui| {
            // ▶︎ Run Workflow
            if ui.add_enabled(!is_running, Button::new("▶︎ Run Workflow")).clicked() {
                // 1) dump Input→JSON to tmp file
                let ts = Local::now().format("%Y%m%d_%H%M%S");
                let mut tmp = env::temp_dir();
                tmp.push(format!("arycal_alignment_{}.json", ts));
                {
                    let mut f = File::create(&tmp).expect("create temp JSON");
                    serde_json::to_writer_pretty(&mut f, &config)
                        .expect("serialize JSON");
                    f.flush().ok();
                }

                // 2) reset console & cancel flag
                self.output_lines.lock().unwrap().clear();
                self.cancel_flag.store(false, Ordering::Relaxed);

                // 3) build the single‐element batch
                let args = vec![tmp.to_string_lossy().into_owned()];
                let batch = vec![args];
                let envs = Some(vec![
                    ("ARYCAL_LOG".to_string(), config.log_level.clone()),
                    ("RAYON_NUM_THREADS".to_string(), config.threads.to_string()),
                ]);

                // 4) spawn batch‐runner with concurrency = 1
                let (tx, rx) = channel();
                let handle = spawn_batch_runner(
                    bin.clone(),
                    batch,
                    Arc::clone(&self.output_lines),
                    tx,
                    Arc::clone(&self.cancel_flag),
                    Arc::clone(&self.current_child),
                    1, 
                    envs,
                );

                self.io_handle = Some(handle);
                self.progress_rx = Some(rx);
                self.start_time = None;
                self.total_processing_time = Duration::ZERO;
            }

            // ❔ Show Help
            if ui.add_enabled(!is_running, Button::new("❔ Show Help")).clicked() {
                let (_h, _c) = super::open_swath_tab::spawn_process_runner(
                    bin.clone(),
                    vec!["--help".into()],
                    Arc::clone(&self.output_lines),
                    None
                );
            }

            // ⏹ Stop
            if ui.add_enabled(is_running, Button::new("⏹ Stop")).clicked() {
                self.cancel_flag.store(true, Ordering::Relaxed);
                for child in self.current_child.lock().unwrap().drain(..) {
                    let _ = child.lock().unwrap().kill();
                }
            }

            // 🗑 Clear
            if ui.add_enabled(!is_running, Button::new("🗑 Clear")).clicked() {
                self.output_lines.lock().unwrap().clear();
            }

            // 💾 Save
            if ui.add_enabled(!is_running, Button::new("💾 Save")).clicked() {
                let ts = Local::now().format("%Y%m%d_%H%M%S");
                let name = format!("arycal_log_{}.txt", ts);
                if let Some(path) = FileDialog::new()
                    .set_title("Save Console Log")
                    .set_file_name(&name)
                    .save_file()
                {
                    let mut f = File::create(&path).unwrap();
                    for line in self.output_lines.lock().unwrap().iter() {
                        let _ = writeln!(f, "{}", line);
                    }
                }
            }
        });

        ui.separator();

        // console output pane
        ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
            for line in self.output_lines.lock().unwrap().iter() {
                ui.label(line);
            }
        });

        // handle Progress messages
        if let Some(rx) = &mut self.progress_rx {
            let mut saw_terminal = false;
    
            for msg in rx.try_iter() {
                match msg {
                    Progress::Started { .. } => {
                        self.start_time = Some(Instant::now());
                        self.output_lines
                            .lock()
                            .unwrap()
                            .push("[Started] alignment".into());
                    }
                    Progress::Finished { duration, .. } => {
                        self.total_processing_time = duration;
                        self.output_lines.lock().unwrap().push(
                            format!("[Finished] Completed in {:.2?}", duration)
                        );
                    }
                    Progress::Failed { code, .. } => {
                        self.output_lines.lock().unwrap().push(
                            format!("[Failed] exit code {}", code)
                        );
                        saw_terminal = true;
                    }
                    Progress::Cancelled { .. } => {
                        self.output_lines.lock().unwrap().push("[Cancelled] by user".into());
                        saw_terminal = true;
                    }
                    Progress::AllDone => {
                        self.output_lines
                            .lock()
                            .unwrap()
                            .push("[AllDone] alignment complete".into());
                        saw_terminal = true;
                    }
                }
            }
    
            if saw_terminal {
                // actually tear down **only** on terminal events
                self.progress_rx = None;
                self.io_handle   = None;
            }
        }

        // repaint while running
        if self.io_handle.is_some() {
            ui.ctx().request_repaint();
        }
    }
}