mod config;

use arycal_cli::input::Input;
use config::AppConfig;

use eframe::App;
use egui::{CentralPanel, ScrollArea, SidePanel, TextEdit, Ui, Window};
use rfd::FileDialog;
use serde_json;
use std::fs;
use std::sync::{Arc, Mutex};

pub struct ArycalApp {
    config: Input,
    config_file_path: String,
    show_log: bool,
    log_messages: Arc<Mutex<Vec<String>>>,
    sidebar_width: f32,
}

impl ArycalApp {
    pub fn new(config_file_path: String) -> Self {
        let config = Input::default();

        Self {
            config,
            config_file_path,
            show_log: true,
            log_messages: Arc::new(Mutex::new(Vec::new())),
            sidebar_width: 250.0,
        }
    }

    pub fn show(&mut self, ctx: &egui::Context) {
        SidePanel::left("config_panel")
            .default_width(self.sidebar_width)
            .resizable(true)
            .show(ctx, |ui| {
                self.edit_config(ui, ctx);
            });

        CentralPanel::default().show(ctx, |ui| {
            self.show_main_panel(ui, ctx);
        });
    }

    fn edit_config(&mut self, ui: &mut Ui, ctx: &egui::Context) {
        ui.heading("XIC Configuration");

        ui.checkbox(&mut self.config.xic.include_precursor, "Include Precursor");
        ui.add(
            egui::Slider::new(&mut self.config.xic.num_isotopes, 1..=10).text("Number of Isotopes"),
        );

        ui.label("Type of XIC files");
        // ui.text_edit_singleline(&mut self.config.xic.r#type);
        // Make a drop-down list for the type of XIC files (currently only sqMass type is supported)
        egui::containers::ComboBox::from_id_salt("xic_type")
            .selected_text(
                self.config
                    .xic
                    .file_type
                    .clone()
                    .unwrap_or_default()
                    .as_str(),
            )
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut self
                        .config
                        .xic
                        .file_type
                        .clone()
                        .unwrap_or_default()
                        .as_str(),
                    "sqMass",
                    "sqMass",
                );
            });

        Self::edit_file_paths(ui, &mut self.config.xic.file_paths, "XIC");

        ui.separator();

        ui.add_space(25.0);

        ui.heading("Features Configuration");

        ui.label("Type of feature files");
        // ui.text_edit_singleline(&mut self.config.features.r#type);
        egui::containers::ComboBox::from_id_salt("features_type")
            .selected_text(
                self.config
                    .features
                    .file_type
                    .clone()
                    .unwrap_or_default()
                    .as_str(),
            )
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut self
                        .config
                        .features
                        .file_type
                        .clone()
                        .unwrap_or_default()
                        .as_str(),
                    "osw",
                    "osw",
                );
            });

        Self::edit_file_paths(ui, &mut self.config.features.file_paths, "Feature");

        ui.add_space(25.0);

        ui.heading("Filters Configuration");
        ui.checkbox(&mut self.config.filters.decoy, "Exclude Decoys");

        ui.separator();

        ui.add_space(25.0);

        ui.heading("Alignment Configuration");

        ui.label("Batch Size");
        let mut batch_size_str = self
            .config
            .alignment
            .batch_size
            .unwrap_or_default()
            .to_string();
        ui.text_edit_singleline(&mut batch_size_str);
        if let Ok(batch_size) = batch_size_str.parse::<usize>() {
            self.config.alignment.batch_size = Some(batch_size);
        }

        ui.label("Alignment Method");
        // ui.text_edit_singleline(&mut self.config.alignment.method);
        // Make a drop-down list for the alignment method (currently DTW or FFT or FFT+DTW
        egui::containers::ComboBox::from_id_salt("alignment_method")
            .selected_text(self.config.alignment.method.clone().to_uppercase())
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut self.config.alignment.method, "DTW".to_string(), "DTW");
                ui.selectable_value(&mut self.config.alignment.method, "FFT".to_string(), "FFT");
                ui.selectable_value(
                    &mut self.config.alignment.method,
                    "FFT+DTW".to_string(),
                    "FFT+DTW",
                );
            });

        ui.label("Reference Type");
        // Make drop-down list for the reference type (currently star, MST, or progressive)
        egui::containers::ComboBox::from_id_salt("reference_type")
            .selected_text(self.config.alignment.reference_type.clone().to_uppercase())
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut self.config.alignment.reference_type,
                    "STAR".to_string(),
                    "star",
                );
                ui.selectable_value(
                    &mut self.config.alignment.reference_type,
                    "MST".to_string(),
                    "mst",
                );
                ui.selectable_value(
                    &mut self.config.alignment.reference_type,
                    "PROGRESSIVE".to_string(),
                    "progressive",
                );
            });

        ui.label("Reference File");
        // Two columns, one checkbox, and one selector for the basename of files in xic file paths
        ui.horizontal(|ui| {
            // Check if reference_run is not none/null
            let mut use_reference = self.config.alignment.reference_run.is_some();
            ui.checkbox(&mut use_reference, "Use Reference Run");

            // Enable or disable the reference file selection box based on the checkbox
            if use_reference {
                egui::containers::ComboBox::from_id_salt("reference_file")
                    .selected_text(
                        self.config
                            .alignment
                            .reference_run
                            .clone()
                            .unwrap_or_default(),
                    )
                    .show_ui(ui, |ui| {
                        for path in self.config.xic.file_paths.iter() {
                            if let Some(file_name) = path.file_name() {
                                if let Some(file_name) = file_name.to_str() {
                                    ui.selectable_value(
                                        &mut self.config.alignment.reference_run,
                                        Some(file_name.to_string()),
                                        file_name,
                                    );
                                }
                            }
                        }
                    });
            } else {
                // Only set reference_run to None if the checkbox is unchecked
                self.config.alignment.reference_run = None;
            }
        });

        ui.checkbox(&mut self.config.alignment.use_tic, "Use TIC");
        ui.add(
            egui::Slider::new(&mut self.config.alignment.smoothing.sgolay_window, 1..=21)
                .text("Savitzky-Golay Window"),
        );
        ui.add(
            egui::Slider::new(&mut self.config.alignment.smoothing.sgolay_order, 1..=5)
                .text("Savitzky-Golay Order"),
        );

        ui.add_space(25.0);

        ui.separator();

        ui.horizontal(|ui| {
            if ui.button("Save Configuration").clicked() {
                if let Err(e) = self.save_config() {
                    if let Ok(mut log) = self.log_messages.lock() {
                        log.push(format!("Failed to save config: {}", e));
                    }
                } else {
                    if let Ok(mut log) = self.log_messages.lock() {
                        log.push("Configuration saved successfully.".to_string());
                    }
                }
            }

            if ui.button("Load Configuration").clicked() {
                if let Err(e) = self.load_config() {
                    if let Ok(mut log) = self.log_messages.lock() {
                        log.push(format!("Failed to load config: {}", e));
                    }
                } else {
                    if let Ok(mut log) = self.log_messages.lock() {
                        log.push("Configuration loaded successfully.".to_string());
                    }
                }
            }
        });

        ui.separator();
    }

    fn edit_file_paths(
        ui: &mut Ui,
        file_paths: &mut Vec<std::path::PathBuf>,
        file_type_name: &str,
    ) {
        // Drag & Drop Zone
        ui.group(|ui| {
            ui.label(format!("Drag and drop {} files here:", file_type_name));

            let dropped_files = ui.input(|i| i.raw.dropped_files.clone());
            if !dropped_files.is_empty() {
                for file in dropped_files {
                    if let Some(path) = file.path {
                        file_paths.push(path);
                    }
                }
            }
        });

        // Display File Paths with Unique ScrollArea ID
        ScrollArea::vertical()
            .id_salt(file_type_name) // Use `file_type_name` to generate a unique ID
            .show(ui, |ui| {
                for path in file_paths.iter() {
                    ui.label(path.to_string_lossy());
                }
            });

        ui.separator();

        // Buttons Below
        ui.horizontal(|ui| {
            if ui
                .button(format!("➕ Add {} File", file_type_name))
                .clicked()
            {
                if let Some(selected_path) = FileDialog::new().pick_file() {
                    file_paths.push(selected_path);
                }
            }

            if !file_paths.is_empty()
                && ui
                    .button(format!("❌ Remove Last {}", file_type_name))
                    .clicked()
            {
                file_paths.pop();
            }
        });
    }

    fn save_config(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Open a file dialog to save the configuration file
        if let Some(file_path) = FileDialog::new().save_file() {
            self.config_file_path = file_path.to_str().unwrap().to_string();
        } else {
            return Ok(());
        }
        let config_json = serde_json::to_string_pretty(&self.config)?;
        fs::write(&self.config_file_path, config_json)?;
        Ok(())
    }

    fn load_config(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Open a file dialog to select the configuration file
        if let Some(file_path) = FileDialog::new().pick_file() {
            self.config_file_path = file_path.to_str().unwrap().to_string();
        } else {
            return Ok(());
        }
        self.config = Input::load(&self.config_file_path)?;
        Ok(())
    }

    fn show_main_panel(&mut self, ui: &mut Ui, ctx: &egui::Context) {
        // Add a thin top horizontal panel for navigation buttons, such as showing/hiding the log panel
        ui.horizontal(|ui| {
            if ui.button("▶ Run Alignment").clicked() {
                match self.run_alignment(ctx) {
                    Ok(_) => {
                        if let Ok(mut log) = self.log_messages.lock() {
                            log.push("Alignment started.".to_string());
                        }
                    }
                    Err(e) => {
                        if let Ok(mut log) = self.log_messages.lock() {
                            log.push(format!("Error starting alignment: {}", e));
                        }
                    }
                }
            }

            if ui.button("Show Log").clicked() {
                self.show_log = true;
            }

            if ui.button("Organize Windows").clicked() {
                ui.ctx().memory_mut(|mem| mem.reset_areas());
            }
        });
        ui.separator();

        ui.heading("Arycal GUI");
        ui.label("Have fun aligning chromatograms!");
        ui.separator();

        self.show_log_panel(ctx);
    }

    fn show_log_panel(&mut self, ctx: &egui::Context) {
        if self.show_log {
            Window::new("Log Panel")
                .pivot(egui::Align2::RIGHT_TOP)
                // .min_size(egui::Vec2::new(200.0, 200.0))
                .default_pos([(self.sidebar_width + 5.0), 100.0]) // Set default position
                .collapsible(true)
                .resizable(true)
                .open(&mut self.show_log)
                .show(ctx, |ui| {
                    ScrollArea::vertical().show(ui, |ui| {
                        let log_messages = self.log_messages.lock().unwrap();
                        for msg in log_messages.iter() {
                            ui.label(msg);
                        }
                    });
                });
        }
    }
}

use egui::Context;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

use arycal_cli::Runner;

use os_pipe::{PipeReader, PipeWriter, pipe};
use std::io;
use std::io::{BufRead, BufReader};
use std::os::unix::io::{AsRawFd, FromRawFd};
use std::process::{Command, Stdio};

impl ArycalApp {
    fn run_alignment(&mut self, ctx: &egui::Context) -> Result<(), Box<dyn std::error::Error>> {
        let log_messages = Arc::clone(&self.log_messages);
        let config = self.config.clone();

        // Wrap `ctx` in an `Arc`
        let ctx = Arc::new(ctx.clone());

        // Create pipes for capturing stdout and stderr
        let (reader, writer) = pipe()?;
        let reader_stderr = reader.try_clone()?; // Clone reader for stderr capture

        // Save the original stdout and stderr
        let original_stdout = io::stdout();
        let original_stderr = io::stderr();

        // Replace stdout and stderr with our pipe
        unsafe {
            libc::dup2(writer.as_raw_fd(), original_stdout.as_raw_fd());
            libc::dup2(writer.as_raw_fd(), original_stderr.as_raw_fd());
        }

        // Spawn a thread to read from the pipe and log messages
        let log_messages_clone = Arc::clone(&log_messages);
        let ctx_clone = Arc::clone(&ctx);
        thread::spawn(move || {
            let reader = BufReader::new(reader);
            for line in reader.lines() {
                if let Ok(line) = line {
                    let mut log = log_messages_clone.lock().unwrap();
                    log.push(line);
                }
                ctx_clone.request_repaint(); // Force UI to update
            }
        });

        // Run the alignment in another thread
        let log_messages_clone = Arc::clone(&log_messages);
        thread::spawn(move || {
            let runner = Runner::new(config).unwrap();
            match runner.run() {
                Ok(_) => {
                    log_messages_clone
                        .lock()
                        .unwrap()
                        .push("Alignment completed successfully.".to_string());
                }
                Err(e) => {
                    log_messages_clone
                        .lock()
                        .unwrap()
                        .push(format!("Error running alignment: {}", e));
                }
            }

            // Restore original stdout and stderr
            unsafe {
                libc::dup2(original_stdout.as_raw_fd(), io::stdout().as_raw_fd());
                libc::dup2(original_stderr.as_raw_fd(), io::stderr().as_raw_fd());
            }
        });

        Ok(())
    }
}

impl App for ArycalApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.show(ctx);
    }
}
