use crate::panels::{config_panel, system_footer, visualization_settings};
use crate::tabs::{
    alignment_tab::AlignmentState,
    open_swath_tab::OpenSwathState,
    spectral_library_tab::PQPState,
    visualization_tab::VisualizationState,
};
use crate::{
    config::AppConfig, panels::peptide_query_parameter_settings::draw_pqp_generation, panels::openswath_settings::draw_open_swath,
    panels::alignment_settings::draw_alignment,
    panels::validation_settings::draw_validation,
    tabs::validation_tab::ValidationState,
};
use arycal_cli::{input::Input, Runner};
use eframe::{
    egui::{CentralPanel, SidePanel, Ui},
    App,
};
use egui::{Align, Layout, TopBottomPanel};
use egui_notify::Toasts;
use log::{error, info};
use rfd::FileDialog;
use std::fmt::format;
use std::time::{Duration, Instant};
use std::{
    fs,
    sync::{Arc, Mutex},
    thread,
};
#[cfg(not(target_arch = "wasm32"))]
use sysinfo::System;

#[derive(Debug, PartialEq)]
pub enum AppTab {
    PQPGeneration,
    OpenSwath,
    Alignment,
    Validation,
    Visualization,
}

pub struct ArycalApp {
    current_tab: AppTab,
    toasts: Toasts,
    /// wraps serde­-loaded settings
    config: Input,
    config_file_path: String,
    /// States for each tab
    pqp_state: PQPState,
    osw_state: OpenSwathState,
    align_state: AlignmentState,
    validation_state: ValidationState,
    vis_state: VisualizationState,
    sidebar_width: f32,
    show_log: bool,
    log_messages: Arc<Mutex<Vec<String>>>,
    alignment_progress: Arc<Mutex<f32>>,

    /// System information
    #[cfg(not(target_arch = "wasm32"))]
    system: System,
    #[cfg(not(target_arch = "wasm32"))]
    last_sys_refresh: Instant,
}

impl ArycalApp {
    pub fn new(config_file_path: String) -> Self {
        // Try to load the config if a path is provided, otherwise use default
        let config = if !config_file_path.is_empty() {
            match Input::load(&config_file_path) {
                Ok(cfg) => {
                    log::info!("Loaded configuration from: {}", config_file_path);
                    cfg
                }
                Err(e) => {
                    log::warn!("Failed to load config '{}': {}. Using defaults.", config_file_path, e);
                    Input::default()
                }
            }
        } else {
            log::info!("No configuration file provided. Using defaults.");
            Input::default()
        };

        #[cfg(not(target_arch = "wasm32"))]
        let mut system = {
            let mut s = System::new_all();
            s.refresh_all();
            s
        };

        Self {
            current_tab: AppTab::Visualization,
            toasts: Toasts::default(),
            config,
            config_file_path,
            pqp_state: PQPState::default(),
            osw_state: OpenSwathState::default(),
            align_state: AlignmentState::default(),
            validation_state: ValidationState::default(),
            vis_state: VisualizationState::default(),
            sidebar_width: 250.0,
            show_log: true,
            log_messages: Arc::new(Mutex::new(Vec::new())),
            alignment_progress: Arc::new(Mutex::new(0.0)),
            #[cfg(not(target_arch = "wasm32"))]
            system,
            #[cfg(not(target_arch = "wasm32"))]
            last_sys_refresh: Instant::now(),
        }
    }
}


impl App for ArycalApp {
    fn update(&mut self, ctx: &eframe::egui::Context, _frame: &mut eframe::Frame) {
        // Load image loaders
        egui_extras::install_image_loaders(ctx);

        // Side panel with nested collapsibles
        SidePanel::left("config_panel")
            .default_width(self.sidebar_width)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {

                    match self.current_tab {
                        AppTab::PQPGeneration  => {
                            ui.add(
                                egui::Image::new(egui::include_image!(
                                    "../../../assets/img/pqp_logo_transparent_bg.png"
                                ))
                                .corner_radius(5),
                            );
                        },
                        AppTab::OpenSwath  => {
                            ui.add(
                                egui::Image::new(egui::include_image!(
                                    "../../../assets/img/openswath.png"
                                ))
                                .corner_radius(5),
                            );
                        },
                        AppTab::Validation  => {
                            ui.add(
                                egui::Image::new(egui::include_image!(
                                    "../../../assets/img/PyProphet_Logo_transparent_bg.png"
                                ))
                                .corner_radius(5),
                            );
                        },
                        _ => {
                            ui.add(
                                egui::Image::new(egui::include_image!(
                                    "../../../assets/img/arycal_logo_submark_transparent.png"
                                ))
                                .corner_radius(5),
                            );
                        }
                    }

                    ui.heading("Configuration");

                    // Horizontal buttons to Load or Save config
                    ui.horizontal(|ui| {
                        if ui.button("Load Config").clicked() {
                            if let Some(path) = FileDialog::new()
                                .set_title("Load Configuration")
                                .add_filter("JSON", &["json"])
                                .pick_file()
                            {
                                match Input::load(&path.to_string_lossy()) {
                                    Ok(config) => {
                                        self.config = config;
                                        self.config_file_path = path.to_string_lossy().to_string();
                                        self.toasts.success("Configuration loaded successfully.");
                                    }
                                    Err(e) => {
                                        error!("Failed to load config: {}", e);
                                        self.toasts.error(format!("Failed to load config: {}", e));
                                    }
                                }
                            }
                        }

                        if ui.button("Save Config").clicked() {
                            // Show a save file dialog with default filename
                            if let Some(path) = rfd::FileDialog::new()
                                .set_title("Save Configuration")
                                .add_filter("Config Files", &["json"])
                                .set_file_name("arycal_gui_config.json")  
                                .save_file()
                            {
                                // Try saving to the chosen path
                                match self.config.save_to_file(&path) {
                                    Ok(_) => {
                                        self.toasts.success("Configuration saved successfully.");
                                    }
                                    Err(e) => {
                                        error!("Failed to save config: {}", e);
                                        self.toasts.error(format!("Failed to save config: {}", e));
                                    }
                                }
                             
                            }
                        }
                        
                    });

                    ui.add_space(10.0);
                    ui.separator();

                    // Shared settings always visible
                    match self.current_tab {
                        AppTab::PQPGeneration | AppTab::OpenSwath | AppTab::Alignment
                        | AppTab::Validation | AppTab::Visualization => {
                            egui::CollapsingHeader::new("File Settings")
                                .default_open(true)
                                .show(ui, |ui| {
                                    config_panel::draw_shared(ui, &self.current_tab, &mut self.config);
                                });
                        },
                        _ => {
                            ui.label("");
                        }
                    }

                    // Tab-specific settings
                    match self.current_tab {
                        AppTab::PQPGeneration => {
                            egui::CollapsingHeader::new("PQP Generation Settings")
                                .default_open(true)
                                .show(ui, |ui| {
                                    draw_pqp_generation(ui, &mut self.config);
                                });
                        }
                        AppTab::OpenSwath => {
                            egui::CollapsingHeader::new("OpenSwath Settings")
                                .default_open(true)
                                .show(ui, |ui| {
                                    draw_open_swath(ui, &mut self.config, &mut self.osw_state);
                                });
                        }
                        AppTab::Alignment => {
                            egui::CollapsingHeader::new("Alignment Settings")
                                .default_open(true)
                                .show(ui, |ui| {
                                    draw_alignment(ui, &mut self.toasts, &mut self.config);
                                });
                        }
                        AppTab::Validation => {
                            egui::CollapsingHeader::new("Validation Settings")
                                .default_open(true)
                                .show(ui, |ui| {
                                    draw_validation(ui, &mut self.config);
                                });
                        }
                        AppTab::Visualization => {
                            egui::CollapsingHeader::new("Visualization Settings")
                                .default_open(true)
                                .show(ui, |ui| {
                                    visualization_settings::draw_visualization_settings(
                                        ui,
                                        &mut self.vis_state,
                                        &mut self.config,
                                    );
                                });
                        }
                    }

                    // System settings
                    egui::CollapsingHeader::new("System Settings")
                        .default_open(true)
                        .show(ui, |ui| {
                            #[cfg(not(target_arch = "wasm32"))]
                            {
                                // the real sysinfo drawer:
                                system_footer::draw_system_settings(
                                    &mut self.system,
                                    &mut self.last_sys_refresh,
                                    ui,
                                    &mut self.config,
                                );
                            }
                            #[cfg(target_arch = "wasm32")]
                            {
                                // stub — do nothing in wasm
                                system_footer::draw_system_settings_stub(ui, &mut self.config);
                            }
                        });

                    // Footer
                    ui.add_space(4.0);
                    ui.separator();
                    ui.horizontal(|ui| {
                        ui.hyperlink_to(
                            format!("{}  ARYCAL on GitHub", egui::special_emojis::GITHUB),
                            "https://github.com/singjc/arycal",
                        );

                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.label(format!("Version {}", clap::crate_version!()));
                        });
                    });
                });
            });

        // Bottom status metrics panel. 
        TopBottomPanel::bottom("status_bar")
            .exact_height(20.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    match self.current_tab {
                        AppTab::PQPGeneration => {
                            let pqp = &self.pqp_state;
        
                            let is_running = pqp.is_running(); // new helper on PQPState
                            // Extract some common timing fields
                            let tavg = pqp.avg_processing_time;
                            let ttotal = pqp.total_processing_time;
        
                            if is_running {
                                // Safe to unwrap current_stage_runner because is_running() == true
                                let stage_runner = pqp.current_stage_runner.as_ref().unwrap();
                                let stage = &stage_runner.stage;
        
                                let stage_idx = pqp.current_stage_index + 1;
                                let n_stages = pqp.pipeline.len();
        
                                // Per-stage progress (fraction inside current stage)
                                let stage_done = pqp.stage_successful_runs + pqp.stage_failed_runs;
                                let stage_total = pqp.total_number_of_runs;
                                let stage_fraction = if stage_total > 0 {
                                    stage_done as f32 / stage_total as f32
                                } else { 0.0 };
        
                                // Overall progress
                                let overall_fraction = pqp.overall_progress_fraction();
                                let overall_done = pqp.successful_runs + pqp.failed_runs;
                                let overall_total = pqp.overall_total_runs;
        
                                ui.label(format!(
                                    "Stage {}/{}: {} ({}/{} runs) | Overall {}/{}",
                                    stage_idx,
                                    n_stages,
                                    stage.display_name,
                                    stage_done,
                                    stage_total,
                                    overall_done,
                                    overall_total
                                ));
        
                                ui.label(format!(
                                    "Stage Avg: {:.1}m Total: {:.1}m",
                                    tavg.as_secs_f64() / 60.0,
                                    ttotal.as_secs_f64() / 60.0
                                ));
        
                                // Right side: two progress bars (overall + stage) aligned right
                                // Some spacing
                                ui.add_space(600.0);
                                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                    // Stage progress
                                    ui.add(
                                        egui::ProgressBar::new(stage_fraction)
                                            .text(format!("Stage {:>3.0}%", stage_fraction * 100.0))
                                            .animate(true)
                                    );
                                    // Some spacing
                                    ui.add_space(2.0);
                                    // Overall progress
                                    ui.add(
                                        egui::ProgressBar::new(overall_fraction)
                                            .text(format!("Overall {:>3.0}%", overall_fraction * 100.0))
                                            .animate(true)
                                    );
                                });
        
                            } else {
                                // Not currently running
                                if pqp.overall_total_runs > 0 && (pqp.successful_runs + pqp.failed_runs) == pqp.overall_total_runs {
                                    // Completed pipeline state
                                    ui.label(format!(
                                        "PQP Generation complete: {} success, {} failed ({} total runs, {} stages). Last Stage Avg: {:.1}m Total: {:.1}m",
                                        pqp.successful_runs,
                                        pqp.failed_runs,
                                        pqp.overall_total_runs,
                                        pqp.pipeline.len(), // after completion pipeline is usually cleared; if cleared use stored count
                                        tavg.as_secs_f64() / 60.0,
                                        ttotal.as_secs_f64() / 60.0
                                    ));
                                } else if pqp.overall_total_runs > 0 && (pqp.successful_runs + pqp.failed_runs) > 0 {
                                    // Partially done (e.g., cancelled)
                                    let overall_fraction = pqp.overall_progress_fraction();
                                    ui.label(format!(
                                        "Pipeline stopped: {} / {} completed ({} success, {} failed).",
                                        pqp.successful_runs + pqp.failed_runs,
                                        pqp.overall_total_runs,
                                        pqp.successful_runs,
                                        pqp.failed_runs
                                    ));
                                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                        ui.add(
                                            egui::ProgressBar::new(overall_fraction)
                                                .text(format!("{:>3.0}%", overall_fraction * 100.0))
                                        );
                                    });
                                } else {
                                    // Idle / never started
                                    ui.label("PQP Generation idle.");
                                }
                            }
                        },

                        AppTab::OpenSwath => {
                            let is_running = self.osw_state.runner.is_some();

                            if is_running {
                                ui.label(format!(
                                    "Processing {}/{} files",
                                    self.osw_state.processing_run_n,
                                    self.osw_state.total_number_of_runs
                                ));
    
                                let tavg = &self.osw_state.avg_processing_time;
                                let ttotal = &self.osw_state.total_processing_time;
                                ui.label(format!(
                                    "Avg Time: {:.2?}min, Total Time: {:.2?}min",
                                    tavg.as_secs() / 60,
                                    ttotal.as_secs() / 60
                                ));
    
                                // Right side: progress bar
                                ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                                    // Get progress based on processing run number and total runs
                                    let progress = self.osw_state.processed_fraction;
                                    ui.add(
                                        egui::ProgressBar::new(progress)
                                            .text(format!("{:.0}%", progress * 100.0))
                                            .animate(true),
                                    );
                                });
                            } else {
                                if self.osw_state.successful_runs > 0 || self.osw_state.failed_runs > 0 {
                                    let tavg = &self.osw_state.avg_processing_time;
                                    let ttotal = &self.osw_state.total_processing_time;
                                    ui.label(format!(
                                        "OpenSwath completed. {} successful, {} failed. Avg Time: {:.2?}min, Total Time: {:.2?}min",
                                        self.osw_state.successful_runs,
                                        self.osw_state.failed_runs,
                                        tavg.as_secs() / 60,
                                        ttotal.as_secs() / 60
                                    ));
                                } else {
                                    ui.label(format!(
                                        "OpenSwath not running. {} files to process.",
                                        self.osw_state.total_number_of_runs
                                    ));
                                }
                            }
                        },
                        AppTab::Alignment => {
                            let a = &self.align_state;
                            let is_running = a.io_handle.is_some();
                            let ttotal = a.total_processing_time;
        
                            if is_running {
                                ui.label(format!(
                                    "Alignment processing {} files with a batch size of {}",
                                    a.n_runs,
                                    a.batch_size
                                ));
                            } else if ttotal > Duration::ZERO {
                                ui.label(format!(
                                    "Alignment complete in {:.1}m",
                                    ttotal.as_secs_f64() / 60.0
                                ));
                            } else {
                                ui.label("Alignment idle.");
                            }
                        },
                        AppTab::Validation => {
                            let vs = &self.validation_state;
                            let is_running = vs.current_stage_runner.is_some();

                            if is_running {
                                // Safe: stage runner exists
                                let runner = vs.current_stage_runner.as_ref().unwrap();
                                let stage = &runner.stage;
                                let stage_idx = vs.current_stage_index + 1;
                                let n_stages = vs.pipeline.len();

                                // per-stage progress
                                let stage_done = vs.stage_successful_runs + vs.stage_failed_runs;
                                let stage_total = vs.total_number_of_runs;
                                let stage_frac = if stage_total > 0 {
                                    stage_done as f32 / stage_total as f32
                                } else {
                                    0.0
                                };

                                // overall progress
                                let overall_done = vs.successful_runs + vs.failed_runs;
                                let overall_total = vs.overall_total_runs;
                                let overall_frac = if overall_total > 0 {
                                    overall_done as f32 / overall_total as f32
                                } else {
                                    0.0
                                };

                                // timings
                                let tavg = vs.avg_processing_time;
                                let ttotal = vs.total_processing_time;

                                ui.label(format!(
                                    "Stage {}/{}: {} ({}/{} runs) | Overall {}/{}",
                                    stage_idx, n_stages, stage.display_name,
                                    stage_done, stage_total,
                                    overall_done, overall_total
                                ));
                                ui.label(format!(
                                    "Stage Avg: {:.1}m Total: {:.1}m",
                                    tavg.as_secs_f64() / 60.0,
                                    ttotal.as_secs_f64() / 60.0
                                ));

                                // push the bars to the right
                                ui.add_space(600.0);
                                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                    ui.add(
                                        egui::ProgressBar::new(stage_frac)
                                            .text(format!("Stage {:>3.0}%", stage_frac * 100.0))
                                            .animate(true),
                                    );
                                    ui.add_space(2.0);
                                    ui.add(
                                        egui::ProgressBar::new(overall_frac)
                                            .text(format!("Overall {:>3.0}%", overall_frac * 100.0))
                                            .animate(true),
                                    );
                                });
                            } else {
                                // completed, cancelled, or idle
                                let done = vs.successful_runs + vs.failed_runs;

                                if vs.overall_total_runs > 0 && done == vs.overall_total_runs {
                                    ui.label(format!(
                                        "Validation complete: {} success, {} failed ({} runs, {} stages). Last Stage Avg: {:.1}m Total: {:.1}m",
                                        vs.successful_runs,
                                        vs.failed_runs,
                                        vs.overall_total_runs,
                                        vs.pipeline.len(),
                                        vs.avg_processing_time.as_secs_f64() / 60.0,
                                        vs.total_processing_time.as_secs_f64() / 60.0
                                    ));
                                } else if vs.overall_total_runs > 0 && done > 0 {
                                    let frac = if vs.overall_total_runs > 0 {
                                        done as f32 / vs.overall_total_runs as f32
                                    } else {
                                        0.0
                                    };
                                    ui.label(format!(
                                        "Validation stopped: {}/{} completed ({} success, {} failed).",
                                        done,
                                        vs.overall_total_runs,
                                        vs.successful_runs,
                                        vs.failed_runs
                                    ));
                                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                        ui.add(
                                            egui::ProgressBar::new(frac)
                                                .text(format!("{:>3.0}%", frac * 100.0)),
                                        );
                                    });
                                } else {
                                    ui.label("Validation idle.");
                                }
                            }
                        }
                        AppTab::Visualization => {
                            // Left side: timings
                            let t = &self.vis_state.precursor_load_time;
                            ui.label(format!("⏱ Loaded Peptides: {}ms, ", t.as_millis()));

                            let t = &self.vis_state.xic_load_time;
                            let t2 = &self.vis_state.plot_render_time;
                            ui.label(format!(
                                "XIC Load: {}µs, Plot Render: {}µs",
                                t.as_micros(),
                                t2.as_micros()
                            ));

                            // Right side: counts
                            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                                ui.label(format!(
                                    "{} precursors   {} peptides",
                                    self.vis_state.num_unique_precursors,
                                    self.vis_state.num_unique_peptides
                                ));
                            });
                        }
                        _ => {
                            // For other tabs, just show the current tab name
                            ui.label(format!("{:?}", self.current_tab));
                        }
                    } // End match
                });
            });

        // Main area with tabs
        CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui
                    .selectable_label(self.current_tab == AppTab::PQPGeneration, "PQP Generation")
                    .clicked()
                {
                    self.current_tab = AppTab::PQPGeneration;
                }
                if ui
                    .selectable_label(self.current_tab == AppTab::OpenSwath, "OpenSwath")
                    .clicked()
                {
                    self.current_tab = AppTab::OpenSwath;
                }
                if ui
                    .selectable_label(self.current_tab == AppTab::Alignment, "Alignment")
                    .clicked()
                {
                    self.current_tab = AppTab::Alignment;
                }
                if ui
                    .selectable_label(self.current_tab == AppTab::Validation, "Validation")
                    .clicked()
                {
                    self.current_tab = AppTab::Validation;
                }
                if ui
                    .selectable_label(self.current_tab == AppTab::Visualization, "Visualization")
                    .clicked()
                {
                    self.current_tab = AppTab::Visualization;
                }
            });
            ui.separator();

            match self.current_tab {
                AppTab::PQPGeneration => {
                    self.pqp_state.ui(ui, &self.config);
                }
                AppTab::OpenSwath => self.osw_state.ui(ui, &self.config),
                AppTab::Alignment => self.align_state.ui(ui, &self.config),
                AppTab::Validation => self.validation_state.ui(ui, &self.config),
                AppTab::Visualization => self.vis_state.ui(ui, &mut self.config),
            }
        });

        self.toasts.show(ctx);
    }
}
