
use std::time::{Duration, Instant};
#[cfg(not(target_arch = "wasm32"))]
use sysinfo::System;
use eframe::{egui::{CentralPanel, SidePanel, Ui}, App};
use egui::{Align, Layout, TopBottomPanel};
use rfd::FileDialog;
use std::{fs, sync::{Arc, Mutex}, thread};
use log::{error, info};
use crate::config::AppConfig;
use arycal_cli::{input::Input, Runner};
use crate::tabs::{visualization_tab::VisualizationState, open_swath_tab::OpenSwathState, alignment_tab::{AlignmentState, AlignmentAction}};
use crate::panels::{config_panel, visualization_settings, system_footer::draw_system_settings};

#[derive(PartialEq)]
enum AppTab {
    Visualization,
    OpenSwath,
    Alignment,
}

pub struct ArycalApp {
    current_tab: AppTab,
    /// wraps serde­-loaded settings
    config: Input,
    config_file_path: String,
    vis_state: VisualizationState,
    osw_state: OpenSwathState,
    align_state: AlignmentState,
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
        let config = Input::default();
        #[cfg(not(target_arch = "wasm32"))]
        let mut system = {
            let mut s = System::new_all();
            s.refresh_all();
            s
        };

        Self {
            current_tab: AppTab::Visualization,
            config,
            config_file_path,
            vis_state: VisualizationState::default(),
            osw_state: OpenSwathState::default(),
            align_state: AlignmentState::default(),
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

    fn run_alignment(&mut self, ctx: &eframe::egui::Context) -> Result<(), Box<dyn std::error::Error>> {
        let logs = Arc::clone(&self.log_messages);
        let cfg = Arc::new(self.config.clone());
        *self.alignment_progress.lock().unwrap() = 0.0;
        let progress = Arc::clone(&self.alignment_progress);

        thread::spawn(move || {
            let mut runner = Runner::new((*cfg).clone(), Some(progress)).unwrap();
            match runner.run() {
                Ok(_) => {
                    logs.lock().unwrap().push("Alignment completed successfully.".into());
                    info!("Alignment done.");
                }
                Err(e) => {
                    let msg = format!("Error: {}", e);
                    logs.lock().unwrap().push(msg.clone());
                    error!("{}", msg);
                }
            }
        });

        Ok(())
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
                ui.add(
                    egui::Image::new(
                        egui::include_image!("../../../assets/img/arycal_logo_submark_transparent.png")
                    )
                    .corner_radius(5)
                );
                ui.heading("Configuration");
                ui.add_space(10.0);
                ui.separator();

                // Shared settings always visible
                egui::CollapsingHeader::new("Shared Settings")
                    .default_open(true)
                    .show(ui, |ui| {
                        config_panel::draw_shared(ui, &mut self.config);
                    });

                // Tab-specific settings
                match self.current_tab {
                    AppTab::Visualization => {
                        egui::CollapsingHeader::new("Visualization Settings")
                            .default_open(true)
                            .show(ui, |ui| {
                                visualization_settings::draw_visualization_settings(ui, &mut self.vis_state, &mut self.config);
                            });
                    }
                    AppTab::OpenSwath => {
                        egui::CollapsingHeader::new("OpenSwath Settings")
                            .default_open(true)
                            .show(ui, |ui| {
                                config_panel::draw_open_swath(ui, &mut self.config);
                            });
                    }
                    AppTab::Alignment => {
                        egui::CollapsingHeader::new("Alignment Settings")
                            .default_open(true)
                            .show(ui, |ui| {
                                config_panel::draw_alignment(ui, &mut self.config);
                            });
                    }
                }

                // System settings
                egui::CollapsingHeader::new("System Settings")
                    .default_open(true)
                    .show(ui, |ui| {
                        #[cfg(not(target_arch = "wasm32"))] {
                            // the real sysinfo drawer:
                            panels::system_footer::draw_system_settings(
                                &mut self.system,
                                &mut self.last_sys_refresh,
                                ui,
                            );
                        }
                        #[cfg(target_arch = "wasm32")] {
                            // stub — do nothing in wasm
                            panels::system_footer::draw_system_settings_stub(ui);
                        }
                    });
            });
        });

        // Main area with tabs
        CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui
                    .selectable_label(self.current_tab == AppTab::Visualization, "Visualization")
                    .clicked()
                {
                    self.current_tab = AppTab::Visualization;
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
            });
            ui.separator();

            match self.current_tab {
                AppTab::Visualization => self.vis_state.ui(ui, &mut self.config),
                AppTab::OpenSwath => self.osw_state.ui(ui, &self.config),
                AppTab::Alignment => {
                    let actions = self.align_state.ui(ui, &self.alignment_progress, &self.log_messages);
                    for action in actions {
                        match action {
                            AlignmentAction::Run => { let _ = self.run_alignment(ctx); }
                            AlignmentAction::ShowLog => { self.show_log = true; }
                        }
                    }
                }
            }
        });

        // Bottom status metrics panel
        TopBottomPanel::bottom("status_bar")
            .exact_height(20.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
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
                });
            });


    }
}



