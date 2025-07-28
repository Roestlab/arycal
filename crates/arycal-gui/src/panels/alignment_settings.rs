use arycal_cli::input::Input;
use arycal_cloudpath::util::find_executable;
use eframe::egui::{Ui, ComboBox, Slider, TextEdit};
use egui_notify::{Toast, Toasts};
use std::{default, env, path::PathBuf};
use arycal_common::config::{self, OpenSwathConfig, PQPConfig, PyProphetConfig};

use crate::utils::{get_install_dir, install_latest_arycal};


/// Draw additional settings for Alignment tab
pub fn draw_alignment(ui: &mut Ui, toasts: &mut Toasts, config: &mut Input) {
    ui.heading("Alignment Settings");

    // — binary path + auto-detect button —
    let default_binary_directory = get_install_dir()
        .expect("Failed to get default binary directory");

    if config.arycal_binary_path.is_none() {
        if let Some(p) = find_executable("arycal", Some(default_binary_directory.as_path())) {
            config.arycal_binary_path = Some(p);
        } 
    }

    ui.horizontal(|ui| {
        ui.add(egui::Label::new("arycal binary:")).on_hover_text(
            "Path to the arycal binary. If not set, it will try to find it in your PATH."
        );

        // Show & edit it as a String
        let mut path_str = config.arycal_binary_path.clone().unwrap_or_default().to_string_lossy().into_owned();
        if ui.add(TextEdit::singleline(&mut path_str)).changed() {
            config.arycal_binary_path = Some(PathBuf::from(path_str.clone()));
        }

        // File dialog button
        if ui.button("…").on_hover_text("Browse for binary").clicked() {
            if let Some(file) = rfd::FileDialog::new()
                .set_title("Select arycal binary")
                // start in current dir or last used
                .set_directory(env::current_dir().unwrap_or_default())
                .pick_file() 
            {
                config.arycal_binary_path = Some(file);
            }
        }

        // Button to install the binary
        if ui.button("Install").on_hover_text("Install arycal binary").clicked() {
            toasts.info("Installing latest arycal binary...").closable(false).duration(None);
            match install_latest_arycal() {
                Ok(_) => {
                    toasts.dismiss_latest_toast();
                    config.arycal_binary_path = Some(get_install_dir().unwrap().join("arycal"));
                    toasts.success(format!("Installed successfully to:\n{}", config.arycal_binary_path.as_ref().unwrap().display())).closable(true).duration(None);
                    log::info!("Installed arycal binary to: {}", config.arycal_binary_path.as_ref().unwrap().display());
                },
                Err(e) => {
                    toasts.dismiss_latest_toast();
                    toasts.error(format!("Failed to install:\n{}", e))
                        .closable(true)
                        .duration(None);
                    log::error!("Failed to install arycal binary: {}", e);
                }
            }
        }
    });

    ui.separator();

    ui.horizontal(|ui| {
        ui.label("Batch Size:");
        let mut bs = config.alignment.batch_size.unwrap_or(0).to_string();
        if ui.add(TextEdit::singleline(&mut bs)).changed() {
            if let Ok(val) = bs.parse::<usize>() {
                config.alignment.batch_size = Some(val);
            }
        }
    });

    ui.horizontal(|ui| {
        ui.label("Method:");
        ComboBox::from_id_salt("alignment_method")
            .selected_text(config.alignment.method.clone())
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut config.alignment.method, "DTW".into(), "DTW");
                ui.selectable_value(&mut config.alignment.method, "FFT".into(), "FFT");
                ui.selectable_value(&mut config.alignment.method, "FFT+DTW".into(), "FFT+DTW");
            });
    });
    ui.label("Reference Type:");
    ComboBox::from_id_salt("reference_type")
        .selected_text(config.alignment.reference_type.clone())
        .show_ui(ui, |ui| {
            ui.selectable_value(&mut config.alignment.reference_type, "star".into(), "STAR");
            ui.selectable_value(&mut config.alignment.reference_type, "mst".into(), "MST");
            ui.selectable_value(&mut config.alignment.reference_type, "progressive".into(), "PROGRESSIVE");
        });
    ui.label("Reference File:");
    ComboBox::from_id_salt("reference_file")
        .selected_text(config.alignment.reference_run.clone().unwrap_or_default())
        .show_ui(ui, |ui| {
            for path in config.xic.file_paths.iter() {
                if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                    ui.selectable_value(&mut config.alignment.reference_run, Some(name.to_string()), name);
                }
            }
        });
    ui.checkbox(&mut config.alignment.use_tic, "Use TIC");
    ui.add(Slider::new(&mut config.alignment.smoothing.sgolay_window, 1..=51).text("Savitzky-Golay Window"));
    ui.add(Slider::new(&mut config.alignment.smoothing.sgolay_order, 1..=21).text("Savitzky-Golay Order"));
    ui.label("RT Mapping Tolerance:");
    let mut rt_tol = config.alignment.rt_mapping_tolerance.unwrap_or_default().to_string();
    if ui.add(TextEdit::singleline(&mut rt_tol)).changed() {
        if let Ok(val) = rt_tol.parse::<f64>() {
            config.alignment.rt_mapping_tolerance = Some(val);
        }
    }
    ui.label("Decoy Peak Mapping Method:");
    ComboBox::from_id_salt("decoy_peak_mapping_method")
        .selected_text(config.alignment.decoy_peak_mapping_method.clone())
        .show_ui(ui, |ui| {
            ui.selectable_value(&mut config.alignment.decoy_peak_mapping_method, "shuffle".into(), "shuffle");
            ui.selectable_value(&mut config.alignment.decoy_peak_mapping_method, "random_regions".into(), "random_regions");
        });
    ui.label("Decoy Window Size:");
    let mut dw = config.alignment.decoy_window_size.unwrap_or_default().to_string();
    if ui.add(TextEdit::singleline(&mut dw)).changed() {
        if let Ok(val) = dw.parse::<usize>() {
            config.alignment.decoy_window_size = Some(val);
        }
    }
    ui.add_space(25.0);
    ui.separator();

    // Add numeric input for the number of threads and a select box for the log level for log::info, debug, trace
    ui.horizontal(|ui| {
        ui.label("Threads:");
        let mut threads = config.threads.to_string();
        if ui.add(TextEdit::singleline(&mut threads)).changed() {
            if let Ok(val) = threads.parse::<usize>() {
                config.threads = val;
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("Log Level:");
        ComboBox::from_id_salt("log_level")
            .selected_text(config.log_level.clone())
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut config.log_level, "info".into(), "Info");
                ui.selectable_value(&mut config.log_level, "debug".into(), "Debug");
                ui.selectable_value(&mut config.log_level, "trace".into(), "Trace");
            });
    });
}