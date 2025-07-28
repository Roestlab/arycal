use std::time::{Duration, Instant};
#[cfg(not(target_arch = "wasm32"))]
use sysinfo::System;
use eframe::egui::Ui;
use egui::ProgressBar;

use arycal_cli::input::Input;


/// Draw the system settings panel.
///
/// - `system` is your cached `sysinfo::System`
/// - `last_refresh` is the Instant you last did a refresh_cpu/refresh_memory
/// - `ui` is the egui Ui to draw into
#[cfg(not(target_arch = "wasm32"))]
pub fn draw_system_settings(
    system: &mut System,
    last_refresh: &mut Instant,
    ui: &mut Ui,
    config: &mut Input,
) {
    // Theme buttons

    use crate::utils::update_python_path;
    egui::widgets::global_theme_preference_buttons(ui);

    // Concurrent processes input
    ui.horizontal(|ui| {
        ui.label("Concurrent Processes:");
        let cap_num_processes = std::thread::available_parallelism().unwrap().get();
        ui.add(
            egui::DragValue::new(&mut config.n_concurrent_processes)
                .speed(1.0)
                .range(1..=cap_num_processes),
        )
        .on_hover_text("Set the number of concurrent processes...");
    });

    ui.add_space(2.0);

    // Python executable selector
    ui.horizontal(|ui| {
        ui.label("Python Executable:");
        let mut py_display = config
            .python_path
            .clone()
            .unwrap_or_else(|| "System Default".to_string());
        if ui.text_edit_singleline(&mut py_display).changed() {
            config.python_path = if py_display.trim().is_empty() {
                None
            } else {
                Some(py_display.trim().to_string())
            };
            update_python_path(&config.python_path, &mut config.last_python_dir);
        }
        if ui.button("…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .set_title("Select Python Executable")
                .pick_file()
            {
                config.python_path = Some(path.display().to_string());
                update_python_path(&config.python_path, &mut config.last_python_dir);
            }
        }
    });

    ui.add_space(4.0);

    // Refresh system info
    if last_refresh.elapsed() >= Duration::from_secs(1) {
        system.refresh_cpu_all();
        system.refresh_memory();
        *last_refresh = Instant::now();
    }
    ui.ctx().request_repaint_after(Duration::from_secs(1));

    // System info panel
    egui::CollapsingHeader::new("System Info")
        .default_open(true)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label("OS:");
                ui.label(format!(
                    "{} {}",
                    System::name().unwrap_or_default(),
                    System::os_version().unwrap_or_default()
                ));
            });

            let cpu = system.global_cpu_usage();
            ui.horizontal(|ui| {
                ui.label("CPU:");
                ui.add(ProgressBar::new(cpu as f32 / 100.0).text(format!("{:.1}%", cpu)));
            });

            let total = system.total_memory() as f64 / 1024.0 / 1024.0 / 1024.0;
            let used = system.used_memory() as f64 / 1024.0 / 1024.0 / 1024.0;
            ui.horizontal(|ui| {
                ui.label("Memory:");
                ui.add(
                    ProgressBar::new((used / total) as f32)
                        .text(format!("{:.1} / {:.1} GB", used, total)),
                );
            });
        });
}



/// Draw the system settings panel.
///
/// - `system` is your cached `sysinfo::System`
/// - `last_refresh` is the Instant you last did a refresh_cpu/refresh_memory
/// - `ui` is the egui Ui to draw into
#[cfg(target_arch = "wasm32")]
pub fn draw_system_settings(
    ui: &mut Ui,
) {
    // theme buttons
    egui::widgets::global_theme_preference_buttons(ui);
    ui.add_space(4.0);
    // ask egui to repaint again in 1s
    ui.ctx().request_repaint_after(Duration::from_secs(1));
    
    ui.add_space(4.0);
    ui.separator();
    ui.hyperlink_to(
        format!("{}  ARYCAL on GitHub", egui::special_emojis::GITHUB),
        "https://github.com/singjc/arycal",
    );
}