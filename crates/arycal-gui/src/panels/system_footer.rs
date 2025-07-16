use std::time::{Duration, Instant};
#[cfg(not(target_arch = "wasm32"))]
use sysinfo::System;
use eframe::egui::Ui;
use egui::ProgressBar;


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
) {
    // theme buttons
    egui::widgets::global_theme_preference_buttons(ui);
    ui.add_space(4.0);

    // only refresh once a second:
    if last_refresh.elapsed() >= Duration::from_secs(1) {
        system.refresh_cpu_all();
        system.refresh_memory();
        *last_refresh = Instant::now();
    }
    // ask egui to repaint again in 1s
    ui.ctx().request_repaint_after(Duration::from_secs(1));

    egui::CollapsingHeader::new("System Info").show(ui, |ui| {
        // OS information
        ui.horizontal(|ui| {
            ui.label("OS:");
            ui.label(format!("{} {}", System::name().unwrap_or_default(), System::os_version().unwrap_or_default()));
        });

        let cpu = system.global_cpu_usage();
        ui.horizontal(|ui| {
            ui.label("CPU:");
            ui.add(
                ProgressBar::new(cpu as f32 / 100.0)
                    .text(format!("{:.1}%", cpu)),
            );
        });

        let total = system.total_memory() as f64 / 1024.0 / 1024.0 / 1024.0;
        let used  = system.used_memory()  as f64 / 1024.0 / 1024.0 / 1024.0;
        ui.horizontal(|ui| {
            ui.label("Memory:");
            ui.add(
                ProgressBar::new((used/total) as f32)
                    .text(format!("{:.1} / {:.1} GB", used, total)),
            );
        });
    });

    ui.add_space(4.0);
    ui.separator();
    ui.hyperlink_to(
        format!("{}  ARYCAL on GitHub", egui::special_emojis::GITHUB),
        "https://github.com/singjc/arycal",
    );
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