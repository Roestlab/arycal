// arycal-gui/src/main.rs
use arycal_gui::ArycalApp;
use eframe::egui;

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "Arycal GUI",
        options,
        Box::new(|_cc| Ok(Box::new(ArycalApp::new("config_arycal.json".to_string())))),
    )
}
