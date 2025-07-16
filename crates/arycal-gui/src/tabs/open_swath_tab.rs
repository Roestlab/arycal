use eframe::egui::Ui;
use arycal_cli::input::Input;

pub struct OpenSwathState {
    // binary path, flags, etc.
}

impl OpenSwathState {
    pub fn default() -> Self { Self {} }
    pub fn ui(&mut self, ui: &mut Ui, _config: &Input) {
        ui.label("OpenSwathWorkflow tab: set binary path & params");
        // ... integrate your runner invocation UI code here
    }
}