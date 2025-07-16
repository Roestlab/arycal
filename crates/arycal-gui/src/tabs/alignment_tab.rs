use std::sync::Arc;
use std::sync::Mutex;
use eframe::{egui::{Ui, ProgressBar}};
use crate::app::ArycalApp;
use arycal_cli::input::Input;


/// Actions that the Alignment tab can request
pub enum AlignmentAction { Run, ShowLog }

pub struct AlignmentState; // no internal mutable state needed for now

impl AlignmentState {
    pub fn default() -> Self { Self }

    /// Renders the Alignment tab.
    /// Returns any requested actions for the host app to execute.
    pub fn ui(
        &self,
        ui: &mut Ui,
        progress: &Arc<Mutex<f32>>,
        logs: &Arc<Mutex<Vec<String>>>
    ) -> Vec<AlignmentAction> {
        let mut actions = Vec::new();
        ui.horizontal(|ui| {
            if ui.button("▶ Run Alignment").clicked() {
                actions.push(AlignmentAction::Run);
            }
            if ui.button("Show Log").clicked() {
                actions.push(AlignmentAction::ShowLog);
            }
        });
        ui.separator();
        ui.heading("Alignment Status");

        // Always show progress bar
        let p = *progress.lock().unwrap();
        ui.add(ProgressBar::new(p).animate(true));

        // Show recent log entries
        for msg in logs.lock().unwrap().iter() {
            ui.label(msg);
        }

        actions
    }

    /// Renders the shared config panel (previously in app.rs)
    pub fn edit_config(&self, ui: &mut Ui, config: &mut Input) {
        // Copy your edit_config and edit_file_paths implementations here
    }
}

