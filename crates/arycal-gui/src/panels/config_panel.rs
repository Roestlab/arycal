use arycal_cli::input::Input;
use eframe::egui::{Ui, ComboBox, Slider, TextEdit};
use std::path::PathBuf;
use arycal_common::config::{PQPConfig, OpenSwathConfig, PyProphetConfig};

use crate::app::AppTab;
use crate::panels::peptide_query_parameter_settings::draw_pqp_file_settings;
use crate::panels::openswath_settings::draw_osw_file_settings;
use crate::panels::validation_settings::draw_validation_file_settings;
use crate::panels::visualization_settings::draw_visualization_file_settings;


/// Draw the shared XIC, Features, and Filters config
pub fn draw_shared(ui: &mut Ui, current_tab: &AppTab, config: &mut Input) {
    ui.heading("File Settings");
    ui.add_space(4.0);

    match current_tab {
        AppTab::Visualization | &AppTab::Alignment => {
            draw_visualization_file_settings(ui, config);
        },
        AppTab::PQPGeneration => {
            let pqp_cfg = config.pqp
                .get_or_insert_with(PQPConfig::default);
            draw_pqp_file_settings(ui, pqp_cfg);
        },
        AppTab::OpenSwath => {
            let osw_cfg: &mut OpenSwathConfig =
                config.openswath.get_or_insert_with(OpenSwathConfig::default);
            draw_osw_file_settings(ui, osw_cfg);
        },
        AppTab::Validation => {
            let val_cfg = config.statistical_validation
                .get_or_insert_with(PyProphetConfig::default);
            draw_validation_file_settings(ui, val_cfg);
        },
        _ => {
            // No additional file settings for PQP Generation tab
            ui.label("No file settings for this tab.");
        }
    }
}




/// Shared file paths drag & drop widget
pub fn edit_file_paths(
    ui: &mut egui::Ui,
    file_paths: &mut Vec<PathBuf>,
    button_name: &str,
    file_type_name: &str,
    file_dialog_title: Option<&str>,
    file_dialog_filter: Option<&[&str]>,
) {
    ui.group(|ui| {
        ui.vertical(|ui| {
            ui.label(format!("Drag & drop {} files here:", file_type_name));

            // 1) Reserve one rect for both the painter and the scrollable list:
            let max_zone_height = 150.0;
            let zone_width  = ui.available_width();
            let zone_height = ui.available_height().min(max_zone_height);
            let zone_size   = egui::vec2(zone_width, zone_height);
            let (resp, painter) = ui.allocate_painter(zone_size, egui::Sense::hover());

            // draw the border
            painter.rect_stroke(
                resp.rect,
                4.0,
                ui.visuals().widgets.active.bg_stroke,
                egui::epaint::StrokeKind::Middle,
            );

            // OS‐drop on the entire zone
            if resp.hovered() {
                let dropped = ui.ctx().input(|i| i.raw.dropped_files.clone());
                for file in dropped {
                    if let Some(path) = file.path {
                        file_paths.push(path);
                    }
                }
            }

            // 2) Now *inside* the same zone.rect*, place a ScrollArea + dnd_drop_zone:
            //    that way your list sits right on top of the painted box.
            ui.allocate_ui_at_rect(resp.rect, |ui| {
                egui::ScrollArea::vertical()
                    .max_height(zone_height)
                    .show(ui, |ui| {
                        // capacity for file_paths.len() rows
                        let needed = file_paths.len() as f32 * 24.0;
                        let rect = ui.available_rect_before_wrap();
                        ui.set_min_size(egui::vec2(rect.width(), needed));

                        // wrap everything in a single drop‐zone so both OS‐drops
                        // and DnD‐reorder work together:
                        let frame = egui::Frame::NONE.inner_margin(4.0);
                        ui.dnd_drop_zone::<usize, ()>(frame, |ui| {
                            for row in 0..file_paths.len() {
                                let id   = ui.make_persistent_id(("files", row));
                                let path = &file_paths[row];
                                let raw  = path.to_string_lossy();
                                let wrapped = raw.replace("/", "/\u{200B}");

                                // drag source carrying its index:
                                let drag_resp =
                                    ui.dnd_drag_source(id, row, |ui| {
                                        ui.horizontal_wrapped(|ui| {
                                            ui.label("☰");
                                            ui.label(wrapped.clone());
                                        });
                                    })
                                    .response;

                                // preview insertion line
                                if let (Some(pos), Some(_)) = (
                                    ui.input(|i| i.pointer.interact_pos()),
                                    drag_resp.dnd_hover_payload::<usize>(),
                                ) {
                                    let r = drag_resp.rect;
                                    let y = if pos.y < r.center().y { r.top() } else { r.bottom() };
                                    ui.painter().hline(r.x_range(), y, egui::Stroke::new(2.0, egui::Color32::WHITE));
                                }

                                // on drop, reorder
                                if let Some(src_idx) = drag_resp.dnd_release_payload::<usize>() {
                                    let r   = drag_resp.rect;
                                    let pos = ui.input(|i| i.pointer.interact_pos()).unwrap();
                                    let to_index = if pos.y < r.center().y { row } else { row + 1 };
                                    let item = file_paths.remove(*src_idx);
                                    let to = to_index.min(file_paths.len());
                                    file_paths.insert(to, item);
                                }
                            }
                        });
                    });
            });

            // 3) Add / remove buttons
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui.button(format!("➕ Add {}…", button_name)).clicked() {
                    let mut dlg = rfd::FileDialog::new();
                    if let Some(t) = file_dialog_title {
                        dlg = dlg.set_title(t);
                    }
                    if let Some(filters) = file_dialog_filter {
                        dlg = dlg.add_filter(file_type_name, filters);
                    }
                    dlg = dlg.add_filter("All Files", &["*"]);
                    if let Some(paths) = dlg.pick_files() {
                        file_paths.extend(paths);
                    }
                }
                if !file_paths.is_empty()
                    && ui.button(format!("❌ Remove Last {}", button_name)).clicked()
                {
                    file_paths.pop();
                }
            });
        });
    });
}
