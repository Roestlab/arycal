use arycal_cli::input::Input;
use eframe::egui::{Ui, ComboBox, Slider, TextEdit, ScrollArea, Sense};
use eframe::egui::{Id, Pos2, Rect};
use egui::{StrokeKind, UiBuilder};
use egui::{
    Frame,  Stroke, Vec2,
    epaint::Color32,
};
use std::path::PathBuf;
use arycal_common::config::{XicFileType, FeaturesFileType};
use rfd::FileDialog;




/// Draw the shared XIC, Features, and Filters config
pub fn draw_shared(ui: &mut Ui, config: &mut Input) {
    ui.heading("File Settings");
    ui.add_space(4.0);

    // XIC config
    egui::CollapsingHeader::new("XIC Files")
    .default_open(true)
    .show(ui, |ui| {
        // if the user hasn’t picked a type yet, but they have dropped in files,
        // try to auto‐detect from the first file’s extension:
        if config.xic.file_type.is_none() {
            if let Some(first_path) = config.xic.file_paths.get(0) {
                if let Some(ext) = first_path.extension().and_then(|e| e.to_str()) {
                    match ext.to_lowercase().as_str() {
                        "sqmass"  => config.xic.file_type = Some(XicFileType::SqMass),
                        "parquet" => config.xic.file_type = Some(XicFileType::parquet),
                        _         => config.xic.file_type = Some(XicFileType::Unknown),
                    }
                }
            }
        }

        // file type combo
        ui.horizontal(|ui| {
            ui.label("XIC Type:");
            let current = config.xic.file_type
                .as_ref()
                .map(|v| match v {
                    XicFileType::SqMass   => "sqMass",
                    XicFileType::parquet  => "parquet",
                    XicFileType::Unknown  => "unknown",
                })
                .unwrap_or("None")
                .to_string();

            ComboBox::from_id_salt("xic_type")
                .selected_text(current)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut config.xic.file_type, Some(XicFileType::SqMass),  "sqMass");
                    ui.selectable_value(&mut config.xic.file_type, Some(XicFileType::parquet), "parquet");
                });
        });

        // file paths drag & drop
        edit_file_paths(ui, &mut config.xic.file_paths, "XIC", Some("Select XIC Files"), Some(&vec!["sqMass", "parquet"]));
    });

    // Features config
    egui::CollapsingHeader::new("Feature Files")
    .default_open(true)
    .show(ui, |ui| {
        // Auto‐detect from the first path’s extension if not set yet
        if config.features.file_type.is_none() {
            if let Some(first) = config.features.file_paths.get(0) {
                if let Some(ext) = first.extension().and_then(|e| e.to_str()) {
                    config.features.file_type = Some(match ext.to_lowercase().as_str() {
                        "osw" => FeaturesFileType::OSW,
                        _     => FeaturesFileType::Unknown,
                    });
                }
            }
        }

        // Now draw the combo box, letting the user override if they like
        ui.horizontal(|ui| {
            ui.label("Features Type:");
            let current = config.features.file_type
                .as_ref()
                .map(|v| match v {
                    FeaturesFileType::OSW     => "osw",
                    FeaturesFileType::Unknown => "unknown",
                })
                .unwrap_or("None")
                .to_string();

            ComboBox::from_id_salt("features_type")
                .selected_text(current)
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut config.features.file_type,
                        Some(FeaturesFileType::OSW),
                        "osw",
                    );
                    ui.selectable_value(
                        &mut config.features.file_type,
                        Some(FeaturesFileType::Unknown),
                        "unknown",
                    );
                });
        });
        // file paths drag & drop
        edit_file_paths(ui, &mut config.features.file_paths, "Feature", Some("Select Feature Files"), Some(&["osw"]));
    });

    // Filters config
    ui.collapsing("Filters Configuration", |ui| {
        ui.checkbox(&mut config.filters.decoy, "Exclude Decoys");
        ui.checkbox(&mut config.filters.include_identifying_transitions.unwrap_or_default(), "Include Identifying Transitions");
        ui.checkbox(&mut config.xic.include_precursor, "Include Precursor");
        ui.add(Slider::new(&mut config.xic.num_isotopes, 1..=10).text("Number of Isotopes"));
        // Add numeric input for config.filters.max_score_ms2_qvalue, range 0 to 1
        ui.horizontal(|ui| {
            ui.label("Max Score MS2 Q-value:");
            let mut max_score = config.filters.max_score_ms2_qvalue.unwrap_or(0.0).to_string();
            if ui.add(TextEdit::singleline(&mut max_score)).changed() {
                if let Ok(val) = max_score.parse::<f64>() {
                    config.filters.max_score_ms2_qvalue = Some(val);
                }
            }
        });
    });
}


/// Draw additional settings for OpenSwath tab
pub fn draw_open_swath(ui: &mut Ui, config: &mut Input) {
    ui.heading("OpenSwath Workflow Settings");
    ui.horizontal(|ui| {
        ui.label("Binary Path:");
        // let mut path = config.open_swath.binary_path.clone().unwrap_or_default();
        let mut path = "";
        if ui.add(TextEdit::singleline(&mut path)).changed() {
            // config.open_swath.binary_path = Some(path.clone());
            todo!("Implement binary path setting");
        }
    });
}

/// Draw additional settings for Alignment tab
pub fn draw_alignment(ui: &mut Ui, config: &mut Input) {
    ui.heading("Alignment Settings");
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
            ui.selectable_value(&mut config.alignment.decoy_peak_mapping_method, "SHUFFLE".into(), "shuffle");
            ui.selectable_value(&mut config.alignment.decoy_peak_mapping_method, "RANDOM_REGION".into(), "random_region");
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
}

/// Shared file paths drag & drop widget
pub fn edit_file_paths(
    ui: &mut egui::Ui,
    file_paths: &mut Vec<PathBuf>,
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
                if ui.button(format!("➕ Add {}…", file_type_name)).clicked() {
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
                    && ui.button(format!("❌ Remove Last {}", file_type_name)).clicked()
                {
                    file_paths.pop();
                }
            });
        });
    });
}
