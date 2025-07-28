use std::collections::BTreeMap;
use arycal_cli::input::Input;
use arycal_cloudpath::osw::PrecursorPeakBoundaries;
use eframe::egui::{Ui, ComboBox, Slider, TextEdit, ScrollArea, Sense};
use arycal_common::config::{FeaturesFileType, PlotMode, VisualizationConfig, XicFileType};
use egui::{FontId, RichText};
use rfd::FileDialog;

use crate::tabs::visualization_tab::VisualizationState;

use super::config_panel::edit_file_paths;

pub fn draw_visualization_file_settings(ui: &mut Ui, config: &mut Input) {
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
                        "parquet" => config.xic.file_type = Some(XicFileType::Parquet),
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
                    XicFileType::Parquet  => "parquet",
                    XicFileType::Unknown  => "unknown",
                })
                .unwrap_or("None")
                .to_string();

            ComboBox::from_id_salt("xic_type")
                .selected_text(current)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut config.xic.file_type, Some(XicFileType::SqMass),  "sqMass");
                    ui.selectable_value(&mut config.xic.file_type, Some(XicFileType::Parquet), "parquet");
                });
        });

        // file paths drag & drop
        edit_file_paths(ui, &mut config.xic.file_paths, "XIC", "XIC files", Some("Select XIC Files"), Some(&vec!["sqMass", "parquet"]));
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
        edit_file_paths(ui, &mut config.features.file_paths, "Feature","Feature files", Some("Select Feature Files"), Some(&["osw"]));
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

pub fn draw_visualization_settings(
    ui: &mut Ui,
    vis: &mut VisualizationState,
    config: &mut Input,
) {
    ui.heading("Peptide Selection");

    // 1) If precursor_table is not loaded, show spinner and exit early
    let table: &BTreeMap<String, Vec<u8>> = match &vis.precursor_table {
        Some(t) => t,
        None => {
            ui.label("⟳ Loading peptide → charge table…");
            return;
        }
    };

    // 2) Grab or initialize the per-UI config bucket
    let viz_cfg: &mut VisualizationConfig =
        config.visualization.get_or_insert_with(VisualizationConfig::default);

    // 0) Add a filter box before the dropdown:
    ui.horizontal(|ui| {
        ui.label("Filter peptides:");
        ui.add(
            TextEdit::singleline(&mut viz_cfg.peptide_filter)
                .hint_text("Type to search…")
                .desired_width(200.0),
        );
    });

    let filter_lower = viz_cfg.peptide_filter.to_lowercase();

    // 1) Peptide dropdown, but only show matching keys:
    ui.horizontal(|ui| {
        ui.label("Peptide:");
        ComboBox::from_id_salt("viz_peptide")
            .wrap()
            .selected_text(
                viz_cfg
                    .selected_peptide
                    .as_deref()
                    .unwrap_or("— select peptide —"),
            )
            .show_ui(ui, |ui| {
                for pep in table.keys().filter(|p| {
                    filter_lower.is_empty() ||
                    p.to_lowercase().contains(&filter_lower)
                }) {
                    ui.selectable_value(
                        &mut viz_cfg.selected_peptide,
                        Some(pep.clone()),
                        pep,
                    );
                }                
            });
    });

    if let Some(peptide) = &viz_cfg.selected_peptide {
        if let Some(charges) = table.get(peptide) {
            // If `selected_charge` is None *or* no longer in the current list, pick the first one.
            let default_charge = charges[0] as usize;
            if viz_cfg
                .selected_charge
                .map(|c| !charges.iter().any(|&x| x as usize == c))
                .unwrap_or(true)
            {
                viz_cfg.selected_charge = Some(default_charge);
            }
    
            ui.horizontal(|ui| {
                ui.label("Charge:");
                ComboBox::from_id_salt("viz_charge")
                    .selected_text(viz_cfg.selected_charge.unwrap().to_string())
                    .show_ui(ui, |ui| {
                        for &c in charges {
                            ui.selectable_value(
                                &mut viz_cfg.selected_charge,
                                Some(c as usize),
                                c.to_string(),
                            );
                        }
                    });
            });
        }
    }

    if let Some(boundaries) = &vis.precursor_boundaries {
        ui.add_space(25.0);
        ui.separator();
        ui.heading("Feature Boundaries Selection");
    
        // Group by run
        let mut runs: std::collections::BTreeMap<String, Vec<&PrecursorPeakBoundaries>> =
            std::collections::BTreeMap::new();
        for b in boundaries {
            runs.entry(b.run_filename.clone()).or_default().push(b);
        }
    
        // --- Top-level buttons for all runs ---
        ui.horizontal(|ui| {
            if ui.button("Select All").clicked() {
                for (run_name, features) in &runs {
                    for b in features.clone() {
                        vis.selected_boundaries.insert((run_name.clone(), b.feature_id));
                    }
                }
            }
            if ui.button("Deselect All").clicked() {
                for (run_name, features) in &runs {
                    for b in features.clone() {
                        vis.selected_boundaries.remove(&(run_name.clone(), b.feature_id));
                    }
                }
            }
            if ui.button("Select First").clicked() {
                vis.selected_boundaries.clear();
                for (run_name, features) in &runs {
                    if let Some(first) = features.iter().min_by_key(|b| b.sorted_feature_id) {
                        vis.selected_boundaries.insert((run_name.clone(), first.feature_id));
                    }
                }
            }
        });
    
        // --- Initialize selected_boundaries ONCE if empty AND not yet initialized ---
        if !vis.boundaries_initialized && vis.selected_boundaries.is_empty() {
            for (run_name, features) in &runs {
                if let Some(first) = features.iter().min_by_key(|b| b.sorted_feature_id) {
                    vis.selected_boundaries.insert((run_name.clone(), first.feature_id));
                }
            }
            vis.boundaries_initialized = true; // mark as initialized
        }

    
        // --- Rendering loop: iterate by reference so we don't move runs ---
        egui::ScrollArea::vertical().max_height(600.0).show(ui, |ui| {
            for (run_name, features) in &mut runs {
                // Sort features for display by sorted_feature_id
                features.sort_by_key(|b| b.sorted_feature_id);
    
                egui::CollapsingHeader::new(run_name.clone())
                    .default_open(false)
                    .show(ui, |ui| {
                        // "Select All / Deselect All" for convenience (per run)
                        ui.horizontal(|ui| {
                            if ui.button("Select All").clicked() {
                                for b in features.clone() {
                                    vis.selected_boundaries.insert((run_name.clone(), b.feature_id));
                                }
                            }
                            if ui.button("Deselect All").clicked() {
                                for b in features.clone() {
                                    vis.selected_boundaries.remove(&(run_name.clone(), b.feature_id));
                                }
                            }
                        });
    
                        egui::Grid::new(format!("boundary_table_{}", run_name))
                            .striped(true)
                            .show(ui, |ui| {
                                // Header
                                ui.label(""); // checkbox
                                ui.label("AlignGrp ID");
                                ui.label("Feature ID");
                                ui.label("Type");
                                ui.label("Rank");
                                ui.label("PG Q");
                                ui.label("Align Q");
                                ui.end_row();
    
                                for b in features.clone() {
                                    let key = (run_name.clone(), b.feature_id);
                                    let mut selected = vis.selected_boundaries.contains(&key);
    
                                    if ui.checkbox(&mut selected, "").changed() {
                                        if selected {
                                            vis.selected_boundaries.insert(key.clone());
                                        } else {
                                            vis.selected_boundaries.remove(&key);
                                        }
                                    }
    
                                    ui.label(b.alignment_group_id.map(|x| x.to_string()).unwrap_or("-".into()));
                                    ui.label(b.feature_id.to_string());
                                    ui.label(b.feature_type.clone());
                                    ui.label(b.peakgroup_rank.to_string());
                                    ui.label(b.peakgroup_qvalue.map(|x| format!("{:.3}", x)).unwrap_or("-".into()));
                                    ui.label(b.alignment_qvalue.map(|x| format!("{:.3}", x)).unwrap_or("-".into()));
                                    ui.end_row();
                                }
                            });
                    });
            }
        });
    }
    

    ui.add_space(25.0);
    ui.separator();
    // Add plot settings
    add_plot_settings(ui, viz_cfg);

    // Keyboard shortcuts
    ui.separator();
    ui.collapsing("Key Bindings", |ui| {
        egui::Grid::new("key_binding_instructions")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let egui::InputOptions {
                    zoom_modifier,
                    horizontal_scroll_modifier,
                    vertical_scroll_modifier,
                    ..
                } = egui::InputOptions::default();

                ui.label("Pan");
                ui.label("Left-drag");
                ui.end_row();

                ui.label("Horizontal pan");
                ui.label(format!(
                    "{} + Scroll",
                    ui.ctx().format_modifiers(horizontal_scroll_modifier)
                ));
                ui.end_row();

                ui.label("Zoom");
                ui.label(format!(
                    "{} + Scroll",
                    ui.ctx().format_modifiers(zoom_modifier)
                ));
                ui.end_row();

                ui.label("Horizontal zoom");
                ui.label(format!(
                    "{} + Scroll",
                    ui.ctx()
                        .format_modifiers(zoom_modifier | horizontal_scroll_modifier)
                ));
                ui.end_row();

                ui.label("Vertical zoom");
                ui.label(format!(
                    "{} + Scroll",
                    ui.ctx()
                        .format_modifiers(zoom_modifier | vertical_scroll_modifier)
                ));
                ui.end_row();

                ui.label("Box zoom");
                ui.label("Right-drag");
                ui.end_row();

                ui.label("Reset");
                ui.label("Double-click");
                ui.end_row();
            });
    });
}

fn add_plot_settings(
    ui: &mut Ui,
    viz_cfg: &mut VisualizationConfig,
) {
    ui.heading("Plot Settings");

    ui.add_space(5.0);

    // Smoothing
    ui.horizontal(|ui| {
        ui.label("Smoothing:");
        ui.radio_value(
            &mut viz_cfg.smoothing_enabled,
            true,
            "On",
        );
        ui.radio_value(
            &mut viz_cfg.smoothing_enabled,
            false,
            "Off",
        );
    });

    if viz_cfg.smoothing_enabled {
        ui.add_space(5.0);
        ui.label("Savitzky-Golay Smoothing Params");

        ui.horizontal(|ui| {
        ui.label("window:");

        ui.add(
            Slider::new(&mut viz_cfg.sgolay_window, 1..=51)
                .step_by(2.0)
                .text("points"),
        );
        });
        ui.horizontal(|ui| {
            ui.label("order:");
            ui.add(Slider::new(&mut viz_cfg.sgolay_order, 1..=21).text("degree"));
        });

        // Ensure that the polynomial order is always less than the window length
        if viz_cfg.sgolay_order >= viz_cfg.sgolay_window {
            viz_cfg.sgolay_order = viz_cfg.sgolay_window - 1;
        }
    }

    ui.horizontal(|ui| {
        ui.checkbox(&mut viz_cfg.link_axis_x, "Link X-axis");
        ui.checkbox(&mut viz_cfg.link_axis_y, "Link Y-axis");
        ui.checkbox(&mut viz_cfg.link_cursor, "Link cursor");
    });
    
    
    ui.add_space(5.0);
    ui.separator();

    ui.label("Plot Aesthetics");
    ui.horizontal(|ui| {
        ui.checkbox(&mut viz_cfg.show_background, "Show background");
        ui.checkbox(&mut viz_cfg.show_grid, "Show grid");
        ui.checkbox(&mut viz_cfg.show_legend, "Show legend");
    });

    ui.add_space(5.0);
    ui.separator();

    ui.label("Plot View Layout");
    ui.horizontal(|ui| {
        ui.label("Plot layout:");
        ui.selectable_value(
            &mut viz_cfg.plot_mode,
            PlotMode::Floating,
            "Floating"
        );
        ui.selectable_value(
            &mut viz_cfg.plot_mode,
            PlotMode::EmbeddedGrid,
            "Embedded"
        );
    });
    if viz_cfg.plot_mode == PlotMode::EmbeddedGrid {
        ui.horizontal(|ui| {
            ui.add(
                egui::Slider::new(&mut viz_cfg.grid_rows, 1..=4)
                    .text("Rows")
            );
            ui.add(
                egui::Slider::new(&mut viz_cfg.grid_cols, 1..=4)
                    .text("Cols")
            );
        });
    }

    if viz_cfg.plot_mode == PlotMode::Floating {
        if ui.button("Reorganize Windows").clicked() {
            // clear saved positions for everything (panels, windows, etc)
            ui.ctx().memory_mut(|mem| mem.reset_areas());
        }
    }

    
    ui.add_space(25.0);
}