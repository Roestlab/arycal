use std::collections::BTreeMap;
use arycal_cli::input::Input;
use eframe::egui::{Ui, ComboBox, Slider, TextEdit, ScrollArea, Sense};
use arycal_common::config::{FeaturesFileType, PlotMode, VisualizationConfig, XicFileType};
use egui::{FontId, RichText};
use rfd::FileDialog;

use crate::tabs::visualization_tab::VisualizationState;


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
    
    ui.add_space(25.0);
    ui.separator();

    // Add plot settings
    add_plot_settings(ui, viz_cfg);
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
    ui.horizontal(|ui| {
        ui.checkbox(&mut viz_cfg.show_title, "Show title");
        ui.checkbox(&mut viz_cfg.show_axis_labels, "Show axis labels");
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