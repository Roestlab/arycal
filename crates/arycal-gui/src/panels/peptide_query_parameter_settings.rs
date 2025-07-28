use std::{env, path::{Path, PathBuf}};

use arycal_cli::input::Input;
use arycal_cloudpath::util::find_executable;
use arycal_common::config::PQPConfig;
use egui::{Checkbox, ComboBox, DragValue, TextEdit, Ui};
use rfd::FileDialog;

use super::config_panel::edit_file_paths;

pub fn draw_pqp_file_settings(ui: &mut Ui, cfg: &mut PQPConfig) {
    // Files
    egui::CollapsingHeader::new("Transition Files")
    .default_open(true)
    .show(ui, |ui| {
        // file paths drag & drop
        edit_file_paths(ui, &mut cfg.file_paths, "Transition list", "Transition Files: tsv, pqp, traML", Some("Select Transition Files"), Some(&vec!["tsv", "pqp", "traML"]));
    });

    // Check if cfg.file_paths is not empty
    // if not empty, and a vector of a single file, then create a single file dialog to setting output file
    // if not empty and a vector of files, then create a file dialog to select output path
    if !cfg.file_paths.is_empty() {
        if cfg.file_paths.len() == 1 {
            // Single file, set output file
            ui.horizontal(|ui| {
                ui.label("Output File:");
                let mut output_path = cfg.output_path.to_string_lossy().into_owned();
                if ui.add(TextEdit::singleline(&mut output_path)).changed() {
                    cfg.output_path = PathBuf::from(output_path);
                }
                if ui.button("…").on_hover_text("Browse to save output file").clicked() {
                    if let Some(file) = FileDialog::new()
                        .set_title("Select Output File")
                        .save_file()
                    {
                        cfg.output_path = file;
                    }
                }
            });
        } else {
            // output type
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Output type:"))
                  .on_hover_text("Force output format, otherwise inferred from extension.");
                ComboBox::from_id_salt("pqp_out_type")
                    .selected_text(&cfg.pqp_out_type)
                    .show_ui(ui, |ui| {
                        for &opt in &["tsv","pqp","TraML"] {
                            ui.selectable_value(&mut cfg.pqp_out_type, opt.to_string(), opt);
                        }
                    });
            });

            // Multiple files, set output directory
            ui.horizontal(|ui| {
                ui.label("Output Directory:");
                let mut output_path = cfg.output_path.to_string_lossy().into_owned();
                if ui.add(TextEdit::singleline(&mut output_path)).changed() {
                    cfg.output_path = PathBuf::from(output_path);
                }
                if ui.button("…").on_hover_text("Browse to save output file").clicked() {
                    if let Some(file) = FileDialog::new()
                        .set_title("Select Output File")
                        .pick_folder()
                    {
                        cfg.output_path = file;
                    }
                }
            });
        }
    }
}

/// Draw additional settings for PQP tab
pub fn draw_pqp_generation(ui: &mut Ui, config: &mut Input) {
    ui.heading("Peptide Query Parameter Generation Settings");

    // Ensure we have a PQPConfig
    let pqp_cfg = config
        .pqp
        .get_or_insert_with(PQPConfig::default);

    // Add single radio toggle to switch between main PQP generation and iRT PQP generation for pqp_cfg.main_mode
    ui.horizontal(|ui| {
        ui.label("Mode:");
        ComboBox::from_id_salt("pqp_mode")
            .selected_text(if pqp_cfg.main_mode { "Main PQP" } else { "iRT PQP" })
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut pqp_cfg.main_mode, true, "Main PQP");
                ui.selectable_value(&mut pqp_cfg.main_mode, false, "iRT PQP");
            });
    });

    // Add Checkboxes to enable/disable each tool
    ui.horizontal(|ui| {
        if pqp_cfg.main_mode {
            ui.add(Checkbox::new(&mut pqp_cfg.tfc_enabled, "TargetedFileConverter"))
            .on_hover_text("Enable/disable TargetedFileConverter tool.");
        }
        ui.add(Checkbox::new(&mut pqp_cfg.osg_enabled, "OpenSwathAssayGenerator"))
          .on_hover_text("Enable/disable OpenSwathAssayGenerator tool.");
        if pqp_cfg.main_mode {
            ui.add(Checkbox::new(&mut pqp_cfg.odg_enabled, "OpenSwathDecoyGenerator"))
              .on_hover_text("Enable/disable OpenSwathDecoyGenerator tool.");
        } else {
            // Switch ODG to False if not in main mode
            pqp_cfg.odg_enabled = false;

            ui.add(Checkbox::new(&mut pqp_cfg.irt_reduce_enabled, "EasyPQP Reduce"))
              .on_hover_text("Enable/disable to create reduced linear iRT PQP from full iRT PQP.");
        }
    });

    // --- TargetedFileConverter panel ---
    if pqp_cfg.tfc_enabled && pqp_cfg.main_mode {
        egui::CollapsingHeader::new("TargetedFileConverter")
        .default_open(true)
        .show(ui, |ui| {
            // — binary path + auto-detect button —
            if pqp_cfg.tfc_binary_path.as_os_str().is_empty() {
                if let Some(p) = find_executable("TargetedFileConverter", None) {
                    pqp_cfg.tfc_binary_path = p;
                }
            }
    
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("TargetedFileConverter binary:")).on_hover_text(
                    "Path to the TargetedFileConverter binary. If not set, it will try to find it in your PATH."
                );
    
                // Show & edit it as a String
                let mut path_str = pqp_cfg.tfc_binary_path.to_string_lossy().into_owned();
                if ui.add(TextEdit::singleline(&mut path_str)).changed() {
                    pqp_cfg.tfc_binary_path = PathBuf::from(path_str.clone());
                }
    
                // File dialog button
                if ui.button("…").on_hover_text("Browse for binary").clicked() {
                    if let Some(file) = rfd::FileDialog::new()
                        .set_title("Select TargetedFileConverter binary")
                        // start in current dir or last used
                        .set_directory(env::current_dir().unwrap_or_default())
                        .pick_file() 
                    {
                        pqp_cfg.tfc_binary_path = file;
                    }
                }
            });
    
            ui.separator();
    
            // // input
            // ui.horizontal(|ui| {
            //     ui.add(egui::Label::new("Input transition list:"))
            //       .on_hover_text("Tab-separated transition list (.tsv, .mrm, .pqp, .TraML).");
            //     let mut s = pqp_cfg.tfc_input.to_string_lossy().to_string();
            //     if ui.add(TextEdit::singleline(&mut s)).changed() {
            //         pqp_cfg.tfc_input = PathBuf::from(&s);
            //     }
            //     if ui.button("…").on_hover_text("Browse input file").clicked() {
            //         if let Some(f) = FileDialog::new()
            //             .add_filter("Transition list", &["tsv","mrm","pqp","TraML"])
            //             .pick_file()
            //         {
            //             pqp_cfg.tfc_input = f;
            //         }
            //     }
            // });
            // // input type
            // ui.horizontal(|ui| {
            //     ui.add(egui::Label::new("Input type:"))
            //       .on_hover_text("Force interpretation, otherwise inferred from extension.");
            //     ComboBox::from_id_salt("tfc_in_type")
            //         .selected_text(&pqp_cfg.tfc_input_type)
            //         .show_ui(ui, |ui| {
            //             for &opt in &["tsv","mrm","pqp","TraML"] {
            //                 ui.selectable_value(&mut pqp_cfg.tfc_input_type, opt.to_string(), opt);
            //             }
            //         });
            // });

            // // output
            // ui.horizontal(|ui| {
            //     ui.add(egui::Label::new("Output file:"))
            //       .on_hover_text("Where to write the converted library.");
            //     let mut s = pqp_cfg.tfc_output.to_string_lossy().to_string();
            //     if ui.add(TextEdit::singleline(&mut s)).changed() {
            //         pqp_cfg.tfc_output = PathBuf::from(&s);
            //     }
            //     if ui.button("…").on_hover_text("Browse output file").clicked() {
            //         if let Some(f) = FileDialog::new()
            //             .add_filter("TraML/TSV/PQP", &["TraML","tsv","pqp"])
            //             .save_file()
            //         {
            //             pqp_cfg.tfc_output = f;
            //         }
            //     }
            // });
            // // output type
            // ui.horizontal(|ui| {
            //     ui.add(egui::Label::new("Output type:"))
            //       .on_hover_text("Force output format, otherwise inferred from extension.");
            //     ComboBox::from_id_salt("tfc_out_type")
            //         .selected_text(&pqp_cfg.tfc_output_type)
            //         .show_ui(ui, |ui| {
            //             for &opt in &["tsv","pqp","TraML"] {
            //                 ui.selectable_value(&mut pqp_cfg.tfc_output_type, opt.to_string(), opt);
            //             }
            //         });
            // });

            // threads & legacy flag
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Threads:"))
                  .on_hover_text("Number of threads to use.");
                ui.add(DragValue::new(&mut pqp_cfg.tfc_threads).speed(1.0));
                ui.add(Checkbox::new(&mut pqp_cfg.tfc_legacy_traml_id, "Legacy TraML IDs"))
                  .on_hover_text("Use old-style TraML identifiers (pqp→TraML).");
            });
            // advanced
            ui.label("Additional flags:");
            ui.add(
                TextEdit::multiline(&mut pqp_cfg.tfc_advanced)
                    .desired_rows(2)
                    .lock_focus(true)
                    .hint_text("e.g. -write_ini my.ini"),
            );
        });
    }
    
    // --- OpenSwathAssayGenerator panel ---
    if pqp_cfg.osg_enabled {
        egui::CollapsingHeader::new("OpenSwathAssayGenerator")
        .default_open(true)
        .show(ui, |ui| {

            // — binary path + auto-detect button —
            if pqp_cfg.osg_binary_path.as_os_str().is_empty() {
                if let Some(p) = find_executable("OpenSwathAssayGenerator", None) {
                    pqp_cfg.osg_binary_path = p;
                }
            }

            ui.horizontal(|ui| {
                ui.add(egui::Label::new("OpenSwathAssayGenerator binary:"))
                .on_hover_text("Path to the OpenSwathAssayGenerator binary. If not set, it will try to find it in your PATH.");
                // Show & edit it as a String
                let mut path_str = pqp_cfg.osg_binary_path.to_string_lossy().into_owned();
                if ui.add(TextEdit::singleline(&mut path_str)).changed() {
                    pqp_cfg.osg_binary_path = PathBuf::from(path_str.clone());
                }
                // File dialog button
                if ui.button("…").on_hover_text("Browse for binary").clicked() {
                    if let Some(file) = rfd::FileDialog::new()
                        .set_title("Select OpenSwathAssayGenerator binary")
                        // start in current dir or last used
                        .set_directory(env::current_dir().unwrap_or_default())
                        .pick_file() 
                    {
                        pqp_cfg.osg_binary_path = file;
                    }
                }
            });

            ui.separator();

            // // input TraML
            // ui.horizontal(|ui| {
            //     ui.add(egui::Label::new("Input:"))
            //     .on_hover_text("Input transition list.");
            //     let mut s = pqp_cfg.osg_input.to_string_lossy().to_string();
            //     if ui.add(TextEdit::singleline(&mut s)).changed() {
            //         pqp_cfg.osg_input = PathBuf::from(&s);
            //     }
            //     if ui.button("…").on_hover_text("Browse TraML").clicked() {
            //         if let Some(f) = FileDialog::new()
            //             .add_filter("TraML", &["TraML"])
            //             .pick_file()
            //         {
            //             pqp_cfg.osg_input = f;
            //         }
            //     }
            // });

            // min/max transitions
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Min transitions:"))
                .on_hover_text("Minimum transitions per assay.");
                ui.add(DragValue::new(&mut pqp_cfg.osg_min_transitions).speed(1.0));
                ui.add(egui::Label::new("Max transitions:"))
                .on_hover_text("Maximum transitions per assay.");
                ui.add(DragValue::new(&mut pqp_cfg.osg_max_transitions).speed(1.0));
            });
            // fragment types & charges
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Allowed frags:"))
                .on_hover_text("Comma-separated, e.g. b,y.");
                ui.text_edit_singleline(&mut pqp_cfg.osg_allowed_fragment_types);
                ui.add(egui::Label::new("Charges:"))
                .on_hover_text("Comma-separated, e.g. 1,2.");
                ui.text_edit_singleline(&mut pqp_cfg.osg_allowed_fragment_charges);
            });

            if pqp_cfg.main_mode {
                // losses & IPF
                ui.horizontal(|ui| {
                    ui.add(Checkbox::new(&mut pqp_cfg.osg_enable_detection_specific_losses, "Spec. losses"))
                    .on_hover_text("Allow specific neutral losses for detection ions.");
                    ui.add(Checkbox::new(&mut pqp_cfg.osg_enable_detection_unspecific_losses, "Unspecific losses"))
                    .on_hover_text("Allow unspecific neutral losses (H2O, NH3, etc).");
                    ui.add(Checkbox::new(&mut pqp_cfg.osg_enable_ipf, "Enable IPF"))
                    .on_hover_text("Generate theoretical ID transitions for IPF scoring.");
                });
                if pqp_cfg.osg_enable_ipf {
                    // unimod file
                    ui.horizontal(|ui| {
                        ui.add(egui::Label::new("Unimod file:"))
                        .on_hover_text("Restricted unimod XML for PTM localization.");
                        let mut s = pqp_cfg.osg_unimod_file.to_string_lossy().to_string();
                        if ui.add(TextEdit::singleline(&mut s)).changed() {
                            pqp_cfg.osg_unimod_file = PathBuf::from(&s);
                        }
                        if ui.button("…").on_hover_text("Browse unimod XML").clicked() {
                            if let Some(f) = FileDialog::new()
                                .add_filter("XML", &["xml"])
                                .pick_file()
                            {
                                pqp_cfg.osg_unimod_file = f;
                            }
                        }
                    });
                }                
            }
            
            // advanced
            ui.label("Additional flags:");
            ui.add(
                TextEdit::multiline(&mut pqp_cfg.osg_advanced)
                    .desired_rows(2)
                    .lock_focus(true)
                    .hint_text("-max_num_alternative_localizations 10000"),
            );
        });
    }

    // --- OpenSwathDecoyGenerator panel ---
    if pqp_cfg.odg_enabled && pqp_cfg.main_mode {
        egui::CollapsingHeader::new("OpenSwathDecoyGenerator")
        .default_open(true)
        .show(ui, |ui| {
            // — binary path + auto-detect button —
            if pqp_cfg.odg_binary_path.as_os_str().is_empty() {
                if let Some(p) = find_executable("OpenSwathDecoyGenerator", None) {
                    pqp_cfg.odg_binary_path = p;
                }
            }

            ui.horizontal(|ui| {
                ui.add(egui::Label::new("OpenSwathDecoyGenerator binary:"))
                .on_hover_text("Path to the OpenSwathDecoyGenerator binary. If not set, it will try to find it in your PATH.");
                // Show & edit it as a String
                let mut path_str = pqp_cfg.odg_binary_path.to_string_lossy().into_owned();
                if ui.add(TextEdit::singleline(&mut path_str)).changed() {
                    pqp_cfg.odg_binary_path = PathBuf::from(path_str.clone());
                }
                // File dialog button
                if ui.button("…").on_hover_text("Browse for binary").clicked() {
                    if let Some(file) = rfd::FileDialog::new()
                        .set_title("Select OpenSwathDecoyGenerator binary")
                        // start in current dir or last used
                        .set_directory(env::current_dir().unwrap_or_default())
                        .pick_file() 
                    {
                        pqp_cfg.odg_binary_path = file;
                    }
                }
            });

            ui.separator();

            // // input TraML
            // ui.horizontal(|ui| {
            //     ui.add(egui::Label::new("Input:"))
            //     .on_hover_text("Transition list of target assays.");
            //     let mut s = pqp_cfg.odg_input.to_string_lossy().to_string();
            //     if ui.add(TextEdit::singleline(&mut s)).changed() {
            //         pqp_cfg.odg_input = PathBuf::from(&s);
            //     }
            //     if ui.button("…").on_hover_text("Browse TraML").clicked() {
            //         if let Some(f) = FileDialog::new().add_filter("TraML", &["TraML"]).pick_file() {
            //             pqp_cfg.odg_input = f;
            //         }
            //     }
            // });

            // method & tag
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Method:"))
                .on_hover_text("reverse, pseudo-reverse, shuffle, shift");
                ComboBox::from_id_salt("odg_method")
                    .selected_text(&pqp_cfg.odg_method)
                    .show_ui(ui, |ui| {
                        for &opt in &["reverse", "pseudo-reverse", "shuffle", "shift"] {
                            ui.selectable_value(&mut pqp_cfg.odg_method, opt.to_string(), opt);
                        }
                    });
                ui.add(egui::Label::new("Decoy tag:"))
                .on_hover_text("Appended to peptide IDs for decoys.");
                ui.text_edit_singleline(&mut pqp_cfg.odg_decoy_tag);
            });
            // fractions
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Min fraction:"))
                .on_hover_text("Minimum decoy/target ratio.");
                ui.add(DragValue::new(&mut pqp_cfg.odg_min_decoy_fraction).speed(0.01));
                ui.add(egui::Label::new("Aim fraction:"))
                .on_hover_text("Target number of decoys per target.");
                ui.add(DragValue::new(&mut pqp_cfg.odg_aim_decoy_fraction).speed(0.01));
            });
            // shuffle & shift params
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Shuffle attempts:"))
                .on_hover_text("Max tries to reduce seq. identity.");
                ui.add(DragValue::new(&mut pqp_cfg.odg_shuffle_max_attempts).speed(1.0));
                ui.add(egui::Label::new("Shuffle identity threshold:"))
                .on_hover_text("Max allowed sequence identity.");
                ui.add(DragValue::new(&mut pqp_cfg.odg_shuffle_sequence_identity_threshold).speed(0.01));
            });
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Shift prec. Δm/z:"))
                .on_hover_text("Shift size in Thomson for precursor.");
                ui.add(DragValue::new(&mut pqp_cfg.odg_shift_precursor_mz_shift).speed(0.01));
                ui.add(egui::Label::new("Shift frag. Δm/z:"))
                .on_hover_text("Shift size in Thomson for fragment.");
                ui.add(DragValue::new(&mut pqp_cfg.odg_shift_product_mz_shift).speed(0.01));
            });
            // advanced
            ui.label("Additional flags:");
            ui.add(
                TextEdit::multiline(&mut pqp_cfg.odg_advanced)
                    .desired_rows(2)
                    .lock_focus(true)
                    .hint_text("-switchKR false"),
            );
        });
    }

    // --- EasyPQP Reduce panel ---
    if pqp_cfg.irt_reduce_enabled && !pqp_cfg.main_mode {
        egui::CollapsingHeader::new("EasyPQP Reduce iRT")
        .default_open(true)
        .show(ui, |ui| {
            // — binary path + auto-detect button —
            if pqp_cfg.easypqp_binary_path.as_os_str().is_empty() {
                if let Some(p) = find_executable("easypqp", None) {
                    pqp_cfg.easypqp_binary_path = p;
                }
            }

            ui.horizontal(|ui| {
                ui.add(egui::Label::new("EasyPQP binary:"))
                .on_hover_text("Path to the EasyPQP binary. If not set, it will try to find it in your PATH.");
                // Show & edit it as a String
                let mut path_str = pqp_cfg.easypqp_binary_path.to_string_lossy().into_owned();
                if ui.add(TextEdit::singleline(&mut path_str)).changed() {
                    pqp_cfg.easypqp_binary_path = PathBuf::from(path_str.clone());
                }
                // File dialog button
                if ui.button("…").on_hover_text("Browse for binary").clicked() {
                    if let Some(file) = rfd::FileDialog::new()
                        .set_title("Select EasyPQP binary")
                        // start in current dir or last used
                        .set_directory(env::current_dir().unwrap_or_default())
                        .pick_file() 
                    {
                        pqp_cfg.easypqp_binary_path = file;
                    }
                }
            });

            ui.separator();

            // numeric input for number of iRT bins and peptides to smaple
            ui.horizontal(|ui| {
                ui.add(egui::Label::new("Number of bins:"))
                  .on_hover_text("Number of bins to fill along gradient.");
                ui.add(DragValue::new(&mut pqp_cfg.irt_bins).speed(1.0));
                ui.add(egui::Label::new("Peptides to sample:"))
                  .on_hover_text("Number of peptides to sample.");
                ui.add(DragValue::new(&mut pqp_cfg.irt_num_peptides).speed(1.0));
            });

        });
    }
}