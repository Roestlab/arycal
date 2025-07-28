use std::path::{Path, PathBuf};

use arycal_cli::input::Input;
use arycal_cloudpath::util::find_executable;
use arycal_common::config::PyProphetConfig;
use egui::{Checkbox, ComboBox, DragValue, TextEdit, Ui};
use rfd::FileDialog;
use serde::de;

use super::config_panel::edit_file_paths;

#[derive(Clone, Debug, Default)]
enum PyPScoreLevel {
    MS1,
    MS2,
    #[default]
    MS1MS2,
    Transition,
    Alignment,
}

impl PyPScoreLevel {
    fn as_str(&self) -> &'static str {
        match self {
            PyPScoreLevel::MS1 => "ms1",
            PyPScoreLevel::MS2 => "ms2",
            PyPScoreLevel::MS1MS2 => "ms1ms2",
            PyPScoreLevel::Transition => "transition",
            PyPScoreLevel::Alignment => "alignment",
        }
    }
}

#[derive(Clone, Debug, Default)]
enum PyPScoreClassifier {
    LDA,

    #[default]
    XGBoost,
    SVM
}

impl PyPScoreClassifier {
    fn as_str(&self) -> &'static str {
        match self {
            PyPScoreClassifier::LDA => "LDA",
            PyPScoreClassifier::XGBoost => "XGBoost",
            PyPScoreClassifier::SVM => "SVM",
        }
    }
    
}

/// Generic struct for main scoring settings
struct ScoreSection<'a> {
    /// Label to show in the header (“MS1”, “MS2”, …)
    label: &'static str,
    /// Which combo-box id salt to use (must be unique per section)
    id_salt: &'static str,
    integrate_ms1: &'a mut bool,
    classifier: &'a mut String,
    ss_num_iter: &'a mut usize,
    xeval_num_iter: &'a mut usize,
    ss_initial_fdr: &'a mut f32,
    ss_iteration_fdr: &'a mut f32,
    advanced_params: &'a mut String,
}

/// Generic struct for main inference settings
struct PyPInferSection<'a> {
    label: &'static str,
    id_salt: &'static str,
    level: &'static str,
    global: &'a mut bool,
    experiment_wide: &'a mut bool,
    run_specific: &'a mut bool,
    ipf_ms1_scoring: &'a mut bool,
    ipf_ms2_scoring: &'a mut bool,
    ipf_max_precursor_pep: &'a mut f32,
    ipf_max_peakgroup_pep: &'a mut f32,
    ipf_max_precursor_peakgroup_pep: &'a mut f32,
    ipf_max_transition_pep: &'a mut f32,
    ipf_propgate_signal: &'a mut bool,
    ipf_max_alignment_pep: &'a mut f32,
    infer_peptidoform_advanced_params: &'a mut String,
}


pub fn draw_validation_file_settings(ui: &mut Ui, cfg: &mut PyProphetConfig) {
    // Files
    egui::CollapsingHeader::new("PyProphet Files")
    .default_open(true)
    .show(ui, |ui| {
        // file paths drag & drop
        edit_file_paths(ui, &mut cfg.file_paths, "PyProphet File", "PyProphet Files: osw, parquet, oswpq, oswpqd, tsv", Some("Select PyProphet File(s)"), Some(&vec!["osw", "parquet", "oswpq", "oswpqd", "tsv"]));
    });

    // // Check if cfg.file_paths is not empty
    // // if not empty, and a vector of a single file, then create a single file dialog to setting output file
    // // if not empty and a vector of files, then create a file dialog to select output path
    // if !cfg.file_paths.is_empty() {
    //     if cfg.file_paths.len() == 1 {
    //         // Single file, set output file
    //         ui.horizontal(|ui| {
    //             ui.label("Output File:");
    //             let mut output_path = cfg.output_path.to_string_lossy().into_owned();
    //             if ui.add(TextEdit::singleline(&mut output_path)).changed() {
    //                 cfg.output_path = PathBuf::from(output_path);
    //             }
    //             if ui.button("…").on_hover_text("Browse to save output file").clicked() {
    //                 if let Some(file) = FileDialog::new()
    //                     .set_title("Select Output File")
    //                     .save_file()
    //                 {
    //                     cfg.output_path = file;
    //                 }
    //             }
    //         });
    //     } else {
    //         // if the user hasn’t picked a type yet, but they have dropped in files,
    //         // try to auto‐detect from the first file’s extension:
    //         if cfg.output_file_type.is_empty() {
    //             if let Some(first_file) = cfg.file_paths.first() {
    //                 if let Some(ext) = first_file.extension() {
    //                     // Convert OsStr to String, stripping leading dot
    //                     cfg.output_file_type = ext.to_string_lossy().trim_start_matches('.').to_string();
    //                 }
    //             }
    //         }

    //         // output type
    //         ui.horizontal(|ui| {
    //             ui.add(egui::Label::new("Output type:"))
    //               .on_hover_text("Force output format, otherwise inferred from extension.");
    //             ComboBox::from_id_salt("pqp_out_type")
    //                 .selected_text(&cfg.output_file_type)
    //                 .show_ui(ui, |ui| {
    //                     for &opt in &["osw", "parquet", "oswpq", "oswpqd", "tsv"] {
    //                         ui.selectable_value(&mut cfg.output_file_type, opt.to_string(), opt);
    //                     }
    //                 });
    //         });

    //         // Multiple files, set output directory
    //         ui.horizontal(|ui| {
    //             ui.label("Output Directory:");
            
    //             // 1. Grab the current directory and file‐name parts
    //             let (mut dir_string, file_name_os) = {
    //                 // Ensure we have a PathBuf for output
    //                 let out = &cfg.output_path;
    //                 let dir = out.parent().unwrap_or(Path::new("."));
    //                 let file_name = out.file_name().unwrap_or_default().to_owned();
    //                 (dir.to_string_lossy().to_string(), file_name)
    //             };
            
    //             // 2. Editable text box for the directory
    //             if ui
    //                 .add(TextEdit::singleline(&mut dir_string).hint_text("directory…"))
    //                 .changed()
    //             {
    //                 // When the user types a new dir, re‐join with the existing file name
    //                 cfg.output_path = Path::new(&dir_string).join(&file_name_os);
    //             }
            
    //             // 3. “Browse” button to pick a folder
    //             if ui.button("…").on_hover_text("Browse output directory").clicked() {
    //                 if let Some(dir) = FileDialog::new()
    //                     .set_title("Select Output Directory")
    //                     .set_directory(&dir_string) // start in the currently shown dir
    //                     .pick_folder()
    //                 {
    //                     // Update the path, preserving the file name
    //                     cfg.output_path = dir.join(&file_name_os);
    //                 }
    //             }
    //         });
        // }
    // }
}

/// Draw PyProphet settings and command toggles
pub fn draw_validation(ui: &mut Ui, config: &mut Input) {
    ui.heading("PyProphet Settings");

    // Ensure we have a PyProphetConfig
    let pyp_cfg = config
        .statistical_validation
        .get_or_insert_with(PyProphetConfig::default);

    // Auto-detect binary on first run
    if pyp_cfg.binary_path.as_os_str().is_empty() {
        if let Some(path) = find_executable("pyprophet", None) {
            pyp_cfg.binary_path = path;
        }
    }

    // Binary path + browse button
    ui.horizontal(|ui| {
        ui.label("PyProphet binary:");
        let mut path_str = pyp_cfg.binary_path.to_string_lossy().into_owned();
        if ui.add(TextEdit::singleline(&mut path_str)).changed() {
            pyp_cfg.binary_path = PathBuf::from(&path_str);
        }
        if ui.button("…").on_hover_text("Browse for pyprophet executable").clicked() {
            if let Some(file) = FileDialog::new()
                .set_title("Select pyprophet executable")
                .pick_file()
            {
                pyp_cfg.binary_path = file;
            }
        }
    });

    ui.separator();

    ui.label("Commands to run:");
    ui.horizontal(|ui| {
        ui.add(Checkbox::new(&mut pyp_cfg.run_merge, "Merge")).on_hover_text("Run `pyprophet merge`");
        ui.add(Checkbox::new(&mut pyp_cfg.run_score, "Score")).on_hover_text("Run `pyprophet score`");
        ui.add(Checkbox::new(&mut pyp_cfg.run_inference, "Inference")).on_hover_text("Run `pyprophet infer` subcommands");
        ui.add(Checkbox::new(&mut pyp_cfg.run_export, "Export")).on_hover_text("Run `pyprophet export` subcommands");
    });

    ui.separator();

    // Draw merge section if merge is enabled
    if pyp_cfg.run_merge {
        egui::CollapsingHeader::new("Merge Settings")
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Merged Output File:");
                    let mut merge_osw_output = pyp_cfg.merge_osw_output.to_string_lossy().into_owned();
                    if ui.add(TextEdit::singleline(&mut merge_osw_output)).changed() {
                        pyp_cfg.merge_osw_output = PathBuf::from(merge_osw_output);
                    }
                    if ui.button("…").on_hover_text("Browse to save output file").clicked() {
                        if let Some(file) = FileDialog::new()
                            .set_title("Select Output File")
                            .save_file()
                        {
                            pyp_cfg.merge_osw_output = file;
                        }
                    }
                });
            });
    }
    
    // Draw scoring section if scoring is enabled
    if pyp_cfg.run_score {
        draw_full_scoring_section(ui, pyp_cfg);
    }

    // Draw inference section if inference is enabled
    if pyp_cfg.run_inference {
        draw_full_inference_section(ui, pyp_cfg);
    }

    // Draw export section if export is enabled
    if pyp_cfg.run_export {
        draw_pyp_export_section(ui, pyp_cfg);
    }
 
}

fn draw_full_scoring_section(ui: &mut Ui, pyp_cfg: &mut PyProphetConfig) {
    if pyp_cfg.run_score {
        
        // Add collapsing header for scoring settings, default open
        egui::CollapsingHeader::new("Scoring Settings")
            .default_open(true)
            .show(ui, |ui| {
                // Checkboxes for level of scoring MS1, MS2, Transition and Alignment scoring
                ui.horizontal(|ui| {
                    ui.add(Checkbox::new(&mut pyp_cfg.score_ms1, "MS1")).on_hover_text("Score MS1 features");
                    ui.add(Checkbox::new(&mut pyp_cfg.score_ms2, "MS2")).on_hover_text("Score MS2 features");
                    ui.add(Checkbox::new(&mut pyp_cfg.score_transition, "Transition")).on_hover_text("Score transitions");
                    ui.add(Checkbox::new(&mut pyp_cfg.score_alignment, "Alignment")).on_hover_text("Score alignment features");
                });

                // Draw individual scoring sections for each level
                if pyp_cfg.score_ms1 {
                    draw_score_section(ui, ScoreSection {
                        label: "MS1",
                        id_salt: "ms1",
                        integrate_ms1: &mut false,
                        classifier: &mut pyp_cfg.classifier_ms1,
                        ss_num_iter: &mut pyp_cfg.ss_num_iter_ms1,
                        xeval_num_iter: &mut pyp_cfg.xeval_num_iter_ms1,
                        ss_initial_fdr: &mut pyp_cfg.ss_initial_fdr_ms1,
                        ss_iteration_fdr: &mut pyp_cfg.ss_iteration_fdr_ms1,
                        advanced_params: &mut pyp_cfg.advanced_params_ms1,
                    });
                }
                
                if pyp_cfg.score_ms2 {
                    draw_score_section(ui, ScoreSection {
                        label: "MS2",
                        id_salt: "ms2",
                        integrate_ms1: &mut pyp_cfg.integrate_ms1,
                        classifier: &mut pyp_cfg.classifier_ms2,
                        ss_num_iter: &mut pyp_cfg.ss_num_iter_ms2,
                        xeval_num_iter: &mut pyp_cfg.xeval_num_iter_ms2,
                        ss_initial_fdr: &mut pyp_cfg.ss_initial_fdr_ms2,
                        ss_iteration_fdr: &mut pyp_cfg.ss_iteration_fdr_ms2,
                        advanced_params: &mut pyp_cfg.advanced_params_ms2,
                    });
                }
                
                if pyp_cfg.score_transition {
                    draw_score_section(ui, ScoreSection {
                        label: "Transition",
                        id_salt: "transition",
                        integrate_ms1: &mut false,
                        classifier: &mut pyp_cfg.classifier_transition,
                        ss_num_iter: &mut pyp_cfg.ss_num_iter_transition,
                        xeval_num_iter: &mut pyp_cfg.xeval_num_iter_transition,
                        ss_initial_fdr: &mut pyp_cfg.ss_initial_fdr_transition,
                        ss_iteration_fdr: &mut pyp_cfg.ss_iteration_fdr_transition,
                        advanced_params: &mut pyp_cfg.advanced_params_transition,
                    });
                }
                
                if pyp_cfg.score_alignment {
                    draw_score_section(ui, ScoreSection {
                        label: "Alignment",
                        id_salt: "alignment",
                        integrate_ms1: &mut false,
                        classifier: &mut pyp_cfg.classifier_alignment,
                        ss_num_iter: &mut pyp_cfg.ss_num_iter_alignment,
                        xeval_num_iter: &mut pyp_cfg.xeval_num_iter_alignment,
                        ss_initial_fdr: &mut pyp_cfg.ss_initial_fdr_alignment,
                        ss_iteration_fdr: &mut pyp_cfg.ss_iteration_fdr_alignment,
                        advanced_params: &mut pyp_cfg.advanced_params_alignment,
                    });
                }
                
                
            });

    }
}

fn draw_score_section(ui: &mut Ui, section: ScoreSection<'_>) {

    ui.separator();
    ui.label(format!("{} Scoring Settings", section.label.to_uppercase()));

    // Checl if ScoreSection label is MS2, if it is, add integrate_ms1 checkbox
    if section.label == "MS2" {
        ui.horizontal(|ui| {
            ui.add(Checkbox::new(section.integrate_ms1, "Integrate MS1 features"))
                .on_hover_text("Integrate MS1 features with MS2 scoring, i.e. `--level=ms1ms2`");
        });
    }

    ComboBox::from_id_salt(section.id_salt)
        .selected_text(section.classifier.as_str())
        .show_ui(ui, |ui| {
            for candidate in [
                PyPScoreClassifier::LDA,
                PyPScoreClassifier::XGBoost,
                PyPScoreClassifier::SVM,
            ] {
                ui.selectable_value(
                    section.classifier,
                    candidate.as_str().to_string(),
                    candidate.as_str(),
                );
            }
        });

    ui.horizontal(|ui| {
        ui.horizontal(|ui| {
            ui.label("ss_num_iter:");
            ui.add(
                DragValue::new(section.ss_num_iter)
                    .speed(0.1)
                    .range(1..=100),
            ).on_hover_text("Number of iterations for semi-supervised learning step.");
        });
        ui.horizontal(|ui| {
            ui.label("xeval_num_iter:");
            ui.add(
                DragValue::new(section.xeval_num_iter)
                    .speed(0.1)
                    .range(1..=100),
            ).on_hover_text("Number of iterations for cross-validation of semi-supervised learning step..");
        });
    });

    ui.horizontal(|ui| {
        ui.horizontal(|ui| {
            ui.label("ss_initial_fdr:");
            ui.add(
                DragValue::new(section.ss_initial_fdr)
                    .speed(0.01)
                    .range(0.0..=1.0),
            ).on_hover_text("Initial FDR cutoff for best scoring targets.");
        });
        ui.horizontal(|ui| {
            ui.label("ss_iteration_fdr:");
            ui.add(
                DragValue::new(section.ss_iteration_fdr)
                    .speed(0.01)
                    .range(0.0..=1.0),
            ).on_hover_text("Iteration FDR cutoff for best scoring targets..");
        });
    });

    ui.label("Additional pyprophet flags:");
    ui.add(
        TextEdit::multiline(section.advanced_params)
            .desired_rows(3)
            .lock_focus(true)
            .hint_text("e.g. --parametric --threads 4"),
    );
}

fn draw_full_inference_section(ui: &mut Ui, pyp_cfg: &mut PyProphetConfig) {
    if pyp_cfg.run_inference {
        egui::CollapsingHeader::new("Inference Settings")
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.add(Checkbox::new(&mut pyp_cfg.run_infer_peptide, "Peptide")).on_hover_text("Run `pyprophet infer peptide`");
                    ui.add(Checkbox::new(&mut pyp_cfg.run_infer_peptidoform, "Peptidoform")).on_hover_text("Run `pyprophet infer peptidoform`");
                    ui.add(Checkbox::new(&mut pyp_cfg.run_infer_protein, "Protein")).on_hover_text("Run `pyprophet infer protein`");
                    ui.add(Checkbox::new(&mut pyp_cfg.run_infer_gene, "Gene")).on_hover_text("Run `pyprophet infer gene`"); 
                });

                if pyp_cfg.run_infer_peptide {
                    draw_pyprophet_infer_section(ui, PyPInferSection {
                        label: "Peptide",
                        id_salt: "peptide",
                        level: "peptide",
                        global: &mut pyp_cfg.infer_peptide_global,
                        experiment_wide: &mut pyp_cfg.infer_peptide_experiment_wide,
                        run_specific: &mut pyp_cfg.infer_peptide_run_specific,
                        ipf_ms1_scoring: &mut false,
                        ipf_ms2_scoring: &mut false,
                        ipf_max_precursor_pep: &mut 0.7,
                        ipf_max_peakgroup_pep: &mut 0.7,
                        ipf_max_precursor_peakgroup_pep: &mut 0.4,
                        ipf_max_transition_pep: &mut 0.6,
                        ipf_propgate_signal: &mut false,
                        ipf_max_alignment_pep: &mut 0.4,
                        infer_peptidoform_advanced_params: &mut String::new(),
                    });
                }
                if pyp_cfg.run_infer_peptidoform {
                    draw_pyprophet_infer_section(ui, PyPInferSection {
                        label: "Peptidoform",
                        id_salt: "peptidoform",
                        level: "peptidoform",
                        global: &mut false,
                        experiment_wide: &mut false,
                        run_specific: &mut false,
                        ipf_ms1_scoring: &mut pyp_cfg.infer_peptidoform_ms1_scoring,
                        ipf_ms2_scoring: &mut pyp_cfg.infer_peptidoform_ms2_scoring,
                        ipf_max_precursor_pep: &mut pyp_cfg.infer_peptidoform_max_precursor_pep,
                        ipf_max_peakgroup_pep: &mut pyp_cfg.infer_peptidoform_max_peakgroup_pep,
                        ipf_max_precursor_peakgroup_pep: &mut pyp_cfg.infer_peptidoform_max_precursor_peakgroup_pep,
                        ipf_max_transition_pep: &mut pyp_cfg.infer_peptidoform_max_transition_pep,
                        ipf_propgate_signal: &mut pyp_cfg.infer_peptidoform_propagate_signal,   
                        ipf_max_alignment_pep: &mut pyp_cfg.infer_peptidoform_max_alignment_pep,
                        infer_peptidoform_advanced_params: &mut pyp_cfg.infer_peptidoform_advanced_params,
                    });
                }
                if pyp_cfg.run_infer_protein {
                    draw_pyprophet_infer_section(ui, PyPInferSection {
                        label: "Protein",
                        id_salt: "protein",
                        level: "protein",
                        global: &mut pyp_cfg.infer_protein_global,
                        experiment_wide: &mut pyp_cfg.infer_protein_experiment_wide,
                        run_specific: &mut pyp_cfg.infer_protein_run_specific,
                        ipf_ms1_scoring: &mut false,
                        ipf_ms2_scoring: &mut false,
                        ipf_max_precursor_pep: &mut 0.7,
                        ipf_max_peakgroup_pep: &mut 0.7,
                        ipf_max_precursor_peakgroup_pep: &mut 0.4,
                        ipf_max_transition_pep: &mut 0.6,
                        ipf_propgate_signal: &mut false,
                        ipf_max_alignment_pep: &mut 0.4,
                        infer_peptidoform_advanced_params: &mut String::new(),
                    });
                }
                if pyp_cfg.run_infer_gene {
                    draw_pyprophet_infer_section(ui, PyPInferSection {
                        label: "Gene",
                        id_salt: "gene",
                        level: "gene",
                        global: &mut pyp_cfg.infer_gene_global,
                        experiment_wide: &mut pyp_cfg.infer_gene_experiment_wide,
                        run_specific: &mut pyp_cfg.infer_gene_run_specific,
                        ipf_ms1_scoring: &mut false,
                        ipf_ms2_scoring: &mut false,
                        ipf_max_precursor_pep: &mut 0.7,
                        ipf_max_peakgroup_pep: &mut 0.7,
                        ipf_max_precursor_peakgroup_pep: &mut 0.4,
                        ipf_max_transition_pep: &mut 0.6,
                        ipf_propgate_signal: &mut false,
                        ipf_max_alignment_pep: &mut 0.4,
                        infer_peptidoform_advanced_params: &mut String::new(),
                    });
                }
            });
    }
}

fn draw_pyprophet_infer_section(ui: &mut Ui, section: PyPInferSection<'_>)  {
    
    ui.separator();
    ui.label(format!("{} Inference Settings", section.label.to_uppercase()));

    if section.label != "Peptidoform" {
        ui.horizontal(|ui| {
            ui.add(Checkbox::new(section.global, "Global"))
                .on_hover_text("Run inference globally across all runs");
            ui.add(Checkbox::new(section.experiment_wide, "Experiment-wide"))
                .on_hover_text("Run inference across all runs in the experiment");
            ui.add(Checkbox::new(section.run_specific, "Run-specific"))
                .on_hover_text("Run inference for each run separately");
        });
    } else {
        ui.horizontal(|ui| {
            ui.add(Checkbox::new(section.ipf_ms1_scoring, "IPF MS1 Scoring"))
                .on_hover_text("Use IPF for MS1 scoring");
            ui.add(Checkbox::new(section.ipf_ms2_scoring, "IPF MS2 Scoring"))
                .on_hover_text("Use IPF for MS2 scoring");
        });
        
        ui.horizontal(|ui| {
            ui.horizontal(|ui| {
                ui.label("IPF Max Precursor PEP:");
                ui.add(
                    DragValue::new(section.ipf_max_precursor_pep)
                        .speed(0.01)
                        .range(0.0..=1.0),
                ).on_hover_text("Maximum precursor PEP for IPF");
            });
            ui.horizontal(|ui| {
                ui.label("IPF Max Peakgroup PEP:");
                ui.add(
                    DragValue::new(section.ipf_max_peakgroup_pep)
                        .speed(0.01)
                        .range(0.0..=1.0),
                ).on_hover_text("Maximum peakgroup PEP for IPF");
            });
        });

        ui.horizontal(|ui| {
            ui.horizontal(|ui| {
                ui.label("IPF Max Precursor Peakgroup PEP:");
                ui.add(
                    DragValue::new(section.ipf_max_precursor_peakgroup_pep)
                        .speed(0.01)
                        .range(0.0..=1.0),
                ).on_hover_text("Maximum precursor peakgroup PEP for IPF");
            });
            ui.horizontal(|ui| {
                ui.label("IPF Max Transition PEP:");
                ui.add(
                    DragValue::new(section.ipf_max_transition_pep)
                        .speed(0.01)
                        .range(0.0..=1.0),
                ).on_hover_text("Maximum transition PEP for IPF");
            });
        });

        ui.horizontal(|ui| {
            ui.horizontal(|ui| {
                ui.add(Checkbox::new(section.ipf_propgate_signal, "IPF Propagate Signal"))
                    .on_hover_text("Propagate signal across runs using IPF");
            });
            if *section.ipf_propgate_signal {
                ui.horizontal(|ui| {
                    ui.label("IPF Max Alignment PEP:");
                    ui.add(
                        DragValue::new(section.ipf_max_alignment_pep)
                            .speed(0.01)
                            .range(0.0..=1.0),
                    ).on_hover_text("Maximum alignment PEP for IPF");
                });
            }
        });

        ui.label("Additional pyprophet flags:");
        ui.add(
            TextEdit::multiline(section.infer_peptidoform_advanced_params)
                .desired_rows(3)
                .lock_focus(true)
                .hint_text("e.g. --ipf_grouped_fdr"),
        );
    }
}

fn draw_pyp_export_section(ui: &mut Ui, pyp_cfg: &mut PyProphetConfig) {
    if pyp_cfg.run_export {
        egui::CollapsingHeader::new("Export Settings")
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.add(Checkbox::new(&mut pyp_cfg.export_tsv, "TSV")).on_hover_text("Run `pyprophet export tsv`");
                    ui.add(Checkbox::new(&mut pyp_cfg.export_precursor_matrix, "Precursor Matrix")).on_hover_text("Run `pyprophet export matrix --level peptide`");
                    ui.add(Checkbox::new(&mut pyp_cfg.export_peptide_matrix, "Peptide Matrix")).on_hover_text("Run `pyprophet export matrix --level peptide`");
                    ui.add(Checkbox::new(&mut pyp_cfg.export_protein_matrix, "Protein Matrix")).on_hover_text("Run `pyprophet export matrix --level protein`");
                    ui.add(Checkbox::new(&mut pyp_cfg.export_parquet, "Parquet")).on_hover_text("Run `pyprophet export parquet`");
                }); 

                ui.horizontal(|ui| {
                    ui.label("Output Directory:");
                    let mut export_output_path = pyp_cfg.export_output_path.to_string_lossy().into_owned();
                    if ui.add(TextEdit::singleline(&mut export_output_path)).changed() {
                        pyp_cfg.export_output_path = PathBuf::from(export_output_path);
                    }
                    if ui.button("…").on_hover_text("Browse to save output file").clicked() {
                        if let Some(file) = FileDialog::new()
                            .set_title("Select Output File")
                            .pick_folder()
                        {
                            pyp_cfg.export_output_path = file;
                        }
                    }
                });
                

                // Add horizontal numeric sliders for ipf_max_peptidoform_pep, max_rs_peakgroup_qvalue
                ui.horizontal(|ui| {
                    ui.horizontal(|ui| {
                        ui.label("Max Peptidoform PEP:");
                        ui.add(
                            DragValue::new(&mut pyp_cfg.ipf_max_peptidoform_pep)
                                .speed(0.01)
                                .range(0.0..=1.0),
                        ).on_hover_text("IPF: Filter results to maximum run-specific peptidoform-level PEP");
                    });
                    ui.horizontal(|ui| {
                        ui.label("Max Peakgroup QValue:");
                        ui.add(
                            DragValue::new(&mut pyp_cfg.max_rs_peakgroup_qvalue)
                                .speed(0.01)
                                .range(0.0..=1.0),
                        ).on_hover_text("Filter results to maximum run-specific peak group-level q-value");
                    });
                });

                ui.horizontal(|ui| {
                    ui.horizontal(|ui| {
                        ui.label("Max Global Peptide QValue:");
                        ui.add(
                            DragValue::new(&mut pyp_cfg.max_global_peptide_qvalue)
                                .speed(0.01)
                                .range(0.0..=1.0),
                        ).on_hover_text("Maximum peptide QValue");
                    });
                    ui.horizontal(|ui| {
                        ui.label("Max Global Protein QValue:");
                        ui.add(
                            DragValue::new(&mut pyp_cfg.max_global_protein_qvalue)
                                .speed(0.01)
                                .range(0.0..=1.0),
                        ).on_hover_text("Maximum protein QValue");
                    });
                });

                // Check if any of the matrix export options are enabled
                if pyp_cfg.export_precursor_matrix || pyp_cfg.export_peptide_matrix || pyp_cfg.export_protein_matrix {
                    ui.add_space(1.0);
                    ui.separator();
                    ui.label("Matrix export settings:");
                    ui.horizontal(|ui| {
                        ui.label("Top N Feature Summarization:");
                        ui.add(
                            DragValue::new(&mut pyp_cfg.top_n)
                                .speed(1.0)
                                .range(1..=100),
                        ).on_hover_text("Number of top intense features to use for summarization in matrix exports");
                        ui.add(Checkbox::new(&mut pyp_cfg.consistent_top_n, "Consistent Feature Summarization:"))
                            .on_hover_text("Whether to use same top features across all runs");
                    });
                    ui.horizontal(|ui| {
                        ui.label("Matrix Normalization Method:");
                        ComboBox::from_id_salt("feature_summarization_normalization_method")
                            .selected_text(&pyp_cfg.normalization_method)
                            .show_ui(ui, |ui| {
                                for method in ["none", "median", "medianmedian", "quantile"] {
                                    ui.selectable_value(&mut pyp_cfg.normalization_method, method.to_string(), method);
                                }
                            });
                    });
                }


                // Check for export_parquet
                if pyp_cfg.export_parquet {
                    ui.add_space(1.0);
                    ui.separator();
                    ui.label("Parquet Export Settings:");

                    // Get extension of first file in pyp_cfg.file_paths
                    let first_file_ext = pyp_cfg.file_paths.first()
                        .and_then(|p| p.extension())
                        .and_then(|s| s.to_str())
                        .unwrap_or("osw");

                    if first_file_ext.to_lowercase() == "sqmass" {
                       // textbox and browse button for pyp_cfg.pqpfile input file
                        ui.horizontal(|ui| {
                            ui.add(egui::Label::new("PQP file (sqMass conversion):")).on_hover_text("PyProphet PQP file. Only required when converting sqMass to parquet.");
                            let mut pqp_input_file = pyp_cfg.pqpfile.to_string_lossy().into_owned();
                            if ui.add(TextEdit::singleline(&mut pqp_input_file)).changed() {
                                pyp_cfg.pqpfile = PathBuf::from(pqp_input_file);
                            }
                            if ui.button("…").on_hover_text("Browse for input file").clicked() {
                                if let Some(file) = FileDialog::new()
                                    .set_title("Select Parquet Input File")
                                    .pick_file()
                                {
                                    pyp_cfg.pqpfile = file;
                                }
                            }
                        }); 
                    }

                    // textbox and browse button for output export_parquet_output_path file
                    ui.horizontal(|ui| {
                        ui.add(egui::Label::new("Parquet Output Path:")).on_hover_text("Output path for parquet file.");
                        let mut parquet_output_path = pyp_cfg.export_parquet_output_path.to_string_lossy().into_owned();
                        if ui.add(TextEdit::singleline(&mut parquet_output_path)).changed() {
                            pyp_cfg.export_parquet_output_path = PathBuf::from(parquet_output_path);
                        }
                        if ui.button("…").on_hover_text("Browse to save output file").clicked() {
                            if let Some(file) = FileDialog::new()
                                .set_title("Select Parquet Output File")
                                .save_file()
                            {
                                pyp_cfg.export_parquet_output_path = file;
                            }
                        }
                    });

                    // checkbox for split_transition_data and split_runs
                    ui.horizontal(|ui| {
                        ui.add(Checkbox::new(&mut pyp_cfg.split_transition_data, "Split Transition Data"))
                            .on_hover_text("Split transition data into separate files");
                        ui.add(Checkbox::new(&mut pyp_cfg.split_runs, "Split Runs"))
                            .on_hover_text("Split runs into separate files");
                    });

                }

            });
    }
}