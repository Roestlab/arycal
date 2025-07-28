use std::process::Child;
use std::thread::JoinHandle;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
use arycal_cloudpath::util::{extract_basename, extract_directory};
use arycal_common::config::PyProphetConfig;
use chrono::Local;
use eframe::egui::Ui;
use arycal_cli::input::Input;
use egui::{Button, Color32, ScrollArea};
use rfd::FileDialog;
use std::{fs::File, path::PathBuf, sync::{
    atomic::{AtomicBool, Ordering}, mpsc::{channel, Receiver, Sender}, Arc, Mutex
}, time::Duration};

use super::open_swath_tab::{spawn_batch_runner, spawn_process_runner, Progress};

// ------------------ Arg list builders ------------------

/// Build the `pyprophet merge` invocation(s).
/// 
/// Returns a Vec of argument‐lists; in practice this will be either:
/// - empty (if no inputs), or
/// - one element: `["merge", "--out", <output>, <in1>, <in2>, …]`
pub fn build_merge_args(cfg: &PyProphetConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();
    // nothing to merge if no files
    if cfg.file_paths.is_empty() {
        return all_args;
    }

    // pyprophet merge --out <output_path> <file1> <file2> …
    let mut args = Vec::new();
    args.push("--no-log-colorize".to_string());
    args.push("merge".to_string());
    args.push("osw".to_string());
    args.push("--out".to_string());
    args.push(cfg.merge_osw_output.to_string_lossy().into_owned());
    // set --template to first input file, if any
    if let Some(first_input) = cfg.file_paths.first() {
        args.push("--template".to_string());
        args.push(first_input.to_string_lossy().into_owned());
    }
    for input in &cfg.file_paths {
        args.push(input.to_string_lossy().into_owned());
    }

    all_args.push(args);
    all_args
}

/// Build the per‐file arg‐lists for one `pyprophet score --level <level>` invocation.
/// Each enabled level becomes one Stage with N runs (one per input file).
pub fn build_score_level_args(
    cfg: &PyProphetConfig,
    level: &str,
    classifier: &str,
    ss_num_iter: usize,
    xeval_num_iter: usize,
    ss_initial_fdr: f32,
    ss_iteration_fdr: f32,
    advanced: &str,
) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();

    let file_paths = if cfg.run_merge {
        // if merging, use the merged output path
        vec![cfg.merge_osw_output.clone()]
    } else {
        // otherwise, use the original file paths
        cfg.file_paths.clone()
    };

    for input in file_paths {
        let mut args = Vec::new();
        // subcommand and level
        args.push("--no-log-colorize".to_string());
        args.push("score".to_string());
        args.push("--level".to_string());
        args.push(level.to_string());

        // classifier
        args.push("--classifier".to_string());
        args.push(classifier.to_string());

        // semi‐supervised iters
        args.push("--ss_num_iter".to_string());
        args.push(ss_num_iter.to_string());

        // cross‐validation iters
        args.push("--xeval_num_iter".to_string());
        args.push(xeval_num_iter.to_string());

        // FDR cutoffs
        args.push("--ss_initial_fdr".to_string());
        args.push(ss_initial_fdr.to_string());
        args.push("--ss_iteration_fdr".to_string());
        args.push(ss_iteration_fdr.to_string());

        // any user‐supplied advanced flags
        if !advanced.trim().is_empty() {
            args.extend(advanced.split_whitespace().map(|s| s.to_string()));
        }

        // finally: the input file
        args.push("--in".to_string());
        args.push(input.to_string_lossy().into_owned());

        all_args.push(args);
    }

    all_args
}

/// Push up to three *separate* stages for peptide inference, one per context flag.
fn push_infer_peptide_stages(pipeline: &mut Vec<Stage>, cfg: &PyProphetConfig) {
    let bin = cfg.binary_path.clone();
    let inputs = if cfg.run_merge {
        // if merging, use the merged output path
        vec![cfg.merge_osw_output.clone()]
    } else {
        // otherwise, use the original file paths
        cfg.file_paths.clone()
    };
    for input in inputs {
        if cfg.infer_peptide_global {
            pipeline.push(Stage {
                tool: ValidationTool::InferPeptideGlobal,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "--no-log-colorize".into(),
                    "infer".into(), "peptide".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "global".into(),
                ]],
                display_name: "pyprophet infer peptide (global)",
            });
        }
        if cfg.infer_peptide_experiment_wide {
            pipeline.push(Stage {
                tool: ValidationTool::InferPeptideExperimentWide,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "--no-log-colorize".into(),
                    "infer".into(), "peptide".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "experiment-wide".into(),
                ]],
                display_name: "pyprophet infer peptide (experiment-wide)",
            });
        }
        if cfg.infer_peptide_run_specific {
            pipeline.push(Stage {
                tool: ValidationTool::InferPeptideRunSpecific,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "--no-log-colorize".into(),
                    "infer".into(), "peptide".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "run-specific".into(),
                ]],
                display_name: "pyprophet infer peptide (run-specific)",
            });
        }
    }
}

fn push_infer_protein_stages(pipeline: &mut Vec<Stage>, cfg: &PyProphetConfig) {
    let bin = cfg.binary_path.clone();
    let inputs = if cfg.run_merge {
        // if merging, use the merged output path
        vec![cfg.merge_osw_output.clone()]
    } else {
        // otherwise, use the original file paths
        cfg.file_paths.clone()
    };
    for input in inputs {
        if cfg.infer_protein_global {
            pipeline.push(Stage {
                tool: ValidationTool::InferProteinGlobal,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "infer".into(), "protein".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "global".into(),
                ]],
                display_name: "pyprophet infer protein (global)",
            });
        }
        if cfg.infer_protein_experiment_wide {
            pipeline.push(Stage {
                tool: ValidationTool::InferProteinExperimentWide,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "infer".into(), "protein".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "experiment-wide".into(),
                ]],
                display_name: "pyprophet infer protein (experiment-wide)",
            });
        }
        if cfg.infer_protein_run_specific {
            pipeline.push(Stage {
                tool: ValidationTool::InferProteinRunSpecific,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "infer".into(), "protein".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "run-specific".into(),
                ]],
                display_name: "pyprophet infer protein (run-specific)",
            });
        }
    }
}

fn push_infer_gene_stages(pipeline: &mut Vec<Stage>, cfg: &PyProphetConfig) {
    let bin = cfg.binary_path.clone();
    let inputs = if cfg.run_merge {
        // if merging, use the merged output path
        vec![cfg.merge_osw_output.clone()]
    } else {
        // otherwise, use the original file paths
        cfg.file_paths.clone()
    };
    for input in inputs {
        if cfg.infer_gene_global {
            pipeline.push(Stage {
                tool: ValidationTool::InferGeneGlobal,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "infer".into(), "gene".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "global".into(),
                ]],
                display_name: "pyprophet infer gene (global)",
            });
        }
        if cfg.infer_gene_experiment_wide {
            pipeline.push(Stage {
                tool: ValidationTool::InferGeneExperimentWide,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "infer".into(), "gene".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "experiment-wide".into(),
                ]],
                display_name: "pyprophet infer gene (experiment-wide)",
            });
        }
        if cfg.infer_gene_run_specific {
            pipeline.push(Stage {
                tool: ValidationTool::InferGeneRunSpecific,
                bin: bin.clone(),
                arg_lists: vec![vec![
                    "infer".into(), "gene".into(),
                    "--in".into(), input.clone().to_string_lossy().into_owned(),
                    "--context".into(), "run-specific".into(),
                ]],
                display_name: "pyprophet infer gene (run-specific)",
            });
        }
    }
}

/// Build the `pyprophet infer peptidoform` invocations, one per input file.
/// Each run will include only the IPF flags the user has enabled/configured.
///
/// # CLI shape
/// ```text
/// pyprophet infer peptidoform --in <input>
///     [--ipf_ms1_scoring]
///     [--ipf_ms2_scoring]
///     --ipf_max_precursor_pep <f32>
///     --ipf_max_peakgroup_pep <f32>
///     --ipf_max_precursor_peakgroup_pep <f32>
///     --ipf_max_transition_pep <f32>
///     [--ipf_propagate_signal --ipf_max_alignment_pep <f32>]
/// ```
///
/// Returns a Vec of argument‐lists; each inner Vec<String> is one process invocation.
pub fn build_infer_peptidoform_args(cfg: &PyProphetConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();

    let file_paths = if cfg.run_merge {
        // if merging, use the merged output path
        vec![cfg.merge_osw_output.clone()]
    } else {
        // otherwise, use the original file paths
        cfg.file_paths.clone()
    };

    for input in file_paths {
        let mut args = Vec::new();
        // subcommand + level
        args.push("--no-log-colorize".to_string());
        args.push("infer".into());
        args.push("peptidoform".into());

        // input
        args.push("--in".into());
        args.push(input.to_string_lossy().into_owned());

        // IPF MS1 scoring?
        if cfg.infer_peptidoform_ms1_scoring {
            args.push("--ipf_ms1_scoring".into());
        }
        else {
            args.push("--no-ipf_ms1_scoring".into());
        }
        // IPF MS2 scoring?
        if cfg.infer_peptidoform_ms2_scoring {
            args.push("--ipf_ms2_scoring".into());
        } else {
            args.push("--no-ipf_ms2_scoring".into());
        }

        // Required PEP thresholds
        args.push("--ipf_max_precursor_pep".into());
        args.push(cfg.infer_peptidoform_max_precursor_pep.to_string());
        args.push("--ipf_max_peakgroup_pep".into());
        args.push(cfg.infer_peptidoform_max_peakgroup_pep.to_string());
        args.push("--ipf_max_precursor_peakgroup_pep".into());
        args.push(cfg.infer_peptidoform_max_precursor_peakgroup_pep.to_string());
        args.push("--ipf_max_transition_pep".into());
        args.push(cfg.infer_peptidoform_max_transition_pep.to_string());

        // Optional propagate‐signal & alignment‐PEP
        if cfg.infer_peptidoform_propagate_signal {
            args.push("--propagate_signal_across_runs".into());
            args.push("--ipf_max_alignment_pep".into());
            args.push(cfg.infer_peptidoform_max_alignment_pep.to_string());
        }

        // any user‐supplied advanced flags
        if !cfg.infer_peptidoform_advanced_params.trim().is_empty() {
            args.extend(cfg.infer_peptidoform_advanced_params.split_whitespace().map(|s| s.to_string()));
        }

        all_args.push(args);
    }

    all_args
}


/// Build the `pyprophet export tsv` invocations, one per input file.
/// If there is exactly one input, writes to `cfg.output_path`.
/// If multiple inputs, treats `cfg.output_path` as a directory and emits
/// `<stem>_export.tsv` into it.
/// Always appends these filters:
///   --ipf_max_peptidoform_pep
///   --max_rs_peakgroup_qvalue
///   --max_global_peptide_qvalue
///   --max_global_protein_qvalue
pub fn build_export_tsv_args(cfg: &PyProphetConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();

    let file_paths = if cfg.run_merge {
        // if merging, use the merged output path
        vec![cfg.merge_osw_output.clone()]
    } else {
        // otherwise, use the original file paths
        cfg.file_paths.clone()
    };

    for input in file_paths {
        let mut args = Vec::new();
        args.push("--no-log-colorize".to_string());
        args.push("export".into());
        args.push("tsv".into());

        // input file
        args.push("--in".into());
        args.push(input.to_string_lossy().into_owned());

        // Build output path:
        // Get base name of input file, and if cfg.output_path is a directory,
        let basename = extract_basename(&input);
        let out_path = if cfg.export_output_path.is_dir() {
            // if output path is a directory, use it
            cfg.export_output_path.join(format!("{basename}_export.tsv"))
        } else {
            // otherwise, replace the extension of the input file with `_export.tsv`
            let dir = extract_directory(&cfg.export_output_path).unwrap_or_else(|| PathBuf::from("."));
            dir.join(format!("{basename}_export.tsv"))
        };

        args.push("--out".into());
        args.push(out_path.to_string_lossy().into_owned());

        // now your additional PEP/Q‐value filters:
        args.push("--ipf_max_peptidoform_pep".into());
        args.push(cfg.ipf_max_peptidoform_pep.to_string());

        args.push("--max_rs_peakgroup_qvalue".into());
        args.push(cfg.max_rs_peakgroup_qvalue.to_string());

        args.push("--max_global_peptide_qvalue".into());
        args.push(cfg.max_global_peptide_qvalue.to_string());

        args.push("--max_global_protein_qvalue".into());
        args.push(cfg.max_global_protein_qvalue.to_string());

        all_args.push(args);
    }

    all_args
}

/// Build one or more `pyprophet export matrix --level <level>` invocations,
/// using either the merged OSW (if `run_merge`) or each original file.
/// Writes into `cfg.export_output_path` (directory or file-root).
pub fn build_export_matrix_level_args(
    cfg: &PyProphetConfig,
    level: &str,
    enabled: bool,
) -> Vec<Vec<String>> {
    if !enabled {
        return Vec::new();
    }

    let mut all_args = Vec::new();

    // decide which inputs to iterate
    let inputs: Vec<PathBuf> = if cfg.run_merge {
        vec![cfg.merge_osw_output.clone()]
    } else {
        cfg.file_paths.clone()
    };

    for input in inputs {
        // turn PathBuf into a string for the --in flag
        let input_str = input.to_string_lossy().into_owned();

        // build output path: if export_output_path is a dir, use it directly;
        // otherwise treat it as a file-root and replace with level matrix name
        let out_path = if cfg.export_output_path.is_dir() {
            let stem = extract_basename(&input);
            cfg.export_output_path.join(format!("{stem}_{level}_matrix.tsv"))
        } else {
            let stem = extract_basename(&input);
            let dir = extract_directory(&cfg.export_output_path)
                .unwrap_or_else(|| PathBuf::from("."));
            dir.join(format!("{stem}_{level}_matrix.tsv"))
        };
        let out_str = out_path.to_string_lossy().into_owned();

        // assemble the args
        let mut args = Vec::new();
        args.push("--no-log-colorize".to_string());
        args.push("export".into());
        args.push("matrix".into());
        args.push("--in".into());
        args.push(input_str);
        args.push("--out".into());
        args.push(out_str);
        args.push("--level".into());
        args.push(level.into());

        // IPF filters
        args.push("--ipf_max_peptidoform_pep".into());
        args.push(cfg.ipf_max_peptidoform_pep.to_string());
        args.push("--max_rs_peakgroup_qvalue".into());
        args.push(cfg.max_rs_peakgroup_qvalue.to_string());
        args.push("--max_global_peptide_qvalue".into());
        args.push(cfg.max_global_peptide_qvalue.to_string());
        args.push("--max_global_protein_qvalue".into());
        args.push(cfg.max_global_protein_qvalue.to_string());

        // summarization flags
        args.push("--top_n".into());
        args.push(cfg.top_n.to_string());
        if cfg.consistent_top_n {
            args.push("--consistent_top".into());
        } else {
            args.push("--no-consistent_top".into());
        }

        // normalization
        args.push("--normalization".into());
        args.push(cfg.normalization_method.clone());

        all_args.push(args);
    }

    all_args
}

/// One-stage precursor-level matrix export
pub fn build_export_matrix_precursor_args(cfg: &PyProphetConfig) -> Vec<Vec<String>> {
    build_export_matrix_level_args(cfg, "precursor", cfg.export_precursor_matrix)
}

/// One-stage peptide-level matrix export
pub fn build_export_matrix_peptide_args(cfg: &PyProphetConfig) -> Vec<Vec<String>> {
    build_export_matrix_level_args(cfg, "peptide", cfg.export_peptide_matrix)
}

/// One-stage protein-level matrix export
pub fn build_export_matrix_protein_args(cfg: &PyProphetConfig) -> Vec<Vec<String>> {
    build_export_matrix_level_args(cfg, "protein", cfg.export_protein_matrix)
}

pub fn build_export_parquet_args(cfg: &PyProphetConfig) -> Vec<Vec<String>> {
    let mut all_args = Vec::new();

    let file_paths = if cfg.run_merge {
        // if merging, use the merged output path
        vec![cfg.merge_osw_output.clone()]
    } else {
        // otherwise, use the original file paths
        cfg.file_paths.clone()
    };

    for input in file_paths {
        let mut args = Vec::new();
        args.push("--no-log-colorize".to_string());
        args.push("export".into());
        args.push("parquet".into());

        // input file
        args.push("--in".into());
        args.push(input.to_string_lossy().into_owned());

        // Build output path:
        args.push("--out".into());
        args.push(cfg.export_parquet_output_path.to_string_lossy().into_owned());


        // if pqpfile is not None or not pqpfile.is_empty():
        if !cfg.pqpfile.to_string_lossy().is_empty() {
            args.push("--pqpfile".into());
            args.push(cfg.pqpfile.to_string_lossy().into_owned());
        }
        
        if cfg.split_transition_data {
            args.push("--split_transition_data".into());
        } else {
            args.push("--no-split_transition_data".into());
        }
        
        if cfg.split_runs {
            args.push("--split_runs".into());
        } else {
            args.push("--no-split_runs".into());
        }

        all_args.push(args);
    }

    all_args
}

// ------------------ Pipeline definitions ------------------


#[derive(Clone, Debug)]
pub enum ValidationTool {
    Merge,
    ScoreMs1,
    ScoreMs2,
    ScoreTransition,
    ScoreAlignment,
    InferPeptideGlobal,
    InferPeptideExperimentWide,
    InferPeptideRunSpecific,
    InferProteinGlobal,
    InferProteinExperimentWide,
    InferProteinRunSpecific,
    InferGeneGlobal,
    InferGeneExperimentWide,
    InferGeneRunSpecific,

    InferPeptidoform,
    ExportTsv,
    ExportMatrixPrecursor,
    ExportMatrixPeptide,
    ExportMatrixProtein,
}

impl ValidationTool {
    pub fn as_str(&self) -> &'static str {
        match self {
            ValidationTool::Merge           => "merge",
            ValidationTool::ScoreMs1      => "score --level=ms1",
            ValidationTool::ScoreMs2      => "score --level=ms2",
            ValidationTool::ScoreTransition => "score --level=transition",
            ValidationTool::ScoreAlignment => "score --level=alignment",
            ValidationTool::InferPeptideGlobal => "infer peptide --context=global",
            ValidationTool::InferPeptideExperimentWide => "infer peptide --context=experiment-wide",
            ValidationTool::InferPeptideRunSpecific => "infer peptide --context=run-specific",
            ValidationTool::InferProteinGlobal => "infer protein --context=global",
            ValidationTool::InferProteinExperimentWide => "infer protein --context=experiment-wide",
            ValidationTool::InferProteinRunSpecific => "infer protein --context=run-specific",
            ValidationTool::InferGeneGlobal => "infer gene --context=global",
            ValidationTool::InferGeneExperimentWide => "infer gene --context=experiment-wide",
            ValidationTool::InferGeneRunSpecific => "infer gene --context=run-specific",
            ValidationTool::InferPeptidoform => "infer peptidoform",
            ValidationTool::ExportTsv       => "export tsv",
            ValidationTool::ExportMatrixPrecursor => "export matrix precursor",
            ValidationTool::ExportMatrixPeptide => "export matrix peptide",
            ValidationTool::ExportMatrixProtein => "export matrix protein",
        }
    }
}

#[derive(Clone)]
pub struct Stage {
    pub tool: ValidationTool,
    pub bin: PathBuf,            // pyp_cfg.binary_path.clone()
    pub arg_lists: Vec<Vec<String>>,
    pub display_name: &'static str, // e.g. "pyprophet merge"
}

pub struct StageRunner {
    pub handle: JoinHandle<()>,
    pub progress_rx: Receiver<Progress>,
    pub stage: Stage,
}

#[derive(PartialEq)]
enum ValidationTabs {
    Console,
    Results,
}

pub struct ValidationState {
    pub active_tab: ValidationTabs,
    pub output_lines: Arc<Mutex<Vec<String>>>,
    pub cancel_flag: Arc<AtomicBool>,
    pub n_concurrent_processes: usize, // Number of concurrent processes to run
    pub current_children: Arc<Mutex<Vec<Arc<Mutex<Child>>>>>,
    pub current_stage_runner: Option<StageRunner>,

    pub pipeline: Vec<Stage>,
    pub current_stage_index: usize,

    // per‐stage metrics
    pub total_number_of_runs: usize,
    pub processing_run_n: usize,
    pub stage_successful_runs: usize,
    pub stage_failed_runs: usize,
    pub processing_time: Vec<Duration>,
    pub avg_processing_time: Duration,
    pub initial_start_time: Instant,
    pub total_processing_time: Duration,
    pub processed_fraction: f32,

    // global progress
    pub successful_runs: usize,
    pub failed_runs: usize,
    pub overall_total_runs: usize,
}

impl ValidationState {
    pub fn default() -> Self { Self {
        active_tab: ValidationTabs::Console,
        output_lines: Arc::new(Mutex::new(Vec::new())),
        cancel_flag: Arc::new(AtomicBool::new(false)),
        n_concurrent_processes: 1, // Default to 1 concurrent process
        current_children: Arc::new(Mutex::new(Vec::new())),
        current_stage_runner: None,
        pipeline: Vec::new(),
        current_stage_index: 0,
        total_number_of_runs: 0,
        processing_run_n: 0,
        stage_successful_runs: 0,
        stage_failed_runs: 0,
        processing_time: Vec::new(),
        avg_processing_time: Duration::ZERO,
        initial_start_time: Instant::now(),
        total_processing_time: Duration::ZERO,
        processed_fraction: 0.0,
        successful_runs: 0,
        failed_runs: 0,
        overall_total_runs: 0,
    } }

    pub fn ui(&mut self, ui: &mut Ui, config: &Input) {
        // update n_concurrent_processes from config
        self.n_concurrent_processes = config.n_concurrent_processes;
        if let Some(pyp_cfg) = config.statistical_validation.as_ref() {
            if pyp_cfg.binary_path.to_string_lossy().is_empty() || !pyp_cfg.binary_path.exists() {
                ui.colored_label(
                    Color32::RED,
                    "PyProphet missing: install with `pip install pyprophet` and set the binary path in the config settings."
                );
            } else {
                ui.horizontal(|ui| {
                    if ui
                        .selectable_label(self.active_tab == ValidationTabs::Console, "Console Log")
                        .clicked()
                    {
                        self.active_tab = ValidationTabs::Console;
                    }
                    // if ui
                    //     .selectable_label(self.active_tab == ValidationTabs::Results, "Results")
                    //     .clicked()
                    // {
                    //     self.active_tab = ValidationTabs::Results;
                    // }
                });
                ui.separator();

                self.console_ui(ui, config);
        
            //     match self.active_tab {
            //         ValidationTabs::Console => self.console_ui(ui, config),
            //         ValidationTabs::Results => self.results_ui(ui, config),
            //     }
            }
        }
    }

    pub fn is_running(&self) -> bool {
        self.current_stage_runner.is_some()
    }

    /// Called when ▶︎ Run Workflow is clicked
    fn start_pipeline(&mut self, cfg: &PyProphetConfig) {
        // 1) Clear any previous cancel signal so new runs actually start
        self.cancel_flag.store(false, Ordering::Relaxed);

        // 2) Reset *all* global & per-stage counters
        self.pipeline.clear();
        self.current_stage_index = 0;
        self.successful_runs = 0;
        self.failed_runs = 0;
        self.stage_successful_runs = 0;
        self.stage_failed_runs = 0;
        self.processing_time.clear();
        self.avg_processing_time = Duration::ZERO;
        self.total_processing_time = Duration::ZERO;
        self.processed_fraction = 0.0;

        if cfg.run_merge {
            self.pipeline.push(Stage {
                tool: ValidationTool::Merge,
                bin: cfg.binary_path.clone(),
                arg_lists: build_merge_args(cfg),
                display_name: "pyprophet merge",
            });
        }
        if cfg.run_score {
            // MS1
            if cfg.score_ms1 {
                let args = build_score_level_args(
                    cfg,
                    "ms1",
                    &cfg.classifier_ms1,
                    cfg.ss_num_iter_ms1,
                    cfg.xeval_num_iter_ms1,
                    cfg.ss_initial_fdr_ms1,
                    cfg.ss_iteration_fdr_ms1,
                    &cfg.advanced_params_ms1,
                );
                self.pipeline.push(Stage {
                    tool: ValidationTool::ScoreMs1,
                    bin: cfg.binary_path.clone(),
                    arg_lists: args,
                    display_name: "pyprophet score ms1",
                });
            }

            // MS2 (or ms1ms2)
            if cfg.score_ms2 {
                let (level, disp) = if cfg.integrate_ms1 {
                    ("ms1ms2", "pyprophet score ms1ms2")
                } else {
                    ("ms2", "pyprophet score ms2")
                };
                let args = build_score_level_args(
                    cfg,
                    level,
                    &cfg.classifier_ms2,
                    cfg.ss_num_iter_ms2,
                    cfg.xeval_num_iter_ms2,
                    cfg.ss_initial_fdr_ms2,
                    cfg.ss_iteration_fdr_ms2,
                    &cfg.advanced_params_ms2,
                );
                self.pipeline.push(Stage {
                    tool: ValidationTool::ScoreMs2,
                    bin: cfg.binary_path.clone(),
                    arg_lists: args,
                    display_name: disp,
                });
            }

            // transition
            if cfg.score_transition {
                let args = build_score_level_args(
                    cfg,
                    "transition",
                    &cfg.classifier_transition,
                    cfg.ss_num_iter_transition,
                    cfg.xeval_num_iter_transition,
                    cfg.ss_initial_fdr_transition,
                    cfg.ss_iteration_fdr_transition,
                    &cfg.advanced_params_transition,
                );
                self.pipeline.push(Stage {
                    tool: ValidationTool::ScoreTransition,
                    bin: cfg.binary_path.clone(),
                    arg_lists: args,
                    display_name: "pyprophet score transition",
                });
            }

            // alignment
            if cfg.score_alignment {
                let args = build_score_level_args(
                    cfg,
                    "alignment",
                    &cfg.classifier_alignment,
                    cfg.ss_num_iter_alignment,
                    cfg.xeval_num_iter_alignment,
                    cfg.ss_initial_fdr_alignment,
                    cfg.ss_iteration_fdr_alignment,
                    &cfg.advanced_params_alignment,
                );
                self.pipeline.push(Stage {
                    tool: ValidationTool::ScoreAlignment,
                    bin: cfg.binary_path.clone(),
                    arg_lists: args,
                    display_name: "pyprophet score alignment",
                });
            }
        }
        
        // Inference stages
        if cfg.run_inference {
            if cfg.run_infer_peptide {
                push_infer_peptide_stages(&mut self.pipeline, cfg);
            }
            if cfg.run_infer_protein {
                push_infer_protein_stages(&mut self.pipeline, cfg);
            }
            if cfg.run_infer_gene {
                push_infer_gene_stages(&mut self.pipeline, cfg);
            }
            if cfg.run_infer_peptidoform {
                let pf_args = build_infer_peptidoform_args(cfg);
                if !pf_args.is_empty() {
                    self.pipeline.push(Stage {
                        tool: ValidationTool::InferPeptidoform,
                        bin: cfg.binary_path.clone(),
                        arg_lists: pf_args,
                        display_name: "pyprophet infer peptidoform",
                    });
                }
            }
        }

        if cfg.run_export {
            if cfg.export_tsv {
                self.pipeline.push(Stage {
                    tool: ValidationTool::ExportTsv,
                    bin: cfg.binary_path.clone(),
                    arg_lists: build_export_tsv_args(cfg),
                    display_name: "pyprophet export tsv",
                });
            }
            
            let prec_args = build_export_matrix_precursor_args(cfg);
            if !prec_args.is_empty() {
                self.pipeline.push(Stage {
                    tool: ValidationTool::ExportMatrixPrecursor,
                    bin: cfg.binary_path.clone(),
                    arg_lists: prec_args,
                    display_name: "pyprophet export matrix (precursor)",
                });
            }
            let pep_args = build_export_matrix_peptide_args(cfg);
            if !pep_args.is_empty() {
                self.pipeline.push(Stage {
                    tool: ValidationTool::ExportMatrixPeptide,
                    bin: cfg.binary_path.clone(),
                    arg_lists: pep_args,
                    display_name: "pyprophet export matrix (peptide)",
                });
            }
            let prot_args = build_export_matrix_protein_args(cfg);
            if !prot_args.is_empty() {
                self.pipeline.push(Stage {
                    tool: ValidationTool::ExportMatrixProtein,
                    bin: cfg.binary_path.clone(),
                    arg_lists: prot_args,
                    display_name: "pyprophet export matrix (protein)",
                });
            }
            if cfg.export_parquet {
                self.pipeline.push(Stage {
                    tool: ValidationTool::ExportTsv,
                    bin: cfg.binary_path.clone(),
                    arg_lists: build_export_parquet_args(cfg),
                    display_name: "pyprophet export parquet",
                });
            }
        }

        self.overall_total_runs = self.pipeline.iter().map(|s| s.arg_lists.len()).sum();

        if self.pipeline.is_empty() {
            self.output_lines.lock().unwrap().push("No steps to run.".into());
        } else {
            self.start_next_stage();
        }
    }

    fn start_next_stage(&mut self) {
        if let Some(stage) = self.pipeline.get(self.current_stage_index).cloned() {
            // --- Reset per‐stage metrics BEFORE spawning ---
            self.total_number_of_runs   = stage.arg_lists.len();
            self.processing_run_n       = 0;
            self.stage_successful_runs  = 0;
            self.stage_failed_runs      = 0;
            self.processing_time.clear();
            self.avg_processing_time    = Duration::ZERO;
            self.total_processing_time  = Duration::ZERO;
            self.processed_fraction     = 0.0;
            // And clear any leftover cancel_flag
            self.cancel_flag.store(false, Ordering::Relaxed);

            // spawn_batch_runner is exactly the same helper you used before:
            let (tx, rx) = channel();
            let handle = spawn_batch_runner(
                stage.bin.clone(),
                stage.arg_lists.clone(),
                Arc::clone(&self.output_lines),
                tx,
                Arc::clone(&self.cancel_flag),
                Arc::clone(&self.current_children),
                self.n_concurrent_processes, 
                None
            );
            self.current_stage_runner = Some(StageRunner { handle, progress_rx: rx, stage });
            // log start…
        } else {
            // finished all stages
            self.output_lines.lock().unwrap()
                .push(format!("All validation steps complete ({} total steps).", self.pipeline.len()));
            self.current_stage_runner = None;
        }
    }

    fn abort_pipeline(&mut self) {
        self.cancel_flag.store(true, Ordering::Relaxed);
        // kill every running child
        for child in self.current_children.lock().unwrap().drain(..) {
            let _ = child.lock().unwrap().kill();
        }
        self.current_stage_runner = None;
        self.pipeline.clear();
        self.output_lines.lock().unwrap().push("Validation cancelled.".into());
    }

    fn console_ui(&mut self, ui: &mut Ui, config: &Input) {
        let is_running = self.current_stage_runner.is_some();

        if let Some(pyp_cfg) = config.statistical_validation.as_ref() {
            ui.horizontal(|ui| {
                // ▶︎ Run Workflow
                if ui.add_enabled(!is_running, Button::new("▶︎ Run Workflow")).clicked() {
                    self.initial_start_time = Instant::now();
                    self.start_pipeline(pyp_cfg);
                }

                // ❔ Show Help
                if ui.add_enabled(!is_running, Button::new("❔ Show Help")).clicked() {
                    if pyp_cfg.run_merge {
                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["merge".into(), "osw".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }

                    if pyp_cfg.run_score {
                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["score".into(), "--helphelp".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }

                    if pyp_cfg.run_inference {
                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["infer".into(), "peptide".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());

                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["infer".into(), "protein".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());

                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["infer".into(), "gene".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());

                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["infer".into(), "peptidoform".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }

                    if pyp_cfg.run_export {
                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["export".into(), "tsv".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());

                        let (_h,_c)=spawn_process_runner(pyp_cfg.binary_path.clone(), vec!["export".into(), "matrix".into(), "--help".into()], Arc::clone(&self.output_lines), None);
                        self.output_lines.lock().unwrap().push("\n".into());
                    }
                }

                // ⏹ Stop
                if ui.add_enabled(is_running, Button::new("⏹ Stop")).clicked() {
                    self.abort_pipeline();
                }

                // 🗑 Clear
                if ui.add_enabled(!is_running, Button::new("🗑 Clear")).clicked() {
                    self.output_lines.lock().unwrap().clear();
                }

                // 💾 Save
                if ui.add_enabled(!is_running, Button::new("💾 Save")).clicked() {
                    let ts = Local::now().format("%Y%m%d_%H%M%S");
                    let default_name = format!("pyprophet_workflow_log_{}.txt", ts);
                    if let Some(path) = FileDialog::new()
                        .set_title("Save Console Log")
                        .set_file_name(&default_name)
                        .save_file()
                    {
                        if let Ok(mut f) = File::create(&path) {
                            for line in self.output_lines.lock().unwrap().iter() {
                                let _ = writeln!(f, "{}", line);
                            }
                        }
                    }
                }
            });

            ui.separator();

            // — Status / progress bars —
            if !self.pipeline.is_empty() {
                // Current stage info
                if let Some(runner) = &self.current_stage_runner {
                    ui.label(format!(
                        "Stage {}/{}: {}",
                        self.current_stage_index + 1,
                        self.pipeline.len(),
                        runner.stage.display_name
                    ));

                    ui.label(format!(
                        "Run {}/{}  (✔ {}  ✖ {})",
                        self.processing_run_n,
                        self.total_number_of_runs,
                        self.stage_successful_runs,
                        self.stage_failed_runs
                    ));

                    // // Stage progress bar
                    // ui.add(
                    //     egui::ProgressBar::new(self.processed_fraction)
                    //         .text(format!("{:.0}%", self.processed_fraction * 100.0))
                    //         .animate(true),
                    // );
                }

                // Overall progress
                let overall_frac = if self.overall_total_runs > 0 {
                    (self.successful_runs + self.failed_runs) as f32
                        / (self.overall_total_runs as f32)
                } else {
                    0.0
                };
                ui.label(format!(
                    "Overall {}/{}  (✔ {}  ✖ {})",
                    self.successful_runs + self.failed_runs,
                    self.overall_total_runs,
                    self.successful_runs,
                    self.failed_runs
                ));
                // ui.add(
                //     egui::ProgressBar::new(overall_frac)
                //         .text(format!("{:.0}%", overall_frac * 100.0))
                //         .animate(true),
                // );
            }

            // — Console Output —
            ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
                for line in self.output_lines.lock().unwrap().iter() {
                    ui.label(line);
                }
            });

            // — Process progress messages for current stage —
            let mut need_repaint = false;
            if let Some(runner) = &self.current_stage_runner {
                for msg in runner.progress_rx.try_iter() {
                    match msg {
                        Progress::Started { index, total } => {
                            self.processing_run_n = index;
                            self.total_number_of_runs = total;
                        }
                        Progress::Finished { duration, .. } => {
                            self.stage_successful_runs += 1;
                            self.successful_runs += 1;
                            self.processing_time.push(duration);
                            let sum: Duration = self.processing_time.iter().sum();
                            // self.total_processing_time = sum;
                            self.total_processing_time = self.initial_start_time.elapsed();
                            self.avg_processing_time = sum / (self.processing_time.len() as u32);
                        }
                        Progress::Failed { index, code } => {
                            self.stage_failed_runs += 1;
                            self.failed_runs += 1;
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(format!(
                                "[{}] Stage {} run {}/{} failed (exit {})",
                                ts,
                                runner.stage.display_name,
                                index,
                                self.total_number_of_runs,
                                code
                            ));
                        }
                        Progress::Cancelled { index } => {
                            self.stage_failed_runs += 1;
                            self.failed_runs += 1;
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(format!(
                                "[{}] Stage {} run {}/{} cancelled",
                                ts,
                                runner.stage.display_name,
                                index,
                                self.total_number_of_runs
                            ));
                            self.abort_pipeline();
                            break;
                        }
                        Progress::AllDone => {
                            let ts = Local::now().format("%Y-%m-%d %H:%M:%S");
                            self.output_lines.lock().unwrap().push(format!(
                                "[{}] Stage {} complete: {} success, {} failed.",
                                ts,
                                runner.stage.display_name,
                                self.stage_successful_runs,
                                self.stage_failed_runs
                            ));
                            // Advance
                            self.current_stage_index += 1;
                            self.current_stage_runner = None;
                            self.start_next_stage();
                            break;
                        }
                    }
                    // update per-stage fraction
                    let done = (self.stage_successful_runs + self.stage_failed_runs) as f32;
                    self.processed_fraction = if self.total_number_of_runs > 0 {
                        done / self.total_number_of_runs as f32
                    } else {
                        0.0
                    };
                }
                need_repaint = true;
            }

            if need_repaint {
                ui.ctx().request_repaint();
            }
        }
    }

    fn results_ui(&mut self, ui: &mut Ui, config: &Input) {
        ui.label("Results tab is not yet implemented.");
    }
}