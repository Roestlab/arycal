use std::collections::{BTreeMap, HashMap};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::{Path, PathBuf};

use arycal_cloudpath::osw::{OswAccess, PrecursorPeakBoundaries};
use arycal_cloudpath::util::extract_basename;
use eframe::egui::{Ui, WidgetText, CollapsingHeader, ScrollArea};
use egui::{Sense, CursorIcon};
use egui_plot::{Plot, Line, PlotPoints, Legend, Corner};
use egui::Color32;
use egui_tiles::{Behavior, TileId, Tile, Tiles, Tree, UiResponse, SimplificationOptions};


use arycal_cli::input::Input;
use arycal_cloudpath::{ChromatogramReader, sqmass::SqMassAccess, xic_parquet::DuckDBParquetChromatogramReader};
use arycal_cloudpath::sqmass::TransitionGroup;
use arycal_common::config::{PlotMode, VisualizationConfig, XicFileType};



// ----------------------------------------------------------------
// 1) The Behavior that drives our tabs:
//
// – we force every pane to have a tab bar by setting
//   `all_panes_must_have_tabs = true`
// – we use the PathBuf list to name each tab
// 1) Our Behavior: gives each pane a tab (even if it's the only tab),
//    names it after the file‐basename, and draws the plot in pane_ui:
struct TabBehavior<'a> {
    paths:      &'a [PathBuf],
    groups:     &'a [TransitionGroup],
    boundaries: &'a [Vec<PrecursorPeakBoundaries>],
    native_id_annotation_map: &'a HashMap<String, String>,
    viz_cfg:    &'a VisualizationConfig,
    selected_boundaries: &'a std::collections::HashSet<(String, i64)>,
    opts:       SimplificationOptions,
}

impl<'a> TabBehavior<'a> {
    fn new(
        paths:   &'a [PathBuf],
        groups:  &'a [TransitionGroup],
        boundaries: &'a [Vec<PrecursorPeakBoundaries>],
        native_id_annotation_map: &'a HashMap<String, String>,
        viz_cfg: &'a VisualizationConfig,
        selected_boundaries: &'a std::collections::HashSet<(String, i64)>,
    ) -> Self {
        let mut opts = SimplificationOptions::default();
        opts.all_panes_must_have_tabs = true;
        opts.join_nested_linear_containers = false;
        Self { paths, groups, boundaries, native_id_annotation_map, viz_cfg, selected_boundaries, opts }
    }
}

impl<'a> Behavior<usize> for TabBehavior<'a> {
    fn tab_title_for_pane(&mut self, pane: &usize) -> WidgetText {
        let name = self.paths[*pane]
            .file_name()
            .and_then(|os| os.to_str())
            .unwrap_or("…");
        name.into()
    }

    fn pane_ui(
        &mut self,
        ui: &mut egui::Ui,
        _tile_id: TileId,
        pane_idx: &mut usize,
    ) -> UiResponse {
        // Just draw your plot here — no custom dragging:
        let group   = &self.groups[*pane_idx];
        let boundaries = &self.boundaries[*pane_idx];
        let native_id_annotation_map = self.native_id_annotation_map;
        let viz_cfg = self.viz_cfg;
        VisualizationState::render_plot_inner(ui, group, native_id_annotation_map, viz_cfg, &boundaries, self.selected_boundaries);
        UiResponse::None
    }

    fn simplification_options(&self) -> SimplificationOptions {
        self.opts.clone() // or copy as needed
    }
}




pub struct VisualizationState {
    /// Generic chromatogram reader (sqMass or Parquet)
    readers: Vec<Option<Box<dyn ChromatogramReader>>>,
    pub run_basenames: Vec<String>,

    show_plot_windows: Vec<bool>,
    /// Loaded transition groups
    groups: Vec<TransitionGroup>,
    /// Loaded peak boundaries for the current peptide+charge
    group_boundaries: Vec<Vec<PrecursorPeakBoundaries>>,
    /// Index of the currently selected group
    selected_group: Option<usize>,
    pub boundaries_initialized: bool,
    pub selected_boundaries: std::collections::HashSet<(String, i64)>,
    /// Remember which feature path we loaded from
    last_feature_path: Option<std::path::PathBuf>,
    /// Remember which peptide+charge we last loaded (to avoid re‐loading every frame).
    last_peptide_charge: Option<(String, usize)>,
    /// A map of native IDs to annotations, e.g. { "native_id": "annotation" }
    pub native_id_annotation_map: Option<HashMap<String, String>>,
    /// Optional precursor table for a map of modified sequence and charge state
    pub precursor_table: Option<BTreeMap<String, Vec<u8>>>,
    /// Optional precursor boundaries
    pub precursor_boundaries: Option<Vec<PrecursorPeakBoundaries>>,

    /// How long the precursor‐table load took
    pub precursor_load_time: std::time::Duration,
    /// How many distinct peptides we saw
    pub num_unique_peptides: usize,
    /// How many peptide–charge precursors total
    pub num_unique_precursors: usize,

    /// timing metrics (in seconds)
    pub xic_load_time:     std::time::Duration,
    pub plot_render_time:  std::time::Duration,

    /// An optional tile‐tree so we don’t re-create every frame:
    tiles: Option<Tree<usize>>,
    /// How many panes we had last time we built the tree
    last_pane_count: Option<usize>,
    /// What grid (rows, cols) we last used
    last_grid: (usize, usize),
}

impl VisualizationState {

    pub fn default() -> Self {
        Self {
            readers: Vec::new(),
            run_basenames: Vec::new(),
            show_plot_windows: Vec::new(),
            groups: Vec::new(),
            group_boundaries: Vec::new(),
            selected_group: None,
            boundaries_initialized: false,
            selected_boundaries: std::collections::HashSet::new(),
            last_feature_path: None,
            last_peptide_charge: None,
            native_id_annotation_map: None,
            precursor_table: None,
            precursor_boundaries: None,
            precursor_load_time: std::time::Duration::ZERO,
            num_unique_peptides: 0,
            num_unique_precursors: 0,
            xic_load_time: std::time::Duration::ZERO,
            plot_render_time: std::time::Duration::ZERO,
            tiles: None,
            last_pane_count: None,
            last_grid: (0, 0), 
        }
    }

    pub fn ui(&mut self, ui: &mut Ui, config: &mut Input) {

        if config.xic.file_paths.is_empty() {
            ui.label("To visualize extracted ion chromatograms, please add XIC files in the configuration panel. Acceptable file formats are sqMass (sqlite-based output from openswathworkflow) or Parquet (produced by pyprophet's `export parquet`).");
            return;
        }
        // reload precursor‐table, reset on feature‐file change…
        self.sync_precursor_table(config);

        // ensure our per-file slots match whatever the user has dropped in:
        self.sync_file_slots(config);

        // if (peptide,charge) chosen, fetch & cache each file’s group
        self.update_chromatogram_groups(config);

        // finally: render one window per file
        let render_start = std::time::Instant::now();
        if let Some(viz_cfg) = &config.visualization {
            match viz_cfg.plot_mode {
                PlotMode::Floating => {
                    self.show_floating_windows(ui, config);
                }
                PlotMode::EmbeddedGrid => {
                    self.show_embedded_with_tabs(ui, config);
                }
            }
        }
        self.plot_render_time = render_start.elapsed();

    }

    fn sync_precursor_table(&mut self, config: &mut Input) {
        // 1) Detect changes in the feature‐file selection
        let current = config.features.file_paths.get(0);
        if current != self.last_feature_path.as_ref() {
            self.precursor_table = None;
            self.last_feature_path = current.cloned();
            self.groups.clear();
            self.selected_group = None;

            if let Some(viz_cfg) = &mut config.visualization {
                viz_cfg.selected_peptide = None;
                viz_cfg.selected_charge  = None;
            }
        }

        // 2) Lazily load the precursor table (peptide→charge map)
        if self.precursor_table.is_none() {
            if let Some(path) = current {
                let _ = self.load_precursor_table(path);
            }
        }
    }

    /// Ensure our per-file vectors (`readers`, `show_plot_windows`, `groups`) match
    /// the number of configured XIC files, opening plot windows only if we already
    /// have a peptide+charge selection.
    fn sync_file_slots(&mut self, config: &Input) {
        let file_count = config.xic.file_paths.len();

        self.run_basenames = config.xic.file_paths
            .iter()
            .map(|p| extract_basename(p))
            .collect::<Vec<_>>();
        
        // Run extract_basename on self.run_basenames again encase the original file had multiple extensions
        self.run_basenames = self.run_basenames
            .iter()
            .map(|b| extract_basename(b))
            .collect::<Vec<_>>();

        println!("Syncing file slots: {} files", file_count);
        println!("Run basenames: {:?}", self.run_basenames);
    
        // should our windows be open *right now*?
        let plots_open = config
            .visualization
            .as_ref()
            .map(|v| v.selected_peptide.is_some() && v.selected_charge.is_some())
            .unwrap_or(false);
    
        // 1) Resize readers & groups to exactly file_count
        self.readers
            .resize_with(file_count, || None);
        self.groups
            .resize_with(file_count, || TransitionGroup::new(String::new()));
        self.group_boundaries
            .resize_with(file_count, || Vec::new());

    
        // 2) Always reset ALL the window‐open flags
        self.show_plot_windows.clear();
        self.show_plot_windows.resize(file_count, plots_open);
    }

    /// If user has selected (pep,charge) *and* there are XIC files,
    ///    make sure each file has a reader, fetch native IDs, then call read_chromatograms.
    fn update_chromatogram_groups(&mut self, config: &Input) {
        // 1) Do we even have a peptide+charge?
        let viz_cfg = match &config.visualization {
            Some(v) if v.selected_peptide.is_some() && v.selected_charge.is_some() => v,
            _ => {
                // no valid selection → clear our “last” so next time someone picks one, we do load
                self.last_peptide_charge = None;
                self.native_id_annotation_map = None;
                return;
            }
        };
        let pep = viz_cfg.selected_peptide.as_ref().unwrap();
        let charge = viz_cfg.selected_charge.unwrap();

        // 2) If it’s the *same* as last time, skip.
        if self
            .last_peptide_charge
            .as_ref()
            .map(|(p, c)| p == pep && *c == charge)
            .unwrap_or(false)
        {
            return;
        }

        // Load peak boundaries if not already loaded'
        self.precursor_boundaries = if let Some(path) = config.features.file_paths.get(0) {
            match self.fetch_precursor_peak_boundaries(path, pep, charge as i32) {
                Ok(boundaries) => Some(boundaries),
                Err(e) => {
                    eprintln!("Error fetching precursor peak boundaries: {}", e);
                    None
                }
            }
        } else {
            None
        };
        self.selected_boundaries.clear();
        self.boundaries_initialized = false;

        // Filter self,precursor_boundaries to only those that match the current self.run_basenames
        if self.precursor_boundaries.is_some() {
            self.precursor_boundaries = self.precursor_boundaries.as_mut().map(|boundaries| {
                boundaries.iter()
                    .filter(|b| self.run_basenames.contains(&b.run_filename))
                    .cloned()
                    .collect()
            });
        }

        println!("precursor_boundaries: {:?}", self.precursor_boundaries);

        // 3) Otherwise update our “last” and actually do the work:
        self.last_peptide_charge = Some((pep.clone(), charge));

        for (i, path) in config.xic.file_paths.iter().enumerate() {
            // ensure reader slot exists
            self.ensure_reader_for_file(i, config);

            // get file basename
            let run_basename = extract_basename(path);
            // run extract_basename again to handle multiple extensions
            // TODO: Need to implement a more robust way to handle multiple extensions
            let run_basename = extract_basename(&run_basename);

            // Filter peak boundaries for this run
            self.group_boundaries[i] = if let Some(boundaries) = &self.precursor_boundaries {
                println!("Filtering boundaries for run: {}", run_basename);
                boundaries.iter().filter(|b| b.run_filename == run_basename).cloned().collect()
            } else {
                Vec::new()
            };            

            // fetch native IDs
            let native_id_annotation_map = self
                .fetch_native_ids(
                    &config.features.file_paths[0],
                    pep,
                    charge as i32,
                    config.xic.include_precursor,
                    config.xic.num_isotopes,
                    config.filters.include_identifying_transitions.unwrap_or(false),
                )
                .unwrap_or_default();

            // store the native_id_annotation_map for later use
            self.native_id_annotation_map = Some(native_id_annotation_map.clone());

            let raw_ids = native_id_annotation_map
                .into_iter()
                .map(|(id, _)| id)
                .collect::<Vec<String>>();

            let native_ids: Vec<&str> = raw_ids.iter().map(String::as_str).collect();

            println!(
                "Fetching data for peptide {} (charge {}) on file #{} → {:?}",
                pep, charge, i, raw_ids
            );

            // build a group_id that actually encodes peptide+charge so it changes
            let group_id = format!("{}-{}-file{}", pep, charge, i);

            // read chromatograms
            if let Some(reader) = &self.readers[i] {
                let start = std::time::Instant::now();
                if let Ok(group) =
                    reader.read_chromatograms("NATIVE_ID", native_ids.clone(), group_id.clone())
                {
                    self.groups[i] = group;
                    self.show_plot_windows[i] = true;
                    self.xic_load_time = start.elapsed();
                }
            }
        }
    } 

    /// Lazily initialize `self.readers[i]` from config.xic.file_paths[i]
    fn ensure_reader_for_file(&mut self, i: usize, config: &Input) {
        if self.readers[i].is_none() {
            self.readers[i] = match config.xic.file_type.as_ref() {
                Some(XicFileType::SqMass) => {
                    config.xic.file_paths[i]
                        .to_str()
                        .and_then(|p| SqMassAccess::new(p).ok())
                        .map(|r| Box::new(r) as _)
                }
                Some(XicFileType::Parquet) => {
                    config.xic.file_paths[i]
                        .to_str()
                        .and_then(|p| DuckDBParquetChromatogramReader::new(p).ok())
                        .map(|r| Box::new(r) as _)
                }
                _ => None,
            };
        }
    }

    /// One floating window PER XIC file, each overlaying all its traces in color
    pub fn show_floating_windows(&mut self, ui: &mut egui::Ui, config: &Input) {
        // 1) Compute the top‐left of the central panel:
        let available = ui.available_rect_before_wrap();
        let offset = available.min;

        // 2) Grid layout math:
        let n    = config.xic.file_paths.len();
        let cols = (n as f32).sqrt().ceil() as usize;
        let dx   = 630.0;
        let dy   = 360.0;

        // 3) Clone out an owned VisualizationConfig so closures
        //    only borrow that and never need `&mut self`.
        let viz_cfg_owned = config
            .visualization
            .clone()
            .unwrap_or_else(VisualizationConfig::default);

        let native_id_annotation_map: HashMap<String,String> = self.native_id_annotation_map
            .clone()              
            .unwrap_or_default(); 

        for i in 0..n {
            let group = &self.groups[i];
            let boundaries = &self.group_boundaries[i];
            let open_flag = &mut self.show_plot_windows[i];
            let path      = &config.xic.file_paths[i];
            let title     = path
                .file_name()
                .and_then(|os| os.to_str())
                .unwrap_or("Chromatograms");

            let col = i % cols;
            let row = i / cols;
            let pos = [
                offset.x + (col as f32) * dx,
                offset.y + (row as f32) * dy,
            ];

            egui::Window::new(title)
                .default_pos(pos)
                // .open(open_flag)
                .resizable(true)
                .default_size([600.0, 400.0])
                .show(ui.ctx(), |ui| {
                    // call your helper with the entire config struct
                    Self::render_plot_inner(ui, group, &native_id_annotation_map, &viz_cfg_owned, &boundaries, &self.selected_boundaries);
                });
        }
    }
    
    /// Renders all loaded chromatogram groups in an N×M grid
    /// inside the main panel.
    pub fn show_embedded_with_tabs(&mut self, ui: &mut Ui, config: &mut Input) {
        let n = self.groups.len();
        if n == 0 { return; }
    
        let viz_cfg = config.visualization
            .get_or_insert_with(VisualizationConfig::default);
        let native_id_annotation_map: HashMap<String,String> = self.native_id_annotation_map
        .clone()              
        .unwrap_or_default(); 

        let rows = viz_cfg.grid_rows as usize;
        let cols = viz_cfg.grid_cols as usize;
        let tree_id = ui.id().with("embedded_tabs");
    
        // rebuild if we have no tree yet, OR pane count changed, OR grid shape changed:
        let rebuild = 
            self.tiles.is_none()
            || self.last_pane_count.map_or(true, |c| c != n)
            || self.last_grid != (rows, cols);

        if rebuild {
            // remember new shape
            self.last_pane_count = Some(n);
            self.last_grid = (rows, cols);

            // build a fresh Tiles<usize>
            let mut tiles = Tiles::<usize>::default();
            let cell_count = rows * cols;

            // distribute pane-indices 0..n into `cell_count` buckets as evenly as possible
            let mut buckets = vec![Vec::new(); cell_count];
            for pane_idx in 0..n {
                let bucket = pane_idx * cell_count / n;
                buckets[bucket].push(pane_idx);
            }

            // for each bucket, insert whatever panes landed there into one tab-tile
            let mut cell_ids = Vec::with_capacity(cell_count);
            for bucket in buckets {
                if bucket.is_empty() {
                    // empty cell → empty tab container
                    cell_ids.push(tiles.insert_tab_tile(Vec::new()));
                } else {
                    // one pane‐tile PER pane_idx in this bucket
                    let pane_tiles: Vec<_> =
                        bucket.into_iter().map(|idx| tiles.insert_pane(idx)).collect();
                    cell_ids.push(tiles.insert_tab_tile(pane_tiles));
                }
            }

            let root = tiles.insert_grid_tile(cell_ids);
            self.tiles = Some(Tree::new(tree_id, root, tiles));
        }
    
        if let Some(tree) = &mut self.tiles {
            // build a behavior that owns references into `self.groups` and `viz_cfg`
            let mut behavior = TabBehavior::new(
                &config.xic.file_paths,
                &self.groups,
                &self.group_boundaries,
                &native_id_annotation_map,
                viz_cfg,
                &self.selected_boundaries,
            );
            // ONE CALL to tree.ui will both layout the grid+tabs and invoke pane_ui for each pane:
            tree.ui(&mut behavior, ui);
        }
    }
      

    /// Helper to draw a single Plot inside whatever UI you give it.
    fn render_plot_inner(
        ui: &mut egui::Ui,
        group: &TransitionGroup,
        native_id_annotation_map: &HashMap<String, String>,
        viz_cfg: &VisualizationConfig,
        boundaries: &[PrecursorPeakBoundaries],
        selected_boundaries: &std::collections::HashSet<(String, i64)>,
    ) {
        // start the builder as before
        let mut plot = Plot::new(group.group_id.clone())
            // .height(350.0)
            .show_background(viz_cfg.show_background)
            .show_grid(viz_cfg.show_grid);

        // legend, axis‐labels, etc...
        if viz_cfg.show_legend {
            plot = plot.legend(Legend::default().position(Corner::RightTop));
        }
        if viz_cfg.show_axis_labels {
            plot = plot
                .x_axis_label("Retention Time (sec)")
                .y_axis_label("Intensity");
        }

        // 1) link axis?
        if viz_cfg.link_axis_x || viz_cfg.link_axis_y {
            // use a shared Id so all plots with this flag join the same group
            let axis_group = egui::Id::new("shared_axis_group");
            plot = plot.link_axis(
                axis_group,
                [viz_cfg.link_axis_x, viz_cfg.link_axis_y],
            );
        }

        // 2) link cursor?
        if viz_cfg.link_cursor {
            let cursor_group = egui::Id::new("shared_cursor_group");
            // [true,true] will share both x and y cursor,
            // but you can choose [true,false], etc.
            plot = plot.link_cursor(cursor_group, [true, true]);
        }

        // finally render
        plot.show(ui, |plot_ui| {
            for (native_id, chrom) in &group.chromatograms {
                let chrom = if viz_cfg.smoothing_enabled {
                    // apply Savitzky–Golay smoothing
                    chrom.smooth_sgolay(viz_cfg.sgolay_window, viz_cfg.sgolay_order).unwrap()
                } else {
                    chrom.clone()
                };

                let mut hasher = std::collections::hash_map::DefaultHasher::new();
                native_id.hash(&mut hasher);
                let hue = (hasher.finish() % 360) as f32;
                let color = Self::hsl_to_rgb(hue, 0.7, 0.5);

                let pts: egui_plot::PlotPoints = chrom
                    .retention_times
                    .iter()
                    .zip(&chrom.intensities)
                    .map(|(&rt, &i)| [rt, i])
                    .collect();

                // Get annotation for this native_id, if available
                let annotation = native_id_annotation_map
                    .get(native_id)
                    .map(|a| a.to_string())
                    .unwrap_or_else(|| native_id.clone());

                plot_ui.line(
                    egui_plot::Line::new(annotation.clone(), pts)
                        .color(color)
                        .width(2.0),
                );
            }

            // Overlay peak boundaries as vertical dashed lines
            if !boundaries.is_empty() {

                let max_y = group.chromatograms
                    .values()
                    .flat_map(|c| c.intensities.iter())
                    .cloned()
                    .fold(0.0, f64::max);
            
                for b in boundaries {
                    let key = (b.run_filename.clone(), b.feature_id);

                    if selected_boundaries.contains(&key) {
                        let left = b.left_width;
                        let right = b.right_width;

                        println!("Drawing boundary for run {} - feature {}: left={}, right={}", b.run_filename, b.feature_id, left, right);
                
                        // Generate a distinct color per feature
                        let mut hasher = std::collections::hash_map::DefaultHasher::new();
                        b.sorted_feature_id.hash(&mut hasher);
                        let hue = (hasher.finish() % 360) as f32;
                        let color = VisualizationState::hsl_to_rgb(hue, 0.8, 0.5);
                
                        // Build legend label
                        let label = format!("F{}", b.sorted_feature_id);
                
                        // // Combine left and right boundaries into one line with a gap
                        // let points = vec![
                        //     [left, 0.0], [left, max_y],
                        //     [f64::NAN, f64::NAN], // break between segments
                        //     [right, 0.0], [right, max_y],
                        // ];
                
                        // let boundary_line = Line::new(
                        //         label.clone(),
                        //         PlotPoints::from(points)
                        //     )
                        //     .color(color)
                        //     .style(egui_plot::LineStyle::dashed_loose())
                        //     .width(1.5);
                
                        // plot_ui.line(boundary_line);
                        let left_line = Line::new(
                            label.clone(),
                            PlotPoints::from(vec![[left, 0.0], [left, max_y]])
                        )
                        .color(color)
                        .style(egui_plot::LineStyle::dashed_loose())
                        .width(1.5);
                        
                        let right_line = Line::new(
                            label.clone(),  // same label
                            PlotPoints::from(vec![[right, 0.0], [right, max_y]])
                        )
                        .color(color)
                        .style(egui_plot::LineStyle::dashed_loose())
                        .width(1.5);
                        // .show_legend(false);  // <-- hide this duplicate in legend
                        
                        plot_ui.line(left_line);
                        plot_ui.line(right_line);
                        
                    }
                } // for loop
            }
            
        });
    }
    

    // Simple HSL→RGB converter (h in [0,360), s,l in [0,1])
    fn hsl_to_rgb(h: f32, s: f32, l: f32) -> Color32 {
        let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
        let h1 = h / 60.0;
        let x = c * (1.0 - (h1 % 2.0 - 1.0).abs());
        let (r1, g1, b1) = match h1 {
            h if h < 1.0 => (c, x, 0.0),
            h if h < 2.0 => (x, c, 0.0),
            h if h < 3.0 => (0.0, c, x),
            h if h < 4.0 => (0.0, x, c),
            h if h < 5.0 => (x, 0.0, c),
            _            => (c, 0.0, x),
        };
        let m = l - c / 2.0;
        let to_u8 = |v: f32| ((v + m).clamp(0.0, 1.0) * 255.0) as u8;
        Color32::from_rgb(to_u8(r1), to_u8(g1), to_u8(b1))
    }

    fn load_precursor_table(&mut self, features_path: &Path) -> anyhow::Result<()> {
        let start = std::time::Instant::now();
    
        let reader = OswAccess::new(features_path.to_str().unwrap(), true)?;
        let table = reader.fetch_full_peptide_precursor_table(true)?;
    
        // compute counts
        let peptide_count = table.len();
        let precursor_count: usize = table.values().map(|v| v.len()).sum();
    
        // store into state
        self.precursor_table = Some(table);
        self.precursor_load_time = start.elapsed();
        self.num_unique_peptides = peptide_count;
        self.num_unique_precursors = precursor_count;
    
        Ok(())
    }

    fn fetch_native_ids(&mut self, features_path: &Path, peptide: &str, charge: i32, include_precursor: bool, max_number_of_isotopes: usize, include_identifying_transitions: bool) -> anyhow::Result<HashMap<String, String>> {
        let reader = OswAccess::new(features_path.to_str().unwrap(), true)?;
        let native_ids = reader.fetch_native_ids(peptide, charge, include_precursor, max_number_of_isotopes, include_identifying_transitions);

        match native_ids {
            Ok(ids) => Ok(ids),
            Err(e) => Err(anyhow::anyhow!("Failed to fetch native IDs: {}", e)),
        }
    }
    
    fn fetch_precursor_peak_boundaries(&mut self, features_path: &Path, peptide: &str, charge: i32) -> anyhow::Result<Vec<PrecursorPeakBoundaries>> {
        let reader = OswAccess::new(features_path.to_str().unwrap(), true)?;
        let boundaries = reader.get_precursor_peak_boundaries(peptide, charge);

        match boundaries {
            Ok(b) => Ok(b),
            Err(e) => Err(anyhow::anyhow!("Failed to fetch precursor peak boundaries: {}", e)),
        }
    }

}
