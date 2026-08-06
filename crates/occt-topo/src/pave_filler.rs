//! `BOPAlgo_PaveFiller` — the Intersection phase of the Boolean Operations
//! algorithm (Phase 19).
//!
//! Source: `BOPAlgo_PaveFiller.hxx/.cxx` + `BOPAlgo_PaveFiller_{1..11}.cxx`
//! (TKBO/BOPAlgo).
//!
//! The class performs the pairwise intersection of the sub-shapes of the
//! arguments in the fixed order:
//!
//! 1. Vertex/Vertex;
//! 2. Vertex/Edge;
//! 3. Edge/Edge;
//! 4. Vertex/Face;
//! 5. Edge/Face;
//! 6. Face/Face.
//!
//! The results of the intersections are stored into the data structure
//! ([`BopdsDS`]) of the algorithm, which later feeds the building phase.
//!
//! This module owns the main class and the pipeline orchestration. The heavy
//! per-pair intersection work lives in the sibling Phase-19 modules
//! (`crate::pave_intersect`, `crate::pave_blocks`, `crate::pave_common`);
//! the pipeline steps delegate to them (and to the pave-common helpers) via
//! [`PaveFillerLike`].

use crate::abs::ShapeType;
use crate::bopds::BopdsDS;
use crate::int_tools_full::IntToolsContext;
use crate::pave_blocks::PaveFillerLike;
use crate::shape::TopoShape;

// ---------------------------------------------------------------------------
// GlueEnum
// ---------------------------------------------------------------------------

/// Gluing option of the algorithm. Source: `BOPAlgo_GlueEnum.hxx`.
///
/// Gluing trades robustness for speed on special inputs in which many
/// sub-shapes coincide:
///
/// - [`GlueEnum::None`] — no gluing (`BOPAlgo_GlueOff`), full general case;
/// - [`GlueEnum::Shift`] — glue coincident faces only after their vertices
///   were shifted together (`BOPAlgo_GlueShift`);
/// - [`GlueEnum::Full`] — treat all coincident faces as a single one, without
///   splitting them (`BOPAlgo_GlueFull`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlueEnum {
    /// No gluing — the general, most robust case.
    None,
    /// Glue shift mode.
    Shift,
    /// Glue full mode.
    Full,
}

impl Default for GlueEnum {
    fn default() -> Self {
        GlueEnum::None
    }
}

// ---------------------------------------------------------------------------
// PaveFiller
// ---------------------------------------------------------------------------

/// The Intersection phase of the Boolean Operations algorithm.
///
/// Owns the [`BopdsDS`] data structure (which accumulates every participating
/// shape, its sub-shapes and the results of their intersections), the shared
/// [`IntToolsContext`] geometry toolkit cache, and the error/warning report of
/// the run.
///
/// The public surface mirrors the OCCT class, reduced to what the Rust
/// pipeline needs. The fixed external interface (relied upon by the sibling
/// Phase-19 modules) is:
/// `new`, `set_arguments`, `perform`, `ds`, `ds_mut`, `has_errors`,
/// `add_error`.
#[derive(Debug, Clone)]
pub struct PaveFiller {
    /// The arguments of the operation.
    arguments: Vec<TopoShape>,
    /// The data structure of the algorithm.
    ds: BopdsDS,
    /// Cached intersection context (geometry tools + 2D classifiers).
    context: IntToolsContext,
    /// Additional tolerance for touching/coinciding detection, floored at
    /// `Precision::Confusion()` (1e-7).
    fuzzy_value: f64,
    /// Fatal alerts — a non-empty list means the algorithm has failed.
    errors: Vec<String>,
    /// Non-fatal alerts.
    warnings: Vec<String>,
    /// Gluing option of the algorithm.
    glue: GlueEnum,
    /// Non-destructive mode: the argument shapes are not modified.
    non_destructive: bool,
    /// Flag that the intersection must be repeated with increased vertex
    /// tolerances (set by the pipeline when new vertices were created).
    repeat_intersection: bool,
    /// Flag that the intersection phase completed successfully.
    intersection_done: bool,
}

impl Default for PaveFiller {
    fn default() -> Self {
        Self::new()
    }
}

impl PaveFiller {
    /// Empty constructor.
    ///
    /// `fuzzy_value` defaults to `Precision::Confusion()` (1e-7), glue is off,
    /// non-destructive mode is off, and the report is empty.
    pub fn new() -> Self {
        Self {
            arguments: Vec::new(),
            ds: BopdsDS::new(),
            context: IntToolsContext::new(),
            fuzzy_value: 1e-7,
            errors: Vec::new(),
            warnings: Vec::new(),
            glue: GlueEnum::None,
            non_destructive: false,
            repeat_intersection: false,
            intersection_done: false,
        }
    }

    // -----------------------------------------------------------------------
    // Options
    // -----------------------------------------------------------------------

    /// Sets the arguments of the operation.
    pub fn set_arguments(&mut self, shapes: &[TopoShape]) {
        self.arguments = shapes.to_vec();
    }

    /// Returns the arguments of the operation.
    pub fn arguments(&self) -> &[TopoShape] {
        &self.arguments
    }

    /// Sets the additional tolerance.
    ///
    /// OCCT clamps the value from below with `Precision::Confusion()`; a NaN
    /// is treated as the default confusion value (same policy as
    /// [`crate::bopalgo_options::BopAlgoOptions::set_fuzzy_value`]).
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.fuzzy_value = if v.is_nan() { 1e-7 } else { v.max(1e-7) };
    }

    /// Returns the additional tolerance.
    pub fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value
    }

    /// Sets the glue option of the algorithm.
    pub fn set_glue(&mut self, glue: GlueEnum) {
        self.glue = glue;
    }

    /// Returns the glue option of the algorithm.
    pub fn glue(&self) -> GlueEnum {
        self.glue
    }

    /// Sets the non-destructive mode: the argument shapes are not modified.
    pub fn set_non_destructive(&mut self, b: bool) {
        self.non_destructive = b;
    }

    /// Returns the non-destructive mode flag.
    pub fn non_destructive(&self) -> bool {
        self.non_destructive
    }

    // -----------------------------------------------------------------------
    // Alert report
    // -----------------------------------------------------------------------

    /// Clears the report (errors and warnings) and the data structure.
    ///
    /// The user-defined options (fuzzy value, glue, non-destructive mode,
    /// arguments) are *not* reset, matching `BOPAlgo_PaveFiller::Clear` +
    /// `BOPAlgo_Options::Clear`.
    pub fn clear(&mut self) {
        self.errors.clear();
        self.warnings.clear();
        self.ds.clear();
        self.context.clear_cached();
        self.intersection_done = false;
        self.repeat_intersection = false;
    }

    /// True if the algorithm has failed (at least one fatal alert).
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// True if the algorithm has generated at least one warning.
    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }

    /// Returns the collected fatal alerts.
    pub fn errors(&self) -> &[String] {
        &self.errors
    }

    /// Returns the collected non-fatal alerts.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Adds a fatal alert (the algorithm has failed).
    pub fn add_error(&mut self, msg: String) {
        self.errors.push(msg);
    }

    /// Adds a non-fatal alert.
    pub fn add_warning(&mut self, msg: String) {
        self.warnings.push(msg);
    }

    // -----------------------------------------------------------------------
    // Data structure access
    // -----------------------------------------------------------------------

    /// Returns the data structure of the algorithm.
    pub fn ds(&self) -> &BopdsDS {
        &self.ds
    }

    /// Returns the mutable data structure of the algorithm.
    pub fn ds_mut(&mut self) -> &mut BopdsDS {
        &mut self.ds
    }

    /// Returns the shared intersection context.
    pub fn context(&self) -> &IntToolsContext {
        &self.context
    }

    /// Returns the mutable intersection context.
    pub fn context_mut(&mut self) -> &mut IntToolsContext {
        &mut self.context
    }

    /// Returns the flag that the intersection must be repeated with increased
    /// vertex tolerances.
    pub fn repeat_intersection(&self) -> bool {
        self.repeat_intersection
    }

    /// Sets the flag that the intersection must be repeated.
    pub fn set_repeat_intersection(&mut self, b: bool) {
        self.repeat_intersection = b;
    }

    /// True if the intersection phase completed successfully.
    pub fn intersection_done(&self) -> bool {
        self.intersection_done
    }

    /// Sets the flag that the intersection phase completed.
    pub fn set_intersection_done(&mut self, b: bool) {
        self.intersection_done = b;
    }

    // -----------------------------------------------------------------------
    // Initialization
    // -----------------------------------------------------------------------

    /// Initializes the algorithm: appends every argument (and its whole
    /// sub-shape subtree) into the data structure and builds the per-argument
    /// index ranges (ranks).
    ///
    /// Mirrors `BOPAlgo_PaveFiller::Init`: the arguments must be non-empty and
    /// non-null; the report is cleared before the data structure is rebuilt.
    pub fn init(&mut self) -> Result<(), String> {
        if self.arguments.is_empty() {
            self.add_error("BOPAlgo_PaveFiller: too few arguments".to_string());
            return Err("BOPAlgo_PaveFiller: too few arguments".to_string());
        }
        for a in &self.arguments {
            if a.is_null() {
                self.add_error("BOPAlgo_PaveFiller: null input shape".to_string());
                return Err("BOPAlgo_PaveFiller: null input shape".to_string());
            }
        }
        self.clear();
        self.ds.set_arguments(self.arguments.clone());
        self.ds.init(&self.arguments);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Main entry points
    // -----------------------------------------------------------------------

    /// Performs the intersection of the sub-shapes of the arguments.
    pub fn perform(&mut self) -> Result<(), String> {
        self.perform_internal()
    }

    /// The pipeline: runs the intersection stages in the OCCT order,
    /// short-circuiting as soon as a stage reports an error.
    ///
    /// Source: `BOPAlgo_PaveFiller::PerformInternal`. The stages:
    ///
    /// `init` → `prepare` → `perform_vv` → `perform_ve` → SD-vertex update →
    /// `perform_ee` → SD-vertex update → `perform_vf` → SD-vertex update →
    /// `perform_ef` → SD-vertex update → `perform_ff` → `make_split_edges` →
    /// SD-vertex update → `make_blocks` → `make_pcurves`.
    ///
    /// The per-stage block *splitting* happens inside each intersection stage
    /// (the OCCT `SplitPaveBlocks`); this step only *redirects* the bound
    /// vertex indices of the blocks to their same-domain representatives
    /// (`BOPDS_DS::UpdatePaveBlocksWithSDVertices`).
    pub fn perform_internal(&mut self) -> Result<(), String> {
        self.init()?;
        self.check_errors()?;
        self.prepare()?;
        self.check_errors()?;
        self.perform_vv()?;
        self.check_errors()?;
        self.perform_ve()?;
        self.check_errors()?;
        self.ds.update_pave_blocks_with_sd_vertices();
        self.perform_ee()?;
        self.check_errors()?;
        self.ds.update_pave_blocks_with_sd_vertices();
        self.perform_vf()?;
        self.check_errors()?;
        self.ds.update_pave_blocks_with_sd_vertices();
        self.perform_ef()?;
        self.check_errors()?;
        self.ds.update_pave_blocks_with_sd_vertices();
        self.perform_ff()?;
        self.check_errors()?;
        // OCCT order: MakeSplitEdges (right after FF) precedes MakeBlocks;
        // MakePCurves follows MakeBlocks.
        self.make_split_edges()?;
        self.check_errors()?;
        self.ds.update_pave_blocks_with_sd_vertices();
        self.make_blocks()?;
        self.check_errors()?;
        self.make_pcurves()?;
        self.check_errors()?;
        self.intersection_done = true;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Pipeline stages (delegated to the sibling Phase-19 modules)
    // -----------------------------------------------------------------------

    /// Preparation stage. Source: `BOPAlgo_PaveFiller::Prepare`
    /// (`BOPAlgo_PaveFiller_7.cxx`).
    ///
    /// The main class keeps no prepare logic of its own: the pave-blocks of
    /// every source edge are initialized (so the later intersection stages see
    /// a splittable block per edge) and the shrunk-range data is computed by
    /// [`crate::pave_common::fill_shrunk_data`].
    fn prepare(&mut self) -> Result<(), String> {
        let n = self.ds().nb_source_shapes();
        for i in 0..n {
            if self.ds().shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Edge) {
                self.ds_mut().init_pave_blocks_for_edge(i);
            }
        }
        crate::pave_common::fill_shrunk_data(self)
    }

    /// Vertex/Vertex intersection. Source: `PerformVV`
    /// (`BOPAlgo_PaveFiller_1.cxx`). Delegates to
    /// [`crate::pave_intersect::perform_vv`].
    fn perform_vv(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_vv(self)
    }

    /// Vertex/Edge intersection. Source: `PerformVE`. Delegates to
    /// [`crate::pave_intersect::perform_ve`].
    fn perform_ve(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ve(self)
    }

    /// Edge/Edge intersection. Source: `PerformEE`. Delegates to
    /// [`crate::pave_intersect::perform_ee`].
    fn perform_ee(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ee(self)
    }

    /// Vertex/Face intersection. Source: `PerformVF`. Delegates to
    /// [`crate::pave_intersect::perform_vf`].
    fn perform_vf(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_vf(self)
    }

    /// Edge/Face intersection. Source: `PerformEF`. Delegates to
    /// [`crate::pave_intersect::perform_ef`].
    fn perform_ef(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ef(self)
    }

    /// Face/Face intersection. Source: `PerformFF`. Delegates to
    /// [`crate::pave_intersect::perform_ff`].
    fn perform_ff(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ff(self)
    }

    /// Groups coincident pave blocks into common blocks. Source: `MakeBlocks`
    /// (`BOPAlgo_PaveFiller_10.cxx`). Delegates to
    /// [`crate::pave_blocks::make_blocks`].
    fn make_blocks(&mut self) -> Result<(), String> {
        crate::pave_blocks::make_blocks(self)
    }

    /// Builds 2D curves of the section edges on the faces. Source: `MakePCurves`.
    /// Delegates to [`crate::pave_blocks::make_pcurves`].
    fn make_pcurves(&mut self) -> Result<(), String> {
        crate::pave_blocks::make_pcurves(self)
    }

    /// Creates the split edges (new edges of the common part). Source:
    /// `MakeSplitEdges`. Delegates to [`crate::pave_blocks::make_split_edges`].
    fn make_split_edges(&mut self) -> Result<(), String> {
        crate::pave_blocks::make_split_edges(self)
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Short-circuits the pipeline when a fatal alert was recorded by a stage.
    fn check_errors(&self) -> Result<(), String> {
        if self.has_errors() {
            let first = self
                .errors
                .first()
                .cloned()
                .unwrap_or_else(|| "BOPAlgo_PaveFiller: unknown error".to_string());
            return Err(first);
        }
        Ok(())
    }
}

/// The minimal contract [`crate::pave_blocks`] needs from its host. The
/// concrete [`PaveFiller`] satisfies it through the DS accessors and the alert
/// report; the p-curve policy uses the trait defaults (build p-curves, both
/// faces of a section).
impl PaveFillerLike for PaveFiller {
    fn ds(&self) -> &BopdsDS {
        &self.ds
    }
    fn ds_mut(&mut self) -> &mut BopdsDS {
        &mut self.ds
    }
    fn add_error(&mut self, msg: String) {
        self.errors.push(msg);
    }
    fn add_warning(&mut self, msg: String) {
        self.warnings.push(msg);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;

    #[test]
    fn defaults() {
        let pf = PaveFiller::new();
        assert!(pf.arguments().is_empty());
        assert_eq!(pf.fuzzy_value(), 1e-7, "fuzzy defaults to Precision::Confusion");
        assert_eq!(pf.glue(), GlueEnum::None, "glue off by default");
        assert!(!pf.non_destructive());
        assert!(!pf.has_errors());
        assert!(!pf.has_warnings());
        assert!(!pf.intersection_done());
        assert_eq!(pf.ds().nb_shapes(), 0);
    }

    #[test]
    fn init_appends_all_subshapes_and_sets_ranges() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        pf.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.init().unwrap();
        // Each unit box contributes 28 shapes (8 V + 12 E + 6 F + 1 Shell + 1 Solid).
        assert_eq!(pf.ds().nb_shapes(), 56);
        assert_eq!(pf.ds().nb_ranges(), 2);
        // A sub-shape of the first argument has rank 0, of the second rank 1.
        let e0 = pf.ds().index(&a.edges[0].0).expect("first edge indexed");
        let e1 = pf.ds().index(&c.edges[0].0).expect("second edge indexed");
        assert_eq!(pf.ds().rank(e0), 0);
        assert_eq!(pf.ds().rank(e1), 1);
        // A sub-shape round-trips through index()/shape().
        let v = pf.ds().index(&a.vertices[0].0).unwrap();
        assert!(pf.ds().shape(v).unwrap().same_tshape(&a.vertices[0].0));
    }

    #[test]
    fn perform_runs_full_pipeline_without_errors() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        pf.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.perform().unwrap();
        assert!(!pf.has_errors(), "errors: {:?}", pf.errors());
        assert!(pf.intersection_done());
        // The intersection created new shapes (SD vertices, crossing vertices,
        // section edges) on top of the 2×28 argument sub-shapes.
        assert!(pf.ds().nb_shapes() >= 56);
    }

    #[test]
    fn errors_and_warnings_accumulate() {
        let mut pf = PaveFiller::new();
        assert!(!pf.has_errors());
        pf.add_error("intersection failed".to_string());
        pf.add_error("builder failed".to_string());
        pf.add_warning("small edges ignored".to_string());
        assert!(pf.has_errors());
        assert!(pf.has_warnings());
        assert_eq!(pf.errors().len(), 2);
        assert_eq!(pf.warnings().len(), 1);
        assert_eq!(pf.errors(), ["intersection failed", "builder failed"]);
        assert_eq!(pf.warnings(), ["small edges ignored"]);
    }

    #[test]
    fn too_few_arguments_fails() {
        let mut pf = PaveFiller::new();
        let err = pf.perform().unwrap_err();
        assert!(err.contains("too few arguments"), "err: {err}");
        assert!(pf.has_errors());
    }

    #[test]
    fn glue_roundtrip() {
        let mut pf = PaveFiller::new();
        assert_eq!(pf.glue(), GlueEnum::None);
        pf.set_glue(GlueEnum::Shift);
        assert_eq!(pf.glue(), GlueEnum::Shift);
        pf.set_glue(GlueEnum::Full);
        assert_eq!(pf.glue(), GlueEnum::Full);
        pf.set_glue(GlueEnum::None);
        assert_eq!(pf.glue(), GlueEnum::None);
    }

    #[test]
    fn fuzzy_value_roundtrip_and_clamp() {
        let mut pf = PaveFiller::new();
        pf.set_fuzzy_value(0.5);
        assert_eq!(pf.fuzzy_value(), 0.5);
        pf.set_fuzzy_value(1e-9);
        assert_eq!(pf.fuzzy_value(), 1e-7, "clamped below to Precision::Confusion");
        pf.set_fuzzy_value(f64::NAN);
        assert_eq!(pf.fuzzy_value(), 1e-7, "NaN degrades to the default");
    }

    #[test]
    fn clear_resets_report_and_ds_keeps_options() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        pf.set_arguments(&[a.solid.0.clone()]);
        pf.set_glue(GlueEnum::Full);
        pf.set_fuzzy_value(0.25);
        pf.perform().unwrap();
        assert!(pf.intersection_done());
        pf.clear();
        assert!(!pf.intersection_done());
        assert!(!pf.has_errors());
        assert!(!pf.has_warnings());
        assert_eq!(pf.ds().nb_shapes(), 0);
        // User options and arguments survive Clear(), like BOPAlgo_PaveFiller::Clear().
        assert_eq!(pf.glue(), GlueEnum::Full);
        assert_eq!(pf.fuzzy_value(), 0.25);
        assert_eq!(pf.arguments().len(), 1);
    }

    #[test]
    fn repeat_intersection_flag_roundtrip() {
        let mut pf = PaveFiller::new();
        assert!(!pf.repeat_intersection());
        pf.set_repeat_intersection(true);
        assert!(pf.repeat_intersection());
    }
}
