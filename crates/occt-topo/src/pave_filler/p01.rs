use super::prelude::*;
use super::*;

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
// EdgeRangeDistance
// ---------------------------------------------------------------------------

/// EF pair with no common part, recorded for later `ProcessExistingPaveBlocks`.
/// Source: `BOPAlgo_PaveFiller::EdgeRangeDistance`.
#[derive(Debug, Clone, Copy)]
pub struct EdgeRangeDistance {
    /// First parameter of the pave block on the edge.
    pub first: f64,
    /// Last parameter of the pave block on the edge.
    pub last: f64,
    /// Minimal edge-face distance over that range.
    pub distance: f64,
}

impl EdgeRangeDistance {
    /// Constructor (`EdgeRangeDistance(first, last, distance)`).
    pub fn new(first: f64, last: f64, distance: f64) -> Self {
        Self { first, last, distance }
    }
}

impl Default for EdgeRangeDistance {
    fn default() -> Self {
        Self {
            first: 0.0,
            last: 0.0,
            distance: f64::MAX,
        }
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
    pub(super) arguments: Vec<TopoShape>,
    /// The data structure of the algorithm.
    pub(super) ds: BopdsDS,
    /// Cached intersection context (geometry tools + 2D classifiers).
    pub(super) context: IntToolsContext,
    /// Additional tolerance for touching/coinciding detection, floored at
    /// `Precision::Confusion()` (1e-7).
    pub(super) fuzzy_value: f64,
    /// Fatal alerts — a non-empty list means the algorithm has failed.
    pub(super) errors: Vec<String>,
    /// Non-fatal alerts.
    pub(super) warnings: Vec<String>,
    /// Gluing option of the algorithm.
    pub(super) glue: GlueEnum,
    /// Non-destructive mode: the argument shapes are not modified.
    pub(super) non_destructive: bool,
    /// Flag that the intersection must be repeated with increased vertex
    /// tolerances (set by the pipeline when new vertices were created).
    pub(super) repeat_intersection: bool,
    /// Flag that the intersection phase completed successfully.
    pub(super) intersection_done: bool,
    /// Primary filler of the Boolean operation (`BOPAlgo_PaveFiller::myIsPrimary`).
    /// Nested `PostTreatFF` fillers set this to false so `ForceInterfEF` is skipped.
    pub(super) is_primary: bool,
    /// `BOPAlgo_PaveFiller::myDistances` — EF pairs with no common part, keyed
    /// by `(original_edge, face)`.
    pub(super) distances: HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    /// `BOPAlgo_PaveFiller::myVertsToAvoidExtension`.
    pub(super) verts_to_avoid_extension: HashSet<usize>,
    /// `BOPAlgo_PaveFiller::myFPBDone` — pave blocks already intersected with
    /// a face, keyed by face index.
    pub(super) fpb_done: HashMap<usize, HashSet<crate::pave_ff_exist::PbKey>>,
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
            is_primary: true,
            distances: HashMap::new(),
            verts_to_avoid_extension: HashSet::new(),
            fpb_done: HashMap::new(),
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

    /// Sets whether this filler is the primary intersection (`SetIsPrimary`).
    pub fn set_is_primary(&mut self, flag: bool) {
        self.is_primary = flag;
    }

    /// True when this filler is the primary intersection (`IsPrimary`).
    pub fn is_primary(&self) -> bool {
        self.is_primary
    }

    /// `BOPAlgo_PaveFiller::myDistances`.
    pub fn distances(&self) -> &HashMap<(usize, usize), Vec<EdgeRangeDistance>> {
        &self.distances
    }

    /// Mutable `myDistances`.
    pub fn distances_mut(&mut self) -> &mut HashMap<(usize, usize), Vec<EdgeRangeDistance>> {
        &mut self.distances
    }

    /// `BOPAlgo_PaveFiller::myVertsToAvoidExtension`.
    pub fn verts_to_avoid_extension(&self) -> &HashSet<usize> {
        &self.verts_to_avoid_extension
    }

    /// Mutable `myVertsToAvoidExtension`.
    pub fn verts_to_avoid_extension_mut(&mut self) -> &mut HashSet<usize> {
        &mut self.verts_to_avoid_extension
    }

    /// `BOPAlgo_PaveFiller::myFPBDone`.
    pub fn fpb_done(&self) -> &HashMap<usize, HashSet<crate::pave_ff_exist::PbKey>> {
        &self.fpb_done
    }

    /// Mutable `myFPBDone`.
    pub fn fpb_done_mut(&mut self) -> &mut HashMap<usize, HashSet<crate::pave_ff_exist::PbKey>> {
        &mut self.fpb_done
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
        self.distances.clear();
        self.verts_to_avoid_extension.clear();
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
    /// `perform_ef` → SD-vertex update → interference-SD update →
    /// `repeat_intersection` → `force_interf_ee` → `force_interf_ef` →
    /// `refine_face_info_in` → `perform_ff` → `update_blocks_with_shared_vertices`
    /// → `make_split_edges` → SD-vertex update → `make_blocks` (EE organizer) →
    /// `make_blocks_ff` (FF section edges) →
    /// `check_self_interference` → interference-SD update →
    /// `release_pave_blocks` → `refine_face_info_on` → `remove_micro_edges` →
    /// `make_pcurves` → `process_de`.
    ///
    /// The per-stage block *splitting* happens inside each intersection stage
    /// (the OCCT `SplitPaveBlocks`); this step only *redirects* the bound
    /// vertex indices of the blocks to their same-domain representatives
    /// (`BOPDS_DS::UpdatePaveBlocksWithSDVertices`). The interference-SD update
    /// (`BOPDS_DS::UpdateInterfsWithSDVertices`) re-points the new-vertex index
    /// of every typed interference at its SD representative — OCCT calls it
    /// after the E/F stage and again after `MakeBlocks`. The post-`MakeBlocks`
    /// tail (`CheckSelfInterference`, `ReleasePaveBlocks`, `RefineFaceInfoOn`,
    /// `RemoveMicroEdges`, `ProcessDE`) mirrors `BOPAlgo_PaveFiller::PerformInternal`
    /// in order.
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
        crate::pave_common::update_interfs_with_sd_vertices(self)?;
        // OCCT `PerformInternal`: after the interference-SD update the
        // intersection is repeated for the vertices whose tolerance was
        // increased, then the edge/edge and edge/face coincidences forced.
        self.repeat_intersection_stage()?;
        self.check_errors()?;
        self.force_interf_ee()?;
        self.check_errors()?;
        self.force_interf_ef()?;
        self.check_errors()?;
        self.perform_ff()?;
        self.check_errors()?;
        // OCCT `PerformInternal`: `UpdateBlocksWithSharedVertices` runs right
        // after FF; in the default destructive mode its gate returns at once.
        self.update_blocks_with_shared_vertices();
        // OCCT `PerformInternal`: after FF the IN face-info pave blocks that
        // are also ON (boundary) blocks are dropped from the face info.
        self.ds.refine_face_info_in();
        // OCCT order: MakeSplitEdges (right after FF) precedes MakeBlocks;
        // MakePCurves follows MakeBlocks.
        self.make_split_edges()?;
        self.check_errors()?;
        self.ds.update_pave_blocks_with_sd_vertices();
        self.make_blocks()?;
        self.check_errors()?;
        // OCCT `PerformInternal`: `MakeBlocks` (`_6.cxx`) builds section edges
        // from the InterfFF curves written by PerformFF.
        crate::pave_ff::make_blocks_ff(self)?;
        self.check_errors()?;
        // OCCT `PerformInternal` (BOPAlgo_PaveFiller.cxx), after `MakeBlocks`:
        // `CheckSelfInterference` (:336) → `UpdateInterfsWithSDVertices` (:338)
        // → `ReleasePaveBlocks` (:339) → `RefineFaceInfoOn` (:340) →
        // `RemoveMicroEdges` (:342) → `MakePCurves` (:344) → `ProcessDE` (:350).
        crate::pave_common::check_self_interference(self)?;
        crate::pave_common::update_interfs_with_sd_vertices(self)?;
        self.ds_mut().release_pave_blocks();
        self.ds_mut().refine_face_info_on();
        crate::pave_common::remove_micro_edges(self);
        self.make_pcurves()?;
        self.check_errors()?;
        crate::pave_de::process_de(self)?;
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
    pub(super) fn prepare(&mut self) -> Result<(), String> {
        let n = self.ds().nb_source_shapes();
        for i in 0..n {
            if self.ds().shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Edge) {
                self.ds_mut().init_pave_blocks_for_edge(i);
            }
        }
        crate::pave_common::fill_shrunk_data(self)?;
        crate::pave_pcurves::prepare_pcurves_on_planes(self)
    }

    /// Vertex/Vertex intersection. Source: `PerformVV`
    /// (`BOPAlgo_PaveFiller_1.cxx`). Delegates to
    /// [`crate::pave_intersect::perform_vv`].
    pub(super) fn perform_vv(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_vv(self)
    }

    /// Vertex/Edge intersection. Source: `PerformVE`. Delegates to
    /// [`crate::pave_intersect::perform_ve`].
    pub(super) fn perform_ve(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ve(self)
    }

    /// Edge/Edge intersection. Source: `PerformEE`. Delegates to
    /// [`crate::pave_intersect::perform_ee`].
    pub(super) fn perform_ee(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ee(self)
    }

    /// Vertex/Face intersection. Source: `PerformVF`. Delegates to
    /// [`crate::pave_intersect::perform_vf`].
    pub(super) fn perform_vf(&mut self) -> Result<(), String> {
        crate::pave_vf::perform_vf(self)
    }

    /// Edge/Face intersection. Source: `PerformEF`. Delegates to
    /// [`crate::pave_intersect::perform_ef`].
    pub(super) fn perform_ef(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ef(self)
    }

    /// Repeats the intersection for the vertices whose tolerance was increased
    /// during the previous stages.
    ///
    /// Source: `BOPAlgo_PaveFiller::RepeatIntersection`
    /// (`BOPAlgo_PaveFiller.cxx`). The vertices whose tolerance grew (and the
    /// source vertices linked to them through the SD map) get an extended
    /// interference pair set (`BOPDS_Iterator::IntersectExt`), then V/V, V/E
    /// and V/F are re-run on those pairs only.
    ///
    /// Named `_stage` to avoid colliding with the
    /// [`repeat_intersection`](Self::repeat_intersection) option accessor.
    pub(super) fn repeat_intersection_stage(&mut self) -> Result<(), String> {
        let extra_map = {
            let ds = self.ds();
            let increased = ds.increased_ss();
            if increased.is_empty() {
                return Ok(());
            }
            let mut extra_map: HashSet<usize> = HashSet::new();
            let n = ds.nb_source_shapes();
            for i in 0..n {
                let is_vertex = ds
                    .shape_info(i)
                    .map(|s| s.shape_type() == ShapeType::Vertex)
                    .unwrap_or(false);
                if !is_vertex {
                    continue;
                }
                // The original vertex had its tolerance increased directly...
                if increased.contains(&i) {
                    extra_map.insert(i);
                    continue;
                }
                // ...or it was linked to a same-domain vertex whose tolerance
                // grew.
                if let Some(n_vsd) = ds.has_shape_sd(i) {
                    if increased.contains(&n_vsd) {
                        extra_map.insert(i);
                    }
                }
            }
            extra_map
        };
        if extra_map.is_empty() {
            return Ok(());
        }
        let buckets = crate::bopds::intersect_ext_pairs(self.ds(), &extra_map);
        // Re-run the vertex stages on the extended pairs.
        crate::pave_intersect::perform_vv_pairs(self, &buckets[0])?;
        self.ds_mut().update_pave_blocks_with_sd_vertices();
        crate::pave_intersect::perform_ve_pairs(self, &buckets[1])?;
        self.ds_mut().update_pave_blocks_with_sd_vertices();
        crate::pave_intersect::perform_vf_pairs(self, &buckets[3])?;
        self.ds_mut().update_pave_blocks_with_sd_vertices();
        crate::pave_common::update_interfs_with_sd_vertices(self)?;
        Ok(())
    }

    /// Force intersection of the edges after the increase of the tolerance
    /// values of their vertices.
    ///
    /// Source: `BOPAlgo_PaveFiller::ForceInterfEE`
    /// (`BOPAlgo_PaveFiller_3.cxx`). Looks for additional edge/edge common
    /// blocks among the pairs of pave blocks bounded by the same vertices.
    pub(super) fn force_interf_ee(&mut self) -> Result<(), String> {
        crate::pave_force_ee::force_interf_ee(self)
    }

    /// Force edge/face intersection after the increase of the tolerance values
    /// of their vertices.
    ///
    /// Source: `BOPAlgo_PaveFiller::ForceInterfEF`
    /// (`BOPAlgo_PaveFiller_5.cxx`). Looks for additional edge/face common
    /// blocks among the pairs of pave blocks whose bounding vertices lie on
    /// the face.
    pub(super) fn force_interf_ef(&mut self) -> Result<(), String> {
        crate::pave_force_ef::force_interf_ef(self)
    }

    /// Face/Face intersection. Source: `PerformFF`. Delegates to
    /// [`crate::pave_intersect::perform_ff`].
    pub(super) fn perform_ff(&mut self) -> Result<(), String> {
        crate::pave_intersect::perform_ff(self)
    }

    /// Updates the pave blocks of the faces with vertices shared by the faces.
    ///
    /// Source: `BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices`
    /// (`BOPAlgo_PaveFiller_6.cxx`). The whole body is gated behind the
    /// non-destructive mode (`if (!myNonDestructive) return;`); the default
    /// destructive mode returns immediately, so only the gate is ported — the
    /// non-destructive-only body (`EstimatePaveOnCurve` + the shared-vertex
    /// updates) is left as a translation boundary.
    pub(super) fn update_blocks_with_shared_vertices(&mut self) {
        if let Err(e) = crate::pave_ff_misc::update_blocks_with_shared_vertices(self) {
            self.add_error(e);
        }
    }

    /// Groups coincident pave blocks into common blocks. The EE organizer
    /// lives in [`crate::pave_blocks::make_blocks`]; FF section edges are built
    /// afterwards by [`crate::pave_ff::make_blocks_ff`] (`BOPAlgo_PaveFiller_6.cxx`).
    pub(super) fn make_blocks(&mut self) -> Result<(), String> {
        crate::pave_blocks::make_blocks(self)
    }

    /// Builds 2D curves of the section edges on the faces. Source: `MakePCurves`.
    /// Delegates to [`crate::pave_blocks::make_pcurves`].
    pub(super) fn make_pcurves(&mut self) -> Result<(), String> {
        crate::pave_pcurves::make_p_curves(self)
    }

    /// Creates the split edges (new edges of the common part). Source:
    /// `MakeSplitEdges`. Delegates to [`crate::pave_blocks::make_split_edges`].
    pub(super) fn make_split_edges(&mut self) -> Result<(), String> {
        crate::pave_split::make_split_edges(self)
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Short-circuits the pipeline when a fatal alert was recorded by a stage.
    pub(super) fn check_errors(&self) -> Result<(), String> {
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
