//! General Fuse builder + Boolean wrapper — the building phase of the exact
//! NURBS boolean (Phase 20).
//!
//! Port of `BOPAlgo_Builder` (`BOPAlgo_Builder.cxx` + `BOPAlgo_Builder_1.cxx`)
//! and `BOPAlgo_BOP` (`BOPAlgo_BOP.hxx/.cxx`) from TKBO/BOPAlgo.
//!
//! [`BopBuilder`] is the *General Fuse* algorithm — the base algorithm of the
//! Boolean Component. It consumes a [`PaveFiller`] whose intersection phase
//! already filled the [`crate::bopds::BopdsDS`] with every participating
//! shape, its pave blocks and the section vertices/edges, then:
//!
//! 1. [`BopBuilder::perform`] runs the intersection phase (`PaveFiller::perform`)
//!    and the building phase;
//! 2. [`BopBuilder::build_bop`] fills the *images* of the sub-shapes — each
//!    source shape maps to the split pieces it expands into after the
//!    intersection ([`crate::bop_hist::BopHistory`] images table + the origins
//!    back-map + the same-domain shapes map);
//! 3. [`BopBuilder::build_result`] assembles the result compound from the
//!    images of the arguments;
//! 4. [`BopBuilder::post_treat`] fixes the vertex tolerances of the result.
//!
//! The heavy face/solid reconstruction is delegated to the sibling Phase-20
//! modules `crate::bop_build_faces` (split faces) and `crate::bop_build_common`
//! (containers, internal shapes, draft solids). Those are filled by parallel
//! agents; this module holds the orchestration and the vertex/edge image
//! filling, with the face/solid stages present as stubs for now.
//!
//! [`BoolOp2`] + [`builder_bop`] / [`fuse`] / [`cut`] / [`common`] provide the
//! `BOPAlgo_BOP` wrapper: the operation type is converted into the states the
//! object/tool faces must have relative to the opposite group to pass into the
//! result, and the states drive the face selection during the building phase.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::bop_build_common::{is_split_to_reverse, BopBuildOps};
use crate::bop_build_faces::BopBuilderLike;
use crate::bop_hist::{is_supported_type, BopHistory};
use crate::bopds::{BopdsDS, BopdsInterf, BopdsPaveBlock};
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::pave_filler::PaveFiller;
use crate::shape::TopoShape;
use crate::topo_tools_full::all_subshapes;

/// Result-shape assembly order — matches the per-type `BuildResult` calls of
/// OCCT `BOPAlgo_Builder::PerformInternal1` (lowest type first).
const RESULT_TYPES: [ShapeType; 8] = [
    ShapeType::Vertex,
    ShapeType::Edge,
    ShapeType::Wire,
    ShapeType::Face,
    ShapeType::Shell,
    ShapeType::Solid,
    ShapeType::CompSolid,
    ShapeType::Compound,
];

/// The General Fuse algorithm — base algorithm of the Boolean Component.
///
/// Source: `BOPAlgo_Builder`.
///
/// The class owns the [`PaveFiller`] (intersection phase), the
/// [`BopHistory`] images table (each old sub-shape → its split pieces), the
/// origins back-map, the same-domain shapes map and the alert report of the
/// run. The public surface mirrors the OCCT class reduced to what the Rust
/// pipeline needs; the sibling Phase-20 modules (`crate::bop_build_faces`,
/// `crate::bop_build_common`) read the images/history and the object/tool
/// state classification from this builder.
#[derive(Debug, Clone)]
pub struct BopBuilder {
    /// Pave Filler — the intersection phase of the algorithm.
    filler: PaveFiller,
    /// Images naming table + history: each old shape → the pieces it was
    /// split/expanded into during the intersection (`myImages`).
    history: BopHistory,
    /// Origins — the back map of the images table: split shape → the source
    /// shapes that expand into it, keyed by the split shape's `TShape` address
    /// (`myOrigins`).
    origins: HashMap<usize, Vec<TopoShape>>,
    /// Same-domain shapes: source index → the coincident shape it maps to
    /// (`myShapesSD`).
    shapes_sd: HashMap<usize, usize>,
    /// Arguments of the operation (deduplicated by TShape identity).
    arguments: Vec<TopoShape>,
    /// Objects group of the current BOP run.
    objects: Vec<TopoShape>,
    /// Tools group of the current BOP run.
    tools: Vec<TopoShape>,
    /// State the object faces must have relative to the tools to pass into the
    /// result.
    obj_state: FaceState,
    /// State the tool faces must have relative to the objects to pass into the
    /// result.
    tools_state: FaceState,
    /// Fatal alerts — a non-empty list means the algorithm has failed.
    errors: Vec<String>,
    /// Non-fatal alerts.
    warnings: Vec<String>,
    /// The shape produced by the last run.
    result_shape: TopoShape,
}

impl Default for BopBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl BopBuilder {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            filler: PaveFiller::new(),
            history: BopHistory::new(),
            origins: HashMap::new(),
            shapes_sd: HashMap::new(),
            arguments: Vec::new(),
            objects: Vec::new(),
            tools: Vec::new(),
            obj_state: FaceState::Out,
            tools_state: FaceState::Out,
            errors: Vec::new(),
            warnings: Vec::new(),
            result_shape: TopoShape::new(ShapeType::Compound),
        }
    }

    // -----------------------------------------------------------------------
    // Arguments
    // -----------------------------------------------------------------------

    /// Sets the list of arguments for the operation, deduplicated by `TShape`
    /// identity (mirrors `BOPAlgo_Builder::SetArguments` + `AddArgument` fence).
    pub fn set_arguments(&mut self, shapes: &[TopoShape]) {
        self.arguments.clear();
        for s in shapes {
            let flat = flatten_location(s);
            if !self.arguments.iter().any(|a| a.same_tshape(&flat)) {
                self.arguments.push(flat);
            }
        }
    }

    /// Adds an argument to the operation, skipping duplicates.
    pub fn add_argument(&mut self, s: &TopoShape) {
        let flat = flatten_location(s);
        if !self.arguments.iter().any(|a| a.same_tshape(&flat)) {
            self.arguments.push(flat);
        }
    }

    /// Returns the arguments of the operation.
    pub fn arguments(&self) -> &[TopoShape] {
        &self.arguments
    }

    // -----------------------------------------------------------------------
    // Alert report
    // -----------------------------------------------------------------------

    /// True if the algorithm has failed (at least one fatal alert).
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// Returns the collected fatal alerts.
    pub fn errors(&self) -> &[String] {
        &self.errors
    }

    /// Returns the collected non-fatal alerts.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Adds a fatal alert.
    pub fn add_error(&mut self, msg: String) {
        self.errors.push(msg);
    }

    /// Adds a non-fatal alert.
    pub fn add_warning(&mut self, msg: String) {
        self.warnings.push(msg);
    }

    // -----------------------------------------------------------------------
    // Data-structure / history access
    // -----------------------------------------------------------------------

    /// Returns the data structure of the algorithm (through the filler).
    pub fn ds(&self) -> &BopdsDS {
        self.filler.ds()
    }

    /// Returns the mutable data structure of the algorithm.
    pub fn ds_mut(&mut self) -> &mut BopdsDS {
        self.filler.ds_mut()
    }

    /// Returns the images/history naming table.
    pub fn history(&self) -> &BopHistory {
        &self.history
    }

    /// Returns the mutable images/history naming table.
    pub fn history_mut(&mut self) -> &mut BopHistory {
        &mut self.history
    }

    /// Returns the origins back-map (`split shape key → source shapes`).
    pub fn origins(&self) -> &HashMap<usize, Vec<TopoShape>> {
        &self.origins
    }

    /// Returns the mutable origins back-map.
    pub fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
        &mut self.origins
    }

    /// Returns the same-domain shapes map of the builder.
    pub fn shapes_sd(&self) -> &HashMap<usize, usize> {
        &self.shapes_sd
    }

    /// The additional tolerance of the operation (from the PaveFiller).
    pub fn fuzzy_value(&self) -> f64 {
        self.filler.fuzzy_value()
    }

    /// True when the PaveFiller runs in non-destructive mode.
    pub fn non_destructive(&self) -> bool {
        self.filler.non_destructive()
    }

    /// Binds `shape` as same-domain with `sd`, mirroring `myShapesSD.Bind`.
    ///
    /// The shapes are resolved to their data-structure indices and the binding
    /// is recorded (identity bindings are skipped).
    pub fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
        if let (Some(i), Some(j)) = (self.filler.ds().index(&shape), self.filler.ds().index(&sd)) {
            if i != j {
                self.shapes_sd.insert(i, j);
            }
        }
    }

    /// Returns the same-domain representative of `shape`, when bound
    /// (`myShapesSD.Seek`), following the SD chain.
    pub fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
        let start = self.filler.ds().index(shape)?;
        let mut cur = start;
        let mut guard = 0;
        while let Some(&j) = self.shapes_sd.get(&cur) {
            cur = j;
            guard += 1;
            if guard > 64 {
                break;
            }
        }
        if cur == start {
            None
        } else {
            self.filler.ds().shape(cur).cloned()
        }
    }

    /// Returns the objects group of the current run.
    pub fn objects(&self) -> &[TopoShape] {
        &self.objects
    }

    /// Returns the tools group of the current run.
    pub fn tools(&self) -> &[TopoShape] {
        &self.tools
    }

    /// Returns the state the object faces must have to pass into the result.
    pub fn obj_state(&self) -> FaceState {
        self.obj_state
    }

    /// Returns the state the tool faces must have to pass into the result.
    pub fn tools_state(&self) -> FaceState {
        self.tools_state
    }

    /// Returns the shape produced by the last run.
    pub fn result(&self) -> &TopoShape {
        &self.result_shape
    }

    /// Returns the PaveFiller of the algorithm.
    pub fn filler(&self) -> &PaveFiller {
        &self.filler
    }

    /// Returns the mutable PaveFiller of the algorithm.
    pub fn filler_mut(&mut self) -> &mut PaveFiller {
        &mut self.filler
    }

    // -----------------------------------------------------------------------
    // Performing the operation
    // -----------------------------------------------------------------------

    /// Performs the General Fuse operation: runs the intersection phase
    /// (`PaveFiller::perform`) followed by the building phase
    /// (`build_bop` → `build_result` → `post_treat`).
    ///
    /// The result of the General Fuse is a compound containing all split parts
    /// of the arguments.
    pub fn perform(&mut self) -> Result<TopoShape, String> {
        let args = self.arguments.clone();
        self.perform_internal(&args, FaceState::Out, &[], FaceState::Out)
    }

    /// The shared pipeline of [`BopBuilder::perform`] and the [`BOPAlgo_BOP`]
    /// wrapper: intersect, fill the images, assemble the result and post-treat.
    ///
    /// `objects`/`tools` are the two groups of the BOP; `obj_state`/
    /// `tools_state` are the states the faces of each group must have relative
    /// to the opposite group to pass into the result. The General Fuse treats
    /// all arguments as objects with state `Out` and no tools.
    fn perform_internal(
        &mut self,
        objects: &[TopoShape],
        obj_state: FaceState,
        tools: &[TopoShape],
        tools_state: FaceState,
    ) -> Result<TopoShape, String> {
        self.reset_run();
        if self.arguments.len() < 2 {
            let msg = "BOPAlgo_Builder: too few arguments".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }

        // Intersection phase.
        self.filler.set_arguments(&self.arguments);
        if let Err(e) = self.filler.perform() {
            let msgs: Vec<String> = self.filler.errors().to_vec();
            for m in msgs {
                self.add_error(m);
            }
            if self.errors.is_empty() {
                self.add_error(e.clone());
            }
            return Err(self.errors.first().cloned().unwrap_or_else(|| e));
        }
        if self.filler.has_errors() {
            let msgs: Vec<String> = self.filler.errors().to_vec();
            for m in msgs {
                self.add_error(m);
            }
            return Err(self.errors.first().cloned().unwrap_or_default());
        }

        // Building phase.
        self.build_bop(objects, obj_state, tools, tools_state)?;
        self.build_result()?;
        self.prepare_history();
        self.post_treat()?;
        Ok(self.result_shape.clone())
    }

    /// Resets the per-run state (report, history, origins, same-domain map,
    /// result). The arguments and the filler options survive.
    fn reset_run(&mut self) {
        self.errors.clear();
        self.warnings.clear();
        self.history.clear();
        self.origins.clear();
        self.shapes_sd.clear();
        self.objects.clear();
        self.tools.clear();
        self.obj_state = FaceState::Out;
        self.tools_state = FaceState::Out;
        self.result_shape = TopoShape::new(ShapeType::Compound);
    }

    /// Clears the content of the algorithm: the report, the history, the
    /// arguments and the filler (`BOPAlgo_Builder::Clear`).
    pub fn clear(&mut self) {
        self.errors.clear();
        self.warnings.clear();
        self.history.clear();
        self.origins.clear();
        self.shapes_sd.clear();
        self.arguments.clear();
        self.objects.clear();
        self.tools.clear();
        self.result_shape = TopoShape::new(ShapeType::Compound);
    }

    // -----------------------------------------------------------------------
    // Building the result
    // -----------------------------------------------------------------------

    /// Builds the result of the Boolean operation of the given type on the
    /// given groups.
    ///
    /// Port of `BOPAlgo_Builder::BuildBOP`, reduced to the orchestration: the
    /// method validates the inputs, records the object/tool groups and their
    /// states, then fills the images of the sub-shapes. The face and solid
    /// stages are delegated to the sibling Phase-20 modules.
    ///
    /// The state filtering itself (which faces pass into the result) is
    /// performed by the face/solid building stages, which read
    /// [`BopBuilder::obj_state`] / [`BopBuilder::tools_state`].
    pub fn build_bop(
        &mut self,
        objects: &[TopoShape],
        obj_state: FaceState,
        tools: &[TopoShape],
        tools_state: FaceState,
    ) -> Result<(), String> {
        if self.has_errors() {
            return Err(self.errors.first().cloned().unwrap_or_default());
        }
        if !matches!(obj_state, FaceState::In | FaceState::Out)
            || !matches!(tools_state, FaceState::In | FaceState::Out)
        {
            let msg = "BOPAlgo_Builder: invalid state for the operation".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        if objects.is_empty() && tools.is_empty() {
            let msg = "BOPAlgo_Builder: too few arguments".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        for s in objects.iter().chain(tools.iter()) {
            if self.filler.ds().index(s).is_none() {
                let msg = "BOPAlgo_Builder: unknown shape for the operation".to_string();
                self.add_error(msg.clone());
                return Err(msg);
            }
        }
        self.objects = objects.to_vec();
        self.tools = tools.to_vec();
        self.obj_state = obj_state;
        self.tools_state = tools_state;

        self.fill_images_vertices()?;
        self.fill_images_edges()?;
        crate::bop_build_common::fill_images_containers(self, ShapeType::Wire)?;
        self.fill_images_faces()?;
        crate::bop_build_common::fill_images_containers(self, ShapeType::Shell)?;
        self.fill_images_solids()?;
        crate::bop_build_common::fill_images_containers(self, ShapeType::CompSolid)?;
        Ok(())
    }

    /// Fills the images of the vertices.
    ///
    /// Port of `BOPAlgo_Builder::FillImagesVertices`: for every same-domain
    /// vertex pair `(v, v_sd)` of the data structure, records `v_sd` as the
    /// image of `v` in the images table, binds the same-domain map and adds
    /// the origin back-map entry `v_sd → [v]`.
    fn fill_images_vertices(&mut self) -> Result<(), String> {
        let sd = self.filler.ds().shapes_sd().clone();
        for (&n_v, &n_vsd) in &sd {
            let Some(v) = self.filler.ds().shape(n_v).cloned() else { continue };
            let Some(vsd) = self.filler.ds().shape(n_vsd).cloned() else { continue };
            self.history.add_image(&v, vsd.clone());
            self.shapes_sd.insert(n_v, n_vsd);
            self.origins.entry(shape_key(&vsd)).or_default().push(v.clone());
        }
        Ok(())
    }

    /// Fills the images of the edges.
    ///
    /// Port of `BOPAlgo_Builder::FillImagesEdges`: every source edge with pave
    /// blocks expands into the split edges of its blocks. The split edges were
    /// built by the intersection phase ([`crate::pave_blocks::make_split_edges`]
    /// → [`crate::algo_tools::AlgoTools::make_split_edge`]) and stored in the
    /// data structure; each block's `edge()` points at its split edge. The
    /// split edge is recorded as the image of the source edge, and the origin
    /// back-map gets the source edge. A block on a common block would record
    /// the same-domain binding of its edge (an identity in this port, since
    /// the block's own edge is already the real split edge).
    fn fill_images_edges(&mut self) -> Result<(), String> {
        let n = self.filler.ds().nb_source_shapes();
        for i in 0..n {
            let Some(si) = self.filler.ds().shape_info(i).cloned() else { continue };
            if si.shape_type() != ShapeType::Edge {
                continue;
            }
            if !self.filler.ds().has_pave_blocks(i) {
                continue;
            }
            let e = si.shape().clone();
            let blocks = self.filler.ds().pave_blocks(i).to_vec();
            for pb in &blocks {
                // A zero-length block (two coincident bound vertices) is not a
                // real interval and must not contribute an image — otherwise the
                // split image list gains a spurious full-edge entry.
                if (pb.last - pb.first).abs() <= occt_core::precision::PCONFUSION {
                    continue;
                }
                let n_sp_r = pb.edge();
                let Some(sp_r) = self.filler.ds().shape(n_sp_r).cloned() else { continue };
                self.history.add_image(&e, sp_r.clone());
                self.origins.entry(shape_key(&sp_r)).or_default().push(e.clone());
                if is_common_block_on_edge(self.filler.ds(), pb) {
                    // The block lies on a common block — a coincident edge shared
                    // by several source edges. The block's own edge is already the
                    // real split edge in this port, so the same-domain binding is
                    // an identity and is not recorded.
                    let n_sp = pb.edge();
                    if n_sp != n_sp_r {
                        self.shapes_sd.insert(n_sp, n_sp_r);
                    }
                }
            }
        }
        Ok(())
    }

    /// Fills the images of the faces — split faces of the arguments.
    ///
    /// Delegates to `crate::bop_build_faces::fill_images_faces`
    /// (`BuildSplitFaces` → `FillSameDomainFaces`), which records the split
    /// faces into the images table with the [`BopBuilder::obj_state`] /
    /// [`BopBuilder::tools_state`] classification.
    fn fill_images_faces(&mut self) -> Result<(), String> {
        crate::bop_build_faces::fill_images_faces(self)
    }

    /// Fills the images of the solids — split solids of the arguments.
    ///
    /// Delegates to `crate::bop_build_solids::build_split_solids_full`, which
    /// splits each interfered solid into its closed-shell pieces and selects the
    /// pieces that belong to the result per the operation states
    /// ([`BopBuilder::obj_state`] / [`BopBuilder::tools_state`]).
    fn fill_images_solids(&mut self) -> Result<(), String> {
        let objects = self.objects.clone();
        let tools = self.tools.clone();
        let obj_state = self.obj_state;
        let tools_state = self.tools_state;
        crate::bop_build_solids::build_split_solids_full(
            self,
            &objects,
            &tools,
            obj_state,
            tools_state,
        )
    }

    /// Assembles the result compound from the images of the arguments.
    ///
    /// Port of `BOPAlgo_Builder::BuildResult(theType)` collapsed into a single
    /// pass over the result types in the OCCT order: for every argument of a
    /// given type, its images (split pieces) are added to the result; an
    /// argument without images is added as-is. A fence map (by `TShape`
    /// identity) keeps the result free of duplicates.
    fn build_result(&mut self) -> Result<TopoShape, String> {
        let b = TopoBuilder::new();
        let mut result = b.make_compound_of(&[]);
        let mut fence: Vec<TopoShape> = Vec::new();
        for t in RESULT_TYPES {
            for arg in &self.arguments {
                if arg.shape_type() != t {
                    continue;
                }
                let is_tool = self.tools.iter().any(|x| x.same_tshape(arg));
                let imgs: Vec<TopoShape> = match self.history.image(arg) {
                    Some(list) => list.to_vec(),
                    None => {
                        // An unsplit argument is returned as-is, except the
                        // tools of a Cut/Common that contributed no piece (they
                        // are entirely removed by the operation — a Cut removes
                        // its tools, a Common keeps only the overlap, which the
                        // object side already covers).
                        if is_tool && self.tools_state == FaceState::In {
                            continue;
                        }
                        vec![arg.clone()]
                    }
                };
                for img in imgs {
                    if !fence.iter().any(|f| f.same_tshape(&img)) {
                        fence.push(img.clone());
                        b.add(&mut result.0, &img);
                    }
                }
            }
        }
        self.result_shape = result.0.clone();
        Ok(result.0)
    }

    /// Post-treats the result shape by correcting the tolerances.
    ///
    /// Port of `BOPAlgo_Builder::PostTreat` (`BOPAlgo_Builder.cxx`): the
    /// `CorrectTolerances` + `CorrectShapeTolerances` pair is covered by
    /// [`AlgoTools::correct_tolerances`]. The OCCT `aMA` avoid-map of the
    /// source V/E/F shapes is only populated in the non-destructive mode,
    /// which this port does not exercise (the tolerance pass runs against the
    /// result shape, matching the default destructive mode).
    fn post_treat(&mut self) -> Result<(), String> {
        AlgoTools::correct_tolerances(&self.result_shape, 0.05);
        Ok(())
    }

    /// Fills the modified/generated/removed relations of the history from the
    /// images table and the result shape.
    ///
    /// Port of `BOPAlgo_Builder::PrepareHistory` (`BOPAlgo_Builder_4.cxx`):
    /// for every source shape of the data structure,
    /// - the split pieces of the shape kept in the result become **modified**
    ///   from it;
    /// - the vertices/edges the intersections created from an EDGE/FACE source
    ///   become **generated** from it ([`BopBuilder::loc_generated`]);
    /// - a shape with no trace in the result (and no surviving splits) is
    ///   marked **removed**.
    fn prepare_history(&mut self) {
        // All shapes of the result (the root and every sub-shape), keyed by
        // (TShape, orientation) — the OCCT `TopTools_MapOfShape` identity.
        let result_keys: HashSet<(usize, u8)> = all_subshapes(&self.result_shape)
            .iter()
            .map(|s| (Arc::as_ptr(&s.tshape) as usize, s.orientation() as u8))
            .collect();
        let in_result = |s: &TopoShape| {
            result_keys.contains(&(Arc::as_ptr(&s.tshape) as usize, s.orientation() as u8))
        };

        let n = self.filler.ds().nb_source_shapes();
        for i in 0..n {
            let Some(si) = self.filler.ds().shape_info(i) else { continue };
            let s = si.shape().clone();
            if !is_supported_type(&s) {
                continue;
            }

            let mut is_modified = false;
            // Modified: the splits of the shape kept in the result, oriented
            // like the source (VERTEX/SOLID take the source orientation; an
            // EDGE/FACE whose direction flips is reversed).
            let splits: Vec<TopoShape> =
                self.history.image(&s).map(|v| v.to_vec()).unwrap_or_default();
            for sp in &splits {
                if !in_result(sp) {
                    continue;
                }
                let mut sp = sp.clone();
                let t = sp.shape_type();
                if t == ShapeType::Vertex || t == ShapeType::Solid {
                    sp.set_orientation(s.orientation());
                } else if is_split_to_reverse(&sp, &s) {
                    sp.set_orientation(sp.orientation().reversed());
                }
                let _ = self.history.add_modified(&s, sp);
                is_modified = true;
            }

            // Generated: the vertices/edges the intersections created from the
            // shape, kept in the result.
            for g in self.loc_generated(&s) {
                if in_result(&g) {
                    let _ = self.history.add_generated(&s, g);
                }
            }

            // Removed: the shape has no trace in the result nor any surviving
            // split.
            if !is_modified && !in_result(&s) {
                let _ = self.history.add_removed(s);
            }
        }
    }

    /// Shapes generated from `s` by the intersections: the vertices created in
    /// E/E and E/F interferences (for an EDGE or FACE source) plus, for a
    /// FACE, the section edges and vertices lying on it.
    ///
    /// Port of `BOPAlgo_Builder::LocGenerated` (`BOPAlgo_Builder_4.cxx`). The
    /// new-vertex part reads the `index_new` records of the E/E and E/F
    /// interferences; the section edges/vertices of a face come from its
    /// [`crate::bopds::BopdsFaceInfo`]. Whether a returned shape actually made
    /// it into the result is checked by the caller
    /// ([`BopBuilder::prepare_history`]).
    fn loc_generated(&self, s: &TopoShape) -> Vec<TopoShape> {
        let mut out: Vec<TopoShape> = Vec::new();
        let a_type = s.shape_type();
        if a_type != ShapeType::Edge && a_type != ShapeType::Face {
            return out;
        }
        let ds = self.filler.ds();
        let Some(n_s) = ds.index(s) else { return out };
        // Untouched shapes carry no generated elements — an edge without pave
        // blocks, a face without a face-info entry (the OCCT `HasReference`
        // guard).
        let is_face = a_type == ShapeType::Face;
        if is_face {
            if !ds.face_info_pool().iter().any(|fi| fi.face_index == n_s) {
                return out;
            }
        } else if !ds.has_pave_blocks(n_s) {
            return out;
        }

        // New vertices of the E/E (edge sources) and E/F interferences
        // containing the shape, deduplicated by their same-domain index.
        let mut fence: Vec<usize> = Vec::new();
        if !is_face {
            Self::collect_interf_vertices(ds, ds.interf_ee(), n_s, &mut fence, &mut out);
        }
        Self::collect_interf_vertices(ds, ds.interf_ef(), n_s, &mut fence, &mut out);
        if !is_face {
            return out;
        }

        // Section edges and section vertices lying on the face.
        if let Some(fi) = ds.face_info_pool().iter().find(|fi| fi.face_index == n_s) {
            for &(e_idx, _, _) in fi.paves() {
                if let Some(e) = ds.shape(e_idx) {
                    out.push(e.clone());
                }
            }
            for &(v_idx, _, _) in fi.verts() {
                if let Some(v) = ds.shape(v_idx) {
                    out.push(v.clone());
                }
            }
        }
        out
    }

    /// Appends the `index_new` vertices of the interferences `ints` containing
    /// `n_s`, deduplicated by their same-domain vertex index (the OCCT
    /// `LocGenerated` loop over `InterfEE` / `InterfEF`).
    fn collect_interf_vertices(
        ds: &BopdsDS,
        ints: &[BopdsInterf],
        n_s: usize,
        fence: &mut Vec<usize>,
        out: &mut Vec<TopoShape>,
    ) {
        for it in ints {
            let Some(n_v) = it.get_index_new() else { continue };
            if !it.contains(n_s) {
                continue;
            }
            // Resolve the vertex through its same-domain twin.
            let n_v = ds.has_shape_sd(n_v).unwrap_or(n_v);
            if fence.contains(&n_v) {
                continue;
            }
            fence.push(n_v);
            if let Some(v) = ds.shape(n_v) {
                out.push(v.clone());
            }
        }
    }
}

/// The face-reconstruction host contract (`crate::bop_build_faces`).
///
/// Every accessor delegates to the corresponding inherent [`BopBuilder`]
/// method; the trait is what lets the sibling Phase-20 modules
/// (`crate::bop_build_faces::fill_images_faces` & friends) run against this
/// builder.
impl BopBuilderLike for BopBuilder {
    fn ds(&self) -> &BopdsDS {
        self.ds()
    }
    fn ds_mut(&mut self) -> &mut BopdsDS {
        self.ds_mut()
    }
    fn history(&self) -> &BopHistory {
        self.history()
    }
    fn history_mut(&mut self) -> &mut BopHistory {
        self.history_mut()
    }
    fn has_errors(&self) -> bool {
        self.has_errors()
    }
    fn add_error(&mut self, msg: String) {
        self.add_error(msg);
    }
    fn add_warning(&mut self, msg: String) {
        self.add_warning(msg);
    }
    fn non_destructive(&self) -> bool {
        self.non_destructive()
    }
    fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value()
    }
    fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
        self.bind_shapes_sd(shape, sd);
    }
    fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
        self.seek_shapes_sd(shape)
    }
}

/// The container/internal-shape host contract (`crate::bop_build_common`).
///
/// All accessors delegate to inherent [`BopBuilder`] methods; the trait lets
/// the solid stage (`crate::bop_build_solids`) and the container stages
/// (`crate::bop_build_common`) run against this builder.
impl BopBuildOps for BopBuilder {
    fn ds(&self) -> &BopdsDS {
        self.ds()
    }
    fn history(&self) -> &BopHistory {
        self.history()
    }
    fn history_mut(&mut self) -> &mut BopHistory {
        self.history_mut()
    }
    fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value()
    }
    fn arguments(&self) -> &[TopoShape] {
        self.arguments()
    }
    fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
        self.origins_mut()
    }
}

/// The Boolean operation to apply. Source: `BOPAlgo_Operation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp2 {
    /// A ∪ B.
    Fuse,
    /// A − B.
    Cut,
    /// A ∩ B.
    Common,
}

impl BoolOp2 {
    /// Converts the operation into the states of the object/tool faces that
    /// pass into the result — the `BOPAlgo_Builder::BuildBOP` conversion:
    ///
    /// - [`BoolOp2::Fuse`] — faces of the objects OUT of the tools + faces of
    ///   the tools OUT of the objects → the union. States `(Out, Out)`.
    /// - [`BoolOp2::Cut`] — faces of the objects OUT of the tools + faces of
    ///   the tools IN the objects → object minus tools. States `(Out, In)`.
    /// - [`BoolOp2::Common`] — faces of the objects IN the tools + faces of
    ///   the tools IN the objects → the intersection. States `(In, In)`.
    pub fn states(&self) -> (FaceState, FaceState) {
        match self {
            BoolOp2::Fuse => (FaceState::Out, FaceState::Out),
            BoolOp2::Cut => (FaceState::Out, FaceState::In),
            BoolOp2::Common => (FaceState::In, FaceState::In),
        }
    }
}

/// Performs the Boolean operation of type `op` on the `objects` and `tools`
/// groups.
///
/// Port of `BOPAlgo_BOP::Perform`: the arguments of both groups are
/// intersected together, then the operation type is converted into the states
/// for the face selection and the result is built.
pub fn builder_bop(
    objects: &[TopoShape],
    tools: &[TopoShape],
    op: BoolOp2,
) -> Result<TopoShape, String> {
    // Bake any non-identity location into the geometry so the intersection
    // phase runs in world coordinates (the DS/intersection kernels read the
    // registered geometry without applying the `TopoShape` location).
    let objects: Vec<TopoShape> = objects.iter().map(flatten_location).collect();
    let tools: Vec<TopoShape> = tools.iter().map(flatten_location).collect();
    let mut all: Vec<TopoShape> = objects.clone();
    for t in &tools {
        if !all.iter().any(|a| a.same_tshape(t)) {
            all.push(t.clone());
        }
    }
    let mut b = BopBuilder::new();
    b.set_arguments(&all);
    let (os, ts) = op.states();
    b.perform_internal(&objects, os, &tools, ts)
}

/// Union of two groups of shapes. Source: `BOPAlgo_BOP` with `BOPAlgo_FUSE`.
pub fn fuse(objects: &[TopoShape], tools: &[TopoShape]) -> Result<TopoShape, String> {
    builder_bop(objects, tools, BoolOp2::Fuse)
}

/// Subtraction of the tools from the objects. Source: `BOPAlgo_CUT`.
pub fn cut(objects: &[TopoShape], tools: &[TopoShape]) -> Result<TopoShape, String> {
    builder_bop(objects, tools, BoolOp2::Cut)
}

/// Intersection of the two groups. Source: `BOPAlgo_COMMON`.
pub fn common(objects: &[TopoShape], tools: &[TopoShape]) -> Result<TopoShape, String> {
    builder_bop(objects, tools, BoolOp2::Common)
}

/// Stable identity key of a shape — the address of its shared `TShape`
/// (the same identity `TopTools_ShapeMapHasher` uses).
fn shape_key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Bakes a shape's `TopoShape` location into its geometry, returning a copy at
/// the world position with an identity location. Shapes without a location are
/// returned unchanged (the shared `TShape` is preserved).
///
/// The Phase 19/20 intersection and rebuild stages read geometry from the
/// registry without applying the `TopoShape` location, so a located argument
/// (e.g. `crate::transform::translated`) would otherwise be intersected in its
/// local frame. Flattening the location at the builder boundary keeps the rest
/// of the pipeline location-agnostic.
fn flatten_location(s: &TopoShape) -> TopoShape {
    if s.location.is_identity() {
        return s.clone();
    }
    let t = s.location.transformation();
    crate::shape_ops::transformed_copy(s, &t).unwrap_or_else(|_| s.clone())
}

/// True if the pave block `pb` lies on a common block of the data structure:
/// its edge (or original edge) is an edge of a common block and its parameter
/// range matches one of the common-block ranges.
fn is_common_block_on_edge(ds: &BopdsDS, pb: &BopdsPaveBlock) -> bool {
    ds.common_blocks().iter().any(|cb| {
        let on_edge = cb.contains_index(pb.edge()) || cb.contains_index(pb.original_edge());
        on_edge && cb.contains_range(pb.first, pb.last, 1e-7)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::GpPnt;

    /// A second unit box at `x ∈ [3, 4]` — disjoint from the `[0, 1]³` box.
    fn far_box() -> TopoShape {
        BRepPrimBox::make_box_corner(&GpPnt::new(3.0, 0.0, 0.0), &GpPnt::new(4.0, 1.0, 1.0)).solid.0
    }

    /// Counts the direct solid children of a compound.
    fn count_solids(shape: &TopoShape) -> usize {
        shape
            .tshape
            .read()
            .unwrap()
            .children
            .iter()
            .filter(|h| h.shape_type() == ShapeType::Solid)
            .count()
    }

    #[test]
    fn defaults() {
        let b = BopBuilder::new();
        assert!(b.arguments().is_empty());
        assert!(!b.has_errors());
        assert!(b.errors().is_empty());
        assert!(b.warnings().is_empty());
        assert!(b.history().is_empty());
        assert!(b.origins().is_empty());
        assert!(b.shapes_sd().is_empty());
        assert_eq!(b.ds().nb_shapes(), 0);
        assert!(b.result().is_compound());
    }

    #[test]
    fn set_arguments_dedupes_by_tshape() {
        let a = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), a.solid.0.clone()]);
        assert_eq!(b.arguments().len(), 1);
        b.add_argument(&a.solid.0);
        assert_eq!(b.arguments().len(), 1);
        // A freshly built box is a distinct TShape and is added.
        let c = unit_box();
        b.add_argument(&c.solid.0);
        assert_eq!(b.arguments().len(), 2);
    }

    #[test]
    fn too_few_arguments_fails() {
        let mut b = BopBuilder::new();
        let err = b.perform().unwrap_err();
        assert!(err.contains("too few"), "err: {err}");
        assert!(b.has_errors());
    }

    #[test]
    fn perform_on_two_overlapping_boxes_succeeds() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        let result = b.perform().unwrap();
        assert!(!b.has_errors(), "errors: {:?}", b.errors());
        assert!(result.is_compound() || result.is_solid(), "type {:?}", result.shape_type());
        // The GF result is a compound of the two (unmodified) solids while the
        // face/solid stages are stubbed; the history still recorded the
        // vertex/edge images produced by the intersection.
        assert!(!b.history().is_empty(), "history should record split images");
        assert!(!b.origins().is_empty(), "origins back-map should be populated");
    }

    #[test]
    fn prepare_history_fills_modified_generated_removed() {
        // Fuse of two identical boxes: the coincident vertex/edge/face images
        // become modified relations, all kept in the result (`PrepareHistory`,
        // `BOPAlgo_Builder_4.cxx`).
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.perform().unwrap();
        let m = b.history().modified_map();
        assert!(!m.is_empty(), "identical boxes must record modified relations");
        let result_sub: Vec<TopoShape> = all_subshapes(b.result());
        for splits in m.values() {
            for sp in splits {
                assert!(
                    result_sub.iter().any(|r| r.same_tshape(sp)),
                    "a modified split must be kept in the result"
                );
            }
        }
        // Any generated relation also references a result shape.
        for gens in b.history().generated_map().values() {
            for gg in gens {
                assert!(result_sub.iter().any(|r| r.same_tshape(gg)));
            }
        }
        // Cut of a disjoint tool: the tool solid has no trace in the result and
        // no surviving splits, so it is marked removed.
        let a = unit_box();
        let bx = far_box();
        let mut b2 = BopBuilder::new();
        b2.set_arguments(&[a.solid.0.clone(), bx.clone()]);
        let (os, ts) = BoolOp2::Cut.states();
        b2.perform_internal(&[a.solid.0.clone()], os, &[bx.clone()], ts).unwrap();
        assert!(b2.history().is_deleted(&bx), "the cut-away tool solid is removed");
    }

    #[test]
    fn fill_images_vertices_maps_same_domain_corners() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // The 8 coincident corners of the two identical boxes are same-domain.
        assert!(!b.ds().shapes_sd().is_empty(), "coincident corners must be same-domain");
        b.fill_images_vertices().unwrap();
        let sd = b.ds().shapes_sd().clone();
        for (&n_v, &n_vsd) in &sd {
            let v = b.ds().shape(n_v).expect("vertex shape");
            let vsd = b.ds().shape(n_vsd).expect("sd vertex shape");
            let img = b.history().image(v).expect("vertex has an image");
            assert!(!img.is_empty());
            assert!(img[0].same_tshape(vsd), "image must be the same-domain vertex");
            let ors = b.origins().get(&shape_key(vsd)).expect("origin recorded");
            assert!(ors.iter().any(|o| o.same_tshape(v)), "origins back-map must point at the source");
        }
    }

    #[test]
    fn fill_images_edges_records_split_edges() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        b.fill_images_edges().unwrap();
        // Every source edge that carries pave blocks has a recorded image
        // (its split piece — the whole edge when no split was needed).
        let n = b.ds().nb_source_shapes();
        let mut n_edge_images = 0;
        for i in 0..n {
            let Some(si) = b.ds().shape_info(i) else { continue };
            if si.shape_type() != ShapeType::Edge {
                continue;
            }
            if !b.ds().has_pave_blocks(i) {
                continue;
            }
            let e = si.shape();
            let img = b.history().image(e).expect("edge has an image");
            assert!(!img.is_empty());
            n_edge_images += 1;
        }
        assert!(n_edge_images >= 24, "12 box edges × 2 boxes = 24 source edges");
    }

    #[test]
    fn builder_bop_fuse_separate_boxes_two_bodies() {
        let a = unit_box();
        let bx = far_box();
        let result = builder_bop(&[a.solid.0.clone()], &[bx], BoolOp2::Fuse).unwrap();
        assert!(result.is_compound() || result.is_solid());
        assert_eq!(count_solids(&result), 2, "fuse of two disjoint boxes keeps both solids");
    }

    #[test]
    fn builder_bop_cut_and_common_run_without_error() {
        let a = unit_box();
        let bx = far_box();
        let cut = builder_bop(&[a.solid.0.clone()], &[bx.clone()], BoolOp2::Cut).unwrap();
        assert!(!cut.is_null());
        let cmn = builder_bop(&[a.solid.0.clone()], &[bx], BoolOp2::Common).unwrap();
        assert!(!cmn.is_null());
    }

    #[test]
    fn history_and_origins_are_accessible() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.perform().unwrap();
        assert!(!b.history().images().is_empty());
        // history_mut allows the sibling Phase-20 modules to add images.
        b.history_mut().add_image(&a.solid.0, c.solid.0.clone());
        assert!(b.history().has_image(&a.solid.0));
    }

    #[test]
    fn trait_support_methods_are_available() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        // Options are forwarded from the PaveFiller.
        assert_eq!(b.fuzzy_value(), 1e-7);
        assert!(!b.non_destructive());
        // origins_mut exposes the back-map for the sibling modules.
        b.origins_mut().entry(shape_key(&c.solid.0)).or_default().push(a.solid.0.clone());
        assert!(b.origins().get(&shape_key(&c.solid.0)).is_some());
        // Same-domain binding round-trips through the DS.
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        let sd = b.ds().shapes_sd().clone();
        if let Some((&n, &m)) = sd.iter().next() {
            let s = b.ds().shape(n).unwrap().clone();
            let rep = b.ds().shape(m).unwrap().clone();
            b.bind_shapes_sd(s.clone(), rep.clone());
            let found = b.seek_shapes_sd(&s).expect("same-domain representative");
            assert!(found.same_tshape(&rep));
        }
    }

    #[test]
    fn bool_op_states_follow_occt_conversion() {
        assert_eq!(BoolOp2::Fuse.states(), (FaceState::Out, FaceState::Out));
        assert_eq!(BoolOp2::Cut.states(), (FaceState::Out, FaceState::In));
        assert_eq!(BoolOp2::Common.states(), (FaceState::In, FaceState::In));
    }

    #[test]
    fn build_bop_rejects_unknown_shapes_and_bad_states() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // A shape that is not an argument of the operation is rejected.
        let err = b
            .build_bop(&[c.solid.0.clone()], FaceState::Out, &[], FaceState::Out)
            .unwrap_err();
        assert!(err.contains("unknown shape"), "err: {err}");

        // An invalid state is rejected. Use a fresh builder — the failure above
        // left the report dirty, and `build_bop` short-circuits on a dirty
        // report before validating the states.
        let mut b2 = BopBuilder::new();
        b2.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b2.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b2.filler_mut().perform().unwrap();
        let err2 = b2
            .build_bop(&[a.solid.0.clone()], FaceState::On, &[], FaceState::Out)
            .unwrap_err();
        assert!(err2.contains("invalid state"), "err: {err2}");
    }

    #[test]
    fn clear_resets_state_keeps_nothing() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.perform().unwrap();
        assert!(!b.history().is_empty());
        b.clear();
        assert!(b.arguments().is_empty());
        assert!(b.history().is_empty());
        assert!(!b.has_errors());
    }

    #[test]
    fn bop_builder_implements_build_ops_and_like() {
        // The two Phase-20 host traits must be implemented and reachable
        // through the generic entry points.
        fn via_ops<B: BopBuildOps>(f: &mut B) -> usize {
            f.ds().nb_shapes() + f.fuzzy_value() as usize
        }
        fn via_like<B: BopBuilderLike>(f: &mut B) -> bool {
            f.has_errors()
        }
        let mut b = BopBuilder::new();
        assert_eq!(via_ops(&mut b), 0, "empty DS, fuzzy 1e-7 floors to 0");
        assert!(!via_like(&mut b));

        // The mutable accessors of both traits delegate to the builder fields.
        let a = unit_box();
        b.history_mut().add_image(&a.solid.0, a.solid.0.clone());
        assert!(b.history().has_image(&a.solid.0));
        b.origins_mut().entry(shape_key(&a.solid.0)).or_default().push(a.solid.0.clone());
        assert!(b.origins().get(&shape_key(&a.solid.0)).is_some());
    }

    #[test]
    fn fill_images_faces_callable_through_trait() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // The face stage runs directly through the trait surface.
        crate::bop_build_faces::fill_images_faces(&mut b).unwrap();
        assert!(!b.has_errors(), "errors: {:?}", b.errors());
    }

    #[test]
    fn fill_images_solids_callable_through_trait() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // The solid stage runs directly through the trait surface.
        crate::bop_build_solids::fill_images_solids(&mut b).unwrap();
        assert!(!b.has_errors(), "errors: {:?}", b.errors());
    }

    #[test]
    fn perform_on_two_make_box_boxes_stays_valid() {
        // The validation case: two overlapping boxes (built with
        // BRepPrimBox::make_box) through the full General Fuse pipeline, with
        // the face and solid stages wired. Two coincident unit boxes: the
        // split-solid stage rebuilds at least one solid image whose union
        // covers the coincident pair.
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let c = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut br = BopBuilder::new();
        br.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        let result = br.perform().unwrap();
        assert!(!br.has_errors(), "errors: {:?}", br.errors());
        assert!(result.is_compound() || result.is_solid(), "type {:?}", result.shape_type());
        assert!(!result.is_null());
        // At least one argument solid was rebuilt into an image by the
        // split-solid stage (coincident faces produce section edges), and the
        // pipeline completes without errors.
        let any_image = br
            .arguments()
            .iter()
            .any(|arg| arg.is_solid() && br.history().has_image(arg));
        assert!(any_image, "a solid image was recorded");
    }
}
