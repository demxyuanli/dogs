use super::prelude::*;
use super::*;

impl BopBuilder {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            filler: PaveFiller::new(),
            history: BopHistory::new(),
            origins: HashMap::new(),
            shapes_sd: HashMap::new(),
            sd_by_key: HashMap::new(),
            arguments: Vec::new(),
            objects: Vec::new(),
            tools: Vec::new(),
            obj_state: FaceState::Out,
            tools_state: FaceState::Out,
            errors: Vec::new(),
            warnings: Vec::new(),
            result_shape: TopoShape::new(ShapeType::Compound),
            in_parts: HashMap::new(),
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

    /// Sets the additional tolerance forwarded to the PaveFiller.
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.filler.set_fuzzy_value(v);
    }

    /// True when the PaveFiller runs in non-destructive mode.
    pub fn non_destructive(&self) -> bool {
        self.filler.non_destructive()
    }

    /// Binds `shape` as same-domain with `sd`, mirroring `myShapesSD.Bind`.
    ///
    /// The shapes are resolved to their data-structure indices and the binding
    /// is recorded. OCCT `myShapesSD.Bind(aF, *pFSD)` also binds the
    /// representative to itself (`IsBound` is true); `BuildDraftSolid` uses
    /// that to take the `IsSplitToReverseWithWarn` branch. Reconstructed faces
    /// that are not in the DS are still recorded by TShape key.
    pub fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
        let k = GeometryRegistry::shape_key(&shape);
        self.sd_by_key.insert(k, sd.clone());
        if let (Some(i), Some(j)) = (self.filler.ds().index(&shape), self.filler.ds().index(&sd)) {
            if i != j {
                self.shapes_sd.insert(i, j);
            }
        }
    }

    /// Returns the same-domain representative of `shape`, when bound
    /// (`myShapesSD.Seek`), following the SD chain.
    pub fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
        if let Some(sd) = self.sd_by_key.get(&GeometryRegistry::shape_key(shape)) {
            return Some(sd.clone());
        }
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

    /// Replaces the result compound (`myShape = ...` after BuildRC / BuildBOP).
    pub fn set_result_shape(&mut self, shape: TopoShape) {
        self.result_shape = shape;
    }

    /// `myInParts` filled by FillIn3DParts.
    pub fn in_parts(&self) -> &HashMap<usize, Vec<TopoShape>> {
        &self.in_parts
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
    /// (`PaveFiller::perform`) followed by FillImages + BuildResult + post-treat.
    ///
    /// The result of the General Fuse is a compound containing all split parts
    /// of the arguments. Object/tool IN/OUT filtering is not applied here
    /// (`BOPAlgo_Builder::PerformInternal1`).
    pub fn perform(&mut self) -> Result<TopoShape, String> {
        let args = self.arguments.clone();
        self.perform_internal(&args, FaceState::Out, &[], FaceState::Out)
    }

    /// Shared pipeline of [`BopBuilder::perform`] and `BOPAlgo_BOP`.
    ///
    /// Mirrors `PerformInternal1`: FillImages (GF, no obj/tool filter) then
    /// [`Self::build_result`] (GF assembly of every argument image). When tools
    /// are present this is a Boolean run, so `BOPAlgo_BOP::BuildShape` applies
    /// Fuse/Cut/Common afterwards (`BOPAlgo_BOP.cxx:563`).
    pub(super) fn perform_internal(
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

        // GF building phase: record the BOP groups for BuildShape, then fill
        // images without IN/OUT selection (`FillImages*` does not read states).
        self.objects = objects.to_vec();
        self.tools = tools.to_vec();
        self.obj_state = obj_state;
        self.tools_state = tools_state;
        self.fill_images()?;
        self.build_result()?;
        if !self.tools.is_empty() {
            crate::bop_bop::build_shape(self)?;
        }
        self.prepare_history();
        self.post_treat()?;
        Ok(self.result_shape.clone())
    }

    /// Resets the per-run state (report, history, origins, same-domain map,
    /// result). The arguments and the filler options survive.
    pub(super) fn reset_run(&mut self) {
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
        self.in_parts.clear();
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
        self.in_parts.clear();
    }

    // -----------------------------------------------------------------------
    // Building the result
    // -----------------------------------------------------------------------

    /// Records object/tool groups and fills GF images.
    ///
    /// This is **not** `BOPAlgo_Builder::BuildBOP` (open-solid face selection,
    /// `crate::bop_build_bop`). It is the FillImages sequence from
    /// `PerformInternal1`. IN/OUT is validated here only so a direct caller
    /// cannot bind an illegal BOP state; FillImages itself does not read the
    /// states. Fuse/Cut/Common filtering is `crate::bop_bop::build_shape`.
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
        self.fill_images()
    }

    /// `FillImagesVertices` ... `FillImagesSolids` (`BOPAlgo_Builder.cxx:336-418`).
    /// No object/tool IN/OUT filter.
    pub(super) fn fill_images(&mut self) -> Result<(), String> {
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
    pub(super) fn fill_images_vertices(&mut self) -> Result<(), String> {
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
    /// back-map gets the source edge. A common-block member uses
    /// `RealPaveBlock` (`BOPDS_DS::RealPaveBlock`) so the image is the shared
    /// split, matching `BOPAlgo_Builder::FillImagesEdges` (`_1.cxx:71-126`).
    pub(super) fn fill_images_edges(&mut self) -> Result<(), String> {
        crate::pave_blocks::make_split_edges(&mut self.filler)?;
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
                let a_pbr = self.filler.ds().real_pave_block(pb);
                let n_sp_r = a_pbr.edge();
                if n_sp_r == 0 {
                    continue;
                }
                let Some(sp_r) = self.filler.ds().shape(n_sp_r).cloned() else { continue };
                self.history.add_image(&e, sp_r.clone());
                self.origins.entry(shape_key(&sp_r)).or_default().push(e.clone());
                if is_common_block_on_edge(self.filler.ds(), pb) {
                    let n_sp = pb.edge();
                    if n_sp != 0 && n_sp != n_sp_r {
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
    /// (`BuildSplitFaces` → `FillSameDomainFaces`). GF only; no obj/tool state.
    pub(super) fn fill_images_faces(&mut self) -> Result<(), String> {
        crate::bop_build_faces::fill_images_faces(self)
    }

    /// Fills the images of the solids — split solids of the arguments.
    ///
    /// FillIn3DParts → BuildSplitSolids → FillInternalShapes
    /// (`BOPAlgo_Builder_3.cxx:70`). GF only; no obj/tool state.
    pub(super) fn fill_images_solids(&mut self) -> Result<(), String> {
        if !crate::bop_images_solids::has_source_solids(self) {
            self.in_parts.clear();
            return Ok(());
        }
        let mut ctx = crate::int_tools_full::IntToolsContext::new();
        let fill = crate::bop_fill_in3d::fill_in_3d_parts_builder(self, &mut ctx);
        self.in_parts = fill.in_parts.clone();
        crate::bop_split_solids_occt::build_split_solids_occt(self, &fill)?;
        crate::bop_fill_internals_occt::fill_internal_shapes_occt(self)?;
        Ok(())
    }

    /// Assembles the result compound from the images of the arguments.
    ///
    /// Port of `BOPAlgo_Builder::BuildResult` (`Builder_1.cxx:130`) and
    /// `BOPAlgo_BOP::BuildResult` (`BOP.cxx:323`): for every argument of a
    /// given type, add its images, or the argument itself if it has none.
    /// Fence by `TShape` identity. No object/tool IN/OUT skip — Cut/Common
    /// membership is `BOPAlgo_BOP::BuildShape` / `BuildRC`.
    pub(super) fn build_result(&mut self) -> Result<TopoShape, String> {
        let b = TopoBuilder::new();
        let mut result = b.make_compound_of(&[]);
        let mut fence: Vec<TopoShape> = Vec::new();
        for t in RESULT_TYPES {
            for arg in &self.arguments {
                if arg.shape_type() != t {
                    continue;
                }
                let imgs: Vec<TopoShape> = match self.history.image(arg) {
                    Some(list) => list.to_vec(),
                    None => vec![arg.clone()],
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
    /// Port of `BOPAlgo_Builder::PostTreat` (`BOPAlgo_Builder.cxx:450`):
    /// `CorrectTolerances` + `CorrectShapeTolerances` with `aMA` filled from
    /// source VERTEX/EDGE/FACE when `NonDestructive` is set.
    pub(super) fn post_treat(&mut self) -> Result<(), String> {
        let mut avoid = HashSet::new();
        if self.non_destructive() {
            let n = self.filler.ds().nb_source_shapes();
            for i in 0..n {
                let Some(si) = self.filler.ds().shape_info(i) else {
                    continue;
                };
                let t = si.shape_type();
                if t == ShapeType::Vertex || t == ShapeType::Edge || t == ShapeType::Face {
                    avoid.insert(GeometryRegistry::shape_key(si.shape()));
                }
            }
        }
        AlgoTools::correct_tolerances_avoid(&self.result_shape, &avoid, 0.05);
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
    pub(super) fn prepare_history(&mut self) {
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
    pub(super) fn loc_generated(&self, s: &TopoShape) -> Vec<TopoShape> {
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
    pub(super) fn collect_interf_vertices(
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
