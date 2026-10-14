//! Resolver plumbing: caches, warnings, and the `ShapeFix_Shape` heal context (`ShapeFix_Shape::myContext` / `ShapeBuild_ReShape`). Source: `ShapeFix_Shape.cxx:75`, `:257`, `ShapeBuild_ReShape.cxx:191`.
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn new(records: &'a HashMap<usize, Record>) -> Self {
        Self {
            records,
            precision: step_precision(records),
            plane_angle_factor: context_plane_angle_factor(records),
            length_factor: context_length_factor(records),
            b: TopoBuilder::new(),
            shape_cache: RefCell::new(HashMap::new()),
            vertex_bind: RefCell::new(HashMap::new()),
            point_cache: RefCell::new(HashMap::new()),
            dir_cache: RefCell::new(HashMap::new()),
            axis_cache: RefCell::new(HashMap::new()),
            curve_cache: RefCell::new(HashMap::new()),
            surface_cache: RefCell::new(HashMap::new()),
            curve2d_cache: RefCell::new(HashMap::new()),
            surface_curve_pcurves: RefCell::new(HashMap::new()),
            edge_curve_ref: RefCell::new(HashMap::new()),
            resolving: RefCell::new(HashSet::new()),
            warnings: RefCell::new(Vec::new()),
            heal_context: RefCell::new(Vec::new()),
        }
    }

    pub(in crate::step) fn record(&self, id: usize) -> Result<&'a Record, String> {
        self.records
            .get(&id)
            .ok_or_else(|| format!("reference to undefined entity #{id}"))
    }

    pub(in crate::step) fn warn(&self, msg: String) {
        self.warnings.borrow_mut().push(msg);
    }

    /// The `ShapeBuild_ReShape` handle of the `FixShape` pass currently running
    /// (`ShapeFix_Shape::myContext`), or `None` when a face is resolved outside
    /// a shell/solid (no shared pass to take part in).
    ///
    /// `ShapeFix_Shape` hands one `ShapeBuild_ReShape` to every face of a pass
    /// (`ShapeFix_Shape.cxx:75`, `ShapeFix_Shell.cxx:108`,
    /// `ShapeFix_Face.cxx:379`) and closes it with `myResult = Context()->Apply(S)`
    /// (`ShapeFix_Shape.cxx:257`), so a face fix that splits/rebuilds an edge is
    /// substituted on the *neighbouring* faces too. Sharing is therefore what
    /// keeps one `TopoDS_Edge` TShape per geometric edge: without it the seam
    /// rebuild of one face leaves a coincident duplicate edge beside the
    /// neighbour's, and BRepMesh discretizes the two independently (T-junctions).
    ///
    /// ALWAYS ON (2026-10-07). The pass used to sit behind an `OCCT_SHARED_HEAL`
    /// switch while it was unverified; it is now the only path, because OCCT has
    /// no such switch. Verified this round on the 23-model gate table plus the
    /// eight `data/occ` models:
    ///  * T0M open edges (edges used by one triangle, welded at 1e-6) 3643 -> 250
    ///    (OCCT reference 4); model edges 5005 -> 4437 with single-pcurve edges
    ///    1298 -> 78 (the remaining pairs are the genuinely unshared STEP edges).
    ///  * a3n00 open edges 571 -> 0 (OCCT 0); acs10 646 -> 4 (OCCT 0).
    ///  * All 15 `data/step` exports are byte-identical to the previous default.
    ///  * `step_obj_gates` 4 passed; mesh areas move toward OCCT (T0M 192654.648 ->
    ///    192654.622 vs OCCT 192658.641).
    pub(in crate::step) fn current_heal_context(
        &self,
    ) -> Option<crate::shape_fix_compose_shell::SharedReShape> {
        self.heal_context.borrow().last().cloned()
    }

    /// `ShapeFix_Shape::Perform`'s shape scope: `FromSTEP.exec.op` = `FixShape`
    /// (`STEPControl_Controller.cxx:201`) runs `ShapeFix_Shape` once per
    /// transferred root with a single `ShapeBuild_ReShape`
    /// (`ShapeFix_Shape.cxx:75`), and the pass ends with
    /// `myResult = Context()->Apply(S)` (`cxx:257`) - the substitution of every
    /// sub-shape a face fix recorded, including on the faces it did not touch.
    /// Returns the owned context to `Apply` when this call opened the pass.
    pub(super) fn open_heal_context(&self) -> Option<crate::shape_fix_compose_shell::SharedReShape> {
        // See `current_heal_context`: the shared `FixShape` context is the only
        // path (OCCT has no switch for it).
        if self.heal_context.borrow().is_empty() {
            let ctx = crate::shape_fix_compose_shell::SharedReShape::new();
            self.heal_context.borrow_mut().push(ctx.clone());
            Some(ctx)
        } else {
            None
        }
    }

    /// Closes the pass opened by `open_heal_context`: pops the handle and
    /// materialises the substitutions on the shape (`ShapeFix_Shape.cxx:257`),
    /// then runs the rest of `ShapeFix_Shape::Perform`'s tail on the root -
    /// `SameParameter` (`cxx:259-261`) and `FixVertexTolerance` (`cxx:267-291`)
    /// - and the closing `Context()->Apply` (`cxx:294`). Both stage gates are
    /// the tool defaults `-1`: STEP registers `FixSameParameterMode`
    /// (`STEPControl_Controller.cxx:211`) and leaves `FixVertexTolMode` at
    /// `ShapeFix_Shape`'s ctor value, so `NeedFix(-1)` is true for both
    /// (`ShapeFix_Root.lxx:101-104`).
    pub(super) fn close_heal_context(
        &self,
        owned: &Option<crate::shape_fix_compose_shell::SharedReShape>,
        shape: TopoShape,
    ) -> TopoShape {
        match self.pop_heal_context(owned) {
            // cxx:257: myResult = Context()->Apply(S)
            Some(mut ctx) => {
                let mut result = ctx.apply(&shape);
                // `ShapeFix_Shape::SameParameter` (`cxx:307-311`) is
                // `ShapeFix::SameParameter(sh, enforce, 0.0)`; `enforce` is
                // false (`cxx:261`).
                crate::shhealing::shape_fix_same_parameter(&result, false, 0.0);
                crate::shhealing::fix_vertex_tolerances_after_same_parameter(&result, &mut ctx);
                // cxx:294.
                result = ctx.apply(&result);
                result
            }
            None => shape,
        }
    }

    pub(super) fn pop_heal_context(
        &self,
        owned: &Option<crate::shape_fix_compose_shell::SharedReShape>,
    ) -> Option<crate::shape_fix_compose_shell::SharedReShape> {
        if owned.is_some() {
            self.heal_context.borrow_mut().pop()
        } else {
            None
        }
    }

}
