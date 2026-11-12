use super::*;

impl GeometryRegistry {
    /// The shared registry. Shapes and geometry live for the whole process,
    /// matching OCCT's reference-counted Handle model.
    pub fn global() -> &'static GeometryRegistry {
        static REGISTRY: OnceLock<GeometryRegistry> = OnceLock::new();
        REGISTRY.get_or_init(|| GeometryRegistry {
            ids: RwLock::new(HashMap::new()),
            face_surfaces: RwLock::new(HashMap::new()),
            surface_by_ptr: RwLock::new(HashMap::new()),
        })
    }

    // ---- vertices ----

    pub fn set_vertex(&self, s: &TopoShape, geom: VertexGeom) {
        let k = key(s);
        self.ids.write().unwrap().insert(k, shape_id(s));
        // T-25 batch 3: the vertex geometry lives on the shape itself now.
        let mut ts = s.tshape.write().unwrap();
        let c = ts.vertex_core_mut();
        c.point = geom.point;
        c.tolerance = geom.tolerance;
        drop(ts);
    }

    pub fn vertex_geom(&self, s: &TopoShape) -> Option<VertexGeom> {
        // T-25 batch 3: read off the vertex's own `TShape`; the side table is
        // only a fallback for shapes registered before this batch.
        // T-25 end phase: the vertex geometry comes off its own `TShape`.
        let c = s.tshape.read().unwrap().vertex_core().cloned()?;
        Some(VertexGeom {
            point: c.point,
            tolerance: c.tolerance,
        })
    }

    /// The vertex point, or the origin when unregistered (placeholder fallback).
    pub fn vertex_point(&self, s: &TopoShape) -> GpPnt {
        self.vertex_geom(s)
            .map(|g| g.point)
            .unwrap_or_else(GpPnt::zero)
    }

    pub fn vertex_tolerance(&self, s: &TopoShape) -> f64 {
        // A shape always carries a `BRep_TVertex`, whose tolerance starts at
        // `RealEpsilon()` (`BRep_TVertex.cxx:27`).
        self.vertex_geom(s)
            .map(|g| g.tolerance)
            .unwrap_or(occt_core::precision::REAL_EPSILON)
    }

    // ---- edges ----

    pub fn set_edge(&self, s: &TopoShape, geom: EdgeGeom) {
        let k = key(s);
        self.ids.write().unwrap().insert(k, shape_id(s));
        // T-25: drain the carrier onto the edge's own `TShape`. A geom that
        // carries no pcurve must NOT clear the slot (clearing is the explicit
        // `remove_pcurves_on_surface`, `ShapeBuild_Edge::RemovePCurve`).
        //
        // One guard covers both slots: they belong to the same `BRep_TEdge`.
        let mut ts = s.tshape.write().unwrap();
        if !(geom.pcurves.is_empty() && geom.pcurve_ranges.is_empty()) {
            let mut locations = OrderedSlots::default();
            for (face_key, _) in geom.pcurves.iter() {
                locations.insert(*face_key, TopLocLocation::identity());
            }
            ts.edge_pcurves = Some(EdgePcurves {
                curves: geom.pcurves.clone(),
                ranges: geom.pcurve_ranges.clone(),
                locations,
            });
        }
        let c = ts.edge_core_mut();
        c.curve = Some(geom.curve.clone());
        c.curve3d_null = geom.curve3d_null;
        c.first = geom.first;
        c.last = geom.last;
        c.tolerance = geom.tolerance;
        c.same_parameter = geom.same_parameter;
        c.same_range = geom.same_range;
        c.degenerated = geom.degenerated;
    }

    /// `BRep_Builder::Range(E, First, Last)`: set the edge's parameter range.
    /// `BRepPrim_OneAxis` relies on it for the degenerate pole edges
    /// (`SetParameters(ETOP/EBOTTOM, …, 0., myAngle)`,
    /// `BRepPrim_OneAxis.cxx:407`, `:418`): the 3D curve is absent but the range
    /// still spans the full period so the edge's pcurve is a full u-isoline.
    pub fn set_edge_range(&self, s: &TopoShape, first: f64, last: f64) {
        let mut ts = s.tshape.write().unwrap();
        let c = ts.edge_core_mut();
        c.first = first;
        c.last = last;
    }
    /// `BRep_Builder::Range(E, First, Last, Only3d=false)`
    /// (`BRep_Builder.cxx:1091-1117`): the range is written to *every*
    /// `BRep_GCurve` representation, i.e. the 3D curve and all pcurves. The
    /// `Only3d = true` variant is [`GeometryRegistry::set_edge_range`].
    pub fn set_ranges_all(&self, s: &TopoShape, first: f64, last: f64) {
        {
            let mut ts = s.tshape.write().unwrap();
            let c = ts.edge_core_mut();
            c.first = first;
            c.last = last;
            // Both the 3D range and every pcurve range of this edge.
            if let Some(pc) = ts.edge_pcurves.as_mut() {
                pc.ranges.for_each_mut(|_, r| *r = (first, last));
            }
        }
    }

    /// `BRep_Builder::UpdateEdge(E, C, Tol)` for the 3D curve
    /// (`BRep_Builder.cxx:1180-1191`): `TE->ChangeCurves()` replaces the
    /// `BRep_Curve3D` and `TE->UpdateTolerance(Tol)` raises the edge tolerance
    /// without lowering it (`BRep_TEdge.lxx:33-37`). The range and the pcurves
    /// are left alone.
    pub fn set_edge_curve3d(&self, s: &TopoShape, curve: Arc<dyn Curve>, tol: f64) {
        let mut ts = s.tshape.write().unwrap();
        let c = ts.edge_core_mut();
        c.curve = Some(curve);
        // `BRep_Builder::UpdateCurves(curves, C, L)` stores the curve in the
        // existing `BRep_Curve3D` representation (`BRep_Builder.cxx:80-86`), so
        // the representation is no longer the null one.
        c.curve3d_null = false;
        if tol > c.tolerance {
            c.tolerance = tol;
        }
    }

    /// `BRep_Builder::UpdateCurves(curves, C, L)` with a *null* curve: the
    /// `BRep_Curve3D` representation stays in the list holding a null curve
    /// (`BRep_Builder.cxx:86-92`), which is what
    /// `BRep_Builder::Degenerated(E, true)` installs (`cxx:1073-1085`) and what
    /// `ShapeBuild_Edge::MakeEdge(edge, C2d, S, L, p1, p2)`
    /// (`ShapeBuild_Edge.cxx:851-870`) leaves behind for a pcurve-only edge.
    /// `BRep_Tool::Curve` then answers null (`BRep_Tool.cxx:182-190`).
    pub fn set_curve3d_null(&self, s: &TopoShape, v: bool) {
        let mut ts = s.tshape.write().unwrap();
        ts.edge_core_mut().curve3d_null = v;
    }

    /// `ShapeBuild_Edge::RemoveCurve3d` (`ShapeBuild_Edge.cxx:704-710`):
    /// `B.UpdateEdge(edge, nullCurve, 0.)` clears the curve of the existing
    /// `BRep_Curve3D` representation and leaves its range in place
    /// (`BRep_Builder::UpdateCurves`, `BRep_Builder.cxx:58-102`).
    pub fn remove_edge_curve3d(&self, s: &TopoShape) {
        let mut ts = s.tshape.write().unwrap();
        let c = ts.edge_core_mut();
        c.curve = None;
        c.curve3d_null = true;
    }

    /// `BRep_TEdge::UpdateTolerance(Tol)` (`BRep_TEdge.lxx:33-37`): raise only.
    pub fn raise_edge_tolerance(&self, s: &TopoShape, tol: f64) {
        let mut ts = s.tshape.write().unwrap();
        let c = ts.edge_core_mut();
        if tol > c.tolerance {
            c.tolerance = tol;
        }
    }

    /// Number of `BRep_CurveRepresentation` entries that are a
    /// `BRep_CurveOnSurface` (`howMuchPCurves`, `ShapeFix_Wire.cxx:2262-2281`).
    /// A seam edge's `BRep_CurveOnClosedSurface` is one entry.
    pub fn pcurve_rep_count(&self, s: &TopoShape) -> usize {
        let ts = s.tshape.read().unwrap();
        ts.edge_pcurves().map(|g| g.curves.len()).unwrap_or(0)
    }

    pub fn edge_geom(&self, s: &TopoShape) -> Option<EdgeGeom> {
        // T-25 end phase: everything comes off the edge's own `TShape`; the
        // side table no longer stores edge geometry.
        let ts = s.tshape.read().unwrap();
        let c = ts.edge_core()?;
        Some(EdgeGeom {
            // An `EdgeGeom` always carries a 3D curve; an edge registered
            // without one is not representable and is reported as absent.
            // On a pcurve-only edge the curve is the `Adaptor3d_CurveOnSurface`
            // `BRepAdaptor_Curve(E, F)` builds, and `curve3d_null` carries the
            // `BRep_Tool::Curve` answer (`tshape.rs` `EdgeGeomCore`).
            curve: c.curve.clone()?,
            curve3d_null: c.curve3d_null,
            first: c.first,
            last: c.last,
            tolerance: c.tolerance,
            same_parameter: c.same_parameter,
            same_range: c.same_range,
            degenerated: c.degenerated,
            pcurves: ts
                .edge_pcurves()
                .map(|p| p.curves.clone())
                .unwrap_or_default(),
            pcurve_ranges: ts
                .edge_pcurves()
                .map(|p| p.ranges.clone())
                .unwrap_or_default(),
        })
    }

    /// The underlying curve handle (clone of the Arc).
    ///
    /// `BRep_Tool::Curve(E, L, First, Last)` (`BRep_Tool.cxx:172-196`): the
    /// curve stored in the edge's `BRep_Curve3D` representation, which is a null
    /// handle on a pcurve-only edge (`EdgeGeomCore::curve3d_null`). Callers that
    /// want what `BRepAdaptor_Curve(E, F)` evaluates read
    /// [`GeometryRegistry::edge_geom`] instead.
    pub fn edge_curve(&self, s: &TopoShape) -> Option<Arc<dyn Curve>> {
        // Read the 3D slot directly: `edge_geom` also clones both pcurve maps
        // and all six `BRep_TEdge` scalars, which this answer does not use.
        let ts = s.tshape.read().unwrap();
        let c = ts.edge_core()?;
        if c.curve3d_null {
            return None;
        }
        c.curve.clone()
    }

    pub fn edge_parameters(&self, s: &TopoShape) -> (f64, f64) {
        // `BRep_TEdge` always carries `myRange`; the 3D curve is a *separate*
        // CurveRepresentation (`BRep_Curve3D`) and may be absent.
        // `ShapeFix_ComposeShell::SplitByLine` builds exactly such an edge
        // (`cxx:2061-2070`: `MakeEdge` + two pcurves + `Range`), so the range
        // must be read from the edge core even when `edge_geom` reports no
        // curve. Returning `(-inf, inf)` here made `set_pcurves` /
        // `set_pcurve_range` consumers store an infinite COS window and the
        // seam edge never meshed.
        let ts = s.tshape.read().expect("poisoned TShape lock");
        match ts.edge_core() {
            Some(c) => (c.first, c.last),
            None => (f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    pub fn edge_tolerance(&self, s: &TopoShape) -> f64 {
        // `BRep_Tool::Tolerance(edge)` -> `BRep_TEdge::Tolerance()`, initial
        // value `RealEpsilon()` (`BRep_TEdge.cxx:33`).
        s.tshape
            .read()
            .unwrap()
            .edge_core()
            .map(|c| c.tolerance)
            .unwrap_or(occt_core::precision::REAL_EPSILON)
    }

    pub fn same_parameter(&self, s: &TopoShape) -> bool {
        let ts = s.tshape.read().unwrap();
        ts.edge_core().map(|c| c.same_parameter).unwrap_or(false)
    }

    /// `BRep_Tool::SameRange(E)`: the flag, not the ranges themselves.
    pub fn same_range(&self, s: &TopoShape) -> bool {
        let ts = s.tshape.read().unwrap();
        ts.edge_core().map(|c| c.same_range).unwrap_or(false)
    }

    pub fn is_degenerated_edge(&self, s: &TopoShape) -> bool {
        let ts = s.tshape.read().unwrap();
        ts.edge_core().map(|c| c.degenerated).unwrap_or(false)
    }

    /// `BRep_TVertex::UpdateTolerance(T)` (`BRep_TVertex.cxx:29-35`): raise
    /// only. Reached from `BRep_Builder::UpdateVertex(V, Par, E, Tol)`
    /// (`BRep_Builder.cxx:1220-1311`).
    pub fn raise_vertex_tolerance(&self, s: &TopoShape, tol: f64) {
        let mut ts = s.tshape.write().unwrap();
        let v = ts.vertex_core_mut();
        if tol > v.tolerance {
            v.tolerance = tol;
        }
    }

    /// `BRep_TVertex::Tolerance(const Standard_Real)` (`BRep_TVertex.hxx:64`):
    /// a direct write, not the raising `UpdateTolerance`. `BRepFill_Sweep`'s
    /// `Filling` uses it to pin the last vertex of an extended seam edge to
    /// `Precision::Confusion()` (`BRepFill_Sweep.cxx:1203`).
    pub fn set_vertex_tolerance(&self, s: &TopoShape, tol: f64) {
        let mut ts = s.tshape.write().unwrap();
        let core = ts.vertex_core_mut();
        core.tolerance = tol;
    }

    /// `BRep_Builder::UpdateVertex(V, Par, E, Tol)` (`BRep_Builder.cxx:1220-1311`)
    /// and its single-parameter write
    /// `BRep_Builder::UpdateVertex(V, Par, E)` (`BRep_Builder.hxx:392-394` ->
    /// `BRep_TVertex::UpdateParameter(V, E, Par)`): bind vertex `v` to
    /// parameter `par` on edge `s`.
    ///
    /// The port keeps the `(vertex, parameter)` pairs on the edge's
    /// `BRep_TEdge` core (`EdgeGeomCore::vertex_params`), which is where
    /// [`crate::brep_tool::BRepTool::parameter_of_vertex`] reads them back
    /// from, so a write here is visible to `BRep_Tool::Parameter(V, E)`.
    pub fn set_vertex_parameter_on_edge(&self, s: &TopoShape, v: &TopoShape, par: f64) {
        let key = key(v);
        let mut ts = s.tshape.write().unwrap();
        let core = ts.edge_core_mut();
        if let Some(slot) = core.vertex_params.iter_mut().find(|(k, _)| *k == key) {
            slot.1 = par;
        } else {
            core.vertex_params.push((key, par));
        }
    }

    /// `BRep_Builder::Degenerated(E, D)` (`BRep_Builder.cxx:1073-1087`): set
    /// the degenerated flag, and when `D` is true drop the 3D curve
    /// (`UpdateCurves(TE->ChangeCurves(), occ::handle<Geom_Curve>(), ...)`,
    /// `cxx:1082-1084`). `UpdateCurves` nulls an existing `BRep_Curve3D`
    /// representation in place and, when the edge has none, appends one holding
    /// the null handle with the range of the last `BRep_GCurve` (`cxx:60-92`);
    /// either way the representation stays present and `BRep_Tool::Curve`
    /// answers null. That presence is load-bearing for
    /// `BRepCheck_Edge::Minimum` (`BRepCheck_Edge.cxx:100-147`): `exist` is set
    /// by the representation, so the null curve raises no `No3DCurve`, and
    /// `myCref` stays null, so a degenerated edge raises no
    /// `InvalidDegeneratedFlag`. `EdgeGeom::curve` is a non-nullable
    /// `Arc<dyn Curve>`, so the null state is carried by
    /// `EdgeGeomCore::curve3d_null`, the flag `edge_curve` already answers
    /// `BRep_Tool::Curve` with; the stored curve keeps whatever a caller that
    /// needs `BRepAdaptor_Curve(edge, face)` would rebuild from the pcurve.
    ///
    /// `D = false` only clears the flag: OCCT restores no curve either.
    pub fn set_degenerated(&self, s: &TopoShape, v: bool) {
        let mut ts = s.tshape.write().unwrap();
        let c = ts.edge_core_mut();
        c.degenerated = v;
        if v {
            c.curve3d_null = true;
        }
    }
}
