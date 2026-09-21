use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// B-spline approximation helpers
// ---------------------------------------------------------------------------

/// Approximate a non-analytic surface with an interpolating B-spline surface.
///
/// Trait objects cannot be downcast, so a genuine `GeomBSplineSurface` is
/// recovered by sampling and re-fitting (`fit_surface_grid`). The fit passes
/// through every grid node, so a sampled B-spline reproduces itself closely;
/// analytic surfaces should never reach this path.
pub(super) fn fit_bspline_surface(s: &dyn Surface) -> Result<GeomBSplineSurface, String> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let nu = 6;
    let nv = 6;
    let mut pts = Vec::with_capacity(nu);
    for i in 0..nu {
        let mut row = Vec::with_capacity(nv);
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            row.push(s.d0(u, v));
        }
        pts.push(row);
    }
    occt_geom::bspline_surface::fit_surface_grid(&pts, 2, 2)
        .or_else(|_| occt_geom::bspline_surface::fit_surface_grid(&pts, 1, 1))
        .map_err(|e| format!("fit_bspline_surface: {e}"))
}

/// `GeomConvert::CurveToBSplineCurve` for a Bezier curve
/// (`GeomToStep_MakeBoundedCurve.cxx:64-75` converts every `Geom_BezierCurve`
/// this way before writing it): a Bezier is exactly a clamped B-spline whose
/// degree is `nb_poles - 1` and whose knots are `0`/`1` with multiplicity
/// `degree + 1`. The port's `GeomBezierCurve` is non-rational, matching
/// `Geom_BezierCurve::IsRational() == false` for this case.
fn bezier_to_bspline(poles: &[GpPnt]) -> Option<GeomBSplineCurve> {
    let degree = poles.len().checked_sub(1)?;
    let mut knots = vec![0.0; degree + 1];
    knots.extend(std::iter::repeat(1.0).take(degree + 1));
    GeomBSplineCurve::new(poles.to_vec(), knots, degree).ok()
}

/// Internal writer state: the entity writer plus identity maps so shared
/// vertices/edges/curves are emitted exactly once.
pub(super) struct WriteCtx {
    pub(super) w: StepWriter,
    pub(super) vertex_ids: HashMap<usize, usize>,
    pub(super) edge_ids: HashMap<usize, usize>,
    pub(super) curve_ids: HashMap<usize, usize>,
    pub(super) geom_ctx: usize,
    pub(super) prod_ctx: usize,
    pub(super) def_ctx: usize,
    /// When true, non-analytic curves/surfaces are written as B-splines
    /// (sampled + fitted) instead of the tangent-line / plane fallbacks.
    pub(super) splines: bool,
}

impl WriteCtx {
    pub(super) fn new() -> Self {
        let mut w = StepWriter::new();
        let ctx = w.emit("APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN')".into());
        let prod_ctx = w.emit(format!("PRODUCT_CONTEXT('',#{ctx},'mechanical')"));
        let def_ctx = w.emit(format!("PRODUCT_DEFINITION_CONTEXT('part definition',#{ctx},'design')"));
        let geom_ctx = w.emit("GEOMETRIC_REPRESENTATION_CONTEXT('','',3)".into());
        Self {
            w,
            vertex_ids: HashMap::new(),
            edge_ids: HashMap::new(),
            curve_ids: HashMap::new(),
            geom_ctx,
            prod_ctx,
            def_ctx,
            splines: false,
        }
    }

    pub(super) fn finish(self) -> String {
        self.w.finish()
    }

    pub(super) fn emit_vertex(&mut self, v: &Vertex) -> usize {
        let key = Arc::as_ptr(&v.0.tshape) as usize;
        if let Some(&id) = self.vertex_ids.get(&key) {
            return id;
        }
        let p = GeometryRegistry::global().vertex_point(&v.0);
        let pid = self.w.add_cartesian_point(&p);
        let id = self.w.emit(format!("VERTEX_POINT('',#{pid})"));
        self.vertex_ids.insert(key, id);
        id
    }

    pub(super) fn emit_edge(&mut self, e: &Edge) -> usize {
        let key = Arc::as_ptr(&e.0.tshape) as usize;
        if let Some(&id) = self.edge_ids.get(&key) {
            return id;
        }
        let (v1, v2) = match edge_vertices(e) {
            (Some(a), Some(b)) => (a, b),
            _ => {
                // Degenerate edge without vertex children: synthesize from the
                // curve endpoints so the STEP file stays well-formed.
                let (a, b) = GeometryRegistry::global().edge_parameters(&e.0);
                let curve = GeometryRegistry::global().edge_curve(&e.0);
                let p1 = curve.as_ref().map(|c| c.d0(a)).unwrap_or_default();
                let p2 = curve.as_ref().map(|c| c.d0(b)).unwrap_or_default();
                let bld = TopoBuilder::new();
                (bld.make_vertex(p1, 0.0), bld.make_vertex(p2, 0.0))
            }
        };
        let vid1 = self.emit_vertex(&v1);
        let vid2 = self.emit_vertex(&v2);
        // `TopoDSToStep_MakeStepEdge.cxx:345` passes `MkCurve.Value()`, which is a
        // null handle when `GeomToStep_MakeCurve` left `done = false`
        // (`MakeCurve.cxx:100-103`); the STEP geometry attribute is then unset,
        // which this writer spells `$`.
        let geom = match self.emit_edge_curve(e) {
            Some(id) => format!("#{id}"),
            None => "$".to_string(),
        };
        let id = self.w.emit(format!("EDGE_CURVE('',#{vid1},#{vid2},{geom},.T.)"));
        self.edge_ids.insert(key, id);
        id
    }

    pub(super) fn emit_edge_curve(&mut self, e: &Edge) -> Option<usize> {
        // `TopoDSToStep_MakeStepEdge.cxx:194-196` branches on
        // `BRepAdaptor_Curve(aEdge).Curve()` being null, which is the case for a
        // **degenerate** edge (`BRep_Tool::Curve` has no 3D curve). This port
        // stores a point-curve stand-in for those, so the degeneracy flag is what
        // routes to the "edge without 3d curve; creating..." branch below.
        let degenerate = crate::brep_tool::BRepTool::is_degenerated(e)
            || GeometryRegistry::global().is_degenerated_edge(&e.0);
        let curve = if degenerate {
            None
        } else {
            GeometryRegistry::global().edge_curve(&e.0)
        };
        let Some(curve) = curve else {
            // `TopoDSToStep_MakeStepEdge.cxx:263-330` ("edge without 3d curve;
            // creating..."): OCCT builds a `Geom_Line` from the pcurve endpoints
            // when the current face is a plane and the pcurve a line, and
            // otherwise fits a B-spline through sampled surface points. This port
            // emits a line through the edge's vertices for both cases — the
            // B-spline fit is UNPORTED.
            let (v1, v2) = edge_vertices(e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                let pid = self.w.add_cartesian_point(&GpPnt::zero());
                let vid = self.w.add_vector(&dir_x(), 1.0);
                return Some(self.w.emit(format!("LINE('',#{pid},#{vid})")));
            };
            let p1 = GeometryRegistry::global().vertex_point(&v1.0);
            let p2 = GeometryRegistry::global().vertex_point(&v2.0);
            let d = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).unwrap_or(dir_x());
            let pid = self.w.add_cartesian_point(&p1);
            let vid = self.w.add_vector(&d, 1.0);
            return Some(self.w.emit(format!("LINE('',#{pid},#{vid})")));
        };
        let ckey = Arc::as_ptr(&curve) as *const () as usize;
        if let Some(&id) = self.curve_ids.get(&ckey) {
            return Some(id);
        }
        let id = self.emit_curve_entity(curve.as_ref());
        if let Some(id) = id {
            self.curve_ids.insert(ckey, id);
        }
        id
    }

    pub(super) fn emit_curve_entity(&mut self, c: &dyn Curve) -> Option<usize> {
        // `GeomToStep_MakeCurve.cxx:50-104` dispatch order:
        // `Geom_Line` → `Geom_Conic` → `Geom_TrimmedCurve` → `Geom_BoundedCurve`
        // → `done = false`. The family comes from the curve's own type
        // (`IsKind`), never from sampling: the previous body classified by six
        // `|d²|` samples with a `(max-min)/max < 0.02` threshold that OCCT has
        // nowhere (audit A3).
        if let Some(l) = c.gp_line() {
            let pid = self.w.add_cartesian_point(&l.location());
            let vid = self.w.add_vector(&l.direction(), 1.0);
            return Some(self.w.emit(format!("LINE('',#{pid},#{vid})")));
        }
        if let Some(circ) = c.gp_circ() {
            return Some(emit_circle_entity(&mut self.w, &circ));
        }
        if let Some(e) = c.gp_ellipse() {
            return Some(emit_ellipse_entity(&mut self.w, &e));
        }
        if let Some(h) = c.gp_hyperbola() {
            return Some(emit_hyperbola_entity(&mut self.w, &h));
        }
        if let Some(p) = c.gp_parabola() {
            return Some(emit_parabola_entity(&mut self.w, &p));
        }
        // `MakeCurve.cxx:66-92`: a `Geom_TrimmedCurve` is written through its
        // **basis** curve. For a conic basis the `gp_*` queries above already
        // returned the basis (the port's `GeomTrimmedCurve` forwards them), so
        // only the remaining non-BSpline/Bezier bases need the recursion here;
        // a BSpline/Bezier basis would be `Segment`-ed by OCCT (`cxx:71-82`) and
        // is represented by this port's own remapped knots instead, so it falls
        // through to the spline arm below.
        if c.is_geom_trimmed() {
            if let Some((basis, _bf, _bl)) = c.untrimmed_basis() {
                if basis.bspline_knots().is_none() && basis.bezier_poles().is_none() {
                    return self.emit_curve_entity(basis.as_ref());
                }
            }
        }
        // `MakeCurve.cxx:94-99`: `Geom_BoundedCurve` → `GeomToStep_MakeBoundedCurve`
        // (`GeomToStep_MakeBoundedCurve.cxx:37-80`): a B-spline goes to
        // `B_SPLINE_CURVE_WITH_KNOTS` (+`_AND_RATIONAL_...` when rational), and a
        // Bezier is first converted with `GeomConvert::CurveToBSplineCurve` and
        // written the same way.
        // UNPORTED: `Geom_BSplineCurve::SetNotPeriodic` (`cxx:46-51`) is not
        // ported, so a periodic B-spline keeps its periodicity here.
        if let (Some(poles), Some(knots), Some(deg)) =
            (c.bspline_poles(), c.bspline_knots(), c.nurbs_degree())
        {
            let bs = match c.bspline_weights() {
                Some(w) if w.len() == poles.len() => {
                    GeomBSplineCurve::rational(poles.to_vec(), w.to_vec(), knots.to_vec(), deg)
                }
                _ => GeomBSplineCurve::new(poles.to_vec(), knots.to_vec(), deg),
            };
            if let Ok(bs) = bs {
                if let Ok(id) = write_bspline_curve(&mut self.w, &bs) {
                    return Some(id);
                }
            }
        }
        if let Some(poles) = c.bezier_poles() {
            if let Some(bs) = bezier_to_bspline(poles) {
                if let Ok(id) = write_bspline_curve(&mut self.w, &bs) {
                    return Some(id);
                }
            }
        }
        // `GeomToStep_MakeCurve.cxx:100-103`: an unrecognised curve (not
        // `Geom_Line`, `Geom_Conic`, `Geom_TrimmedCurve` or `Geom_BoundedCurve`)
        // leaves `done = false` and `theCurve` null. The caller
        // (`TopoDSToStep_MakeStepEdge.cxx:260-261,345`) then builds the
        // `EDGE_CURVE` with that null handle, i.e. the geometry attribute is
        // simply unset — see `emit_edge`, which writes `$` for `None`. The
        // previous body emitted a spline fit or a tangent line instead, so the
        // file never showed the deviation (task T-71).
        None
    }

    pub(super) fn emit_wire(&mut self, w: &Wire) -> usize {
        let edges = edges_of_wire(w);
        let mut refs = Vec::with_capacity(edges.len());
        for e in &edges {
            let edge_id = self.emit_edge(e);
            // `TopoDSToStep_MakeStepWire.cxx:258-263` builds each EDGE_LOOP entry
            // as `StepShape_OrientedEdge` with a derived edgeStart/edgeEnd
            // (`RWStepShape_RWOrientedEdge.cxx:78-84` writes them as `*`) and the
            // boolean occurrence flag `(anEdge.Orientation() == TopAbs_FORWARD)`.
            // The edge's own vertices belong to the EDGE_CURVE emitted above, so
            // only the occurrence orientation is written here.
            let occ = if e.0.orientation().is_forward() { ".T." } else { ".F." };
            refs.push(self.w.emit(format!("ORIENTED_EDGE('',*,*,#{edge_id},{occ})")));
        }
        self.w.emit(format!("EDGE_LOOP('',({}))", join_refs(&refs)))
    }

    pub(super) fn emit_face(&mut self, f: &Face) -> usize {
        let surf_ref = self.emit_surface(f);
        let wires = wires_of_face(f);
        let mut bounds = Vec::with_capacity(wires.len());
        // `TopoDSToStep_MakeStepFace.cxx:297-327`: each wire of the FORWARD face
        // becomes a plain `StepShape_FaceBound`. OCCT's writer never builds a
        // `StepShape_FaceOuterBound` (the only `new StepShape_FaceOuterBound`
        // in the OCCT sources is the reader's `NewEntity`,
        // `RWStepAP214_GeneralModule.cxx:5756`); a DRAWEXE 8.0.0
        // `testwritestep` of `data/HoledPlate.step` confirms 38 `FACE_BOUND`
        // and 0 `FACE_OUTER_BOUND`.
        // The boolean is the wire orientation relative to the FORWARD face,
        // negated when the face itself is REVERSED (`cxx:319-325`):
        //   face FORWARD  -> `(CurrentWire.Orientation() == TopAbs_FORWARD)`
        //   face REVERSED -> `(CurrentWire.Orientation() == TopAbs_REVERSED)`
        let face_forward = f.0.orientation().is_forward();
        for w in &wires {
            let loop_ref = self.emit_wire(w);
            // `wires_of_face` composes the face orientation into the child
            // (`TopoDS_Iterator`), so undo it to recover the wire's stored
            // orientation relative to the face.
            let rel_forward = w.0.orientation().is_forward() == face_forward;
            let bound_forward = if face_forward { rel_forward } else { !rel_forward };
            let b = if bound_forward { ".T." } else { ".F." };
            bounds.push(self.w.emit(format!("FACE_BOUND('',#{loop_ref},{b})")));
        }
        // ISO 10303-42 `advanced_face(name, bounds, face_geometry, same_sense)`:
        // the surface is the third argument. OCCT's reader
        // (`RWStepShape_RWAdvancedFace`) and `StepToTopoDS_TranslateFace` expect
        // this order; an OCCT-written cylinder reads back only with it.
        // `TopoDSToStep_MakeStepFace.cxx:479`:
        // `Fpms->Init(aName, aBounds, Spms, aFace.Orientation() == TopAbs_FORWARD)`
        // -- the face's relative orientation in the shell, not a constant.
        let same_sense = if face_forward { ".T." } else { ".F." };
        self.w
            .emit(format!("ADVANCED_FACE('',({}),#{surf_ref},{same_sense})", join_refs(&bounds)))
    }

    pub(super) fn emit_surface(&mut self, f: &Face) -> usize {
        let Some(surf) = GeometryRegistry::global().face_surface(&f.0) else {
            return self.emit_plane_fallback();
        };
        match surface_kind(surf.as_ref()) {
            SurfKind::Plane => {
                let pln = face_plane(f).unwrap_or_else(|| GpPln::new(GpAx3::standard()));
                let ax = self.w.add_axis2_placement_3d(
                    &pln.location(),
                    &pln.axis().direction(),
                    &pln.x_axis().direction(),
                );
                self.w.emit(format!("PLANE('',#{ax})"))
            }
            SurfKind::Sphere => {
                let (u0, _, v0, v1) = surf_bounds(surf.as_ref());
                let vm = 0.5 * (v0 + v1);
                let center = sphere_center(surf.as_ref()).unwrap_or_else(GpPnt::zero);
                let r = surf.d0(u0, vm).distance(&center);
                let ax = self.w.add_axis2_placement_3d(&center, &dir_z(), &dir_x());
                self.w.emit(format!("SPHERICAL_SURFACE('',#{ax},{})", step_real(r)))
            }
            SurfKind::Cylinder => {
                if let Some((ax2, r)) = cylinder_params(surf.as_ref()) {
                    let ax = self.emit_ax2(&ax2);
                    self.w.emit(format!("CYLINDRICAL_SURFACE('',#{ax},{})", step_real(r)))
                } else {
                    self.emit_plane_fallback()
                }
            }
            SurfKind::Cone => {
                if let Some((ax2, r, a)) = cone_params(surf.as_ref()) {
                    let ax = self.emit_ax2(&ax2);
                    self.w.emit(format!(
                        "CONICAL_SURFACE('',#{ax},{},{})",
                        step_real(r),
                        step_real(a)
                    ))
                } else {
                    self.emit_plane_fallback()
                }
            }
            SurfKind::Torus => {
                if let Some((ax2, maj, min)) = torus_params(surf.as_ref()) {
                    let ax = self.emit_ax2(&ax2);
                    self.w.emit(format!(
                        "TOROIDAL_SURFACE('',#{ax},{},{})",
                        step_real(maj),
                        step_real(min)
                    ))
                } else {
                    self.emit_plane_fallback()
                }
            }
            SurfKind::Other => {
                // B-spline or otherwise non-analytic surface: emit a fitted
                // B-spline when spline output is requested, else a unit plane.
                if self.splines {
                    if let Ok(bs) = fit_bspline_surface(surf.as_ref()) {
                        if let Ok(id) = write_bspline_surface(&mut self.w, &bs) {
                            return id;
                        }
                    }
                }
                self.emit_plane_fallback()
            }
        }
    }

    pub(super) fn emit_ax2(&mut self, ax2: &GpAx2) -> usize {
        self.w
            .add_axis2_placement_3d(&ax2.location(), &ax2.direction(), &ax2.x_direction())
    }

    pub(super) fn emit_plane_fallback(&mut self) -> usize {
        let ax = self.w.add_axis2_placement_3d(&GpPnt::zero(), &dir_z(), &dir_x());
        self.w.emit(format!("PLANE('',#{ax})"))
    }

    pub(super) fn emit_shell(&mut self, sh: &Shell) -> usize {
        let faces = children_of_type(&sh.0, ShapeType::Face);
        let refs: Vec<usize> = faces
            .iter()
            .map(|f| self.emit_face(&Face(f.clone())))
            .collect();
        self.w
            .emit(format!("CLOSED_SHELL('',({}))", join_refs(&refs)))
    }

    pub(super) fn emit_solid(&mut self, s: &Solid) -> usize {
        let shells = children_of_type(&s.0, ShapeType::Shell);
        let refs: Vec<usize> = shells
            .iter()
            .map(|sh| self.emit_shell(&Shell(sh.clone())))
            .collect();
        let outer = refs[0];
        self.w.emit(format!("MANIFOLD_SOLID_BREP('',#{outer})"))
    }

    /// Emit the shape (or, for a compound, each child) and return the entity
    /// ids to place in the representation's item list.
    pub(super) fn emit_top(&mut self, shape: &TopoShape) -> Vec<usize> {
        match shape.shape_type() {
            ShapeType::Compound => {
                let kids: Vec<TopoShape> = shape
                    .tshape
                    .read()
                    .unwrap()
                    .children
                    .clone();
                let mut out = Vec::new();
                for k in kids {
                    out.extend(self.emit_top(&k));
                }
                out
            }
            ShapeType::Solid => vec![self.emit_solid(&Solid(shape.clone()))],
            ShapeType::Shell => vec![self.emit_shell(&Shell(shape.clone()))],
            ShapeType::Face => vec![self.emit_face(&Face(shape.clone()))],
            ShapeType::Wire => vec![self.emit_wire(&Wire(shape.clone()))],
            ShapeType::Edge => vec![self.emit_edge(&Edge(shape.clone()))],
            ShapeType::Vertex => vec![self.emit_vertex(&Vertex(shape.clone()))],
            _ => Vec::new(),
        }
    }

    pub(super) fn write_shape_named(&mut self, name: &str, shape: &TopoShape) -> Option<usize> {
        let items = self.emit_top(shape);
        if items.is_empty() {
            return None;
        }
        let n = esc_str(name);
        let rep = self.w.emit(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('{n}',({}),#{})",
            join_refs(&items),
            self.geom_ctx
        ));
        // Product scaffolding, matching `STEPControl_Writer`'s graph (compare an
        // OCCT-written file): `PRODUCT`'s context field is a set of entity refs
        // (`(#n)`), `PRODUCT_DEFINITION_SHAPE` characterizes the
        // `PRODUCT_DEFINITION` (not the `PRODUCT`), and the representation is
        // bound by `SHAPE_DEFINITION_REPRESENTATION`
        // (`STEPConstruct_Styles` / `STEPControl_Writer`; schema
        // `shape_definition_representation.definition : characterized_definition`).
        let product = self
            .w
            .emit(format!("PRODUCT('{n}','{n}','',(#{}))", self.prod_ctx));
        let formation = self
            .w
            .emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        let def = self.w.emit(format!(
            "PRODUCT_DEFINITION('','','',#{formation},#{})",
            self.def_ctx
        ));
        let pds = self.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{def})"));
        self.w
            .emit(format!("SHAPE_DEFINITION_REPRESENTATION(#{pds},#{rep})"));
        Some(rep)
    }
}

/// Direct children of `s` whose type is `t`.
pub(super) fn children_of_type(s: &TopoShape, t: ShapeType) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == t)
        .cloned()
        .collect()
}

/// Extract (axis2, radius) for a cylinder surface by sampling two rings.
pub(super) fn cylinder_params(s: &dyn Surface) -> Option<(GpAx2, f64)> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let vm = 0.5 * (v0 + v1);
    let (c0, r) = ring_center_radius(s, u0, u1, vm)?;
    let axis = ring_axis(s, u0, u1, vm).unwrap_or(dir_z());
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&c0, &s.d0(u0, vm))).unwrap_or(dir_x());
    let ax2 = GpAx2::new(c0, axis, xdir).ok()?;
    Some((ax2, r))
}

/// Extract (axis2, radius, semi_angle) for a cone surface.
pub(super) fn cone_params(s: &dyn Surface) -> Option<(GpAx2, f64, f64)> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let vm = 0.5 * (v0 + v1);
    let v2 = vm + 0.25 * (v1 - v0);
    let (c0, r0) = ring_center_radius(s, u0, u1, vm)?;
    let (c1, r1) = ring_center_radius(s, u0, u1, v2)?;
    let dc = c0.distance(&c1);
    if dc < 1e-12 {
        return None;
    }
    let mut axis = GpDir::from_vec(&GpVec::from_pnts(&c0, &c1)).unwrap_or(dir_z());
    let mut dr = r1 - r0;
    if dr < 0.0 {
        axis = axis.reversed();
        dr = -dr;
    }
    let semi_angle = (dr / dc).atan();
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&c0, &s.d0(u0, vm))).unwrap_or(dir_x());
    let ax2 = GpAx2::new(c0, axis, xdir).ok()?;
    Some((ax2, r0, semi_angle))
}

/// Extract (axis2, major_radius, minor_radius) for a torus surface.
pub(super) fn torus_params(s: &dyn Surface) -> Option<(GpAx2, f64, f64)> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let v_half = v0 + 0.5 * (v1 - v0);
    let (c0, r0) = ring_center_radius(s, u0, u1, v0)?;
    let (_, r1) = ring_center_radius(s, u0, u1, v_half)?;
    let axis = ring_axis(s, u0, u1, v0).unwrap_or(dir_z());
    let major = (r0 + r1) / 2.0;
    let minor = (r0 - r1).abs() / 2.0;
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&c0, &s.d0(u0, v0))).unwrap_or(dir_x());
    let ax2 = GpAx2::new(c0, axis, xdir).ok()?;
    Some((ax2, major, minor))
}

/// Serialize the whole model to a STEP physical file.
pub fn write_step(model: &BRepModel) -> String {
    let mut ctx = WriteCtx::new();
    for ms in &model.shapes {
        ctx.write_shape_named(&ms.name, &ms.shape);
    }
    ctx.finish()
}

/// Serialize the model to a STEP file on disk.
pub fn write_step_file(path: &str, model: &BRepModel) -> std::io::Result<()> {
    std::fs::write(path, write_step(model))
}

/// Serialize a single shape (wrapped in a one-shape model).
pub fn write_shape_step(shape: &TopoShape) -> String {
    let mut model = BRepModel::new();
    model.add("Shape", shape.clone());
    write_step(&model)
}

/// Serialize a single shape, emitting B-spline geometry for non-analytic
/// curves and surfaces (full `STEPControl_Writer` spline coverage).
///
/// Analytic geometry — planes/cylinders/cones/spheres/tori and
/// lines/circles/ellipses/parabolas — is written exactly. Anything else is
/// approximated by an interpolating B-spline and written as
/// `B_SPLINE_CURVE_WITH_KNOTS` / `B_SPLINE_SURFACE_WITH_KNOTS`, which the
/// reader reconstructs as a valid B-spline face/edge.
pub fn write_step_with_splines(shape: &TopoShape) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    ctx.write_shape_named("Shape", shape)
        .ok_or_else(|| "write_step_with_splines: nothing to write".to_string())?;
    Ok(ctx.finish())
}

/// Serialize a single shape, overriding the `PRODUCT` name attribute.
pub fn write_step_with_name(shape: &TopoShape, name: &str) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_name: nothing to write".to_string())?;
    Ok(ctx.finish())
}

/// Serialize a single shape with a minimal STEP style block.
///
/// A `COLOUR_RGB` and the `SURFACE_STYLE_FILL_AREA` / `SURFACE_STYLE_USAGE`
/// style chain are attached to the shape's representation through a
/// `STYLED_ITEM` record. The style is decorative: STEP readers that do not
/// interpret presentation styles simply skip these records.
pub fn write_step_with_color(shape: &TopoShape, color: (f64, f64, f64)) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    let rep = ctx
        .write_shape_named("Shape", shape)
        .ok_or_else(|| "write_step_with_color: nothing to write".to_string())?;
    let colour = ctx.w.emit(format!(
        "COLOUR_RGB('',{},{},{})",
        step_real(color.0),
        step_real(color.1),
        step_real(color.2)
    ));
    let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
    let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
    ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    Ok(ctx.finish())
}

/// SI unit definitions for [`write_step_with_units`].
///
/// The length unit is given as the number of metres per STEP length unit
/// (1.0 = metre, 0.001 = millimetre); the angle unit as radians per STEP
/// angle unit (1.0 = radian, `π/180` ≈ 0.01745 = degree). The emitted
/// `SI_UNIT` records attach these to a `GEOMETRIC_REPRESENTATION_CONTEXT`.
#[derive(Debug, Clone)]
pub struct StepUnits {
    pub length_unit_m: f64,
    pub angle_unit_rad: f64,
}

impl Default for StepUnits {
    fn default() -> Self {
        Self {
            length_unit_m: 1.0,
            angle_unit_rad: 1.0,
        }
    }
}

/// Emit the SI unit block (`DIMENSIONAL_EXPONENTS` + `SI_UNIT` records).
///
/// A dimensional-exponent vector of a length measure (metre^1) and one of a
/// plane angle (radian, dimensionless in the length slots) are shared by the
/// two `SI_UNIT` records. Readers that do not interpret units skip them.
pub(super) fn emit_si_units(w: &mut StepWriter, units: &StepUnits) {
    let _ = units;
    let len_dim = w.emit("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)".into());
    let ang_dim = w.emit("DIMENSIONAL_EXPONENTS(0.,0.,0.,1.,0.,0.,0.)".into());
    w.emit(format!("SI_UNIT(#{len_dim},*,.METRE.)"));
    w.emit(format!("SI_UNIT(#{ang_dim},*,.RADIAN.)"));
}

/// Serialize a single shape with an explicit PRODUCT name and SI unit block.
///
/// The unit records (`DIMENSIONAL_EXPONENTS` + `SI_UNIT`) are emitted after
/// the shape representation, mirroring the unit scaffolding a full
/// `STEPControl_Writer` places in the DATA section. The length/angle scale
/// factors from `units` are recorded so a downstream reader can convert to
/// metres; readers that do not interpret units simply skip the records.
pub fn write_step_with_units(
    shape: &TopoShape,
    name: &str,
    units: &StepUnits,
) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_units: nothing to write".to_string())?;
    emit_si_units(&mut ctx.w, units);
    Ok(ctx.finish())
}

/// Consolidated options for [`write_step_with_options`].
///
/// Combines the PRODUCT name, an optional colour style, SI units and the
/// spline-capable writer into one call — the port of `STEPControl_Writer`'s
/// per-shape configuration (name attribute, presentation colour and unit
/// system).
#[derive(Debug, Clone)]
pub struct StepWriteOptions {
    pub name: String,
    pub color: Option<(f64, f64, f64)>,
    pub units: StepUnits,
    pub splines: bool,
}

impl Default for StepWriteOptions {
    fn default() -> Self {
        Self {
            name: "Shape".into(),
            color: None,
            units: StepUnits::default(),
            splines: false,
        }
    }
}

/// Serialize a single shape with the full [`StepWriteOptions`] configuration.
///
/// This is the one-stop entry point: it writes the shape representation and
/// product scaffolding, then decorates it with a colour style (when
/// `color` is set), an SI unit block, and (when `splines` is set) B-spline
/// output for non-analytic geometry.
pub fn write_step_with_options(shape: &TopoShape, opts: &StepWriteOptions) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = opts.splines;
    let rep = ctx
        .write_shape_named(&opts.name, shape)
        .ok_or_else(|| "write_step_with_options: nothing to write".to_string())?;
    if let Some((r, g, b)) = opts.color {
        let colour = ctx.w.emit(format!(
            "COLOUR_RGB('',{},{},{})",
            step_real(r),
            step_real(g),
            step_real(b)
        ));
        let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
        let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
        ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    }
    emit_si_units(&mut ctx.w, &opts.units);
    Ok(ctx.finish())
}

/// Serialize a single shape with an explicit name and a `COLOUR_RGB` style.
///
/// Combines [`write_step_with_name`] and [`write_step_with_color`]: the
/// product carries `name` and the shape representation is decorated with a
/// `SURFACE_STYLE_FILL_AREA` colour.
pub fn write_step_with_name_and_color(
    shape: &TopoShape,
    name: &str,
    color: (f64, f64, f64),
) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    let rep = ctx
        .write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_name_and_color: nothing to write".to_string())?;
    let colour = ctx.w.emit(format!(
        "COLOUR_RGB('',{},{},{})",
        step_real(color.0),
        step_real(color.1),
        step_real(color.2)
    ));
    let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
    let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
    ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    Ok(ctx.finish())
}

/// Write a shape to a STEP file on disk, using the spline-capable writer.
pub fn write_step_file_with_splines(path: &str, shape: &TopoShape) -> std::io::Result<()> {
    let content = write_step_with_splines(shape).map_err(|e| std::io::Error::other(e))?;
    std::fs::write(path, content)
}

/// Write an assembly to a STEP file on disk.
pub fn write_step_assembly_file(path: &str, a: &StepAssembly) -> std::io::Result<()> {
    let content = write_step_assembly(a).map_err(|e| std::io::Error::other(e))?;
    std::fs::write(path, content)
}

/// A lightweight assembly description: named parts plus a parent→children tree.
///
/// The top-level `name` is emitted as the assembly's own product definition,
/// so assembly trees may reference it as a parent in `children`.
#[derive(Debug, Clone)]
pub struct StepAssembly {
    pub name: String,
    pub products: Vec<(String, TopoShape)>,
    pub children: Vec<(String, Vec<String>)>,
}

/// Serialize an assembly: a `PRODUCT` / `PRODUCT_DEFINITION` pair per part
/// (with the part's shape representation), plus `NEXT_ASSEMBLY_USAGE_OCCURRENCE`
/// records linking the product definitions along the assembly tree.
///
/// Every name in `children` (both parents and children) must be either the
/// assembly `name` or a key of `products`; unknown names return an error.
pub fn write_step_assembly(a: &StepAssembly) -> Result<String, String> {
    write_assembly_inner(a, false)
}

/// Serialize an assembly with B-spline output for non-analytic part geometry.
///
/// Equivalent to [`write_step_assembly`] with the spline-capable geometry
/// classification enabled, so parts whose faces/edges are B-spline (or
/// otherwise non-analytic) round-trip through `B_SPLINE_SURFACE_WITH_KNOTS` /
/// `B_SPLINE_CURVE_WITH_KNOTS`.
pub fn write_step_assembly_with_splines(a: &StepAssembly) -> Result<String, String> {
    write_assembly_inner(a, true)
}
