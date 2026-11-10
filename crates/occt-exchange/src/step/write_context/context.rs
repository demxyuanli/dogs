use super::super::prelude::*;
use super::super::*;

use super::bspline_convert::{bounded_bspline_of_curve, bounded_bspline_surface, curve_to_bspline};
use super::surface_params::{children_of_type, cone_params, cylinder_params, torus_params};


/// Internal writer state: the entity writer plus identity maps so shared
/// vertices/edges/curves are emitted exactly once.
pub(crate) struct WriteCtx {
    pub(crate) w: StepWriter,
    pub(crate) vertex_ids: HashMap<usize, usize>,
    pub(crate) edge_ids: HashMap<usize, usize>,
    pub(crate) curve_ids: HashMap<usize, usize>,
    pub(crate) geom_ctx: usize,
    pub(crate) prod_ctx: usize,
    pub(crate) def_ctx: usize,
    /// When true, non-analytic curves/surfaces are written as B-splines
    /// (sampled + fitted) instead of the tangent-line / plane fallbacks.
    pub(crate) splines: bool,
}

impl WriteCtx {
    pub(crate) fn new() -> Self {
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

    pub(crate) fn finish(self) -> String {
        self.w.finish()
    }

    pub(crate) fn emit_vertex(&mut self, v: &Vertex) -> usize {
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

    pub(crate) fn emit_edge(&mut self, e: &Edge) -> usize {
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

    pub(crate) fn emit_edge_curve(&mut self, e: &Edge) -> Option<usize> {
        // `TopoDSToStep_MakeStepEdge.cxx:194-196` branches on
        // `BRepAdaptor_Curve(aEdge).Curve()` being null, which is the case for a
        // **degenerate** edge (`BRep_Tool::Curve` has no 3D curve). This port
        // stores a point-curve stand-in for those, so the degeneracy flag is what
        // routes to the "edge without 3d curve; creating..." branch below.
        let degenerate = occt_topo::brep_tool::BRepTool::is_degenerated(e)
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

    pub(crate) fn emit_curve_entity(&mut self, c: &dyn Curve) -> Option<usize> {
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
        // only the remaining non-BSpline/Bezier bases need the recursion here.
        // A B-spline basis is cut first (`cxx:71-76`):
        // `BS = B->Copy(); BS->Segment(T->FirstParameter(), T->LastParameter())`.
        if c.is_geom_trimmed() {
            if let Some((basis, _bf, _bl)) = c.untrimmed_basis() {
                if let (Some(poles), Some(knots), Some(deg)) = (
                    basis.bspline_poles(),
                    basis.bspline_knots(),
                    basis.nurbs_degree(),
                ) {
                    let bs = match basis.bspline_weights() {
                        Some(w) if w.len() == poles.len() => GeomBSplineCurve::rational(
                            poles.to_vec(),
                            w.to_vec(),
                            knots.to_vec(),
                            deg,
                        ),
                        _ => GeomBSplineCurve::new(poles.to_vec(), knots.to_vec(), deg),
                    };
                    if let Ok(mut bs) = bs {
                        // `MakeCurve.cxx:73`: the basis is copied before it is
                        // cut, so the copy keeps `B->IsPeriodic()` and
                        // `Geom_BSplineCurve::Segment` runs its periodic branch
                        // (`Geom_BSplineCurve.cxx:769-779`: `SetOrigin` +
                        // `SetNotPeriodic`). The port's constructors do not
                        // carry the flag, so it is restored on the copy.
                        bs.set_periodic_flag(basis.is_periodic());
                        // `theTolerance` defaults to `Precision::PConfusion()`
                        // (`Geom_BSplineCurve.hxx:322-324`). OCCT lets a
                        // failure escape; the port keeps the previous
                        // behaviour (the trimmed curve's own remapped knots,
                        // written by the arm below) instead of dropping the
                        // edge geometry.
                        if bs
                            .segment(
                                c.first_parameter(),
                                c.last_parameter(),
                                occt_core::precision::PCONFUSION,
                            )
                            .is_ok()
                        {
                            if let Ok(id) = write_bspline_curve(&mut self.w, &bs) {
                                return Some(id);
                            }
                        }
                    }
                } else if basis.bezier_poles().is_none() {
                    return self.emit_curve_entity(basis.as_ref());
                }
            }
        }
        // `MakeCurve.cxx:94-99`: `Geom_BoundedCurve` → `GeomToStep_MakeBoundedCurve`
        // (`GeomToStep_MakeBoundedCurve.cxx:37-80`): a B-spline goes to
        // `B_SPLINE_CURVE_WITH_KNOTS` (+`_AND_RATIONAL_...` when rational) after
        // the periodic arm (`cxx:44-51`) has opened it, and a Bezier is first
        // converted with `GeomConvert::CurveToBSplineCurve` and written the same
        // way.
        if c.bspline_poles().is_some() {
            if let Some(bs) = bounded_bspline_of_curve(c) {
                if let Ok(id) = write_bspline_curve(&mut self.w, &bs) {
                    return Some(id);
                }
            }
        }
        if c.bezier_poles().is_some() {
            // `GeomToStep_MakeCurve.cxx:77-82` trims a Bezier basis with
            // `Geom_BezierCurve::Segment` before converting it; `Segment` is
            // UNPORTED in this port, so a trimmed Bezier keeps the previous
            // behaviour (convert its **basis**, the port's documented
            // divergence) instead of dropping the edge geometry.
            let basis = c.untrimmed_basis().map(|(b, _, _)| b);
            let target: &dyn Curve = basis.as_deref().unwrap_or(c);
            if let Some(bs) = curve_to_bspline(target) {
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

    pub(crate) fn emit_wire(&mut self, w: &Wire) -> usize {
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

    pub(crate) fn emit_face(&mut self, f: &Face) -> usize {
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

    pub(crate) fn emit_surface(&mut self, f: &Face) -> usize {
        let Some(surf) = GeometryRegistry::global().face_surface(&f.0) else {
            return self.emit_plane_fallback();
        };
        // `GeomToStep_MakeSurface.cxx:44-48` transfers a `Geom_BoundedSurface`
        // (`Geom_BSplineSurface` / `Geom_BezierSurface` /
        // `Geom_RectangularTrimmedSurface`) by dynamic type **before** the
        // elementary, swept and offset arms, and
        // `GeomToStep_MakeBoundedSurface.cxx:39-95` writes the surface's own
        // poles and knots. The geometric classifier below cannot identify those
        // classes, so it must not see them first. Only the spline writer takes
        // this arm: the plain writer keeps its previous entity choice.
        if self.splines {
            if let Some(bs) = bounded_bspline_surface(surf.as_ref()) {
                if let Ok(id) = write_bspline_surface(&mut self.w, &bs) {
                    return id;
                }
            }
        }
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
                // `GeomToStep_MakeSurface.cxx:78-84`: an unrecognised surface
                // leaves `done = false`; the caller then builds the
                // `ADVANCED_FACE` with the surface attribute unset. The port
                // writes its plane placeholder instead, and never resamples the
                // surface (a fit is not in OCCT's writer).
                self.emit_plane_fallback()
            }
        }
    }

    pub(crate) fn emit_ax2(&mut self, ax2: &GpAx2) -> usize {
        self.w
            .add_axis2_placement_3d(&ax2.location(), &ax2.direction(), &ax2.x_direction())
    }

    pub(crate) fn emit_plane_fallback(&mut self) -> usize {
        let ax = self.w.add_axis2_placement_3d(&GpPnt::zero(), &dir_z(), &dir_x());
        self.w.emit(format!("PLANE('',#{ax})"))
    }

    pub(crate) fn emit_shell(&mut self, sh: &Shell) -> usize {
        let faces = children_of_type(&sh.0, ShapeType::Face);
        let refs: Vec<usize> = faces
            .iter()
            .map(|f| self.emit_face(&Face(f.clone())))
            .collect();
        self.w
            .emit(format!("CLOSED_SHELL('',({}))", join_refs(&refs)))
    }

    pub(crate) fn emit_solid(&mut self, s: &Solid) -> usize {
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
    pub(crate) fn emit_top(&mut self, shape: &TopoShape) -> Vec<usize> {
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

    pub(crate) fn write_shape_named(&mut self, name: &str, shape: &TopoShape) -> Option<usize> {
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
            "PRODUCT_DEFINITION('','',#{formation},#{})",
            self.def_ctx
        ));
        let pds = self.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{def})"));
        self.w
            .emit(format!("SHAPE_DEFINITION_REPRESENTATION(#{pds},#{rep})"));
        Some(rep)
    }
}
