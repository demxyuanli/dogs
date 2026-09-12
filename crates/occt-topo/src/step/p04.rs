use super::prelude::*;
use super::*;

impl<'a> Resolver<'a> {

    pub(super) fn new(records: &'a HashMap<usize, Record>) -> Self {
        Self {
            records,
            b: TopoBuilder::new(),
            shape_cache: RefCell::new(HashMap::new()),
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
        }
    }

    pub(super) fn record(&self, id: usize) -> Result<&'a Record, String> {
        self.records
            .get(&id)
            .ok_or_else(|| format!("reference to undefined entity #{id}"))
    }

    pub(super) fn warn(&self, msg: String) {
        self.warnings.borrow_mut().push(msg);
    }

    pub(super) fn resolve_shape(&self, id: usize) -> Result<TopoShape, String> {
        if let Some(s) = self.shape_cache.borrow().get(&id) {
            return Ok(s.clone());
        }
        let rec = self.record(id)?;
        if !self.resolving.borrow_mut().insert(id) {
            return Err(format!("cyclic entity reference #{id}"));
        }
        let shape = match rec.type_name.as_str() {
            "VERTEX_POINT" => self.resolve_vertex(rec),
            "EDGE_CURVE" => self.resolve_edge(rec),
            "ORIENTED_EDGE" => self.resolve_oriented_edge(rec),
            "EDGE_LOOP" => self.resolve_loop(rec),
            "VERTEX_LOOP" => self.resolve_vertex_loop(rec),
            "FACE_OUTER_BOUND" => self.resolve_outer_bound(rec),
            "FACE_BOUND" => self.resolve_outer_bound(rec),
            "ADVANCED_FACE" => self.resolve_face(rec),
            "CLOSED_SHELL" => self.resolve_shell(rec),
            "MANIFOLD_SOLID_BREP" => self.resolve_solid(rec),
            other => {
                self.warn(format!("unsupported topological entity {other} (#{id})"));
                Err(format!("unsupported entity {other} (#{id})"))
            }
        };
        self.resolving.borrow_mut().remove(&id);
        match shape {
            Ok(s) => {
                self.shape_cache.borrow_mut().insert(id, s.clone());
                Ok(s)
            }
            Err(e) => Err(e),
        }
    }

    pub(super) fn resolve_vertex(&self, rec: &'a Record) -> Result<TopoShape, String> {
        if rec.args.len() < 2 {
            return Err("VERTEX_POINT: bad args".into());
        }
        let pid = parse_ref(&rec.args[1]).ok_or("VERTEX_POINT: bad point ref")?;
        let p = self.resolve_point(pid)?;
        Ok(self.b.make_vertex(p, 0.0).0)
    }

    pub(super) fn resolve_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let start_ref = parse_ref(&rec.args[1]).ok_or("EDGE_CURVE: bad start ref")?;
        let end_ref = parse_ref(&rec.args[2]).ok_or("EDGE_CURVE: bad end ref")?;
        let curve_ref = parse_ref(&rec.args[3]).ok_or("EDGE_CURVE: bad curve ref")?;
        let v1 = self.resolve_shape(start_ref)?;
        let v2 = self.resolve_shape(end_ref)?;
        if !v1.is_vertex() || !v2.is_vertex() {
            return Err("EDGE_CURVE: endpoints are not vertices".into());
        }
        let curve = self.resolve_curve(curve_ref)?;
        let p1 = GeometryRegistry::global().vertex_point(&v1);
        let p2 = GeometryRegistry::global().vertex_point(&v2);
        let (curve, first, last) = edge_from_curve3d(curve, &p1, &p2);
        let mut e = self.b.make_edge(curve, first, last);
        // `BRep_Builder` / `BRepLib_MakeEdge`: first vertex FORWARD, last REVERSED
        // (`TopoDS_Builder::Add`). Needed so `TopExp::LastVertex` exists on a
        // closed EDGE_CURVE (same TVertex stored twice) and
        // `ShapeAnalysis_Edge::FirstVertex` works on REVERSED seam uses.
        self.b.add_edge_vertices(&mut e, &Vertex(v1), &Vertex(v2));
        // `MakeFromCurve3D` (`cxx:478-479`) `UpdateVertex(1.000001 * dist)` is
        // skipped until `ShapeAnalysis_Curve::Project` is Extrema-accurate:
        // sampler residuals were written into vertex tolerance and
        // `MaxFaceTolerance` / edge discret over-tessellated Shape.step.
        // Remember the edge's curve entity for the face-level pcurve association
        // (a SURFACE_CURVE's pcurve is matched to the face's surface by ref).
        self.edge_curve_ref
            .borrow_mut()
            .insert(Arc::as_ptr(&e.0.tshape) as usize, curve_ref);
        Ok(e.0)
    }

    pub(super) fn resolve_oriented_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let edge_ref = parse_ref(&rec.args[3]).ok_or("ORIENTED_EDGE: bad edge ref")?;
        let mut s = self.resolve_shape(edge_ref)?;
        // ORIENTED_EDGE(name, *, *, edge, orientation): a .F. reverses the edge
        // in the wire. Without this, two ORIENTED_EDGEs referencing the same
        // EDGE_CURVE (e.g. a cylinder side wall's two seam generatrices) both
        // return the same forward edge, and the wire loses one — the surface
        // never closes.
        if let Some(o) = rec.args.get(4).map(|s| s.trim().to_string()) {
            if o == ".F." {
                s.set_orientation(Orientation::Reversed);
            }
        }
        Ok(s)
    }

    pub(super) fn resolve_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let items = parse_ref_list(&rec.args[1]);
        let mut edges = Vec::with_capacity(items.len());
        for &it in &items {
            let s = self.resolve_shape(it)?;
            if !s.is_edge() {
                return Err(format!("#{it}: expected EDGE in EDGE_LOOP"));
            }
            edges.push(Edge(s));
        }
        let mut wire = self.b.make_wire(&edges).0;
        // The EDGE_LOOP lists its ORIENTED_EDGEs in any cyclic order; the file
        // is not guaranteed to place consecutive edges next to each other (a
        // CAD writer may list them out of sequence). Reorder the wire's edges
        // into a connected chain by following shared vertices, as the OCCT STEP
        // reader does when assembling the loop.
        crate::algo_tools::AlgoTools::orient_edges_on_wire(&mut wire);
        Ok(wire)
    }

    /// A `VERTEX_LOOP(name, vertex)` is the boundary of a degenerate face — a
    /// loop reduced to a single vertex (sphere pole, cone apex). It maps to a
    /// wire containing no edges (the face still references it via FACE_BOUND).
    pub(super) fn resolve_vertex_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let _ = parse_ref(&rec.args[1]).ok_or("VERTEX_LOOP: bad vertex ref")?;
        Ok(self.b.make_wire(&[]).0)
    }

    pub(super) fn resolve_outer_bound(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let loop_ref = parse_ref(&rec.args[1]).ok_or("FACE_OUTER_BOUND: bad loop ref")?;
        let mut s = self.resolve_shape(loop_ref)?;
        if !s.is_wire() {
            return Err("FACE_OUTER_BOUND: loop is not a wire".into());
        }
        // `StepShape_FaceBound::Orientation`: `.F.` reverses the wire
        // (`StepToTopoDS_TranslateFace.cxx` PolyLoop / EdgeLoop bind).
        if !parse_logical(rec.args.get(2).map(String::as_str), true) {
            s.set_orientation(Orientation::Reversed);
        }
        Ok(s)
    }

    pub(super) fn resolve_face(&self, rec: &'a Record) -> Result<TopoShape, String> {
        // Two argument layouts occur in the wild:
        //  * STEP-214 (ISO standard, written by OCCT/FreeCAD): surface is the
        //    third argument — `ADVANCED_FACE(name, bounds, surface, same_sense)`.
        //  * this port's own writer (step.rs write path): surface is the second
        //    argument — `ADVANCED_FACE('', #surface, (bounds), .T.)`.
        // Detect by checking which argument holds a surface reference.
        let surf_ref = if let Some(r) = parse_ref(&rec.args[2]) {
            r
        } else {
            parse_ref(&rec.args[1]).ok_or("ADVANCED_FACE: bad surface ref")?
        };
        let surface = self.resolve_surface(surf_ref)?;
        // Bounds are the other of the two leading arguments.
        let bounds = if rec.args.get(2).map(|s| s.starts_with('(')).unwrap_or(false) {
            parse_ref_list(&rec.args[2])
        } else {
            parse_ref_list(&rec.args[1])
        };
        // `StepToTopoDS_TranslateFace.cxx:572-573` / `716-752`.
        let same_sense = parse_logical(rec.args.last().map(String::as_str), true);
        let a_same_sense = if self.step_surface_is_reversed(surf_ref) {
            !same_sense
        } else {
            same_sense
        };
        // `StepToTopoDS_TranslateFace.cxx:608-634` — a lone VertexLoop on a
        // sphere / BSpline / revolution is the whole closed surface; add
        // `BRepLib_MakeFace` natural bounds and skip the vertex loop itself.
        if bounds.len() == 1
            && self.bound_loop_is_vertex_loop(bounds[0])
            && (surface.gp_sphere().is_some()
                || surface.is_bspline_surface()
                || surface.is_surface_of_revolution())
        {
            let mut face = crate::brep_lib_make_face::make_face_from_surface(
                surface,
                occt_core::precision::Precision::CONFUSION,
            );
            face.0.set_orientation(if same_sense {
                Orientation::Forward
            } else {
                Orientation::Reversed
            });
            return Ok(face.0);
        }
        let mut wires = Vec::with_capacity(bounds.len());
        for &b in &bounds {
            let mut s = self.resolve_shape(b)?;
            if !s.is_wire() {
                return Err(format!("#{b}: expected wire in face bounds"));
            }
            // Bound orientation is already on the wire from `resolve_outer_bound`.
            // A reversed `Face_Surface` (or negative-major torus) flips it again
            // so the stored wire matches CAS.CADE (`cxx:714-723`).
            if !a_same_sense {
                s.reverse();
            }
            wires.push(Wire(s));
        }
        let mut face = self.b.make_face(surface, &wires);
        face.0.set_orientation(if same_sense {
            Orientation::Forward
        } else {
            Orientation::Reversed
        });
        let face_key = GeometryRegistry::shape_key(&face.0);
        // Associate each wire edge's SURFACE_CURVE pcurve with this face's
        // surface (`BRep_Builder::UpdateEdge(edge, pcurve, face, tol)`).
        let face_ori = face.0.orientation();
        for w in &wires {
            for e in edges_of_wire(w) {
                self.associate_edge_pcurve(&e, surf_ref, face_key, face_ori)?;
            }
            crate::shhealing::check_pcurves_and_shift(w, &face);
        }
        Ok(face.0)
    }

    /// True when `bound_id` is a `FACE_BOUND` / `FACE_OUTER_BOUND` whose loop
    /// is a `VERTEX_LOOP` (`StepToTopoDS_TranslateFace.cxx:608`).
    fn bound_loop_is_vertex_loop(&self, bound_id: usize) -> bool {
        let Ok(rec) = self.record(bound_id) else {
            return false;
        };
        match rec.type_name.as_str() {
            "FACE_BOUND" | "FACE_OUTER_BOUND" => parse_ref(&rec.args[1])
                .and_then(|id| self.record(id).ok())
                .map(|r| r.type_name == "VERTEX_LOOP")
                .unwrap_or(false),
            "VERTEX_LOOP" => true,
            _ => false,
        }
    }

    /// `StepToTopoDS_TranslateFace.cxx:479-491` — SolidWorks torus with a
    /// negative major radius is treated as a reversed `Face_Surface`.
    pub(super) fn step_surface_is_reversed(&self, surf_id: usize) -> bool {
        let Ok(rec) = self.record(surf_id) else {
            return false;
        };
        match rec.type_name.as_str() {
            "RECTANGULAR_TRIMMED_SURFACE" => rec
                .args
                .get(1)
                .and_then(|s| parse_ref(s))
                .map(|basis| self.step_surface_is_reversed(basis))
                .unwrap_or(false),
            "TOROIDAL_SURFACE" => rec
                .args
                .get(2)
                .and_then(|s| parse_f64(s).ok())
                .map(|maj| maj < 0.0)
                .unwrap_or(false),
            _ => false,
        }
    }

    /// Match a wire edge's SURFACE_CURVE pcurve to `surf_ref` and attach its 2D
    /// curve to the edge for `face_key`. Source:
    /// `StepToTopoDS_GeometricTool::PCurve` + `StepToTopoDS_TranslateEdge::MakePCurve`.
    pub(super) fn associate_edge_pcurve(
        &self,
        edge: &Edge,
        surf_ref: usize,
        face_key: usize,
        face_ori: Orientation,
    ) -> Result<(), String> {
        let key = Arc::as_ptr(&edge.0.tshape) as usize;
        let Some(&curve_ref) = self.edge_curve_ref.borrow().get(&key) else {
            return Ok(());
        };
        let Some(pcurves) = self.surface_curve_pcurves.borrow().get(&curve_ref).cloned() else {
            return Ok(());
        };
        // Each associated_geometry entry is a pcurve_or_surface: a PCURVE (with
        // a 2D curve) or a SURFACE (intersection curve). Only PCURVEs carry the
        // pcurve; a SURFACE entry fails to resolve and is skipped. The raw 2D
        // curve is stored with its own parameterization; the meshing maps the
        // edge's 3D parameter onto it (`BRepMesh_EdgeParameterProvider`).
        let mut matched: Vec<Arc<dyn Curve2d>> = Vec::new();
        for pc in pcurves {
            if let Ok((basis_surf, c2d)) = self.resolve_pcurve(pc) {
                if basis_surf == surf_ref {
                    matched.push(c2d);
                }
            }
        }
        // A seam edge carries two pcurves on the same face (one per side of the
        // seam, e.g. the cone's `u = 0` and `u = 2π`). Order them forward-then-
        // reversed (`ShapeAnalysis_Curve::SelectForwardSeam`) so the wire assembly
        // can hand each traversal its own side.
        if matched.len() == 2 {
            // `StepToTopoDS_TranslateEdgeLoop.cxx:721-734`: SelectForwardSeam,
            // then invert 1↔2 when Compose(Edge, Wire, Face) is not FORWARD.
            // `edges_of_wire` already composes the parent wire onto the edge.
            let mut fwd =
                crate::pcurve_full::select_forward_seam(matched[0].as_ref(), matched[1].as_ref());
            if fwd != 0 {
                let cumul = Orientation::compose(face_ori, edge.0.orientation());
                if cumul != Orientation::Forward {
                    fwd = 3 - fwd;
                }
                if fwd == 2 {
                    matched.swap(0, 1);
                }
            }
        }
        if !matched.is_empty() {
            GeometryRegistry::global().set_edge_pcurves(&edge.0, face_key, matched);
        }
        Ok(())
    }

    pub(super) fn resolve_shell(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let faces = parse_ref_list(&rec.args[1]);
        let mut fs = Vec::with_capacity(faces.len());
        for &it in &faces {
            let s = self.resolve_shape(it)?;
            if !s.is_face() {
                return Err(format!("#{it}: expected FACE in CLOSED_SHELL"));
            }
            fs.push(Face(s));
        }
        Ok(self.b.make_shell(&fs).0)
    }

    pub(super) fn resolve_solid(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let outer = parse_ref(&rec.args[1]).ok_or("MANIFOLD_SOLID_BREP: bad shell ref")?;
        let s = self.resolve_shape(outer)?;
        if !s.is_shell() {
            return Err("MANIFOLD_SOLID_BREP: outer is not a shell".into());
        }
        Ok(self.b.make_solid(&[Shell(s)]).0)
    }

    pub(super) fn resolve_point(&self, id: usize) -> Result<GpPnt, String> {
        if let Some(p) = self.point_cache.borrow().get(&id) {
            return Ok(*p);
        }
        let rec = self.record(id)?;
        if rec.type_name != "CARTESIAN_POINT" {
            self.warn(format!("expected CARTESIAN_POINT, got {} (#{id})", rec.type_name));
            return Err(format!("expected CARTESIAN_POINT at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("CARTESIAN_POINT #{id}: {e}"))?;
        let p = GpPnt::from_xyz(&v);
        self.point_cache.borrow_mut().insert(id, p);
        Ok(p)
    }

    pub(super) fn resolve_direction(&self, id: usize) -> Result<GpDir, String> {
        if let Some(d) = self.dir_cache.borrow().get(&id) {
            return Ok(*d);
        }
        let rec = self.record(id)?;
        if rec.type_name != "DIRECTION" {
            self.warn(format!("expected DIRECTION, got {} (#{id})", rec.type_name));
            return Err(format!("expected DIRECTION at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        let d = GpDir::new(v.x, v.y, v.z).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        self.dir_cache.borrow_mut().insert(id, d);
        Ok(d)
    }

    pub(super) fn resolve_vector(&self, id: usize) -> Result<GpVec, String> {
        let rec = self.record(id)?;
        if rec.type_name != "VECTOR" {
            self.warn(format!("expected VECTOR, got {} (#{id})", rec.type_name));
            return Err(format!("expected VECTOR at #{id}"));
        }
        let dir_ref = parse_ref(&rec.args[1]).ok_or("VECTOR: bad direction ref")?;
        let mag = parse_f64(&rec.args[2])?;
        let dir = self.resolve_direction(dir_ref)?;
        Ok(GpVec::from_xyz(&dir.xyz().multiplied(mag)))
    }

    /// `AXIS1_PLACEMENT(name, location, axis_direction)` → `GpAx1`.
    pub(super) fn resolve_axis1(&self, id: usize) -> Result<GpAx1, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS1_PLACEMENT" {
            self.warn(format!("expected AXIS1_PLACEMENT, got {} (#{id})", rec.type_name));
            return Err(format!("expected AXIS1_PLACEMENT at #{id}"));
        }
        let loc_ref = parse_ref(&rec.args[1]).ok_or("AXIS1: bad location ref")?;
        let dir_ref = parse_ref(&rec.args[2]).ok_or("AXIS1: bad direction ref")?;
        let loc = self.resolve_point(loc_ref)?;
        let dir = self.resolve_direction(dir_ref)?;
        Ok(GpAx1::new(loc, dir))
    }

    pub(super) fn resolve_axis2(&self, id: usize) -> Result<GpAx2, String> {
        if let Some(a) = self.axis_cache.borrow().get(&id) {
            return Ok(*a);
        }
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_3D" {
            self.warn(format!(
                "expected AXIS2_PLACEMENT_3D, got {} (#{id})",
                rec.type_name
            ));
            return Err(format!("expected AXIS2_PLACEMENT_3D at #{id}"));
        }
        let loc = self
            .resolve_point(parse_ref(&rec.args[1]).ok_or("AXIS2: bad location ref")?)?;
        let axis = match parse_ref(&rec.args[2]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_z(),
        };
        let refd = match parse_ref(&rec.args[3]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_x(),
        };
        let ax2 = GpAx2::new(loc, axis, refd).map_err(|e| format!("AXIS2_PLACEMENT_3D: {e}"))?;
        self.axis_cache.borrow_mut().insert(id, ax2);
        Ok(ax2)
    }

    pub(super) fn resolve_curve(&self, id: usize) -> Result<Arc<dyn Curve>, String> {
        if let Some(c) = self.curve_cache.borrow().get(&id) {
            return Ok(c.clone());
        }
        let rec = self.record(id)?;
        let curve: Arc<dyn Curve> = match rec.type_name.as_str() {
            "LINE" => {
                let pnt = parse_ref(&rec.args[1]).ok_or("LINE: bad point ref")?;
                let vec = parse_ref(&rec.args[2]).ok_or("LINE: bad vector ref")?;
                let p = self.resolve_point(pnt)?;
                let v = self.resolve_vector(vec)?;
                let d = GpDir::from_vec(&v).map_err(|e| format!("LINE: {e}"))?;
                Arc::new(GeomLine::new(GpLin::from_pnt_dir(p, d)))
            }
            "CIRCLE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CIRCLE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomCircle::new(GpCirc::new(self.resolve_axis2(ax)?, r)))
            }
            "ELLIPSE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("ELLIPSE: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                Arc::new(GeomEllipse::new(GpElips::new(
                    self.resolve_axis2(ax)?,
                    maj,
                    min,
                )))
            }
            "HYPERBOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("HYPERBOLA: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                Arc::new(GeomHyperbola::new(GpHypr::new(
                    self.resolve_axis2(ax)?,
                    maj,
                    min,
                )))
            }
            "PARABOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("PARABOLA: bad axis ref")?;
                let f = parse_f64(&rec.args[2])?;
                Arc::new(GeomParabola::new(GpParab::new(self.resolve_axis2(ax)?, f)))
            }
            "SURFACE_CURVE" | "SEAM_CURVE" => {
                // SURFACE_CURVE/SEAM_CURVE(name, curve_3d, pcurves, master_rep):
                // the 3D curve is the second argument; the pcurve list
                // (associated_geometry) is carried for the face-level pcurve
                // association. SEAM_CURVE is the seam of a closed surface.
                let c3d = parse_ref(&rec.args[1]).ok_or("SURFACE_CURVE: bad 3D curve ref")?;
                self.surface_curve_pcurves
                    .borrow_mut()
                    .insert(id, parse_ref_list(&rec.args[2]));
                self.resolve_curve(c3d)?
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                // Layout (10 args): name, degree, control_points, weights|SELF,
                // curve_form, closed, self_intersect, knots, multiplicities, knot_spec.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let weights_arg = rec.args.get(3).map(|s| s.trim().to_string());
                // B_SPLINE_CURVE_WITH_KNOTS layout:
                // (name, degree, control_points, curve_form, closed,
                //  self_intersect, knot_multiplicities, knots, knot_spec).
                let knots = expand_knots(
                    &parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()")),
                );
                let curve = match weights_arg.as_deref() {
                    // No weights: non-rational curve (curve_form is UNSPECIFIED
                    // or a non-rational flag like CIRCULAR/LINEAR).
                    None | Some("SELF") | Some(".UNSPECIFIED.") | Some(".CIRCULAR.")
                    | Some(".LINEAR.") => GeomBSplineCurve::new(poles, knots, degree),
                    Some(w) if w.starts_with('.') => {
                        // Another non-rational curve form marker.
                        GeomBSplineCurve::new(poles, knots, degree)
                    }
                    Some(w) => {
                        let weights = parse_real_list(w);
                        if weights.len() != poles.len() {
                            return Err("B_SPLINE_CURVE: weight count mismatch".into());
                        }
                        GeomBSplineCurve::rational(poles, weights, knots, degree)
                    }
                };
                Arc::new(curve.map_err(|e| format!("B_SPLINE_CURVE: {e}"))?)
            }
            "POLYLINE" => {
                // A connected sequence of CARTESIAN_POINTs; represent it as a
                // degree-1 B-spline that passes through every vertex.
                let pts: Vec<GpPnt> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                if pts.len() < 2 {
                    return Err("POLYLINE: need at least 2 points".into());
                }
                let knots = uniform_knots_for(pts.len(), 1);
                Arc::new(
                    GeomBSplineCurve::new(pts, knots, 1)
                        .map_err(|e| format!("POLYLINE: {e}"))?,
                )
            }
            "B_SPLINE_CURVE" => {
                // Plain B-spline (no explicit knots): the knot vector is the
                // clamped uniform one implied by the pole count and degree.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let knots = uniform_knots_for(poles.len(), degree);
                Arc::new(
                    GeomBSplineCurve::new(poles, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE: {e}"))?,
                )
            }
            "TRIMMED_CURVE" => {
                // Layout: name, basis_curve, trim_1, trim_2, sense_agreement,
                // master_representation. The bounds are the 4th/5th attributes.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("TRIMMED_CURVE: bad basis ref")?;
                let basis = self.resolve_curve(basis_ref)?;
                let a = parse_f64(&rec.args[3])?;
                let b = parse_f64(&rec.args[4])?;
                Arc::new(GeomTrimmedCurve::new(basis, a, b))
            }
            "OFFSET_CURVE_3D" => {
                // Layout: name, basis_curve, direction, distance, self_intersect,
                // curve_form.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("OFFSET_CURVE_3D: bad curve ref")?;
                let dir_ref = parse_ref(&rec.args[2]).ok_or("OFFSET_CURVE_3D: bad dir ref")?;
                let offset = parse_f64(&rec.args[3])?;
                let basis = self.resolve_curve(basis_ref)?;
                let dir = self.resolve_direction(dir_ref)?;
                Arc::new(GeomOffsetCurve::new(basis, offset, dir))
            }
            other => {
                self.warn(format!("unsupported curve entity {other} (#{id})"));
                return Err(format!("unsupported curve entity {other} (#{id})"));
            }
        };
        self.curve_cache.borrow_mut().insert(id, curve.clone());
        Ok(curve)
    }

    /// 2D `CARTESIAN_POINT` → `GpPnt2d` (STEP pcurves use two-component points).
    pub(super) fn resolve_point_2d(&self, id: usize) -> Result<GpPnt2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "CARTESIAN_POINT" {
            return Err(format!("expected CARTESIAN_POINT at #{id}"));
        }
        let (x, y) = parse_xy(&rec.args[1]).map_err(|e| format!("CARTESIAN_POINT #{id}: {e}"))?;
        Ok(GpPnt2d::new(x, y))
    }

    /// 2D `DIRECTION` → `GpDir2d`.
    pub(super) fn resolve_direction_2d(&self, id: usize) -> Result<GpDir2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "DIRECTION" {
            return Err(format!("expected DIRECTION at #{id}"));
        }
        let (x, y) = parse_xy(&rec.args[1]).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        GpDir2d::new(x, y).map_err(|e| format!("DIRECTION #{id}: {e}"))
    }

    /// 2D `VECTOR` → `GpVec2d`.
    pub(super) fn resolve_vector_2d(&self, id: usize) -> Result<GpVec2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "VECTOR" {
            return Err(format!("expected VECTOR at #{id}"));
        }
        let dir_ref = parse_ref(&rec.args[1]).ok_or("VECTOR: bad direction ref")?;
        let mag = parse_f64(&rec.args[2])?;
        let dir = self.resolve_direction_2d(dir_ref)?;
        Ok(GpVec2d::new(dir.x * mag, dir.y * mag))
    }

    /// 2D `AXIS2_PLACEMENT_2D` → `GpAx22d` (the ref_direction is X; Y is its
    /// counter-clockwise normal).
    pub(super) fn resolve_axis22d(&self, id: usize) -> Result<GpAx22d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_2D" {
            return Err(format!("expected AXIS2_PLACEMENT_2D at #{id}"));
        }
        let loc = self.resolve_point_2d(parse_ref(&rec.args[1]).ok_or("AXIS2_2D: bad loc ref")?)?;
        let vx = match parse_ref(&rec.args[2]) {
            Some(r) => self.resolve_direction_2d(r)?,
            None => GpDir2d::new(1.0, 0.0).unwrap(),
        };
        let vy = GpDir2d::new(-vx.y, vx.x).map_err(|e| format!("AXIS2_2D: {e}"))?;
        GpAx22d::new(loc, vx, vy).map_err(|e| format!("AXIS2_PLACEMENT_2D: {e}"))
    }

    /// 2D curve → `Arc<dyn Curve2d>` (STEP pcurve geometry). Source:
    /// `StepToGeom::MakeCurve2d`.
    pub(super) fn resolve_curve_2d(&self, id: usize) -> Result<(Arc<dyn Curve2d>, (f64, f64)), String> {
        if let Some(c) = self.curve2d_cache.borrow().get(&id) {
            return Ok(c.clone());
        }
        let rec = self.record(id)?;
        let (curve, range): (Arc<dyn Curve2d>, (f64, f64)) = match rec.type_name.as_str() {
            "LINE" => {
                let pnt = parse_ref(&rec.args[1]).ok_or("LINE: bad point ref")?;
                let vec = parse_ref(&rec.args[2]).ok_or("LINE: bad vector ref")?;
                let p = self.resolve_point_2d(pnt)?;
                let v = self.resolve_vector_2d(vec)?;
                let d = GpDir2d::new(v.x(), v.y()).map_err(|e| format!("LINE: {e}"))?;
                let mag = (v.x() * v.x() + v.y() * v.y()).sqrt();
                (Arc::new(Geom2dLine::new(GpAx2d::new(p, d))), (0.0, mag))
            }
            "CIRCLE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CIRCLE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                let ax22 = self.resolve_axis22d(ax)?;
                let c: Arc<dyn Curve2d> = Arc::new(Geom2dCircle::new(GpCirc2d::new(ax22, r)));
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "ELLIPSE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("ELLIPSE: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                let ax22 = self.resolve_axis22d(ax)?;
                let c: Arc<dyn Curve2d> =
                    Arc::new(Geom2dEllipse::new(GpElips2d::new(ax22, maj, min)));
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = parse_f64(&rec.args[1])? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let knots = expand_knots(
                    &parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()")),
                );
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(
                    Geom2dBSplineCurve::new(xs, ys, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE_2D: {e}"))?,
                );
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "POLYLINE" => {
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                if pts.len() < 2 {
                    return Err("POLYLINE: need at least 2 points".into());
                }
                let knots = uniform_knots_for(pts.len(), 1);
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(
                    Geom2dBSplineCurve::new(xs, ys, knots, 1)
                        .map_err(|e| format!("POLYLINE: {e}"))?,
                );
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "B_SPLINE_CURVE" => {
                let degree = parse_f64(&rec.args[1])? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let knots = uniform_knots_for(pts.len(), degree);
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(
                    Geom2dBSplineCurve::new(xs, ys, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE_2D: {e}"))?,
                );
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "TRIMMED_CURVE" => {
                let basis_ref = parse_ref(&rec.args[1]).ok_or("TRIMMED_CURVE: bad basis ref")?;
                let (basis, _) = self.resolve_curve_2d(basis_ref)?;
                let a = parse_f64(&rec.args[3])?;
                let b = parse_f64(&rec.args[4])?;
                (Arc::new(Geom2dTrimmedCurve::new(basis, a, b)), (a, b))
            }
            other => {
                self.warn(format!("unsupported 2D curve entity {other} (#{id})"));
                return Err(format!("unsupported 2D curve entity {other} (#{id})"));
            }
        };
        self.curve2d_cache.borrow_mut().insert(id, (curve.clone(), range));
        Ok((curve, range))
    }

    /// `PCURVE(name, basis_surface, reference_to_curve)` → the pcurve's basis
    /// surface reference and its 2D curve (the first item of the
    /// DEFINITIONAL_REPRESENTATION). Source: `StepToTopoDS_TranslateEdge::MakePCurve`.
    pub(super) fn resolve_pcurve(&self, id: usize) -> Result<(usize, Arc<dyn Curve2d>), String> {
        let rec = self.record(id)?;
        if rec.type_name != "PCURVE" {
            return Err(format!("expected PCURVE at #{id}"));
        }
        let basis_surf = parse_ref(&rec.args[1]).ok_or("PCURVE: bad basis surface ref")?;
        let dri = parse_ref(&rec.args[2]).ok_or("PCURVE: bad reference_to_curve ref")?;
        let dri_rec = self.record(dri)?;
        if dri_rec.type_name != "DEFINITIONAL_REPRESENTATION" {
            return Err(format!("expected DEFINITIONAL_REPRESENTATION at #{dri}"));
        }
        let items = parse_ref_list(&dri_rec.args[1]);
        let curve_ref = *items.first().ok_or("PCURVE: empty reference_to_curve")?;
        let (c2d, _range) = self.resolve_curve_2d(curve_ref)?;
        // `StepToTopoDS_TranslateEdge::MakePCurve` then
        // `GeomConvert_Units::DegreeToRadian` (`GeomConvert_Units.cxx:173-302`).
        let surf = self.resolve_surface(basis_surf)?;
        let c2d = self.degree_to_radian(c2d, curve_ref, surf.as_ref());
        Ok((basis_surf, c2d))
    }

    /// `GeomConvert_Units::DegreeToRadian`. STEP stores cylinder/cone/sphere/torus
    /// pcurve `U` in the file angle unit and cone `V` as axis length; OCCT's
    /// `ElSLib` cone uses `V` along the generatrix (`Length / cos(semiAngle)`).
    /// Length/angle factors stay 1: this reader does not rescale 3D coordinates
    /// from `SI_UNIT` records (`step/p02.rs` write-side units are skipped on read).
    ///
    /// Line location/direction are scaled then the direction is renormalized
    /// (`GeomConvert_Units.cxx:260-276`), so a V-iso generatrix keeps parameter
    /// speed 1 while a U-iso circle at written `V = height` moves to ElSLib `V`.
    /// BSpline poles are scaled in place (`cxx:278-294`). Plane Circle/Ellipse
    /// use `SetScale(origin, LengthFact)` (`cxx:215-222`). Circle/Ellipse on a
    /// revolution surface go through `Geom2dConvert::CurveToBSplineCurve` first
    /// (`cxx:229-236`) — that convert is not in this crate; those entities keep
    /// the raw 2D curve. Parabola/Hyperbola return unchanged (`cxx:237-252`).
    pub(super) fn degree_to_radian(
        &self,
        c2d: Arc<dyn Curve2d>,
        curve_ref: usize,
        surf: &dyn Surface,
    ) -> Arc<dyn Curve2d> {
        const LENGTH_FACT: f64 = 1.0;
        const ANGLE_FACT: f64 = 1.0;
        let is_plane = classify_surface(surf) == SurfaceKind::Plane;
        // `GeomConvert_Units.cxx:191-227`. Offset (and other non-analytic
        // kinds) return the pcurve unchanged. Revolution is U-angle only.
        let (u_fact, v_fact) = if let Some((_, alpha)) = surf.cone_ref() {
            (ANGLE_FACT, LENGTH_FACT / alpha.cos())
        } else if surf.is_surface_of_revolution() {
            (ANGLE_FACT, LENGTH_FACT)
        } else {
            match classify_surface(surf) {
                SurfaceKind::Sphere | SurfaceKind::Torus => (ANGLE_FACT, ANGLE_FACT),
                SurfaceKind::Cylinder => (ANGLE_FACT, LENGTH_FACT),
                SurfaceKind::Plane => (LENGTH_FACT, LENGTH_FACT),
                SurfaceKind::Cone | SurfaceKind::Other => return c2d,
            }
        };
        let Ok(rec) = self.record(curve_ref) else {
            return c2d;
        };
        match rec.type_name.as_str() {
            "CIRCLE" | "ELLIPSE" if is_plane => {
                let mut t = occt_core::gp::GpTrsf2d::default();
                if t.set_scale(&GpPnt2d::new(0.0, 0.0), LENGTH_FACT).is_err() {
                    return c2d;
                }
                let mut scaled = c2d.clone_dyn();
                scaled.transform(&t);
                Arc::from(scaled)
            }
            "PARABOLA" | "HYPERBOLA" => c2d,
            _ if (u_fact - 1.0).abs() <= 1e-15 && (v_fact - 1.0).abs() <= 1e-15 => c2d,
            "LINE" => {
                let loc = c2d.d0(0.0);
                let new_loc = GpPnt2d::new(loc.x() * u_fact, loc.y() * v_fact);
                let (_, tan) = c2d.d1(0.0);
                let Ok(new_dir) = GpDir2d::new(tan.x() * u_fact, tan.y() * v_fact) else {
                    return c2d;
                };
                Arc::new(Geom2dLine::from_pnt_dir(new_loc, new_dir))
            }
            "B_SPLINE_CURVE_WITH_KNOTS" | "B_SPLINE_CURVE" | "POLYLINE" => self
                .scale_pcurve_bspline_poles(curve_ref, u_fact, v_fact)
                .unwrap_or(c2d),
            _ => c2d,
        }
    }

    /// `GeomConvert_Units.cxx:278-294` — affinity on BSpline poles.
    fn scale_pcurve_bspline_poles(
        &self,
        curve_ref: usize,
        u_fact: f64,
        v_fact: f64,
    ) -> Option<Arc<dyn Curve2d>> {
        let rec = self.record(curve_ref).ok()?;
        let (pts, knots, degree): (Vec<GpPnt2d>, Vec<f64>, usize) = match rec.type_name.as_str() {
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = parse_f64(&rec.args[1]).ok()? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()
                    .ok()?;
                let knots = expand_knots(
                    &parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()")),
                );
                (pts, knots, degree)
            }
            "B_SPLINE_CURVE" => {
                let degree = parse_f64(&rec.args[1]).ok()? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()
                    .ok()?;
                let knots = uniform_knots_for(pts.len(), degree);
                (pts, knots, degree)
            }
            "POLYLINE" => {
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()
                    .ok()?;
                if pts.len() < 2 {
                    return None;
                }
                let knots = uniform_knots_for(pts.len(), 1);
                (pts, knots, 1)
            }
            _ => return None,
        };
        let xs: Vec<f64> = pts.iter().map(|p| p.x() * u_fact).collect();
        let ys: Vec<f64> = pts.iter().map(|p| p.y() * v_fact).collect();
        Geom2dBSplineCurve::new(xs, ys, knots, degree)
            .ok()
            .map(|c| Arc::new(c) as Arc<dyn Curve2d>)
    }
}

