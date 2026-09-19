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

/// Approximate a non-analytic curve with an interpolating B-spline curve.
///
/// Sampled at `n` parameter values across the edge's `[a, b]` range and
/// interpolated with a cubic (degree 1 on failure) clamped B-spline.
pub(super) fn fit_bspline_curve(c: &dyn Curve, a: f64, b: f64) -> Result<GeomBSplineCurve, String> {
    let (lo, hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    };
    let n = 8;
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let u = lo + (hi - lo) * i as f64 / (n - 1) as f64;
        pts.push(c.d0(u));
    }
    crate::loft::interpolate_bspline(&pts, 3)
        .or_else(|_| crate::loft::interpolate_bspline(&pts, 1))
        .map_err(|e| format!("fit_bspline_curve: {e}"))
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
        let curve_ref = self.emit_edge_curve(e);
        let id = self.w.emit(format!("EDGE_CURVE('',#{vid1},#{vid2},#{curve_ref},.T.)"));
        self.edge_ids.insert(key, id);
        id
    }

    pub(super) fn emit_edge_curve(&mut self, e: &Edge) -> usize {
        let Some(curve) = GeometryRegistry::global().edge_curve(&e.0) else {
            // No registered curve: emit a line through the endpoint vertices.
            let (v1, v2) = edge_vertices(e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                let pid = self.w.add_cartesian_point(&GpPnt::zero());
                let vid = self.w.add_vector(&dir_x(), 1.0);
                return self.w.emit(format!("LINE('',#{pid},#{vid})"));
            };
            let p1 = GeometryRegistry::global().vertex_point(&v1.0);
            let p2 = GeometryRegistry::global().vertex_point(&v2.0);
            let d = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).unwrap_or(dir_x());
            let pid = self.w.add_cartesian_point(&p1);
            let vid = self.w.add_vector(&d, 1.0);
            return self.w.emit(format!("LINE('',#{pid},#{vid})"));
        };
        let ckey = Arc::as_ptr(&curve) as *const () as usize;
        if let Some(&id) = self.curve_ids.get(&ckey) {
            return id;
        }
        let (a, b) = GeometryRegistry::global().edge_parameters(&e.0);
        let id = self.emit_curve_entity(curve.as_ref(), a, b);
        self.curve_ids.insert(ckey, id);
        id
    }

    pub(super) fn emit_curve_entity(&mut self, c: &dyn Curve, a: f64, b: f64) -> usize {
        match classify_curve(c, a, b) {
            CurveKind::Line => {
                let origin = c.d0(0.0);
                let d1 = c.d1(0.0).1;
                let dir = GpDir::from_vec(&d1).unwrap_or(dir_x());
                let pid = self.w.add_cartesian_point(&origin);
                let vid = self.w.add_vector(&dir, 1.0);
                self.w.emit(format!("LINE('',#{pid},#{vid})"))
            }
            CurveKind::Circle => emit_circle_entity(&mut self.w, c, a),
            CurveKind::Ellipse => emit_ellipse_entity(&mut self.w, c, a),
            CurveKind::Parabola => emit_parabola_entity(&mut self.w, c),
            CurveKind::Other => {
                // B-spline / trimmed / offset / hyperbola curve: emit a real
                // B-spline when spline output is requested, otherwise fall back
                // to a tangent line at the start parameter so the file stays
                // valid (a linear approximation of the geometry).
                if self.splines {
                    if let Ok(bs) = fit_bspline_curve(c, a, b) {
                        if let Ok(id) = write_bspline_curve(&mut self.w, &bs) {
                            return id;
                        }
                    }
                }
                let p0 = c.d0(a);
                let d1 = c.d1(a).1;
                let dir = GpDir::from_vec(&d1).unwrap_or(dir_x());
                let pid = self.w.add_cartesian_point(&p0);
                let vid = self.w.add_vector(&dir, 1.0);
                self.w.emit(format!("LINE('',#{pid},#{vid})"))
            }
        }
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
