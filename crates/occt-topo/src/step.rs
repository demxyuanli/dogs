//! Phase 4 module: step — STEP (ISO 10303-21) physical-file exchange.
//!
//! Ports `STEPControl_Writer` / `STEPControl_Reader` for the BRep port.
//! Writes and reads the classic EXPRESS entity set for B-rep solids:
//! `CARTESIAN_POINT`, `DIRECTION`, `VECTOR`, `AXIS2_PLACEMENT_3D`,
//! `LINE`, `CIRCLE`, `ELLIPSE`, `PARABOLA`, `PLANE`, `CYLINDRICAL_SURFACE`,
//! `CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`,
//! `VERTEX_POINT`, `EDGE_CURVE`, `ORIENTED_EDGE`, `EDGE_LOOP`,
//! `FACE_OUTER_BOUND`, `ADVANCED_FACE`, `CLOSED_SHELL`,
//! `MANIFOLD_SOLID_BREP`, plus the product/representation scaffolding
//! (`PRODUCT`, `PRODUCT_DEFINITION`, `PRODUCT_DEFINITION_SHAPE`,
//! `SHAPE_REPRESENTATION`, `PRODUCT_DEFINITION_SHAPE_REPRESENTATION`,
//! `ADVANCED_BREP_SHAPE_REPRESENTATION`).
//!
//! Curves and surfaces cannot be downcast from `Arc<dyn Curve>` /
//! `Arc<dyn Surface>`, so geometry is classified by sampling invariants
//! (constant zero second derivative ⇒ line, constant curvature + periodic ⇒
//! circle, planar + unbounded ⇒ plane, equidistant samples ⇒ sphere, ...).
//! This mirrors `GeomAdaptor`'s type tag at a slightly higher cost.

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{
    GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpElips, GpHypr, GpLin, GpParab, GpPln, GpPnt,
    GpSphere, GpTorus, GpVec, GpXyz,
};
use occt_geom::{
    Curve, GeomBSplineCurve, GeomCircle, GeomCone, GeomCylinder, GeomEllipse, GeomHyperbola,
    GeomLine, GeomParabola, GeomPlane, GeomSphere, GeomTorus, Surface,
};

use crate::abs::ShapeType;
use crate::brep_surface::{classify_surface, face_plane, sphere_center, SurfaceKind};
use crate::builder::TopoBuilder;
use crate::model::BRepModel;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, edge_vertices, wires_of_face};

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Render an f64 in STEP real syntax (always carries a decimal point in the
/// mantissa; the shortest round-tripping form).
fn step_real(x: f64) -> String {
    if !x.is_finite() {
        return "0.".into();
    }
    let s = format!("{:?}", x);
    // Rust's `{:?}` gives e.g. "0.0", "2.5", "1e-9", "1e20", "inf".
    if let Some(pos) = s.find(['e', 'E']) {
        let mant = &s[..pos];
        if !mant.contains('.') {
            return format!("{mant}.{}", &s[pos..]);
        }
        s
    } else if s.contains('.') {
        s
    } else {
        format!("{s}.")
    }
}

/// Escape a label for a STEP string literal (single quotes doubled).
fn esc_str(s: &str) -> String {
    s.replace('\'', "''")
}

/// `#1,#2,#3` from a slice of ids.
fn join_refs(ids: &[usize]) -> String {
    ids.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",")
}

fn dir_x() -> GpDir {
    GpDir::new(1.0, 0.0, 0.0).unwrap()
}
fn dir_y() -> GpDir {
    GpDir::new(0.0, 1.0, 0.0).unwrap()
}
fn dir_z() -> GpDir {
    GpDir::new(0.0, 0.0, 1.0).unwrap()
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        0.5 * (a.x() + b.x()),
        0.5 * (a.y() + b.y()),
        0.5 * (a.z() + b.z()),
    )
}

/// 3×3 determinant.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Circumcenter of three non-collinear 3D points (perpendicular-bisector
/// system solved by Cramer's rule).
fn circle_center3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.square_modulus() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.x() * p.x() + p.y() * p.y() + p.z() * p.z();
    let mat = [
        [d1.x(), d1.y(), d1.z()],
        [d2.x(), d2.y(), d2.z()],
        [n.x, n.y, n.z],
    ];
    let rhs = [
        0.5 * (n2(b) - n2(a)),
        0.5 * (n2(c) - n2(a)),
        a.x() * n.x + a.y() * n.y + a.z() * n.z,
    ];
    let d = det3(&mat);
    if d.abs() < 1e-30 {
        return None;
    }
    let mut o = [0.0; 3];
    for k in 0..3 {
        let mut m = mat;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        o[k] = det3(&m) / d;
    }
    Some(GpPnt::new(o[0], o[1], o[2]))
}

// ---------------------------------------------------------------------------
// Geometry classification
// ---------------------------------------------------------------------------

/// Analytic curve family, classified by sampling (GeomAdaptor substitute).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurveKind {
    Line,
    Circle,
    Ellipse,
    Parabola,
    Other, // hyperbola / B-spline / unsupported
}

fn classify_curve(c: &dyn Curve, a: f64, b: f64) -> CurveKind {
    let (lo, hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    };
    let span = hi - lo;
    let n = 6;
    let mut d2s = Vec::with_capacity(n);
    for i in 0..n {
        let u = lo + span * i as f64 / (n - 1) as f64;
        d2s.push(c.d2(u).2.magnitude());
    }
    let max_d2 = d2s.iter().cloned().fold(0.0_f64, f64::max);
    if max_d2 < 1e-9 {
        return CurveKind::Line;
    }
    let min_d2 = d2s.iter().cloned().fold(f64::INFINITY, f64::min);
    if (max_d2 - min_d2) / max_d2 < 0.02 {
        // constant |d²|: circle (closed) or parabola (open)
        if c.is_periodic() {
            CurveKind::Circle
        } else {
            CurveKind::Parabola
        }
    } else if c.is_periodic() {
        CurveKind::Ellipse
    } else {
        CurveKind::Other
    }
}

/// Analytic surface family. Plane/sphere reuse `brep_surface`'s tested
/// classifiers; cylinder/cone/torus are detected from constant-v rings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

fn surf_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| {
        if a.is_finite() && b.is_finite() && b > a {
            (a, b)
        } else {
            (-1.0, 1.0)
        }
    };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// Center + radius of the ring obtained by sweeping `u` at fixed `v`.
fn ring_center_radius(s: &dyn Surface, u0: f64, u1: f64, v: f64) -> Option<(GpPnt, f64)> {
    let span = u1 - u0;
    let p0 = s.d0(u0, v);
    let p1 = s.d0(u0 + 0.25 * span, v);
    let p2 = s.d0(u0 + 0.5 * span, v);
    let c = circle_center3(&p0, &p1, &p2)?;
    Some((c, c.distance(&p0)))
}

/// Axis direction of a ruled-of-revolution surface from a constant-v ring's
/// plane normal.
fn ring_axis(s: &dyn Surface, u0: f64, u1: f64, v: f64) -> Option<GpDir> {
    let span = u1 - u0;
    let p0 = s.d0(u0, v);
    let p1 = s.d0(u0 + 0.25 * span, v);
    let p2 = s.d0(u0 + 0.5 * span, v);
    let n = GpVec::from_pnts(&p0, &p1)
        .xyz()
        .crossed(GpVec::from_pnts(&p0, &p2).xyz());
    GpDir::from_vec(&GpVec::from_xyz(&n)).ok()
}

/// Whether the v-lines (constant u) are straight — cone generators are,
/// torus minor circles are not.
fn v_line_straight(s: &dyn Surface, u: f64, v0: f64, v1: f64) -> bool {
    let span = v1 - v0;
    let p0 = s.d0(u, v0);
    let p1 = s.d0(u, v0 + 0.25 * span);
    let p2 = s.d0(u, v0 + 0.5 * span);
    let a = GpVec::from_pnts(&p0, &p1);
    let b = GpVec::from_pnts(&p0, &p2);
    let cross = a.xyz().crossed(b.xyz()).modulus();
    cross < 1e-6 * (a.magnitude() * b.magnitude()).max(1e-12)
}

fn surface_kind(s: &dyn Surface) -> SurfKind {
    match classify_surface(s) {
        SurfaceKind::Plane => return SurfKind::Plane,
        SurfaceKind::Sphere => return SurfKind::Sphere,
        _ => {}
    }
    let (u0, u1, v0, v1) = surf_bounds(s);
    let vm = 0.5 * (v0 + v1);
    let v2 = vm + 0.25 * (v1 - v0);
    if let (Some((c0, r0)), Some((c1, r1))) = (
        ring_center_radius(s, u0, u1, vm),
        ring_center_radius(s, u0, u1, v2),
    ) {
        if (r1 - r0).abs() <= 1e-6 * r0.abs().max(1.0) && c0.distance(&c1) > 1e-9 {
            return SurfKind::Cylinder;
        }
        if v_line_straight(s, u0, v0, v1) {
            return SurfKind::Cone;
        }
        return SurfKind::Torus;
    }
    SurfKind::Other
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Incremental STEP writer — emits `#N=...;` records with an internal counter.
pub struct StepWriter {
    next_id: usize,
    lines: Vec<String>,
}

impl StepWriter {
    pub fn new() -> Self {
        Self {
            next_id: 1,
            lines: Vec::new(),
        }
    }

    /// Emit one entity body (`TYPE(a,b,...)`) and return its record id.
    pub fn emit(&mut self, entity: String) -> usize {
        let id = self.next_id;
        self.next_id += 1;
        self.lines.push(format!("#{id}={entity};"));
        id
    }

    pub fn add_cartesian_point(&mut self, p: &GpPnt) -> usize {
        self.emit(format!(
            "CARTESIAN_POINT('',({},{},{}))",
            step_real(p.x()),
            step_real(p.y()),
            step_real(p.z())
        ))
    }

    pub fn add_direction(&mut self, d: &GpDir) -> usize {
        self.emit(format!(
            "DIRECTION('',({},{},{}))",
            step_real(d.x()),
            step_real(d.y()),
            step_real(d.z())
        ))
    }

    /// VECTOR(name, direction, magnitude).
    pub fn add_vector(&mut self, d: &GpDir, magnitude: f64) -> usize {
        let did = self.add_direction(d);
        self.emit(format!("VECTOR('',#{did},{})", step_real(magnitude)))
    }

    /// AXIS2_PLACEMENT_3D(name, location, axis, ref_direction).
    pub fn add_axis2_placement_3d(&mut self, loc: &GpPnt, axis: &GpDir, ref_dir: &GpDir) -> usize {
        let lid = self.add_cartesian_point(loc);
        let aid = self.add_direction(axis);
        let rid = self.add_direction(ref_dir);
        self.emit(format!("AXIS2_PLACEMENT_3D('',#{lid},#{aid},#{rid})"))
    }

    /// Assemble the complete physical file.
    pub fn finish(self) -> String {
        let mut out = String::new();
        out.push_str("ISO-10303-21;\n");
        out.push_str("HEADER;\n");
        out.push_str("FILE_DESCRIPTION(('BRep'),'2;1');\n");
        out.push_str("FILE_NAME('model.step','2026-07-31T00:00:00',('',''),(''),'','rust','');\n");
        out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n");
        out.push_str("ENDSEC;\n");
        out.push_str("DATA;\n");
        for l in &self.lines {
            out.push_str(l);
            out.push('\n');
        }
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        out
    }
}

/// Internal writer state: the entity writer plus identity maps so shared
/// vertices/edges/curves are emitted exactly once.
struct WriteCtx {
    w: StepWriter,
    vertex_ids: HashMap<usize, usize>,
    edge_ids: HashMap<usize, usize>,
    curve_ids: HashMap<usize, usize>,
    geom_ctx: usize,
    prod_ctx: usize,
    def_ctx: usize,
}

impl WriteCtx {
    fn new() -> Self {
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
        }
    }

    fn finish(self) -> String {
        self.w.finish()
    }

    fn emit_vertex(&mut self, v: &Vertex) -> usize {
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

    fn emit_edge(&mut self, e: &Edge) -> usize {
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

    fn emit_edge_curve(&mut self, e: &Edge) -> usize {
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

    fn emit_curve_entity(&mut self, c: &dyn Curve, a: f64, b: f64) -> usize {
        match classify_curve(c, a, b) {
            CurveKind::Line => {
                let origin = c.d0(0.0);
                let d1 = c.d1(0.0).1;
                let dir = GpDir::from_vec(&d1).unwrap_or(dir_x());
                let pid = self.w.add_cartesian_point(&origin);
                let vid = self.w.add_vector(&dir, 1.0);
                self.w.emit(format!("LINE('',#{pid},#{vid})"))
            }
            CurveKind::Circle => {
                let p0 = c.d0(a);
                let p1 = c.d0(a + PI / 2.0);
                let p2 = c.d0(a + PI);
                let center = circle_center3(&p0, &p1, &p2).unwrap_or_else(GpPnt::zero);
                let r = center.distance(&p0);
                let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &p0)).unwrap_or(dir_x());
                let ydir = GpDir::from_vec(&GpVec::from_pnts(&center, &p1)).unwrap_or(dir_y());
                let axis = xdir.crossed(&ydir).unwrap_or(dir_z());
                let ax = self.w.add_axis2_placement_3d(&center, &axis, &xdir);
                self.w.emit(format!("CIRCLE('',#{ax},{})", step_real(r)))
            }
            CurveKind::Ellipse => {
                let p0 = c.d0(a);
                let p_half = c.d0(a + PI / 2.0);
                let p_pi = c.d0(a + PI);
                let center = midpoint(&p0, &p_pi);
                let major = center.distance(&p0);
                let minor = center.distance(&p_half);
                let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &p0)).unwrap_or(dir_x());
                let ydir = GpDir::from_vec(&GpVec::from_pnts(&p_half, &center)).unwrap_or(dir_y());
                let axis = xdir.crossed(&ydir).unwrap_or(dir_z());
                let ax = self.w.add_axis2_placement_3d(&center, &axis, &xdir);
                self.w.emit(format!(
                    "ELLIPSE('',#{ax},{},{})",
                    step_real(major),
                    step_real(minor)
                ))
            }
            CurveKind::Parabola => {
                let vertex = c.d0(0.0);
                let d1 = c.d1(0.0).1;
                let d2 = c.d2(0.0).2;
                let f = 0.5 / d2.magnitude().max(1e-30);
                let xdir = GpDir::from_vec(&d2).unwrap_or(dir_x());
                let ydir = GpDir::from_vec(&d1).unwrap_or(dir_y());
                let axis = xdir.crossed(&ydir).unwrap_or(dir_z());
                let ax = self.w.add_axis2_placement_3d(&vertex, &axis, &xdir);
                self.w.emit(format!("PARABOLA('',#{ax},{})", step_real(f)))
            }
            CurveKind::Other => {
                // Fallback: a tangent line at the start parameter keeps the
                // file valid; the geometry is a linear approximation.
                let p0 = c.d0(a);
                let d1 = c.d1(a).1;
                let dir = GpDir::from_vec(&d1).unwrap_or(dir_x());
                let pid = self.w.add_cartesian_point(&p0);
                let vid = self.w.add_vector(&dir, 1.0);
                self.w.emit(format!("LINE('',#{pid},#{vid})"))
            }
        }
    }

    fn emit_wire(&mut self, w: &Wire) -> usize {
        let edges = edges_of_wire(w);
        let mut refs = Vec::with_capacity(edges.len());
        for e in &edges {
            let edge_id = self.emit_edge(e);
            let (v1, v2) = edge_vertices(e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                continue;
            };
            let k1 = Arc::as_ptr(&v1.0.tshape) as usize;
            let k2 = Arc::as_ptr(&v2.0.tshape) as usize;
            let (Some(&v1_id), Some(&v2_id)) = (self.vertex_ids.get(&k1), self.vertex_ids.get(&k2))
            else {
                continue;
            };
            refs.push(self.w.emit(format!("ORIENTED_EDGE('',#{v1_id},#{v2_id},#{edge_id},.T.)")));
        }
        self.w.emit(format!("EDGE_LOOP('',({}))", join_refs(&refs)))
    }

    fn emit_face(&mut self, f: &Face) -> usize {
        let surf_ref = self.emit_surface(f);
        let wires = wires_of_face(f);
        let mut bounds = Vec::with_capacity(wires.len());
        for w in &wires {
            let loop_ref = self.emit_wire(w);
            bounds.push(self.w.emit(format!("FACE_OUTER_BOUND('',#{loop_ref},.T.)")));
        }
        self.w.emit(format!("ADVANCED_FACE('',#{surf_ref},({}),.T.)", join_refs(&bounds)))
    }

    fn emit_surface(&mut self, f: &Face) -> usize {
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
            SurfKind::Other => self.emit_plane_fallback(),
        }
    }

    fn emit_ax2(&mut self, ax2: &GpAx2) -> usize {
        self.w
            .add_axis2_placement_3d(&ax2.location(), &ax2.direction(), &ax2.x_direction())
    }

    fn emit_plane_fallback(&mut self) -> usize {
        let ax = self.w.add_axis2_placement_3d(&GpPnt::zero(), &dir_z(), &dir_x());
        self.w.emit(format!("PLANE('',#{ax})"))
    }

    fn emit_shell(&mut self, sh: &Shell) -> usize {
        let faces = children_of_type(&sh.0, ShapeType::Face);
        let refs: Vec<usize> = faces
            .iter()
            .map(|f| self.emit_face(&Face(f.clone())))
            .collect();
        self.w
            .emit(format!("CLOSED_SHELL('',({}))", join_refs(&refs)))
    }

    fn emit_solid(&mut self, s: &Solid) -> usize {
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
    fn emit_top(&mut self, shape: &TopoShape) -> Vec<usize> {
        match shape.shape_type() {
            ShapeType::Compound => {
                let kids: Vec<TopoShape> = shape
                    .tshape
                    .read()
                    .unwrap()
                    .children
                    .iter()
                    .map(|h| TopoShape::from_handle(h.clone()))
                    .collect();
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

    fn write_shape_named(&mut self, name: &str, shape: &TopoShape) {
        let items = self.emit_top(shape);
        if items.is_empty() {
            return;
        }
        let n = esc_str(name);
        let rep = self.w.emit(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('{n}',({}),#{})",
            join_refs(&items),
            self.geom_ctx
        ));
        let product = self
            .w
            .emit(format!("PRODUCT('{n}','{n}','',({}))", self.prod_ctx));
        let formation = self
            .w
            .emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        self.w.emit(format!(
            "PRODUCT_DEFINITION('','','',#{formation},#{})",
            self.def_ctx
        ));
        let pds = self.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{product})"));
        self.w
            .emit(format!("PRODUCT_DEFINITION_SHAPE_REPRESENTATION('',#{pds},#{rep})"));
    }
}

/// Direct children of `s` whose type is `t`.
fn children_of_type(s: &TopoShape, t: ShapeType) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.read().unwrap().shape_type == t)
        .map(|h| TopoShape::from_handle(h.clone()))
        .collect()
}

/// Extract (axis2, radius) for a cylinder surface by sampling two rings.
fn cylinder_params(s: &dyn Surface) -> Option<(GpAx2, f64)> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let vm = 0.5 * (v0 + v1);
    let (c0, r) = ring_center_radius(s, u0, u1, vm)?;
    let axis = ring_axis(s, u0, u1, vm).unwrap_or(dir_z());
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&c0, &s.d0(u0, vm))).unwrap_or(dir_x());
    let ax2 = GpAx2::new(c0, axis, xdir).ok()?;
    Some((ax2, r))
}

/// Extract (axis2, radius, semi_angle) for a cone surface.
fn cone_params(s: &dyn Surface) -> Option<(GpAx2, f64, f64)> {
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
fn torus_params(s: &dyn Surface) -> Option<(GpAx2, f64, f64)> {
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

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// A parsed `#N=TYPE(...)` data record.
#[derive(Debug, Clone)]
struct Record {
    type_name: String,
    args: Vec<String>,
}

/// Split a comma-separated argument list at the top nesting level, respecting
/// strings (with `''` escapes), parens and brackets.
fn split_top(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0usize;
    let mut cur = String::new();
    let mut chars = s.chars().peekable();
    let mut in_str = false;
    while let Some(c) = chars.next() {
        if in_str {
            cur.push(c);
            if c == '\'' {
                if chars.peek() == Some(&'\'') {
                    cur.push(chars.next().unwrap());
                } else {
                    in_str = false;
                }
            }
        } else {
            match c {
                '\'' => {
                    in_str = true;
                    cur.push(c);
                }
                '(' | '[' => {
                    depth += 1;
                    cur.push(c);
                }
                ')' | ']' => {
                    depth -= 1;
                    cur.push(c);
                }
                ',' if depth == 0 => {
                    out.push(cur.trim().to_string());
                    cur.clear();
                }
                _ => cur.push(c),
            }
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_string());
    }
    out
}

/// Parse the entity body `TYPE(a,b,...)` into its type name and argument list.
fn parse_entity_body(body: &str) -> (String, Vec<String>) {
    let body = body.trim();
    let open = body.find('(').unwrap_or(body.len());
    let type_name = body[..open].trim().to_string();
    let close = body.rfind(')').unwrap_or(body.len());
    let inner = if close > open + 1 {
        &body[open + 1..close]
    } else {
        ""
    };
    (type_name, split_top(inner))
}

/// Split the DATA section into `#id=TYPE(...)` records.
fn parse_records(data: &str) -> Result<HashMap<usize, Record>, String> {
    let bytes = data.as_bytes();
    let mut records = HashMap::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j == i + 1 {
                i += 1;
                continue;
            }
            let id: usize = data[i + 1..j]
                .parse()
                .map_err(|_| format!("bad entity id at offset {i}"))?;
            let mut k = j;
            while k < bytes.len() && bytes[k].is_ascii_whitespace() {
                k += 1;
            }
            if k >= bytes.len() || bytes[k] != b'=' {
                return Err(format!("malformed record #{id}: expected '='"));
            }
            let (body, next) = parse_entity_body_text(data, k + 1)?;
            let (type_name, args) = parse_entity_body(&body);
            records.insert(id, Record { type_name, args });
            i = next;
        } else {
            i += 1;
        }
    }
    Ok(records)
}

/// Scan one entity body starting after the `=`, stopping at the terminating
/// `;`. Returns the body text and the index just past the `;`.
fn parse_entity_body_text(data: &str, start: usize) -> Result<(String, usize), String> {
    let bytes = data.as_bytes();
    let mut i = start;
    let mut depth = 0usize;
    let mut in_str = false;
    while i < bytes.len() {
        let c = bytes[i] as char;
        if in_str {
            if c == '\'' {
                if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                    i += 2;
                    continue;
                }
                in_str = false;
            }
        } else {
            match c {
                '\'' => in_str = true,
                '(' => depth += 1,
                ')' => {
                    if depth == 0 {
                        return Err("unbalanced parens in entity body".into());
                    }
                    depth -= 1;
                    if depth == 0 {
                        let body = data[start..=i].to_string();
                        let mut j = i + 1;
                        while j < bytes.len() && bytes[j].is_ascii_whitespace() {
                            j += 1;
                        }
                        if j >= bytes.len() || bytes[j] != b';' {
                            return Err("missing ';' after entity body".into());
                        }
                        return Ok((body, j + 1));
                    }
                }
                _ => {}
            }
        }
        i += 1;
    }
    Err("unterminated entity body".into())
}

fn parse_ref(s: &str) -> Option<usize> {
    s.trim().strip_prefix('#')?.trim().parse().ok()
}

fn parse_ref_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| parse_ref(&a))
        .collect()
}

fn parse_f64(s: &str) -> Result<f64, String> {
    s.trim()
        .parse()
        .map_err(|_| format!("bad real literal '{s}'"))
}

fn parse_xyz(s: &str) -> Result<GpXyz, String> {
    let s = s.trim();
    let inner = s.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<String> = split_top(inner)
        .into_iter()
        .map(|p| p.trim().to_string())
        .collect();
    if parts.len() != 3 {
        return Err(format!("expected 3-component tuple, got '{s}'"));
    }
    Ok(GpXyz::new(
        parse_f64(&parts[0])?,
        parse_f64(&parts[1])?,
        parse_f64(&parts[2])?,
    ))
}

fn parse_str(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        s[1..s.len() - 1].replace("''", "'")
    } else {
        s.to_string()
    }
}

/// Reference resolver with memoization; unknown/unsupported entities are
/// recorded as warnings and skipped.
struct Resolver<'a> {
    records: &'a HashMap<usize, Record>,
    b: TopoBuilder,
    shape_cache: RefCell<HashMap<usize, TopoShape>>,
    point_cache: RefCell<HashMap<usize, GpPnt>>,
    dir_cache: RefCell<HashMap<usize, GpDir>>,
    axis_cache: RefCell<HashMap<usize, GpAx2>>,
    curve_cache: RefCell<HashMap<usize, Arc<dyn Curve>>>,
    surface_cache: RefCell<HashMap<usize, Arc<dyn Surface>>>,
    resolving: RefCell<HashSet<usize>>,
    warnings: RefCell<Vec<String>>,
}

impl<'a> Resolver<'a> {
    fn new(records: &'a HashMap<usize, Record>) -> Self {
        Self {
            records,
            b: TopoBuilder::new(),
            shape_cache: RefCell::new(HashMap::new()),
            point_cache: RefCell::new(HashMap::new()),
            dir_cache: RefCell::new(HashMap::new()),
            axis_cache: RefCell::new(HashMap::new()),
            curve_cache: RefCell::new(HashMap::new()),
            surface_cache: RefCell::new(HashMap::new()),
            resolving: RefCell::new(HashSet::new()),
            warnings: RefCell::new(Vec::new()),
        }
    }

    fn record(&self, id: usize) -> Result<&'a Record, String> {
        self.records
            .get(&id)
            .ok_or_else(|| format!("reference to undefined entity #{id}"))
    }

    fn warn(&self, msg: String) {
        self.warnings.borrow_mut().push(msg);
    }

    fn resolve_shape(&self, id: usize) -> Result<TopoShape, String> {
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
            "FACE_OUTER_BOUND" => self.resolve_outer_bound(rec),
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

    fn resolve_vertex(&self, rec: &'a Record) -> Result<TopoShape, String> {
        if rec.args.len() < 2 {
            return Err("VERTEX_POINT: bad args".into());
        }
        let pid = parse_ref(&rec.args[1]).ok_or("VERTEX_POINT: bad point ref")?;
        let p = self.resolve_point(pid)?;
        Ok(self.b.make_vertex(p, 0.0).0)
    }

    fn resolve_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
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
        let (first, last) = edge_params_for_curve(curve.as_ref(), &p1, &p2);
        let mut e = self.b.make_edge(curve, first, last);
        self.b.add(&mut e.0, &v1);
        self.b.add(&mut e.0, &v2);
        Ok(e.0)
    }

    fn resolve_oriented_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let edge_ref = parse_ref(&rec.args[3]).ok_or("ORIENTED_EDGE: bad edge ref")?;
        self.resolve_shape(edge_ref)
    }

    fn resolve_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let items = parse_ref_list(&rec.args[1]);
        let mut edges = Vec::with_capacity(items.len());
        for &it in &items {
            let s = self.resolve_shape(it)?;
            if !s.is_edge() {
                return Err(format!("#{it}: expected EDGE in EDGE_LOOP"));
            }
            edges.push(Edge(s));
        }
        Ok(self.b.make_wire(&edges).0)
    }

    fn resolve_outer_bound(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let loop_ref = parse_ref(&rec.args[1]).ok_or("FACE_OUTER_BOUND: bad loop ref")?;
        let s = self.resolve_shape(loop_ref)?;
        if !s.is_wire() {
            return Err("FACE_OUTER_BOUND: loop is not a wire".into());
        }
        Ok(s)
    }

    fn resolve_face(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let surf_ref = parse_ref(&rec.args[1]).ok_or("ADVANCED_FACE: bad surface ref")?;
        let surface = self.resolve_surface(surf_ref)?;
        let bounds = parse_ref_list(&rec.args[2]);
        let mut wires = Vec::with_capacity(bounds.len());
        for &b in &bounds {
            let s = self.resolve_shape(b)?;
            if !s.is_wire() {
                return Err(format!("#{b}: expected wire in face bounds"));
            }
            wires.push(Wire(s));
        }
        Ok(self.b.make_face(surface, &wires).0)
    }

    fn resolve_shell(&self, rec: &'a Record) -> Result<TopoShape, String> {
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

    fn resolve_solid(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let outer = parse_ref(&rec.args[1]).ok_or("MANIFOLD_SOLID_BREP: bad shell ref")?;
        let s = self.resolve_shape(outer)?;
        if !s.is_shell() {
            return Err("MANIFOLD_SOLID_BREP: outer is not a shell".into());
        }
        Ok(self.b.make_solid(&[Shell(s)]).0)
    }

    fn resolve_point(&self, id: usize) -> Result<GpPnt, String> {
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

    fn resolve_direction(&self, id: usize) -> Result<GpDir, String> {
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

    fn resolve_vector(&self, id: usize) -> Result<GpVec, String> {
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

    fn resolve_axis2(&self, id: usize) -> Result<GpAx2, String> {
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

    fn resolve_curve(&self, id: usize) -> Result<Arc<dyn Curve>, String> {
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
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let mults = parse_usize_list(&rec.args[3]);
                let knot_vals = parse_real_list(&rec.args[4]);
                let mut knots = Vec::new();
                for (i, &m) in mults.iter().enumerate() {
                    let kv = knot_vals.get(i).copied().unwrap_or(0.0);
                    for _ in 0..m {
                        knots.push(kv);
                    }
                }
                Arc::new(
                    GeomBSplineCurve::new(poles, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE: {e}"))?,
                )
            }
            other => {
                self.warn(format!("unsupported curve entity {other} (#{id})"));
                return Err(format!("unsupported curve entity {other} (#{id})"));
            }
        };
        self.curve_cache.borrow_mut().insert(id, curve.clone());
        Ok(curve)
    }

    fn resolve_surface(&self, id: usize) -> Result<Arc<dyn Surface>, String> {
        if let Some(s) = self.surface_cache.borrow().get(&id) {
            return Ok(s.clone());
        }
        let rec = self.record(id)?;
        let surface: Arc<dyn Surface> = match rec.type_name.as_str() {
            "PLANE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("PLANE: bad axis ref")?;
                Arc::new(GeomPlane::new(GpPln::new(self.resolve_axis2(ax)?.to_ax3())))
            }
            "CYLINDRICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CYLINDRICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomCylinder::new(
                    GpCylinder::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("CYLINDRICAL_SURFACE: {e}"))?,
                ))
            }
            "CONICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CONICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                let a = parse_f64(&rec.args[3])?;
                Arc::new(GeomCone::new(
                    GpCone::new(self.resolve_axis2(ax)?.to_ax3(), r, a)
                        .map_err(|e| format!("CONICAL_SURFACE: {e}"))?,
                ))
            }
            "SPHERICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("SPHERICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomSphere::new(
                    GpSphere::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("SPHERICAL_SURFACE: {e}"))?,
                ))
            }
            "TOROIDAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("TOROIDAL_SURFACE: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                Arc::new(GeomTorus::new(
                    GpTorus::new(self.resolve_axis2(ax)?.to_ax3(), maj, min)
                        .map_err(|e| format!("TOROIDAL_SURFACE: {e}"))?,
                ))
            }
            other => {
                self.warn(format!("unsupported surface entity {other} (#{id})"));
                return Err(format!("unsupported surface entity {other} (#{id})"));
            }
        };
        self.surface_cache.borrow_mut().insert(id, surface.clone());
        Ok(surface)
    }

    /// Resolve a representation record into (name, shapes).
    fn resolve_representation(&self, id: usize) -> Result<(String, Vec<TopoShape>), String> {
        let rec = self.record(id)?;
        let name = parse_str(&rec.args[0]);
        let items = parse_ref_list(&rec.args[1]);
        let mut shapes = Vec::with_capacity(items.len());
        for &it in &items {
            shapes.push(self.resolve_shape(it)?);
        }
        Ok((name, shapes))
    }
}

/// Compute an edge's parameter range from its endpoint vertex points, based on
/// the reconstructed curve's analytic type.
fn edge_params_for_curve(curve: &dyn Curve, p1: &GpPnt, p2: &GpPnt) -> (f64, f64) {
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    let (lo, hi) = if f.is_finite() && l.is_finite() && l > f {
        (f, l)
    } else {
        (0.0, 1.0)
    };
    match classify_curve(curve, lo, hi) {
        CurveKind::Line => {
            let origin = curve.d0(0.0);
            let d1 = curve.d1(0.0).1;
            let dir = GpVec::from_xyz(&d1.xyz().normalized());
            (
                GpVec::from_pnts(&origin, p1).dot(&dir),
                GpVec::from_pnts(&origin, p2).dot(&dir),
            )
        }
        CurveKind::Circle => {
            if p1.distance(p2) < 1e-9 {
                return (0.0, 2.0 * PI);
            }
            let pa = curve.d0(lo);
            let pb = curve.d0(lo + (hi - lo) / 4.0);
            let pc = curve.d0(lo + (hi - lo) / 2.0);
            let center = circle_center3(&pa, &pb, &pc).unwrap_or_else(GpPnt::zero);
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &pa)).unwrap_or(dir_x());
            let ydir = GpDir::from_vec(&GpVec::from_pnts(&center, &pb)).unwrap_or(dir_y());
            let ang = |p: &GpPnt| {
                let v = GpVec::from_pnts(&center, p);
                v.xyz().dot(ydir.xyz()).atan2(v.xyz().dot(xdir.xyz()))
            };
            let t1 = ang(p1);
            let mut t2 = ang(p2);
            if t2 <= t1 {
                t2 += 2.0 * PI;
            }
            (t1, t2)
        }
        CurveKind::Ellipse => {
            if p1.distance(p2) < 1e-9 {
                return (0.0, 2.0 * PI);
            }
            let pa = curve.d0(lo);
            let pb = curve.d0(lo + (hi - lo) / 4.0);
            let pc = curve.d0(lo + (hi - lo) / 2.0);
            let center = midpoint(&pa, &pc);
            let a_major = center.distance(&pa).max(1e-30);
            let b_minor = center.distance(&pb).max(1e-30);
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &pa)).unwrap_or(dir_x());
            let ydir = GpDir::from_vec(&GpVec::from_pnts(&pb, &center)).unwrap_or(dir_y());
            let ang = |p: &GpPnt| {
                let v = GpVec::from_pnts(&center, p);
                (-v.xyz().dot(ydir.xyz()) / b_minor)
                    .atan2(v.xyz().dot(xdir.xyz()) / a_major)
            };
            let t1 = ang(p1);
            let mut t2 = ang(p2);
            if t2 <= t1 {
                t2 += 2.0 * PI;
            }
            (t1, t2)
        }
        CurveKind::Parabola => {
            let vertex = curve.d0(0.0);
            let d1 = curve.d1(0.0).1;
            let ydir = GpDir::from_vec(&d1).unwrap_or(dir_y());
            let dir = GpVec::from_xyz(ydir.xyz());
            (
                GpVec::from_pnts(&vertex, p1).dot(&dir),
                GpVec::from_pnts(&vertex, p2).dot(&dir),
            )
        }
        CurveKind::Other => {
            if f.is_finite() && l.is_finite() && l > f {
                (f, l)
            } else {
                (0.0, 1.0)
            }
        }
    }
}

fn parse_usize_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

fn parse_real_list(s: &str) -> Vec<f64> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

/// Parse a STEP physical file, returning the model plus collected warnings.
fn read_step_impl(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    if !content.contains("END-ISO-10303-21") {
        return Err("STEP file: missing END-ISO-10303-21 terminator".into());
    }
    let data_start = content
        .find("DATA;")
        .ok_or("STEP file: missing DATA section")?;
    let after_data = &content[data_start + 5..];
    let data_end = after_data
        .find("ENDSEC;")
        .ok_or("STEP file: DATA section not closed by ENDSEC")?;
    let data = &after_data[..data_end];
    let records = parse_records(data)?;
    let resolver = Resolver::new(&records);

    let mut rep_ids: Vec<usize> = records
        .iter()
        .filter(|(_, r)| {
            r.type_name == "ADVANCED_BREP_SHAPE_REPRESENTATION" || r.type_name == "SHAPE_REPRESENTATION"
        })
        .map(|(id, _)| *id)
        .collect();
    rep_ids.sort_unstable();

    let mut model = BRepModel::new();
    if !rep_ids.is_empty() {
        for id in rep_ids {
            match resolver.resolve_representation(id) {
                Ok((name, shapes)) => {
                    if shapes.is_empty() {
                        continue;
                    }
                    let shape = if shapes.len() == 1 {
                        shapes.into_iter().next().unwrap()
                    } else {
                        let b = TopoBuilder::new();
                        b.make_compound_of(&shapes).0
                    };
                    model.add(&name, shape);
                }
                Err(e) => resolver.warn(format!("representation #{id}: {e}")),
            }
        }
    } else {
        // Fallback for bare files without a representation layer: treat any
        // MANIFOLD_SOLID_BREP as a root shape.
        let mut roots: Vec<usize> = records
            .iter()
            .filter(|(_, r)| r.type_name == "MANIFOLD_SOLID_BREP")
            .map(|(id, _)| *id)
            .collect();
        roots.sort_unstable();
        for id in roots {
            match resolver.resolve_shape(id) {
                Ok(s) => {
                    model.add("", s);
                }
                Err(e) => resolver.warn(format!("shape #{id}: {e}")),
            }
        }
    }

    let warnings = resolver.warnings.into_inner();
    Ok((model, warnings))
}

/// Parse a STEP physical file into a `BRepModel`.
pub fn read_step(content: &str) -> Result<BRepModel, String> {
    Ok(read_step_impl(content)?.0)
}

/// Parse a STEP physical file, also returning collected warnings for skipped
/// or unsupported entities.
pub fn read_step_with_warnings(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    read_step_impl(content)
}

/// Read a STEP physical file from disk.
pub fn read_step_file(path: &str) -> Result<BRepModel, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step(&content)
}

/// Round-trip helper: write `shape` to STEP, read it back, and count the
/// distinct vertices / edges / faces in the reconstructed model.
pub fn step_roundtrip_counts(shape: &TopoShape) -> Result<(usize, usize, usize), String> {
    let step = write_shape_step(shape);
    let model = read_step(&step)?;
    let mut nv = 0usize;
    let mut ne = 0usize;
    let mut nf = 0usize;
    for ms in &model.shapes {
        let c = crate::topo_tools_full::shape_counts(&ms.shape);
        nv += c.get(&ShapeType::Vertex).copied().unwrap_or(0);
        ne += c.get(&ShapeType::Edge).copied().unwrap_or(0);
        nf += c.get(&ShapeType::Face).copied().unwrap_or(0);
    }
    Ok((nv, ne, nf))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_surface::SurfaceKind;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::{faces_of, vertices_of};

    #[test]
    fn box_roundtrip_counts_points_and_faces() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let step = write_shape_step(&b.solid.0);
        let model = read_step(&step).expect("parse step");
        assert_eq!(model.len(), 1);

        let (nv, ne, nf) = step_roundtrip_counts(&b.solid.0).expect("roundtrip counts");
        assert_eq!((nv, ne, nf), (8, 12, 6));

        let shape = &model.shapes[0].shape;
        let verts = vertices_of(shape);
        assert_eq!(verts.len(), 8);
        let pts: Vec<GpPnt> = verts
            .iter()
            .map(|v| GeometryRegistry::global().vertex_point(&v.0))
            .collect();
        for &(x, y, z) in &[
            (0., 0., 0.),
            (2., 0., 0.),
            (2., 3., 0.),
            (0., 3., 0.),
            (0., 0., 4.),
            (2., 0., 4.),
            (2., 3., 4.),
            (0., 3., 4.),
        ] {
            assert!(
                pts.iter().any(|p| p.is_equal(&GpPnt::new(x, y, z))),
                "missing corner ({x},{y},{z})"
            );
        }

        let faces = faces_of(shape);
        assert_eq!(faces.len(), 6);
        for f in &faces {
            let surf = GeometryRegistry::global().face_surface(&f.0);
            assert!(surf.is_some(), "face has a surface");
            if let Some(s) = surf {
                assert_eq!(classify_surface(s.as_ref()), SurfaceKind::Plane);
            }
        }
    }

    #[test]
    fn sphere_roundtrip_surface() {
        let s = BRepPrimSphere::make_sphere(2.5);
        let step = write_shape_step(&s.solid.0);
        let model = read_step(&step).expect("parse sphere step");
        assert_eq!(model.len(), 1);
        let shape = &model.shapes[0].shape;
        let faces = faces_of(shape);
        assert_eq!(faces.len(), 1);
        let surf = GeometryRegistry::global()
            .face_surface(&faces[0].0)
            .expect("sphere face has surface");
        assert_eq!(classify_surface(surf.as_ref()), SurfaceKind::Sphere);
    }

    #[test]
    fn model_two_shapes_names_preserved() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let s = BRepPrimSphere::make_sphere(1.0);
        let mut model = BRepModel::new();
        model.add("Box", b.solid.0.clone());
        model.add("Sphere", s.solid.0.clone());
        let step = write_step(&model);
        let got = read_step(&step).expect("parse model step");
        assert_eq!(got.len(), 2);
        let mut names: Vec<&str> = got.names();
        names.sort();
        assert_eq!(names, vec!["Box", "Sphere"]);
    }

    #[test]
    fn malformed_missing_end_iso_errors() {
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=DIRECTION('',(1.,0.,0.));\nENDSEC;\n";
        assert!(read_step(s).is_err());
    }

    #[test]
    fn reader_skips_unknown_entity_with_warning() {
        // A minimal file whose representation references an unsupported entity;
        // the reader must collect a warning and skip the shape gracefully.
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN');\n\
#2=WIDGET_FROBNICATOR('',42.);\n\
#3=GEOMETRIC_REPRESENTATION_CONTEXT('','',3);\n\
#4=ADVANCED_BREP_SHAPE_REPRESENTATION('widget',(#2),#3);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let (model, warnings) = read_step_with_warnings(s).expect("parse");
        assert!(model.is_empty(), "unsupported shape is skipped");
        assert!(
            warnings.iter().any(|w| w.contains("WIDGET_FROBNICATOR")),
            "expected a warning about the unknown entity, got {warnings:?}"
        );
    }

    #[test]
    fn step_real_formats_decimal_point() {
        assert_eq!(step_real(0.0), "0.0");
        assert_eq!(step_real(2.5), "2.5");
        assert_eq!(step_real(-3.0), "-3.0");
        assert!(step_real(1e20).contains('.'));
    }

    #[test]
    fn split_top_handles_nested_lists() {
        let args = split_top("'',(1.,0.,0.),(#2,#3)");
        assert_eq!(args.len(), 3);
        assert_eq!(args[0], "''");
        assert_eq!(args[1], "(1.,0.,0.)");
        assert_eq!(args[2], "(#2,#3)");
    }

    #[test]
    fn cylinder_roundtrip() {
        use crate::primitives::BRepPrimCylinder;
        let c = BRepPrimCylinder::make_cylinder(2.0, 10.0);
        let (nv, ne, nf) = step_roundtrip_counts(&c.solid.0).expect("counts");
        assert_eq!((nv, ne, nf), (2, 3, 3));
        let s = write_shape_step(&c.solid.0);
        let m = read_step(&s).expect("read cylinder");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert_eq!(fs.len(), 3);
        // Two cap faces should classify as planes.
        let planar = fs
            .iter()
            .filter(|f| {
                GeometryRegistry::global()
                    .face_surface(&f.0)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Plane)
                    .unwrap_or(false)
            })
            .count();
        assert_eq!(planar, 2);
    }

    #[test]
    fn step_file_io_roundtrip() {
        use crate::primitives::BRepPrimSphere;
        let s = BRepPrimSphere::make_sphere(1.5);
        let mut model = BRepModel::new();
        model.add("Ball", s.solid.0.clone());
        let path = std::env::temp_dir().join("occt_step_test.step");
        let p = path.to_str().unwrap();
        write_step_file(p, &model).expect("write file");
        let got = read_step_file(p).expect("read file");
        assert_eq!(got.len(), 1);
        assert_eq!(got.names(), vec!["Ball"]);
        std::fs::remove_file(p).ok();
    }
}
