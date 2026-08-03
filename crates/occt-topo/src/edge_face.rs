//! Edge/face full intersection — Phase 17b.
//!
//! Port of `IntTools_EdgeFace` (TKBO): computes the common parts between an
//! edge and a face in 3D space.
//!
//! A common part is either:
//!
//! * an **edge** — a (sub-)range of the edge that lies on the face surface
//!   (coincidence), encoded as [`CommonPartType::Edge`] with a non-degenerate
//!   range;
//! * a **point** — a single parameter at which the edge touches/crosses the
//!   face, encoded as [`CommonPartType::Edge`] with a degenerate range
//!   `[t, t]`.
//!
//! The point/edge distinction mirrors OCCT's `TopAbs_VERTEX` / `TopAbs_EDGE`.
//! The port's [`CommonPartType`] collapses both to `Edge` (per its own doc: "a
//! vertex or coincident arc"), so callers distinguish them via the range length
//! or [`EdgeFace::point_parameters`].
//!
//! Flow (mirrors `IntTools_EdgeFace::Perform`):
//!
//! 1. `check_data` — edge geometry accessibility (degenerated / non-geometric).
//! 2. tolerance preparation (`tol_e + tol_f`, with the B-spline special case);
//! 3. quick coincidence check ([`EdgeFace::is_coincident`]) when requested;
//! 4. [`BeanFaceIntersector`] over the bean range → candidate ranges;
//! 5. `IsProjectable` filter (distance + face 2D restriction);
//! 6. `MakeType` + `CheckTouch` classification per range (point vs edge);
//! 7. line/cylinder and circle/plane touch refinement.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::{Curve, Surface};

use crate::bean_face::BeanFaceIntersector;
use crate::brep_surface::{is_planar, sphere_center, surface_closest_params, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FaceState, FClass2d};
use crate::intcurvesurface::perform_curve_surface;
use crate::inttools_data::{CommonPartType, CommonPrt, IntRange};
use crate::shape::{Edge, Face};

const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Curve classification helpers
// ---------------------------------------------------------------------------

/// Coarse analytic kind of a 3D curve, mirroring the `GeomAbs_CurveType` subset
/// the algorithm branches on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurveKind {
    Line,
    Circle,
    BSpline,
    Other,
}

/// Whether the curve is geometrically a straight line (unbounded, or all
/// samples collinear with the first–last chord).
fn is_line_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return true;
    }
    if (b - a).abs() <= 1e-15 {
        return false;
    }
    let p0 = c.d0(a);
    let pl = c.d0(b);
    let size = p0.distance(&pl);
    if size <= 1e-30 {
        return false;
    }
    let d0 = GpVec::from_pnts(&p0, &pl);
    let tol = 1e-6 * size;
    for i in 1..8 {
        let p = c.d0(a + (b - a) * i as f64 / 8.0);
        if GpVec::from_pnts(&p0, &p).crossed(&d0).magnitude() > tol * size {
            return false;
        }
    }
    true
}

/// First three non-collinear samples of a point set, if they exist.
fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
    let p0 = pts[0];
    let mut i1 = None;
    for (i, p) in pts.iter().enumerate().skip(1) {
        if GpVec::from_pnts(&p0, p).magnitude() > 1e-9 {
            i1 = Some(i);
            break;
        }
    }
    let i1 = i1?;
    let p1 = pts[i1];
    let d0 = GpVec::from_pnts(&p0, &p1);
    for p in pts.iter().skip(i1 + 1) {
        if GpVec::from_pnts(&p0, p).crossed(&d0).magnitude() > 1e-9 * d0.magnitude().max(1e-9) {
            return Some((p0, p1, *p));
        }
    }
    None
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-20 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear 3D points, if it exists.
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.crossed(&d2);
    if n.magnitude() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.xyz().x, n.xyz().y, n.xyz().z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n.xyz())];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

/// Whether the curve is geometrically a (planar) circle: every sample is
/// coplanar and equidistant from a common center.
fn is_circle_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let (p0, p1, p2) = match first_three_spanning(&pts) {
        Some(x) => x,
        None => return false,
    };
    let center = match circumcenter(&p0, &p1, &p2) {
        Some(c) => c,
        None => return false,
    };
    let radius = p0.distance(&center);
    if radius <= 1e-30 {
        return false;
    }
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return false;
    }
    let nv = nrm.divided(m);
    let scale = radius.max(1.0);
    let tol = 1e-6 * scale;
    for p in pts {
        let v = GpVec::from_pnts(&center, &p);
        if v.dot(&nv).abs() > tol {
            return false;
        }
        if (v.magnitude() - radius).abs() > tol {
            return false;
        }
    }
    true
}

/// `(center, radius, unit plane normal)` of a circle-like curve.
fn circle_geometry(c: &dyn Curve) -> Option<(GpPnt, f64, GpVec)> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let (p0, p1, p2) = first_three_spanning(&pts)?;
    let center = circumcenter(&p0, &p1, &p2)?;
    let radius = p0.distance(&center);
    if radius <= 1e-30 {
        return None;
    }
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return None;
    }
    Some((center, radius, nrm.divided(m)))
}

/// Classify a 3D curve by sampling geometric invariants.
fn curve_kind(c: &dyn Curve) -> CurveKind {
    if is_line_like(c) {
        return CurveKind::Line;
    }
    if is_circle_like(c) {
        return CurveKind::Circle;
    }
    if c.is_periodic() && (c.period() - 2.0 * PI).abs() < 1e-9 {
        // Ellipse / other periodic conic: analytic (not B-spline).
        CurveKind::Other
    } else {
        CurveKind::BSpline
    }
}

/// Approximate `BRepAdaptor_Curve::Resolution(tol)`: the parameter increment
/// over which the curve moves at most ~`tol`. Used for the "whole range" tests.
fn curve_resolution(c: &dyn Curve, tol: f64) -> f64 {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return tol.max(PCONFUSION);
    }
    let span = b - a;
    if span.abs() <= 1e-15 {
        return 1.0;
    }
    let n = 64;
    let h = span / n as f64;
    let mut max_speed = 0.0;
    let mut prev = c.d0(a);
    for i in 1..=n {
        let u = a + h * i as f64;
        let p = c.d0(u);
        let d = p.distance(&prev);
        let speed = d / h;
        if speed > max_speed {
            max_speed = speed;
        }
        prev = p;
    }
    if max_speed <= 1e-30 {
        return tol.max(PCONFUSION);
    }
    (tol / max_speed).max(PCONFUSION)
}

// ---------------------------------------------------------------------------
// Analytic surface geometry recovery (dyn Surface cannot be downcast)
// ---------------------------------------------------------------------------

fn surface_sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

fn quadric_u0(s: &dyn Surface) -> f64 {
    let (u0, _) = s.u_range();
    if u0.is_finite() {
        u0
    } else {
        0.0
    }
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// `(origin, unit-u, unit-v)` frame of a planar surface (natural parameterization).
fn plane_frame(s: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let (o, x, y) = s.d1(0.0, 0.0);
    let mx = x.magnitude();
    let my = y.magnitude();
    if mx < 1e-30 || my < 1e-30 {
        return None;
    }
    Some((o, x.divided(mx), y.divided(my)))
}

/// Exact projection of `p` onto a planar surface: `(u, v, distance)`.
fn plane_projection(s: &dyn Surface, p: &GpPnt) -> (f64, f64, f64) {
    let (o, x, y) = match plane_frame(s) {
        Some(f) => f,
        None => {
            let (u, v) = surface_closest_params(s, p, 24, 24);
            let q = s.d0(u, v);
            return (u, v, p.distance(&q));
        }
    };
    let n = x.crossed(&y);
    let m = n.magnitude();
    if m < 1e-30 {
        let (u, v) = surface_closest_params(s, p, 24, 24);
        let q = s.d0(u, v);
        return (u, v, p.distance(&q));
    }
    let n = n.divided(m);
    let d = GpVec::from_pnts(&o, p);
    let dist = d.dot(&n).abs();
    let proj = p.translated_vec(&n.multiplied_scalar(-d.dot(&n)));
    let w = GpVec::from_pnts(&o, &proj);
    (w.dot(&x), w.dot(&y), dist)
}

/// Plane normal of a planar surface.
fn plane_normal(s: &dyn Surface) -> GpVec {
    match plane_frame(s) {
        Some((_, x, y)) => x.crossed(&y).normalized(),
        None => GpVec::zero(),
    }
}

/// `(axis point, axis direction, radius)` of a cylindrical surface.
fn cylinder_geometry(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let u0 = quadric_u0(s);
    let (_, _, v0, _) = surface_sample_bounds(s);
    let c0 = midpoint(&s.d0(u0, v0), &s.d0(u0 + PI, v0));
    let c1 = midpoint(&s.d0(u0, v0 + 1.0), &s.d0(u0 + PI, v0 + 1.0));
    let zvec = GpVec::from_pnts(&c0, &c1);
    let zm = zvec.magnitude();
    if zm < 1e-12 {
        return None;
    }
    let z = zvec.divided(zm);
    let r = s.d0(u0, v0).distance(&c0);
    if r < 1e-12 {
        return None;
    }
    Some((c0, z, r))
}

fn cylinder_matches(s: &dyn Surface, a: &GpPnt, z: &GpVec, r: f64) -> bool {
    let (u0, u1, v0, v1) = surface_sample_bounds(s);
    let tol = 1e-4 * r.abs().max(1.0);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let d = GpVec::from_pnts(a, &s.d0(u, v));
            let rho = d.coord.subtracted(&z.xyz().multiplied(d.dot(z))).modulus();
            if (rho - r).abs() > tol {
                return false;
            }
        }
    }
    true
}

/// `(apex, axis direction, tan(semi-angle))` of a conical surface.
fn cone_geometry(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let u0 = quadric_u0(s);
    let (_, _, v0, _) = surface_sample_bounds(s);
    let c0 = midpoint(&s.d0(u0, v0), &s.d0(u0 + PI, v0));
    let c1 = midpoint(&s.d0(u0, v0 + 1.0), &s.d0(u0 + PI, v0 + 1.0));
    let zvec = GpVec::from_pnts(&c0, &c1);
    let zm = zvec.magnitude();
    if zm < 1e-12 {
        return None;
    }
    let z = zvec.divided(zm);
    let r0 = s.d0(u0, v0).distance(&c0);
    let r1 = s.d0(u0, v0 + 1.0).distance(&c1);
    let tan_alpha = (r1 - r0) / zm;
    if tan_alpha.abs() < 1e-12 {
        return None;
    }
    let s0 = r0 / tan_alpha;
    let apex = c0.translated_vec(&z.multiplied_scalar(-s0));
    Some((apex, z, tan_alpha))
}

fn cone_matches(s: &dyn Surface, apex: &GpPnt, z: &GpVec, tan_alpha: f64) -> bool {
    let cosa = 1.0 / (1.0 + tan_alpha * tan_alpha).sqrt();
    let (u0, u1, v0, v1) = surface_sample_bounds(s);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let d = GpVec::from_pnts(apex, &s.d0(u, v));
            let err = (d.dot(z) * d.dot(z) - d.square_magnitude() * cosa * cosa).abs();
            if err > 1e-4 * d.square_magnitude().max(1.0) {
                return false;
            }
        }
    }
    true
}

/// `(center, axis direction, major radius, minor radius)` of a torus.
fn torus_geometry(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64, f64)> {
    let u0 = quadric_u0(s);
    let o = midpoint(&s.d0(u0, 0.0), &s.d0(u0 + PI, 0.0));
    let r_outer = s.d0(u0, 0.0).distance(&o);
    let o2 = midpoint(&s.d0(u0, PI), &s.d0(u0 + PI, PI));
    let r_inner = s.d0(u0, PI).distance(&o2);
    if o.distance(&o2) > 1e-6 * r_outer.max(1.0) {
        return None;
    }
    let major = 0.5 * (r_outer + r_inner);
    let minor = 0.5 * (r_outer - r_inner);
    if major <= 1e-12 || minor <= 1e-12 {
        return None;
    }
    let x = GpVec::from_pnts(&o, &s.d0(u0, 0.0)).divided(r_outer);
    let y = GpVec::from_pnts(&o, &s.d0(u0 + PI / 2.0, 0.0)).divided(r_outer);
    let z = x.crossed(&y).normalized();
    Some((o, z, major, minor))
}

fn torus_matches(s: &dyn Surface, o: &GpPnt, z: &GpVec, major: f64, minor: f64) -> bool {
    let (u0, u1, v0, v1) = surface_sample_bounds(s);
    let tol = 1e-4 * major.max(minor).max(1.0);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let d = GpVec::from_pnts(o, &s.d0(u, v));
            let dz = d.dot(z);
            let rho = d.coord.subtracted(&z.xyz().multiplied(dz)).modulus();
            let err = ((rho - major).powi(2) + dz * dz).sqrt() - minor;
            if err.abs() > tol {
                return false;
            }
        }
    }
    true
}

/// Full analytic surface classification (plane, sphere, cylinder, cone, torus).
fn surface_kind(s: &dyn Surface) -> SurfaceKind {
    if is_planar(s, 6, 6, 1e-6) {
        return SurfaceKind::Plane;
    }
    if let Some((a, z, r)) = cylinder_geometry(s) {
        if cylinder_matches(s, &a, &z, r) {
            return SurfaceKind::Cylinder;
        }
    }
    if let Some((apex, z, ta)) = cone_geometry(s) {
        if cone_matches(s, &apex, &z, ta) {
            return SurfaceKind::Cone;
        }
    }
    if let Some((o, z, maj, min)) = torus_geometry(s) {
        if torus_matches(s, &o, &z, maj, min) {
            return SurfaceKind::Torus;
        }
    }
    if sphere_center(s).is_some() {
        return SurfaceKind::Sphere;
    }
    SurfaceKind::Other
}

/// Shortcut distance for a point lying on the axis of an analytic surface.
///
/// Port of `IntTools_EdgeFace::IsEqDistance`: when the point is within `TOL` of
/// the cylinder axis / cone axis / torus major circle, the surface distance is
/// the analytic radius directly (projection is ill-defined there).
fn is_eq_distance(p: &GpPnt, s: &dyn Surface) -> Option<f64> {
    const TOL: f64 = 1e-7;
    match surface_kind(s) {
        SurfaceKind::Cylinder => {
            let (a, z, r) = cylinder_geometry(s)?;
            let dc = GpVec::from_pnts(&a, p).crossed(&z).magnitude();
            if dc < TOL {
                Some(r)
            } else {
                None
            }
        }
        SurfaceKind::Cone => {
            let (apex, z, ta) = cone_geometry(s)?;
            let dc = GpVec::from_pnts(&apex, p).crossed(&z).magnitude();
            if dc < TOL {
                Some(p.distance(&apex) * ta)
            } else {
                None
            }
        }
        SurfaceKind::Torus => {
            let (o, _, maj, min) = torus_geometry(s)?;
            let dc = (o.distance(p) - maj).abs();
            if dc < TOL {
                Some(min)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// `IntTools_Tools::IsDirsCoinside`: unit directions equal or opposite within
/// `0.0002` (measured as the distance between the unit-sphere points).
fn dirs_coinside(d1: &GpVec, d2: &GpVec) -> bool {
    let d = d1.subtracted(d2).magnitude();
    d < 0.0002 || (2.0 - d).abs() < 0.0002
}

/// Whether a circle-like curve is coplanar with a planar surface.
fn is_coplanar(curve: &dyn Curve, surface: &dyn Surface) -> bool {
    if !is_circle_like(curve) || !is_planar(surface, 6, 6, 1e-6) {
        return false;
    }
    let (_, _, n) = match circle_geometry(curve) {
        Some(g) => g,
        None => return false,
    };
    dirs_coinside(&n, &plane_normal(surface))
}

/// Whether a circle-like curve is tangent to a planar surface (its center is at
/// distance ≈ radius from the plane).
fn is_radius(curve: &dyn Curve, surface: &dyn Surface, criteria: f64) -> bool {
    if !is_circle_like(curve) || !is_planar(surface, 6, 6, 1e-6) {
        return false;
    }
    let (center, r, _) = match circle_geometry(curve) {
        Some(g) => g,
        None => return false,
    };
    let pn = plane_normal(surface);
    let po = surface.d0(0.0, 0.0);
    let d = GpVec::from_pnts(&po, &center).dot(&pn).abs();
    (d - r).abs() < criteria
}

// ---------------------------------------------------------------------------
// EdgeFace
// ---------------------------------------------------------------------------

/// Edge/face intersection algorithm. Port of `IntTools_EdgeFace`.
///
/// Not `Debug`-derived: the cached `Arc<dyn Curve>`/`Arc<dyn Surface>` handles
/// are not `Debug`. `Clone` is provided (Arc handles clone cheaply).
#[derive(Clone)]
pub struct EdgeFace {
    edge: Edge,
    face: Face,
    range: IntRange,
    fuzzy_value: f64,
    quick_coincidence_check: bool,
    curve: Option<Arc<dyn Curve>>,
    surface: Option<Arc<dyn Surface>>,
    criteria: f64,
    is_done: bool,
    error_status: i32,
    common_parts: Vec<CommonPrt>,
    face_classifier: Option<FClass2d>,
}

impl Default for EdgeFace {
    fn default() -> Self {
        Self::new()
    }
}

impl EdgeFace {
    /// Empty constructor (`IntTools_EdgeFace()`).
    pub fn new() -> Self {
        Self {
            edge: Edge::new(),
            face: Face::new(),
            range: IntRange::new_unchecked(f64::NEG_INFINITY, f64::INFINITY),
            fuzzy_value: CONFUSION,
            quick_coincidence_check: false,
            curve: None,
            surface: None,
            criteria: CONFUSION,
            is_done: false,
            error_status: 1,
            common_parts: Vec::new(),
            face_classifier: None,
        }
    }

    // ---- setters / getters --------------------------------------------------

    /// Sets the edge for intersection.
    pub fn set_edge(&mut self, edge: Edge) {
        self.edge = edge;
    }

    /// Returns the edge.
    pub fn edge(&self) -> &Edge {
        &self.edge
    }

    /// Sets the face for intersection.
    pub fn set_face(&mut self, face: Face) {
        self.face = face;
    }

    /// Returns the face.
    pub fn face(&self) -> &Face {
        &self.face
    }

    /// Sets the boundaries of the edge to process.
    pub fn set_range(&mut self, first: f64, last: f64) {
        self.range = IntRange::new_unchecked(first, last);
    }

    /// Returns the processing range.
    pub fn range(&self) -> IntRange {
        self.range
    }

    /// Sets the fuzzy value (clamped to `Precision::Confusion`).
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.fuzzy_value = v.max(CONFUSION);
    }

    /// Returns the fuzzy value.
    pub fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value
    }

    /// Sets the quick coincidence check flag.
    pub fn set_quick_coincidence_check(&mut self, b: bool) {
        self.quick_coincidence_check = b;
    }

    /// Returns the quick coincidence check flag.
    pub fn is_coincidence_checked_quickly(&self) -> bool {
        self.quick_coincidence_check
    }

    // ---- performing ----------------------------------------------------------

    /// Launches the intersection.
    ///
    /// `Err` is returned only for hard setup failures (missing geometry, invalid
    /// range). Algorithmic failures set [`error_status`](Self::error_status)
    /// (`2`/`3`/`4` per OCCT) and still return `Ok`.
    pub fn perform(&mut self) -> Result<(), String> {
        self.common_parts.clear();
        self.error_status = 0;
        self.check_data();
        if self.error_status != 0 {
            return Ok(());
        }

        let curve = BRepTool::edge_curve(&self.edge)
            .ok_or_else(|| "EdgeFace::perform: edge has no curve".to_string())?;
        let surface = BRepTool::face_surface(&self.face)
            .ok_or_else(|| "EdgeFace::perform: face has no surface".to_string())?;
        let (ef, el) = BRepTool::edge_parameters(&self.edge);
        if !ef.is_finite() || !el.is_finite() || el - ef <= 1e-15 {
            return Err("EdgeFace::perform: edge has an empty or unbounded parameter range".into());
        }
        // Default the processing range to the whole edge when not explicitly set.
        if !self.range.is_valid() || !self.range.first.is_finite() || !self.range.last.is_finite() {
            self.range = IntRange::new_unchecked(ef, el);
        }

        self.is_done = false;
        self.curve = Some(curve.clone());
        self.surface = Some(surface.clone());

        let c_kind = curve_kind(curve.as_ref());
        let s_kind = surface_kind(surface.as_ref());

        // Prepare myCriteria.
        let fuzz = self.fuzzy_value * 0.5;
        let tol_f = BRepTool::face_tolerance(&self.face) + fuzz;
        let tol_e = BRepTool::edge_tolerance(&self.edge) + fuzz;
        self.criteria = match c_kind {
            CurveKind::BSpline => {
                let diff1 = tol_e / tol_f.max(1e-300);
                let diff2 = tol_f / tol_e.max(1e-300);
                if diff1 > 100.0 || diff2 > 100.0 {
                    tol_e.max(tol_f)
                } else {
                    1.5 * tol_e + tol_f
                }
            }
            _ => tol_e + tol_f,
        };

        // 2D classifier for the face's UV restriction (used by coincidence and
        // projectability checks).
        let cl_tol = BRepTool::face_tolerance(&self.face).max(PCONFUSION);
        self.face_classifier = Some(FClass2d::new(&self.face, cl_tol)?);

        if self.quick_coincidence_check && self.is_coincident() {
            let mut cp = CommonPrt::new();
            cp.part_type = CommonPartType::Edge;
            cp.range = self.range;
            cp.face = Some(self.face.0.clone());
            self.common_parts.push(cp);
            self.is_done = true;
            return Ok(());
        }

        let mut intersector = BeanFaceIntersector::new();
        intersector.initialize(curve.clone(), surface.clone(), tol_e, tol_f);
        intersector.set_bean_parameters(self.range.first, self.range.last);
        // The surface parameter window must be the face's *trimmed* UV bounds
        // (from the boundary wires), not the unbounded surface range — otherwise
        // a coplanar edge extending past the face is reported as fully on-face.
        let (u0, u1, v0, v1) = crate::wireframe::face_uv_bounds(&self.face, surface.as_ref());
        intersector.set_surface_parameters(u0, u1, v0, v1);
        intersector.perform()?;
        if !intersector.is_done() {
            return Ok(());
        }

        for r in intersector.result() {
            let mid = 0.5 * (r.first + r.last);
            if self.is_projectable(mid) {
                let mut cp = CommonPrt::new();
                cp.range = r;
                cp.face = Some(self.face.0.clone());
                self.common_parts.push(cp);
            }
        }

        let nb = self.common_parts.len();
        for i in 0..nb {
            let mut cp = self.common_parts[i].clone();
            self.make_type(&mut cp);
            self.common_parts[i] = cp;
        }

        // Line/Cylinder and Circle/Plane common-part refinement.
        let special = (c_kind == CurveKind::Line && s_kind == SurfaceKind::Cylinder)
            || (c_kind == CurveKind::Circle
                && s_kind == SurfaceKind::Plane
                && !is_coplanar(curve.as_ref(), surface.as_ref())
                && !is_radius(curve.as_ref(), surface.as_ref(), self.criteria));
        if special {
            self.refine_touch_parts();
        }

        self.is_done = true;
        Ok(())
    }

    /// Whether the computation was successful.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// Completion code: `0` success; `1` not started; `2`/`3` invalid input;
    /// `4` projection failed.
    pub fn error_status(&self) -> i32 {
        self.error_status
    }

    /// The resulting common parts.
    pub fn common_parts(&self) -> &[CommonPrt] {
        &self.common_parts
    }

    /// For each common part, `Some(t)` when the part is a single point at edge
    /// parameter `t`, `None` when it is a coincident edge sub-range.
    pub fn point_parameters(&self) -> Vec<Option<f64>> {
        self.common_parts
            .iter()
            .map(|cp| {
                if cp.range.length() <= self.criteria * 2.0 {
                    Some(cp.range.first)
                } else {
                    None
                }
            })
            .collect()
    }

    // ---- internals -----------------------------------------------------------

    /// `CheckData`: sets `error_status` to `2` (degenerated edge) or `3`
    /// (non-geometric edge).
    fn check_data(&mut self) {
        if BRepTool::is_degenerated(&self.edge) {
            self.error_status = 2;
        }
        if BRepTool::edge_curve(&self.edge).is_none() {
            self.error_status = 3;
        }
    }

    /// Whether the edge is entirely on the face, sampled over the range.
    /// Port of `IntTools_EdgeFace::IsCoincident`.
    fn is_coincident(&self) -> bool {
        let curve = self.curve.clone().expect("curve set");
        let surface = self.surface.clone().expect("surface set");
        let a_nb_seg = if curve_kind(curve.as_ref()) == CurveKind::Line
            && surface_kind(surface.as_ref()) == SurfaceKind::Plane
        {
            2
        } else {
            23
        };
        let a_tresh = 0.5;
        let a_tresh_idx_f = ((a_nb_seg + 1) as f64 * 0.25) as i32;
        let a_tresh_idx_l = ((a_nb_seg + 1) as f64 * 0.75) as i32;

        let (mut a_t1, mut a_t2) = (self.range.first, self.range.last);
        if a_t2 - a_t1 <= 1e-12 {
            return false;
        }
        let a_bnd_shift = 0.01 * (a_t2 - a_t1);
        a_t1 += a_bnd_shift;
        a_t2 -= a_bnd_shift;
        if a_t2 <= a_t1 {
            return false;
        }
        let d_t = (a_t2 - a_t1) / a_nb_seg as f64;

        let mut is_classified = false;
        let mut i_cnt = 0;
        for i in 0..=a_nb_seg {
            let a_t = a_t1 + i as f64 * d_t;
            let a_p = curve.d0(a_t);
            let (u, v, a_d) = self.project_point(&a_p);
            if a_d > self.criteria {
                if a_d > 100.0 * self.criteria {
                    return false;
                }
                continue;
            }
            i_cnt += 1;
            if ((0 < i) && (i < a_tresh_idx_f)) || ((a_tresh_idx_l < i) && (i < a_nb_seg)) {
                continue;
            }
            if is_classified && (i != a_nb_seg) {
                continue;
            }
            let state = self
                .face_classifier
                .as_ref()
                .map(|c| c.perform(GpPnt2d::new(u, v)))
                .unwrap_or(FaceState::In);
            if state == FaceState::Out {
                return false;
            }
            if i != 0 {
                is_classified = true;
            }
        }
        let a_coeff = i_cnt as f64 / (a_nb_seg + 1) as f64;
        a_coeff > a_tresh
    }

    /// Whether the curve point at parameter `t` is on the face within
    /// `myCriteria` and inside its 2D restriction.
    /// Port of `IntTools_EdgeFace::IsProjectable`.
    fn is_projectable(&self, t: f64) -> bool {
        let curve = self.curve.clone().expect("curve set");
        let p = curve.d0(t);
        let (u, v, dist) = self.project_point(&p);
        if dist > self.criteria {
            return false;
        }
        match &self.face_classifier {
            Some(cl) => cl.perform(GpPnt2d::new(u, v)) != FaceState::Out,
            None => true,
        }
    }

    /// Signed distance from the curve point at `t` to the surface, minus
    /// `myCriteria`. Port of `IntTools_EdgeFace::DistanceFunction`.
    fn distance_function(&self, t: f64) -> f64 {
        let curve = self.curve.clone().expect("curve set");
        let surface = self.surface.clone().expect("surface set");
        let p = curve.d0(t);
        if let Some(d) = is_eq_distance(&p, surface.as_ref()) {
            return d - self.criteria;
        }
        let (_, _, dist) = self.project_point(&p);
        dist - self.criteria
    }

    /// Project a point onto the face surface: `(u, v, distance)`. Planar faces
    /// are projected analytically (exact distance); others use the grid+refine
    /// projector.
    fn project_point(&self, p: &GpPnt) -> (f64, f64, f64) {
        let surface = self.surface.clone().expect("surface set");
        if is_planar(surface.as_ref(), 6, 6, 1e-6) {
            plane_projection(surface.as_ref(), p)
        } else {
            let (u, v) = surface_closest_params(surface.as_ref(), p, 24, 24);
            let q = surface.d0(u, v);
            (u, v, p.distance(&q))
        }
    }

    /// Distance from the curve point at `t` to the face surface (unsigned).
    fn surface_distance_at(&self, t: f64) -> f64 {
        let curve = self.curve.clone().expect("curve set");
        let p = curve.d0(t);
        let (_, _, d) = self.project_point(&p);
        d
    }

    /// Sample the curve–surface distance over `[t0, t1]`. Returns
    /// `(min, max, param_at_min)` with golden-section refinement around the min.
    fn distance_profile(&self, t0: f64, t1: f64, n: usize) -> (f64, f64, f64) {
        let curve = self.curve.clone().expect("curve set");
        if t1 <= t0 {
            let p = curve.d0(t0);
            let (_, _, d) = self.project_point(&p);
            return (d, d, t0);
        }
        let mut min_d = f64::INFINITY;
        let mut max_d = 0.0;
        let mut min_t = t0;
        for i in 0..=n {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            let p = curve.d0(t);
            let (_, _, d) = self.project_point(&p);
            if d < min_d {
                min_d = d;
                min_t = t;
            }
            if d > max_d {
                max_d = d;
            }
        }
        let h = (t1 - t0) / n as f64;
        let lo = (min_t - h).max(t0);
        let hi = (min_t + h).min(t1);
        if hi > lo {
            let (rt, rd) = golden_1d(&|t: f64| self.surface_distance_at(t), lo, hi, 1e-10);
            if rd < min_d {
                min_d = rd;
                min_t = rt;
            }
        }
        (min_d, max_d, min_t)
    }

    /// Classify a range as an edge (coincident sub-range) or a point.
    /// Port of `IntTools_EdgeFace::MakeType`.
    ///
    /// A point collapses the range to `[t, t]` (the touch/mid parameter); an
    /// edge keeps the original non-degenerate range. Both carry
    /// [`CommonPartType::Edge`].
    fn make_type(&mut self, cp: &mut CommonPrt) -> i32 {
        let af1 = cp.range.first;
        let al1 = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let a_pf = curve.d0(af1);
        let a_pl = curve.d0(al1);
        let df1 = a_pf.distance(&a_pl);
        let a_cr = curve_resolution(curve.as_ref(), self.criteria);
        let is_whole_range =
            (af1 - self.range.first).abs() < a_cr && (al1 - self.range.last).abs() < a_cr;

        if df1 > self.criteria * 2.0 {
            // A long common part is an EDGE (coincident on-face range) when its
            // interior lies on the face — the whole edge, or a partial on-face
            // sub-range such as a coplanar edge trimmed to the face boundary.
            // It collapses to a touch POINT only when the middle is off the
            // face (a tangency at the range end).
            let tm = 0.5 * (af1 + al1);
            if is_whole_range || self.is_projectable(tm) {
                cp.part_type = CommonPartType::Edge;
                return 0;
            }
        }

        let mut tm = 0.5 * (af1 + al1);
        if !self.check_touch(cp, &mut tm) {
            tm = 0.5 * (af1 + al1);
        }
        cp.part_type = CommonPartType::Edge;
        cp.range = IntRange::new_unchecked(tm, tm);
        0
    }

    /// Whether the range contains a touch point within `myCriteria`, and the
    /// touch parameter. Port of `IntTools_EdgeFace::CheckTouch`.
    fn check_touch(&self, cp: &CommonPrt, tx: &mut f64) -> bool {
        let a_tf = cp.range.first;
        let a_tl = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let a_cr = curve_resolution(curve.as_ref(), self.criteria);
        if (a_tf - self.range.first).abs() < a_cr && (a_tl - self.range.last).abs() < a_cr {
            return false; // whole range: keep EDGE
        }

        let (min_d, max_d, min_t) = self.distance_profile(a_tf, a_tl, 32);
        // Extrema parallel case: the distance is nearly constant over the range.
        if max_d - min_d <= 0.05 * (min_d + self.criteria.max(1e-9)) {
            return false;
        }
        let mut a_dist2 = min_d * min_d;
        let mut a_tx = min_t;

        // Exact curve–surface intersection fallback (the `Extrema` aNbExt == 0
        // branch of `IntTools_EdgeFace::CheckTouch`): when the sampled profile
        // found no near-surface minimum, consult the exact intersector.
        if a_dist2 > self.criteria * self.criteria && a_tl > a_tf {
            let surface = self.surface.clone().expect("surface set");
            let (u0, u1, v0, v1) = surface_sample_bounds(surface.as_ref());
            if let Ok(hr) = perform_curve_surface(
                curve.as_ref(),
                surface.as_ref(),
                (a_tf, a_tl),
                (u0, u1, v0, v1),
            ) {
                for i in 0..hr.nb_points() {
                    let p = hr.point(i);
                    if p.param() >= a_tf && p.param() <= a_tl {
                        a_dist2 = 0.0;
                        a_tx = p.param();
                        break;
                    }
                }
            }
        }

        let b1 = self.distance_function(a_tf) + self.criteria;
        if b1 * b1 < a_dist2 {
            a_dist2 = b1 * b1;
            a_tx = a_tf;
        }
        let b2 = self.distance_function(a_tl) + self.criteria;
        if b2 * b2 < a_dist2 {
            a_dist2 = b2 * b2;
            a_tx = a_tl;
        }
        let bm = self.distance_function(0.5 * (a_tf + a_tl)) + self.criteria;
        if bm * bm < a_dist2 {
            a_dist2 = bm * bm;
            a_tx = 0.5 * (a_tf + a_tl);
        }

        if a_dist2 > self.criteria * self.criteria {
            return false;
        }
        *tx = a_tx;
        if (a_tx - a_tf).abs() < PCONFUSION {
            return true;
        }
        if (a_tx - a_tl).abs() < PCONFUSION {
            return true;
        }
        if a_tx > a_tf && a_tx < a_tl {
            return true;
        }
        false
    }

    /// Vertex-specific touch refinement. Port of `IntTools_EdgeFace::CheckTouchVertex`.
    fn check_touch_vertex(&self, cp: &CommonPrt, tx: &mut f64) -> bool {
        let a_tf = cp.range.first;
        let a_tl = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let a_type = curve_kind(curve.as_ref());
        let a_eps_t = if a_type == CurveKind::Line { 9e-5 } else { 8e-5 };
        let a_tm = 0.5 * (a_tf + a_tl);
        let a_dist2 = {
            let d = self.distance_function(a_tm);
            d * d
        };
        if a_tl <= a_tf {
            return false;
        }
        let (min_d, max_d, min_t) = self.distance_profile(a_tf, a_tl, 32);
        if max_d - min_d <= 0.05 * (min_d + self.criteria.max(1e-9)) {
            return false;
        }
        let a_dist2_new = min_d * min_d;
        if a_dist2_new > a_dist2 {
            *tx = a_tm;
            return true;
        }
        if a_dist2_new > self.criteria * self.criteria {
            return false;
        }
        let a_tx = min_t;
        if (a_tx - a_tf).abs() < a_eps_t {
            return false;
        }
        if (a_tx - a_tl).abs() < a_eps_t {
            return false;
        }
        if a_tx > a_tf && a_tx < a_tl {
            *tx = a_tx;
            return true;
        }
        false
    }

    /// The line/cylinder and circle/plane special treatment: refine EDGE/VERTEX
    /// common parts into touch points when the range is only tangent.
    fn refine_touch_parts(&mut self) {
        for i in 0..self.common_parts.len() {
            if self.common_parts[i].part_type != CommonPartType::Edge {
                continue;
            }
            let is_point = self.common_parts[i].range.length() <= self.criteria * 2.0;
            let cp = self.common_parts[i].clone();
            let mut tx = 0.0;
            let touched = if is_point {
                self.check_touch_vertex(&cp, &mut tx)
            } else {
                self.check_touch(&cp, &mut tx)
            };
            if touched {
                self.common_parts[i].range = IntRange::new_unchecked(tx, tx);
            }
        }
    }
}

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    const GOLD: f64 = 0.6180339887498949;
    let mut a = lo;
    let mut b = hi;
    let mut c = b - GOLD * (b - a);
    let mut d = a + GOLD * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    while (b - a) > eps {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - GOLD * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + GOLD * (b - a);
            fd = f(d);
        }
    }
    let x = 0.5 * (a + b);
    (x, f(x))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::inttools::edge_face_intersections;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn box_face(idx: usize) -> (crate::brep_extrema::test_box::UnitBox, Face) {
        let ub = unit_box();
        let face = ub.faces[idx].clone();
        (ub, face)
    }

    #[test]
    fn default_state_is_not_started() {
        let ef = EdgeFace::new();
        assert!(!ef.is_done());
        assert_eq!(ef.error_status(), 1);
        assert!(ef.common_parts().is_empty());
        assert_eq!(ef.fuzzy_value(), CONFUSION);
        assert!(!ef.is_coincidence_checked_quickly());
    }

    #[test]
    fn line_in_face_plane_is_edge_part() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Edge lying in the face plane, crossing the face interior.
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.0), &GpPnt::new(1.0, 0.5, 0.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);

        let cps = ef.common_parts();
        assert_eq!(cps.len(), 1, "common parts: {cps:?}");
        assert_eq!(cps[0].part_type, CommonPartType::Edge);
        assert!(cps[0].range.length() > 0.5, "coincident range: {cps:?}");
        let params = ef.point_parameters();
        assert_eq!(params.len(), 1);
        assert!(params[0].is_none(), "in-plane edge is not a point: {params:?}");

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn line_piercing_face_is_point_and_matches_approximate() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 2.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);

        let cps = ef.common_parts();
        assert_eq!(cps.len(), 1, "common parts: {cps:?}");
        let params = ef.point_parameters();
        assert_eq!(params.len(), 1);
        let t = params[0].expect("piercing line is a point");
        let curve = BRepTool::edge_curve(&e).expect("edge curve");
        let p = curve.d0(t);
        assert!(
            p.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-4,
            "intersection point {p:?} at t={t}"
        );

        // Cross-check with the existing approximate solver.
        let hits = edge_face_intersections(&e, &face, 1e-9);
        assert_eq!(hits.len(), 1, "approximate hits: {hits:?}");
        assert!(p.distance(&hits[0].1) < 1e-4, "{p:?} vs {:?}", hits[0].1);

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn line_outside_face_is_empty() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Pierces the plane but outside the face's UV domain.
        let e = b.make_edge_segment(&GpPnt::new(2.0, 2.0, -1.0), &GpPnt::new(2.0, 2.0, 1.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 2.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);
        assert!(
            ef.common_parts().is_empty(),
            "common parts: {:?}",
            ef.common_parts()
        );

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn edge_endpoint_on_face_is_point() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Starts exactly on the face, then leaves along +Z.
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, 0.0), &GpPnt::new(0.5, 0.5, 1.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        assert_eq!(ef.error_status(), 0);

        let params = ef.point_parameters();
        assert_eq!(params.len(), 1, "parts: {:?}", ef.common_parts());
        let t = params[0].expect("endpoint touch is a point");
        assert!(t.abs() < 1e-4, "touch parameter {t}, want ~0");

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn quick_coincidence_detects_in_plane_edge() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.0), &GpPnt::new(1.0, 0.5, 0.0));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.set_quick_coincidence_check(true);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        let cps = ef.common_parts();
        assert_eq!(cps.len(), 1, "common parts: {cps:?}");
        assert_eq!(cps[0].part_type, CommonPartType::Edge);
        assert!(cps[0].range.length() > 0.5);

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn quick_coincidence_rejects_offset_line() {
        let (ub, face) = box_face(0); // bottom, z = 0
        let b = TopoBuilder::new();
        // Parallel to the plane, 0.5 above it — not coincident.
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.5), &GpPnt::new(1.0, 0.5, 0.5));

        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.set_quick_coincidence_check(true);
        ef.perform().expect("perform succeeds");
        assert!(ef.is_done());
        // Offset line never touches the surface: no common parts.
        assert!(
            ef.common_parts().is_empty(),
            "common parts: {:?}",
            ef.common_parts()
        );

        clear_tree(&e.0);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn non_geometric_edge_sets_error_status() {
        let (ub, face) = box_face(0);
        // An edge with no registered curve is "non-geometric" → error status 3.
        let e = Edge::new();
        let mut ef = EdgeFace::new();
        ef.set_edge(e.clone());
        ef.set_face(face.clone());
        ef.set_range(0.0, 1.0);
        ef.perform().expect("perform does not hard-fail");
        assert!(!ef.is_done());
        assert_eq!(ef.error_status(), 3);
        clear_tree(&ub.solid.0);
    }

    #[test]
    fn fuzzy_value_clamped_to_confusion() {
        let mut ef = EdgeFace::new();
        ef.set_fuzzy_value(0.0);
        assert_eq!(ef.fuzzy_value(), CONFUSION);
        ef.set_fuzzy_value(1e-3);
        assert_eq!(ef.fuzzy_value(), 1e-3);
    }
}
