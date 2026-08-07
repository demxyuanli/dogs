//! Exact curve–surface intersection.
//!
//! Port of `IntCurveSurface_HInter` (TKGeomAlgo) at the level the precise
//! boolean needs: given a parametric curve and a surface, compute every point
//! (and, for the coincidence degeneracy, segment) where the curve meets the
//! surface.
//!
//! Dispatch mirrors the OCCT `PerformBounds` → `PerformConicSurf` /
//! `InternalPerform` split:
//!
//! - **Analytic path** for conic curves against quadric surfaces: a line
//!   against a plane/sphere/cylinder/cone (substitute the line's affine
//!   parametrization into the surface's implicit equation → linear/quadratic
//!   polynomial), a line against a torus (`IntAna_IntLinTorus` quartic), and a
//!   circle/ellipse against a plane or sphere (`a·cos + b·sin = c` via
//!   `math_TrigonometricFunctionRoots`). Tangency / coincidence produces a
//!   zero-length tangent point or an `On` segment.
//! - **General path** for everything else (B-spline curves, non-quadric
//!   surfaces): dense curve sampling, signed-distance sign-change bisection for
//!   transverse crossings, and golden-section refinement of near-zero distance
//!   minima for tangencies.
//!
//! Trait objects (`Arc<dyn Curve>` / `Arc<dyn Surface>`) cannot be downcast, so
//! analytic dispatch classifies curves and surfaces by geometric invariants —
//! the established pattern in this port (see `brep_surface::classify_surface`).

use std::cmp::Ordering;

use occt_core::gp::{GpAx3, GpDir, GpLin, GpPnt, GpTorus, GpVec};
use occt_geom::intana::line_torus_intersect;
use occt_geom::{Curve, Surface};

use crate::brep_surface::{is_planar, sphere_center, surface_closest_params, surface_normal};

/// Tolerance for the analytic path (exact solves).
const ANALYTIC_TOL: f64 = 1e-9;
/// Tolerance for the general sampling path.
const GENERAL_TOL: f64 = 1e-6;
/// Parameter-merge window: points closer than this are merged as a tangency.
const MERGE_TOL: f64 = 1e-7;
/// Angular tolerance for transition classification (mirrors `Precision::Angular`).
const ANG_TOL: f64 = 1e-12;

// ---------------------------------------------------------------------------
// Data classes
// ---------------------------------------------------------------------------

/// Transition of the curve across the surface at an intersection point.
/// Mirrors `TopAbs_State`: `In` = entering the material, `Out` = leaving,
/// `On` = tangent (the curve only touches the surface), `Unknown` = undecided.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    In,
    Out,
    On,
    Unknown,
}

/// A point common to the curve and the surface.
///
/// Mirrors `IntCurveSurface_IntersectionPoint`: the parameter on the curve
/// (`param`, OCCT `W`), the surface parameters (`u`, `v`), the 3D point and
/// the curve transition at the crossing.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionPoint {
    pub param: f64,
    pub u: f64,
    pub v: f64,
    pub pnt: GpPnt,
    pub state: State,
}

impl IntersectionPoint {
    pub fn new(param: f64, u: f64, v: f64, pnt: GpPnt, state: State) -> Self {
        Self { param, u, v, pnt, state }
    }
    /// Parameter on the curve.
    pub fn param(&self) -> f64 {
        self.param
    }
    /// U parameter on the surface.
    pub fn u(&self) -> f64 {
        self.u
    }
    /// V parameter on the surface.
    pub fn v(&self) -> f64 {
        self.v
    }
    /// The 3D intersection point.
    pub fn pnt(&self) -> GpPnt {
        self.pnt
    }
    /// Transition of the curve at the point.
    pub fn state(&self) -> State {
        self.state
    }
}

/// A curve span `[first, last]` lying on the surface (the coincidence case,
/// e.g. a line in a plane or a circle on a sphere).
///
/// Mirrors `IntCurveSurface_IntersectionSegment`.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionSegment {
    first: IntersectionPoint,
    last: IntersectionPoint,
}

impl IntersectionSegment {
    pub fn new(first: IntersectionPoint, last: IntersectionPoint) -> Self {
        Self { first, last }
    }
    pub fn first_point(&self) -> &IntersectionPoint {
        &self.first
    }
    pub fn second_point(&self) -> &IntersectionPoint {
        &self.last
    }
}

/// Result of a curve–surface intersection: discrete points plus (degenerate)
/// coincident segments.
#[derive(Debug, Clone, Default)]
pub struct HInterResult {
    points: Vec<IntersectionPoint>,
    segments: Vec<IntersectionSegment>,
}

impl HInterResult {
    /// Number of isolated intersection points.
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }
    /// The `i`-th intersection point (0-based).
    pub fn point(&self, i: usize) -> &IntersectionPoint {
        &self.points[i]
    }
    /// All intersection points.
    pub fn points(&self) -> &[IntersectionPoint] {
        &self.points
    }
    /// Number of coincidence segments.
    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }
    /// The `i`-th coincidence segment (0-based).
    pub fn segment(&self, i: usize) -> &IntersectionSegment {
        &self.segments[i]
    }
    /// All coincidence segments.
    pub fn segments(&self) -> &[IntersectionSegment] {
        &self.segments
    }
}

// ---------------------------------------------------------------------------
// Curve classification
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurveKind {
    Line,
    Circle,
    Ellipse,
    Other,
}

fn classify_curve(c: &dyn Curve) -> CurveKind {
    if is_line_like(c) {
        return CurveKind::Line;
    }
    if c.is_periodic() && (c.period() - 2.0 * std::f64::consts::PI).abs() < 1e-9 {
        if is_circle_like(c) {
            CurveKind::Circle
        } else {
            CurveKind::Ellipse
        }
    } else {
        CurveKind::Other
    }
}

/// Whether the curve is geometrically a straight line: either unbounded (the
/// port's `GeomLine`) or with collinear samples.
fn is_line_like(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() && (b - a).abs() <= 1e-15 {
        return false;
    }
    // Choose a sampling window: the natural range, or [-2, 2] for unbounded
    // curves (enough to expose a parabola/hyperbola's curvature).
    let (lo, hi) = if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        (-2.0, 2.0)
    };
    let p0 = c.d0(lo);
    let pl = c.d0(hi);
    let size = p0.distance(&pl);
    if size <= 1e-30 {
        return false;
    }
    let d0 = GpVec::from_pnts(&p0, &pl);
    let tol = 1e-7 * size;
    for i in 1..8 {
        let u = lo + (hi - lo) * i as f64 / 8.0;
        let p = c.d0(u);
        if GpVec::from_pnts(&p0, &p).crossed(&d0).magnitude() > tol * size {
            return false;
        }
    }
    true
}

/// Whether the curve is a circle: all samples are coplanar and equidistant
/// from one center.
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
    let nvec = nrm.divided(m);
    let scale = radius.max(1.0);
    let tol = 1e-6 * scale;
    for p in &pts {
        let v = GpVec::from_pnts(&center, p);
        if v.dot(&nvec).abs() > tol {
            return false;
        }
        if (v.magnitude() - radius).abs() > tol {
            return false;
        }
    }
    true
}

/// First three non-collinear points of a sample set.
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

/// Circumcenter of three non-collinear points (perpendicular-bisector system).
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

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// Reconstructed conic geometry. For a natural-parameter periodic 2π curve,
/// `d0(u) = o + x·r·cos(u) + y·r·sin(u)` (circle) or
/// `d0(u) = o + x·a·cos(u) − y·b·sin(u)` (ellipse).
enum ConicGeom {
    Circle { o: GpPnt, r: f64, x: GpVec, y: GpVec },
    Ellipse { o: GpPnt, a: f64, b: f64, x: GpVec, y: GpVec },
}

/// Reconstruct a natural-parameter conic from the curve's `d0`, verifying the
/// parameterization so a non-natural curve (e.g. a trimmed reparametrization)
/// falls back to the general path.
fn conic_geometry(c: &dyn Curve, kind: CurveKind) -> Option<ConicGeom> {
    match kind {
        CurveKind::Circle => {
            let p0 = c.d0(0.0);
            let p1 = c.d0(std::f64::consts::FRAC_PI_2);
            let p2 = c.d0(std::f64::consts::PI);
            let o = circumcenter(&p0, &p1, &p2)?;
            let r = p0.distance(&o);
            if r <= 1e-30 {
                return None;
            }
            let x = GpVec::from_pnts(&o, &p0).divided(r);
            let y = GpVec::from_pnts(&o, &p1).divided(r);
            // Verify natural parameterization at a few samples.
            for &u in &[0.3f64, 1.1f64, 5.0f64] {
                let expect = o.translated_vec(&x.multiplied_scalar(r * u.cos()))
                    .translated_vec(&y.multiplied_scalar(r * u.sin()));
                if c.d0(u).distance(&expect) > 1e-6 * r.max(1.0) {
                    return None;
                }
            }
            Some(ConicGeom::Circle { o, r, x, y })
        }
        CurveKind::Ellipse => {
            let p0 = c.d0(0.0);
            let pm = c.d0(std::f64::consts::PI);
            let pq = c.d0(std::f64::consts::FRAC_PI_2);
            let o = midpoint(&p0, &pm);
            let a = p0.distance(&o);
            let b = pq.distance(&o);
            if a <= 1e-30 || b <= 1e-30 {
                return None;
            }
            let x = GpVec::from_pnts(&o, &p0).divided(a);
            let y = GpVec::from_pnts(&pq, &o).divided(b); // P(π/2) = o − y·b
            for &u in &[0.3f64, 1.1f64, 5.0f64] {
                let expect = o.translated_vec(&x.multiplied_scalar(a * u.cos()))
                    .translated_vec(&y.multiplied_scalar(-b * u.sin()));
                if c.d0(u).distance(&expect) > 1e-6 * a.max(b).max(1.0) {
                    return None;
                }
            }
            Some(ConicGeom::Ellipse { o, a, b, x, y })
        }
        _ => None,
    }
}

/// Reconstruct a line as `d0(u) = loc + v·u` (v need not be unit: for a
/// `GeomLine` it is the unit direction; for a trimmed/other line-like curve it
/// is the chord direction scaled to the parameter range).
fn line_geometry(c: &dyn Curve) -> Option<(GpPnt, GpVec)> {
    let p0 = c.d0(0.0);
    let p1 = c.d0(1.0);
    let v = GpVec::from_pnts(&p0, &p1);
    if v.magnitude() < 1e-30 {
        return None;
    }
    Some((p0, v))
}

// ---------------------------------------------------------------------------
// Surface classification and reconstruction
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

/// Reconstructed quadric surface geometry. Each variant carries the data needed
/// for the implicit-equation solves and for exact `(u, v)` parameter recovery.
#[derive(Debug, Clone)]
enum SurfaceGeom {
    Plane { o: GpPnt, x: GpVec, y: GpVec, n: GpVec },
    Sphere { o: GpPnt, r: f64, x: GpVec, y: GpVec, z: GpVec },
    Cylinder { a: GpPnt, z: GpVec, r: f64 },
    Cone { apex: GpPnt, z: GpVec, cosa: f64 },
    Torus { o: GpPnt, major: f64, minor: f64, x: GpVec, z: GpVec },
}

/// Finite, sane sampling bounds (unbounded directions clamp to ±1).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// The natural `(0, 0)`-anchored sampling point for a quadric's frame axes
/// (all ported quadrics start their periodic U at 0).
fn quadric_origin_u(s: &dyn Surface) -> f64 {
    let (u0, _) = s.u_range();
    if u0.is_finite() { u0 } else { 0.0 }
}

fn classify_surface_kind(s: &dyn Surface) -> SurfaceKind {
    if is_planar(s, 8, 8, 1e-6) {
        return SurfaceKind::Plane;
    }
    if let Some(geom) = build_sphere_geom(s) {
        if sphere_matches(s, &geom) {
            return SurfaceKind::Sphere;
        }
    }
    if let Some(geom) = build_cylinder_geom(s) {
        if cylinder_matches(s, &geom) {
            return SurfaceKind::Cylinder;
        }
    }
    if let Some(geom) = build_cone_geom(s) {
        if cone_matches(s, &geom) {
            return SurfaceKind::Cone;
        }
    }
    if let Some(geom) = build_torus_geom(s) {
        if torus_matches(s, &geom) {
            return SurfaceKind::Torus;
        }
    }
    SurfaceKind::Other
}

fn build_surface_geom(s: &dyn Surface, kind: SurfaceKind) -> Option<SurfaceGeom> {
    match kind {
        SurfaceKind::Plane => {
            let u0 = match s.u_range().0 {
                a if a.is_finite() => a,
                _ => 0.0,
            };
            let v0 = match s.v_range().0 {
                a if a.is_finite() => a,
                _ => 0.0,
            };
            let o = s.d0(u0, v0);
            let x = GpVec::from_pnts(&o, &s.d0(u0 + 1.0, v0));
            let y = GpVec::from_pnts(&o, &s.d0(u0, v0 + 1.0));
            let n = x.crossed(&y);
            if n.magnitude() < 1e-12 {
                return None;
            }
            Some(SurfaceGeom::Plane { o, x, y, n })
        }
        SurfaceKind::Sphere => build_sphere_geom(s),
        SurfaceKind::Cylinder => build_cylinder_geom(s),
        SurfaceKind::Cone => build_cone_geom(s),
        SurfaceKind::Torus => build_torus_geom(s),
        SurfaceKind::Other => None,
    }
}

/// Sphere: center from `sphere_center`, radius from a sample, orthonormal
/// frame from the natural (0, 0), (π/2, 0), (0, π/2) samples.
fn build_sphere_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let o = sphere_center(s)?;
    let u0 = quadric_origin_u(s);
    let r = s.d0(u0, 0.0).distance(&o);
    if r < 1e-12 {
        return None;
    }
    let x = GpVec::from_pnts(&o, &s.d0(u0, 0.0)).divided(r);
    let y = GpVec::from_pnts(&o, &s.d0(u0 + std::f64::consts::FRAC_PI_2, 0.0)).divided(r);
    let z = GpVec::from_pnts(&o, &s.d0(u0, std::f64::consts::FRAC_PI_2)).divided(r);
    if !frame_ok(&x, &y, &z) {
        return None;
    }
    Some(SurfaceGeom::Sphere { o, r, x, y, z })
}

/// Verify a candidate sphere reconstruction against a sample grid.
fn sphere_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Sphere { o, r, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
    let tol = 1e-4 * r.abs().max(1.0);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            if (s.d0(u, v).distance(o) - r).abs() > tol {
                return false;
            }
        }
    }
    true
}

/// Cylinder: axis point = circle-center at the v start, axis direction from the
/// v advance, radius = half the opposite-sample distance.
fn build_cylinder_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let u0 = quadric_origin_u(s);
    let (_, _, v0, _) = sample_bounds(s);
    let c0 = midpoint(&s.d0(u0, v0), &s.d0(u0 + std::f64::consts::PI, v0));
    let c1 = midpoint(&s.d0(u0, v0 + 1.0), &s.d0(u0 + std::f64::consts::PI, v0 + 1.0));
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
    let x = GpVec::from_pnts(&c0, &s.d0(u0, v0)).divided(r);
    let y = GpVec::from_pnts(&c0, &s.d0(u0 + std::f64::consts::FRAC_PI_2, v0)).divided(r);
    // Reject a degenerate frame (x, y not perpendicular to z).
    if x.crossed(&y).dot(&z).abs() < 1e-6 {
        return None;
    }
    Some(SurfaceGeom::Cylinder { a: c0, z, r })
}

/// Verify a candidate cylinder reconstruction against a sample grid.
fn cylinder_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Cylinder { a, z, r, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
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

/// Cone: axis from the circle-center line, apex at the extrapolated zero
/// radius, semi-angle from the radius growth rate.
fn build_cone_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let u0 = quadric_origin_u(s);
    let (_, _, v0, _) = sample_bounds(s);
    let c0 = midpoint(&s.d0(u0, v0), &s.d0(u0 + std::f64::consts::PI, v0));
    let c1 = midpoint(&s.d0(u0, v0 + 1.0), &s.d0(u0 + std::f64::consts::PI, v0 + 1.0));
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
    let cosa = 1.0 / (1.0 + tan_alpha * tan_alpha).sqrt();
    // Distance from apex to the v0 circle plane along the axis: r0 / tanα.
    let s0 = r0 / tan_alpha;
    let apex = c0.translated_vec(&z.multiplied_scalar(-s0));
    Some(SurfaceGeom::Cone { apex, z, cosa })
}

/// Verify a candidate cone reconstruction against a sample grid.
fn cone_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Cone { apex, z, cosa, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
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

/// Torus: center from the v = 0 circle, major/minor radii from the v = 0 and
/// v = π circle radii, frame from the natural parameterization.
fn build_torus_geom(s: &dyn Surface) -> Option<SurfaceGeom> {
    let u0 = quadric_origin_u(s);
    let pi = std::f64::consts::PI;
    let o = midpoint(&s.d0(u0, 0.0), &s.d0(u0 + pi, 0.0));
    let r_outer = s.d0(u0, 0.0).distance(&o);
    let o2 = midpoint(&s.d0(u0, pi), &s.d0(u0 + pi, pi));
    let r_inner = s.d0(u0, pi).distance(&o2);
    if o.distance(&o2) > 1e-6 * r_outer.max(1.0) {
        return None;
    }
    let major = 0.5 * (r_outer + r_inner);
    let minor = 0.5 * (r_outer - r_inner);
    if major <= 1e-12 || minor <= 1e-12 {
        return None;
    }
    let x = GpVec::from_pnts(&o, &s.d0(u0, 0.0)).divided(r_outer);
    let y = GpVec::from_pnts(&o, &s.d0(u0 + pi / 2.0, 0.0)).divided(r_outer);
    let z = GpVec::from_pnts(
        &o,
        &s.d0(u0, pi / 2.0),
    )
    .subtracted(&x.multiplied_scalar(major))
    .divided(minor);
    if !frame_ok(&x, &y, &z) {
        return None;
    }
    Some(SurfaceGeom::Torus { o, major, minor, x, z })
}

/// Verify a candidate torus reconstruction against a sample grid.
fn torus_matches(s: &dyn Surface, geom: &SurfaceGeom) -> bool {
    let SurfaceGeom::Torus { o, major, minor, z, .. } = geom else { return false; };
    let (u0, u1, v0, v1) = sample_bounds(s);
    let tol = 1e-4 * (*major).max(*minor).max(1.0);
    for i in 0..8 {
        for j in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 7.0;
            let v = v0 + (v1 - v0) * j as f64 / 7.0;
            let d = GpVec::from_pnts(o, &s.d0(u, v));
            let dz = d.dot(z);
            let rho = d.coord.subtracted(&z.xyz().multiplied(dz)).modulus();
            let err = ((rho - major) * (rho - major) + dz * dz - minor * minor).abs();
            if err > tol * (rho + major + minor).max(1.0) {
                return false;
            }
        }
    }
    true
}

/// Whether x, y, z form an orthonormal frame (up to scale/reflection noise).
fn frame_ok(x: &GpVec, y: &GpVec, z: &GpVec) -> bool {
    (x.magnitude() - 1.0).abs() < 1e-6
        && (y.magnitude() - 1.0).abs() < 1e-6
        && (z.magnitude() - 1.0).abs() < 1e-6
        && x.dot(y).abs() < 1e-6
        && x.dot(z).abs() < 1e-6
        && y.dot(z).abs() < 1e-6
}

/// Wrap `x` into `[start, start + period)` (positive modulo).
fn wrap_periodic(x: f64, start: f64, period: f64) -> f64 {
    let p = period.abs();
    if p <= 1e-30 {
        return x;
    }
    start + (x - start).rem_euclid(p)
}

/// Snap a cosine/sine value to ±1 when within `eps` of a pole, so `asin` does
/// not amplify a ~1e-16 frame error into a ~1e-8 angular error.
fn snap_pole(s: f64) -> f64 {
    if s > 1.0 - 1e-9 {
        1.0
    } else if s < -1.0 + 1e-9 {
        -1.0
    } else {
        s
    }
}

/// Surface `(u, v)` parameters of a point, exact for reconstructed quadrics
/// (so the reported parameters are precise and the UV-bounds validation is
/// meaningful), falling back to the sampling projector for other surfaces.
/// Periodic U/V directions are wrapped into the surface's natural range.
fn surface_params(surface: &dyn Surface, geom: Option<&SurfaceGeom>, p: &GpPnt) -> (f64, f64) {
    let (u0, _, _, _) = sample_bounds(surface);
    if let Some(g) = geom {
        match g {
            // Plane and sphere frames are reconstructed from the surface's own
            // natural parameterization, so the recovered (u, v) are exact.
            SurfaceGeom::Plane { o, x, y, .. } => {
                let d = GpVec::from_pnts(o, p);
                return (d.dot(x), d.dot(y));
            }
            SurfaceGeom::Sphere { o, r, x, y, z } => {
                let d = GpVec::from_pnts(o, p);
                let u = wrap_periodic(d.dot(y).atan2(d.dot(x)), u0, 2.0 * std::f64::consts::PI);
                let s = snap_pole(d.dot(z) / r);
                return (u, s.asin());
            }
            // Cylinder/cone/torus frames use a reconstructed axis/reference that
            // can differ by a constant offset from the surface's natural
            // parameterization; the sampling projector keeps (u, v) consistent.
            _ => return surface_closest_params(surface, p, 24, 24),
        }
    }
    surface_closest_params(surface, p, 24, 24)
}

// ---------------------------------------------------------------------------
// Analytic solves
// ---------------------------------------------------------------------------

/// Result of an analytic curve-vs-quadric solve in the curve parameter.
enum SolveResult {
    /// Discrete roots in the curve parameter.
    Roots(Vec<f64>),
    /// The curve lies on the surface (coincidence) — an `On` segment.
    OnSurface,
    /// The curve is disjoint (including parallel-and-offset).
    None,
}

/// Real roots of a·x² + b·x + c = 0 (stable), with the linear/parallel
/// degeneracy folded in.
fn solve_quadratic(a: f64, b: f64, c: f64, surface_eps: f64) -> SolveResult {
    if a.abs() < 1e-300 {
        if b.abs() < 1e-300 {
            if c.abs() <= surface_eps {
                return SolveResult::OnSurface;
            }
            return SolveResult::None;
        }
        return SolveResult::Roots(vec![-c / b]);
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return SolveResult::None;
    }
    if disc.abs() <= 1e-12 * a.abs().max(1.0) {
        return SolveResult::Roots(vec![-b / (2.0 * a)]);
    }
    let sq = disc.sqrt();
    let q = -0.5 * (b + b.signum() * sq);
    let r1 = q / a;
    let r2 = if q.abs() > 1e-300 { c / q } else { (-b - sq) / (2.0 * a) };
    SolveResult::Roots(vec![r1, r2])
}

fn line_plane(loc: &GpPnt, v: &GpVec, o: &GpPnt, n: &GpVec) -> SolveResult {
    let denom = n.dot(v);
    if denom.abs() < ANG_TOL {
        // Parallel to the plane: intersects only if it lies in it.
        let d = GpVec::from_pnts(o, loc);
        if n.dot(&d).abs() <= ANALYTIC_TOL {
            SolveResult::OnSurface
        } else {
            SolveResult::None
        }
    } else {
        let t = n.dot(&GpVec::from_pnts(loc, o)) / denom;
        SolveResult::Roots(vec![t])
    }
}

fn line_sphere(loc: &GpPnt, v: &GpVec, c: &GpPnt, r: f64) -> SolveResult {
    let d = GpVec::from_pnts(c, loc); // loc − c
    let a = v.square_magnitude();
    if a < 1e-30 {
        return SolveResult::None;
    }
    let b = 2.0 * d.dot(v);
    let cc = d.square_magnitude() - r * r;
    solve_quadratic(a, b, cc, 1e-9)
}

fn line_cylinder(loc: &GpPnt, v: &GpVec, a: &GpPnt, z: &GpVec, r: f64) -> SolveResult {
    let d0 = GpVec::from_pnts(a, loc); // loc − a
    let d0z = d0.dot(z);
    let dz = v.dot(z);
    let a2 = v.square_magnitude() - dz * dz;
    let b = 2.0 * (d0.dot(v) - d0z * dz);
    let c = d0.square_magnitude() - d0z * d0z - r * r;
    solve_quadratic(a2, b, c, 1e-9)
}

fn line_cone(loc: &GpPnt, v: &GpVec, apex: &GpPnt, z: &GpVec, cosa: f64) -> SolveResult {
    let d0 = GpVec::from_pnts(apex, loc); // loc − apex
    let d0z = d0.dot(z);
    let dz = v.dot(z);
    let v2 = v.square_magnitude();
    let a = dz * dz - v2 * cosa * cosa;
    let b = 2.0 * (d0z * dz - d0.dot(v) * cosa * cosa);
    let c = d0z * d0z - d0.square_magnitude() * cosa * cosa;
    solve_quadratic(a, b, c, 1e-9)
}

/// `a·cos(u) + b·sin(u) = c` roots in `[inf, sup]` via the trigonometric
/// reduction (`math_TrigonometricFunctionRoots`).
fn trig_solve(a: f64, b: f64, c: f64, inf: f64, sup: f64, surface_eps: f64) -> SolveResult {
    if a.abs() < 1e-12 && b.abs() < 1e-12 {
        if c.abs() <= surface_eps {
            return SolveResult::OnSurface;
        }
        return SolveResult::None;
    }
    SolveResult::Roots(occt_math::trig_roots(a, b, -c, inf, sup))
}

// ---------------------------------------------------------------------------
// Main entry
// ---------------------------------------------------------------------------

/// Compute the intersection of `curve` with `surface`.
///
/// * `cu_range` — the curve parameter window to search (`(first, last)`); for
///   unbounded curves (lines) this must be finite.
/// * `uv_range` — the surface `(u0, v0, u1, v1)` window to search; unbounded
///   surface directions are clamped to it.
///
/// Returns discrete intersection points (sorted by curve parameter) and, for
/// the coincidence degeneracy (e.g. a line lying in a plane), `On` segments.
pub fn perform_curve_surface(
    curve: &dyn Curve,
    surface: &dyn Surface,
    cu_range: (f64, f64),
    uv_range: (f64, f64, f64, f64),
) -> Result<HInterResult, String> {
    let (t0, t1) = effective_curve_range(curve, cu_range)?;
    let uv = effective_uv_range(surface, uv_range)?;

    let curve_kind = classify_curve(curve);
    let surface_kind = classify_surface_kind(surface);
    let geom = build_surface_geom(surface, surface_kind);

    let mut points: Vec<IntersectionPoint> = Vec::new();
    let mut segments: Vec<IntersectionSegment> = Vec::new();

    // Analytic path for conic curves against supported quadrics.
    if curve_kind != CurveKind::Other {
        if let Some((mut pts, mut segs)) = geom
            .as_ref()
            .and_then(|g| perform_conic_surf(curve, curve_kind, surface, g, t0, t1, uv))
        {
            points.append(&mut pts);
            segments.append(&mut segs);
        } else {
            let mut pts = general_curve_surface(curve, surface, geom.as_ref(), t0, t1, uv, GENERAL_TOL);
            points.append(&mut pts);
        }
    } else {
        let mut pts = general_curve_surface(curve, surface, geom.as_ref(), t0, t1, uv, GENERAL_TOL);
        points.append(&mut pts);
    }

    // Sort by curve parameter and merge near-duplicate (tangent) points.
    points.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(Ordering::Equal));
    let mut deduped: Vec<IntersectionPoint> = Vec::with_capacity(points.len());
    for p in points {
        if let Some(last) = deduped.last_mut() {
            if (p.param - last.param).abs() < MERGE_TOL && last.pnt.distance(&p.pnt) <= 1e-6 {
                last.state = State::On;
                continue;
            }
        }
        deduped.push(p);
    }

    Ok(HInterResult { points: deduped, segments })
}

/// Clamp the curve parameter window to the curve's own range. Unbounded curves
/// (lines) rely on the caller-provided finite window.
fn effective_curve_range(curve: &dyn Curve, cu_range: (f64, f64)) -> Result<(f64, f64), String> {
    let (ca, cb) = cu_range;
    let (lo, hi) = if ca <= cb { (ca, cb) } else { (cb, ca) };
    let fa = curve.first_parameter();
    let fb = curve.last_parameter();
    let a = if fa.is_finite() { fa.max(lo) } else { lo };
    let b = if fb.is_finite() { fb.min(hi) } else { hi };
    if !a.is_finite() || !b.is_finite() {
        return Err("perform_curve_surface: unbounded curve parameter range".into());
    }
    if b - a <= 1e-15 {
        return Err("perform_curve_surface: empty curve range".into());
    }
    Ok((a, b))
}

/// Clamp the surface `(u0, v0, u1, v1)` window to the surface's own ranges.
/// Unbounded directions take the caller's bounds directly.
fn effective_uv_range(surface: &dyn Surface, uv: (f64, f64, f64, f64)) -> Result<(f64, f64, f64, f64), String> {
    let (u0, u1, v0, v1) = uv;
    let (su0, su1) = surface.u_range();
    let (sv0, sv1) = surface.v_range();
    let clamp_pair = |a: f64, b: f64, sa: f64, sb: f64| -> (f64, f64) {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let a = if sa.is_finite() { a.max(sa) } else { a };
        let b = if sb.is_finite() { b.min(sb) } else { b };
        (a, b)
    };
    let (u0, u1) = clamp_pair(u0, u1, su0, su1);
    let (v0, v1) = clamp_pair(v0, v1, sv0, sv1);
    if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
        return Err("perform_curve_surface: unbounded surface UV range".into());
    }
    if u1 - u0 <= 1e-15 || v1 - v0 <= 1e-15 {
        return Err("perform_curve_surface: empty surface UV range".into());
    }
    Ok((u0, u1, v0, v1))
}

/// Analytic dispatch for conic curves against quadric surfaces. Returns
/// `None` when the combination has no analytic solve (caller falls back to the
/// general path).
fn perform_conic_surf(
    curve: &dyn Curve,
    kind: CurveKind,
    surface: &dyn Surface,
    geom: &SurfaceGeom,
    t0: f64,
    t1: f64,
    uv: (f64, f64, f64, f64),
) -> Option<(Vec<IntersectionPoint>, Vec<IntersectionSegment>)> {
    let empty = (Vec::new(), Vec::new());
    let result = match (kind, geom) {
        (CurveKind::Line, SurfaceGeom::Plane { o, n, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_plane(&loc, &v, o, n)
        }
        (CurveKind::Line, SurfaceGeom::Sphere { o, r, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_sphere(&loc, &v, o, *r)
        }
        (CurveKind::Line, SurfaceGeom::Cylinder { a, z, r, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_cylinder(&loc, &v, a, z, *r)
        }
        (CurveKind::Line, SurfaceGeom::Cone { apex, z, cosa, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            line_cone(&loc, &v, apex, z, *cosa)
        }
        (CurveKind::Line, SurfaceGeom::Torus { o, major, minor, x, z, .. }) => {
            let (loc, v) = match line_geometry(curve) {
                Some(x) => x,
                None => return Some(empty),
            };
            let dir = match GpDir::from_vec(&v) {
                Ok(d) => d,
                Err(_) => return Some(empty),
            };
            let lin = GpLin::from_pnt_dir(loc, dir);
            let ax3 = match (
                GpDir::from_vec(z),
                GpDir::from_vec(x),
            ) {
                (Ok(zd), Ok(xd)) => match GpAx3::new(*o, zd, &xd) {
                    Ok(a) => a,
                    Err(_) => return Some(empty),
                },
                _ => return Some(empty),
            };
            let torus = match GpTorus::new(ax3, *major, *minor) {
                Ok(t) => t,
                Err(_) => return Some(empty),
            };
            let pts = line_torus_intersect(&lin, &torus);
            let v2 = v.square_magnitude();
            SolveResult::Roots(pts.into_iter().map(|p| GpVec::from_pnts(&loc, &p).dot(&v) / v2).collect())
        }
        (CurveKind::Circle, SurfaceGeom::Plane { o, n, .. }) => {
            let Some(ConicGeom::Circle { o: co, r, x, y }) = conic_geometry(curve, kind) else {
                return None;
            };
            let dq = GpVec::from_pnts(&co, o);
            trig_solve(r * n.dot(&x), r * n.dot(&y), n.dot(&dq), t0, t1, ANALYTIC_TOL)
        }
        (CurveKind::Circle, SurfaceGeom::Sphere { o, r: sph_r, .. }) => {
            let Some(ConicGeom::Circle { o: co, r, x, y }) = conic_geometry(curve, kind) else {
                return None;
            };
            let d = GpVec::from_pnts(o, &co); // circle center − sphere center
            let a = 2.0 * r * d.dot(&x);
            let b = 2.0 * r * d.dot(&y);
            let cc = sph_r * sph_r - r * r - d.square_magnitude();
            trig_solve(a, b, cc, t0, t1, ANALYTIC_TOL)
        }
        (CurveKind::Ellipse, SurfaceGeom::Plane { o, n, .. }) => {
            let Some(ConicGeom::Ellipse { o: eo, a, b, x, y }) = conic_geometry(curve, kind) else {
                return None;
            };
            let dq = GpVec::from_pnts(&eo, o);
            trig_solve(a * n.dot(&x), -b * n.dot(&y), n.dot(&dq), t0, t1, ANALYTIC_TOL)
        }
        _ => return None,
    };

    Some(match result {
        SolveResult::Roots(roots) => {
            let pts = roots
                .into_iter()
                .filter_map(|w| {
                    if w < t0 - ANALYTIC_TOL || w > t1 + ANALYTIC_TOL {
                        return None;
                    }
                    compute_append_point(curve, surface, Some(geom), w, uv, ANALYTIC_TOL)
                })
                .collect();
            (pts, Vec::new())
        }
        SolveResult::OnSurface => (Vec::new(), on_surface_segment(curve, surface, geom, t0, t1)),
        SolveResult::None => (Vec::new(), Vec::new()),
    })
}

/// Build an `On` coincidence segment spanning `[t0, t1]`.
fn on_surface_segment(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: &SurfaceGeom,
    t0: f64,
    t1: f64,
) -> Vec<IntersectionSegment> {
    let p0 = curve.d0(t0);
    let p1 = curve.d0(t1);
    let (u0, v0) = surface_params(surface, Some(geom), &p0);
    let (u1, v1) = surface_params(surface, Some(geom), &p1);
    let a = IntersectionPoint::new(t0, u0, v0, p0, State::On);
    let b = IntersectionPoint::new(t1, u1, v1, p1, State::On);
    vec![IntersectionSegment::new(a, b)]
}

/// Validate a candidate parameter and build an [`IntersectionPoint`]: the
/// curve point must lie on the surface (within `tol`) and its surface
/// parameters must fall inside `uv`.
fn compute_append_point(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    w: f64,
    uv: (f64, f64, f64, f64),
    tol: f64,
) -> Option<IntersectionPoint> {
    let pnt = curve.d0(w);
    let (su, sv) = surface_params(surface, geom, &pnt);
    let (u0, u1, v0, v1) = uv;
    if su < u0 - tol || su > u1 + tol || sv < v0 - tol || sv > v1 + tol {
        return None;
    }
    // Validate the point lies on the surface. Analytic points are checked
    // against the quadric's implicit equation (exact, and independent of the
    // reported (u, v)); the general path uses the normal (signed) distance.
    let on_surface = match geom {
        Some(g) => quadric_on_surface(g, &pnt, tol.max(1e-6)),
        None => signed_dist(surface, None, &pnt).abs() <= tol.max(1e-6),
    };
    if !on_surface {
        return None;
    }
    let state = compute_state(curve, surface, geom, w, su, sv);
    Some(IntersectionPoint::new(w, su, sv, pnt, state))
}

/// Exact outward normal of a reconstructed quadric at a point on it. Used for
/// the transition classification so tangencies are not corrupted by the
/// first-order error of the finite-difference normal (which tilts the normal by
/// a curvature term).
fn quadric_normal(geom: &SurfaceGeom, p: &GpPnt) -> Option<GpVec> {
    match geom {
        SurfaceGeom::Plane { n, .. } => {
            let m = n.magnitude();
            if m < 1e-12 {
                None
            } else {
                Some(n.divided(m))
            }
        }
        SurfaceGeom::Sphere { o, .. } => {
            let d = GpVec::from_pnts(o, p);
            let m = d.magnitude();
            if m < 1e-12 {
                None
            } else {
                Some(d.divided(m))
            }
        }
        SurfaceGeom::Cylinder { a, z, .. } => {
            let d = GpVec::from_pnts(a, p);
            let radial = d.coord.subtracted(&z.xyz().multiplied(d.dot(z)));
            let m = radial.modulus();
            if m < 1e-12 {
                None
            } else {
                Some(GpVec::from_xyz(&radial.divided(m)))
            }
        }
        SurfaceGeom::Cone { apex, z, cosa, .. } => {
            let d = GpVec::from_pnts(apex, p);
            // f = |d|²·cos²α − (d·Z)²; ∇f ∝ cos²α·d − (d·Z)·Z.
            let grad = d.coord.multiplied(cosa * cosa).subtracted(&z.xyz().multiplied(d.dot(z)));
            let m = grad.modulus();
            if m < 1e-12 {
                None
            } else {
                Some(GpVec::from_xyz(&grad.divided(m)))
            }
        }
        SurfaceGeom::Torus { o, major, z, .. } => {
            let d = GpVec::from_pnts(o, p);
            let dz = d.dot(z);
            let radial = d.coord.subtracted(&z.xyz().multiplied(dz));
            let rho = radial.modulus();
            if rho < 1e-12 {
                return None;
            }
            // f = (ρ − R)² + z² − r²; ∇f ∝ (ρ − R)/ρ·d_perp + z·Z.
            let g = radial.multiplied((rho - major) / rho).added(&z.xyz().multiplied(dz));
            let m = g.modulus();
            if m < 1e-12 {
                None
            } else {
                Some(GpVec::from_xyz(&g.divided(m)))
            }
        }
    }
}

/// Whether a point lies on a reconstructed quadric, via the quadric's implicit
/// equation (relative tolerance). Independent of the surface parameterization,
/// so it stays exact for quadrics whose reconstructed reference differs from
/// the natural `(u, v)` origin.
fn quadric_on_surface(geom: &SurfaceGeom, p: &GpPnt, tol: f64) -> bool {
    match geom {
        SurfaceGeom::Plane { o, n, .. } => {
            let d = GpVec::from_pnts(o, p);
            d.dot(n).abs() / n.magnitude() <= tol
        }
        SurfaceGeom::Sphere { o, r, .. } => (p.distance(o) - r).abs() <= tol,
        SurfaceGeom::Cylinder { a, z, r, .. } => {
            let d = GpVec::from_pnts(a, p);
            let rho = d.coord.subtracted(&z.xyz().multiplied(d.dot(z))).modulus();
            (rho - r).abs() <= tol
        }
        SurfaceGeom::Cone { apex, z, cosa, .. } => {
            let d = GpVec::from_pnts(apex, p);
            let lhs = d.dot(z) * d.dot(z);
            let rhs = d.square_magnitude() * cosa * cosa;
            (lhs - rhs).abs() <= tol * d.square_magnitude().max(1.0)
        }
        SurfaceGeom::Torus { o, major, minor, z, .. } => {
            let d = GpVec::from_pnts(o, p);
            let dz = d.dot(z);
            let rho = d.coord.subtracted(&z.xyz().multiplied(dz)).modulus();
            let lhs = (rho - major) * (rho - major) + dz * dz;
            (lhs - minor * minor).abs() <= tol * (rho + major + minor).max(1.0)
        }
    }
}

/// Transition of the curve at the crossing: aligned with the surface normal →
/// `Out`, opposed → `In`, perpendicular → `On`.
fn compute_state(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    w: f64,
    su: f64,
    sv: f64,
) -> State {
    let (_, d1) = curve.d1(w);
    let dm = d1.magnitude();
    if dm < 1e-12 {
        return State::Unknown;
    }
    let n = match geom.and_then(|g| quadric_normal(g, &curve.d0(w))) {
        Some(n) => n,
        None => surface_normal(surface, su, sv),
    };
    let nm = n.magnitude();
    if nm < 1e-12 {
        return State::Unknown;
    }
    let cos = d1.dot(&n) / (dm * nm);
    if -cos > ANG_TOL {
        State::In
    } else if cos > ANG_TOL {
        State::Out
    } else {
        State::On
    }
}

// ---------------------------------------------------------------------------
// General path (sampling + bisection)
// ---------------------------------------------------------------------------

/// Signed distance from a point to `surface`: positive when the point is on the
/// side of the surface normal, negative otherwise. When the reconstructed
/// quadric geometry is available its exact normal is used (the trait-level
/// `surface_normal` finite-difference fallback is ill-defined on unbounded
/// parameter directions); otherwise the sampling normal is used.
fn signed_dist(surface: &dyn Surface, geom: Option<&SurfaceGeom>, p: &GpPnt) -> f64 {
    let (su, sv) = surface_closest_params(surface, p, 16, 16);
    let q = surface.d0(su, sv);
    let n = match geom.and_then(|g| quadric_normal(g, p)) {
        Some(n) => n,
        None => surface_normal(surface, su, sv),
    };
    let d = GpVec::from_pnts(&q, p);
    let nm = n.magnitude();
    if nm < 1e-12 {
        d.magnitude()
    } else {
        d.dot(&n) / nm
    }
}

/// General sampling path: uniform curve sampling, sign-change bisection for
/// transverse crossings, golden-section refinement of near-zero distance
/// minima for tangencies.
fn general_curve_surface(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    t0: f64,
    t1: f64,
    uv: (f64, f64, f64, f64),
    tol: f64,
) -> Vec<IntersectionPoint> {
    let n = 512usize;
    let mut points: Vec<IntersectionPoint> = Vec::new();
    let mut prev_u = t0;
    let mut prev_s = signed_dist(surface, geom, &curve.d0(t0));
    for i in 1..=n {
        let u = t0 + (t1 - t0) * i as f64 / n as f64;
        let s = signed_dist(surface, geom, &curve.d0(u));
        if prev_s * s < 0.0 {
            // Transverse crossing.
            if let Some(w) = bisect_root(curve, surface, geom, prev_u, u, prev_s, tol) {
                if let Some(pt) = compute_append_point(curve, surface, geom, w, uv, tol) {
                    points.push(pt);
                }
            }
        } else if s.abs() <= tol.max(1e-6) && prev_s.abs() > s.abs() {
            // Near-zero minimum → likely tangency.
            if let Some(w) = refine_min(curve, surface, geom, prev_u, u, tol) {
                if let Some(mut pt) = compute_append_point(curve, surface, geom, w, uv, tol) {
                    pt.state = State::On;
                    points.push(pt);
                }
            }
        }
        prev_u = u;
        prev_s = s;
    }
    points
}

/// Bisect a sign change of the signed distance over `[ua, ub]`.
fn bisect_root(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    ua: f64,
    ub: f64,
    sa: f64,
    tol: f64,
) -> Option<f64> {
    let mut lo = ua;
    let mut hi = ub;
    let mut flo = sa;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let sm = signed_dist(surface, geom, &curve.d0(mid));
        if flo * sm <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo = sm;
        }
        if (hi - lo).abs() < 1e-10 {
            break;
        }
        if sm.abs() < 1e-10 {
            break;
        }
    }
    let w = 0.5 * (lo + hi);
    let p = curve.d0(w);
    if signed_dist(surface, geom, &p).abs() <= tol.max(1e-6) {
        Some(w)
    } else {
        None
    }
}

/// Golden-section refinement of the unsigned distance minimum over `[lo, hi]`.
fn refine_min(
    curve: &dyn Curve,
    surface: &dyn Surface,
    geom: Option<&SurfaceGeom>,
    lo: f64,
    hi: f64,
    tol: f64,
) -> Option<f64> {
    let f = |u: f64| signed_dist(surface, geom, &curve.d0(u)).abs();
    let (u, d) = golden_1d(&f, lo, hi, 1e-10);
    if d <= tol.max(1e-6) {
        Some(u)
    } else {
        None
    }
}

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
fn golden_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64, eps: f64) -> (f64, f64) {
    const GOLD: f64 = 0.618_033_988_749_894_9;
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
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpPln, GpSphere, GpTorus};
    use occt_geom::{GeomBSplineCurve, GeomCircle, GeomCone, GeomCylinder, GeomLine, GeomPlane, GeomSphere, GeomTorus};

    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_tool::BRepTool;

    const PI: f64 = std::f64::consts::PI;

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    fn sphere_at(origin: GpPnt, r: f64) -> GeomSphere {
        let ax3 = GpAx3::new(origin, dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        GeomSphere::new(GpSphere::new(ax3, r).unwrap())
    }

    #[test]
    fn line_through_unit_sphere_two_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, -2.0), dir(0.0, 0.0, 1.0));
        // z = −1 at t = 1, z = +1 at t = 3.
        let res = perform_curve_surface(&line, &sphere, (-2.0, 4.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        let p0 = res.point(0).pnt();
        let p1 = res.point(1).pnt();
        assert!((p0.coord.modulus() - 1.0).abs() < 1e-9, "p0 {:?}", p0);
        assert!((p1.coord.modulus() - 1.0).abs() < 1e-9, "p1 {:?}", p1);
        assert!((p0.z() + 1.0).abs() < 1e-9 || (p1.z() + 1.0).abs() < 1e-9);
        assert!(res.point(0).param() < res.point(1).param());
    }

    #[test]
    fn line_outside_sphere_no_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let line = GeomLine::from_pnt_dir(GpPnt::new(3.0, 0.0, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, &sphere, (-1.0, 1.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
    }

    #[test]
    fn line_through_plane_single_point() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, &plane, (-1.0, 1.0), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!(p.distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-9, "p {:?}", p);
    }

    #[test]
    fn line_in_plane_is_on_segment() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        // Line lying in z = 0.
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, 0.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &plane, (-2.0, 2.0), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 1, "segments: {:?}", res.segments());
        assert_eq!(res.segment(0).first_point().state(), State::On);
    }

    #[test]
    fn circle_in_plane_is_on_segment() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let res = perform_curve_surface(&circle, &plane, (0.0, 2.0 * PI), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 1, "segments: {:?}", res.segments());
        assert_eq!(res.segment(0).first_point().state(), State::On);
    }

    #[test]
    fn circle_parallel_to_plane_no_intersection() {
        // Circle in z = 1, plane z = 0.
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 1.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &plane, (0.0, 2.0 * PI), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 0);
    }

    #[test]
    fn circle_crosses_sphere_two_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Circle radius 1 centered (1.5, 0, 0) in the XY plane.
        let ax2 = GpAx2::new(GpPnt::new(1.5, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &sphere, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            assert!((p.coord.modulus() - 1.0).abs() < 1e-9, "p {:?}", p);
        }
    }

    #[test]
    fn circle_inside_sphere_no_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 0.5));
        let res = perform_curve_surface(&circle, &sphere, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 0);
    }

    #[test]
    fn bspline_through_sphere_general_path() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // A B-spline through (−2,0,0) and (2,0,0) with a bow toward z = 0.5;
        // it crosses the unit sphere near x = ±1.
        let poles = vec![
            GpPnt::new(-2.0, 0.0, 0.0),
            GpPnt::new(-0.5, 0.0, 0.6),
            GpPnt::new(0.5, 0.0, 0.6),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        // Degree-2 clamped B-spline: 4 poles → 4 + 2 + 1 = 7 knots.
        let curve = Arc::new(GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap());
        let res = perform_curve_surface(curve.as_ref(), &sphere, (0.0, 1.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            assert!((p.coord.modulus() - 1.0).abs() < 1e-5, "p {:?}", p);
        }
    }

    #[test]
    fn box_face_surface_crossed_by_line() {
        let b = unit_box();
        // Bottom face (z = 0).
        let face = &b.faces[0];
        let surf = BRepTool::face_surface(face).expect("face surface");
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.5, 0.5, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, surf.as_ref(), (-1.0, 1.0), (0.0, 1.0, 0.0, 1.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!(p.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6, "p {:?}", p);
    }

    #[test]
    fn box_face_surface_missed_by_line() {
        let b = unit_box();
        let face = &b.faces[0]; // bottom, z = 0
        let surf = BRepTool::face_surface(face).expect("face surface");
        // Line outside the face's UV domain (x beyond the box).
        let line = GeomLine::from_pnt_dir(GpPnt::new(5.0, 5.0, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, surf.as_ref(), (-1.0, 1.0), (0.0, 1.0, 0.0, 1.0)).unwrap();
        assert_eq!(res.nb_points(), 0, "points: {:?}", res.points());
    }

    #[test]
    fn tangent_line_touches_sphere_single_on_point() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Line at y = 1, tangent to the unit sphere at (0, 1, 0).
        let line = GeomLine::from_pnt_dir(GpPnt::new(-2.0, 1.0, 0.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &sphere, (-2.0, 2.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!(p.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-6, "p {:?}", p);
        assert_eq!(res.point(0).state(), State::On);
    }

    #[test]
    fn ellipse_crosses_plane_two_points() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        // Ellipse centered (0,0,0.5), semi-axes 2 (X) and 1 (Z), lying in the
        // XZ plane: d0(u) = (2·cos u, 0, 0.5 + sin u). It crosses the plane
        // z = 0 where sin u = −0.5 → u = 7π/6, 11π/6, giving two points.
        let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 0.5), dir(0.0, 1.0, 0.0), dir(1.0, 0.0, 0.0)).unwrap();
        let elips = occt_core::gp::GpElips::new(ax2, 2.0, 1.0);
        let ellipse = occt_geom::GeomEllipse::new(elips);
        let res = perform_curve_surface(&ellipse, &plane, (0.0, 2.0 * PI), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            assert!(p.z().abs() < 1e-9, "p {:?}", p);
        }
    }

    #[test]
    fn line_through_cylinder_two_points() {
        // Cylinder radius 1 along Z; line at y = 0.5, z = 0 crosses at
        // x = ±√(1 − 0.25) = ±0.866.
        let cyl = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        let surface = GeomCylinder::new(cyl);
        let line = GeomLine::from_pnt_dir(GpPnt::new(-2.0, 0.5, 0.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &surface, (-2.0, 4.0), (0.0, 2.0 * PI, -10.0, 10.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            let rho = (p.x() * p.x() + p.y() * p.y()).sqrt();
            assert!((rho - 1.0).abs() < 1e-6, "p {:?}", p);
        }
    }

    #[test]
    fn line_through_cone_two_points() {
        // Cone (geometric tip at (0,0,−1), radius 1 at z = 0, half-angle 45°):
        // radius at z = 1 is 2, so the line along X at z = 1 crosses at x = ±2.
        // Placement at the origin (radius 1 at z = 0); the tip is then at
        // z = −RefRadius/tan(45°) = −1.
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        let cone = GpCone::new(ax3, 1.0, PI / 4.0).unwrap();
        let surface = GeomCone::new(cone);
        let line = GeomLine::from_pnt_dir(GpPnt::new(-3.0, 0.0, 1.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &surface, (-2.0, 6.0), (0.0, 2.0 * PI, -10.0, 10.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            // On the cone: radius = z + 1 (tip at z = −1, slope 1).
            let rho = (p.x() * p.x() + p.y() * p.y()).sqrt();
            assert!((rho - (p.z() + 1.0)).abs() < 1e-6, "p {:?}", p);
        }
    }

    #[test]
    fn line_through_torus_four_points() {
        // Torus at origin R = 3, r = 1; line along X through the centre cuts at
        // x = ±2, ±4 (port of IntAna_IntLinTorus).
        let torus = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let surface = GeomTorus::new(torus);
        let line = GeomLine::from_pnt_dir(GpPnt::new(-6.0, 0.0, 0.0), dir(1.0, 0.0, 0.0));
        // x = −6 + t cuts the torus at t = 2, 4, 8, 10 (x = −4, −2, 2, 4).
        let res = perform_curve_surface(&line, &surface, (0.0, 12.0), (0.0, 2.0 * PI, 0.0, 2.0 * PI)).unwrap();
        assert_eq!(res.nb_points(), 4, "points: {:?}", res.points());
        let mut xs: Vec<f64> = res.points().iter().map(|p| p.pnt().x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (x, expected) in xs.iter().zip([-4.0, -2.0, 2.0, 4.0].iter()) {
            assert!((x - expected).abs() < 1e-5, "x {x} expected {expected}");
        }
    }

    #[test]
    fn line_entering_sphere_state_in_out() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, -2.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, &sphere, (-2.0, 4.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2);
        assert_eq!(res.point(0).state(), State::In, "first crossing enters");
        assert_eq!(res.point(1).state(), State::Out, "second crossing exits");
    }

    #[test]
    fn circle_on_sphere_is_on_segment() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Unit circle in the XY plane = the unit sphere's equator.
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let res = perform_curve_surface(&circle, &sphere, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 1, "segments: {:?}", res.segments());
        assert_eq!(res.segment(0).first_point().state(), State::On);
    }

    #[test]
    fn circle_crosses_cylinder_general_path() {
        let cyl = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        let surface = GeomCylinder::new(cyl);
        // Circle radius 1 centred (1.5, 0, 0) in the XY plane: points satisfy
        // x = 1.5 + cos u, y = sin u; distance to the Z axis is √(3.25 + 3 cos u),
        // equal to 1 when cos u = −0.75 → two points.
        let ax2 = GpAx2::new(GpPnt::new(1.5, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &surface, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -10.0, 10.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            let rho = (p.x() * p.x() + p.y() * p.y()).sqrt();
            assert!((rho - 1.0).abs() < 1e-5, "p {:?}", p);
        }
    }

    #[test]
    fn trimmed_circle_range_filters_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Circle radius 1 centred (1.5, 0, 0) in the XY plane crosses the unit
        // sphere at u = acos(−0.75) ≈ 2.42 and u = 2π − 2.42 ≈ 3.86. Restricting
        // to [2.0, 3.0] keeps only the first.
        let ax2 = GpAx2::new(GpPnt::new(1.5, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &sphere, (2.0, 3.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!((p.coord.modulus() - 1.0).abs() < 1e-9, "p {:?}", p);
    }
}
