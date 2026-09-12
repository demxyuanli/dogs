//! Curve/surface sampling classifiers used by EdgeFace.
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::PCONFUSION;
use occt_geom::{Curve, Surface};

use crate::brep_surface::{is_planar, sphere_center, surface_closest_params, SurfaceKind};

const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Curve classification helpers
// ---------------------------------------------------------------------------

/// Coarse analytic kind of a 3D curve, mirroring the `GeomAbs_CurveType` subset
/// the algorithm branches on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CurveKind {
    Line,
    Circle,
    BSpline,
    Other,
}

/// Whether the curve is geometrically a straight line (unbounded, or all
/// samples collinear with the first–last chord).
pub(crate) fn is_line_like(c: &dyn Curve) -> bool {
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
pub(crate) fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
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

pub(crate) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

pub(crate) fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
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
pub(crate) fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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
pub(crate) fn is_circle_like(c: &dyn Curve) -> bool {
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
pub(crate) fn circle_geometry(c: &dyn Curve) -> Option<(GpPnt, f64, GpVec)> {
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
pub(crate) fn curve_kind(c: &dyn Curve) -> CurveKind {
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
pub(crate) fn curve_resolution(c: &dyn Curve, tol: f64) -> f64 {
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

pub(crate) fn surface_sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

pub(crate) fn quadric_u0(s: &dyn Surface) -> f64 {
    let (u0, _) = s.u_range();
    if u0.is_finite() {
        u0
    } else {
        0.0
    }
}

pub(crate) fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// `(origin, unit-u, unit-v)` frame of a planar surface (natural parameterization).
pub(crate) fn plane_frame(s: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let (o, x, y) = s.d1(0.0, 0.0);
    let mx = x.magnitude();
    let my = y.magnitude();
    if mx < 1e-30 || my < 1e-30 {
        return None;
    }
    Some((o, x.divided(mx), y.divided(my)))
}

/// Exact projection of `p` onto a planar surface: `(u, v, distance)`.
pub(crate) fn plane_projection(s: &dyn Surface, p: &GpPnt) -> (f64, f64, f64) {
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
pub(crate) fn plane_normal(s: &dyn Surface) -> GpVec {
    match plane_frame(s) {
        Some((_, x, y)) => x.crossed(&y).normalized(),
        None => GpVec::zero(),
    }
}

/// `(axis point, axis direction, radius)` of a cylindrical surface.
pub(crate) fn cylinder_geometry(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
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

pub(crate) fn cylinder_matches(s: &dyn Surface, a: &GpPnt, z: &GpVec, r: f64) -> bool {
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
pub(crate) fn cone_geometry(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
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

pub(crate) fn cone_matches(s: &dyn Surface, apex: &GpPnt, z: &GpVec, tan_alpha: f64) -> bool {
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
pub(crate) fn torus_geometry(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64, f64)> {
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

pub(crate) fn torus_matches(s: &dyn Surface, o: &GpPnt, z: &GpVec, major: f64, minor: f64) -> bool {
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
pub(crate) fn surface_kind(s: &dyn Surface) -> SurfaceKind {
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
pub(crate) fn is_eq_distance(p: &GpPnt, s: &dyn Surface) -> Option<f64> {
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
pub(crate) fn dirs_coinside(d1: &GpVec, d2: &GpVec) -> bool {
    let d = d1.subtracted(d2).magnitude();
    d < 0.0002 || (2.0 - d).abs() < 0.0002
}

/// Whether a circle-like curve is coplanar with a planar surface.
pub(crate) fn is_coplanar(curve: &dyn Curve, surface: &dyn Surface) -> bool {
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
pub(crate) fn is_radius(curve: &dyn Curve, surface: &dyn Surface, criteria: f64) -> bool {
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
