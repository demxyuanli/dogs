//! Curve/surface sampling classifiers used by BeanFace analytic paths.
use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_geom::{Curve, Surface};



const PI: f64 = std::f64::consts::PI;


/// Coarse analytic kind of a curve (sampling classification).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FastCurveKind {
    Line,
    Circle,
    Ellipse,
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

/// First three non-collinear samples of a curve.
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

/// Circumcenter of three non-collinear points.
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

pub(crate) fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// Whether the curve is a full circle-like closed curve: all samples coplanar
/// and equidistant from a common centre.
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
pub(crate) fn circle_geometry(c: &dyn Curve) -> Option<(GpPnt, f64, GpDir)> {
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
    let d = GpDir::from_vec(&nrm.divided(m)).ok()?;
    Some((center, radius, d))
}

/// Whether the curve is an ellipse-like closed curve (coplanar, periodic 2π,
/// not a circle).
pub(crate) fn is_ellipse_like(c: &dyn Curve) -> bool {
    if is_circle_like(c) {
        return false;
    }
    if !c.is_periodic() {
        return false;
    }
    let period = c.period();
    if (period - 2.0 * PI).abs() > 1e-6 {
        return false;
    }
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    // All samples coplanar.
    let p0 = c.d0(a);
    let p1 = c.d0(a + (b - a) / 6.0);
    let p2 = c.d0(a + 2.0 * (b - a) / 6.0);
    let (p0, p1, p2) = match first_three_spanning(&[p0, p1, p2]) {
        Some(x) => x,
        None => return false,
    };
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return false;
    }
    let nv = nrm.divided(m);
    for i in 0..24 {
        let u = a + (b - a) * i as f64 / 24.0;
        let v = GpVec::from_pnts(&p0, &c.d0(u));
        if v.dot(&nv).abs() > 1e-5 * v.magnitude().max(1.0) {
            return false;
        }
    }
    true
}

/// `(center, axis direction)` of an ellipse-like curve (natural 2π
/// parametrisation: opposite samples share the centre, the plane normal is the
/// axis).
pub(crate) fn ellipse_geometry(c: &dyn Curve) -> Option<(GpPnt, GpDir)> {
    let center = midpoint(&c.d0(0.0), &c.d0(PI));
    let p0 = c.d0(0.0);
    let pq = c.d0(PI / 2.0);
    let nrm = GpVec::from_pnts(&center, &p0).crossed(&GpVec::from_pnts(&center, &pq));
    let m = nrm.magnitude();
    if m <= 1e-30 {
        return None;
    }
    let d = GpDir::from_vec(&nrm.divided(m)).ok()?;
    // Validate the centre with the perpendicular pair.
    let c2 = midpoint(&c.d0(PI / 2.0), &c.d0(3.0 * PI / 2.0));
    if c2.distance(&center) > 1e-4 * center.distance(&p0).max(1.0) {
        return None;
    }
    Some((center, d))
}

/// `(location, unit direction)` of a line-like curve.
pub(crate) fn line_geometry(c: &dyn Curve) -> Option<(GpPnt, GpDir)> {
    let p0 = c.d0(0.0);
    let p1 = c.d0(1.0);
    let v = GpVec::from_pnts(&p0, &p1);
    let m = v.magnitude();
    if m <= 1e-30 {
        return None;
    }
    let d = GpDir::from_vec(&v.divided(m)).ok()?;
    Some((p0, d))
}

pub(crate) fn classify_fast_curve(c: &dyn Curve) -> FastCurveKind {
    if is_line_like(c) {
        return FastCurveKind::Line;
    }
    if is_circle_like(c) {
        return FastCurveKind::Circle;
    }
    if is_ellipse_like(c) {
        return FastCurveKind::Ellipse;
    }
    FastCurveKind::Other
}

/// `(location, x_dir, y_dir, normal)` of a planar surface (natural frame).
pub(crate) fn plane_geometry(s: &dyn Surface) -> Option<(GpPnt, GpDir, GpDir, GpDir)> {
    let (u0, _) = s.u_range();
    let (v0, _) = s.v_range();
    let u0 = if u0.is_finite() { u0 } else { 0.0 };
    let v0 = if v0.is_finite() { v0 } else { 0.0 };
    let o = s.d0(u0, v0);
    let xv = GpVec::from_pnts(&o, &s.d0(u0 + 1.0, v0));
    let yv = GpVec::from_pnts(&o, &s.d0(u0, v0 + 1.0));
    let nv = xv.crossed(&yv);
    let n = GpDir::from_vec(&nv).ok()?;
    let x = GpDir::from_vec(&xv).ok()?;
    let y = GpDir::from_vec(&yv).ok()?;
    Some((o, x, y, n))
}

/// `(center, radius)` of a spherical surface.
pub(crate) fn sphere_geometry(s: &dyn Surface) -> Option<(GpPnt, f64)> {
    let center = crate::brep_surface::sphere_center(s)?;
    let (u0, _) = s.u_range();
    let u0 = if u0.is_finite() { u0 } else { 0.0 };
    let r = s.d0(u0, 0.0).distance(&center);
    if r < 1e-12 {
        return None;
    }
    Some((center, r))
}

/// `(axis, radius)` of a cylindrical surface (axis direction from the V
/// advance, axis point = circle centre at V start, radius = half the opposite
/// sample distance).
pub(crate) fn cylinder_geometry(s: &dyn Surface) -> Option<(occt_core::gp::GpAx1, f64)> {
    let u0 = 0.0;
    let c0 = midpoint(&s.d0(u0, 0.0), &s.d0(u0 + PI, 0.0));
    let c1 = midpoint(&s.d0(u0, 1.0), &s.d0(u0 + PI, 1.0));
    let z = GpVec::from_pnts(&c0, &c1);
    let zm = z.magnitude();
    if zm < 1e-12 {
        return None;
    }
    let zd = GpDir::from_vec(&z.divided(zm)).ok()?;
    let r = s.d0(u0, 0.0).distance(&c0);
    if r < 1e-12 {
        return None;
    }
    Some((occt_core::gp::GpAx1::new(c0, zd), r))
}

/// Signed distance from `p` to a plane (`plane_geometry` frame).
pub(crate) fn plane_distance(ploc: &GpPnt, nrm: &GpDir, p: &GpPnt) -> f64 {
    GpVec::from_pnts(ploc, p).dot(&GpVec::from_xyz(nrm.xyz())).abs()
}
