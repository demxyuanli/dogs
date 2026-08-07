//! Analytic intersection of primitives. Port of the `IntAna` package
//! (TKGeomBase): `IntAna_Int3Pln`, `IntAna_QuadQuadGeo` (the analytic
//! quadric-quadric cases: plane-plane, plane-sphere, sphere-sphere,
//! plane-cylinder, plane-cone, cylinder-cylinder, cylinder-sphere,
//! sphere-cone, cone-cone, cylinder-cone) and `IntAna_IntLinTorus`.
//!
//! Shape-level intersection (`IntCurvesFace`, TKTopAlgo) is deliberately out
//! of scope; this module is the geometric kernel only.

use std::cmp::Ordering;
use std::f64::consts::PI;

use occt_core::elib::clib;
use occt_core::gp::dir::DirAxis;
use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpDir2d, GpElips, GpHypr, GpLin,
    GpParab, GpPln, GpPnt, GpPnt2d, GpSphere, GpTorus, GpVec, GpVec2d,
};
use occt_core::precision::{ANGULAR, CONFUSION};

// ---------------------------------------------------------------------------
// Polynomial root helpers (self-contained; `occt-geom` does not depend on
// `occt-math`). Ports `math_DirectPolynomialRoots`.
// ---------------------------------------------------------------------------

fn quadratic_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < 1e-300 {
        if b.abs() < 1e-300 {
            return Vec::new();
        }
        return vec![-c / b];
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return Vec::new();
    }
    let sq = disc.sqrt();
    let q = -0.5 * (b + b.signum() * sq);
    let r1 = q / a;
    let r2 = if q.abs() > 1e-300 { c / q } else { (-b - sq) / (2.0 * a) };
    let mut v = vec![r1, r2];
    v.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    v
}

fn cubic_roots(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < 1e-300 {
        return quadratic_roots(b, c, d);
    }
    let b = b / a;
    let c = c / a;
    let d = d / a;
    let p = c - b * b / 3.0;
    let q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
    let disc = q * q / 4.0 + p * p * p / 27.0;
    let shift = b / 3.0;
    let mut roots = Vec::new();
    if disc > 1e-14 {
        let sq = disc.sqrt();
        let u = (-q / 2.0 + sq).cbrt();
        let v = (-q / 2.0 - sq).cbrt();
        roots.push(u + v - shift);
    } else if disc >= -1e-14 {
        if q.abs() < 1e-14 {
            roots.push(-shift);
        } else {
            let u = (-q / 2.0).cbrt();
            roots.push(2.0 * u - shift);
            roots.push(-u - shift);
        }
    } else {
        let r = 2.0 * (-p / 3.0).sqrt();
        let theta = ((3.0 * q / (2.0 * p)) * (-3.0 / p).sqrt()).acos() / 3.0;
        let two_pi = 2.0 * std::f64::consts::PI;
        roots.push(r * theta.cos() - shift);
        roots.push(r * (theta - two_pi / 3.0).cos() - shift);
        roots.push(r * (theta + two_pi / 3.0).cos() - shift);
    }
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    let mut out: Vec<f64> = Vec::new();
    for x in roots {
        if out.last().map_or(true, |&l| (x - l).abs() > 1e-9 * (1.0 + x.abs())) {
            out.push(x);
        }
    }
    out
}

fn poly4(a: f64, b: f64, c: f64, d: f64, e: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x * x + d * x + e
}

fn polish4(a: f64, b: f64, c: f64, d: f64, e: f64, mut x: f64) -> f64 {
    for _ in 0..10 {
        let fp = (4.0 * a * x + 3.0 * b) * x * x + 2.0 * c * x + d;
        if fp.abs() < 1e-300 {
            break;
        }
        let xn = x - poly4(a, b, c, d, e, x) / fp;
        if !xn.is_finite() || (xn - x).abs() > 1.0 + x.abs() {
            break;
        }
        if (xn - x).abs() < 1e-13 * (1.0 + x.abs()) {
            return xn;
        }
        x = xn;
    }
    x
}

/// Real roots of a·x⁴ + b·x³ + c·x² + d·x + e = 0 (Ferrari).
fn quartic_roots(a4: f64, a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    if a4.abs() < 1e-300 {
        return cubic_roots(a3, a2, a1, a0);
    }
    let a = a3 / a4;
    let b = a2 / a4;
    let c = a1 / a4;
    let d = a0 / a4;
    let p = b - 3.0 * a * a / 8.0;
    let q = c - a * b / 2.0 + a * a * a / 8.0;
    let r = d - a * c / 4.0 + a * a * b / 16.0 - 3.0 * a * a * a * a / 256.0;
    let ms = cubic_roots(8.0, -4.0 * p, -8.0 * r, 4.0 * p * r - q * q);
    let mut m: Option<f64> = None;
    for mm in ms {
        if 2.0 * mm - p >= 0.0 && m.map_or(true, |best| mm > best) {
            m = Some(mm);
        }
    }
    let m = match m {
        Some(m) => m,
        None => return Vec::new(),
    };
    let two_m_p = 2.0 * m - p;
    let mut yroots: Vec<f64> = Vec::new();
    if two_m_p < 1e-12 {
        let rad = m * m - r;
        if rad < 0.0 {
            return Vec::new();
        }
        let sr = rad.sqrt();
        for ym in [sr, -sr] {
            let val = -m + ym;
            if val >= 0.0 {
                let y = val.sqrt();
                yroots.push(y);
                if y > 1e-12 {
                    yroots.push(-y);
                }
            }
        }
    } else {
        let s = two_m_p.sqrt();
        let t = q / (2.0 * s);
        yroots.extend(quadratic_roots(1.0, -s, m + t));
        yroots.extend(quadratic_roots(1.0, s, m - t));
    }
    let mut out: Vec<f64> = Vec::new();
    for y in yroots {
        let x = y - a / 4.0;
        if poly4(1.0, a, b, c, d, x).abs() < 1e-6 {
            out.push(polish4(1.0, a, b, c, d, x));
        }
    }
    out.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    let mut dedup: Vec<f64> = Vec::new();
    for x in out {
        if dedup.last().map_or(true, |&l| (x - l).abs() > 1e-8 * (1.0 + x.abs())) {
            dedup.push(x);
        }
    }
    dedup
}

// ---------------------------------------------------------------------------
// Plane helpers.
// ---------------------------------------------------------------------------

/// Plane coefficients in OCCT form: `A·x + B·y + C·z + D = 0` with
/// (A,B,C) the unit normal (`gp_Pln::Coefficients`).
fn plane_coeffs(p: &GpPln) -> (f64, f64, f64, f64) {
    let n = p.pos.direction();
    let loc = p.location();
    let (a, b, c) = (n.x(), n.y(), n.z());
    let d = -(a * loc.x() + b * loc.y() + c * loc.z());
    (a, b, c, d)
}

/// A direction perpendicular to `z` (needed to build a coordinate frame).
fn perp_x_dir(z: &GpDir) -> GpDir {
    let base = if z.x().abs() < 0.9 {
        GpDir::from_axis(DirAxis::X)
    } else {
        GpDir::from_axis(DirAxis::Y)
    };
    base.cross(z).unwrap_or_else(|_| base)
}

/// Build an `Ax2` whose Z is `z` and X is `x` (falling back to a
/// perpendicular if `x` is degenerate).
fn ax2_from_dirs(origin: GpPnt, z: GpDir, x: GpDir) -> GpAx2 {
    GpAx2::new(origin, z, x).unwrap_or_else(|_| GpAx2::new(origin, z, perp_x_dir(&z)).unwrap())
}

// ---------------------------------------------------------------------------
// IntAna_Int3Pln
// ---------------------------------------------------------------------------

/// Result of intersecting three planes.
#[derive(Debug, Clone, PartialEq)]
pub enum Intersect3Pln {
    /// The three planes meet at a single point.
    Point(GpPnt),
    /// Two independent planes coincide and the third cuts them: common line.
    Line(GpLin),
    /// No common point/line (parallel or coincident degeneracy).
    NoResult,
}

fn solve3(m: [[f64; 3]; 3], rhs: [f64; 3]) -> Option<[f64; 3]> {
    let mut a = m;
    let mut b = rhs;
    for col in 0..3 {
        let mut piv = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(piv, col);
        b.swap(piv, col);
        let d = a[col][col];
        for r in (col + 1)..3 {
            let f = a[r][col] / d;
            for c in col..3 {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut s = b[i];
        for j in (i + 1)..3 {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

/// The intersection line of two non-parallel planes (port of the
/// least-norm solution `p = (d1·(n2×v) + d2·(v×n1))/|v|²`, `v = n1×n2`,
/// for `n·p = d`).
fn plane_plane_line(p1: &GpPln, p2: &GpPln) -> Option<GpLin> {
    let (a1, b1, c1, d1) = plane_coeffs(p1);
    let (a2, b2, c2, d2) = plane_coeffs(p2);
    let n1 = GpVec::new(a1, b1, c1);
    let n2 = GpVec::new(a2, b2, c2);
    // Work with the offset form n·p = d.
    let (d1, d2) = (-d1, -d2);
    let v = n1.crossed(&n2);
    let v2 = v.square_magnitude();
    if v2 < 1e-20 {
        return None;
    }
    let p0 = n2
        .crossed(&v)
        .multiplied_scalar(d1)
        .added(&v.crossed(&n1).multiplied_scalar(d2))
        .divided(v2);
    let dir = GpDir::from_xyz(&v.xyz()).ok()?;
    Some(GpLin::from_pnt_dir(GpPnt::from_xyz(&p0.coord), dir))
}

/// Intersection of three planes. Mirrors `IntAna_Int3Pln` (3×3 linear solve;
/// `NoResult` when the system is singular) and extends it with the common-line
/// case the OCCT class leaves as empty.
pub fn three_planes_intersect(p1: &GpPln, p2: &GpPln, p3: &GpPln) -> Result<Intersect3Pln, String> {
    let (a1, b1, c1, d1) = plane_coeffs(p1);
    let (a2, b2, c2, d2) = plane_coeffs(p2);
    let (a3, b3, c3, d3) = plane_coeffs(p3);
    let m = [
        [a1, b1, c1],
        [a2, b2, c2],
        [a3, b3, c3],
    ];
    // n·p + D = 0  ⇒  M·p = -D
    let rhs = [-d1, -d2, -d3];
    if let Some(x) = solve3(m, rhs) {
        return Ok(Intersect3Pln::Point(GpPnt::new(x[0], x[1], x[2])));
    }
    // Singular: try the pair with two independent normals.
    let pairs = [(p1, p2, p3), (p1, p3, p2), (p2, p3, p1)];
    for (pa, pb, pc) in pairs {
        if let Some(line) = plane_plane_line(pa, pb) {
            // Third plane must contain a representative point of the line.
            let q = clib::line_value(&line, 0.0);
            let (a, b, c, d) = plane_coeffs(pc);
            if (a * q.x() + b * q.y() + c * q.z() + d).abs() < 1e-9 {
                return Ok(Intersect3Pln::Line(line));
            }
        }
    }
    Ok(Intersect3Pln::NoResult)
}

// ---------------------------------------------------------------------------
// IntAna_QuadQuadGeo: plane-plane, plane-sphere, sphere-sphere,
// plane-cylinder, plane-cone.
// ---------------------------------------------------------------------------

/// Analytic intersection curve of two quadrics.
#[derive(Debug, Clone)]
pub enum QuadricIntersection {
    Point(GpPnt),
    Line(GpLin),
    TwoLines(GpLin, GpLin),
    Circle(GpCirc),
    /// Two circles (cylinder∩sphere, sphere∩cone, cone×cone on one axis…).
    TwoCircles(GpCirc, GpCirc),
    Ellipse(GpElips),
    /// Two ellipses (equal-radius cylinders with intersecting axes).
    TwoEllipses(GpElips, GpElips),
    Parabola(GpParab),
    Hyperbola(GpHypr),
    /// The two quadrics coincide (same plane/sphere/cone).
    Same,
    /// No intersection.
    None,
}

impl PartialEq for QuadricIntersection {
    fn eq(&self, other: &Self) -> bool {
        use QuadricIntersection::*;
        match (self, other) {
            (None, None) | (Same, Same) => true,
            (Point(a), Point(b)) => a == b,
            (Line(a), Line(b)) => a == b,
            (TwoLines(a1, a2), TwoLines(b1, b2)) => a1 == b1 && a2 == b2,
            (Circle(a), Circle(b)) => a.location() == b.location() && a.radius() == b.radius(),
            (TwoCircles(a1, a2), TwoCircles(b1, b2)) => {
                a1.location() == b1.location() && a1.radius() == b1.radius()
                    && a2.location() == b2.location() && a2.radius() == b2.radius()
            }
            (Ellipse(a), Ellipse(b)) => a == b,
            (TwoEllipses(a1, a2), TwoEllipses(b1, b2)) => a1 == b1 && a2 == b2,
            (Parabola(a), Parabola(b)) => a == b,
            (Hyperbola(a), Hyperbola(b)) => a == b,
            _ => false,
        }
    }
}

/// A quadric accepted by [`quadric_quadric`].
#[derive(Debug, Clone)]
pub enum Quadric {
    Plane(GpPln),
    Sphere(GpSphere),
    Cylinder(GpCylinder),
    Cone(GpCone),
}

/// Plane ∩ plane. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Pln)`.
pub fn quadric_quadric_planes(p1: &GpPln, p2: &GpPln, tol_ang: f64, tol: f64) -> QuadricIntersection {
    let (a1, b1, c1, d1) = plane_coeffs(p1);
    let (a2, b2, c2, d2) = plane_coeffs(p2);
    let n1 = GpVec::new(a1, b1, c1);
    let n2 = GpVec::new(a2, b2, c2);
    let vd = n1.crossed(&n2);
    let loc1 = p1.location();
    let loc2 = p2.location();
    let dist1 = a2 * loc1.x() + b2 * loc1.y() + c2 * loc1.z() + d2;
    let dist2 = a1 * loc2.x() + b1 * loc2.y() + c1 * loc2.z() + d1;
    let mvd = vd.magnitude();
    if mvd <= tol_ang {
        // Normals collinear: planes identical or parallel.
        if dist1.abs() <= tol && dist2.abs() <= tol {
            QuadricIntersection::Same
        } else {
            QuadricIntersection::None
        }
    } else {
        match plane_plane_line(p1, p2) {
            Some(l) => QuadricIntersection::Line(l),
            None => QuadricIntersection::None,
        }
    }
}

/// Plane ∩ sphere. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Sphere)`.
pub fn quadric_quadric_plane_sphere(p: &GpPln, s: &GpSphere) -> QuadricIntersection {
    let (a, b, c, d) = plane_coeffs(p);
    let loc = s.location();
    let radius = s.radius();
    let dist = a * loc.x() + b * loc.y() + c * loc.z() + d;
    if (dist.abs() - radius).abs() < radius.abs() * f64::EPSILON {
        // Tangent: single point, the projection of the center onto the plane.
        let pt = GpPnt::new(loc.x() - dist * a, loc.y() - dist * b, loc.z() - dist * c);
        QuadricIntersection::Point(pt)
    } else if dist.abs() < radius {
        let center = GpPnt::new(loc.x() - dist * a, loc.y() - dist * b, loc.z() - dist * c);
        let mut dir1 = *p.axis().direction();
        if !p.is_direct() {
            dir1.reverse();
        }
        let dir2 = *p.position().x_direction();
        let r = (radius * radius - dist * dist).sqrt();
        let pos = ax2_from_dirs(center, dir1, dir2);
        QuadricIntersection::Circle(GpCirc::new(pos, r))
    } else {
        QuadricIntersection::None
    }
}

/// Sphere ∩ sphere. Port of `IntAna_QuadQuadGeo::Perform(gp_Sphere, gp_Sphere)`.
pub fn quadric_quadric_sphere_sphere(s1: &GpSphere, s2: &GpSphere, tol: f64) -> QuadricIntersection {
    let o1 = s1.location();
    let o2 = s2.location();
    let d = o1.distance(&o2);
    let r1 = s1.radius();
    let r2 = s2.radius();
    let (rmin, rmax) = if r1 > r2 { (r2, r1) } else { (r1, r2) };

    if d <= tol && (r1 - r2).abs() <= tol {
        return QuadricIntersection::Same;
    }
    if d <= tol {
        return QuadricIntersection::None;
    }
    let dir = GpDir::from_vec(&GpVec::from_pnts(&o1, &o2)).unwrap_or_default();
    let t = rmax - d - rmin;
    if t >= 0.0 && t <= tol {
        let t2 = if r1 == rmax {
            (r1 + (r2 + d)) * 0.5
        } else {
            (-r1 + (d - r2)) * 0.5
        };
        let pt = o1.translated_vec(&GpVec::from_xyz(dir.xyz()).multiplied_scalar(t2));
        return QuadricIntersection::Point(pt);
    }
    if d > (r1 + r2 + tol) || rmax > (d + rmin + tol) {
        return QuadricIntersection::None;
    }
    let mut alpha = 0.5 * (r1 * r1 - r2 * r2 + d * d) / d;
    let mut beta = (r1 * r1 - alpha * alpha).max(0.0).sqrt();
    let eps_circle = 0.01 * CONFUSION;
    if beta <= eps_circle {
        alpha = (r1 + (d - r2)) * 0.5;
        let pt = o1.translated_vec(&GpVec::from_xyz(dir.xyz()).multiplied_scalar(alpha));
        QuadricIntersection::Point(pt)
    } else {
        if beta <= 0.0 {
            beta = 0.0;
        }
        let center = o1.translated_vec(&GpVec::from_xyz(dir.xyz()).multiplied_scalar(alpha));
        let pos = ax2_from_dirs(center, dir, perp_x_dir(&dir));
        QuadricIntersection::Circle(GpCirc::new(pos, beta))
    }
}

/// Result of a line-plane test (port of `IntAna_IntConicQuad::Perform(L, P)`).
enum LinePlaneHit {
    /// The line crosses the plane: parameter on the line and the point.
    Cross { param: f64, point: GpPnt },
    /// The line is parallel to the plane.
    Parallel,
}

fn line_plane(l: &GpLin, p: &GpPln, tol_ang: f64) -> LinePlaneHit {
    let (a, b, c, d) = plane_coeffs(p);
    let o = l.location();
    let dir = l.direction();
    let (al, bl, cl) = (dir.x(), dir.y(), dir.z());
    let direc = a * al + b * bl + c * cl;
    let dis = a * o.x() + b * o.y() + c * o.z() + d;
    if direc.abs() < tol_ang {
        LinePlaneHit::Parallel
    } else {
        let param = -dis / direc;
        LinePlaneHit::Cross {
            param,
            point: o.translated_vec(&GpVec::new(al * param, bl * param, cl * param)),
        }
    }
}

/// Plane ∩ cylinder. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cylinder)`
/// (H = 0, so the near-parallel direction-refinement branch is skipped).
pub fn quadric_quadric_plane_cylinder(
    p: &GpPln,
    cyl: &GpCylinder,
    tol_ang: f64,
    tol: f64,
) -> QuadricIntersection {
    let (a, b, c, d) = plane_coeffs(p);
    let radius = cyl.radius();
    let axec = cyl.axis();
    let axis_origin = *axec.location();
    let axis_dir = *axec.direction();
    let (x, y, z) = (axis_origin.x(), axis_origin.y(), axis_origin.z());
    // Signed distance from the axis origin to the plane.
    let dist = a * x + b * y + c * z + d;
    let normp = p.pos.direction();

    // Parallel detection: |axis · plane normal| below the angular tolerance.
    let direc = (a * axis_dir.x() + b * axis_dir.y() + c * axis_dir.z()).abs();
    if direc < tol_ang {
        // Axis parallel to the plane → the section is one or two lines.
        let omega = GpPnt::new(x - dist * a, y - dist * b, z - dist * c);
        if (dist.abs() - radius).abs() < tol {
            // Tangent: a single line.
            QuadricIntersection::Line(GpLin::from_pnt_dir(omega, axis_dir))
        } else if dist.abs() < radius {
            let h = (radius * radius - dist * dist).sqrt();
            let axey = axis_dir.xyz().crossed(&normp.xyz());
            let p1 = omega.translated_vec(&GpVec::from_xyz(&axey.multiplied(-h)));
            let p2 = omega.translated_vec(&GpVec::from_xyz(&axey.multiplied(h)));
            QuadricIntersection::TwoLines(
                GpLin::from_pnt_dir(p1, axis_dir),
                GpLin::from_pnt_dir(p2, axis_dir),
            )
        } else {
            QuadricIntersection::None
        }
    } else {
        // The axis crosses the plane; the section is a circle (axis ⊥ plane)
        // or an ellipse.
        let lin = GpLin::from_pnt_dir(axis_origin, axis_dir);
        let center = match line_plane(&lin, p, tol_ang) {
            LinePlaneHit::Cross { point, .. } => point,
            LinePlaneHit::Parallel => return QuadricIntersection::None,
        };
        let axey = normp.xyz().crossed(&axis_dir.xyz());
        let sint = axey.modulus();
        if sint < tol / radius {
            // Perpendicular axis → circle (axes of the cylinder).
            let cdir = cyl.position().direction();
            let xdir = *cyl.position().x_direction();
            let pos = ax2_from_dirs(center, cdir, xdir);
            QuadricIntersection::Circle(GpCirc::new(pos, radius))
        } else {
            // Oblique cut → ellipse.
            let cost = axis_dir.dot(&normp).abs();
            let axex = axey.crossed(&normp.xyz());
            let pos = ax2_from_dirs(center, normp, GpDir::from_xyz(&axex).unwrap_or_default());
            QuadricIntersection::Ellipse(GpElips::new(pos, radius / cost, radius))
        }
    }
}

/// Plane ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cone)`.
pub fn quadric_quadric_plane_cone(
    p: &GpPln,
    cone: &GpCone,
    tol_ang: f64,
    tol: f64,
) -> QuadricIntersection {
    let (a, b, c, d) = plane_coeffs(p);
    let axec = cone.axis();
    let axis_dir = *axec.direction();
    let apex = cone.apex();
    let dist = a * apex.x() + b * apex.y() + c * apex.z() + d;
    let mut normp = p.pos.direction();
    if !p.is_direct() {
        normp.reverse();
    }
    let axey = normp.xyz().crossed(&axis_dir.xyz());
    let axex = axey.crossed(&normp.xyz());
    let angl = cone.semi_angle();
    let cosa = angl.cos();
    let sina = angl.sin().abs();
    let sint = axey.modulus();
    let cost = axis_dir.dot(&normp).abs();
    // sin((π/2 − t) − angl) = cos(t + angl); zero ⇒ plane contains a generatrix.
    let costa = cost * cosa - sint * sina;

    if dist.abs() < tol {
        // Plane contains the apex: point, one or two lines.
        if costa.abs() < tol_ang {
            // Parallel to a generatrix → one line through the apex.
            let mut ptonaxe = apex.coord.added(&axis_dir.xyz().multiplied(10.0));
            let d2 = a * ptonaxe.x + b * ptonaxe.y + c * ptonaxe.z + d;
            ptonaxe = ptonaxe.subtracted(&normp.xyz().multiplied(d2));
            let dir = GpDir::from_xyz(&ptonaxe.subtracted(&apex.coord)).unwrap_or_default();
            QuadricIntersection::Line(GpLin::from_pnt_dir(apex, dir))
        } else if cost < sina {
            // Plane interior to the cone → two generatrices.
            let dh = (sina * sina - cost * cost).sqrt() / cosa;
            let d1 = GpDir::from_xyz(&axex.added(&axey.multiplied(dh))).unwrap_or_default();
            let d2 = GpDir::from_xyz(&axex.subtracted(&axey.multiplied(dh))).unwrap_or_default();
            QuadricIntersection::TwoLines(
                GpLin::from_pnt_dir(apex, d1),
                GpLin::from_pnt_dir(apex, d2),
            )
        } else {
            // Plane exterior → apex only.
            QuadricIntersection::Point(apex)
        }
    } else if cost < tol_ang {
        // Plane contains the axis direction → hyperbola.
        let pt1 = apex.coord.subtracted(&normp.xyz().multiplied(dist));
        let major = (dist / angl.tan()).abs();
        let minor = dist.abs();
        let pos = ax2_from_dirs(
            GpPnt::from_xyz(&pt1),
            normp,
            GpDir::from_xyz(&axex).unwrap_or_default(),
        );
        QuadricIntersection::Hyperbola(GpHypr::new(pos, major, minor))
    } else {
        // The cone axis crosses the plane at `center`.
        let lin = GpLin::from_pnt_dir(*axec.location(), axis_dir);
        let (center, param_on_axis) = match line_plane(&lin, p, tol_ang) {
            LinePlaneHit::Cross { point, param } => (point, param),
            LinePlaneHit::Parallel => return QuadricIntersection::None,
        };
        let distance = apex.distance(&center);
        let mut axex = axex;
        let mut axey = axey;
        // The apex parameter on the axis is negative in the OCCT cone frame;
        // flip the section axes when the center is on the apex side.
        if param_on_axis + cone.radius() / angl.tan() < 0.0 {
            axex.reverse();
            axey.reverse();
        }

        if costa.abs() < tol_ang {
            // Parallel to a generatrix → parabola.
            let deltacenter = distance / 2.0 / cosa;
            let axex_n = GpDir::from_xyz(&axex.normalized()).unwrap_or_default();
            let pt1 = center.coord.subtracted(&axex_n.xyz().multiplied(deltacenter));
            let focal = deltacenter * sina * sina;
            let pos = ax2_from_dirs(GpPnt::from_xyz(&pt1), normp, axex_n);
            QuadricIntersection::Parabola(GpParab::new(pos, focal))
        } else if sint < tol_ang {
            // Plane perpendicular to the axis → circle.
            let r = distance * angl.tan().abs();
            let cdir = cone.position().direction();
            let xdir = *cone.position().x_direction();
            let pos = ax2_from_dirs(center, cdir, xdir);
            QuadricIntersection::Circle(GpCirc::new(pos, r))
        } else if cost < sina {
            // Hyperbola.
            let deltacenter = sint * sina * sina * distance / (sina * sina - cost * cost);
            let axex_n = GpDir::from_xyz(&axex.normalized()).unwrap_or_default();
            let pt1 = center.coord.subtracted(&axex_n.xyz().multiplied(deltacenter));
            let major = cost * sina * cosa * distance / (sina * sina - cost * cost);
            let minor = cost * sina * distance / (sina * sina - cost * cost).sqrt();
            let pos = ax2_from_dirs(GpPnt::from_xyz(&pt1), normp, axex_n);
            QuadricIntersection::Hyperbola(GpHypr::new(pos, major, minor))
        } else {
            // cost > sina → ellipse.
            let radius = cost * sina * cosa * distance / (cost * cost - sina * sina);
            let deltacenter = sint * sina * sina * distance / (cost * cost - sina * sina);
            let axex_n = GpDir::from_xyz(&axex.normalized()).unwrap_or_default();
            let pt1 = center.coord.added(&axex_n.xyz().multiplied(deltacenter));
            let minor = cost * sina * distance / (cost * cost - sina * sina).sqrt();
            let pos = ax2_from_dirs(GpPnt::from_xyz(&pt1), normp, axex_n);
            QuadricIntersection::Ellipse(GpElips::new(pos, radius, minor))
        }
    }
}

// ---------------------------------------------------------------------------
// IntAna_QuadQuadGeo: cylinder×cylinder, cylinder×sphere, sphere×cone,
// cone×cone, cylinder×cone.
// ---------------------------------------------------------------------------

/// Snap an axis direction that is nearly axis-aligned to the exact axis.
/// Port of `RefineDir` (IntAna_QuadQuadGeo.cxx): it keeps the cross products
/// used by the axis relation code away from exact degeneracies.
fn refine_dir(d: &mut GpDir) {
    let mut c = [d.x(), d.y(), d.z()];
    let (mut m, mut n) = (0, 0);
    for &v in &c {
        if v == 1.0 || v == -1.0 {
            m += 1;
        } else if v != 0.0 {
            n += 1;
        }
    }
    if m > 0 && n > 0 {
        let eps = f64::EPSILON;
        let (r1, r2) = (1.0 - eps, 1.0 + eps);
        for k in 0..3 {
            let num = c[k].abs();
            if num > r1 && num < r2 {
                c[k] = if c[k] > 0.0 { 1.0 } else { -1.0 };
                c[(k + 1) % 3] = 0.0;
                c[(k + 2) % 3] = 0.0;
                break;
            }
        }
        if let Ok(dd) = GpDir::new(c[0], c[1], c[2]) {
            *d = dd;
        }
    }
}

/// Scalar triple product `a·(b×c)` — OCCT `Det33` in AxeOperator.
fn det33(a: &GpVec, b: &GpVec, c: &GpVec) -> f64 {
    a.dot(&b.crossed(c))
}

/// Relation between two axes — port of `AxeOperator`
/// (IntAna_QuadQuadGeo.cxx). Returns `(parallel, distance, coplanar,
/// intersection point)`; `coplanar` needs the axes to (nearly) meet
/// (`distance < eps_dist` and the triple product within tolerance), and the
/// point is only computed for concurrent non-parallel axes.
fn axe_operator(a1: &GpAx1, a2: &GpAx1, eps_dist: f64, eps_para: f64) -> (bool, f64, bool, Option<GpPnt>) {
    let mut v1 = *a1.direction();
    let mut v2 = *a2.direction();
    refine_dir(&mut v1);
    refine_dir(&mut v2);
    let (p1, p2) = (*a1.location(), *a2.location());
    let (w1, w2) = (GpVec::from_xyz(v1.xyz()), GpVec::from_xyz(v2.xyz()));
    let parallel = w1.cross_magnitude(&w2) <= eps_para;
    let distance = if parallel {
        let rel = GpVec::from_pnts(&p1, &p2);
        rel.subtracted(&w1.multiplied_scalar(rel.dot(&w1))).magnitude()
    } else {
        w1.crossed(&w2).normalized().dot(&GpVec::from_pnts(&p1, &p2)).abs()
    };
    let mut coplanar = false;
    let mut pt_intersect = None;
    if distance < eps_dist {
        let det = det33(&w1, &w2, &GpVec::from_pnts(&p2, &p1)); // rows V1, V2, P1−P2
        if det.abs() <= eps_dist {
            coplanar = true;
            if !parallel {
                // Concurrent axes: intersection point P1 + A·V1.
                let sm = GpVec::from_pnts(&p1, &p2);
                let d1 = w1.y() * w2.x() - w1.x() * w2.y();
                let d2 = w1.z() * w2.y() - w1.y() * w2.z();
                let d3 = w1.z() * w2.x() - w1.x() * w2.z();
                let a = if d1 != 0.0 && d1.abs() >= d2.abs() && d1.abs() >= d3.abs() {
                    (sm.y() * w2.x() - sm.x() * w2.y()) / d1
                } else if d2 != 0.0 && d2.abs() >= d1.abs() && d2.abs() >= d3.abs() {
                    (sm.z() * w2.y() - sm.y() * w2.z()) / d2
                } else {
                    (sm.z() * w2.x() - sm.x() * w2.z()) / d3
                };
                pt_intersect = Some(p1.translated_vec(&w1.multiplied_scalar(a)));
            }
        }
    }
    (parallel, distance, coplanar, pt_intersect)
}

/// The perpendicular segment between two non-parallel axes: signed distance and
/// the parameters of the closest points on each axis. Port of
/// `AxeOperator::Distance`.
fn axe_distance(a1: &GpAx1, a2: &GpAx1) -> (f64, f64, f64) {
    let w1 = GpVec::from_xyz(a1.direction().xyz());
    let w2 = GpVec::from_xyz(a2.direction().xyz());
    let o1o2 = GpVec::from_pnts(a1.location(), a2.location());
    let n = w1.crossed(&w2);
    if n.magnitude() < 1e-12 {
        return (0.0, 0.0, 0.0);
    }
    let n = n.normalized();
    let d = det33(&w1, &w2, &n);
    if d != 0.0 {
        let dist = det33(&w1, &w2, &o1o2) / d;
        let p1 = det33(&o1o2, &w2, &n) / (-d);
        let p2 = det33(&w1, &o1o2, &n) / d;
        (dist, p1, p2)
    } else {
        (0.0, 0.0, 0.0)
    }
}

/// Distance from `p` to the axis line.
fn dist_point_axis(p: &GpPnt, ax: &GpAx1) -> f64 {
    let dir = GpVec::from_xyz(ax.direction().xyz());
    let rel = GpVec::from_pnts(ax.location(), p);
    rel.subtracted(&dir.multiplied_scalar(rel.dot(&dir))).magnitude()
}

fn mid_pnt(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// A circle with the given center / plane normal / radius (the frame X
/// direction is arbitrary — it does not affect the curve).
fn circle_with_normal(center: GpPnt, normal: GpDir, radius: f64) -> GpCirc {
    GpCirc::new(ax2_from_dirs(center, normal, perp_x_dir(&normal)), radius)
}

/// A plane through `pt` with the given normal.
fn plane_normal_at(pt: GpPnt, normal: GpDir) -> Option<GpPln> {
    let ax3 = GpAx3::new(pt, normal, &perp_x_dir(&normal)).ok()?;
    Some(GpPln::new(ax3))
}

/// Cylinder ∩ cylinder. Port of `IntAna_QuadQuadGeo::Perform(gp_Cylinder,
/// gp_Cylinder)`: parallel axes give two generatrix lines (or one when tangent,
/// `Same`/`None` for coincident/disjoint), intersecting equal-radius axes give
/// the two bisector-plane ellipses, external tangency a point, and everything
/// else (`NoGeometricSolution`) → `None` for the numeric walker.
pub fn quadric_quadric_cylinder_cylinder(
    c1: &GpCylinder,
    c2: &GpCylinder,
    tol: f64,
) -> QuadricIntersection {
    let (parallel, dist, coplanar, pt_inter) =
        axe_operator(&c1.axis(), &c2.axis(), CONFUSION, ANGULAR);
    let r1 = c1.radius();
    let r2 = c2.radius();
    let rmr = (r1 - r2).abs();
    let rmr_relative = rmr / r1.max(r2);
    let dir_cyl = c1.position().direction();
    let wdir = GpVec::from_xyz(dir_cyl.xyz());

    if parallel {
        if dist <= tol {
            return if rmr <= tol {
                QuadricIntersection::Same
            } else {
                QuadricIntersection::None
            };
        }
        // Parallel axes, strictly separated. Project the 2nd location onto the
        // 1st cylinder base plane and intersect the two base circles.
        let p1 = c1.location();
        let p2t = c2.location();
        let proj = wdir.dot(&GpVec::from_pnts(&p1, &p2t));
        let p2 = p2t.translated_vec(&wdir.multiplied_scalar(-proj));
        let r1p2 = r1 + r2;
        if dist > r1p2 + tol {
            QuadricIntersection::None
        } else if (r1p2 - dist) <= f64::EPSILON {
            // External tangency: one generatrix line.
            let pt1 = p1.translated_vec(&GpVec::from_pnts(&p1, &p2).multiplied_scalar(r1 / r1p2));
            QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
        } else if dist > rmr {
            // Two generatrix lines (or one when the base circles are tangent).
            let a_r1r1 = r1 * r1;
            let a_cos = 0.5 * (a_r1r1 - r2 * r2 + dist * dist) / (r1 * dist);
            let a_sin2 = 1.0 - a_cos * a_cos;
            let is_tangent = 4.0 * a_r1r1 * a_sin2 < tol * tol;
            let dir_a1a2 = GpVec::from_pnts(&p1, &p2).divided(dist);
            if is_tangent {
                let pt1 = p1.translated_vec(&dir_a1a2.multiplied_scalar(r1 * a_cos));
                QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
            } else {
                let a_sin = a_sin2.sqrt();
                let axd = *c1.position().x_direction();
                let ayd = *c1.position().y_direction();
                let r1x = GpVec::from_xyz(axd.xyz()).multiplied_scalar(r1);
                let r1y = GpVec::from_xyz(ayd.xyz()).multiplied_scalar(r1);
                let adx = dir_a1a2.dot(&GpVec::from_xyz(axd.xyz()));
                let ady = dir_a1a2.dot(&GpVec::from_xyz(ayd.xyz()));
                let (ndx, ndy) = (adx * a_cos - ady * a_sin, ady * a_cos + adx * a_sin);
                let pt1 = p1
                    .translated_vec(&r1x.multiplied_scalar(ndx))
                    .translated_vec(&r1y.multiplied_scalar(ndy));
                let (ndx, ndy) = (adx * a_cos + ady * a_sin, ady * a_cos - adx * a_sin);
                let pt2 = p1
                    .translated_vec(&r1x.multiplied_scalar(ndx))
                    .translated_vec(&r1y.multiplied_scalar(ndy));
                QuadricIntersection::TwoLines(
                    GpLin::from_pnt_dir(pt1, dir_cyl),
                    GpLin::from_pnt_dir(pt2, dir_cyl),
                )
            }
        } else if dist > rmr - tol {
            // Internal tangency: one generatrix line.
            let mut r1_rmr = r1 / rmr;
            if r1 < r2 {
                r1_rmr = -r1_rmr;
            }
            let pt1 = p1.translated_vec(&GpVec::from_pnts(&p1, &p2).multiplied_scalar(r1_rmr));
            QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
        } else {
            QuadricIntersection::None
        }
    } else if rmr_relative <= 1e-13 && coplanar {
        // Equal-radius cylinders with intersecting axes → two ellipses in the
        // bisector planes. Frame: `dir1`/`dir2` are the bisectors (perpendicular),
        // each used as the plane normal of one ellipse.
        let Some(pt) = pt_inter else { return QuadricIntersection::None };
        let wd2 = GpVec::from_xyz(c2.position().direction().xyz());
        let ang = wdir.angle(&wd2);
        let b = (0.5 * (PI - ang)).sin().abs();
        let a = (0.5 * ang).sin().abs();
        if a == 0.0 || b == 0.0 {
            return QuadricIntersection::Same;
        }
        let Ok(d1) = GpDir::from_vec(&wdir.added(&wd2)) else { return QuadricIntersection::None };
        let Ok(d2) = GpDir::from_vec(&wdir.subtracted(&wd2)) else { return QuadricIntersection::None };
        let (mut p1, mut p1bis) = (r1 / b, r1);
        let (mut p2, mut p2bis) = (r1 / a, r1);
        if p1 < p1bis {
            std::mem::swap(&mut p1, &mut p1bis);
        }
        if p2 < p2bis {
            std::mem::swap(&mut p2, &mut p2bis);
        }
        let Ok(ax2_1) = GpAx2::new(pt, d1, d2) else { return QuadricIntersection::None };
        let Ok(ax2_2) = GpAx2::new(pt, d2, d1) else { return QuadricIntersection::None };
        QuadricIntersection::TwoEllipses(
            GpElips::new(ax2_1, p1, p1bis),
            GpElips::new(ax2_2, p2, p2bis),
        )
    } else if (dist - r1 - r2).abs() < tol {
        // External tangency with intersecting (non-parallel) axes: a point on
        // the common perpendicular.
        let d1 = *c1.axis().direction();
        let d2 = *c2.axis().direction();
        let (_, p1p, p2p) = axe_distance(&c1.axis(), &c2.axis());
        let p1 = c1
            .axis()
            .location()
            .translated_vec(&GpVec::from_xyz(d1.xyz()).multiplied_scalar(-p1p));
        let p2 = c2
            .axis()
            .location()
            .translated_vec(&GpVec::from_xyz(d2.xyz()).multiplied_scalar(-p2p));
        let Ok(dir) = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)) else { return QuadricIntersection::None };
        let pt = p1.translated_vec(&GpVec::from_xyz(dir.xyz()).multiplied_scalar(r1));
        QuadricIntersection::Point(pt)
    } else {
        QuadricIntersection::None
    }
}

/// Cylinder ∩ sphere. Port of `IntAna_QuadQuadGeo::Perform(gp_Cylinder,
/// gp_Sphere)`: when the sphere center lies on the cylinder axis the section is
/// one or two circles (radius = cylinder radius, centered at
/// `center ± √(r_sph² − r_cyl²)·axis`); otherwise `NoGeometricSolution` → `None`.
pub fn quadric_quadric_cylinder_sphere(
    cyl: &GpCylinder,
    sph: &GpSphere,
    _tol: f64,
) -> QuadricIntersection {
    let pt = sph.location();
    // OCCT tests the axes to intersect at the sphere center exactly; a small
    // tolerance keeps the closed form on near-axis configurations.
    if dist_point_axis(&pt, &cyl.axis()) > 1e-9 {
        return QuadricIntersection::None;
    }
    let r_cyl = cyl.radius();
    let r_sph = sph.radius();
    if r_sph < r_cyl {
        return QuadricIntersection::None; // IntAna_Empty
    }
    let dist = (r_sph * r_sph - r_cyl * r_cyl).sqrt();
    let dir = cyl.position().direction();
    let w = GpVec::from_xyz(dir.xyz());
    let c1 = pt.translated_vec(&w.multiplied_scalar(dist));
    let circ1 = circle_with_normal(c1, dir, r_cyl);
    if dist > f64::EPSILON {
        let c2 = pt.translated_vec(&w.multiplied_scalar(-dist));
        QuadricIntersection::TwoCircles(circ1, circle_with_normal(c2, dir, r_cyl))
    } else {
        QuadricIntersection::Circle(circ1)
    }
}

/// Sphere ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Sphere, gp_Cone)`:
/// when the sphere center lies on the cone axis the section is one or two
/// circles — the roots of the 2D cross-section quadratic
/// `(1+tg²)x² + 2·tg²·d·x + tg²·d² − r² = 0`, with `d` the apex→center
/// distance. Otherwise `NoGeometricSolution` → `None`.
pub fn quadric_quadric_sphere_cone(
    sph: &GpSphere,
    cone: &GpCone,
    _tol: f64,
) -> QuadricIntersection {
    let pt = sph.location();
    if dist_point_axis(&pt, &cone.axis()) > 1e-9 {
        return QuadricIntersection::None;
    }
    let apex = cone.apex();
    let d = pt.distance(&apex);
    let condir = if d > f64::EPSILON {
        let Ok(c) = GpDir::from_vec(&GpVec::from_pnts(&apex, &pt)) else {
            return QuadricIntersection::None;
        };
        c
    } else {
        cone.position().direction()
    };
    let rad = sph.radius();
    let tga = cone.semi_angle().tan();
    let tgatga = tga * tga;
    let roots = quadratic_roots(1.0 + tgatga, 2.0 * tgatga * d, -rad * rad + d * d * tgatga);
    if roots.is_empty() {
        return QuadricIntersection::None; // IntAna_Empty
    }
    let w = GpVec::from_xyz(condir.xyz());
    let mut circles: Vec<GpCirc> = Vec::new();
    for x in roots {
        let dpx = d + x;
        let center = apex.translated_vec(&w.multiplied_scalar(dpx));
        let r = (tga * dpx).abs();
        if r <= 0.01 * CONFUSION {
            continue; // IntAna_PointAndCircle: degenerate radius → point
        }
        circles.push(circle_with_normal(center, condir, r));
    }
    match circles.len() {
        1 => QuadricIntersection::Circle(circles.pop().unwrap()),
        2 => {
            let c2 = circles.pop().unwrap();
            let c1 = circles.pop().unwrap();
            QuadricIntersection::TwoCircles(c1, c2)
        }
        _ => QuadricIntersection::None,
    }
}

/// Cone ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cone)` for the
/// tractable branches: coincident axes (two circles / a point / `Same`),
/// parallel axes with equal semi-angle (a conic in the plane through the two
/// apexes, from the plane∩cone closed form), and coincident apexes (one or two
/// generatrix lines). The common-generatrix case and everything else return
/// `None` for the numeric walker.
pub fn quadric_quadric_cone_cone(
    c1: &GpCone,
    c2: &GpCone,
    tol_ang: f64,
    tol: f64,
) -> QuadricIntersection {
    let tg1 = c1.semi_angle().tan();
    let mut tg2 = c2.semi_angle().tan();
    if tg1 * tg2 < 0.0 {
        tg2 = -tg2;
    }
    let tol2 = tol * tol;
    let ap1 = c1.apex();
    let ap2 = c2.apex();
    let d_a1a2 = ap1.square_distance(&ap2);
    let (parallel, dist_axes, _coplanar, _pt_inter) =
        axe_operator(&c1.axis(), &c2.axis(), 1e-14, ANGULAR);

    // 1 — coincident axes: two circles where the cone radii match (or the two
    // cones coincide / touch at the apex).
    if parallel && dist_axes < 1e-14 {
        let p = c1.apex();
        let d = c1.position().direction();
        let w = GpVec::from_xyz(d.xyz());
        let offset = w.dot(&GpVec::from_pnts(&p, &ap2));
        if (tg1 - tg2).abs() > ANGULAR {
            if offset.abs() < 1e-10 {
                return QuadricIntersection::Point(p);
            }
            let x1 = offset * tg2 / (tg1 + tg2);
            let x2 = offset * tg2 / (tg2 - tg1);
            let c1c = p.translated_vec(&w.multiplied_scalar(x1));
            let c2c = p.translated_vec(&w.multiplied_scalar(x2));
            QuadricIntersection::TwoCircles(
                circle_with_normal(c1c, d, (x1 * tg1).abs()),
                circle_with_normal(c2c, d, (x2 * tg1).abs()),
            )
        } else if offset.abs() < 1e-10 {
            QuadricIntersection::Same
        } else {
            let x = 0.5 * offset;
            QuadricIntersection::Circle(circle_with_normal(
                p.translated_vec(&w.multiplied_scalar(x)),
                d,
                (x * tg1).abs(),
            ))
        }
    }
    // 2 — parallel axes with (nearly) equal semi-angles: the intersection lies
    // in the plane through the two apexes; reduce to the plane∩cone conic.
    else if (tg1 - tg2).abs() < tol_ang && parallel {
        let da1 = c1.position().direction();
        let o1o2 = GpVec::from_pnts(&ap1, &ap2);
        let o1o2n = o1o2.normalized();
        let o1o2_da1 = GpVec::from_xyz(da1.xyz()).dot(&o1o2n);
        let o1_proj = o1o2n.subtracted(&GpVec::from_xyz(da1.xyz()).multiplied_scalar(o1o2_da1));
        let Ok(db1) = GpDir::from_vec(&o1_proj) else { return QuadricIntersection::None };
        let y_o1o2 = o1o2.dot(&GpVec::from_xyz(da1.xyz()));
        let abstg1 = tg1.abs();
        let x2 = (dist_axes / abstg1 - y_o1o2) * 0.5;
        let x1 = x2 + y_o1o2;
        let p1 = ap1
            .translated_vec(&GpVec::from_xyz(da1.xyz()).multiplied_scalar(x1))
            .translated_vec(&GpVec::from_xyz(db1.xyz()).multiplied_scalar(x1 * abstg1));
        let p1_m = GpVec::from_pnts(&p1, &mid_pnt(&ap1, &ap2));
        let da1_x_db1 = GpVec::from_xyz(da1.xyz()).crossed(&GpVec::from_xyz(db1.xyz()));
        let ortho = da1_x_db1.crossed(&p1_m);
        let Ok(n) = GpDir::from_vec(&ortho) else { return QuadricIntersection::None };
        let Some(pln) = plane_normal_at(p1, n) else { return QuadricIntersection::None };
        quadric_quadric_plane_cone(&pln, c1, tol_ang, tol)
    }
    // 3 — coincident apexes: one or two generatrix lines (or `None`).
    else if d_a1a2 < tol2 {
        cone_cone_common_apex(c1, c2, tg1, tg2, tol)
    } else {
        // 4/5 — common generatrix / general: no analytic closed form.
        QuadricIntersection::None
    }
}

/// Coincident-apex cones: `IntAna_QuadQuadGeo::Perform` branch 3 — a 2D
/// section analysis determines touch/intersection, then the one or two
/// generatrix lines are built through the shared apex.
fn cone_cone_common_apex(c1: &GpCone, c2: &GpCone, tg1: f64, tg2: f64, tol: f64) -> QuadricIntersection {
    let half_pi = 0.5 * PI;
    let d1 = 1.0;
    let p0 = GpPnt2d::new(0.0, 0.0);
    let ax1 = c1.axis();
    let ax2 = c2.axis();
    let mut gamma = ax1.direction().angle(ax2.direction());
    if gamma > half_pi {
        gamma = PI - gamma;
    }
    let (cos_g, sin_g) = (gamma.cos(), gamma.sin());
    let tg_beta1 = tg1.abs();
    let tg_beta2 = tg2.abs();
    let r1 = d1 * tg_beta1;
    let p1 = GpPnt2d::new(d1, r1);
    // Project P1 onto the 2nd axis line (in the plane of the two axes) to find
    // whether the section circles overlap, touch or miss.
    let v_ax2 = GpVec2d::new(cos_g, sin_g);
    let Ok(_) = GpDir2d::from_vec2d(&v_ax2) else { return QuadricIntersection::None };
    let v = GpVec2d::new(p1.x() - p0.x(), p1.y() - p0.y());
    let mut dx = v_ax2.dot(&v);
    let pa2 = p0.translated_vec(&v_ax2.multiplied_scalar(dx));
    dx = pa2.distance(&p0);
    let r2 = dx * tg_beta2;
    let rd2 = pa2.distance(&p1);
    if rd2 > r2 + tol {
        return QuadricIntersection::None; // IntAna_Empty
    }
    let i_ret = if rd2 < r2 - tol { 2 } else { 1 };
    // 3D construction: two planes perpendicular to the axes through the ring
    // points Q1/Q2 intersect in the line through the section mid-point QX.
    let q_apex1 = c1.apex();
    let d3_ax1 = *ax1.direction();
    let w1 = GpVec::from_xyz(d3_ax1.xyz());
    let qa1 = q_apex1.translated_vec(&w1.multiplied_scalar(d1));
    let dx = w1.dot(&GpVec::from_xyz(ax2.direction().xyz()));
    let d3_ax2 = if dx < 0.0 {
        ax2.direction().reversed()
    } else {
        *ax2.direction()
    };
    let w2 = GpVec::from_xyz(d3_ax2.xyz());
    let d2 = d1 * ((1.0 + tg_beta1 * tg_beta1) / (1.0 + tg_beta2 * tg_beta2)).sqrt();
    let qa2 = q_apex1.translated_vec(&w2.multiplied_scalar(d2));
    let Some(pln1) = plane_normal_at(qa1, d3_ax1) else { return QuadricIntersection::None };
    let Some(pln2) = plane_normal_at(qa2, d3_ax2) else { return QuadricIntersection::None };
    let Some(lin) = plane_plane_line(&pln1, &pln2) else { return QuadricIntersection::None };
    let wl = GpVec::from_xyz(lin.direction().xyz());
    let orig = lin.location();
    let vr = GpVec::from_pnts(&qa1, &orig);
    let dx = wl.dot(&vr);
    let qx = orig.translated_vec(&wl.multiplied_scalar(dx));
    if i_ret == 1 {
        // One tangency line.
        let Ok(dir) = GpDir::from_vec(&GpVec::from_pnts(&q_apex1, &qx)) else {
            return QuadricIntersection::None;
        };
        QuadricIntersection::Line(GpLin::from_pnt_dir(q_apex1, dir))
    } else {
        // Two intersection lines.
        let da = qa1.distance(&qx);
        let ddx = (r1 * r1 - da * da).sqrt();
        let qx1 = qx.translated_vec(&wl.multiplied_scalar(ddx));
        let qx2 = qx.translated_vec(&wl.multiplied_scalar(-ddx));
        let Ok(dir1) = GpDir::from_vec(&GpVec::from_pnts(&q_apex1, &qx1)) else {
            return QuadricIntersection::None;
        };
        let Ok(dir2) = GpDir::from_vec(&GpVec::from_pnts(&q_apex1, &qx2)) else {
            return QuadricIntersection::None;
        };
        QuadricIntersection::TwoLines(
            GpLin::from_pnt_dir(q_apex1, dir1),
            GpLin::from_pnt_dir(q_apex1, dir2),
        )
    }
}

/// Cylinder ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Cylinder,
/// gp_Cone)`: only the coincident-axis case has a closed form (two circles at
/// `apex ± r_cyl/tan(angle)` along the axis); otherwise `NoGeometricSolution`.
pub fn quadric_quadric_cylinder_cone(
    cyl: &GpCylinder,
    cone: &GpCone,
    _tol: f64,
) -> QuadricIntersection {
    let (parallel, dist, _, _) = axe_operator(&cyl.axis(), &cone.axis(), 1e-14, ANGULAR);
    if !(parallel && dist < 1e-14) {
        return QuadricIntersection::None;
    }
    let pt = cone.apex();
    let dist = cyl.radius() / cone.semi_angle().tan();
    let dir = cyl.position().direction();
    let w = GpVec::from_xyz(dir.xyz());
    let r = cyl.radius();
    QuadricIntersection::TwoCircles(
        circle_with_normal(pt.translated_vec(&w.multiplied_scalar(dist)), dir, r),
        circle_with_normal(pt.translated_vec(&w.multiplied_scalar(-dist)), dir, r),
    )
}

/// Analytic quadric-quadric intersection dispatcher for the tractable pairs.
pub fn quadric_quadric(q1: &Quadric, q2: &Quadric, tol_ang: f64, tol: f64) -> QuadricIntersection {
    use Quadric::*;
    match (q1, q2) {
        (Plane(a), Plane(b)) => quadric_quadric_planes(a, b, tol_ang, tol),
        (Plane(p), Sphere(s)) | (Sphere(s), Plane(p)) => quadric_quadric_plane_sphere(p, s),
        (Sphere(a), Sphere(b)) => quadric_quadric_sphere_sphere(a, b, tol),
        (Plane(p), Cylinder(c)) | (Cylinder(c), Plane(p)) => {
            quadric_quadric_plane_cylinder(p, c, tol_ang, tol)
        }
        (Plane(p), Cone(c)) | (Cone(c), Plane(p)) => quadric_quadric_plane_cone(p, c, tol_ang, tol),
        (Cylinder(a), Cylinder(b)) => quadric_quadric_cylinder_cylinder(a, b, tol),
        (Cylinder(c), Sphere(s)) | (Sphere(s), Cylinder(c)) => quadric_quadric_cylinder_sphere(c, s, tol),
        (Sphere(s), Cone(c)) | (Cone(c), Sphere(s)) => quadric_quadric_sphere_cone(s, c, tol),
        (Cone(a), Cone(b)) => quadric_quadric_cone_cone(a, b, tol_ang, tol),
        (Cylinder(c), Cone(k)) | (Cone(k), Cylinder(c)) => quadric_quadric_cylinder_cone(c, k, tol),
    }
}

// ---------------------------------------------------------------------------
// IntAna_IntLinTorus
// ---------------------------------------------------------------------------

/// Intersection of a line with a torus. Port of `IntAna_IntLinTorus::Perform`:
/// a quartic in the line parameter in the torus reference frame, with each
/// root verified by re-evaluation on the torus.
pub fn line_torus_intersect(l: &GpLin, t: &GpTorus) -> Vec<GpPnt> {
    let pl = l.location();
    let dl = l.direction();
    let tor_loc = t.location();
    // Reparametrize so the line location is nearest the torus location.
    let param_of_new_pl = GpVec::from_pnts(&pl, &tor_loc).dot(&GpVec::from_xyz(dl.xyz()));
    let new_pl = pl.translated_vec(&GpVec::from_xyz(dl.xyz()).multiplied_scalar(param_of_new_pl));

    // Express the line in the torus reference frame.
    let pos = t.position();
    let (xd, yd, zd) = (
        GpVec::from_xyz(pos.x_direction().xyz()),
        GpVec::from_xyz(pos.y_direction().xyz()),
        GpVec::from_xyz(pos.direction().xyz()),
    );
    let v = GpVec::from_pnts(&pos.location(), &new_pl);
    let (x0, y0, z0) = (v.dot(&xd), v.dot(&yd), v.dot(&zd));
    let dv = GpVec::from_xyz(dl.xyz());
    let (x1, y1, z1) = (dv.dot(&xd), dv.dot(&yd), dv.dot(&zd));

    let r = t.major_radius();
    let r2 = r * r;
    let rr = t.minor_radius();
    let rr2 = rr * rr;

    let a = x1 * x1 + y1 * y1 + z1 * z1;
    let b = 2.0 * (x1 * x0 + y1 * y0 + z1 * z0);
    let c = x0 * x0 + y0 * y0 + z0 * z0 - (r2 + rr2);

    let a4 = a * a;
    let a3 = 2.0 * a * b;
    let a2 = 2.0 * a * c + 4.0 * r2 * z1 * z1 + b * b;
    let a1 = 2.0 * b * c + 8.0 * r2 * z1 * z0;
    let a0 = c * c + 4.0 * r2 * (z0 * z0 - rr2);

    let mut out = Vec::new();
    let mut seen: Vec<GpPnt> = Vec::new();
    for mut tt in quartic_roots(a4, a3, a2, a1, a0) {
        tt += param_of_new_pl;
        let p = clib::line_value(l, tt);
        // Verify the point lies on the torus (OCCT re-checks the square distance).
        let v2 = GpVec::from_pnts(&tor_loc, &p);
        let dz = v2.dot(&zd);
        let perp = v2.coord.subtracted(&zd.xyz().multiplied(dz));
        let rho = perp.modulus();
        let err = ((rho - r) * (rho - r) + dz * dz - rr2).abs();
        if err < 1e-7 {
            if seen.iter().all(|s| s.distance(&p) > 1e-7) {
                seen.push(p);
                out.push(p);
            }
        }
    }
    out.sort_by(|a, b| a.x().partial_cmp(&b.x()).unwrap_or(Ordering::Equal));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpAx3;
    use occt_core::precision::ANGULAR;

    const PI: f64 = std::f64::consts::PI;

    fn plane(origin: GpPnt, normal: GpDir) -> GpPln {
        GpPln::new(
            GpAx3::new(origin, normal, &perp_x_dir(&normal)).unwrap_or_default(),
        )
    }

    fn sphere_at(origin: GpPnt, r: f64) -> GpSphere {
        let ax3 = GpAx3::new(origin, GpDir::from_axis(DirAxis::Z), &GpDir::from_axis(DirAxis::X))
            .unwrap_or_default();
        GpSphere::new(ax3, r).unwrap()
    }

    fn on_plane(p: &GpPnt, pl: &GpPln, tol: f64) -> bool {
        let (a, b, c, d) = plane_coeffs(pl);
        (a * p.x() + b * p.y() + c * p.z() + d).abs() < tol
    }

    fn on_cone(p: &GpPnt, cone: &GpCone, tol: f64) -> bool {
        let apex = cone.apex();
        let axis = *cone.axis().direction();
        let v = GpVec::from_pnts(&apex, p);
        let h = v.dot(&GpVec::from_xyz(axis.xyz()));
        let perp = v.coord.subtracted(&axis.xyz().multiplied(h));
        let rho = perp.modulus();
        (rho - h.abs() * cone.semi_angle().tan()).abs() < tol
    }

    fn on_cylinder(p: &GpPnt, cyl: &GpCylinder, tol: f64) -> bool {
        let axis = *cyl.axis().direction();
        let v = GpVec::from_pnts(cyl.axis().location(), p);
        let h = v.dot(&GpVec::from_xyz(axis.xyz()));
        let perp = v.coord.subtracted(&axis.xyz().multiplied(h));
        (perp.modulus() - cyl.radius()).abs() < tol
    }

    fn on_sphere(p: &GpPnt, sph: &GpSphere, tol: f64) -> bool {
        (p.distance(&sph.location()) - sph.radius()).abs() < tol
    }

    fn cylinder_axis(center: GpPnt, axis: GpDir, r: f64) -> GpCylinder {
        let x = perp_x_dir(&axis);
        GpCylinder::new(GpAx3::new(center, axis, &x).unwrap_or_default(), r).unwrap()
    }

    /// A cone whose geometric apex is `apex` (radius 0 anchor).
    fn cone_apex_at(apex: GpPnt, axis: GpDir, semi: f64) -> GpCone {
        let x = perp_x_dir(&axis);
        GpCone::new(GpAx3::new(apex, axis, &x).unwrap_or_default(), 0.0, semi).unwrap()
    }

    #[test]
    fn three_planes_known_point() {
        // x=1, y=2, z=3 meet at (1,2,3).
        let p1 = plane(GpPnt::new(1.0, 0.0, 0.0), GpDir::from_axis(DirAxis::X));
        let p2 = plane(GpPnt::new(0.0, 2.0, 0.0), GpDir::from_axis(DirAxis::Y));
        let p3 = plane(GpPnt::new(0.0, 0.0, 3.0), GpDir::from_axis(DirAxis::Z));
        match three_planes_intersect(&p1, &p2, &p3).unwrap() {
            Intersect3Pln::Point(p) => {
                assert!((p.x() - 1.0).abs() < 1e-9);
                assert!((p.y() - 2.0).abs() < 1e-9);
                assert!((p.z() - 3.0).abs() < 1e-9, "got {p:?}");
            }
            other => panic!("expected point, got {other:?}"),
        }
    }

    #[test]
    fn three_planes_parallel_no_result() {
        let p1 = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let p2 = plane(GpPnt::new(0.0, 0.0, 1.0), GpDir::from_axis(DirAxis::Z));
        let p3 = plane(GpPnt::new(0.0, 0.0, 2.0), GpDir::from_axis(DirAxis::Z));
        assert_eq!(three_planes_intersect(&p1, &p2, &p3).unwrap(), Intersect3Pln::NoResult);
    }

    #[test]
    fn plane_sphere_circle_z0() {
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let s = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        match quadric_quadric_plane_sphere(&pl, &s) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9, "center {:?}", c.location());
                assert!((c.radius() - 1.0).abs() < 1e-9, "radius {}", c.radius());
                // Every sampled point lies in z=0 at distance 1 from the origin.
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!(on_plane(&p, &pl, 1e-9), "{p:?} not in plane");
                    assert!((p.distance(&GpPnt::new(0., 0., 0.)) - 1.0).abs() < 1e-9);
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_sphere_no_intersection() {
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let s = sphere_at(GpPnt::new(0.0, 0.0, 5.0), 1.0);
        assert_eq!(quadric_quadric_plane_sphere(&pl, &s), QuadricIntersection::None);
    }

    #[test]
    fn sphere_sphere_circle() {
        // Two unit spheres, centers 1 apart: circle radius √(1 − (d/2)²) = √0.75
        // at the midpoint (0.5, 0, 0).
        let s1 = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        let s2 = sphere_at(GpPnt::new(1.0, 0.0, 0.0), 1.0);
        match quadric_quadric_sphere_sphere(&s1, &s2, 1e-7) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0.5, 0.0, 0.0)) < 1e-9, "center {:?}", c.location());
                assert!((c.radius() - (0.75f64).sqrt()).abs() < 1e-9, "radius {}", c.radius());
                // Every point on the circle is on both spheres.
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!((p.distance(&GpPnt::new(0., 0., 0.)) - 1.0).abs() < 1e-9);
                    assert!((p.distance(&GpPnt::new(1., 0., 0.)) - 1.0).abs() < 1e-9);
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn sphere_sphere_disjoint() {
        let s1 = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        let s2 = sphere_at(GpPnt::new(5.0, 0.0, 0.0), 1.0);
        assert_eq!(quadric_quadric_sphere_sphere(&s1, &s2, 1e-7), QuadricIntersection::None);
    }

    #[test]
    fn plane_cylinder_circle_perpendicular_axis() {
        // Cylinder along Z at origin r=1, plane z=0 → circle radius 1 at origin.
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let cyl = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        match quadric_quadric_plane_cylinder(&pl, &cyl, ANGULAR, CONFUSION) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
                assert!((c.radius() - 1.0).abs() < 1e-9);
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!(on_plane(&p, &pl, 1e-9));
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_cylinder_two_lines_parallel_axis() {
        // Cylinder along X, base at (0,1,0.5), r=1; plane z=0 → two lines at
        // y = 1 ± √(1 − 0.5²) in z=0.
        let xdir = GpDir::from_axis(DirAxis::X);
        let zdir = GpDir::from_axis(DirAxis::Z);
        let ax3 = GpAx3::new(GpPnt::new(0.0, 1.0, 0.5), xdir, &zdir).unwrap();
        let cyl = GpCylinder::new(ax3, 1.0).unwrap();
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), zdir);
        match quadric_quadric_plane_cylinder(&pl, &cyl, ANGULAR, CONFUSION) {
            QuadricIntersection::TwoLines(l1, l2) => {
                let h = (0.75f64).sqrt();
                let p1 = clib::line_value(&l1, 0.0);
                let p2 = clib::line_value(&l2, 0.0);
                for p in [p1, p2] {
                    assert!(on_plane(&p, &pl, 1e-9), "{p:?} not in plane");
                    assert!((p.y() - (1.0 + h)).abs() < 1e-9 || (p.y() - (1.0 - h)).abs() < 1e-9, "y {}", p.y());
                }
            }
            other => panic!("expected two lines, got {other:?}"),
        }
    }

    #[test]
    fn plane_cone_circle_perpendicular_axis() {
        // Cone apex at origin, axis +Z, semi-angle 30°; plane z=1 → circle of
        // radius tan(30°) at (0,0,1).
        let zdir = GpDir::from_axis(DirAxis::Z);
        let h = 1.0;
        let radius = h * (PI / 6.0).tan();
        let loc = GpPnt::new(0.0, 0.0, -h);
        let ax3 = GpAx3::new(loc, zdir, &GpDir::from_axis(DirAxis::X)).unwrap();
        let cone = GpCone::new(ax3, radius, PI / 6.0).unwrap();
        assert!(cone.apex().distance(&GpPnt::new(0., 0., 0.)) < 1e-9, "apex {:?}", cone.apex());
        let pl = plane(GpPnt::new(0.0, 0.0, 1.0), zdir);
        match quadric_quadric_plane_cone(&pl, &cone, ANGULAR, CONFUSION) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0., 0., 1.)) < 1e-9, "center {:?}", c.location());
                assert!((c.radius() - radius).abs() < 1e-9, "radius {} vs {}", c.radius(), radius);
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!(on_plane(&p, &pl, 1e-9));
                    assert!(on_cone(&p, &cone, 1e-6), "{p:?} not on cone");
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_cone_ellipse_lies_on_quadrics() {
        // Oblique plane (not through the apex, not perpendicular, not parallel
        // to a generatrix) cutting a cone → ellipse; verify membership by
        // sampling.
        let cone = GpCone::new(GpAx3::standard(), 1.0, PI / 6.0).unwrap();
        let pl = plane(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.3, 0.954).unwrap());
        match quadric_quadric_plane_cone(&pl, &cone, ANGULAR, CONFUSION) {
            QuadricIntersection::Ellipse(e) => {
                for k in 0..16 {
                    let p = clib::ellipse_value(&e, 2.0 * PI * k as f64 / 16.0);
                    assert!(on_plane(&p, &pl, 1e-6), "{p:?} not in plane");
                    assert!(on_cone(&p, &cone, 1e-5), "{p:?} not on cone");
                }
            }
            other => panic!("expected ellipse, got {other:?}"),
        }
    }

    #[test]
    fn line_torus_x_axis_four_points() {
        // Torus at origin R=3 r=1; line along X through origin cuts at x=±2, ±4.
        let torus = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let line = GpLin::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::from_axis(DirAxis::X));
        let pts = line_torus_intersect(&line, &torus);
        assert_eq!(pts.len(), 4, "got {pts:?}");
        let expected = [-4.0, -2.0, 2.0, 4.0];
        let mut got: Vec<f64> = pts.iter().map(|p| p.x()).collect();
        got.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (g, e) in got.iter().zip(expected.iter()) {
            assert!((g - e).abs() < 1e-6, "got {g} expected {e}: {pts:?}");
        }
    }

    #[test]
    fn line_torus_miss() {
        // Line along Z through the center misses (the minor circle is off-axis).
        let torus = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let line = GpLin::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::from_axis(DirAxis::Z));
        assert!(line_torus_intersect(&line, &torus).is_empty());
    }

    #[test]
    fn cylinder_cylinder_parallel_two_lines() {
        // Two unit cylinders, parallel Z axes through (0,0,0) and (0.5,0,0):
        // base circles intersect at x=0.25, y=±√(1−0.25²).
        let z = GpDir::from_axis(DirAxis::Z);
        let c1 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 1.0);
        let c2 = cylinder_axis(GpPnt::new(0.5, 0.0, 0.0), z, 1.0);
        match quadric_quadric_cylinder_cylinder(&c1, &c2, 1e-7) {
            QuadricIntersection::TwoLines(l1, l2) => {
                let h = (0.9375f64).sqrt();
                for l in [l1, l2] {
                    let p = clib::line_value(&l, 0.0);
                    assert!((p.y().abs() - h).abs() < 1e-9, "y {}", p.y());
                    assert!((p.x() - 0.25).abs() < 1e-9, "x {}", p.x());
                    assert!(on_cylinder(&p, &c1, 1e-9), "{p:?} not on cyl1");
                    assert!(on_cylinder(&p, &c2, 1e-9), "{p:?} not on cyl2");
                }
            }
            other => panic!("expected two lines, got {other:?}"),
        }
    }

    #[test]
    fn cylinder_cylinder_nested_empty() {
        // Concentric cylinders of different radii never meet.
        let z = GpDir::from_axis(DirAxis::Z);
        let c1 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 1.0);
        let c2 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 2.0);
        assert_eq!(
            quadric_quadric_cylinder_cylinder(&c1, &c2, 1e-7),
            QuadricIntersection::None
        );
    }

    #[test]
    fn cylinder_cylinder_intersecting_ellipses() {
        // Equal unit cylinders, perpendicular axes through the origin: two
        // bisector-plane ellipses (x²+z²=1 ∧ x²+y²=1 → y=±z).
        let z = GpDir::from_axis(DirAxis::Z);
        let x = GpDir::from_axis(DirAxis::X);
        let c1 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 1.0);
        let c2 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), x, 1.0);
        match quadric_quadric_cylinder_cylinder(&c1, &c2, 1e-7) {
            QuadricIntersection::TwoEllipses(e1, e2) => {
                assert!(e1.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
                assert!(e2.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
                for e in [e1, e2] {
                    for k in 0..16 {
                        let p = clib::ellipse_value(&e, 2.0 * PI * k as f64 / 16.0);
                        assert!(on_cylinder(&p, &c1, 1e-6), "{p:?} not on cyl1");
                        assert!(on_cylinder(&p, &c2, 1e-6), "{p:?} not on cyl2");
                    }
                }
            }
            other => panic!("expected two ellipses, got {other:?}"),
        }
    }

    #[test]
    fn cylinder_sphere_two_circles() {
        // Unit sphere at the origin, Z-axis cylinder r=0.5 through it: circles
        // at z=±√(1−0.25)=±0.866, radius 0.5.
        let z = GpDir::from_axis(DirAxis::Z);
        let cyl = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 0.5);
        let sph = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        match quadric_quadric_cylinder_sphere(&cyl, &sph, 1e-7) {
            QuadricIntersection::TwoCircles(c1, c2) => {
                for c in [c1, c2] {
                    assert!((c.radius() - 0.5).abs() < 1e-9, "radius {}", c.radius());
                    assert!(c.location().x().abs() < 1e-9 && c.location().y().abs() < 1e-9);
                    assert!((c.location().z().abs() - (0.75f64).sqrt()).abs() < 1e-9);
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cylinder(&p, &cyl, 1e-6), "{p:?} not on cylinder");
                        assert!(on_sphere(&p, &sph, 1e-6), "{p:?} not on sphere");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }

    #[test]
    fn sphere_cone_two_circles() {
        // Cone apex at the origin, axis +Z, 30°; sphere center (0,0,3), r=2:
        // two circles at z≈3.40 (r≈1.96) and z≈1.10 (r≈0.64).
        let z = GpDir::from_axis(DirAxis::Z);
        let cone = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 6.0);
        let sph = sphere_at(GpPnt::new(0.0, 0.0, 3.0), 2.0);
        match quadric_quadric_sphere_cone(&sph, &cone, 1e-7) {
            QuadricIntersection::TwoCircles(c1, c2) => {
                for c in [c1, c2] {
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cone(&p, &cone, 1e-5), "{p:?} not on cone");
                        assert!(on_sphere(&p, &sph, 1e-5), "{p:?} not on sphere");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }

    #[test]
    fn cone_cone_same_axis_two_circles() {
        // Same axis +Z; cone1 apex origin 30°, cone2 apex (0,0,1) 45°.
        let z = GpDir::from_axis(DirAxis::Z);
        let c1 = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 6.0);
        let c2 = cone_apex_at(GpPnt::new(0.0, 0.0, 1.0), z, PI / 4.0);
        match quadric_quadric_cone_cone(&c1, &c2, ANGULAR, 1e-7) {
            QuadricIntersection::TwoCircles(a, b) => {
                for c in [a, b] {
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cone(&p, &c1, 1e-5), "{p:?} not on cone1");
                        assert!(on_cone(&p, &c2, 1e-5), "{p:?} not on cone2");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }

    #[test]
    fn cone_cone_common_apex_two_lines() {
        // Two 30° cones sharing the apex, axes differing by 40°: two generatrix
        // lines through the apex.
        let z = GpDir::from_axis(DirAxis::Z);
        let axis2 = GpDir::new((40.0f64).to_radians().sin(), 0.0, (40.0f64).to_radians().cos()).unwrap();
        let c1 = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 6.0);
        let c2 = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), axis2, PI / 6.0);
        match quadric_quadric_cone_cone(&c1, &c2, ANGULAR, 1e-7) {
            QuadricIntersection::TwoLines(l1, l2) => {
                for l in [l1, l2] {
                    for t in [-2.0, -1.0, 1.0, 2.0] {
                        let p = clib::line_value(&l, t);
                        assert!(on_cone(&p, &c1, 1e-5), "{p:?} not on cone1");
                        assert!(on_cone(&p, &c2, 1e-5), "{p:?} not on cone2");
                    }
                }
            }
            other => panic!("expected two lines, got {other:?}"),
        }
    }

    #[test]
    fn cylinder_cone_same_axis_two_circles() {
        // Cylinder r=0.5 axis Z through origin; cone apex origin 45°: circles at
        // z=±0.5/tan(45°)=±0.5, radius 0.5.
        let z = GpDir::from_axis(DirAxis::Z);
        let cyl = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 0.5);
        let cone = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 4.0);
        match quadric_quadric_cylinder_cone(&cyl, &cone, 1e-7) {
            QuadricIntersection::TwoCircles(c1, c2) => {
                for c in [c1, c2] {
                    assert!((c.radius() - 0.5).abs() < 1e-9);
                    assert!(c.location().x().abs() < 1e-9 && c.location().y().abs() < 1e-9);
                    assert!((c.location().z().abs() - 0.5).abs() < 1e-9);
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cylinder(&p, &cyl, 1e-6), "{p:?} not on cylinder");
                        assert!(on_cone(&p, &cone, 1e-6), "{p:?} not on cone");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }
}
