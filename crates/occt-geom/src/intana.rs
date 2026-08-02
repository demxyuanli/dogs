//! Analytic intersection of primitives. Port of the `IntAna` package
//! (TKGeomBase): `IntAna_Int3Pln`, `IntAna_QuadQuadGeo` (the analytic
//! quadric-quadric cases: plane-plane, plane-sphere, sphere-sphere,
//! plane-cylinder, plane-cone) and `IntAna_IntLinTorus`.
//!
//! Shape-level intersection (`IntCurvesFace`, TKTopAlgo) is deliberately out
//! of scope; this module is the geometric kernel only.

use std::cmp::Ordering;

use occt_core::elib::clib;
use occt_core::gp::dir::DirAxis;
use occt_core::gp::{
    GpAx2, GpCirc, GpCone, GpCylinder, GpDir, GpElips, GpHypr, GpLin, GpParab, GpPln, GpPnt,
    GpSphere, GpTorus, GpVec,
};
use occt_core::precision::CONFUSION;

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
    Ellipse(GpElips),
    Parabola(GpParab),
    Hyperbola(GpHypr),
    /// The two quadrics coincide (same plane/sphere).
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
            (Ellipse(a), Ellipse(b)) => a == b,
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
        _ => QuadricIntersection::None,
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
}
