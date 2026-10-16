use super::prelude::*;
use super::IntAnaIntConicQuad;

use occt_math::{Gauss, MathMatrix, MathVector};


// ---------------------------------------------------------------------------
// Polynomial root helpers.
//
// These are the ported `math_DirectPolynomialRoots`
// (`occt-core/src/math_direct_poly_roots.rs`) consumed as the OCCT
// `Value(1 .. NbSolutions)` list. The closed-form quadratic/cubic and the
// Ferrari quartic that used to live here were a second, invented solver and
// are gone.
// ---------------------------------------------------------------------------

pub(crate) use occt_core::math_direct_poly_roots::{quadratic_roots, quartic_roots};

// ---------------------------------------------------------------------------
// Plane helpers.
// ---------------------------------------------------------------------------

/// Plane coefficients in OCCT form: `A·x + B·y + C·z + D = 0` with
/// (A,B,C) the unit normal (`gp_Pln::Coefficients`).
pub(super) fn plane_coeffs(p: &GpPln) -> (f64, f64, f64, f64) {
    let n = p.pos.direction();
    let loc = p.location();
    let (a, b, c) = (n.x(), n.y(), n.z());
    let d = -(a * loc.x() + b * loc.y() + c * loc.z());
    (a, b, c, d)
}

/// A direction perpendicular to `z` (needed to build a coordinate frame).
pub(super) fn perp_x_dir(z: &GpDir) -> GpDir {
    let base = if z.x().abs() < 0.9 {
        GpDir::from_axis(DirAxis::X)
    } else {
        GpDir::from_axis(DirAxis::Y)
    };
    base.cross(z).unwrap_or_else(|_| base)
}

/// Build an `Ax2` whose Z is `z` and X is `x` (falling back to a
/// perpendicular if `x` is degenerate).
pub(super) fn ax2_from_dirs(origin: GpPnt, z: GpDir, x: GpDir) -> GpAx2 {
    GpAx2::new(origin, z, x).unwrap_or_else(|_| GpAx2::new(origin, z, perp_x_dir(&z)).unwrap())
}

// ---------------------------------------------------------------------------
// IntAna_Int3Pln
// ---------------------------------------------------------------------------

/// Result of intersecting three planes. `IntAna_Int3Pln` answers a single
/// point or nothing at all (`IntAna_Int3Pln.hxx:29-70`): `NoResult` is
/// `IsEmpty()`.
#[derive(Debug, Clone, PartialEq)]
pub enum Intersect3Pln {
    /// The three planes meet at a single point.
    Point(GpPnt),
    /// `empt = true` (`IntAna_Int3Pln.cxx:53-56`): parallel/coincident planes.
    NoResult,
}

/// `IntAna_Int3Pln::Perform(P1, P2, P3)` (`IntAna_Int3Pln.cxx:40-64`).
///
/// `m` is `M(1..3, 1..3)` built from `P.Coefficients` and `d` is `V(1..3)`,
/// the `D` part of the same call (`cxx:44-46`). The system is solved by
/// `math_Gauss Resol(M, gp::Resolution())` (`cxx:51`): `IsDone()` false means
/// `empt = true`, which is the `None` returned here. Otherwise `V = -V` and
/// `Resol.Solve(V)` give the coordinates (`cxx:58-61`).
pub(super) fn solve3(m: [[f64; 3]; 3], d: [f64; 3]) -> Option<[f64; 3]> {
    let mut matrix = MathMatrix::new(1, 3, 1, 3);
    for r in 0..3 {
        for c in 0..3 {
            matrix.set_value(r + 1, c + 1, m[r][c]);
        }
    }
    let resol = Gauss::new(&matrix, REAL_SMALL);
    if !resol.is_done() {
        return None;
    }
    let mut v = MathVector::new(1, 3);
    for i in 0..3 {
        v.set_value(i + 1, -d[i]);
    }
    let x = resol.solve(&v)?;
    Some([x.value(1), x.value(2), x.value(3)])
}

/// `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Pln, TolAng, Tol)`
/// (`IntAna_QuadQuadGeo.cxx:389-511`), the line arm.
///
/// `None` covers `IntAna_Empty` (collinear normals, `cxx:415-418`) and every
/// early exit that leaves `done = false` (`cxx:463-509`); `Some` is the
/// `typeres = IntAna_Line` / `nbint = 1` result of `cxx:444-445` (`pt1` the
/// midpoint of the two normal-projected origins, `dir1 = gp_Dir(vd)`).
pub(super) fn plane_plane_line(p1: &GpPln, p2: &GpPln, tol_ang: f64, tol: f64) -> Option<GpLin> {
    let (a1, b1, c1, d1) = plane_coeffs(p1);
    let (a2, b2, c2, d2) = plane_coeffs(p2);
    let n1 = GpVec::new(a1, b1, c1);
    let n2 = GpVec::new(a2, b2, c2);
    let vd = n1.crossed(&n2);
    let loc1 = p1.location();
    let loc2 = p2.location();
    let dist1 = a2 * loc1.x() + b2 * loc1.y() + c2 * loc1.z() + d2;
    let dist2 = a1 * loc2.x() + b1 * loc2.y() + c1 * loc2.z() + d1;
    let a_mvd = vd.magnitude();
    if a_mvd <= tol_ang {
        // `cxx:415-418`: normals are collinear, the planes are same or parallel.
        return None;
    }
    // `cxx:423-431`: the `aEps = 1.e-16` clamp and the two parameters along
    // the normals.
    let denom = a1 * a2 + b1 * b2 + c1 * c2;
    let ddenom = 1.0 - denom * denom;
    let denom = if ddenom.abs() <= A_EPS { A_EPS } else { ddenom };
    let par1 = dist1 / denom;
    let par2 = -dist2 / denom;
    let inter1 = n1.crossed(&vd);
    let inter2 = n2.crossed(&vd);
    let x1 = loc1.x() + par1 * inter1.x();
    let y1 = loc1.y() + par1 * inter1.y();
    let z1 = loc1.z() + par1 * inter1.z();
    let x2 = loc2.x() + par2 * inter2.x();
    let y2 = loc2.y() + par2 * inter2.y();
    let z2 = loc2.z() + par2 * inter2.z();
    let mut pt1 = GpPnt::new((x1 + x2) * 0.5, (y1 + y2) * 0.5, (z1 + z2) * 0.5);
    let dir1 = GpDir::from_xyz(&vd.xyz()).ok()?;
    if a_mvd < A_TRESH_ANG {
        // `cxx:449-508`: with a near-zero angle between the planes the origin
        // above carries ~1e-4 of error, so OCCT re-derives it from two
        // line/plane solves (`IntAna_IntConicQuad`).
        let a_dist1 = a1 * pt1.x() + b1 * pt1.y() + c1 * pt1.z() + d1;
        let a_dist2 = a2 * pt1.x() + b2 * pt1.y() + c2 * pt1.z() + d2;
        if a_dist1.abs() > A_TRESH_DIST || a_dist2.abs() > A_TRESH_DIST {
            let a_dn1 = GpDir::from_xyz(&n1.xyz()).ok()?;
            let a_l1 = GpLin::from_pnt_dir(pt1, a_dn1);
            let icq = IntAnaIntConicQuad::line_plane(&a_l1, p1, tol_ang, tol, 0.0);
            if !icq.is_done() {
                return None;
            }
            let a_pnt1 = *icq.point(1);
            let a_dl2 = dir1.crossed(&a_dn1).ok()?;
            let a_l2 = GpLin::from_pnt_dir(a_pnt1, a_dl2);
            let icq = IntAnaIntConicQuad::line_plane(&a_l2, p2, tol_ang, tol, 0.0);
            if !icq.is_done() || icq.is_parallel() {
                return None;
            }
            pt1 = *icq.point(1);
        }
    }
    Some(GpLin::from_pnt_dir(pt1, dir1))
}

/// `aEps = 1.e-16` (`IntAna_QuadQuadGeo.cxx:423`).
const A_EPS: f64 = 1.0e-16;
/// `aTreshAng = 2.e-6` (1.e-4 deg) (`IntAna_QuadQuadGeo.cxx:460`).
const A_TRESH_ANG: f64 = 2.0e-6;
/// `aTreshDist = 1.e-12` (`IntAna_QuadQuadGeo.cxx:461`).
const A_TRESH_DIST: f64 = 1.0e-12;

/// Intersection of three planes. Port of `IntAna_Int3Pln::Perform`
/// (`IntAna_Int3Pln.cxx:40-64`): a single point, or `NoResult` when the
/// `math_Gauss` decomposition is singular (`empt = true`).
pub fn three_planes_intersect(p1: &GpPln, p2: &GpPln, p3: &GpPln) -> Result<Intersect3Pln, String> {
    let (a1, b1, c1, d1) = plane_coeffs(p1);
    let (a2, b2, c2, d2) = plane_coeffs(p2);
    let (a3, b3, c3, d3) = plane_coeffs(p3);
    let m = [
        [a1, b1, c1],
        [a2, b2, c2],
        [a3, b3, c3],
    ];
    let d = [d1, d2, d3];
    match solve3(m, d) {
        Some(x) => Ok(Intersect3Pln::Point(GpPnt::new(x[0], x[1], x[2]))),
        None => Ok(Intersect3Pln::NoResult),
    }
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
    /// `done = false`: the overload declined the result it had built and the
    /// caller must fall back to the bi-parametric method instead of reading a
    /// curve. `IntAna_QuadQuadGeo.cxx:929-950` (plane/cone) is the only
    /// producer so far, and every consumer tests `IsDone()` before using the
    /// result (`IntPatch_ImpImpIntersection.cxx:3480-3483`).
    NotDone,
}

impl PartialEq for QuadricIntersection {
    fn eq(&self, other: &Self) -> bool {
        use QuadricIntersection::*;
        match (self, other) {
            (None, None) | (Same, Same) => true,
            (NotDone, NotDone) => true,
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

/// Plane ∩ plane. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Pln)`
/// (`IntAna_QuadQuadGeo.cxx:389-511`): collinear normals give `IntAna_Same`
/// or `IntAna_Empty` (`cxx:414-418`), otherwise the single `IntAna_Line`
/// built by [`plane_plane_line`].
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
        match plane_plane_line(p1, p2, tol_ang, tol) {
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
    if (dist.abs() - radius).abs() < epsilon(radius) {
        // `cxx:998-1003`: tangency, the sphere center projected onto the plane.
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
pub(super) enum LinePlaneHit {
    /// The line crosses the plane: parameter on the line and the point.
    Cross { param: f64, point: GpPnt },
    /// The line is parallel to the plane.
    Parallel,
}

pub(super) fn line_plane(l: &GpLin, p: &GpPln, tol_ang: f64) -> LinePlaneHit {
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

/// Plane ∩ cylinder. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cylinder,
/// Tolang, Tol, H)` (`IntAna_QuadQuadGeo.cxx:543-723`).
///
/// `h` is the `H` argument (`cxx:546`, default `0` in `hxx:94`): the height of
/// the cylinder face along its axis. `IntPCy`
/// (`IntPatch_ImpImpIntersection.cxx:2579-2588`) passes `VMax - VMin`, and `0`
/// when either V bound is infinite. It feeds only the near-parallel
/// refinement (`cxx:568-598`), which is what turns a strongly eccentric
/// ellipse into a line pair when the axis deviates from the plane by less than
/// `2*Tol` over the whole face height.
pub fn quadric_quadric_plane_cylinder(
    p: &GpPln,
    cyl: &GpCylinder,
    tol_ang: f64,
    tol: f64,
    h: f64,
) -> QuadricIntersection {
    let (a, b, c, d) = plane_coeffs(p);
    let radius = cyl.radius();
    let axec = cyl.axis();
    let axis_origin = *axec.location();
    let axis_dir = *axec.direction();
    let (x, y, z) = (axis_origin.x(), axis_origin.y(), axis_origin.z());
    // Signed distance from the axis origin to the plane (`cxx:567`).
    let dist = a * x + b * y + c * z + d;
    let normp = p.pos.direction();

    // `cxx:568-598`: with the axis nearly perpendicular to the plane normal,
    // decide from the deviation over the face height whether the section has
    // to be reported as lines; `tolang`/`toltang` are relaxed when it is.
    let mut tolang = tol_ang;
    let mut toltang = tol;
    let mut newparams = false;
    let ldv = GpVec::from_xyz(axis_dir.xyz());
    let npv = GpVec::from_xyz(normp.xyz());
    let d_a = ldv.angle(&npv).abs();
    if d_a > PI / 4.0 {
        let dangle = (d_a - PI / 2.0).abs();
        if dangle > tol_ang {
            let sinda = dangle.sin().abs();
            let dif = (sinda - tol).abs();
            if dif < tol || (h > 0.0 && sinda * h < 2.0 * tol) {
                tolang = sinda * 2.0;
                // Keep the endpoint-distance check of IntAna_IntConicQuad from
                // rejecting the parallel status for a small deviation over H.
                toltang = tol.max(sinda * h * 1.01);
                newparams = true;
            }
        }
    }

    // `cxx:559` `gp_Lin axec(Cl.Axis())`. `cxx:601`: the parallel test is
    // `IntAna_IntConicQuad::IsParallel()` with those tolerances and the face
    // height as the bounding length.
    let axec_lin = GpLin::from_pnt_dir(axis_origin, axis_dir);
    let inter = IntAnaIntConicQuad::line_plane(&axec_lin, p, tolang, toltang, h);

    if inter.is_parallel() {
        // `cxx:604-684`: line type, one or two generatrices.
        let omega = GpPnt::new(x - dist * a, y - dist * b, z - dist * c);
        if (dist.abs() - radius).abs() < tol {
            // `cxx:610-634`: tangent plane → a single line.
            let dir1 = if newparams {
                // `cxx:617-631`: rebuild the direction from the point 100 units
                // along the axis, re-projected onto the plane.
                let t = GpPnt::new(
                    x + 100.0 * axis_dir.x(),
                    y + 100.0 * axis_dir.y(),
                    z + 100.0 * axis_dir.z(),
                );
                let distt = a * t.x() + b * t.y() + c * t.z() + d;
                let p1 = GpPnt::new(
                    t.x() - distt * a,
                    t.y() - distt * b,
                    t.z() - distt * c,
                );
                // `gp_Dir dd1(vv1)` (`cxx:629`) raises on a null vector.
                let Ok(dd1) = GpDir::from_vec(&GpVec::from_pnts(&omega, &p1)) else {
                    return QuadricIntersection::None;
                };
                dd1
            } else {
                axis_dir
            };
            QuadricIntersection::Line(GpLin::from_pnt_dir(omega, dir1))
        } else if dist.abs() < radius {
            // `cxx:636-675`: two lines.
            let hh = (radius * radius - dist * dist).sqrt();
            let axey = axis_dir.xyz().crossed(&normp.xyz());
            let p1 = omega.translated_vec(&GpVec::from_xyz(&axey.multiplied(-hh)));
            let p2 = omega.translated_vec(&GpVec::from_xyz(&axey.multiplied(hh)));
            if newparams {
                // `cxx:647-668`: same re-projection, with the offset recomputed
                // at the far point (`ht`, clamped to 0 when the argument is
                // negative — `cxx:654-655`).
                let t = GpPnt::new(
                    x + 100.0 * axis_dir.x(),
                    y + 100.0 * axis_dir.y(),
                    z + 100.0 * axis_dir.z(),
                );
                let distt = a * t.x() + b * t.y() + c * t.z() + d;
                let arg = radius * radius - distt * distt;
                let ht = if arg > 0.0 { arg.sqrt() } else { 0.0 };
                let omega1 = GpPnt::new(
                    t.x() - distt * a,
                    t.y() - distt * b,
                    t.z() - distt * c,
                );
                let ppt1 = omega1.translated_vec(&GpVec::from_xyz(&axey.multiplied(-ht)));
                let ppt2 = omega1.translated_vec(&GpVec::from_xyz(&axey.multiplied(ht)));
                let Ok(dd1) = GpDir::from_vec(&GpVec::from_pnts(&p1, &ppt1)) else {
                    return QuadricIntersection::None;
                };
                let Ok(dd2) = GpDir::from_vec(&GpVec::from_pnts(&p2, &ppt2)) else {
                    return QuadricIntersection::None;
                };
                QuadricIntersection::TwoLines(
                    GpLin::from_pnt_dir(p1, dd1),
                    GpLin::from_pnt_dir(p2, dd2),
                )
            } else {
                QuadricIntersection::TwoLines(
                    GpLin::from_pnt_dir(p1, axis_dir),
                    GpLin::from_pnt_dir(p2, axis_dir),
                )
            }
        } else {
            // `cxx:680-683`: IntAna_Empty.
            QuadricIntersection::None
        }
    } else {
        // `cxx:685-718`: the axis crosses the plane in one point; the section
        // is a circle (axis ⊥ plane) or an ellipse.
        let axey = normp.xyz().crossed(&axis_dir.xyz());
        let sint = axey.modulus();
        let center = *inter.point(1);
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

/// `EllipseLimit = 1.0E+9` (function-static in `IntAna_QuadQuadGeo.cxx:927`).
const ELLIPSE_LIMIT: f64 = 1.0e9;
/// `HyperbolaLimit = 2.0E+6` (function-static in `IntAna_QuadQuadGeo.cxx:928`).
const HYPERBOLA_LIMIT: f64 = 2.0e6;

/// Plane ∩ cone. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cone,
/// Tolang, Tol)` (`IntAna_QuadQuadGeo.cxx:752-972`).
///
/// The tail `cxx:922-951` rejects results the cxx calls "too extreme" — an
/// `IntAna_Ellipse` whose `param1` or `param1bis` exceeds
/// `EllipseLimit = 1.0E+9` (`cxx:926-934`, any `nbint`) and an
/// `IntAna_Hyperbola` whose parameters exceed `HyperbolaLimit = 2.0E+6`
/// (`cxx:935-942` on `param2`/`param2bis` for `nbint >= 2`, `cxx:943-951` on
/// `param1`/`param1bis` for `nbint >= 1`) — by setting `done = false` and
/// returning, so the caller treats that pair with the bi-parametric method
/// instead. In this port the two parameter pairs are exactly the
/// `major_radius`/`minor_radius` carried by the returned `GpElips`/`GpHypr`
/// (`param1 = MajorRadius()`, `param1bis = MinorRadius()`, and for both
/// hyperbola arms `param2 = param1`, `param2bis = param1bis`, `nbint = 2`),
/// so the check runs on the built curve and yields
/// [`QuadricIntersection::NotDone`].
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

    let res = if dist.abs() < tol {
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
    };

    // `cxx:922-951`: the tail rejects ellipses and hyperbolas it calls "too
    // extreme" by setting `done = false` and returning, which tells the caller
    // to treat the pair with the bi-parametric method instead. `param1` /
    // `param1bis` are the major / minor radii stored above, and both hyperbola
    // arms also carry them in `param2` / `param2bis` with `nbint = 2`, so the
    // two cxx blocks (`cxx:935-942`, `cxx:943-951`) reduce to one test here.
    match &res {
        QuadricIntersection::Ellipse(e) => {
            if e.major_radius().abs() > ELLIPSE_LIMIT || e.minor_radius().abs() > ELLIPSE_LIMIT {
                return QuadricIntersection::NotDone;
            }
        }
        QuadricIntersection::Hyperbola(h) => {
            if h.major_radius.abs() > HYPERBOLA_LIMIT
                || h.minor_radius.abs() > HYPERBOLA_LIMIT
            {
                return QuadricIntersection::NotDone;
            }
        }
        _ => {}
    }
    res
}
// ---------------------------------------------------------------------------
// IntAna_QuadQuadGeo: cylinder×cylinder, cylinder×sphere, sphere×cone,
// cone×cone, cylinder×cone.
// ---------------------------------------------------------------------------

/// Snap an axis direction that is nearly axis-aligned to the exact axis.
/// Port of `RefineDir` (IntAna_QuadQuadGeo.cxx): it keeps the cross products
/// used by the axis relation code away from exact degeneracies.
pub(super) fn refine_dir(d: &mut GpDir) {
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
pub(super) fn det33(a: &GpVec, b: &GpVec, c: &GpVec) -> f64 {
    a.dot(&b.crossed(c))
}
