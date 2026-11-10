//! `IntAna_QuadQuadGeo::Perform` arms for cylinder/cylinder and cylinder/sphere.

use super::super::prelude::*;
use super::*;

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
        } else if (r1p2 - dist) <= REAL_SMALL {
            // External tangency: one generatrix line.
            let pt1 = p1.translated_vec(&GpVec::from_pnts(&p1, &p2).multiplied_scalar(r1 / r1p2));
            QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
        } else if dist > rmr {
            // Two generatrix lines (or one when the base circles are tangent).
            let r1r1 = r1 * r1;
            let cos = 0.5 * (r1r1 - r2 * r2 + dist * dist) / (r1 * dist);
            let sin2 = 1.0 - cos * cos;
            let is_tangent = 4.0 * r1r1 * sin2 < tol * tol;
            let dir_a1a2 = GpVec::from_pnts(&p1, &p2).divided(dist);
            if is_tangent {
                let pt1 = p1.translated_vec(&dir_a1a2.multiplied_scalar(r1 * cos));
                QuadricIntersection::Line(GpLin::from_pnt_dir(pt1, dir_cyl))
            } else {
                let sin_val = sin2.sqrt();
                let axd = *c1.position().x_direction();
                let ayd = *c1.position().y_direction();
                let r1x = GpVec::from_xyz(axd.xyz()).multiplied_scalar(r1);
                let r1y = GpVec::from_xyz(ayd.xyz()).multiplied_scalar(r1);
                let adx = dir_a1a2.dot(&GpVec::from_xyz(axd.xyz()));
                let ady = dir_a1a2.dot(&GpVec::from_xyz(ayd.xyz()));
                let (ndx, ndy) = (adx * cos - ady * sin_val, ady * cos + adx * sin_val);
                let pt1 = p1
                    .translated_vec(&r1x.multiplied_scalar(ndx))
                    .translated_vec(&r1y.multiplied_scalar(ndy));
                let (ndx, ndy) = (adx * cos + ady * sin_val, ady * cos - adx * sin_val);
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
        // `cxx:1224` guards this arm with `RmR_Relative <=
        // myEPSILON_CYLINDER_DELTA_RADIUS` (`InitTolerances`, `cxx:356`) and
        // `A1A2.Intersect()`; in the non-parallel branch that is exactly
        // `coplanar`.
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
    // `cxx:1378-1379`: the axes must meet exactly at the sphere center, or be
    // coincident parallel axes. The exact `AxeOperator` predicate is the OCCT
    // control flow -- no extra tolerance is added here.
    if !axes_meet_at_point(&cyl.axis(), sph.position().axis(), &pt) {
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
