//! `IntAna_QuadQuadGeo::Perform(gp_Sphere, gp_Cone)`.

use super::super::prelude::*;
use super::*;

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
    // `cxx:1924-1925`: same exact-axis predicate as the cylinder/sphere arm.
    if !axes_meet_at_point(&cone.axis(), sph.position().axis(), &pt) {
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
