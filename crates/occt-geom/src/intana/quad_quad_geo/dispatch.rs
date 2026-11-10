//! `quadric_quadric` dispatcher plus the cylinder/cone arm it reaches.

use super::super::prelude::*;
use super::*;

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
            // `H` defaults to 0 (`IntAna_QuadQuadGeo.hxx:94`): only `IntPCy`
            // knows the face height (`IntPatch_ImpImpIntersection.cxx:2584-2588`),
            // and it calls the arm directly.
            quadric_quadric_plane_cylinder(p, c, tol_ang, tol, 0.0)
        }
        (Plane(p), Cone(c)) | (Cone(c), Plane(p)) => quadric_quadric_plane_cone(p, c, tol_ang, tol),
        (Cylinder(a), Cylinder(b)) => quadric_quadric_cylinder_cylinder(a, b, tol),
        (Cylinder(c), Sphere(s)) | (Sphere(s), Cylinder(c)) => quadric_quadric_cylinder_sphere(c, s, tol),
        (Sphere(s), Cone(c)) | (Cone(c), Sphere(s)) => quadric_quadric_sphere_cone(s, c, tol),
        (Cone(a), Cone(b)) => quadric_quadric_cone_cone(a, b, tol),
        (Cylinder(c), Cone(k)) | (Cone(k), Cylinder(c)) => quadric_quadric_cylinder_cone(c, k, tol),
    }
}
