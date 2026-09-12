//! Surface–surface intersection curves.
//! Source: `IntPatch_Intersection`, `IntAna_IntConicQuad`, `IntWalk_PWalking` (TKMath / TKGeomAlgo).
//!
//! Ports OCCT's surface-intersection machinery to exact analytic closed forms
//! where the pair is amenable (plane∩sphere, sphere∩sphere, plane∩cylinder)
//! and a marching-squares grid tracer for the general case. The result of
//! intersecting two parametric surfaces is one or more 3D curves together with
//! the `(u, v)` parameters of each sampled curve point on both surfaces — the
//! data `BRepAlgoAPI_Section` needs to build p-curves on the trimmed faces.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
use occt_geom::{
    Curve, GeomCylinder, GeomLine, GeomPlane, GeomSphere, GeomTorus, Surface,
};

use crate::brep_surface::{classify_surface, SurfaceKind};

/// A single 3D intersection curve between two surfaces, sampled into points
/// and per-surface parameters.
#[derive(Clone)]
pub struct IntersectionCurve {
    /// The 3D curve of intersection.
    pub curve: Arc<dyn Curve>,
    /// Sample points on the curve (curve evaluated on a uniform parameter grid).
    pub points: Vec<GpPnt>,
    /// The `(u, v)` parameters of each sample point on surface A.
    pub on_a: Vec<(f64, f64)>,
    /// The `(u, v)` parameters of each sample point on surface B.
    pub on_b: Vec<(f64, f64)>,
}

impl std::fmt::Debug for IntersectionCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IntersectionCurve")
            .field("points", &self.points)
            .field("on_a", &self.on_a)
            .field("on_b", &self.on_b)
            .finish()
    }
}

/// Result of intersecting two surfaces.
#[derive(Clone)]
pub enum SurfaceIntersection {
    /// One or more distinct intersection curves.
    Curves(Vec<IntersectionCurve>),
    /// The surfaces are geometrically coincident over a region.
    Coincident,
    /// No intersection found.
    None,
}

impl std::fmt::Debug for SurfaceIntersection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SurfaceIntersection::Curves(c) => f.debug_tuple("Curves").field(c).finish(),
            SurfaceIntersection::Coincident => f.write_str("Coincident"),
            SurfaceIntersection::None => f.write_str("None"),
        }
    }
}

#[path = "intpatch_geom.rs"]
mod geom;
#[path = "intpatch_bspline.rs"]
mod bspline;
#[path = "intpatch_analytic.rs"]
mod analytic;
#[path = "intpatch_trace.rs"]
mod trace;

pub use geom::*;
pub use bspline::*;
pub use analytic::{
    intersect_plane_cone, intersect_plane_cylinder, intersect_plane_sphere, intersect_plane_torus,
    intersect_sphere_sphere,
};
pub use trace::{
    chain_intersection_points, intersect_general_surfaces, intersection_curve_points,
    polyline_to_curve, surfaces_intersect_general, trace_surface_curve, PolylineCurve,
};

use analytic::sample_curve_on;

#[path = "intpatch_impimp.rs"]
pub mod impimp;
pub use impimp::{ImpImpIntersection, IntStatus as ImpImpStatus};

#[path = "intpatch_impprm.rs"]
mod impprm;
#[path = "intpatch_prmprm.rs"]
mod prmprm;
#[path = "intpatch_special_points.rs"]
pub mod special_points;
#[path = "intpatch_aline_to_wline.rs"]
pub mod aline_to_wline;
#[path = "intpatch_wline_tool.rs"]
pub mod wline_tool;
#[path = "intpatch_intersection.rs"]
pub mod intersection;
pub use intersection::PatchIntersection;

/// Sample points from the intersection curves of `a` and `b`, retaining only
/// those within `tol` of both surfaces (verification helper).
pub fn points_on_both(a: &dyn Surface, b: &dyn Surface, tol: f64, nsamples: usize) -> Vec<GpPnt> {
    let mut out = Vec::new();
    match surface_surface_intersection(a, b, tol) {
        SurfaceIntersection::Curves(curves) => {
            for ic in &curves {
                let (f0, f1) = (ic.curve.first_parameter(), ic.curve.last_parameter());
                let n = nsamples.max(4);
                for i in 0..n {
                    let t = if (f1 - f0).is_finite() {
                        f0 + (f1 - f0) * i as f64 / (n - 1) as f64
                    } else {
                        -2.0 + 4.0 * i as f64 / (n - 1) as f64
                    };
                    let p = ic.curve.d0(t);
                    if distance_to_surface(&p, a) <= tol && distance_to_surface(&p, b) <= tol {
                        out.push(p);
                    }
                }
            }
        }
        _ => {}
    }
    out
}

// ---------------------------------------------------------------------------
// Dispatcher
// ---------------------------------------------------------------------------

/// Intersect two surfaces, dispatching on their analytic types and falling back
/// to the grid tracer for general pairs.
pub fn surface_surface_intersection(a: &dyn Surface, b: &dyn Surface, tol: f64) -> SurfaceIntersection {
    let kind_a = classify_surface(a);
    let kind_b = classify_surface(b);

    match (kind_a, kind_b) {
        (SurfaceKind::Plane, SurfaceKind::Plane) => {
            let pa = match plane_from_surface(a) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let pb = match plane_from_surface(b) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let n1 = GpVec::from_xyz(pa.axis().direction().xyz());
            let n2 = GpVec::from_xyz(pb.axis().direction().xyz());
            if n1.cross_magnitude(&n2) <= tol {
                let dd = n1
                    .normalized()
                    .dot(&GpVec::from_pnts(&pa.location(), &pb.location()))
                    .abs();
                return if dd <= tol { SurfaceIntersection::Coincident } else { SurfaceIntersection::None };
            }
            return match crate::face_face::plane_plane_intersection(&pa, &pb) {
                Some((origin, dir)) => {
                    let lin = occt_core::gp::GpLin::from_pnt_dir(origin, GpDir::from_vec(&dir).unwrap_or_default());
                    let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(lin));
                    SurfaceIntersection::Curves(vec![sample_curve_on(curve, a, b, 32)])
                }
                None => SurfaceIntersection::None,
            };
        }
        (SurfaceKind::Plane, SurfaceKind::Sphere) => {
            let pa = match plane_from_surface(a) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let (c, r) = match sphere_params(b) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            return match intersect_plane_sphere(&pa, c, r) {
                Some(ic) => SurfaceIntersection::Curves(vec![ic]),
                None => SurfaceIntersection::None,
            };
        }
        (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let pb = match plane_from_surface(b) {
                Some(p) => p,
                None => return SurfaceIntersection::None,
            };
            let (c, r) = match sphere_params(a) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            return match intersect_plane_sphere(&pb, c, r) {
                Some(ic) => SurfaceIntersection::Curves(vec![ic]),
                None => SurfaceIntersection::None,
            };
        }
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let (c1, r1) = match sphere_params(a) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            let (c2, r2) = match sphere_params(b) {
                Some(x) => x,
                None => return SurfaceIntersection::None,
            };
            return match intersect_sphere_sphere(c1, r1, c2, r2) {
                Some(ic) => SurfaceIntersection::Curves(vec![ic]),
                None => SurfaceIntersection::None,
            };
        }
        _ => {}
    }

    // General fallback: grid tracer chained into per-curve polylines.
    intersect_general_surfaces(a, b, tol)
}
// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::GpAx1;

    const TOL: f64 = 1e-6;

    fn unit_sphere(center: GpPnt) -> Arc<dyn Surface> {
        let ax3 = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        Arc::new(GeomSphere::new(GpSphere::new(ax3, 1.0).unwrap()))
    }

    fn plane_z(z: f64) -> GpPln {
        GpPln::new(GpAx3::new(GpPnt::new(0.0, 0.0, z), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap())
    }

    #[test]
    fn plane_sphere_circle_radius_and_on_surface() {
        // Plane z = 0.5 cutting the unit sphere at origin: circle radius sqrt(1-0.25) ≈ 0.866.
        let pln = plane_z(0.5);
        let ic = intersect_plane_sphere(&pln, GpPnt::zero(), 1.0).expect("circle");
        assert_eq!(ic.points.len(), 48);
        for p in &ic.points {
            let r = p.distance(&GpPnt::zero());
            assert!((r - 1.0).abs() < 1e-9, "on sphere: r {r}");
            assert!((p.z() - 0.5).abs() < 1e-9, "in plane: z {}", p.z());
        }
        // Radius of the circle ≈ sqrt(3/4).
        let rad = ic.points[0].distance(&GpPnt::new(0.0, 0.0, 0.5));
        assert!((rad - (0.75f64.sqrt())).abs() < 1e-9, "circle radius {rad}");
        // Every sampled point is within tol of both surfaces.
        let s = unit_sphere(GpPnt::zero());
        let a: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        for p in &ic.points {
            assert!(distance_to_surface(p, a.as_ref()) < TOL);
            assert!(distance_to_surface(p, s.as_ref()) < TOL);
        }
    }

    #[test]
    fn plane_sphere_no_intersection() {
        let pln = plane_z(2.0); // above the unit sphere
        assert!(intersect_plane_sphere(&pln, GpPnt::zero(), 1.0).is_none());
    }

    #[test]
    fn plane_sphere_tangent() {
        let pln = plane_z(1.0);
        let ic = intersect_plane_sphere(&pln, GpPnt::zero(), 1.0).expect("tangent circle");
        // Tangent: radius 0, all points collapse to the tangent point.
        let p = ic.points[0];
        assert!(p.distance(&GpPnt::new(0.0, 0.0, 1.0)) < 1e-9, "tangent point {p:?}");
    }

    #[test]
    fn sphere_sphere_circle_on_both() {
        // Two unit spheres, centers 1.5 apart.
        let ic = intersect_sphere_sphere(GpPnt::zero(), 1.0, GpPnt::new(1.5, 0.0, 0.0), 1.0).expect("circle");
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(1.5, 0.0, 0.0));
        // Intersection plane at x = 0.75, radius sqrt(1 - 0.75²).
        let expected_r = (1.0 - 0.75f64 * 0.75).sqrt();
        for p in &ic.points {
            assert!((p.x() - 0.75).abs() < 1e-9, "plane x {}", p.x());
            let r = GpPnt::new(0.0, p.y(), p.z()).distance(&GpPnt::zero());
            assert!((r - expected_r).abs() < 1e-9, "radius {r} (expected {expected_r})");
            assert!(distance_to_surface(p, s1.as_ref()) < TOL);
            assert!(distance_to_surface(p, s2.as_ref()) < TOL);
        }
    }

    #[test]
    fn sphere_sphere_tangent_or_none() {
        // Tangent spheres (d = r1 + r2) → degenerate single point (radius 0).
        let ic = intersect_sphere_sphere(GpPnt::zero(), 1.0, GpPnt::new(2.0, 0.0, 0.0), 1.0).expect("tangent circle");
        let p = ic.points[0];
        assert!(p.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-6, "tangent point {p:?}");

        // Disjoint spheres → None.
        assert!(intersect_sphere_sphere(GpPnt::zero(), 1.0, GpPnt::new(5.0, 0.0, 0.0), 1.0).is_none());
    }

    #[test]
    fn plane_cylinder_circle() {
        // Plane z = 1 perpendicular to the Z-axis cylinder radius 2 → circle radius 2 at z=1.
        let pln = plane_z(1.0);
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let ic = intersect_plane_cylinder(&pln, &ax, 2.0).expect("circle");
        for p in &ic.points {
            assert!((p.z() - 1.0).abs() < 1e-9);
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - 2.0).abs() < 1e-9, "radius {r}");
        }
    }

    #[test]
    fn plane_cylinder_parallel_lines() {
        // Plane x = 0.5 parallel to the Z-axis cylinder radius 1 → two lines at y = ±sqrt(0.75).
        let pln = GpPln::new(GpAx3::new(
            GpPnt::new(0.5, 0.0, 0.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let ic = intersect_plane_cylinder(&pln, &ax, 1.0).expect("lines");
        assert!(ic.points.len() >= 2);
        let mut ys: Vec<f64> = ic.points.iter().map(|p| p.y().abs()).collect();
        ys.sort_by(f64::total_cmp);
        let target = (0.75f64).sqrt();
        assert!((ys[0] - target).abs() < 1e-6, "line y {target}, got {}", ys[0]);
    }

    #[test]
    fn sphere_cylinder_via_trace() {
        // Unit sphere at origin ∩ cylinder radius 0.5 axis Z → circle at |z| = sqrt(0.75).
        let s = unit_sphere(GpPnt::zero());
        let ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let c: Arc<dyn Surface> = Arc::new(GeomCylinder::new(GpCylinder::new(ax3, 0.5).unwrap()));
        let pts = trace_surface_curve(s.as_ref(), c.as_ref(), 0.03);
        assert!(!pts.is_empty(), "traced points");
        for p in pts.iter().take(20) {
            assert!(distance_to_surface(p, s.as_ref()) < 0.05);
            assert!(distance_to_surface(p, c.as_ref()) < 0.05);
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - 0.5).abs() < 0.05, "on cylinder r {r}");
        }
    }

    #[test]
    fn refine_point_on_sphere_converges() {
        let s = unit_sphere(GpPnt::zero());
        let target = GpPnt::new(0.0, 0.0, 1.0);
        let (u, v, p) = refine_point_on_surface(s.as_ref(), target, 0.5, 0.5, 20);
        assert!(p.distance(&target) < 1e-6, "refined {p:?} vs {target:?}");
        assert!((v - std::f64::consts::FRAC_PI_2).abs() < 1e-4, "v {v}");
        let _ = u;
    }

    #[test]
    fn points_on_both_filters() {
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(1.5, 0.0, 0.0));
        let pts = points_on_both(s1.as_ref(), s2.as_ref(), 1e-6, 32);
        assert!(!pts.is_empty());
        for p in &pts {
            assert!(distance_to_surface(p, s1.as_ref()) < 1e-6);
            assert!(distance_to_surface(p, s2.as_ref()) < 1e-6);
        }
    }

    #[test]
    fn surface_surface_intersection_sphere_sphere() {
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(1.5, 0.0, 0.0));
        match surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-6) {
            SurfaceIntersection::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                assert!(!curves[0].points.is_empty());
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn surface_surface_intersection_disjoint_none() {
        let s1 = unit_sphere(GpPnt::zero());
        let s2 = unit_sphere(GpPnt::new(5.0, 0.0, 0.0));
        assert!(matches!(
            surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-6),
            SurfaceIntersection::None
        ));
    }

    #[test]
    fn plane_cone_circle() {
        // Cone apex at origin, semi-angle 30°, plane z = 2 → circle radius 2·tan30 ≈ 1.155.
        let pln = plane_z(2.0);
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let semi = 30f64.to_radians();
        let ic = intersect_plane_cone(&pln, GpPnt::zero(), &ax, semi).expect("circle");
        let expected_r = 2.0 * semi.tan();
        for p in &ic.points {
            assert!((p.z() - 2.0).abs() < 1e-9);
            let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((r - expected_r).abs() < 1e-9, "radius {r} (expected {expected_r})");
        }
    }

    #[test]
    fn plane_torus_two_circles() {
        // Torus major 3 minor 1; plane y=0 (contains axis) → two circles radius 1 at x=±3.
        let pln = GpPln::new(GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ).unwrap());
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let ic = intersect_plane_torus(&pln, GpPnt::zero(), &ax, 3.0, 1.0).expect("two circles");
        assert!(!ic.points.is_empty());
        let torus_ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let t: Arc<dyn Surface> = Arc::new(GeomTorus::new(GpTorus::new(torus_ax3, 3.0, 1.0).unwrap()));
        for p in &ic.points {
            assert!((p.y()).abs() < 1e-9, "in plane y {}", p.y());
            assert!(distance_to_surface(p, t.as_ref()) < 1e-6);
            // Each point lies on one of the two circles (radius 1) centered at x = ±3.
            let d1 = p.distance(&GpPnt::new(3.0, 0.0, 0.0));
            let d2 = p.distance(&GpPnt::new(-3.0, 0.0, 0.0));
            assert!(
                (d1 - 1.0).abs() < 1e-9 || (d2 - 1.0).abs() < 1e-9,
                "on a minor circle: d1 {d1}, d2 {d2}"
            );
        }
    }

    #[test]
    fn trace_and_polyline_curve() {
        let c = PolylineCurve {
            pts: vec![GpPnt::zero(), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(2.0, 0.0, 0.0)],
        };
        assert!(c.d0(0.0).distance(&GpPnt::zero()) < 1e-12);
        assert!(c.d0(0.5).distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-12);
        assert!(c.d0(1.0).distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-12);
    }

    #[test]
    fn plane_plane_intersection_line() {
        let xy = plane_z(0.0);
        let yz = GpPln::new(GpAx3::new(
            GpPnt::zero(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let s1: Arc<dyn Surface> = Arc::new(GeomPlane::new(xy.clone()));
        let s2: Arc<dyn Surface> = Arc::new(GeomPlane::new(yz.clone()));
        match surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-9) {
            SurfaceIntersection::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                for p in &curves[0].points {
                    assert!(p.z().abs() < 1e-9);
                    assert!(p.x().abs() < 1e-9);
                }
            }
            other => panic!("expected line, got {other:?}"),
        }
    }

    // -- general (non-analytic) surface intersection extensions --

    fn sphere_r(center: GpPnt, r: f64) -> Arc<dyn Surface> {
        let ax3 = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        Arc::new(GeomSphere::new(GpSphere::new(ax3, r).unwrap()))
    }

    fn cylinder_z(r: f64) -> Arc<dyn Surface> {
        let ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        Arc::new(GeomCylinder::new(GpCylinder::new(ax3, r).unwrap()))
    }

    /// A degree-1 B-spline patch over a tilted plane `z = u + 0.5·v` on a
    /// 4×4 grid spanning `[0, 1]²`.
    fn tilted_plane_patch() -> Arc<dyn Surface> {
        let mut grid: Vec<Vec<GpPnt>> = Vec::new();
        for i in 0..=3usize {
            let mut row = Vec::new();
            for j in 0..=3usize {
                let u = i as f64 / 3.0;
                let v = j as f64 / 3.0;
                row.push(GpPnt::new(u, v, u + 0.5 * v));
            }
            grid.push(row);
        }
        make_bspline_surface_from_grid(&grid, 1, 1).expect("patch fit")
    }

    #[test]
    fn general_plane_plane_line() {
        let s1: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(0.0)));
        let yz = GpPln::new(GpAx3::new(
            GpPnt::zero(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let s2: Arc<dyn Surface> = Arc::new(GeomPlane::new(yz));
        // The analytic dispatcher still returns a line for two planes.
        match surface_surface_intersection(s1.as_ref(), s2.as_ref(), 1e-6) {
            SurfaceIntersection::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                for p in &curves[0].points {
                    assert!(p.z().abs() < 1e-6);
                    assert!(p.x().abs() < 1e-6);
                }
            }
            other => panic!("expected line, got {other:?}"),
        }
    }

    #[test]
    fn general_sphere_sphere_via_trace() {
        // Two radius-2 spheres with centers 3 apart intersect in a circle;
        // the general tracer (bypassing the analytic sphere-sphere path) finds it.
        let a = sphere_r(GpPnt::zero(), 2.0);
        let b = sphere_r(GpPnt::new(3.0, 0.0, 0.0), 2.0);
        let tol = 0.05;
        match intersect_general_surfaces(a.as_ref(), b.as_ref(), tol) {
            SurfaceIntersection::Curves(curves) => {
                assert!(!curves.is_empty(), "expected at least one traced curve");
                for ic in &curves {
                    assert!(ic.points.len() >= 2, "chained points");
                    for p in &ic.points {
                        assert!(distance_to_surface(p, a.as_ref()) < 4.0 * tol, "on sphere a: {p:?}");
                        assert!(distance_to_surface(p, b.as_ref()) < 4.0 * tol, "on sphere b: {p:?}");
                    }
                }
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn general_bspline_plane_intersection() {
        let patch = tilted_plane_patch();
        let pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(0.5)));
        let tol = 0.02;
        match intersect_general_surfaces(patch.as_ref(), pln.as_ref(), tol) {
            SurfaceIntersection::Curves(curves) => {
                assert!(!curves.is_empty(), "expected a traced line");
                let mut found = 0;
                for ic in &curves {
                    for p in &ic.points {
                        assert!(distance_to_surface(p, patch.as_ref()) < 6.0 * tol);
                        assert!(distance_to_surface(p, pln.as_ref()) < 6.0 * tol);
                        found += 1;
                    }
                }
                assert!(found >= 2, "at least two points on the intersection line");
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn general_cylinder_plane_parallel_lines() {
        let cyl = cylinder_z(1.0);
        let pln_x = GpPln::new(GpAx3::new(
            GpPnt::new(0.5, 0.0, 0.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        ).unwrap());
        let pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln_x));
        let tol = 0.05;
        match intersect_general_surfaces(cyl.as_ref(), pln.as_ref(), tol) {
            SurfaceIntersection::Curves(curves) => {
                assert!(!curves.is_empty(), "expected generatrix lines");
                for ic in &curves {
                    for p in &ic.points {
                        assert!(distance_to_surface(p, cyl.as_ref()) < 4.0 * tol, "on cylinder: {p:?}");
                        assert!(distance_to_surface(p, pln.as_ref()) < 4.0 * tol, "on plane: {p:?}");
                    }
                }
            }
            other => panic!("expected curves, got {other:?}"),
        }
    }

    #[test]
    fn chain_intersection_merges() {
        // Two polylines whose endpoint regions overlap within tol merge to one chain.
        let chain_a: Vec<GpPnt> = (0..=10).map(|i| GpPnt::new(i as f64 * 0.1, 0.0, 0.0)).collect();
        let chain_b: Vec<GpPnt> = (0..=10).map(|i| GpPnt::new(0.95 + i as f64 * 0.1, 0.0, 0.0)).collect();
        let mut pts = chain_a.clone();
        pts.extend(chain_b);
        let chains = chain_intersection_points(&pts, 0.12);
        assert_eq!(chains.len(), 1, "overlapping polylines merge into one chain");
        assert_eq!(chains[0].len(), 22, "all 22 points chained in order");
        // Well-separated curves stay separate.
        let far: Vec<GpPnt> = (0..=5).map(|i| GpPnt::new(5.0 + i as f64 * 0.1, 0.0, 0.0)).collect();
        let mut pts2 = pts;
        pts2.extend(far);
        let chains2 = chain_intersection_points(&pts2, 0.12);
        assert_eq!(chains2.len(), 2, "a separated curve starts a new chain");
    }

    #[test]
    fn polyline_to_curve_on_params() {
        let s = unit_sphere(GpPnt::zero());
        let pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(0.5)));
        // Circle of intersection: radius sqrt(0.75) at z = 0.5.
        let r = (0.75f64).sqrt();
        let poly: Vec<GpPnt> = (0..16)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 16.0;
                GpPnt::new(r * a.cos(), r * a.sin(), 0.5)
            })
            .collect();
        let ic = polyline_to_curve(&poly, s.as_ref(), pln.as_ref(), 1e-4).expect("curve");
        assert_eq!(ic.points.len(), poly.len());
        assert_eq!(ic.on_a.len(), poly.len());
        assert_eq!(ic.on_b.len(), poly.len());
        for (p, (ua, va)) in ic.points.iter().zip(&ic.on_a) {
            let q = s.d0(*ua, *va);
            assert!(q.distance(p) < 1e-5, "on_a reconstructs {p:?} as {q:?}");
        }
        for (p, (ub, vb)) in ic.points.iter().zip(&ic.on_b) {
            let q = pln.d0(*ub, *vb);
            assert!(q.distance(p) < 1e-5, "on_b reconstructs {p:?} as {q:?}");
        }
    }

    #[test]
    fn general_nonintersecting_empty() {
        let s = unit_sphere(GpPnt::zero());
        let far_pln: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_z(10.0)));
        assert!(matches!(
            intersect_general_surfaces(s.as_ref(), far_pln.as_ref(), 0.05),
            SurfaceIntersection::None
        ));
        assert!(!surfaces_intersect_general(s.as_ref(), far_pln.as_ref(), 0.05));
        assert!(surfaces_intersect_general(s.as_ref(), unit_sphere(GpPnt::new(1.0, 0.0, 0.0)).as_ref(), 0.05));
    }

    #[test]
    fn fit_bspline_grid_corners() {
        let mut grid: Vec<Vec<GpPnt>> = Vec::new();
        for i in 0..=3usize {
            let mut row = Vec::new();
            for j in 0..=3usize {
                let u = i as f64 / 3.0;
                let v = j as f64 / 3.0;
                row.push(GpPnt::new(u * 2.0 - 1.0, v * 2.0 - 1.0, u * u + v));
            }
            grid.push(row);
        }
        let patch = fit_bspline_grid(&grid, 1, 1).expect("fit");
        let corners = [
            (0.0, 0.0, &grid[0][0]),
            (1.0, 0.0, &grid[3][0]),
            (0.0, 1.0, &grid[0][3]),
            (1.0, 1.0, &grid[3][3]),
        ];
        for (u, v, expected) in corners {
            let got = patch.d0(u, v);
            assert!(got.distance(expected) < 1e-9, "corner ({u},{v}): got {got:?} expected {expected:?}");
        }
    }
}
