//! Regression tests for `Extrema_ExtSS` / `Extrema_ExtElSS` / `Extrema_GenExtSS`.
//!
//! Port-local: the assertions pin the distances and point sets that the ported
//! `.cxx` control flow produces. OCCT itself has no unit test for these pairs.

use super::*;
use crate::bspline_surface::{bspline_surface_uniform_knots, fit_surface_grid, GeomBSplineSurface};
use crate::{GeomPlane, GeomSphere};
use occt_core::gp::{GpAx3, GpPnt, GpSphere as GpSphereT};
use occt_math::{MathFunctionSetWithDerivatives, MathVector};

fn unit_sphere_at(p: GpPnt) -> GeomSphere {
    let mut ax = GpAx3::standard();
    ax.set_location(p);
    GeomSphere::new(GpSphereT::new(ax, 1.0).unwrap())
}

fn plane_z(z: f64) -> GeomPlane {
    let mut pl = occt_core::gp::GpPln::new(GpAx3::standard());
    pl.set_location(&GpPnt::new(0.0, 0.0, z));
    GeomPlane::new(pl)
}

#[test]
fn sphere_sphere_extrema_min_and_max() {
    let s1 = unit_sphere_at(GpPnt::new(0.0, 0.0, 0.0));
    let s2 = unit_sphere_at(GpPnt::new(4.0, 0.0, 0.0));
    // Sphere x sphere is not a Plane x Plane pair, so this is the
    // `Extrema_GenExtSS` arm (`Extrema_ExtSS.cxx:180-221`). OCCT 8.0.0
    // (`specs\occt_probe\run_here.bat dummy --ssx`) reports
    // `SSX sph-sph nbExt=2 sqDist=4.0000000000003126 / 36`.
    let ext = ExtremaExtSS::new(&s1, &s2, 0.0, 0.0);
    assert!(ext.is_done(), "Extrema_ExtSS not done");
    assert!(!ext.is_parallel(), "surfaces are not planes");
    let all = surface_surface_extrema_all(&s1, &s2);
    assert_eq!(all.len(), 2, "two GenExtSS solutions");
    let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
    let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
    assert!((min.distance - 2.0).abs() < 1e-6, "min {}", min.distance);
    assert!((min.p1.x() - 1.0).abs() < 1e-6 && (min.p2.x() - 3.0).abs() < 1e-6, "closest {min:?}");
    assert!((max.distance - 6.0).abs() < 1e-6, "max {}", max.distance);
    assert!((max.p1.x() + 1.0).abs() < 1e-6 && (max.p2.x() - 5.0).abs() < 1e-6, "farthest {max:?}");
}

/// Relative comparison against an OCCT 8.0.0 reading.
fn rel_close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs().max(1e-300)
}

/// Sphere x an *unbounded* plane adaptor: the `Extrema_GenExtSS` arm, and the
/// answer is not the geometric 4.0.
///
/// `GeomAdaptor_Surface` takes `Geom_Plane::Bounds` verbatim
/// (`GeomAdaptor_Surface.hxx:150-160` -> `load`, `cxx:401-434`), i.e.
/// `±Precision::Infinite()` = `2e100`, so `Extrema_GenExtSS` samples the plane
/// over `[-2e100, 2e100]`. Its `math_BFGS` start step is
/// `0.015 / sqrt(|dir|^2)` (`math_BFGS.cxx:115-135`) -- ulp-sized against that
/// span -- so the quasi-Newton iterates stay in the `2e100` neighbourhood, and
/// the pair is then dropped by the `Extrema_ExtSS` window filter because its
/// `V1` walked to `-1.5786`, outside the sphere's `v` window `[-pi/2, pi/2]`
/// (`Extrema_ExtSS.cxx:162-166`). Only the grid maximum survives.
///
/// OCCT 8.0.0 ground truth (`specs\occt_probe\run_here.bat dummy --ssx`):
/// ```text
/// SSX sph-pln done=1 isPar=0 nbExt=1
/// SSX  sph-pln[1] sqDist=7.992002e+200
///      uv1=(0.0015707963267948967,-1.5700109286314992)
///      uv2=(1.999e+100,1.999e+100)
///      q1=(0.00078539711370610363,1.2336999159625383e-06,-0.99999969157487834)
///      q2=(1.999e+100,1.999e+100,5)
/// SSX gen-sph-pln nbExt=2 [1] sqDist=2.2138509695290867e+198 [2] sqDist=7.992002e+200
/// ```
///
/// The removed analytic `plane_sphere_extrema` arm answered 4.0 here. OCCT has
/// no analytic plane/sphere surface extrema at all:
/// `Extrema_ExtElSS::Perform(gp_Pln, gp_Sphere)` throws
/// `Standard_NotImplemented` (`Extrema_ExtElSS.cxx:62-69`).
#[test]
fn surface_surface_sphere_plane_is_gen_ext_ss() {
    let s1 = unit_sphere_at(GpPnt::new(0.0, 0.0, 0.0));
    let s2 = plane_z(5.0);
    let tol = occt_core::precision::PCONFUSION;

    // The raw engine stores both solutions: the grid/BFGS minimum and the grid
    // maximum.
    let gen = ExtremaGenExtSS::with_surfaces(&s1, &s2, 20, 20, tol, tol);
    assert!(gen.is_done());
    assert_eq!(gen.nb_ext(), 2, "GenExtSS stores two solutions");
    assert!(rel_close(gen.square_distance(1), 2.2138509695290867e198));
    assert!(rel_close(gen.square_distance(2), 7.992002e200));

    let ext = ExtremaExtSS::new(&s1, &s2, tol, tol);
    assert!(ext.is_done(), "Extrema_ExtSS not done");
    assert!(!ext.is_parallel(), "a sphere is not a plane");
    assert_eq!(ext.nb_ext(), 1, "the BFGS pair leaves the sphere's v window");
    assert!(rel_close(ext.square_distance(1), 7.992002e200));
    let (p1, p2) = ext.points(1);
    assert!(rel_close(p1.value().x(), 0.00078539711370610363));
    assert!(rel_close(p1.value().y(), 1.2336999159625383e-06));
    assert!(rel_close(p1.value().z(), -0.99999969157487834));
    assert!(rel_close(p2.value().x(), 1.999e100));
    assert!(rel_close(p2.value().y(), 1.999e100));
    assert_eq!(p2.value().z(), 5.0);

    let all = surface_surface_extrema_all(&s1, &s2);
    assert_eq!(all.len(), 1);
    assert!(rel_close(all[0].distance, 7.992002e200f64.sqrt()));
}

#[test]
fn surface_surface_planes_parallel() {
    let p1 = GeomPlane::new(occt_core::gp::GpPln::new(GpAx3::standard()));
    let p2 = plane_z(3.0);
    // Plane x plane is the analytic `Extrema_ExtElSS` arm; a parallel pair
    // records only the square distance (`Extrema_ExtElSS.cxx:46-53`).
    let ext = ExtremaExtSS::new(&p1, &p2, 0.0, 0.0);
    assert!(ext.is_done(), "Extrema_ExtSS not done");
    assert!(ext.is_parallel(), "parallel planes");
    assert_eq!(ext.nb_ext(), 1, "parallel planes report one square distance");
    assert!(
        (ext.square_distance(1) - 9.0).abs() < 1e-9,
        "parallel planes {}",
        ext.square_distance(1)
    );
}

#[test]
fn surface_surface_planes_coincident() {
    let p1 = GeomPlane::new(occt_core::gp::GpPln::new(GpAx3::standard()));
    let p2 = GeomPlane::new(occt_core::gp::GpPln::new(GpAx3::standard()));
    let ext = ExtremaExtSS::new(&p1, &p2, 0.0, 0.0);
    assert!(ext.is_done(), "Extrema_ExtSS not done");
    assert!(ext.is_parallel(), "coincident planes are parallel");
    assert_eq!(ext.nb_ext(), 1, "coincident planes report one square distance");
    assert!(ext.square_distance(1).abs() < 1e-9, "coincident planes");
}

fn paraboloid_plus_one() -> GeomBSplineSurface {
    let (nu, nv) = (3, 3);
    let points: Vec<Vec<GpPnt>> = (0..nu)
        .map(|i| {
            (0..nv)
                .map(|j| {
                    let u = i as f64 / (nu - 1) as f64;
                    let v = j as f64 / (nv - 1) as f64;
                    GpPnt::new(u, v, u * u + v * v + 1.0)
                })
                .collect()
        })
        .collect();
    fit_surface_grid(&points, 2, 2).unwrap()
}

#[test]
fn surface_surface_newton_paraboloid_plane_min() {
    // S1(u,v) = (u, v, u²+v²+1) (exact degree-2 B-spline), S2 = plane
    // z = 0.5 over [0,1]². Minimum distance 0.5 at (u1,v1)=(0,0).
    let s1 = paraboloid_plus_one();
    let s2 = GeomBSplineSurface::new(
        vec![
            vec![GpPnt::new(0.0, 0.0, 0.5), GpPnt::new(0.0, 1.0, 0.5)],
            vec![GpPnt::new(1.0, 0.0, 0.5), GpPnt::new(1.0, 1.0, 0.5)],
        ],
        bspline_surface_uniform_knots(2, 2, 1, 1).0,
        bspline_surface_uniform_knots(2, 2, 1, 1).1,
        1,
        1,
    )
    .unwrap();
    let all = surface_surface_extrema_all(&s1, &s2);
    assert!(!all.is_empty(), "no surface-surface extrema");
    let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
    assert!((min.distance - 0.5).abs() < 1e-6, "min {}", min.distance);
    // Orthogonality holds at the stored solution (`Extrema_FuncExtSS::Value`).
    let mut func = ExtremaFuncExtSS::with_surfaces(&s1, &s2);
    let mut uv = MathVector::new(1, 4);
    uv.set_value(1, min.u1);
    uv.set_value(2, min.v1.unwrap());
    uv.set_value(3, min.u2);
    uv.set_value(4, min.v2.unwrap());
    let mut f = MathVector::new(1, 4);
    func.value(&uv, &mut f);
    let residual: f64 = (1..=4).map(|i| f.value(i).abs()).sum();
    assert!(residual < 1e-3, "orthogonality F={residual}");
}

#[test]
fn sphere_sphere_gen_ext_ss_points() {
    let s1 = unit_sphere_at(GpPnt::new(0.0, 0.0, 0.0));
    let s2 = unit_sphere_at(GpPnt::new(4.0, 0.0, 0.0));
    let ext = ExtremaExtSS::new(&s1, &s2, 0.0, 0.0);
    assert!(ext.is_done(), "Extrema_ExtSS not done");
    // `Extrema_ExtSS::Points(N)` reads the `Extrema_GenExtSS` solution N.
    let (p1, p2) = ext.points(1);
    assert!((p1.value().square_distance(&p2.value()) - ext.square_distance(1)).abs() < 1e-9);
    assert!((p1.value().distance(&p2.value()) - 2.0).abs() < 1e-6, "{:?}", p1.value());
}
