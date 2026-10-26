//! `IntPolyh_Tools` (`IntPolyh_Tools.hxx`, `.cxx:24-194`).
//!
//! Sampling helpers: whether the sampling zone may be enlarged, the UV
//! sampling itself, the deflection of that sampling, and the point/normal grid
//! used by the shifted (advanced) intersection.

use occt_core::gp::GpPnt;
use occt_core::gp::GpVec;
use occt_core::precision::{Precision, REAL_SMALL};
use occt_geom::Surface;

use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};

use super::arrays::ArrayOfPointNormal;

/// `IsEnlargePossible(theSurf, theUEnlarge, theVEnlarge)` (`cxx:26-57`):
/// only BSpline and Bezier surfaces with finite, non-closed and non-periodic
/// parameters may be enlarged.
pub fn is_enlarge_possible(surface: &dyn Surface) -> (bool, bool) {
    let mut u_enlarge = false;
    let mut v_enlarge = false;

    let ty = surface_kind(surface);
    if ty == GeomAbsSurfaceType::BSplineSurface || ty == GeomAbsSurfaceType::BezierSurface {
        if !surface.is_u_closed() && !surface.is_u_periodic() {
            let (u0, u1) = surface.u_range();
            if !Precision::is_infinite(u0) && !Precision::is_infinite(u1) {
                u_enlarge = true;
            }
        }

        if !surface.is_v_closed() && !surface.is_v_periodic() {
            let (v0, v1) = surface.v_range();
            if !Precision::is_infinite(v0) && !Precision::is_infinite(v1) {
                v_enlarge = true;
            }
        }
    }

    (u_enlarge, v_enlarge)
}

/// `EnlargeZone(theSurf, u0, u1, v0, v1)` (`cxx:63-83`): grow by 1% of the
/// current span on each side.
fn enlarge_zone(surface: &dyn Surface, u0: &mut f64, u1: &mut f64, v0: &mut f64, v1: &mut f64) {
    let (to_enlarge_u, to_enlarge_v) = is_enlarge_possible(surface);
    if to_enlarge_u {
        let delta_u = 0.01 * (*u1 - *u0).abs();
        *u0 -= delta_u;
        *u1 += delta_u;
    }
    if to_enlarge_v {
        let delta_v = 0.01 * (*v1 - *v0).abs();
        *v0 -= delta_v;
        *v1 += delta_v;
    }
}

/// `MakeSampling(theSurf, theNbSU, theNbSV, theEnlargeZone, theUPars,
/// theVPars)` (`cxx:88-139`).
pub fn make_sampling(
    surface: &dyn Surface,
    nb_su: i32,
    nb_sv: i32,
    enlarge: bool,
) -> (Vec<f64>, Vec<f64>) {
    let (mut u0, mut u1) = surface.u_range();
    let (mut v0, mut v1) = surface.v_range();

    if enlarge {
        enlarge_zone(surface, &mut u0, &mut u1, &mut v0, &mut v1);
    }

    let nb_samples_u1 = nb_su - 1;
    let nb_samples_v1 = nb_sv - 1;

    let du = (u1 - u0) / f64::from(nb_samples_u1);
    let dv = (v1 - v0) / f64::from(nb_samples_v1);

    let mut u_pars = Vec::with_capacity(nb_su.max(0) as usize);
    for i in 0..nb_su {
        let mut a_u = u0 + f64::from(i) * du;
        if i == nb_samples_u1 {
            a_u = u1;
        }
        u_pars.push(a_u);
    }

    let mut v_pars = Vec::with_capacity(nb_sv.max(0) as usize);
    for i in 0..nb_sv {
        let mut a_v = v0 + f64::from(i) * dv;
        if i == nb_samples_v1 {
            a_v = v1;
        }
        v_pars.push(a_v);
    }

    (u_pars, v_pars)
}

/// `ComputeDeflection(theSurf, theUPars, theVPars)` (`cxx:144-151`): the
/// deflection over-estimation of `IntCurveSurface_ThePolyhedronOfHInter`
/// built on the given sampling.
pub fn compute_deflection(surface: &dyn Surface, u_pars: &[f64], v_pars: &[f64]) -> f64 {
    let polyhedron =
        crate::intcurvesurface::ThePolyhedronOfHInter::with_params(surface, u_pars, v_pars);
    polyhedron.deflection_over_estimation()
}

/// `FillArrayOfPointNormal(theSurf, theUPars, theVPars, thePoints)`
/// (`cxx:155-193`): the grid of sample points with their unit normals, in
/// U-major order.
pub fn fill_array_of_point_normal(
    surface: &dyn Surface,
    u_pars: &[f64],
    v_pars: &[f64],
    points: &mut ArrayOfPointNormal,
) {
    let nb_u = u_pars.len() as i32;
    let nb_v = v_pars.len() as i32;
    let mut i_cnt: i32 = 0;
    points.init(nb_u * nb_v);
    for i in 0..nb_u {
        let a_u = u_pars[i as usize];
        for j in 0..nb_v {
            let a_v = v_pars[j as usize];
            let (a_p, a_du, a_dv) = surface.d1(a_u, a_v);
            let mut a_v_norm = a_du.crossed(&a_dv);
            let a_length = a_v_norm.magnitude();
            if a_length > REAL_SMALL {
                a_v_norm = a_v_norm.divided(a_length);
            } else {
                a_v_norm = GpVec::new(0.0, 0.0, 0.0);
            }

            let pn = points.get_mut(i_cnt);
            pn.point = GpPnt::new(a_p.x(), a_p.y(), a_p.z());
            pn.normal = a_v_norm;

            i_cnt += 1;
        }
    }
    points.set_nb_items(i_cnt);
}
