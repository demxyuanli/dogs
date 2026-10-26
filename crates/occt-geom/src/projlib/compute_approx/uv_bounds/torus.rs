//! `Function_SetUVBounds`, `GeomAbs_Torus` arm.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx:890-1012`
//! (OCCT 8.0.0).
//!
//! The torus is the only surface with two independent seam windows: `U` and
//! `V` are both sampled (`cxx:902-962`) and both opened if they wrapped.

use occt_core::elib::slib;
use occt_core::gp::GpPnt;
use occt_core::precision;

use super::UvBounds;
use crate::curve::Curve;
use crate::projlib::compute_approx::function::function_compute_step;
use crate::surface::Surface;

const PI: f64 = std::f64::consts::PI;

/// `Function_SetUVBounds` torus arm (`cxx:890-1012`).
pub(super) fn set_uv_bounds(
    b: &mut UvBounds,
    curve: &dyn Curve,
    surface: &dyn Surface,
    w1: f64,
    w2: f64,
    p1: &GpPnt,
) {
    let tr = surface.gp_torus().expect("GeomAbs_Torus");
    let (u1, v1) = slib::torus_parameters(&tr.pos, tr.major_radius, tr.minor_radius, p1);
    let r = tr.minor_radius;
    let step = function_compute_step(curve, r);
    b.u1 = u1;
    b.u2 = u1;
    b.v1 = v1;
    b.v2 = v1;

    let mut u1 = u1;
    let mut v1 = v1;
    let mut delta_u = 0.0f64;
    let mut delta_v = 0.0f64;
    let mut dmax_u = 0.0f64;
    let mut dmax_v = 0.0f64;
    let mut pmin_u = w1;
    let mut pmax_u = w1;
    let mut pmin_v = w1;
    let mut pmax_v = w1;
    let plim = w2 + 0.1 * step;

    let mut par = w1 + step;
    while par <= plim {
        let q = curve.d0(par);
        let (mut u, mut v) = slib::torus_parameters(&tr.pos, tr.major_radius, tr.minor_radius, &q);
        u += delta_u;
        v += delta_v;
        let mut du = u - u1;
        let mut dv = v - v1;
        if du > PI {
            u -= delta_u;
            delta_u -= 2.0 * PI;
            u += delta_u;
            du = u - u1;
        } else if du < -PI {
            u -= delta_u;
            delta_u += 2.0 * PI;
            u += delta_u;
            du = u - u1;
        }
        if dv > PI {
            v -= delta_v;
            delta_v -= 2.0 * PI;
            v += delta_v;
            dv = v - v1;
        } else if dv < -PI {
            v -= delta_v;
            delta_v += 2.0 * PI;
            v += delta_v;
            dv = v - v1;
        }
        dmax_u = dmax_u.max(du.abs());
        dmax_v = dmax_v.max(dv.abs());
        if u < b.u1 {
            b.u1 = u;
            pmin_u = par;
        }
        if u > b.u2 {
            b.u2 = u;
            pmax_u = par;
        }
        if v < b.v1 {
            b.v1 = v;
            pmin_v = par;
        }
        if v > b.v2 {
            b.v2 = v;
            pmax_v = par;
        }
        u1 = u;
        v1 = v;
        par += step;
    }

    if (pmin_u - w1).abs() > precision::PCONFUSION && (pmin_u - w2).abs() > precision::PCONFUSION {
        b.u1 -= dmax_u * 0.5;
    }
    if (pmax_u - w1).abs() > precision::PCONFUSION && (pmax_u - w2).abs() > precision::PCONFUSION {
        b.u2 += dmax_u * 0.5;
    }
    if (pmin_v - w1).abs() > precision::PCONFUSION && (pmin_v - w2).abs() > precision::PCONFUSION {
        b.v1 -= dmax_v * 0.5;
    }
    if (pmax_v - w1).abs() > precision::PCONFUSION && (pmax_v - w2).abs() > precision::PCONFUSION {
        b.v2 += dmax_v * 0.5;
    }

    if (0.0..=2.0 * PI).contains(&b.u1) && (0.0..=2.0 * PI).contains(&b.u2) {
        b.u1 = 0.0;
        b.u2 = 2.0 * PI;
        b.u_couture = false;
    } else {
        let u = (b.u1 + b.u2) / 2.0;
        b.u1 = u - PI;
        b.u2 = u + PI;
        b.u_couture = true;
    }
    if (0.0..=2.0 * PI).contains(&b.v1) && (0.0..=2.0 * PI).contains(&b.v2) {
        b.v_couture = false;
    } else {
        let v = (b.v1 + b.v2) / 2.0;
        b.v1 = v - PI;
        b.v2 = v + PI;
        b.v_couture = true;
    }
}
