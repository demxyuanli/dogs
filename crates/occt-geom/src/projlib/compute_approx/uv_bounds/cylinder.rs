//! `Function_SetUVBounds`, `GeomAbs_Cylinder` arm.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx:463-591`
//! (OCCT 8.0.0).
//!
//! An ellipse is read at its two ends and midpoint (`cxx:468-521`); every other
//! curve type is sampled across `[W1, W2 + 0.1 * Step]` (`cxx:522-590`). The
//! window is a pure `U` window — `VCouture` is always cleared.

use occt_core::elib::slib;
use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;

use super::UvBounds;
use crate::curve::Curve;
use crate::projlib::compute_approx::function::function_compute_step;
use crate::projlib::projected_curve::curve_type;
use crate::surface::Surface;

const PI: f64 = std::f64::consts::PI;

/// `Function_SetUVBounds` cylinder arm (`cxx:463-591`).
pub(super) fn set_uv_bounds(
    b: &mut UvBounds,
    curve: &dyn Curve,
    surface: &dyn Surface,
    w1: f64,
    w2: f64,
    p1: &GpPnt,
    p2: &GpPnt,
    p: &GpPnt,
) {
    let cyl = surface.gp_cylinder().expect("GeomAbs_Cylinder");
    b.v_couture = false;

    if curve_type(curve) == CurveType::Ellipse {
        // `cxx:468-521`.
        let (u1, _v1) = slib::cylinder_parameters(&cyl.pos, p1);
        let (u2, _v2) = slib::cylinder_parameters(&cyl.pos, p2);
        let (u_mid, _v) = slib::cylinder_parameters(&cyl.pos, p);
        b.u1 = u1.min(u2);
        b.u2 = u1.max(u2);

        if !curve.is_closed() {
            if b.u1 < u_mid && u_mid < b.u2 {
                let u = (b.u1 + b.u2) / 2.0;
                b.u1 = u - PI;
                b.u2 = u + PI;
            } else {
                let u = (b.u1 + b.u2) / 2.0;
                if b.u1 < u {
                    b.u1 = u - 2.0 * PI;
                    b.u2 = u;
                } else {
                    b.u1 = u;
                    b.u2 = u + 2.0 * PI;
                }
            }
            b.u_couture = true;
        } else {
            // `cxx:502-520`: the derivative at `W1` against `D1U` at the
            // `(U1, U2)` pair the ellipse ends produced tells the direction.
            let (_, t) = curve.d1(w1);
            let (_, d1u, _d1v) = surface.d1(u1, u2);
            let d_u = t.dot(&d1u);
            b.u_couture = true;
            if d_u > 0.0 {
                b.u2 = b.u1 + 2.0 * PI;
            } else {
                b.u2 = b.u1;
                b.u1 -= 2.0 * PI;
            }
        }
        return;
    }

    // `cxx:522-590`.
    let (u1, _v1) = slib::cylinder_parameters(&cyl.pos, p1);
    let r = cyl.radius;
    let step = function_compute_step(curve, r);
    b.u1 = u1;
    b.u2 = u1;
    let mut u1 = u1;

    let mut delta = 0.0f64;
    let mut dmax = 0.0f64;
    let mut pmin = w1;
    let mut pmax = w1;
    let plim = w2 + 0.1 * step;

    let mut par = w1 + step;
    while par <= plim {
        let q = curve.d0(par);
        let (mut u, _v) = slib::cylinder_parameters(&cyl.pos, &q);
        u += delta;
        let mut d = u - u1;
        if d > PI {
            u -= delta;
            delta -= 2.0 * PI;
            u += delta;
            d = u - u1;
        } else if d < -PI {
            u -= delta;
            delta += 2.0 * PI;
            u += delta;
            d = u - u1;
        }
        dmax = dmax.max(d.abs());
        if u < b.u1 {
            b.u1 = u;
            pmin = par;
        }
        if u > b.u2 {
            b.u2 = u;
            pmax = par;
        }
        u1 = u;
        par += step;
    }

    if (pmin - w1).abs() > precision::PCONFUSION && (pmin - w2).abs() > precision::PCONFUSION {
        b.u1 -= dmax * 0.5;
    }
    if (pmax - w1).abs() > precision::PCONFUSION && (pmax - w2).abs() > precision::PCONFUSION {
        b.u2 += dmax * 0.5;
    }

    if (0.0..=2.0 * PI).contains(&b.u1) && (0.0..=2.0 * PI).contains(&b.u2) {
        b.u_couture = false;
    } else {
        let u = (b.u1 + b.u2) / 2.0;
        b.u1 = u - PI;
        b.u2 = u + PI;
        b.u_couture = true;
    }
}
