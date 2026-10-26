//! `Function_SetUVBounds`, `GeomAbs_Cone` arm.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx:287-461`
//! (OCCT 8.0.0).
//!
//! Two branches: a parabola / hyperbola / ellipse is read at its two ends and
//! its midpoint only (`cxx:312-331`), every other curve type is sampled along
//! its range (`cxx:332-460`) and the seam crossing is detected from the
//! running `U` deltas.

use occt_core::elib::slib;
use occt_core::gp::{GpLin, GpPnt};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;

use super::UvBounds;
use crate::curve::Curve;
use crate::projlib::compute_approx::function::{function_compute_step, is_equal};
use crate::projlib::projected_curve::curve_type;
use crate::surface::Surface;

const PI: f64 = std::f64::consts::PI;

/// `Function_SetUVBounds` cone arm (`cxx:287-461`).
///
/// `w1` / `w2` / `w` / `p1` / `p2` / `p` are the curve range and the three
/// sampled points; the apex correction (`cxx:295-309`) rewrites them in place.
#[allow(clippy::too_many_arguments)]
pub(super) fn set_uv_bounds(
    b: &mut UvBounds,
    curve: &dyn Curve,
    surface: &dyn Surface,
    w1: &mut f64,
    w2: &mut f64,
    w: &mut f64,
    p1: &mut GpPnt,
    p2: &mut GpPnt,
    p: &mut GpPnt,
) {
    let tol = precision::epsilon(1.0);
    let ptol = precision::PCONFUSION;
    let cone = surface.gp_cone().expect("GeomAbs_Cone");
    b.v_couture = false;

    // Calculation of cone parameters for P == ConeApex often produces wrong
    // values of U (`cxx:293-295`).
    let apex = cone.apex();
    if xyz_is_equal(&apex, p1, tol) {
        *w1 += ptol;
        *p1 = curve.d0(*w1);
    }
    if xyz_is_equal(&apex, p2, tol) {
        *w2 -= ptol;
        *p2 = curve.d0(*w2);
    }
    if xyz_is_equal(&apex, p, tol) {
        *w += ptol;
        *p = curve.d0(*w);
    }

    match curve_type(curve) {
        CurveType::Parabola | CurveType::Hyperbola | CurveType::Ellipse => {
            // `cxx:312-331`.
            let (u1, _v1) = slib::cone_parameters(&cone.pos, cone.radius, cone.semi_angle, p1);
            let (u2, _v2) = slib::cone_parameters(&cone.pos, cone.radius, cone.semi_angle, p2);
            let (u, _v) = slib::cone_parameters(&cone.pos, cone.radius, cone.semi_angle, p);
            b.u1 = u1.min(u2);
            b.u2 = u1.max(u2);
            if (u1 < u && u < u2) && !curve.is_closed() {
                b.u_couture = false;
            } else {
                b.u_couture = true;
                b.u2 = b.u1 + 2.0 * PI;
            }
        }
        _ => sample_uv_bounds(b, curve, surface, *w1, *w2, p1, p2, p),
    }
}

/// The `default:` curve arm (`cxx:332-460`): walk the range in `Step`
/// increments, unwrap the seam from the running `U1`, and open the window on
/// the seam when it wrapped.
fn sample_uv_bounds(
    b: &mut UvBounds,
    curve: &dyn Curve,
    surface: &dyn Surface,
    w1: f64,
    w2: f64,
    p1: &GpPnt,
    p2: &GpPnt,
    p: &GpPnt,
) {
    let cone = surface.gp_cone().expect("GeomAbs_Cone");
    let (u1, _v1) = slib::cone_parameters(&cone.pos, cone.radius, cone.semi_angle, p1);
    let (ul, _v2) = slib::cone_parameters(&cone.pos, cone.radius, cone.semi_angle, p2);
    let a_lin = GpLin::new(cone.axis());
    let r = (a_lin.distance(p1) + a_lin.distance(p2) + a_lin.distance(p)) / 3.0;

    b.u1 = u1;
    b.u2 = u1;
    let uf = u1;
    let mut u1 = u1;

    let step = if curve_type(curve) == CurveType::Line {
        let nbp = 3;
        (w2 - w1) / (nbp - 1) as f64
    } else {
        function_compute_step(curve, r)
    };

    let mut delta = 0.0f64;
    let mut pmin = w1;
    let mut pmax = w1;
    let mut dmax = 0.0f64;
    let isclandper = !curve.is_closed() && !curve.is_periodic();
    let mut is_first = true;

    let mut par = w1 + step;
    while par <= w2 {
        if !isclandper {
            par += step;
        }
        let q = curve.d0(par);
        let (mut u, _v) = slib::cone_parameters(&cone.pos, cone.radius, cone.semi_angle, &q);
        u += delta;
        let mut d = u - u1;
        if d > PI {
            if is_equal(u, 2.0 * PI, 1.0e-10)
                && (0.0..=PI).contains(&u1)
                && is_equal(u, ul, 1.0e-10)
                && !is_equal(uf, 0.0, 1.0e-10)
                && isclandper
            {
                u = 0.0;
            } else if is_first {
                // Protection against first-last point on seam.
                u1 = 2.0 * PI;
            } else if par + step >= w2 {
                u = 0.0;
            } else {
                delta -= 2.0 * PI;
            }
            u += delta;
            d = u - u1;
        } else if d < -PI {
            if is_equal(u, 0.0, 1.0e-10)
                && (PI..=2.0 * PI).contains(&u1)
                && is_equal(u, ul, 1.0e-10)
                && !is_equal(uf, 2.0 * PI, 1.0e-10)
                && isclandper
            {
                u = 2.0 * PI;
            } else if is_first {
                u1 = 0.0;
            } else if par + step >= w2 {
                u = 2.0 * PI;
            } else {
                delta += 2.0 * PI;
            }
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
        is_first = false;
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

/// `gp_XYZ::IsEqual` (`gp_XYZ.hxx`): component-wise, strict `<`. The cone arm
/// calls it on `Cone.Apex().XYZ()` and the sampled point's `XYZ()`.
fn xyz_is_equal(a: &GpPnt, b: &GpPnt, tol: f64) -> bool {
    (a.x() - b.x()).abs() < tol && (a.y() - b.y()).abs() < tol && (a.z() - b.z()).abs() < tol
}
