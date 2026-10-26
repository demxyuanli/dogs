//! `Function_SetUVBounds`, `GeomAbs_Sphere` arm.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx:593-889`
//! (OCCT 8.0.0).
//!
//! A `GeomAbs_Circle` is analysed analytically against the sphere's seam
//! meridian (`cxx:596-810`): the number of intersections of the circle with
//! the seam decides whether the 2d result is closed, and a circle parallel to
//! the equator at the seam height raises `VCouture`. Every other curve type is
//! sampled (`cxx:812-889`).
//!
//! `myV1` / `myV2` are overwritten with `-1e100` / `1e100` right before the
//! arm returns (`cxx:806-808`, the `box+sphere >>` marker), which makes the
//! sphere window unbounded in `V` for `Function_Value`.

use occt_core::elib::slib;
use occt_core::gp::{GpAx3, GpPln, GpPnt, GpTrsf};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;

use super::UvBounds;
use crate::curve::Curve;
use crate::projlib::compute_approx::function::function_compute_step;
use crate::projlib::projected_curve::curve_type;
use crate::surface::Surface;

const PI: f64 = std::f64::consts::PI;

/// `Function_SetUVBounds` sphere arm (`cxx:593-889`).
pub(super) fn set_uv_bounds(
    b: &mut UvBounds,
    curve: &dyn Curve,
    surface: &dyn Surface,
    w1: f64,
    w2: f64,
    p1: &GpPnt,
    p2: &mut GpPnt,
    p: &GpPnt,
) {
    b.v_couture = false;
    let sp = surface.gp_sphere().expect("GeomAbs_Sphere");

    if curve_type(curve) == CurveType::Circle {
        circle_arm(b, curve, &sp, w1, w2, p1, p2, p);
        // `box+sphere >>` (`cxx:806-808`).
        b.v1 = -1.0e+100;
        b.v2 = 1.0e+100;
        return;
    }

    // `cxx:812-889`: same sampling as the cylinder arm, `[0, 2*PI]` window
    // normalized to the full period instead of opened at the seam.
    let (u1, _v1) = slib::sphere_parameters(&sp.pos, p1);
    let r = sp.radius;
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
        let (mut u, _v) = slib::sphere_parameters(&sp.pos, &q);
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
        b.u1 = 0.0;
        b.u2 = 2.0 * PI;
        b.u_couture = false;
    } else {
        let u = (b.u1 + b.u2) / 2.0;
        b.u1 = u - PI;
        b.u2 = u + PI;
        b.u_couture = true;
    }

    b.v_couture = false;
}

/// `cxx:596-810`: the `GeomAbs_Circle` analysis.
#[allow(clippy::too_many_arguments)]
fn circle_arm(
    b: &mut UvBounds,
    curve: &dyn Curve,
    sp: &occt_core::gp::GpSphere,
    w1: f64,
    w2: f64,
    p1: &GpPnt,
    p2: &mut GpPnt,
    p: &GpPnt,
) {
    b.u_couture = true;

    // `cxx:619-625`: the circle carried into the sphere's frame, then its
    // plane's coefficients `A x + B y + C z + D = 0`.
    let mut circle = curve.gp_circ().expect("GeomAbs_Circle");
    let mut trsf = GpTrsf::identity();
    trsf.set_transformation(&sp.pos);
    circle.transform(&trsf);
    let r = sp.radius;
    let pln = GpPln::new(GpAx3::from_ax2(&circle.position()));
    let (ca, cb, cc, cd) = pln.coefficients();
    let _ = cb;
    let tol = 1.0e-10;

    // `cxx:617-677`: how many times the seam is crossed. 0 or 2 solutions ->
    // the pcurve is closed and lives in `[Uc - PI, Uc + PI]`; 1 -> the pcurve
    // is open and lives in `[U1, U1 +/- 2*PI]`.
    //
    // The system solved is
    //   x^2 + y^2 + z^2     = R^2   (1)
    //   A x + B y + C z + D = 0     (2)
    //   x > 0                       (3)
    //   y = 0                       (4)
    // where (1)(2) is the circle and (1)(3)(4) the seam.
    let mut nb_solutions = 0i32;
    if cc.abs() < tol {
        if ca.abs() > tol && (cd / ca) < 0.0 {
            if (r - (cd / ca).abs()) > tol {
                nb_solutions = 2;
            } else if (r - (cd / ca).abs()).abs() < tol {
                nb_solutions = 1;
            } else {
                nb_solutions = 0;
            }
        }
    } else {
        let mut delta = r * r * (ca * ca + cc * cc) - cd * cd;
        delta *= cc * cc;
        if delta.abs() < tol * tol {
            if ca * cd > 0.0 {
                nb_solutions = 1;
            }
        } else if delta > 0.0 {
            delta = delta.sqrt();
            if -ca * cd + delta > tol {
                nb_solutions += 1;
            }
            if -ca * cd - delta > tol {
                nb_solutions += 1;
            }
        }
    }

    // `cxx:679-705`: `U` of the first end, with the seam ambiguity resolved
    // by looking one `dt` inside the range.
    let eps = 10.0 * precision::epsilon(1.0);
    let dt = precision::PCONFUSION.max(0.01 * (w2 - w1));
    let (mut u1, _v1) = slib::sphere_parameters(&sp.pos, p1);
    if u1.abs() < eps {
        let pd = curve.d0(w1 + dt);
        let (ud, _vd) = slib::sphere_parameters(&sp.pos, &pd);
        if (u1 - ud).abs() > PI {
            u1 = 2.0 * PI;
        }
    } else if (2.0 * PI - u1).abs() < eps {
        let pd = curve.d0(w1 + dt);
        let (ud, _vd) = slib::sphere_parameters(&sp.pos, &pd);
        if (u1 - ud).abs() > PI {
            u1 = 0.0;
        }
    }

    // `cxx:706-729`: `U2` of the second end, seam-adjusted the same way. OCCT
    // then overwrites `U2` at `cxx:733` without reading the adjusted value, so
    // this block is dead; it is kept for control-flow parity.
    let (mut u2, _v2) = slib::sphere_parameters(&sp.pos, p2);
    if u2.abs() < eps {
        let pd = curve.d0(w2 - dt);
        let (ud, _vd) = slib::sphere_parameters(&sp.pos, &pd);
        if (u2 - ud).abs() > PI {
            u2 = 2.0 * PI;
        }
    } else if (2.0 * PI - u2).abs() < eps {
        let pd = curve.d0(w2 - dt);
        let (ud, _vd) = slib::sphere_parameters(&sp.pos, &pd);
        if (u2 - ud).abs() > PI {
            u2 = 0.0;
        }
    }
    let _ = u2;

    // `cxx:730-733`: `UU` of the midpoint is written but never read.
    let (_uu, _v3) = slib::sphere_parameters(&sp.pos, p);
    *p2 = curve.d0(w1 + PI / 8.0);
    let (u2, _v4) = slib::sphere_parameters(&sp.pos, p2);

    if nb_solutions == 1 {
        // `cxx:734-767`.
        if (u1 - u2).abs() > PI {
            // crossing the seam
            if u1 > PI {
                b.u1 = u1;
                b.u2 = u1 + 2.0 * PI;
            } else {
                b.u2 = u1;
                b.u1 = u1 - 2.0 * PI;
            }
        } else if u1 > u2 {
            // not crossing the seam
            b.u2 = u1;
            b.u1 = u1 - 2.0 * PI;
        } else {
            b.u1 = u1;
            b.u2 = u1 + 2.0 * PI;
        }
    } else {
        // 0 or 2 solutions: `cxx:764-768`.
        let center = circle.location();
        let (u, _v) = slib::sphere_parameters(&GpAx3::standard(), &center);
        b.u1 = u - PI;
        b.u2 = u + PI;
    }

    // `cxx:771-804`: evaluate the `VCouture`.
    if cc == 0.0 || ((cd / cc).abs() - r).abs() > 1.0e-10 {
        b.v_couture = false;
    } else {
        b.v_couture = true;
        b.u_couture = true;
        if cd / cc < 0.0 {
            b.v1 = -PI / 2.0;
            b.v2 = 3.0 * PI / 2.0;
        } else {
            b.v1 = -3.0 * PI / 2.0;
            b.v2 = PI / 2.0;
        }

        // If `P1` is at a vertex (`+/- R` on the axis) the test is not
        // significant.
        let pp = p1.transformed(&trsf);
        if (pp.x() * pp.x() + pp.y() * pp.y() + pp.z() * pp.z() - r * r).abs() < tol {
            let center = circle.location();
            let (u, _v) = slib::sphere_parameters(&GpAx3::standard(), &center);
            b.u1 = u - PI;
            b.u2 = u + PI;
            b.v_couture = false;
        }
    }
}
