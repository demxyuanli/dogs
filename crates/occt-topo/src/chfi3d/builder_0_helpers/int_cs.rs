//! `ChFi3d_IntCS` - the fast curve / surface intersection of the corner code.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4296-4380`
//! (OCCT 8.0.0).
//!
//! `IntCurveSurface_HInter::Perform(C, S)` answers a set of points; the body
//! keeps the one whose curve parameter lies inside the curve's own range, whose
//! `(u, v)` lies inside the adaptor's parameter box after folding by whole
//! periods, and then picks it by the requested end:
//!
//! - `wc < -1.e100` (`keepfirst`) keeps the smallest curve parameter,
//! - `wc > 1.e100` (`keeplast`) keeps the largest,
//! - otherwise the point nearest to the incoming `wc`.
//!
//! Translation notes:
//! - the port's `perform_curve_surface` is `IntCurveSurface_HInter`; an `Err`
//!   there is the not-done state, answered as `false` like `!IsDone()`,
//! - `ChFi3d_InPeriod(up, u1, u1 + S->UPeriod(), 1.e-8)` uses the same helper
//!   as the `u` / `v` folds of [`crate::chfi3d::builder_helpers::in_period`],
//! - `p2dS.SetCoord(up, vp)` and `wc = pint.W()` are the two out-parameters; the
//!   port keeps OCCT's in/out use of `wc` without touching it on the `false`
//!   paths.

use occt_core::gp::GpPnt2d;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;

use crate::chfi3d::builder_helpers::in_period;
use crate::intcurvesurface::perform_curve_surface;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_IntCS(S, C, p2dS, wc)` (`cxx:4302-4380`).
///
/// `s` is the trimmed adaptor OCCT passes as `S`: its `FirstUParameter` /
/// `LastUParameter` box is the acceptance window, and the surface behind it is
/// what `Perform` intersects the curve with.
pub fn ch_fi3d_int_cs(
    s: &SurfaceAdaptor,
    c: &dyn Curve,
    p2d_s: &mut GpPnt2d,
    wc: &mut f64,
) -> bool {
    // `double uf = C->FirstParameter(), ul = C->LastParameter(); double u1 =
    // S->FirstUParameter() ...;` (`cxx:4309-4313`).
    let uf = c.first_parameter();
    let ul = c.last_parameter();
    let u1 = s.first_u_parameter();
    let u2 = s.last_u_parameter();
    let v1 = s.first_v_parameter();
    let v2 = s.last_v_parameter();

    // `Intersection.Perform(C, S);` (`cxx:4314`).
    let Some(surface) = s.surface() else {
        return false;
    };
    let Ok(intersection) = perform_curve_surface(
        c,
        surface.as_ref(),
        (uf, ul),
        (u1, v1, u2, v2),
    ) else {
        return false;
    };

    // `bool keepfirst = (wc < -1.e100), keeplast = (wc > 1.e100); double temp
    // = 0.;` (`cxx:4315-4317`).
    let keepfirst = *wc < -1.0e100;
    let keeplast = *wc > 1.0e100;
    let mut temp = 0.0;
    if keepfirst {
        temp = 1.0e100;
    }
    if keeplast {
        temp = -1.0e100;
    }
    // `double dist = 2.e100;` (`cxx:4324`).
    let mut dist = 2.0e100;

    // `int nbp = Intersection.NbPoints(), i, isol = 0;` (`cxx:4327`); the port
    // keeps the 1-based index and maps it to 0-based only at the read.
    let nbp = intersection.nb_points();
    let mut isol = 0usize;
    for i in 1..=nbp {
        let pint = intersection.point(i - 1);
        // `double up = pint.U(); double vp = pint.V();` (`cxx:4331-4332`).
        let mut up = pint.u;
        let mut vp = pint.v;
        // `if (S->IsUPeriodic()) { up = ChFi3d_InPeriod(up, u1, u1 +
        // S->UPeriod(), 1.e-8); }` and the `v` twin (`cxx:4333-4340`).
        if surface.is_u_periodic() {
            up = in_period(up, u1, u1 + surface.u_period(), 1.0e-8);
        }
        if surface.is_v_periodic() {
            vp = in_period(vp, v1, v1 + surface.v_period(), 1.0e-8);
        }
        // `if (uf <= pint.W() && ul >= pint.W() && u1 <= up && u2 >= up &&
        // v1 <= vp && v2 >= vp)` (`cxx:4341-4342`).
        if uf <= pint.param && ul >= pint.param && u1 <= up && u2 >= up && v1 <= vp && v2 >= vp
        {
            if keepfirst && pint.param < temp {
                temp = pint.param;
                isol = i;
            } else if keeplast && pint.param > temp {
                temp = pint.param;
                isol = i;
            } else if (pint.param - *wc).abs() < dist {
                dist = (pint.param - *wc).abs();
                isol = i;
            }
        }
    }

    // `if (isol == 0) { return false; }` (`cxx:4360-4363`).
    if isol == 0 {
        return false;
    }

    // `pint = Intersection.Point(isol);` and the two folds again
    // (`cxx:4364-4375`).
    let pint = intersection.point(isol - 1);
    let mut up = pint.u;
    let mut vp = pint.v;
    if surface.is_u_periodic() {
        up = in_period(up, u1, u1 + surface.u_period(), 1.0e-8);
    }
    if surface.is_v_periodic() {
        vp = in_period(vp, v1, v1 + surface.v_period(), 1.0e-8);
    }
    // `p2dS.SetCoord(up, vp); wc = pint.W(); return true;`
    // (`cxx:4376-4378`).
    *p2d_s = GpPnt2d::new(up, vp);
    *wc = pint.param;
    true
}
