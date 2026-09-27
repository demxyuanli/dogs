//! Port of `Extrema_GCurveLocator` (`Extrema_GCurveLocator.hxx:32-130`) as
//! instantiated by `Geom2dInt_TheCurveLocatorOfTheProjPCurOfGInter`
//! (`Geom2dInt_TheCurveLocatorOfTheProjPCurOfGInter.hxx:26-28`:
//! `Extrema_GCurveLocator<Adaptor2d_Curve2d, Geom2dInt_Geom2dCurveTool,
//! Extrema_POnCurv2d, gp_Pnt2d>`). Returns `(parameter, point)` in place of
//! `Extrema_POnCurv2d::SetValues(parameter, point)`.

use occt_core::gp::GpPnt2d;

use crate::curve::Curve2d;
use super::curve_tool;

/// `Locate(P, C, NbU)` (`hxx:45-72`): the sample `C(u_i)`, i=1..NbU-1, closest
/// to `P`. OCCT throws `Standard_OutOfRange` for `NbU < 2`; the port requires
/// the same.
pub fn locate(p: &GpPnt2d, c: &dyn Curve2d, nb_u: i32) -> (f64, GpPnt2d) {
    assert!(nb_u >= 2, "Extrema_GCurveLocator::Locate: NbU < 2");
    let mut u = curve_tool::first_parameter(c);
    let pas_u = (curve_tool::last_parameter(c) - u) / f64::from(nb_u - 1);
    let mut dist2_min = f64::MAX;
    let mut u_min = 0.0f64;
    let mut pnt_min = GpPnt2d::new(0.0, 0.0);
    let mut no_sample = 1;
    while no_sample < nb_u {
        let pt = curve_tool::value(c, u);
        let dist2 = pt.square_distance(p);
        if dist2 < dist2_min {
            dist2_min = dist2;
            u_min = u;
            pnt_min = pt;
        }
        no_sample += 1;
        u += pas_u;
    }
    (u_min, pnt_min)
}

/// `Locate(P, C, NbU, Umin, Usup)` (`hxx:84-129`): same as [`locate`] but the
/// sampling window is intersected with `[Umin, Usup]`.
pub fn locate_range(
    p: &GpPnt2d,
    c: &dyn Curve2d,
    nb_u: i32,
    u_min_in: f64,
    u_sup_in: f64,
) -> (f64, GpPnt2d) {
    assert!(nb_u >= 2, "Extrema_GCurveLocator::Locate: NbU < 2");
    let u_inf = curve_tool::first_parameter(c);
    let u_last = curve_tool::last_parameter(c);
    let u1 = u_inf.min(u_last);
    let u2 = u_inf.max(u_last);
    let mut u11 = u_min_in.min(u_sup_in);
    let mut u12 = u_min_in.max(u_sup_in);
    let real_epsilon = f64::EPSILON; // RealEpsilon()
    if u11 < u1 - real_epsilon {
        u11 = u1;
    }
    if u12 > u2 + real_epsilon {
        u12 = u2;
    }
    let mut u = u11;
    let pas_u = (u12 - u) / f64::from(nb_u - 1);
    let mut dist2_min = f64::MAX;
    let mut u_min = 0.0f64;
    let mut pnt_min = GpPnt2d::new(0.0, 0.0);
    let mut no_sample = 1;
    while no_sample < nb_u {
        let pt = curve_tool::value(c, u);
        let dist2 = pt.square_distance(p);
        if dist2 < dist2_min {
            dist2_min = dist2;
            u_min = u;
            pnt_min = pt;
        }
        no_sample += 1;
        u += pas_u;
    }
    (u_min, pnt_min)
}
