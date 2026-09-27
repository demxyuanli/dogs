//! Port of `Geom2dInt_TheProjPCurOfGInter`
//! (`Geom2dInt_TheProjPCurOfGInter_0.cxx:27-79`): the `ProjectOnPCurveTool`
//! of the `Geom2dInt` instantiation of `IntImpParGen_Intersector`.
//!
//! `FindParameter` locates the closest sample (`Extrema_GCurveLocator`) and
//! refines it with `Extrema_GenLocateExtPC`.

use occt_core::gp::GpPnt2d;

use crate::curve::Curve2d;
use super::curve_locator;
use super::curve_tool;
use super::gen_locate_ext_pc::GenLocateExtPC;

/// `FindParameter(C, P, LowParameter, HighParameter, Tol)` (`cxx:27-65`).
pub fn find_parameter_range(
    c: &dyn Curve2d,
    p: &GpPnt2d,
    low_parameter: f64,
    high_parameter: f64,
    _tol: f64,
) -> f64 {
    // cxx:34-35.
    let nb_pts = curve_tool::nb_samples_curve(c) as i32;
    let the_eps_x = curve_tool::eps_x();
    // cxx:38-44.
    let (defaultparam, _p_on_c) =
        curve_locator::locate_range(p, c, nb_pts, low_parameter, high_parameter);
    // cxx:45: Geom2dInt_TheLocateExtPCOfTheProjPCurOfGInter Loc(P, C, defaultparam, theEpsX).
    let mut loc = GenLocateExtPC::default();
    loc.initialize(
        c,
        curve_tool::first_parameter(c),
        curve_tool::last_parameter(c),
        the_eps_x,
    );
    loc.perform(p, defaultparam);
    if !loc.is_done() {
        return defaultparam; // cxx:47-51
    }
    if !loc.is_min() {
        return defaultparam; // cxx:54-58
    }
    loc.point(1).parameter() // cxx:61
}

/// `FindParameter(C, P, Tol)` (`cxx:67-79`).
pub fn find_parameter(c: &dyn Curve2d, p: &GpPnt2d, tol: f64) -> f64 {
    find_parameter_range(
        c,
        p,
        curve_tool::first_parameter(c),
        curve_tool::last_parameter(c),
        tol,
    )
}
