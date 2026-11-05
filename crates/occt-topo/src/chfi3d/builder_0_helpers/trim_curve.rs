//! `ChFi3d_TrimCurve` - the trimmed curve between two 3D points.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4713-4768`
//! (OCCT 8.0.0). The analytic arms are `ElCLib::Parameter`; every other curve
//! goes through `GeomAPI_ProjectPointOnCurve tool(P, gc)` and requires exactly
//! one solution (`cxx:4749-4764`, `StdFail_NotDone`), which is the `Err` here.
//! `GeomAPI_ProjectPointOnCurve(P, C)` runs `Extrema_ExtPC` over
//! `C.FirstParameter() .. C.LastParameter()` (`GeomAPI_ProjectPointOnCurve.cxx:
//! 26-33`, `:58`), i.e. [`point_curve_extrema_all`].
//!
//! The result is `new Geom_TrimmedCurve(gc, uf, ul)` (`cxx:4767`); OCCT does not
//! reorder `uf` / `ul`, so a `ul < uf` trim is left as it is.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::GpPnt;
use occt_geom::curve::Curve;
use occt_geom::extrema_pc::{extrema_curve_type, point_curve_extrema_all, GeomAbsCurveType};
use occt_geom::trimmed::GeomTrimmedCurve;

/// `void ChFi3d_TrimCurve(const occ::handle<Geom_Curve>& gc, const gp_Pnt&
/// FirstP, const gp_Pnt& LastP, occ::handle<Geom_TrimmedCurve>& gtc)`
/// (`cxx:4713-4768`).
pub fn ch_fi3d_trim_curve(
    gc: &Arc<dyn Curve>,
    first_p: &GpPnt,
    last_p: &GpPnt,
) -> Result<GeomTrimmedCurve, String> {
    Ok(ch_fi3d_trim_curve_range(gc, first_p, last_p)?.0)
}

/// [`ch_fi3d_trim_curve`] plus the `uf` / `ul` it trims with, which are the
/// `FirstParameter()` / `LastParameter()` of the `Geom_TrimmedCurve` OCCT
/// returns (it stores them as given, without reordering).
pub fn ch_fi3d_trim_curve_range(
    gc: &Arc<dyn Curve>,
    first_p: &GpPnt,
    last_p: &GpPnt,
) -> Result<(GeomTrimmedCurve, f64, f64), String> {
    // `GeomAdaptor_Curve gac(gc); switch (gac.GetType())` (`cxx:4718-4719`).
    let (uf, ul);
    match extrema_curve_type(gc.as_ref()) {
        GeomAbsCurveType::Line => {
            let line = gc
                .gp_line()
                .ok_or_else(|| "ChFi3d_TrimCurve: the line curve has no gp_Lin".to_string())?;
            uf = clib::line_parameter(&line, first_p);
            ul = clib::line_parameter(&line, last_p);
        }
        GeomAbsCurveType::Circle => {
            let circ = gc
                .gp_circ()
                .ok_or_else(|| "ChFi3d_TrimCurve: the circle curve has no gp_Circ".to_string())?;
            uf = clib::circle_parameter(&circ.pos, first_p);
            ul = clib::circle_parameter(&circ.pos, last_p);
        }
        GeomAbsCurveType::Ellipse => {
            let elips = gc.gp_ellipse().ok_or_else(|| {
                "ChFi3d_TrimCurve: the ellipse curve has no gp_Elips".to_string()
            })?;
            uf = clib::parameter_elips(&elips, first_p);
            ul = clib::parameter_elips(&elips, last_p);
        }
        GeomAbsCurveType::Hyperbola => {
            let hypr = gc.gp_hyperbola().ok_or_else(|| {
                "ChFi3d_TrimCurve: the hyperbola curve has no gp_Hypr".to_string()
            })?;
            uf = clib::parameter_hypr(&hypr, first_p);
            ul = clib::parameter_hypr(&hypr, last_p);
        }
        GeomAbsCurveType::Parabola => {
            let parab = gc.gp_parabola().ok_or_else(|| {
                "ChFi3d_TrimCurve: the parabola curve has no gp_Parab".to_string()
            })?;
            uf = clib::parameter_parab(&parab, first_p);
            ul = clib::parameter_parab(&parab, last_p);
        }
        // `default:` (`cxx:4746-4765`).
        _ => {
            let all = point_curve_extrema_all(gc.as_ref(), first_p);
            if all.len() != 1 {
                return Err(
                    "ChFi3d_TrimCurve() - no projection results for the first point".to_string()
                );
            }
            uf = all[0].u1;
            let all = point_curve_extrema_all(gc.as_ref(), last_p);
            if all.len() != 1 {
                return Err(
                    "ChFi3d_TrimCurve() - no projection results for the second point".to_string()
                );
            }
            ul = all[0].u1;
        }
    }
    // `gtc = new Geom_TrimmedCurve(gc, uf, ul);` (`cxx:4767`).
    Ok((GeomTrimmedCurve::new(Arc::clone(gc), uf, ul), uf, ul))
}
