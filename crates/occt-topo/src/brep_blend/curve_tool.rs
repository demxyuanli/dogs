//! `BRepBlend_HCurveTool` / `BRepBlend_HCurve2dTool` — the sampling counts the
//! walking builders ask of a curve adaptor.
//!
//! Source (OCCT 8.0.0), `ModelingAlgorithms/TKFillet/BRepBlend/`:
//! `BRepBlend_HCurveTool.hxx:41-121` / `.cxx:25-57` and
//! `BRepBlend_HCurve2dTool.hxx:43-125` / `.cxx:25-57`.
//!
//! Only `NbSamples` has a body in either class (`cxx:25-57` each); every other
//! member is `C->Method(...)` with no branch, which in this port is the
//! `occt_geom::curve::Curve` / `occt_geom2d::curve::Curve2d` trait call itself,
//! so no forwarding class exists. This module keeps the one translated body.
//!
//! `NbSamples` translates `Geom_BSplineCurve::NbKnots()` (`hxx` of
//! `Geom_BSplineCurve`: the count of *distinct* knots) through
//! [`occt_core::bspl::unperiodize::distinct_knots_and_mults`], because the
//! adaptor trait exposes the flat knot vector (`Curve::bspline_knots`).

use occt_core::bspl::unperiodize::distinct_knots_and_mults;
use occt_geom::curve::Curve;
use occt_geom2d::Curve2d;

/// `static double nbsOther = 10.0` (`cxx:31` in both classes).
const NBS_OTHER: f64 = 10.0;
/// `nbs = 2` for a line (`cxx:35`).
const NBS_LINE: f64 = 2.0;
/// `nbs = 3 + C->Bezier()->NbPoles()` (`cxx:39`).
const NBS_BEZIER_BASE: f64 = 3.0;
/// `if (nbs < 2.0) nbs = 2.0;` (`cxx:48-51`).
const NBS_MIN: f64 = 2.0;
/// `if (nbs > 50) nbs = 50;` (`cxx:53-56`).
const NBS_MAX: f64 = 50.0;

/// `BRepBlend_HCurveTool::NbSamples(C, U0, U1)` (`cxx:25-57`).
pub fn nb_samples(c: &dyn Curve, u0: f64, u1: f64) -> i32 {
    let mut nbs = NBS_OTHER;

    if c.is_line() {
        nbs = NBS_LINE;
    } else if let Some(poles) = c.bezier_poles() {
        nbs = NBS_BEZIER_BASE + poles.len() as f64;
    } else if let Some(knots) = c.bspline_knots() {
        let (distinct, _) = distinct_knots_and_mults(knots);
        nbs = distinct.len() as f64;
        nbs *= c.nurbs_degree().unwrap_or(0) as f64;
        // `BSC->LastParameter() - BSC->FirstParameter()`.
        let (first, last) = match (knots.first(), knots.last()) {
            (Some(f), Some(l)) => (*f, *l),
            _ => (0.0, 0.0),
        };
        nbs *= last - first;
        nbs /= u1 - u0;
        if nbs < NBS_MIN {
            nbs = NBS_MIN;
        }
    }
    if nbs > NBS_MAX {
        nbs = NBS_MAX;
    }
    nbs as i32
}

/// `BRepBlend_HCurve2dTool::NbSamples(C, U0, U1)` (`cxx:25-57`).
pub fn nb_samples2d(c: &dyn Curve2d, u0: f64, u1: f64) -> i32 {
    let mut nbs = NBS_OTHER;

    if c.is_line() {
        nbs = NBS_LINE;
    } else if c.is_bezier2d() {
        nbs = NBS_BEZIER_BASE + c.bezier_nb_poles().unwrap_or(0) as f64;
    } else if c.is_bspline2d() {
        nbs = c.bspline_nb_knots().unwrap_or(0) as f64;
        nbs *= c.bspline_degree().unwrap_or(0) as f64;
        // `BSC->LastParameter() - BSC->FirstParameter()`.
        let (first, last) = match c.bspline_knots2d() {
            Some(knots) => match (knots.first(), knots.last()) {
                (Some(f), Some(l)) => (*f, *l),
                _ => (0.0, 0.0),
            },
            None => (0.0, 0.0),
        };
        nbs *= last - first;
        nbs /= u1 - u0;
        if nbs < NBS_MIN {
            nbs = NBS_MIN;
        }
    }
    if nbs > NBS_MAX {
        nbs = NBS_MAX;
    }
    nbs as i32
}
