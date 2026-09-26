//! `Extrema_Curve2dTool` (`Extrema_Curve2dTool.hxx:38-157`) and the
//! `Geom2dAdaptor_Curve` queries it forwards to, as free functions over the
//! port's `dyn Curve2d`.
//!
//! OCCT's `Extrema_Curve2dTool` is a static facade whose methods are one-line
//! calls on `Adaptor2d_Curve2d` (`Extrema_Curve2dTool.hxx:43-156`); the port's
//! `Curve2d` trait already covers `FirstParameter`/`LastParameter`/
//! `IsPeriodic`/`Period`/`Value`/`D0`/`D1`/`D2`/`GetType`/Line/Circle/
//! Ellipse/Hyperbola/Parabola/NbIntervals/Intervals.
//! The methods that are *not* a plain delegation are implemented here:
//! `Resolution` (`Geom2dAdaptor_Curve.cxx:1186-1219`), `IsClosed`
//! (`Geom2dAdaptor_Curve.cxx:588-600`), the type table
//! (`Geom2dAdaptor_Curve.cxx:96-120`) and the interval grading
//! (`Geom2dAdaptor_Curve.cxx:409-573`).

use super::prelude::*;

use occt_core::kernel::geomabs::CurveType;

use super::analytic_solvers::PI;

/// Unwrap `Geom2d_TrimmedCurve`s to their basis:
/// `Geom2dAdaptor_Curve::load` (`Geom2dAdaptor_Curve.cxx:285-288`) replaces a
/// trimmed curve with its basis and keeps `myFirst`/`myLast` at the trimmed
/// range, so the *type* resolves on the basis while the range stays trimmed.
pub(super) fn adaptor_basis2d(mut c: &dyn Curve2d) -> &dyn Curve2d {
    while let Some(b) = c.trimmed_basis() {
        c = b;
    }
    c
}

/// `Extrema_Curve2dTool::Resolution` (`Extrema_Curve2dTool.hxx:121-124`) ->
/// `Geom2dAdaptor_Curve::Resolution` (`Geom2dAdaptor_Curve.cxx:1186-1219`).
///
/// **UNPORTED**: the `GeomAbs_BezierCurve` (`:1206-1210`) and
/// `GeomAbs_BSplineCurve` (`:1211-1215`) arms call
/// `Geom2d_BezierCurve::Resolution` / `Geom2d_BSplineCurve::Resolution`
/// (`Geom2d_BezierCurve.cxx:679-690`, `Geom2d_BSplineCurve_1.cxx:764-800`,
/// `BSplCLib::Resolution`), which this port does not have; those two types take
/// the `default:` arm `Precision::Parametric(Ruv)` (`:1216-1218`), i.e.
/// `Ruv * Precision::PConfusion() / Precision::Confusion()`.
pub(super) fn adaptor_resolution2d(c: &dyn Curve2d, ruv: f64) -> f64 {
    let base = adaptor_basis2d(c);
    if base.is_line() {
        return ruv; // `:1190-1191`
    }
    if let Some(a_circ) = base.gp_circ2d() {
        let r = a_circ.radius();
        return if r > ruv / 2.0 {
            2.0 * (ruv / (2.0 * r)).asin() // `:1194-1197`
        } else {
            2.0 * PI // `:1199-1201`
        };
    }
    if let Some(an_elips) = base.gp_elips2d() {
        return ruv / an_elips.major_radius; // `:1203-1205`
    }
    ruv * PCONFUSION / CONFUSION // `:1216-1218`
}

/// `Extrema_Curve2dTool::IsClosed` (`Extrema_Curve2dTool.hxx:69`) ->
/// `Geom2dAdaptor_Curve::IsClosed` (`Geom2dAdaptor_Curve.cxx:588-600`): an
/// infinite range is never closed; otherwise the two ends coincide within
/// `Precision::Confusion()`.
pub(super) fn adaptor_is_closed2d(c: &dyn Curve2d) -> bool {
    let (first, last) = (c.first_parameter(), c.last_parameter());
    if Precision::is_positive_infinite(last) || Precision::is_negative_infinite(first) {
        return false; // `:596-599`
    }
    c.d0(first).distance(&c.d0(last)) <= CONFUSION // `:592-594`
}

/// `Geom2dAdaptor_Curve::GetType` (`Geom2dAdaptor_Curve.cxx:96-120`, the
/// `load` type table `:285-346`), resolved on the basis of a trimmed curve.
pub(super) fn adaptor_type2d(c: &dyn Curve2d) -> CurveType {
    let base = adaptor_basis2d(c);
    if base.is_line() {
        CurveType::Line
    } else if base.gp_circ2d().is_some() {
        CurveType::Circle
    } else if base.gp_elips2d().is_some() {
        CurveType::Ellipse
    } else if base.gp_parab2d().is_some() {
        CurveType::Parabola
    } else if base.gp_hypr2d().is_some() {
        CurveType::Hyperbola
    } else if base.offset_basis().is_some() {
        CurveType::OffsetCurve
    } else if base.bezier_nb_poles().is_some() {
        CurveType::BezierCurve
    } else if base.bspline_nb_knots().is_some() || base.bspline_degree().is_some() {
        CurveType::BSplineCurve
    } else {
        CurveType::OtherCurve
    }
}

/// `Geom2dAdaptor_Curve::IsPeriodic` (`Geom2dAdaptor_Curve.cxx:604-607`) =
/// `myCurve->IsPeriodic()` on the **basis** (a trimmed curve is unwrapped by
/// `load`, `:285-288`).
pub(super) fn adaptor_is_periodic2d(c: &dyn Curve2d) -> bool {
    adaptor_basis2d(c).is_periodic()
}

/// `Geom2dAdaptor_Curve::Period` (`Geom2dAdaptor_Curve.cxx:611-614`):
/// `myCurve->LastParameter() - myCurve->FirstParameter()` on the basis.
pub(super) fn adaptor_period2d(c: &dyn Curve2d) -> f64 {
    let base = adaptor_basis2d(c);
    base.last_parameter() - base.first_parameter()
}

/// `Extrema_Curve2dTool::NbIntervals`/`Intervals`
/// (`Extrema_Curve2dTool.hxx:51-62`) -> `Geom2dAdaptor_Curve::{NbIntervals,Intervals}`
/// (`Geom2dAdaptor_Curve.cxx:409-573`).
///
/// The B-spline arm (`:411-451`, `:492-534`) is the port's
/// `Curve2d::parameter_intervals` (`occt_core::bspl::adaptor_intervals`); the
/// offset arm (`:453-480`, `:536-566`) asks a basis adaptor for the *shape
/// below* (`C0->C1, C1->C2, C2->C3, else CN`) and then forces its own ends;
/// every other type is the single span `[First, Last]` (`:482-485`, `:568-572`).
pub(super) fn adaptor_intervals2d(c: &dyn Curve2d, s: u8) -> Vec<f64> {
    let base = adaptor_basis2d(c);
    if let Some(basis) = base.offset_basis() {
        let base_s = match s {
            0 => 2, // GeomAbs_C0 -> GeomAbs_C1 (`:463-465`)
            2 => 4, // GeomAbs_C1 -> GeomAbs_C2 (`:466-468`)
            4 => 5, // GeomAbs_C2 -> GeomAbs_C3 (`:469-471`)
            _ => 6, // GeomAbs_CN (`:472-474`)
        };
        let (lo, hi) = (c.first_parameter(), c.last_parameter());
        let mut iv = adaptor_intervals2d(basis, base_s);
        if iv.len() >= 2 {
            iv[0] = lo; // `:564`
            let n = iv.len();
            iv[n - 1] = hi; // `:565`
        }
        return iv;
    }
    base.parameter_intervals(s)
}