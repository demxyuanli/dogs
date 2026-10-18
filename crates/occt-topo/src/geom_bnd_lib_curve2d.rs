//! Dispatcher for 2D curve bounding boxes.
//!
//! Source: `GeomBndLib_Curve2d.cxx` (`Box(U1,U2,Tol)` / `Add`) and the
//! `BndLib_Add2dCurve::Add(Geom2d_Curve, ...)` entry point (`:76-83`) that
//! builds it. The handle constructor (`:97-171`) classifies the curve with the
//! `Geom2d_Curve` downcast chain Line / Circle / Ellipse / Hyperbola /
//! Parabola / Bezier / BSpline / Offset — unwrapping a `Geom2d_TrimmedCurve`
//! to its `BasisCurve()` first (`:150-157`) — and everything else becomes
//! `GeomBndLib_OtherCurve2d`.
//!
//! This port has no `down_cast`, so [`curve2d_abs_type`] reads the equivalent
//! `Geom2dAdaptor_Curve::GetType` table (`occt_geom2d::adaptor_curve::
//! get_type2d`) and each arm rebuilds the concrete geometry from the `Curve2d`
//! trait queries (`gp_lin2d` / `gp_circ2d` / `gp_elips2d` / `gp_hypr2d` /
//! `gp_parab2d` / `poles2d` / `bspline_copy2d` / `offset_basis` +
//! `offset_value`). Line and Circle keep a sample-based reconstruction
//! (`lin2d_from_curve` / `circ2d_from_curve`) as the stand-in for the
//! unported `Adaptor2d_Curve2d::Line()/Circle()` of the port's own wrapper
//! classes (`pcurve::PlaneKeepParam2d` deliberately does not forward
//! `gp_circ2d`).
//!
//! `BoxOptimal` is not dispatched: PerformAreas uses `Add(..., 0.)` → `Box()`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d};
use occt_core::kernel::geomabs::CurveType;
use occt_geom2d::adaptor_curve::type2d;
use occt_geom2d::curve::Curve2d;

use crate::geom_bnd_lib_bezier2d::box_bezier_as_curve;
use crate::geom_bnd_lib_bspline2d::box_bspline_as_curve;
use crate::geom_bnd_lib_circle2d::box_circ_range_2d;
use crate::geom_bnd_lib_ellipse2d::box_elips_range;
use crate::geom_bnd_lib_hyperbola2d::box_hypr;
use crate::geom_bnd_lib_line2d::box_lin_2d;
use crate::geom_bnd_lib_offset2d::box_offset_as_curve;
use crate::geom_bnd_lib_other_curve2d::box_other_2d;
use crate::geom_bnd_lib_parabola2d::box_parab;

/// Analytic 2D curve type matching `GeomAbs_CurveType` as used by the dispatcher.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsCurve2d {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    Bezier,
    BSpline,
    Offset,
    Other,
}

/// Map a pcurve to the dispatcher kind (`GeomBndLib_Curve2d` ctor from handle).
pub fn curve2d_abs_type(curve: &dyn Curve2d) -> GeomAbsCurve2d {
    match type2d(curve) {
        CurveType::Line => GeomAbsCurve2d::Line,
        CurveType::Circle => GeomAbsCurve2d::Circle,
        CurveType::Ellipse => GeomAbsCurve2d::Ellipse,
        CurveType::Hyperbola => GeomAbsCurve2d::Hyperbola,
        CurveType::Parabola => GeomAbsCurve2d::Parabola,
        CurveType::BezierCurve => GeomAbsCurve2d::Bezier,
        CurveType::BSplineCurve => GeomAbsCurve2d::BSpline,
        CurveType::OffsetCurve => GeomAbsCurve2d::Offset,
        CurveType::OtherCurve => GeomAbsCurve2d::Other,
    }
}

/// Reconstruct a `gp_Lin2d` from an unbounded 2D curve (`d0(0)`, `d1(0)`).
///
/// Stand-in for `Adaptor2d_Curve2d::Line()` on the port's wrapper classes that
/// report `is_line()` without forwarding `gp_lin2d` (`pcurve::PlaneKeepParam2d`
/// over a projected line, whose `gp_Lin2d` this port does not build).
pub fn lin2d_from_curve(curve: &dyn Curve2d) -> Option<GpLin2d> {
    let p = curve.d0(0.0);
    let (_, t) = curve.d1(0.0);
    let d = occt_core::gp::GpDir2d::from_vec2d(&t).ok()?;
    Some(GpLin2d::from_pnt_dir(p, d))
}

/// Reconstruct a `gp_Circ2d` from a 2-pi periodic curve (midpoint of `d0(0)`
/// and `d0(pi)` is the centre for the standard parametrisation).
///
/// Stand-in for `Adaptor2d_Curve2d::Circle()` on the port's wrapper classes
/// that do not forward `gp_circ2d` (see [`lin2d_from_curve`]).
pub fn circ2d_from_curve(curve: &dyn Curve2d) -> Option<GpCirc2d> {
    let p0 = curve.d0(0.0);
    let ppi = curve.d0(std::f64::consts::PI);
    let centre = occt_core::gp::GpPnt2d::new(0.5 * (p0.x() + ppi.x()), 0.5 * (p0.y() + ppi.y()));
    let r = centre.distance(&p0);
    if r <= occt_core::precision::CONFUSION {
        return None;
    }
    let mut c = GpCirc2d::default();
    c.set_location(centre);
    c.set_radius(r);
    Some(c)
}

/// `GeomBndLib_Curve2d::Box(theU1, theU2, theTol)` for a `Curve2d` handle.
pub fn box_curve2d(curve: &dyn Curve2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    match curve2d_abs_type(curve) {
        GeomAbsCurve2d::Line => {
            if let Some(lin) = curve.gp_lin2d().or_else(|| lin2d_from_curve(curve)) {
                if let Ok(b) = box_lin_2d(&lin, u1, u2, tol) {
                    return b;
                }
            }
            box_other_2d(curve, u1, u2, tol)
        }
        GeomAbsCurve2d::Circle => {
            if let Some(c) = curve.gp_circ2d().or_else(|| circ2d_from_curve(curve)) {
                return box_circ_range_2d(&c, u1, u2, tol);
            }
            box_other_2d(curve, u1, u2, tol)
        }
        GeomAbsCurve2d::Ellipse => {
            if let Some(e) = curve.gp_elips2d() {
                return box_elips_range(&e, u1, u2, tol);
            }
            box_other_2d(curve, u1, u2, tol)
        }
        GeomAbsCurve2d::Hyperbola => {
            if let Some(h) = curve.gp_hypr2d() {
                if let Ok(b) = box_hypr(&h, u1, u2, tol) {
                    return b;
                }
            }
            box_other_2d(curve, u1, u2, tol)
        }
        GeomAbsCurve2d::Parabola => {
            if let Some(p) = curve.gp_parab2d() {
                if let Ok(b) = box_parab(&p, u1, u2, tol) {
                    return b;
                }
            }
            box_other_2d(curve, u1, u2, tol)
        }
        GeomAbsCurve2d::Bezier => box_bezier_as_curve(curve, u1, u2, tol),
        GeomAbsCurve2d::BSpline => box_bspline_as_curve(curve, u1, u2, tol),
        GeomAbsCurve2d::Offset => box_offset_as_curve(curve, u1, u2, tol),
        GeomAbsCurve2d::Other => box_other_2d(curve, u1, u2, tol),
    }
}

/// `GeomBndLib_Curve2d::Add(theU1, theU2, theTol, theBox)`.
pub fn add_curve2d(curve: &dyn Curve2d, u1: f64, u2: f64, tol: f64, bnd_box: &mut BndBox2d) {
    bnd_box.add_box(&box_curve2d(curve, u1, u2, tol));
}

/// Direct analytic boxes when the caller already has the `gp_*` (the static
/// overloads of `GeomBndLib_*2d::Box`).
pub fn box_lin_gp(lin: &GpLin2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    box_lin_2d(lin, u1, u2, tol).unwrap_or_else(|_| BndBox2d::new())
}

pub fn box_circ_gp(c: &GpCirc2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    box_circ_range_2d(c, u1, u2, tol)
}

pub fn box_elips_gp(e: &GpElips2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    box_elips_range(e, u1, u2, tol)
}

pub fn box_hypr_gp(h: &GpHypr2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    box_hypr(h, u1, u2, tol).unwrap_or_else(|_| BndBox2d::new())
}

pub fn box_parab_gp(p: &GpParab2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    box_parab(p, u1, u2, tol).unwrap_or_else(|_| BndBox2d::new())
}
