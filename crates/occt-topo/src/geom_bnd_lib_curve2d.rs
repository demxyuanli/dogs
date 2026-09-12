//! Dispatcher for 2D curve bounding boxes.
//!
//! Source: `GeomBndLib_Curve2d.cxx` (`Box(U1,U2,Tol)` / `Add`). The Geom2d
//! handle constructor downcasts Line / Circle / Ellipse / Hyperbola /
//! Parabola / Bezier / BSpline / Offset; everything else is
//! `GeomBndLib_OtherCurve2d`. This port classifies via [`pc_curve_kind`]
//! because `Arc<dyn Curve2d>` has no `down_cast`. Ellipse/hyperbola/parabola
//! fall through to OtherCurve2d sampling unless the caller already holds the
//! concrete `gp_*2d` (the analytic `box_*` functions remain available).
//!
//! `BoxOptimal` is not dispatched: PerformAreas uses `Add(..., 0.)` → `Box()`.

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d};
use occt_geom2d::curve::Curve2d;

use crate::geom_bnd_lib_bspline2d::box_bspline_as_curve;
use crate::geom_bnd_lib_circle2d::box_circ_range;
use crate::geom_bnd_lib_ellipse2d::box_elips_range;
use crate::geom_bnd_lib_hyperbola2d::box_hypr;
use crate::geom_bnd_lib_line2d::box_lin;
use crate::geom_bnd_lib_other2d::box_other;
use crate::geom_bnd_lib_parabola2d::box_parab;
use crate::pcurve::{pc_curve_kind, CurveKind};

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
    match pc_curve_kind(curve) {
        CurveKind::Line => GeomAbsCurve2d::Line,
        CurveKind::Circle => GeomAbsCurve2d::Circle,
        CurveKind::BSpline => GeomAbsCurve2d::BSpline,
        CurveKind::Other => GeomAbsCurve2d::Other,
    }
}

/// Reconstruct a `gp_Lin2d` from an unbounded 2D curve (`d0(0)`, `d1(0)`).
pub fn lin2d_from_curve(curve: &dyn Curve2d) -> Option<GpLin2d> {
    let p = curve.d0(0.0);
    let (_, t) = curve.d1(0.0);
    let d = occt_core::gp::GpDir2d::from_vec2d(&t).ok()?;
    Some(GpLin2d::from_pnt_dir(p, d))
}

/// Reconstruct a `gp_Circ2d` from a 2-pi periodic curve (midpoint of `d0(0)`
/// and `d0(pi)` is the centre for the standard parametrisation).
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
pub fn box_curve2d(curve: &dyn Curve2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    match curve2d_abs_type(curve) {
        GeomAbsCurve2d::Line => {
            if let Some(lin) = lin2d_from_curve(curve) {
                if let Ok(b) = box_lin(&lin, the_u1, the_u2, the_tol) {
                    return b;
                }
            }
            box_other(curve, the_u1, the_u2, the_tol)
        }
        GeomAbsCurve2d::Circle => {
            if let Some(c) = circ2d_from_curve(curve) {
                return box_circ_range(&c, the_u1, the_u2, the_tol);
            }
            box_other(curve, the_u1, the_u2, the_tol)
        }
        GeomAbsCurve2d::BSpline => box_bspline_as_curve(curve, the_u1, the_u2, the_tol),
        _ => box_other(curve, the_u1, the_u2, the_tol),
    }
}

/// `GeomBndLib_Curve2d::Add(theU1, theU2, theTol, theBox)`.
pub fn add_curve2d(curve: &dyn Curve2d, the_u1: f64, the_u2: f64, the_tol: f64, the_box: &mut BndBox2d) {
    the_box.add_box(&box_curve2d(curve, the_u1, the_u2, the_tol));
}

/// Direct analytic boxes when the caller already has the `gp_*` (the static
/// overloads of `GeomBndLib_*2d::Box`).
pub fn box_lin_gp(lin: &GpLin2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    box_lin(lin, u1, u2, tol).unwrap_or_else(|_| BndBox2d::new())
}

pub fn box_circ_gp(c: &GpCirc2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    box_circ_range(c, u1, u2, tol)
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
