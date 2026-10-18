//! `Geom2dAdaptor_Curve` queries that are not a plain delegation.
//!
//! Sources: `Geom2dAdaptor_Curve.cxx` (the `load` type table `:279-346`,
//! `GetType` `:96-120`). These are the queries an `Adaptor2d_Curve2d` facade is
//! expected to answer on the *stored* curve, so they live on their own rather
//! than inside one consumer (`Extrema_Curve2dTool`, `GeomBndLib_Curve2d` and
//! `BndLib_Add2dCurve` all read the same table).

use occt_core::kernel::geomabs::CurveType;

use crate::curve::Curve2d;

/// Unwrap `Geom2d_TrimmedCurve`s to their basis:
/// `Geom2dAdaptor_Curve::load` (`Geom2dAdaptor_Curve.cxx:285-288`) replaces a
/// trimmed curve with its basis and keeps `myFirst`/`myLast` at the trimmed
/// range, so the *type* resolves on the basis while the range stays trimmed.
pub fn basis2d(mut c: &dyn Curve2d) -> &dyn Curve2d {
    while let Some(b) = c.trimmed_basis() {
        c = b;
    }
    c
}

/// `Geom2dAdaptor_Curve::GetType` (`Geom2dAdaptor_Curve.cxx:96-120`), i.e. the
/// `load` type table (`:279-346`) resolved on the basis of a trimmed curve.
///
/// The type queries are the port's stand-ins for OCCT's
/// `DynamicType() == STANDARD_TYPE(Geom2d_*)` branches; each is answered only
/// by the matching concrete curve, and a `Geom2d_TrimmedCurve` forwards to its
/// basis because `load` recurses into it (`:285-288`).
pub fn type2d(c: &dyn Curve2d) -> CurveType {
    let base = basis2d(c);
    if base.is_line() {
        CurveType::Line // `:292-296`
    } else if base.gp_circ2d().is_some() {
        CurveType::Circle // `:289-292`
    } else if base.gp_elips2d().is_some() {
        CurveType::Ellipse // `:299-303`
    } else if base.gp_parab2d().is_some() {
        CurveType::Parabola // `:304-308`
    } else if base.gp_hypr2d().is_some() {
        CurveType::Hyperbola // `:309-313`
    } else if base.offset_basis().is_some() {
        CurveType::OffsetCurve // `:330-342`
    } else if base.bezier_nb_poles().is_some() {
        CurveType::BezierCurve // `:314-321`
    } else if base.bspline_nb_knots().is_some() || base.bspline_degree().is_some() {
        CurveType::BSplineCurve // `:322-329`
    } else {
        CurveType::OtherCurve // `:343-345`
    }
}
