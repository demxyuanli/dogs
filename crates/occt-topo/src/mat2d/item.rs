//! The `Geom2d_Geometry` item a `MAT2d` contour is made of: a
//! `Geom2d_TrimmedCurve` or a `Geom2d_CartesianPoint`.
//!
//! `MAT2d_Circuit` / `MAT2d_MiniPath` take their input as
//! `NCollection_Sequence<NCollection_Sequence<occ::handle<Geom2d_Geometry>>>`
//! (`MAT2d_Circuit.hxx:44-48`, `MAT2d_MiniPath.hxx:53-57`) and every curve of
//! that figure is down-cast to `Geom2d_TrimmedCurve`
//! (`MAT2d_Circuit.cxx:103-105`, `:236-238`), so the two variants below are the
//! whole of the geometry the port has to carry.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

/// `occ::handle(Geom2d_Geometry)` restricted to the circuit items.
#[derive(Clone)]
pub enum Mat2dItem {
    /// `Geom2d_TrimmedCurve`.
    Curve(Geom2dTrimmedCurve),
    /// `Geom2d_CartesianPoint`.
    Point(GpPnt2d),
}

impl std::fmt::Debug for Mat2dItem {
    /// `DynamicType()` plus the trim range, for the diagnostics that OCCT
    /// prints under `OCCT_DEBUG`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Curve(c) => write!(
                f,
                "Mat2dItem::Curve[{:?}, {:?}]",
                c.first_parameter(),
                c.last_parameter()
            ),
            Self::Point(p) => write!(f, "Mat2dItem::Point[{:?}, {:?}]", p.x(), p.y()),
        }
    }
}

impl Mat2dItem {
    /// `DynamicType() == STANDARD_TYPE(Geom2d_CartesianPoint)`.
    pub fn is_point(&self) -> bool {
        matches!(self, Self::Point(_))
    }

    /// `occ::down_cast<Geom2d_Point>(...)->Pnt2d()`.
    pub fn point(&self) -> Option<GpPnt2d> {
        match self {
            Self::Point(p) => Some(*p),
            Self::Curve(_) => None,
        }
    }

    /// `occ::down_cast<Geom2d_Curve>(...)`.
    pub fn curve(&self) -> Option<&Geom2dTrimmedCurve> {
        match self {
            Self::Curve(c) => Some(c),
            Self::Point(_) => None,
        }
    }

    /// The `Geom2d_TrimmedCurve` of an item that OCCT assumes is a curve; panics
    /// where the OCCT down-cast would yield a null handle.
    pub fn expect_curve(&self) -> &Geom2dTrimmedCurve {
        self.curve()
            .unwrap_or_else(|| panic!("MAT2d item is not a Geom2d_TrimmedCurve"))
    }

    /// `Curve->StartPoint()` (`Geom2d_TrimmedCurve::StartPoint`,
    /// `Geom2d_TrimmedCurve.cxx:63-66`): `Value(FirstParameter())`.
    pub fn start_point(&self) -> GpPnt2d {
        let c = self.expect_curve();
        c.d0(c.first_parameter())
    }

    /// `Curve->EndPoint()` (`Geom2d_TrimmedCurve::EndPoint`, `cxx:68-71`).
    pub fn end_point(&self) -> GpPnt2d {
        let c = self.expect_curve();
        c.d0(c.last_parameter())
    }

    /// `Geom2d_Geometry::Copy()`: a point is immutable and shared, a curve is
    /// duplicated so that `Reverse()` on the copy leaves the item untouched
    /// (`MAT2d_Circuit.cxx:533-536`).
    pub fn copy_item(&self) -> Self {
        match self {
            Self::Curve(c) => Self::Curve(c.clone()),
            Self::Point(p) => Self::Point(*p),
        }
    }

    /// `Geom2d_Curve::DN(U, N)` on a curve item.
    pub fn dn(&self, u: f64, n: i32) -> Option<GpVec2d> {
        self.curve().map(|c| c.eval_dn(u, n))
    }

    /// `Geom2d_Curve::ReversedParameter(U)` on a curve item.
    pub fn reversed_parameter(&self, u: f64) -> Option<f64> {
        self.curve().map(|c| c.reversed_parameter(u))
    }
}

/// `NCollection_Sequence<occ::handle(Geom2d_Geometry)>`: one line of the figure.
pub type Mat2dLine = Vec<Mat2dItem>;

/// `NCollection_Sequence<NCollection_Sequence<occ::handle(Geom2d_Geometry)>>`:
/// the figure, one entry per contour line.
pub type Mat2dFigure = Vec<Mat2dLine>;
