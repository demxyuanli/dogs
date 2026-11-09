//! `Contap_IType` and `Contap_TFunction` - the two enumerations that classify
//! a contour line and the function that produced it.
//!
//! Sources: `Contap_IType.hxx:22-27`, `Contap_TFunction.hxx:22-28`.

/// `enum Contap_IType` (`Contap_IType.hxx:22-27`): the kind of contour line a
/// [`super::ContapLine`] carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash)]
pub enum ContapIType {
    /// `Contap_Lin`: a straight line.
    Lin,
    /// `Contap_Circle`: a circle.
    Circle,
    /// `Contap_Walking`: a walked intersection line stored in the
    /// `IntSurf_LineOn2S` buffer.
    #[default]
    Walking,
    /// `Contap_Restriction`: a piece of the surface's restriction.
    Restriction,
}

/// `enum Contap_TFunction` (`Contap_TFunction.hxx:22-28`): the form of the
/// contour function `F(u, v) = 0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash)]
pub enum ContapTFunction {
    /// `Contap_ContourStd`: `N . Dir` (parallel projection).
    #[default]
    ContourStd,
    /// `Contap_ContourPrs`: `N . (Eye - P)` (perspective).
    ContourPrs,
    /// `Contap_DraftStd`: `N . Dir - cosAng * |N|` (draft, parallel).
    DraftStd,
    /// `Contap_DraftPrs`: draft with a perspective eye. Both functions leave
    /// the value untouched in this arm (`Contap_SurfFunction.cxx:102-105`,
    /// `Contap_ArcFunction.cxx:79-82`); it is only produced by the `Set(Eye,
    /// Angle)` overloads.
    DraftPrs,
}
