//! `Blend_CSFunction` — the deferred function for a blending surface between a
//! surface and a curve, driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_CSFunction.hxx:44-199`,
//! `Blend_CSFunction.cxx:18-52` (OCCT 8.0.0).
//!
//! The vector `X` may be `(U, V, W)`: the parametric coordinates on the surface
//! and the parameter on the curve (`hxx:47-53`).
//!
//! Translation notes:
//! - `NbVariables` returns `3` (`cxx:20-23`) and `Pnt1` / `Pnt2` return
//!   `PointOnC()` / `PointOnS()` (`cxx:25-33`) — bodies at this level that
//!   redefine required base methods, so an implementor writes
//!   `fn nb_variables(&self) -> usize { 3 }`,
//!   `fn pnt1(&self) -> &GpPnt { self.point_on_c() }`,
//!   `fn pnt2(&self) -> &GpPnt { self.point_on_s() }` and
//!   `fn section_with_derivatives2(...) -> bool { false }` (`cxx:35-47`). See
//!   [`super::app_function`] for the convention.
//! - `Pnt1` returns the point on the **curve** and `Pnt2` the point on the
//!   **surface** here, the opposite of `Blend_Function`; the pair is kept as
//!   written.
//! - `GetMinimalDistance` (`cxx:49-52`) throws `Standard_NotImplemented`; the
//!   port reaches that through `BlendAppFunction::get_minimal_distance`.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use super::app_function::BlendAppFunction;

/// `class Blend_CSFunction : public Blend_AppFunction` (`hxx:44-199`).
pub trait BlendCsFunction: BlendAppFunction {
    /// `PointOnS()` (`hxx:118`): the point on the surface.
    fn point_on_s(&self) -> &GpPnt;

    /// `PointOnC()` (`hxx:121`): the point on the curve.
    fn point_on_c(&self) -> &GpPnt;

    /// `Pnt2d()` (`hxx:124`): the `U, V` coordinates of the point on the
    /// surface.
    fn pnt2d(&self) -> &GpPnt2d;

    /// `ParameterOnC()` (`hxx:127`): the parameter of the point on the curve.
    fn parameter_on_c(&self) -> f64;

    /// `IsTangencyPoint()` (`hxx:131`).
    fn is_tangency_point(&self) -> bool;

    /// `TangentOnS()` (`hxx:134`).
    fn tangent_on_s(&self) -> &GpVec;

    /// `Tangent2d()` (`hxx:138`): the tangent at `PointOnS`, in the parametric
    /// space of the surface.
    fn tangent2d(&self) -> &GpVec2d;

    /// `TangentOnC()` (`hxx:141`).
    fn tangent_on_c(&self) -> &GpVec;

    /// `Tangent(U, V, TgS, NormS)` (`hxx:146-149`): the tangent at the section
    /// and the normal of the surface at `(U, V)`.
    fn tangent(&self, u: f64, v: f64, tg_s: &mut GpVec, norm_s: &mut GpVec);
}
