//! `Blend_SurfRstFunction` — the deferred function for a blending surface
//! between a surface and a pcurve on another surface, driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_SurfRstFunction.hxx:38-223`,
//! `Blend_SurfRstFunction.cxx:19-35` (OCCT 8.0.0).
//!
//! The vector `X` may be `(U, V, W)`: the parametric coordinates of the section
//! extremity on the surface, then on the curve (`hxx:41-48`).
//!
//! Translation notes:
//! - `NbVariables` (`hxx:52`) and `NbEquations` (`hxx:55`) are pure here; the
//!   comment "Returns 3 (default value). Can be redefined." describes the
//!   intent, not a body — every concrete `BRepBlend_SurfRst*` writes its `3`.
//! - `Pnt1` / `Pnt2` have bodies at this level (`cxx:22-30`): `PointOnS()` and
//!   `PointOnRst()`. They redefine required base methods, so an implementor
//!   writes `fn pnt1(&self) -> &GpPnt { self.point_on_s() }` and
//!   `fn pnt2(&self) -> &GpPnt { self.point_on_rst() }`. See
//!   [`super::app_function`] for the convention.
//! - `GetMinimalDistance` (`cxx:32-34`) throws `Standard_NotImplemented`; the
//!   port reaches that through `BlendAppFunction::get_minimal_distance` and the
//!   implementor answers `Err(..)` with the same message.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use super::app_function::BlendAppFunction;

/// `class Blend_SurfRstFunction : public Blend_AppFunction` (`hxx:38-223`).
pub trait BlendSurfRstFunction: BlendAppFunction {
    /// `PointOnS()` (`hxx:120`): the point on the surface.
    fn point_on_s(&self) -> &GpPnt;

    /// `PointOnRst()` (`hxx:123`): the point on the curve.
    fn point_on_rst(&self) -> &GpPnt;

    /// `Pnt2dOnS()` (`hxx:126`): the `U, V` coordinates of the point on the
    /// surface.
    fn pnt2d_on_s(&self) -> &GpPnt2d;

    /// `Pnt2dOnRst()` (`hxx:130`): the `U, V` coordinates of the point on the
    /// curve on surface.
    fn pnt2d_on_rst(&self) -> &GpPnt2d;

    /// `ParameterOnRst()` (`hxx:133`): the parameter of the point on the curve.
    fn parameter_on_rst(&self) -> f64;

    /// `IsTangencyPoint()` (`hxx:137`).
    fn is_tangency_point(&self) -> bool;

    /// `TangentOnS()` (`hxx:140`).
    fn tangent_on_s(&self) -> &GpVec;

    /// `Tangent2dOnS()` (`hxx:144`): the tangent at `PointOnS`, in the
    /// parametric space of the surface.
    fn tangent2d_on_s(&self) -> &GpVec2d;

    /// `TangentOnRst()` (`hxx:147`).
    fn tangent_on_rst(&self) -> &GpVec;

    /// `Tangent2dOnRst()` (`hxx:151`): the tangent at `PointOnRst`, in the
    /// parametric space of the second surface.
    fn tangent2d_on_rst(&self) -> &GpVec2d;

    /// `Decroch(Sol, NS, TgS)` (`hxx:157`): a function-specific "leaving the
    /// restriction" criterion.
    fn decroch(&self, sol: &occt_math::MathVector, ns: &mut GpVec, tg_s: &mut GpVec) -> bool;
}
