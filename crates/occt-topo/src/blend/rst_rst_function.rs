//! `Blend_RstRstFunction` — the deferred function for a blending surface
//! between two restrictions (pcurves), driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_RstRstFunction.hxx:42-233`,
//! `Blend_RstRstFunction.cxx:19-35` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NbVariables` (`hxx:52`) / `NbEquations` (`hxx:55`) are pure; the comment
//!   "Returns 2 (default value). Can be redefined." describes the intent, and
//!   every concrete `BRepBlend_RstRst*` returns `2`.
//! - `Pnt1` / `Pnt2` have bodies here (`cxx:22-30`): `PointOnRst1()` and
//!   `PointOnRst2()`. They redefine required base methods, so an implementor
//!   writes those two one-liners. See [`super::app_function`] for the
//!   convention.
//! - `GetMinimalDistance` (`cxx:32-34`) throws `Standard_NotImplemented`; the
//!   port reaches that through `BlendAppFunction::get_minimal_distance`.
//! - `Decroch` returns a [`BlendDecrochStatus`] here, unlike the `bool` of
//!   `Blend_SurfRstFunction` (`hxx:154-157` vs `hxx:157`).

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::MathVector;

use super::app_function::BlendAppFunction;
use super::status::BlendDecrochStatus;

/// `class Blend_RstRstFunction : public Blend_AppFunction` (`hxx:42-233`).
pub trait BlendRstRstFunction: BlendAppFunction {
    /// `PointOnRst1()` (`hxx:114`).
    fn point_on_rst1(&self) -> &GpPnt;

    /// `PointOnRst2()` (`hxx:117`).
    fn point_on_rst2(&self) -> &GpPnt;

    /// `Pnt2dOnRst1()` (`hxx:120`).
    fn pnt2d_on_rst1(&self) -> &GpPnt2d;

    /// `Pnt2dOnRst2()` (`hxx:124`).
    fn pnt2d_on_rst2(&self) -> &GpPnt2d;

    /// `ParameterOnRst1()` (`hxx:127`).
    fn parameter_on_rst1(&self) -> f64;

    /// `ParameterOnRst2()` (`hxx:130`).
    fn parameter_on_rst2(&self) -> f64;

    /// `IsTangencyPoint()` (`hxx:134`).
    fn is_tangency_point(&self) -> bool;

    /// `TangentOnRst1()` (`hxx:137`).
    fn tangent_on_rst1(&self) -> &GpVec;

    /// `Tangent2dOnRst1()` (`hxx:141`).
    fn tangent2d_on_rst1(&self) -> &GpVec2d;

    /// `TangentOnRst2()` (`hxx:144`).
    fn tangent_on_rst2(&self) -> &GpVec;

    /// `Tangent2dOnRst2()` (`hxx:148`).
    fn tangent2d_on_rst2(&self) -> &GpVec2d;

    /// `Decroch(Sol, NRst1, TgRst1, NRst2, TgRst2)` (`hxx:154-157`): the
    /// function-specific decrochage criterion. OCCT warns it may be called
    /// without a previous `IsSolution`, in which case the computed values can
    /// be meaningless.
    fn decroch(
        &self,
        sol: &MathVector,
        n_rst1: &mut GpVec,
        tg_rst1: &mut GpVec,
        n_rst2: &mut GpVec,
        tg_rst2: &mut GpVec,
    ) -> BlendDecrochStatus;
}
