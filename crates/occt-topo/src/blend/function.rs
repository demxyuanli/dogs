//! `Blend_Function` — the deferred function for a blending surface between two
//! surfaces, driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_Function.hxx:44-121`,
//! `Blend_Function.cxx:18-58` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Values` / `Derivatives` is `(U1, V1, U2, V2)`:
//! the parametric coordinates of the two extremities of a section on the first
//! and second surface (`hxx:37-42`).
//!
//! Translation notes:
//! - `TwistOnS1` / `TwistOnS2` (`cxx:32-42`) and the 10-array `Section`
//!   (`cxx:44-58`) are bodies written at this level and are defaults here.
//!   `NbSection` / `Pnt1` / `Pnt2` are also bodies at this level in OCCT
//!   (`cxx:19-32`) but they redefine required base methods, which Rust
//!   sub-traits cannot supply; an implementor writes
//!   `fn nb_variables(&self) -> usize { 4 }`,
//!   `fn pnt1(&self) -> &GpPnt { self.point_on_s1() }`,
//!   `fn pnt2(&self) -> &GpPnt { self.point_on_s2() }`,
//!   `fn section_with_derivatives2(...) -> bool { false }` — exactly those
//!   bodies. See [`super::app_function`] for the port-wide convention.
//! - `Tangent(...)` (`hxx:87-96`) returns four vectors through references; the
//!   port keeps the four `&mut` outputs.

use occt_core::gp::{GpPnt, GpVec, GpVec2d};

use super::app_function::BlendAppFunction;

/// `class Blend_Function : public Blend_AppFunction` (`hxx:44-121`).
pub trait BlendFunction: BlendAppFunction {
    /// `PointOnS1()` (`hxx:59`): the point on the first surface, at
    /// `Sol(1), Sol(2)`.
    fn point_on_s1(&self) -> &GpPnt;

    /// `PointOnS2()` (`hxx:63`): the point on the second surface, at
    /// `Sol(3), Sol(4)`.
    fn point_on_s2(&self) -> &GpPnt;

    /// `IsTangencyPoint()` (`hxx:68`): `true` when the tangent vectors at
    /// `PointOnS1` and/or `PointOnS2` cannot be computed.
    fn is_tangency_point(&self) -> bool;

    /// `TangentOnS1()` (`hxx:71`).
    fn tangent_on_s1(&self) -> &GpVec;

    /// `Tangent2dOnS1()` (`hxx:75`).
    fn tangent2d_on_s1(&self) -> &GpVec2d;

    /// `TangentOnS2()` (`hxx:78`).
    fn tangent_on_s2(&self) -> &GpVec;

    /// `Tangent2dOnS2()` (`hxx:82`).
    fn tangent2d_on_s2(&self) -> &GpVec2d;

    /// `Tangent(U1, V1, U2, V2, TgFirst, TgLast, NormFirst, NormLast)`
    /// (`hxx:87-96`): the tangents at the beginning and the end of the section,
    /// and the normals of the surfaces at those points.
    #[allow(clippy::too_many_arguments)]
    fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut GpVec,
        tg_last: &mut GpVec,
        norm_first: &mut GpVec,
        norm_last: &mut GpVec,
    );

    /// `TwistOnS1()` (`cxx:34-37`): `false`.
    fn twist_on_s1(&self) -> bool {
        false
    }

    /// `TwistOnS2()` (`cxx:39-42`): `false`.
    fn twist_on_s2(&self) -> bool {
        false
    }
}
