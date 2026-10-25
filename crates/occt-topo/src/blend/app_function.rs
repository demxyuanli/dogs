//! `Blend_AppFunction` — the deferred function an approximation driver
//! (`AppBlend_Approx` / `Blend_AppSurface`) consumes to build one blending
//! surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_AppFunction.hxx:39-192`,
//! `Blend_AppFunction.cxx:20-23` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The `math_FunctionSetWithDerivatives` base (`hxx:44`, pure virtuals at
//!   `:50-73`) is the crate trait of the same name, so `NbVariables` /
//!   `NbEquations` / `Value` / `Derivatives` / `Values` are *not* declared
//!   again here; `Blend_AppFunction` keeps only what that base does not have.
//! - C++ overloading: the two `Set` (`hxx:78`, `:84`), the two `GetTolerance`
//!   (`:87`, `:145`) and the three `Section` (`:158`, `:166`, `:173`) overloads
//!   each get one Rust name — `set` / `set_interval`, `get_tolerance` /
//!   `get_tolerance_3d`, `section` / `section_poles` /
//!   `section_with_derivatives2`.
//! - `NCollection_Array1<T>&` out-parameters become `&mut [T]`; OCCT's
//!   `Standard_OutOfRange` for a too-short array is the caller's slice length
//!   here, which is checked by the callee when it indexes.
//! - `Parameter(P)` (`cxx:20-23`) is the only method with a body at this level
//!   and keeps it as a default.
//!
//! Port-wide convention for the `Blend` interface family: Rust subtraits have
//! no inherited method bodies, so the bodies OCCT writes in an intermediate
//! class (`Blend_Function::NbVariables` returning `4`, `Blend_SurfRstFunction::Pnt1`
//! returning `PointOnS()`, ...) cannot be attached to the sub-trait. Every
//! implementor writes them explicitly; each sub-trait's doc comment names the
//! OCCT `.cxx` line that supplies the required body.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_math::{MathFunctionSetWithDerivatives, MathVector};

use super::point::BlendPoint;

/// `class Blend_AppFunction : public math_FunctionSetWithDerivatives`
/// (`hxx:39-192`).
pub trait BlendAppFunction: MathFunctionSetWithDerivatives {
    /// `Set(Param)` (`hxx:78`): the parameter along the guide line, which
    /// selects the plane the solution is searched in.
    fn set(&mut self, param: f64);

    /// `Set(First, Last)` (`hxx:84`): the bounds of the parametric interval on
    /// the guide line.
    fn set_interval(&mut self, first: f64, last: f64);

    /// `GetTolerance(Tolerance, Tol)` (`hxx:87`): the parametric tolerance of
    /// each variable, for the 3d tolerance `tol`.
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64);

    /// `GetBounds(InfBound, SupBound)` (`hxx:92`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector);

    /// `IsSolution(Sol, Tol)` (`hxx:99`): `true` when `sol` is a zero of the
    /// function at the current guide parameter.
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool;

    /// `GetMinimalDistance()` (`hxx:103`): the minimal distance between the two
    /// extremities of the computed sections.
    ///
    /// `Blend_SurfRstFunction`, `Blend_RstRstFunction` and `Blend_CSFunction`
    /// override this with `throw Standard_NotImplemented` (`cxx:32-34`,
    /// `:49-51`); the port carries those as `Err`, so a leaf that does not
    /// compute a distance returns the same failure OCCT raises.
    fn get_minimal_distance(&self) -> Result<f64, String>;

    /// `Pnt1()` (`hxx:106`): the point on the first support.
    fn pnt1(&self) -> &GpPnt;

    /// `Pnt2()` (`hxx:110`): the point on the second support.
    fn pnt2(&self) -> &GpPnt;

    /// `IsRational()` (`hxx:114`).
    fn is_rational(&self) -> bool;

    /// `GetSectionSize()` (`hxx:117`): the length of the maximum section.
    fn get_section_size(&self) -> f64;

    /// `GetMinimalWeight(Weigths)` (`hxx:119`).
    fn get_minimal_weight(&self, weights: &mut [f64]);

    /// `NbIntervals(S)` (`hxx:124`): the number of continuity intervals for
    /// `s`; may be one when the continuity is at least `s`.
    fn nb_intervals(&self, s: Shape) -> usize;

    /// `Intervals(T, S)` (`hxx:131`): the parameters bounding the continuity
    /// intervals for `s`.
    fn intervals(&self, t: &mut [f64], s: Shape);

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`hxx:135`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    );

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1D)` (`hxx:145`):
    /// the tolerance to reach in approximation.
    #[allow(clippy::too_many_arguments)]
    fn get_tolerance_3d(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut MathVector,
        tol1d: &mut MathVector,
    );

    /// `Knots(TKnots)` (`hxx:152`).
    fn knots(&mut self, tknots: &mut [f64]);

    /// `Mults(TMults)` (`hxx:154`).
    fn mults(&mut self, tmults: &mut [i32]);

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`
    /// (`hxx:158`): the first/last section, with first derivatives. Returns
    /// `true` when the derivatives were computed.
    #[allow(clippy::too_many_arguments)]
    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        dpoles: &mut [GpVec],
        poles2d: &mut [GpPnt2d],
        dpoles2d: &mut [GpVec2d],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool;

    /// `Section(P, Poles, Poles2d, Weigths)` (`hxx:166`): the section without
    /// derivatives.
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    );

    /// `Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weigths, DWeigths, D2Weigths)` (`hxx:173`): the first/last section, with
    /// first and second derivatives.
    ///
    /// This is required (OCCT `= 0`) at this level; `Blend_Function` and
    /// `Blend_CSFunction` override it with `return false` (`Blend_Function.cxx:44-58`,
    /// `Blend_CSFunction.cxx:25-41`), so their implementors answer `false`.
    #[allow(clippy::too_many_arguments)]
    fn section_with_derivatives2(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        dpoles: &mut [GpVec],
        d2poles: &mut [GpVec],
        poles2d: &mut [GpPnt2d],
        dpoles2d: &mut [GpVec2d],
        d2poles2d: &mut [GpVec2d],
        weights: &mut [f64],
        dweights: &mut [f64],
        d2weights: &mut [f64],
    ) -> bool;

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`hxx:186`): the parametric
    /// resolution on the `IC2d`-th 2d curve.
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64);

    /// `Parameter(P)` (`cxx:20-23`): the parameter of `p`, used to impose the
    /// parameters in the approximation.
    fn parameter(&self, p: &BlendPoint) -> f64 {
        p.parameter()
    }
}
