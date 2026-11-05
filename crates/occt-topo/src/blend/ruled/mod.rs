//! `BlendFunc_Ruled` — the ruled blending surface between two surfaces: the
//! section is the straight segment joining a point of the first surface to a
//! point of the second, both lying in the section plane of the guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_Ruled.hxx:34-207`,
//! `BlendFunc_Ruled.cxx:31-744` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(U1, V1, U2, V2)`.
//! With `nplan` the unit tangent of the guide line at the current parameter and
//! `theD = -(nplan . ptgui)`, the four equations are
//! 1. `nplan . pts1 + theD` — the first extremity lies in the section plane,
//! 2. `nplan . pts2 + theD` — the second extremity lies in the same plane,
//! 3. `(pts2 - pts1) . nor1` — the segment is tangent to the first surface,
//! 4. `(pts2 - pts1) . nor2` — the segment is tangent to the second surface,
//!
//! where `nor1` / `nor2` are the surface normals projected into the section
//! plane (`ns1 - (nplan . ns1) nplan`, normalised).
//!
//! `BlendFunc_Ruled` is the degenerate member of the `BlendFunc_*` family: it
//! has no radius, so
//! - no `Eps` guard around the projected-normal norms (`cxx:189-194`,
//!   `:307-318`): a singular projection divides by zero exactly as OCCT does,
//! - `IsSolution` tests all four residuals against `Tol` one by one
//!   (`cxx:106-109`) instead of the `|f1|` plus squared-norm test of
//!   `ConstRad` / `EvolRad`,
//! - `Set(First, Last)` raises `Standard_NotImplemented` (`cxx:60-63`); the
//!   port panics with that message,
//! - `GetSectionSize` raises `Standard_NotImplemented` (`cxx:580-583`),
//! - `distmin` is initialised to `RealLast()` and never written (`cxx:42`,
//!   `:175`), so `GetMinimalDistance` answers it unchanged.
//!
//! Translation notes:
//! - `GetSection` (`cxx:469-568`) and `AxeRot` (`cxx:713-730`) are public in
//!   OCCT but have no caller anywhere in the OCCT tree; they are ported as
//!   inherent methods so the class is complete.
//! - `math_Gauss Resol(gradsol)` (`cxx:149`, `:543`) uses the default
//!   `MinPivot = 1.0e-20` (`math_Gauss.hxx:49`); the port passes it
//!   explicitly.
//! - `TangentOnS1` / `Tangent2dOnS1` / `TangentOnS2` / `Tangent2dOnS2` raise
//!   `Standard_DomainError` when `istangent` (`cxx:432-466`); the trait
//!   returns a value, and every OCCT caller guards the call with
//!   `!IsTangencyPoint()`, so the stored value is returned.
//! - `TwistOnS1` / `TwistOnS2` are not overridden by this class
//!   (`hxx:34-207`), so the `Blend_Function` bodies apply: both answer
//!   `false` (`Blend_Function.cxx:32-42`).

mod equations;
mod ruled_inv;
mod section;
mod solution;

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{BlendAppFunction, BlendFunction, BlendPoint};

pub use ruled_inv::RuledInv;

/// `math_Gauss`'s default `MinPivot` (`math_Gauss.hxx:49`), the pivot threshold
/// of the two linear systems of this class (`cxx:149`, `:543`).
pub(super) const GAUSS_MIN_PIVOT: f64 = 1.0e-20;

/// `p.XYZ()` as a `gp_Vec` — `nplan.XYZ().Dot(ptgui.XYZ())` (`cxx:57`) and the
/// mid-point expressions of the section files.
pub(super) fn pnt_vec(p: &GpPnt) -> GpVec {
    GpVec::from_xyz(p.xyz())
}

/// `class BlendFunc_Ruled : public Blend_Function` (`hxx:34-207`).
pub struct Ruled {
    /// `surf1` (`hxx:188`).
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:189`).
    surf2: Arc<dyn Surface>,
    /// `curv` (`hxx:190`): the guide line.
    curv: Arc<dyn Curve>,
    /// `pts1` (`hxx:191`): the extremity on the first surface, written by
    /// `Value` / `Values` / `IsSolution`.
    pts1: GpPnt,
    /// `pts2` (`hxx:192`).
    pts2: GpPnt,
    /// `istangent` (`hxx:193`).
    istangent: bool,
    /// `tg1` (`hxx:194`).
    tg1: GpVec,
    /// `tg12d` (`hxx:195`).
    tg12d: GpVec2d,
    /// `tg2` (`hxx:196`).
    tg2: GpVec,
    /// `tg22d` (`hxx:197`).
    tg22d: GpVec2d,
    /// `ptgui` (`hxx:198`): the guide point at the current parameter, written
    /// by `Set(Param)`.
    ptgui: GpPnt,
    /// `d1gui` (`hxx:199`).
    d1gui: GpVec,
    /// `d2gui` (`hxx:200`).
    d2gui: GpVec,
    /// `nplan` (`hxx:201`): the unit tangent of the guide line.
    nplan: GpVec,
    /// `normtg` (`hxx:202`): `|d1gui|`.
    normtg: f64,
    /// `theD` (`hxx:203`): `-(nplan . ptgui)`.
    the_d: f64,
    /// `distmin` (`hxx:204`): `RealLast()` at construction and never written
    /// (`cxx:42`, `:175`).
    distmin: f64,
}

impl Ruled {
    /// `BlendFunc_Ruled(S1, S2, C)` (`cxx:31-43`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf1,
            surf2,
            curv,
            pts1: GpPnt::zero(),
            pts2: GpPnt::zero(),
            istangent: true,
            tg1: GpVec::zero(),
            tg12d: GpVec2d::zero(),
            tg2: GpVec::zero(),
            tg22d: GpVec2d::zero(),
            ptgui: GpPnt::zero(),
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            nplan: GpVec::zero(),
            normtg: 0.0,
            the_d: 0.0,
            // `distmin(RealLast())` (`cxx:42`).
            distmin: f64::MAX,
        }
    }

    /// `Set(Param)` (`cxx:51-58`): `curv->D2(Param, ptgui, d1gui, d2gui)`, then
    /// the section frame `normtg` / `nplan` / `theD`. `istangent` is reset.
    pub fn set_param(&mut self, param: f64) {
        let (ptgui, d1gui, d2gui) = self.curv.d2(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        self.the_d = -self.nplan.dot(&pnt_vec(&self.ptgui));
        self.istangent = true;
    }
}

impl MathFunctionSetWithDerivatives for Ruled {
    /// `Blend_Function::NbVariables()` = `4` (`Blend_Function.cxx:19-22`).
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`cxx:45-48`): `4`.
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:178-201`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:203-294`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, None, d)
    }

    /// `Values(X, F, D)` (`cxx:296-393`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, Some(f), d)
    }
}

impl BlendAppFunction for Ruled {
    /// `Set(Param)` (`cxx:51-58`).
    fn set(&mut self, param: f64) {
        self.set_param(param);
    }

    /// `Set(First, Last)` (`cxx:60-63`): `throw Standard_NotImplemented`.
    fn set_interval(&mut self, _first: f64, _last: f64) {
        panic!("BlendFunc_Ruled::Set(First, Last): Standard_NotImplemented");
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:65-72`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:74-93`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:95-168`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`cxx:173-176`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Ok(self.distmin)
    }

    /// `Pnt1()`: `Blend_Function::Pnt1()` = `PointOnS1()`
    /// (`Blend_Function.cxx:24-27`).
    fn pnt1(&self) -> &GpPnt {
        &self.pts1
    }

    /// `Pnt2()`: `Blend_Function::Pnt2()` = `PointOnS2()`
    /// (`Blend_Function.cxx:29-32`).
    fn pnt2(&self) -> &GpPnt {
        &self.pts2
    }

    /// `IsRational()` (`cxx:572-575`): `false`.
    fn is_rational(&self) -> bool {
        false
    }

    /// `GetSectionSize()` (`cxx:580-583`): `throw Standard_NotImplemented`.
    fn get_section_size(&self) -> f64 {
        panic!("BlendFunc_Ruled::GetSectionSize(): Standard_NotImplemented");
    }

    /// `GetMinimalWeight(Weights)` (`cxx:586-589`): every weight is `1.`.
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        for w in weights.iter_mut() {
            *w = 1.0;
        }
    }

    /// `NbIntervals(S)` (`cxx:593-596`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.guide_nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:600-603`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.guide_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:606-612`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        self.section_shape(nb_poles, nb_knots, degree, nb_poles2d);
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1D)`
    /// (`cxx:618-626`): only `Tol3d` is written.
    fn get_tolerance_3d(
        &self,
        bound_tol: f64,
        _surf_tol: f64,
        _angle_tol: f64,
        tol3d: &mut MathVector,
        _tol1d: &mut MathVector,
    ) {
        tol3d.init(bound_tol);
    }

    /// `Knots(TKnots)` (`cxx:628-632`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:634-637`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:652-690`).
    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        dpoles: &mut [GpVec],
        poles2d: &mut [GpPnt2d],
        dpoles2d: &mut [GpVec2d],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool {
        self.section_with_derivatives(p, poles, dpoles, poles2d, dpoles2d, weights, dweights)
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:692-711`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, Poles, DPoles, D2Poles, ...)` (`cxx:638-650`): `false`.
    #[allow(clippy::too_many_arguments)]
    fn section_with_derivatives2(
        &mut self,
        _p: &BlendPoint,
        _poles: &mut [GpPnt],
        _dpoles: &mut [GpVec],
        _d2poles: &mut [GpVec],
        _poles2d: &mut [GpPnt2d],
        _dpoles2d: &mut [GpVec2d],
        _d2poles2d: &mut [GpVec2d],
        _weights: &mut [f64],
        _dweights: &mut [f64],
        _d2weights: &mut [f64],
    ) -> bool {
        false
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:731-744`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendFunction for Ruled {
    /// `PointOnS1()` (`cxx:416-421`).
    fn point_on_s1(&self) -> &GpPnt {
        &self.pts1
    }

    /// `PointOnS2()` (`cxx:423-428`).
    fn point_on_s2(&self) -> &GpPnt {
        &self.pts2
    }

    /// `IsTangencyPoint()` (`cxx:430-433`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS1()` (`cxx:435-441`).
    fn tangent_on_s1(&self) -> &GpVec {
        &self.tg1
    }

    /// `Tangent2dOnS1()` (`cxx:451-458`).
    fn tangent2d_on_s1(&self) -> &GpVec2d {
        &self.tg12d
    }

    /// `TangentOnS2()` (`cxx:443-449`).
    fn tangent_on_s2(&self) -> &GpVec {
        &self.tg2
    }

    /// `Tangent2dOnS2()` (`cxx:460-466`).
    fn tangent2d_on_s2(&self) -> &GpVec2d {
        &self.tg22d
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:396-414`).
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
    ) {
        self.compute_tangent(u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last);
    }

    // `TwistOnS1` / `TwistOnS2` keep the `Blend_Function` defaults (`false`);
    // `BlendFunc_Ruled.hxx:34-207` does not override them.
}
