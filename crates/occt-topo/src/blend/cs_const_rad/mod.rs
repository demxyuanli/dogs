//! `BlendFunc_CSConstRad` — the constant-radius blending function between a
//! curve and a surface, driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSConstRad.hxx:36-222`,
//! `BlendFunc_CSConstRad.cxx:37-989` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(U, V, W)`: the
//! parameters of the point on the surface and the parameter of the point on the
//! curve. With `nplan` the unit tangent of the guide line at the current
//! parameter, `theD = -(nplan . ptgui)` and `ns` the surface normal projected
//! into the section plane, the three equations are
//! 1. `nplan . pts + theD` — the surface point lies in the section plane,
//! 2. `nplan . ptc + theD` — the curve point lies in the section plane,
//! 3. `|ray * ns + (pts - ptc)|^2 - ray^2` — the curve point is at the distance
//!    `|ray|` from the centre of the arc.
//!
//! Translation notes:
//! - `Value` (`cxx:228-247`), `Derivatives` (`cxx:251-306`) and `Values`
//!   (`cxx:312-375`) are three copies in OCCT and stay three bodies here: the
//!   projected normal is normalised term by term in `Value`
//!   (`nplan.Dot(ns) / norm` then `-1. / norm`) but built first and divided
//!   afterwards in `Derivatives` / `Values` (`SetLinearForm` then `Divide`), so
//!   the rounding differs and is kept as written.
//! - There is no `Eps` guard around `|nplan x ns|` in this class (unlike
//!   `BlendFunc_ConstRad` / `BlendFunc_CSCircular`): a degenerate projected
//!   normal divides by zero exactly where OCCT does.
//! - `Set(First, Last)` throws `Standard_NotImplemented` (`cxx:95-98`); the
//!   port panics with the same message.
//! - `Section(P, ..., D2Poles, ...)` throws `Standard_DomainError`
//!   (`cxx:506-517`); the port panics with the same message.
//! - `TangentOnS` / `TangentOnC` / `Tangent2d` throw `Standard_DomainError`
//!   when `istangent` (`cxx:414-443`); the trait returns a value and every OCCT
//!   caller guards the call with `!IsTangencyPoint()`, so the stored value is
//!   returned.
//! - `GetSection` (`cxx:523-670`) and the analytic `Section(Param, ...)`
//!   (`cxx:467-504`) are kept as inherent methods: the former has no virtual
//!   slot, the latter is called by the visualization code only.
//! - `GeomFill::GetCircle` answers `Result` here because a pole array of the
//!   wrong length is OCCT's `Standard_OutOfRange`; every caller sizes the array
//!   from `GetShape`, so the failure is discarded the way the rest of the blend
//!   package discards it.

pub mod constants;
mod equations;
mod section;
mod solution;

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{BlendAppFunction, BlendCsFunction, BlendFuncSectionShape, BlendPoint};

/// `class BlendFunc_CSConstRad : public Blend_CSFunction` (`hxx:36-222`).
pub struct CsConstRad {
    /// `surf` (`hxx:199`).
    surf: Arc<dyn Surface>,
    /// `curv` (`hxx:200`).
    curv: Arc<dyn Curve>,
    /// `guide` (`hxx:201`).
    guide: Arc<dyn Curve>,
    /// `pts` (`hxx:202`): the point on the surface.
    pts: GpPnt,
    /// `ptc` (`hxx:203`): the point on the curve.
    ptc: GpPnt,
    /// `pt2d` (`hxx:204`).
    pt2d: GpPnt2d,
    /// `prmc` (`hxx:205`): the parameter on the curve.
    prmc: f64,
    /// `istangent` (`hxx:206`).
    istangent: bool,
    /// `tgs` (`hxx:207`).
    tgs: GpVec,
    /// `tg2d` (`hxx:208`).
    tg2d: GpVec2d,
    /// `tgc` (`hxx:209`).
    tgc: GpVec,
    /// `ray` (`hxx:210`): the signed radius of the section.
    ray: f64,
    /// `choix` (`hxx:211`): the fillet quadrant.
    choix: i32,
    /// `ptgui` (`hxx:212`): the guide-line point.
    ptgui: GpPnt,
    /// `d1gui` (`hxx:213`).
    d1gui: GpVec,
    /// `d2gui` (`hxx:214`).
    d2gui: GpVec,
    /// `nplan` (`hxx:215`): the unit tangent of the guide line.
    nplan: GpVec,
    /// `normtg` (`hxx:216`).
    normtg: f64,
    /// `theD` (`hxx:217`).
    the_d: f64,
    /// `maxang` (`hxx:218`): the widest section angle seen by `IsSolution`.
    maxang: f64,
    /// `minang` (`hxx:219`): the narrowest section angle.
    minang: f64,
    /// `mySShape` (`hxx:220`).
    my_sshape: BlendFuncSectionShape,
    /// `myTConv` (`hxx:221`): `Convert_TgtThetaOver2` after construction
    /// (`cxx:57`).
    my_tconv: ParameterisationType,
    /// Scratch storage for `Value` / `Derivatives` / `Values`: `hxx` declares
    /// no `E` / `DEDX` for this class (the three bodies use local `gp_Vec`s and
    /// `IsSolution` its own `valsol` / `gradsol`, `cxx:145-147`), so this is
    /// only the port's output buffer, never read across calls.
    e: MathVector,
    /// Output buffer for the jacobian; see `e`.
    dedx: MathMatrix,
}

impl CsConstRad {
    /// `BlendFunc_CSConstRad(S, C, CG)` (`cxx:39-57`).
    pub fn new(surf: Arc<dyn Surface>, curv: Arc<dyn Curve>, guide: Arc<dyn Curve>) -> Self {
        Self {
            surf,
            curv,
            guide,
            pts: GpPnt::zero(),
            ptc: GpPnt::zero(),
            pt2d: GpPnt2d::new(0.0, 0.0),
            prmc: 0.0,
            istangent: true,
            tgs: GpVec::zero(),
            tg2d: GpVec2d::new(0.0, 0.0),
            tgc: GpVec::zero(),
            ray: 0.0,
            choix: 0,
            ptgui: GpPnt::zero(),
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            nplan: GpVec::zero(),
            normtg: 0.0,
            the_d: 0.0,
            maxang: f64::MIN,
            minang: f64::MAX,
            my_sshape: BlendFuncSectionShape::Rational,
            my_tconv: ParameterisationType::TgtThetaOver2,
            e: MathVector::new(1, 3),
            dedx: MathMatrix::new(1, 3, 1, 3),
        }
    }

    /// `Set(Radius, Choix)` (`cxx:68-73`): the radius is always negative; the
    /// quadrant only changes `choix`.
    pub fn set_radius(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        self.ray = -radius.abs();
    }

    /// `Set(TypeSection)` (`cxx:78-81`).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_sshape = type_section;
    }

    /// `Set(Param)` (`cxx:85-91`): the section plane.
    pub fn set_param(&mut self, param: f64) {
        let (ptgui, d1gui, d2gui) = self.guide.d2(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        self.the_d = -self.nplan.dot(&pnt_vec(&self.ptgui));
    }

    /// `ray` (`hxx:210`).
    pub fn ray(&self) -> f64 {
        self.ray
    }

    /// `choix` (`hxx:211`).
    pub fn choix(&self) -> i32 {
        self.choix
    }

    /// `maxang` (`hxx:218`).
    pub fn max_ang(&self) -> f64 {
        self.maxang
    }

    /// `minang` (`hxx:219`).
    pub fn min_ang(&self) -> f64 {
        self.minang
    }
}

/// `p.XYZ()` as a `gp_Vec` — `nplan.XYZ().Dot(pts.XYZ())` (`cxx:234-235`).
pub(super) fn pnt_vec(p: &GpPnt) -> GpVec {
    GpVec::from_xyz(p.xyz())
}

impl MathFunctionSetWithDerivatives for CsConstRad {
    /// `Blend_CSFunction::NbVariables()` = `3` (`Blend_CSFunction.cxx:20-23`);
    /// this class does not override it (`hxx:36-222`).
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:63-66`): `3`.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:228-247`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let ok = self.evaluate_value(x);
        f.initialized(&self.e);
        ok
    }

    /// `Derivatives(X, D)` (`cxx:251-306`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let ok = self.evaluate_derivatives(x);
        for i in 1..=3 {
            for j in 1..=3 {
                d.set_value(i, j, self.dedx.value(i, j));
            }
        }
        ok
    }

    /// `Values(X, F, D)` (`cxx:312-375`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let ok = self.evaluate_values(x);
        f.initialized(&self.e);
        for i in 1..=3 {
            for j in 1..=3 {
                d.set_value(i, j, self.dedx.value(i, j));
            }
        }
        ok
    }
}

impl BlendAppFunction for CsConstRad {
    /// `Set(Param)` (`cxx:85-91`).
    fn set(&mut self, param: f64) {
        self.set_param(param);
    }

    /// `Set(First, Last)` (`cxx:95-98`) throws `Standard_NotImplemented`.
    fn set_interval(&mut self, _first: f64, _last: f64) {
        panic!("Standard_NotImplemented: BlendFunc_CSConstRad::Set");
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:101-106`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:111-133`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:135-224`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()`: `Blend_CSFunction::GetMinimalDistance` throws
    /// `Standard_NotImplemented` (`Blend_CSFunction.cxx:49-52`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Err("Standard_NotImplemented: Blend_CSFunction::GetMinimalDistance".into())
    }

    /// `Pnt1()`: `Blend_CSFunction::Pnt1()` = `PointOnC()`
    /// (`Blend_CSFunction.cxx:25-28`).
    fn pnt1(&self) -> &GpPnt {
        &self.ptc
    }

    /// `Pnt2()`: `Blend_CSFunction::Pnt2()` = `PointOnS()`
    /// (`Blend_CSFunction.cxx:30-33`).
    fn pnt2(&self) -> &GpPnt {
        &self.pts
    }

    /// `IsRational()` (`cxx:673-676`).
    fn is_rational(&self) -> bool {
        self.shape_is_rational()
    }

    /// `GetSectionSize()` (`cxx:680-683`): `maxang * |ray|`.
    fn get_section_size(&self) -> f64 {
        self.section_size()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:687-691`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`cxx:695-698`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.curve_nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:702-705`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.curve_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:710-714`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        self.section_shape(nb_poles, nb_knots, degree, nb_poles2d);
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)`
    /// (`cxx:718-733`).
    fn get_tolerance_3d(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut MathVector,
        tol1d: &mut MathVector,
    ) {
        self.fill_section_tolerance(bound_tol, surf_tol, angle_tol, tol3d, tol1d);
    }

    /// `Knots(TKnots)` (`cxx:736-739`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:742-745`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:804-981`).
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

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:750-802`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` throws
    /// `Standard_DomainError("BlendFunc_CSConstRad::Section : Not implemented")`
    /// (`cxx:506-517`).
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
        panic!("Standard_DomainError: BlendFunc_CSConstRad::Section : Not implemented");
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:983-988`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendCsFunction for CsConstRad {
    /// `PointOnS()` (`cxx:377-380`).
    fn point_on_s(&self) -> &GpPnt {
        &self.pts
    }

    /// `PointOnC()` (`cxx:384-387`).
    fn point_on_c(&self) -> &GpPnt {
        &self.ptc
    }

    /// `Pnt2d()` (`cxx:391-394`).
    fn pnt2d(&self) -> &GpPnt2d {
        &self.pt2d
    }

    /// `ParameterOnC()` (`cxx:398-401`).
    fn parameter_on_c(&self) -> f64 {
        self.prmc
    }

    /// `IsTangencyPoint()` (`cxx:405-408`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS()` (`cxx:412-420`).
    fn tangent_on_s(&self) -> &GpVec {
        &self.tgs
    }

    /// `Tangent2d()` (`cxx:435-443`).
    fn tangent2d(&self) -> &GpVec2d {
        &self.tg2d
    }

    /// `TangentOnC()` (`cxx:424-433`).
    fn tangent_on_c(&self) -> &GpVec {
        &self.tgc
    }

    /// `Tangent(U, V, TgS, NormS)` (`cxx:447-463`).
    fn tangent(&self, u: f64, v: f64, tg_s: &mut GpVec, norm_s: &mut GpVec) {
        self.compute_tangent(u, v, tg_s, norm_s);
    }
}
