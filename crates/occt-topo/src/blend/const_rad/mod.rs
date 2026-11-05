//! `BlendFunc_ConstRad` — the constant-radius blending function between two
//! surfaces, driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRad.hxx:37-279`,
//! `BlendFunc_ConstRad.cxx:44-1956` (OCCT 8.0.0).
//!
//! The vector `X` of `Value` / `Derivatives` / `Values` is `(U1, V1, U2, V2)`:
//! the parameters of the two extremities of the section, one on each surface.
//! With `nplan` the unit tangent of the guide line at the current parameter,
//! `theD = -(nplan . ptgui)` and `ns1` / `ns2` the surface normals projected
//! into the section plane, the four equations are
//! 1. `(nplan . (pts1 + pts2)) / 2 + theD` — the mid-point of the two
//!    extremities lies in the section plane,
//! 2. `ray1 * ns1 - ray2 * ns2 + (pts1 - pts2) = 0` — the two extremities are
//!    the endpoints of a circular arc of radius `|ray1|` centred on the
//!    guide-line point.
//!
//! Translation notes:
//! - `ComputeValues` (`cxx:136-777`) is the mandatory entry point of every
//!   equation: it caches the last `(X, T)` and the derivative order reached
//!   (`myXOrder` / `myTOrder`), then writes `E`, `DEDX`, `DEDT`, `D2EDX2`,
//!   `D2EDXDT` and `D2EDT2`. The port keeps the cache and its two sentinels
//!   (`xval` / `tval`, `CACHE_EMPTY`).
//! - OCCT keeps `ptgui`, `d1gui`, `d2gui`, `d3gui`, `invnormtg`, `dinvnormtg`
//!   as function-local `static`s (`cxx:142-147`) so that they survive a cache
//!   hit. `d3gui` and `dinvnormtg` are written and read inside the same
//!   `case 2` arm, so they are locals here; the other four are fields. A cache
//!   hit on `t` requires `T == tval`, so the values carried over are exactly
//!   the ones the `static`s would hold.
//! - `surf1->D3` / `surf2->D3` (`cxx:240-248`) is [`Surface::d3`], which
//!   answers `None` for a surface whose `EvalD3` is not translated. OCCT
//!   raises `Geom_UndefinedDerivative` there, so the port panics at the same
//!   place instead of inventing a zero third derivative.
//! - `BlendFunc::ComputeNormal` / `ComputeDNormal` (`cxx:260-278`) return
//!   `false` when `CSLib` cannot define the normal, and OCCT then leaves
//!   `nsurf1` / `nsurf2` untouched; the port keeps that by writing only the
//!   `Some` answer. A throwing `Surf->DN` is the port's `Err`, which OCCT lets
//!   escape, so the panic message carries it.
//! - `BlendFunc_Tensor` (`hxx:23`, `hxx:269`) is structurally the same class as
//!   `GeomFill_Tensor`; the port reuses [`occt_geom::geom_fill::Tensor`].
//! - `dnsurf1` / `dnsurf2` (`hxx:250-251`) are never written nor read in OCCT
//!   (the same holds for the commented-out `gp_Vec normal` locals at
//!   `cxx:257-269`), so they are not carried.
//! - `TangentOnS1` / `Tangent2dOnS1` / `TangentOnS2` / `Tangent2dOnS2` /
//!   `TwistOnS1` / `TwistOnS2` raise `Standard_DomainError` when `istangent`
//!   (`cxx:1027-1137`); the trait returns a value, and every OCCT caller
//!   guards the call with `!IsTangencyPoint()`, so the stored value is
//!   returned.
//! - `Section(First, Last)` keeps `tcurv = curv->Trim(First, Last, 1.e-12)`
//!   (`cxx:790-793`). `GeomAdaptor_Curve::Trim` keeps the basis
//!   parameterisation and this class never asks `tcurv` for its bounds, so the
//!   port evaluates `curv` and keeps the interval as documentation.

pub mod constants;
mod compute_values;
mod compute_values_order2;
mod const_rad_inv;
mod section;
mod section_circ;
mod section_order2;
mod solution;

pub use const_rad_inv::ConstRadInv;
pub use constants::{CACHE_EMPTY, EPS, GAUSS_SOLUTION_EPS, POINTS_EQUAL_TOL, SVD_EPS};

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::geom_fill::Tensor;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{BlendAppFunction, BlendFuncSectionShape, BlendFunction, BlendPoint};

/// `class BlendFunc_ConstRad : public Blend_Function` (`hxx:37-279`).
pub struct ConstRad {
    /// `surf1` (`hxx:208`).
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:209`).
    surf2: Arc<dyn Surface>,
    /// `curv` (`hxx:210`): the guide line.
    curv: Arc<dyn Curve>,
    /// `tcurv` (`hxx:211`): `curv->Trim(First, Last, 1.e-12)`; see the module
    /// note -- evaluation is identical, so it aliases `curv`.
    tcurv: Arc<dyn Curve>,
    /// `pts1` (`hxx:212`).
    pts1: GpPnt,
    /// `pts2` (`hxx:213`).
    pts2: GpPnt,
    /// `istangent` (`hxx:214`).
    istangent: bool,
    /// `tg1` (`hxx:215`).
    tg1: GpVec,
    /// `tg12d` (`hxx:216`).
    tg12d: GpVec2d,
    /// `tg2` (`hxx:217`).
    tg2: GpVec,
    /// `tg22d` (`hxx:218`).
    tg22d: GpVec2d,
    /// `param` (`hxx:219`): the current guide parameter, set by `Set(Param)`.
    param: f64,
    /// `ray1` (`hxx:220`): the signed radius on the first surface.
    ray1: f64,
    /// `ray2` (`hxx:221`): the signed radius on the second surface.
    ray2: f64,
    /// `choix` (`hxx:222`): the fillet quadrant.
    choix: i32,
    /// `myXOrder` (`hxx:223`): the derivative order cached in `xval`.
    my_x_order: i32,
    /// `myTOrder` (`hxx:224`): the derivative order cached in `tval`.
    my_t_order: i32,
    /// `xval` (`hxx:225`).
    x_val: MathVector,
    /// `tval` (`hxx:226`).
    tval: f64,
    /// `ptgui` (`cxx:144`, a function-local `static`): the guide point.
    ptgui: GpPnt,
    /// `d1gui` (`cxx:145`): the first derivative of the guide line.
    d1gui: GpVec,
    /// `d2gui` (`cxx:145`): the second derivative of the guide line.
    d2gui: GpVec,
    /// `invnormtg` (`cxx:145`): `1 / |d1gui|`.
    invnormtg: f64,
    /// `d1u1` (`hxx:227`).
    d1u1: GpVec,
    /// `d1u2` (`hxx:228`).
    d1u2: GpVec,
    /// `d1v1` (`hxx:229`).
    d1v1: GpVec,
    /// `d1v2` (`hxx:230`).
    d1v2: GpVec,
    /// `d2u1` (`hxx:231`).
    d2u1: GpVec,
    /// `d2v1` (`hxx:232`).
    d2v1: GpVec,
    /// `d2uv1` (`hxx:233`).
    d2uv1: GpVec,
    /// `d2u2` (`hxx:234`).
    d2u2: GpVec,
    /// `d2v2` (`hxx:235`).
    d2v2: GpVec,
    /// `d2uv2` (`hxx:236`).
    d2uv2: GpVec,
    /// `d3u1` (`cxx:142`, a function-local `static`): a third derivative of
    /// `surf1`, written by the `Order == 2` arm of the `X` positioning and read
    /// by [`Self::compute_second_order`].
    d3u1: GpVec,
    /// `d3v1` (`cxx:142`).
    d3v1: GpVec,
    /// `d3uuv1` (`cxx:142`).
    d3uuv1: GpVec,
    /// `d3uvv1` (`cxx:142`).
    d3uvv1: GpVec,
    /// `d3u2` (`cxx:142`).
    d3u2: GpVec,
    /// `d3v2` (`cxx:142`).
    d3v2: GpVec,
    /// `d3uuv2` (`cxx:142`).
    d3uuv2: GpVec,
    /// `d3uvv2` (`cxx:142`).
    d3uvv2: GpVec,
    /// `dn1w` (`hxx:237`): the derivative of the projected normal of `surf1`.
    dn1w: GpVec,
    /// `dn2w` (`hxx:238`).
    dn2w: GpVec,
    /// `d2n1w` (`hxx:239`).
    d2n1w: GpVec,
    /// `d2n2w` (`hxx:240`).
    d2n2w: GpVec,
    /// `nplan` (`hxx:241`): the unit tangent of the guide line.
    nplan: GpVec,
    /// `nsurf1` (`hxx:242`): the (unnormalised) normal of the first surface.
    nsurf1: GpVec,
    /// `dns1u1` (`hxx:243`).
    dns1u1: GpVec,
    /// `dns1u2` (`hxx:244`).
    dns1u2: GpVec,
    /// `dns1v1` (`hxx:245`).
    dns1v1: GpVec,
    /// `dns1v2` (`hxx:246`).
    dns1v2: GpVec,
    /// `nsurf2` (`hxx:247`).
    nsurf2: GpVec,
    /// `dnplan` (`hxx:248`).
    dnplan: GpVec,
    /// `d2nplan` (`hxx:249`).
    d2nplan: GpVec,
    /// `dndu1` (`hxx:252`).
    dndu1: GpVec,
    /// `dndu2` (`hxx:253`).
    dndu2: GpVec,
    /// `dndv1` (`hxx:254`).
    dndv1: GpVec,
    /// `dndv2` (`hxx:255`).
    dndv2: GpVec,
    /// `d2ndu1` (`hxx:256`).
    d2ndu1: GpVec,
    /// `d2ndu2` (`hxx:257`).
    d2ndu2: GpVec,
    /// `d2ndv1` (`hxx:258`).
    d2ndv1: GpVec,
    /// `d2ndv2` (`hxx:259`).
    d2ndv2: GpVec,
    /// `d2nduv1` (`hxx:260`).
    d2nduv1: GpVec,
    /// `d2nduv2` (`hxx:261`).
    d2nduv2: GpVec,
    /// `d2ndtu1` (`hxx:262`).
    d2ndtu1: GpVec,
    /// `d2ndtu2` (`hxx:263`).
    d2ndtu2: GpVec,
    /// `d2ndtv1` (`hxx:264`).
    d2ndtv1: GpVec,
    /// `d2ndtv2` (`hxx:264`).
    d2ndtv2: GpVec,
    /// `E` (`hxx:266`): the four equation values.
    e: MathVector,
    /// `DEDX` (`hxx:267`).
    dedx: MathMatrix,
    /// `DEDT` (`hxx:268`).
    dedt: MathVector,
    /// `D2EDX2` (`hxx:269`).
    d2edx2: Tensor,
    /// `D2EDXDT` (`hxx:270`).
    d2edxdt: MathMatrix,
    /// `D2EDT2` (`hxx:271`).
    d2edt2: MathVector,
    /// `maxang` (`hxx:272`): the widest section angle seen by `IsSolution`.
    maxang: f64,
    /// `minang` (`hxx:273`): the narrowest section angle.
    minang: f64,
    /// `distmin` (`hxx:274`).
    distmin: f64,
    /// `mySShape` (`hxx:275`).
    my_sshape: BlendFuncSectionShape,
    /// `myTConv` (`hxx:276`): `Convert_TgtThetaOver2` after construction
    /// (`cxx:75`).
    my_tconv: ParameterisationType,
}

impl ConstRad {
    /// `BlendFunc_ConstRad(S1, S2, C)` (`cxx:46-77`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        let mut x_val = MathVector::new(1, 4);
        // `xval.Init(-9.876e100)` (`cxx:71`).
        x_val.init(CACHE_EMPTY);
        Self {
            tcurv: curv.clone(),
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
            param: 0.0,
            ray1: 0.0,
            ray2: 0.0,
            choix: 0,
            my_x_order: -1,
            my_t_order: -1,
            x_val,
            tval: CACHE_EMPTY,
            ptgui: GpPnt::zero(),
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            invnormtg: 0.0,
            d1u1: GpVec::zero(),
            d1u2: GpVec::zero(),
            d1v1: GpVec::zero(),
            d1v2: GpVec::zero(),
            d2u1: GpVec::zero(),
            d2v1: GpVec::zero(),
            d2uv1: GpVec::zero(),
            d2u2: GpVec::zero(),
            d2v2: GpVec::zero(),
            d2uv2: GpVec::zero(),
            d3u1: GpVec::zero(),
            d3v1: GpVec::zero(),
            d3uuv1: GpVec::zero(),
            d3uvv1: GpVec::zero(),
            d3u2: GpVec::zero(),
            d3v2: GpVec::zero(),
            d3uuv2: GpVec::zero(),
            d3uvv2: GpVec::zero(),
            dn1w: GpVec::zero(),
            dn2w: GpVec::zero(),
            d2n1w: GpVec::zero(),
            d2n2w: GpVec::zero(),
            nplan: GpVec::zero(),
            nsurf1: GpVec::zero(),
            dns1u1: GpVec::zero(),
            dns1u2: GpVec::zero(),
            dns1v1: GpVec::zero(),
            dns1v2: GpVec::zero(),
            nsurf2: GpVec::zero(),
            dnplan: GpVec::zero(),
            d2nplan: GpVec::zero(),
            dndu1: GpVec::zero(),
            dndu2: GpVec::zero(),
            dndv1: GpVec::zero(),
            dndv2: GpVec::zero(),
            d2ndu1: GpVec::zero(),
            d2ndu2: GpVec::zero(),
            d2ndv1: GpVec::zero(),
            d2ndv2: GpVec::zero(),
            d2nduv1: GpVec::zero(),
            d2nduv2: GpVec::zero(),
            d2ndtu1: GpVec::zero(),
            d2ndtu2: GpVec::zero(),
            d2ndtv1: GpVec::zero(),
            d2ndtv2: GpVec::zero(),
            e: MathVector::new(1, 4),
            dedx: MathMatrix::new(1, 4, 1, 4),
            dedt: MathVector::new(1, 4),
            d2edx2: Tensor::new(4, 4, 4),
            d2edxdt: MathMatrix::new(1, 4, 1, 4),
            d2edt2: MathVector::new(1, 4),
            maxang: f64::MIN,
            minang: f64::MAX,
            distmin: f64::MAX,
            my_sshape: BlendFuncSectionShape::Rational,
            my_tconv: ParameterisationType::TgtThetaOver2,
        }
    }

    /// `Set(Radius, Choix)` (`cxx:86-118`): the sign of each radius follows the
    /// fillet quadrant.
    pub fn set_radius(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => {
                self.ray1 = -radius;
                self.ray2 = -radius;
            }
            3 | 4 => {
                self.ray1 = radius;
                self.ray2 = -radius;
            }
            5 | 6 => {
                self.ray1 = radius;
                self.ray2 = radius;
            }
            7 | 8 => {
                self.ray1 = -radius;
                self.ray2 = radius;
            }
            _ => {
                self.ray1 = -radius;
                self.ray2 = -radius;
            }
        }
    }

    /// `Set(TypeSection)` (`cxx:122-125`).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_sshape = type_section;
    }

    /// `ray1` (`hxx:220`).
    pub fn ray1(&self) -> f64 {
        self.ray1
    }

    /// `ray2` (`hxx:221`).
    pub fn ray2(&self) -> f64 {
        self.ray2
    }

    /// `choix` (`hxx:222`).
    pub fn choix(&self) -> i32 {
        self.choix
    }

    /// `maxang` (`hxx:272`).
    pub fn max_ang(&self) -> f64 {
        self.maxang
    }

    /// `minang` (`hxx:273`).
    pub fn min_ang(&self) -> f64 {
        self.minang
    }
}

/// `nplan.XYZ().Dot(ptgui.XYZ())` and the two-point mid-point expressions of
/// `ComputeValues` (`cxx:287-291`).
pub(super) fn pnt_vec(p: &GpPnt) -> GpVec {
    GpVec::from_xyz(p.xyz())
}

impl MathFunctionSetWithDerivatives for ConstRad {
    /// `Blend_Function::NbVariables()` = `4` (`Blend_Function.cxx:19-22`).
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`cxx:79-82`): `4`.
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:974-979`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let ok = self.compute_values(x, 0, false, 0.0);
        f.initialized(&self.e);
        ok
    }

    /// `Derivatives(X, D)` (`cxx:984-988`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let ok = self.compute_values(x, 1, false, 0.0);
        for i in 1..=4 {
            for j in 1..=4 {
                d.set_value(i, j, self.dedx.value(i, j));
            }
        }
        ok
    }

    /// `Values(X, F, D)` (`cxx:992-997`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let ok = self.compute_values(x, 1, false, 0.0);
        f.initialized(&self.e);
        for i in 1..=4 {
            for j in 1..=4 {
                d.set_value(i, j, self.dedx.value(i, j));
            }
        }
        ok
    }
}

impl BlendAppFunction for ConstRad {
    /// `Set(Param)` (`cxx:781-784`): the current guide parameter.
    fn set(&mut self, param: f64) {
        self.param = param;
    }

    /// `Set(First, Last)` (`cxx:790-793`): `tcurv = curv->Trim(First, Last,
    /// 1.e-12)`. See the module note: `Trim` does not change the evaluated
    /// values, so `tcurv` keeps aliasing `curv`.
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// `GetTolerance(Tolerance, Tol)` (`cxx:797-804`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.fill_parametric_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:808-828`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:833-962`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`cxx:967-970`).
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

    /// `IsRational()` (`cxx:1204-1207`).
    fn is_rational(&self) -> bool {
        self.shape_is_rational()
    }

    /// `GetSectionSize()` (`cxx:1211-1214`): `maxang * |ray1|`.
    fn get_section_size(&self) -> f64 {
        self.section_size()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:1218-1221`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`cxx:1225-1228`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.guide_nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`cxx:1232-1235`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.guide_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:1239-1243`).
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
    /// (`cxx:1249-1262`).
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

    /// `Knots(TKnots)` (`cxx:1266-1269`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.section_knots(tknots);
    }

    /// `Mults(TMults)` (`cxx:1273-1276`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:1356-1524`).
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

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:1284-1352`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` (`cxx:1527-1905`).
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
    ) -> bool {
        self.section_with_second_derivatives(
            p, poles, dpoles, d2poles, poles2d, dpoles2d, d2poles2d, weights, dweights, d2weights,
        )
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1942-1956`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.compute_resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendFunction for ConstRad {
    /// `PointOnS1()` (`cxx:1003-1006`).
    fn point_on_s1(&self) -> &GpPnt {
        &self.pts1
    }

    /// `PointOnS2()` (`cxx:1012-1015`).
    fn point_on_s2(&self) -> &GpPnt {
        &self.pts2
    }

    /// `IsTangencyPoint()` (`cxx:1019-1022`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS1()` (`cxx:1026-1033`).
    fn tangent_on_s1(&self) -> &GpVec {
        &self.tg1
    }

    /// `Tangent2dOnS1()` (`cxx:1047-1055`).
    fn tangent2d_on_s1(&self) -> &GpVec2d {
        &self.tg12d
    }

    /// `TangentOnS2()` (`cxx:1037-1043`).
    fn tangent_on_s2(&self) -> &GpVec {
        &self.tg2
    }

    /// `Tangent2dOnS2()` (`cxx:1060-1067`).
    fn tangent2d_on_s2(&self) -> &GpVec2d {
        &self.tg22d
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:1064-1114`).
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

    /// `TwistOnS1()` (`cxx:1119-1126`).
    fn twist_on_s1(&self) -> bool {
        self.tg1.dot(&self.nplan) < 0.0
    }

    /// `TwistOnS2()` (`cxx:1131-1138`).
    fn twist_on_s2(&self) -> bool {
        self.tg2.dot(&self.nplan) < 0.0
    }
}
