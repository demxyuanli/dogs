//! `BlendFunc_ConstThroat` — the symmetric chamfer whose constant throat is
//! the height of the isosceles triangle of the section.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstThroat.hxx:27-110`,
//! `BlendFunc_ConstThroat.cxx:27-281` (OCCT 8.0.0).
//!
//! The four equations, with `ptgui` the guide-line point, `nplan` its unit
//! tangent and `theD = -(nplan . ptgui)`:
//! 1. `nplan . pts1 + theD` — the first point lies in the section plane,
//! 2. `nplan . pts2 + theD` — the second point lies in the section plane,
//! 3. `|mid - ptgui|^2 - Throat^2` — the mid-point is at the throat distance,
//! 4. `|pts1 - ptgui|^2 - |pts2 - ptgui|^2` — the two are symmetric about the
//!    guide line.
//!
//! Translation notes:
//! - The trailing `#define Eps 1.e-15` of the `.cxx` is not used by any body
//!   of this class; it is not carried.
//! - `TangentOnS*` / `Tangent2dOnS*` throw `Standard_DomainError` when
//!   `istangent` (`cxx:183-224`); the trait returns a value and every OCCT
//!   caller guards the call with `!IsTangencyPoint()`, so the stored value is
//!   returned.
//! - `GetSectionSize()` throws `Standard_NotImplemented` (`cxx:276-280`); the
//!   port panics with the same message.
//! - The inherited `Set(First, Last)` is empty (`GenChamfer.cxx:50`), so
//!   `set_interval` is a no-op.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{Gauss, MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{BlendAppFunction, BlendFunction, BlendPoint};

use super::gen_chamfer::GenChamfer;
use super::shared::{self, GAUSS_MIN_PIVOT, mid_pnt, pnt_vec};

/// `class BlendFunc_ConstThroat : public BlendFunc_GenChamfer`
/// (`hxx:27-110`).
pub struct ConstThroat {
    /// The `BlendFunc_GenChamfer` base sub-object.
    pub(super) gen: GenChamfer,
    /// `pts1` (`hxx:88`): the point on the first support.
    pub(super) pts1: GpPnt,
    /// `pts2` (`hxx:89`): the point on the second support.
    pub(super) pts2: GpPnt,
    /// `d1u1` (`hxx:90`).
    pub(super) d1u1: GpVec,
    /// `d1v1` (`hxx:91`).
    pub(super) d1v1: GpVec,
    /// `d1u2` (`hxx:92`).
    pub(super) d1u2: GpVec,
    /// `d1v2` (`hxx:93`).
    pub(super) d1v2: GpVec,
    /// `istangent` (`hxx:94`).
    pub(super) istangent: bool,
    /// `tg1` (`hxx:95`).
    pub(super) tg1: GpVec,
    /// `tg12d` (`hxx:96`).
    pub(super) tg12d: GpVec2d,
    /// `tg2` (`hxx:97`).
    pub(super) tg2: GpVec,
    /// `tg22d` (`hxx:98`).
    pub(super) tg22d: GpVec2d,
    /// `param` (`hxx:99`).
    pub(super) param: f64,
    /// `Throat` (`hxx:100`).
    pub(super) throat: f64,
    /// `ptgui` (`hxx:102`).
    pub(super) ptgui: GpPnt,
    /// `nplan` (`hxx:103`).
    pub(super) nplan: GpVec,
    /// `normtg` (`hxx:104`).
    pub(super) normtg: f64,
    /// `theD` (`hxx:105`).
    pub(super) the_d: f64,
    /// `d1gui` (`hxx:106`).
    pub(super) d1gui: GpVec,
    /// `d2gui` (`hxx:107`).
    pub(super) d2gui: GpVec,
}

impl ConstThroat {
    /// `BlendFunc_ConstThroat(S1, S2, C)` (`cxx:27-37`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            gen: GenChamfer::new(surf1, surf2, curv),
            pts1: GpPnt::zero(),
            pts2: GpPnt::zero(),
            d1u1: GpVec::zero(),
            d1v1: GpVec::zero(),
            d1u2: GpVec::zero(),
            d1v2: GpVec::zero(),
            istangent: false,
            tg1: GpVec::zero(),
            tg12d: GpVec2d::new(0.0, 0.0),
            tg2: GpVec::zero(),
            tg22d: GpVec2d::new(0.0, 0.0),
            param: 0.0,
            throat: 0.0,
            ptgui: GpPnt::zero(),
            nplan: GpVec::zero(),
            normtg: 0.0,
            the_d: 0.0,
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
        }
    }

    /// `Set(aThroat, dis1, Choix)` (`cxx:42-46`): the second argument is
    /// unused in OCCT; it is only kept for the caller's `ChFi3d` interface.
    pub fn set_throat(&mut self, a_throat: f64, _dis1: f64, choix: i32) {
        self.throat = a_throat;
        self.gen.choix = choix;
    }

    /// `Set(Param)` (`cxx:48-56`): the section plane.
    pub fn set_param(&mut self, param: f64) {
        self.param = param;
        let (ptgui, d1gui, d2gui) = self.gen.curv.d2(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        self.the_d = -self.nplan.dot(&pnt_vec(&self.ptgui));
    }

    /// `Throat` (`hxx:100`).
    pub fn throat(&self) -> f64 {
        self.throat
    }

    /// `Value(X, F)` (`cxx:114-133`).
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) {
        self.pts1 = self.gen.surf1.value(x.value(1), x.value(2));
        self.pts2 = self.gen.surf2.value(x.value(3), x.value(4));

        f.set_value(1, self.nplan.dot(&pnt_vec(&self.pts1)) + self.the_d);
        f.set_value(2, self.nplan.dot(&pnt_vec(&self.pts2)) + self.the_d);

        let mid = mid_pnt(&self.pts1, &self.pts2);
        let vmid = GpVec::from_pnts(&self.ptgui, &mid);
        f.set_value(3, vmid.square_magnitude() - self.throat * self.throat);

        let vref1 = GpVec::from_pnts(&self.ptgui, &self.pts1);
        let vref2 = GpVec::from_pnts(&self.ptgui, &self.pts2);
        f.set_value(
            4,
            vref1.square_magnitude() - vref2.square_magnitude(),
        );
    }

    /// `Derivatives(X, D)` (`cxx:137-159`).
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) {
        let (p1, u1, v1) = self.gen.surf1.d1(x.value(1), x.value(2));
        self.pts1 = p1;
        self.d1u1 = u1;
        self.d1v1 = v1;
        let (p2, u2, v2) = self.gen.surf2.d1(x.value(3), x.value(4));
        self.pts2 = p2;
        self.d1u2 = u2;
        self.d1v2 = v2;

        let mid = mid_pnt(&self.pts1, &self.pts2);
        let vmid = GpVec::from_pnts(&self.ptgui, &mid);
        let v1ref = GpVec::from_pnts(&self.ptgui, &self.pts1);
        let v2ref = GpVec::from_pnts(&self.ptgui, &self.pts2);

        d.set_value(1, 1, self.nplan.dot(&self.d1u1));
        d.set_value(1, 2, self.nplan.dot(&self.d1v1));
        d.set_value(1, 3, 0.0);
        d.set_value(1, 4, 0.0);
        d.set_value(2, 1, 0.0);
        d.set_value(2, 2, 0.0);
        d.set_value(2, 3, self.nplan.dot(&self.d1u2));
        d.set_value(2, 4, self.nplan.dot(&self.d1v2));
        d.set_value(3, 1, vmid.dot(&self.d1u1));
        d.set_value(3, 2, vmid.dot(&self.d1v1));
        d.set_value(3, 3, vmid.dot(&self.d1u2));
        d.set_value(3, 4, vmid.dot(&self.d1v2));
        d.set_value(4, 1, 2.0 * v1ref.dot(&self.d1u1));
        d.set_value(4, 2, 2.0 * v1ref.dot(&self.d1v1));
        d.set_value(4, 3, -2.0 * v2ref.dot(&self.d1u2));
        d.set_value(4, 4, -2.0 * v2ref.dot(&self.d1v2));
    }

    /// `IsSolution(Sol, Tol)` (`cxx:60-108`): the four values are within
    /// `Tol` (`Tol^2` for the last two), then the section tangents are solved
    /// for from the jacobian with `math_Gauss` (pivot `1.e-20`, no fallback).
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut secmember = MathVector::new(1, 4);
        let mut valsol = MathVector::new(1, 4);
        let mut gradsol = MathMatrix::new(1, 4, 1, 4);

        self.compute_value(sol, &mut valsol);
        self.compute_derivatives(sol, &mut gradsol);

        self.gen.tol = tol;

        if valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= tol * tol
            && valsol.value(4).abs() <= tol * tol
        {
            let dnplan = GpVec::zero().set_linear_form_2(
                1.0 / self.normtg,
                &self.d2gui,
                -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
                &self.nplan,
            );

            let temp1 = GpVec::from_pnts(&self.ptgui, &self.pts1);
            let temp2 = GpVec::from_pnts(&self.ptgui, &self.pts2);
            let mid = mid_pnt(&self.pts1, &self.pts2);
            let tempmid = GpVec::from_pnts(&self.ptgui, &mid);

            let ng = self.nplan.dot(&self.d1gui);
            secmember.set_value(1, ng - dnplan.dot(&temp1));
            secmember.set_value(2, ng - dnplan.dot(&temp2));
            secmember.set_value(3, 2.0 * self.d1gui.dot(&tempmid));
            secmember.set_value(4, 2.0 * self.d1gui.dot(&temp2) - 2.0 * self.d1gui.dot(&temp1));

            let (p1, u1, v1) = self.gen.surf1.d1(sol.value(1), sol.value(2));
            self.pts1 = p1;
            self.d1u1 = u1;
            self.d1v1 = v1;
            let (p2, u2, v2) = self.gen.surf2.d1(sol.value(3), sol.value(4));
            self.pts2 = p2;
            self.d1u2 = u2;
            self.d1v2 = v2;

            let resol = Gauss::new(&gradsol, GAUSS_MIN_PIVOT);
            if resol.is_done() {
                if let Some(x) = resol.solve(&secmember) {
                    secmember = x;
                }
                self.tg1 = GpVec::zero().set_linear_form_2(
                    secmember.value(1),
                    &self.d1u1,
                    secmember.value(2),
                    &self.d1v1,
                );
                self.tg2 = GpVec::zero().set_linear_form_2(
                    secmember.value(3),
                    &self.d1u2,
                    secmember.value(4),
                    &self.d1v2,
                );
                self.tg12d = GpVec2d::new(secmember.value(1), secmember.value(2));
                self.tg22d = GpVec2d::new(secmember.value(3), secmember.value(4));
                self.istangent = false;
            } else {
                self.istangent = true;
            }

            self.gen.distmin = self.gen.distmin.min(self.pts1.distance(&self.pts2));

            return true;
        }

        false
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:228-274`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn compute_tangent(
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
        let (_, du1, dv1) = self.gen.surf1.d1(u1, v1);
        *norm_first = du1.crossed(&dv1);

        let (_, du2, dv2) = self.gen.surf2.d1(u2, v2);
        *norm_last = du2.crossed(&dv2);

        *tg_first = self.nplan.crossed(norm_first).normalized();
        *tg_last = self.nplan.crossed(norm_last).normalized();

        let mut rev_first = false;
        let mut rev_last = false;

        let choix = self.gen.choix;
        if choix == 2 || choix == 5 {
            rev_first = true;
            rev_last = true;
        }
        if choix == 4 || choix == 7 {
            rev_last = true;
        }
        if choix == 3 || choix == 8 {
            rev_first = true;
        }

        if rev_first {
            tg_first.reverse();
        }
        if rev_last {
            tg_last.reverse();
        }
    }
}

impl MathFunctionSetWithDerivatives for ConstThroat {
    /// `NbVariables()` = `4` (`Blend_Function.cxx:18-21`); not overridden by
    /// this class.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`GenChamfer.cxx:44-47`).
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:114-133`): always `true`.
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f);
        true
    }

    /// `Derivatives(X, D)` (`cxx:137-159`): always `true`.
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d);
        true
    }

    /// `Values(X, F, D)` (`GenChamfer.cxx:95-99`): `Value(X, F)` then
    /// `Derivatives(X, D)`, short-circuiting on the former.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let val = self.value(x, f);
        val && self.derivatives(x, d)
    }
}

impl BlendAppFunction for ConstThroat {
    /// `Set(Param)` (`cxx:48-56`).
    fn set(&mut self, param: f64) {
        self.set_param(param);
    }

    /// `Set(First, Last)` (`GenChamfer.cxx:50`): empty.
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// `GetTolerance(Tolerance, Tol)` (`GenChamfer.cxx:54-61`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.gen.fill_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`GenChamfer.cxx:65-84`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.gen.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:60-108`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`GenChamfer.cxx:88-91`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Ok(self.gen.minimal_distance())
    }

    /// `Pnt1()` = `PointOnS1()` (`Blend_Function.cxx:25-28`).
    fn pnt1(&self) -> &GpPnt {
        &self.pts1
    }

    /// `Pnt2()` = `PointOnS2()` (`Blend_Function.cxx:30-33`).
    fn pnt2(&self) -> &GpPnt {
        &self.pts2
    }

    /// `IsRational()` (`GenChamfer.cxx:125-128`).
    fn is_rational(&self) -> bool {
        self.gen.is_rational()
    }

    /// `GetSectionSize()` throws
    /// `Standard_NotImplemented("BlendFunc_ConstThroat::GetSectionSize()")`
    /// (`cxx:276-280`).
    fn get_section_size(&self) -> f64 {
        panic!("Standard_NotImplemented: BlendFunc_ConstThroat::GetSectionSize()");
    }

    /// `GetMinimalWeight(Weights)` (`GenChamfer.cxx:132-135`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.gen.fill_minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`GenChamfer.cxx:139-142`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.gen.nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`GenChamfer.cxx:146-149`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.gen.guide_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`GenChamfer.cxx:153-159`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        self.gen
            .section_shape(nb_poles, nb_knots, degree, nb_poles2d);
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)`
    /// (`GenChamfer.cxx:165-172`).
    fn get_tolerance_3d(
        &self,
        bound_tol: f64,
        _surf_tol: f64,
        _angle_tol: f64,
        tol3d: &mut MathVector,
        _tol1d: &mut MathVector,
    ) {
        self.gen.section_tolerance_3d(bound_tol, tol3d);
    }

    /// `Knots(TKnots)` (`GenChamfer.cxx:176-180`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.gen.section_knots(tknots);
    }

    /// `Mults(TMults)` (`GenChamfer.cxx:184-188`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.gen.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`GenChamfer.cxx:205-254`).
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
    ) -> bool {
        let tol = self.gen.tol;
        shared::gen_section_with_derivatives(
            self, tol, p, poles, dpoles, poles2d, dpoles2d, weights, dweights,
        )
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`GenChamfer.cxx:256-284`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        shared::gen_section_poles(self, p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` returns `false`
    /// (`Blend_Function.cxx:44-58`).
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

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`GenChamfer.cxx:286-300`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.gen.resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendFunction for ConstThroat {
    /// `PointOnS1()` (`cxx:162-165`).
    fn point_on_s1(&self) -> &GpPnt {
        &self.pts1
    }

    /// `PointOnS2()` (`cxx:169-172`).
    fn point_on_s2(&self) -> &GpPnt {
        &self.pts2
    }

    /// `IsTangencyPoint()` (`cxx:176-179`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS1()` (`cxx:183-192`).
    fn tangent_on_s1(&self) -> &GpVec {
        &self.tg1
    }

    /// `Tangent2dOnS1()` (`cxx:205-213`).
    fn tangent2d_on_s1(&self) -> &GpVec2d {
        &self.tg12d
    }

    /// `TangentOnS2()` (`cxx:196-204`).
    fn tangent_on_s2(&self) -> &GpVec {
        &self.tg2
    }

    /// `Tangent2dOnS2()` (`cxx:217-225`).
    fn tangent2d_on_s2(&self) -> &GpVec2d {
        &self.tg22d
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:228-274`).
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
    ) {
        self.compute_tangent(
            u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last,
        );
    }
}
