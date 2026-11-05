//! `BlendFunc_ConstThroatWithPenetration` — the chamfer whose section is a
//! right-angled triangle with a constant leg, so the two faces penetrate.
//!
//! Source:
//! `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstThroatWithPenetration.hxx:28-100`,
//! `BlendFunc_ConstThroatWithPenetration.cxx:27-185` (OCCT 8.0.0).
//!
//! Same `X = (U1, V1, U2, V2)` and section plane as
//! [`super::ConstThroat`], with the last two equations replaced by
//! 3. `|pts1 - ptgui|^2 - Throat^2`,
//! 4. `(pts1 - ptgui) . (pts2 - pts1)` — the section is right-angled at the
//!    first point.
//!
//! Translation notes:
//! - Only `Value`, `Derivatives`, `IsSolution`, the four tangent getters and
//!   `GetSectionSize` are overridden by this class; everything else is the
//!   inherited `ConstThroat` body and is delegated to `parent`.
//! - The inherited `Tangent(U1, V1, U2, V2, ...)` (with its `choix` sign
//!   rules) is *not* overridden.
//! - The tangent getters throw `Standard_DomainError` when `istangent`
//!   (`cxx:138-177`); the trait returns a value and every OCCT caller guards
//!   the call, so the stored value is returned.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{Gauss, MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{BlendAppFunction, BlendFunction, BlendPoint};

use super::const_throat::ConstThroat;
use super::shared::{GAUSS_MIN_PIVOT, pnt_vec};

/// `class BlendFunc_ConstThroatWithPenetration : public BlendFunc_ConstThroat`
/// (`hxx:28-100`).
pub struct ConstThroatWithPenetration {
    /// The `BlendFunc_ConstThroat` base sub-object.
    pub(super) parent: ConstThroat,
}

impl ConstThroatWithPenetration {
    /// `BlendFunc_ConstThroatWithPenetration(S1, S2, C)` (`cxx:27-34`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            parent: ConstThroat::new(surf1, surf2, curv),
        }
    }

    /// `Set(aThroat, dis1, Choix)` (`ConstThroat.cxx:42-46`).
    pub fn set_throat(&mut self, a_throat: f64, dis1: f64, choix: i32) {
        self.parent.set_throat(a_throat, dis1, choix);
    }

    /// `Set(Param)` (`ConstThroat.cxx:48-56`).
    pub fn set_param(&mut self, param: f64) {
        self.parent.set_param(param);
    }

    /// `Throat` (`ConstThroat.hxx:100`).
    pub fn throat(&self) -> f64 {
        self.parent.throat()
    }

    /// `Value(X, F)` (`cxx:91-107`).
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) {
        self.parent.pts1 = self.parent.gen.surf1.value(x.value(1), x.value(2));
        self.parent.pts2 = self.parent.gen.surf2.value(x.value(3), x.value(4));

        f.set_value(
            1,
            self.parent.nplan.dot(&pnt_vec(&self.parent.pts1)) + self.parent.the_d,
        );
        f.set_value(
            2,
            self.parent.nplan.dot(&pnt_vec(&self.parent.pts2)) + self.parent.the_d,
        );

        let vref = GpVec::from_pnts(&self.parent.ptgui, &self.parent.pts1);
        f.set_value(3, vref.square_magnitude() - self.parent.throat * self.parent.throat);

        let vec12 = GpVec::from_pnts(&self.parent.pts1, &self.parent.pts2);
        f.set_value(4, vref.dot(&vec12));
    }

    /// `Derivatives(X, D)` (`cxx:112-135`).
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) {
        let (p1, u1, v1) = self.parent.gen.surf1.d1(x.value(1), x.value(2));
        self.parent.pts1 = p1;
        self.parent.d1u1 = u1;
        self.parent.d1v1 = v1;
        let (p2, u2, v2) = self.parent.gen.surf2.d1(x.value(3), x.value(4));
        self.parent.pts2 = p2;
        self.parent.d1u2 = u2;
        self.parent.d1v2 = v2;

        let vref = GpVec::from_pnts(&self.parent.ptgui, &self.parent.pts1);
        let vec12 = GpVec::from_pnts(&self.parent.pts1, &self.parent.pts2);

        d.set_value(1, 1, self.parent.nplan.dot(&self.parent.d1u1));
        d.set_value(1, 2, self.parent.nplan.dot(&self.parent.d1v1));
        d.set_value(1, 3, 0.0);
        d.set_value(1, 4, 0.0);
        d.set_value(2, 1, 0.0);
        d.set_value(2, 2, 0.0);
        d.set_value(2, 3, self.parent.nplan.dot(&self.parent.d1u2));
        d.set_value(2, 4, self.parent.nplan.dot(&self.parent.d1v2));
        d.set_value(3, 1, 2.0 * vref.dot(&self.parent.d1u1));
        d.set_value(3, 2, 2.0 * vref.dot(&self.parent.d1v1));
        d.set_value(3, 3, 0.0);
        d.set_value(3, 4, 0.0);
        d.set_value(4, 1, self.parent.d1u1.dot(&vec12) - vref.dot(&self.parent.d1u1));
        d.set_value(4, 2, self.parent.d1v1.dot(&vec12) - vref.dot(&self.parent.d1v1));
        d.set_value(4, 3, vref.dot(&self.parent.d1u2));
        d.set_value(4, 4, vref.dot(&self.parent.d1v2));
    }

    /// `IsSolution(Sol, Tol)` (`cxx:36-85`): the last value is only compared
    /// with `Tol` (not `Tol^2` as in `ConstThroat`).
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut secmember = MathVector::new(1, 4);
        let mut valsol = MathVector::new(1, 4);
        let mut gradsol = MathMatrix::new(1, 4, 1, 4);

        self.compute_value(sol, &mut valsol);
        self.compute_derivatives(sol, &mut gradsol);

        self.parent.gen.tol = tol;

        if valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= tol * tol
            && valsol.value(4).abs() <= tol
        {
            let dnplan = GpVec::zero().set_linear_form_2(
                1.0 / self.parent.normtg,
                &self.parent.d2gui,
                -1.0 / self.parent.normtg * self.parent.nplan.dot(&self.parent.d2gui),
                &self.parent.nplan,
            );

            let temp1 = GpVec::from_pnts(&self.parent.ptgui, &self.parent.pts1);
            let temp2 = GpVec::from_pnts(&self.parent.ptgui, &self.parent.pts2);
            let temp3 = GpVec::from_pnts(&self.parent.pts1, &self.parent.pts2);

            let ng = self.parent.nplan.dot(&self.parent.d1gui);
            secmember.set_value(1, ng - dnplan.dot(&temp1));
            secmember.set_value(2, ng - dnplan.dot(&temp2));
            secmember.set_value(3, 2.0 * self.parent.d1gui.dot(&temp1));
            secmember.set_value(4, self.parent.d1gui.dot(&temp3));

            let (p1, u1, v1) = self.parent.gen.surf1.d1(sol.value(1), sol.value(2));
            self.parent.pts1 = p1;
            self.parent.d1u1 = u1;
            self.parent.d1v1 = v1;
            let (p2, u2, v2) = self.parent.gen.surf2.d1(sol.value(3), sol.value(4));
            self.parent.pts2 = p2;
            self.parent.d1u2 = u2;
            self.parent.d1v2 = v2;

            let resol = Gauss::new(&gradsol, GAUSS_MIN_PIVOT);
            if resol.is_done() {
                if let Some(x) = resol.solve(&secmember) {
                    secmember = x;
                }
                self.parent.tg1 = GpVec::zero().set_linear_form_2(
                    secmember.value(1),
                    &self.parent.d1u1,
                    secmember.value(2),
                    &self.parent.d1v1,
                );
                self.parent.tg2 = GpVec::zero().set_linear_form_2(
                    secmember.value(3),
                    &self.parent.d1u2,
                    secmember.value(4),
                    &self.parent.d1v2,
                );
                self.parent.tg12d = GpVec2d::new(secmember.value(1), secmember.value(2));
                self.parent.tg22d = GpVec2d::new(secmember.value(3), secmember.value(4));
                self.parent.istangent = false;
            } else {
                self.parent.istangent = true;
            }

            self.parent.gen.distmin = self
                .parent
                .gen
                .distmin
                .min(self.parent.pts1.distance(&self.parent.pts2));

            return true;
        }

        false
    }
}

impl MathFunctionSetWithDerivatives for ConstThroatWithPenetration {
    /// `NbVariables()` = `4`, not overridden.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` = `4`, not overridden.
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:91-107`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f);
        true
    }

    /// `Derivatives(X, D)` (`cxx:112-135`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d);
        true
    }

    /// `Values(X, F, D)` (`GenChamfer.cxx:95-99`), not overridden.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let val = self.value(x, f);
        val && self.derivatives(x, d)
    }
}

impl BlendAppFunction for ConstThroatWithPenetration {
    /// `Set(Param)` (`ConstThroat.cxx:48-56`), not overridden.
    fn set(&mut self, param: f64) {
        self.parent.set_param(param);
    }

    /// `Set(First, Last)` (`GenChamfer.cxx:50`): empty.
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// `GetTolerance(Tolerance, Tol)` (`GenChamfer.cxx:54-61`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.parent.gen.fill_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`GenChamfer.cxx:65-84`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.parent.gen.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:36-85`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`GenChamfer.cxx:88-91`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Ok(self.parent.gen.minimal_distance())
    }

    /// `Pnt1()` = `PointOnS1()`.
    fn pnt1(&self) -> &GpPnt {
        &self.parent.pts1
    }

    /// `Pnt2()` = `PointOnS2()`.
    fn pnt2(&self) -> &GpPnt {
        &self.parent.pts2
    }

    /// `IsRational()` (`GenChamfer.cxx:125-128`).
    fn is_rational(&self) -> bool {
        self.parent.gen.is_rational()
    }

    /// `GetSectionSize()` throws
    /// `Standard_NotImplemented("BlendFunc_ConstThroatWithPenetration::GetSectionSize()")`
    /// (`cxx:182-185`).
    fn get_section_size(&self) -> f64 {
        panic!(
            "Standard_NotImplemented: BlendFunc_ConstThroatWithPenetration::GetSectionSize()"
        );
    }

    /// `GetMinimalWeight(Weights)` (`GenChamfer.cxx:132-135`).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        self.parent.gen.fill_minimal_weights(weights);
    }

    /// `NbIntervals(S)` (`GenChamfer.cxx:139-142`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.parent.gen.nb_intervals(s) as usize
    }

    /// `Intervals(T, S)` (`GenChamfer.cxx:146-149`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        self.parent.gen.guide_intervals(t, s);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`GenChamfer.cxx:153-159`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        self.parent
            .gen
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
        self.parent.gen.section_tolerance_3d(bound_tol, tol3d);
    }

    /// `Knots(TKnots)` (`GenChamfer.cxx:176-180`).
    fn knots(&mut self, tknots: &mut [f64]) {
        self.parent.gen.section_knots(tknots);
    }

    /// `Mults(TMults)` (`GenChamfer.cxx:184-188`).
    fn mults(&mut self, tmults: &mut [i32]) {
        self.parent.gen.section_mults(tmults);
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`GenChamfer.cxx:205-254`), with the overridden `Value` / `IsSolution`
    /// reached through the traits.
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
        let tol = self.parent.gen.tol;
        super::shared::gen_section_with_derivatives(
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
        super::shared::gen_section_poles(self, p, poles, poles2d, weights);
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
        self.parent.gen.resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendFunction for ConstThroatWithPenetration {
    /// `PointOnS1()` (`ConstThroat.cxx:162-165`), not overridden.
    fn point_on_s1(&self) -> &GpPnt {
        &self.parent.pts1
    }

    /// `PointOnS2()` (`ConstThroat.cxx:169-172`), not overridden.
    fn point_on_s2(&self) -> &GpPnt {
        &self.parent.pts2
    }

    /// `IsTangencyPoint()` (`ConstThroat.cxx:176-179`), not overridden.
    fn is_tangency_point(&self) -> bool {
        self.parent.istangent
    }

    /// `TangentOnS1()` (`cxx:138-146`).
    fn tangent_on_s1(&self) -> &GpVec {
        &self.parent.tg1
    }

    /// `Tangent2dOnS1()` (`cxx:160-168`).
    fn tangent2d_on_s1(&self) -> &GpVec2d {
        &self.parent.tg12d
    }

    /// `TangentOnS2()` (`cxx:148-158`).
    fn tangent_on_s2(&self) -> &GpVec {
        &self.parent.tg2
    }

    /// `Tangent2dOnS2()` (`cxx:170-177`).
    fn tangent2d_on_s2(&self) -> &GpVec2d {
        &self.parent.tg22d
    }

    /// `Tangent(...)` (`ConstThroat.cxx:228-274`), not overridden.
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
        self.parent
            .compute_tangent(u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last);
    }
}
