//! `BlendFunc_Chamfer` — the "ordinary" chamfer: the distances from the
//! guide line (spine) to the two supports are constant.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_Chamfer.hxx:32-93`,
//! `BlendFunc_Chamfer.cxx:28-246` (OCCT 8.0.0).
//!
//! The vector `X` is `(U1, V1, U2, V2)`. The four equations are the two
//! [`super::Corde`] systems, one per support: each point lies in the section
//! plane and at its own distance from the guide-line point.
//!
//! Translation notes:
//! - `Value` / `Derivatives` copy the `2x2` blocks of the two `Corde` systems
//!   and write the cross blocks as `0.` (`cxx:107-129`).
//! - `IsSolution` short-circuits on the first `Corde` (`cxx:65-66`), so the
//!   second one is not evaluated when the first fails.
//! - `Tangent` takes the guide-line frame from `corde1` (`cxx:203-204`).
//! - `GetSectionSize()` throws `Standard_NotImplemented`
//!   (`cxx:242-246`); the port panics with the same message.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{BlendAppFunction, BlendFunction, BlendPoint};

use super::corde::Corde;
use super::gen_chamfer::GenChamfer;

/// `class BlendFunc_Chamfer : public BlendFunc_GenChamfer` (`hxx:32-93`).
pub struct Chamfer {
    /// The `BlendFunc_GenChamfer` base sub-object.
    pub(super) gen: GenChamfer,
    /// `corde1` (`hxx:91`): the support-surface / guide-line system of `surf1`.
    corde1: Corde,
    /// `corde2` (`hxx:92`): the support-surface / guide-line system of `surf2`.
    corde2: Corde,
}

impl Chamfer {
    /// `BlendFunc_Chamfer(S1, S2, CG)` (`cxx:28-36`): the two `Corde` helpers
    /// are built on the same surfaces and guide line.
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        let corde1 = Corde::new(surf1.clone(), curv.clone());
        let corde2 = Corde::new(surf2.clone(), curv.clone());
        Self {
            gen: GenChamfer::new(surf1, surf2, curv),
            corde1,
            corde2,
        }
    }

    /// `Set(Dist1, Dist2, Choix)` (`cxx:40-45`): the distances are used as
    /// given; only the "quadrant" is stored by the base class.
    pub fn set_distances(&mut self, dist1: f64, dist2: f64, choix: i32) {
        self.corde1.set_dist(dist1);
        self.corde2.set_dist(dist2);
        self.gen.choix = choix;
    }

    /// `Set(Param)` (`cxx:48-52`): both `Corde` frames are refreshed (the
    /// guide line is evaluated twice, as in OCCT).
    pub fn set_param(&mut self, param: f64) {
        self.corde1.set_param(param);
        self.corde2.set_param(param);
    }

    /// `Value(X, F)` (`cxx:78-96`): always `true`.
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) {
        let mut x2 = MathVector::new(1, 2);
        let mut f2 = MathVector::new(1, 2);

        x2.set_value(1, x.value(1));
        x2.set_value(2, x.value(2));
        self.corde1.value(&x2, &mut f2);
        f.set_value(1, f2.value(1));
        f.set_value(2, f2.value(2));

        x2.set_value(1, x.value(3));
        x2.set_value(2, x.value(4));
        self.corde2.value(&x2, &mut f2);
        f.set_value(3, f2.value(1));
        f.set_value(4, f2.value(2));
    }

    /// `Derivatives(X, D)` (`cxx:100-131`): always `true`.
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) {
        let mut x2 = MathVector::new(1, 2);
        let mut d2 = MathMatrix::new(1, 2, 1, 2);

        x2.set_value(1, x.value(1));
        x2.set_value(2, x.value(2));
        self.corde1.derivatives(&x2, &mut d2);
        for i in 1..=2 {
            for j in 1..=2 {
                d.set_value(i, j, d2.value(i, j));
                d.set_value(i, j + 2, 0.0);
            }
        }

        x2.set_value(1, x.value(3));
        x2.set_value(2, x.value(4));
        self.corde2.derivatives(&x2, &mut d2);
        for i in 1..=2 {
            for j in 1..=2 {
                d.set_value(i + 2, j + 2, d2.value(i, j));
                d.set_value(i + 2, j, 0.0);
            }
        }
    }

    /// `IsSolution(Sol, Tol)` (`cxx:56-75`).
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut sol1 = MathVector::new(1, 2);
        let mut sol2 = MathVector::new(1, 2);
        sol1.set_value(1, sol.value(1));
        sol1.set_value(2, sol.value(2));
        sol2.set_value(1, sol.value(3));
        sol2.set_value(2, sol.value(4));

        let issol =
            self.corde1.is_solution(&sol1, tol) && self.corde2.is_solution(&sol2, tol);
        self.gen.tol = tol;
        if issol {
            self.gen.distmin = self
                .gen
                .distmin
                .min(self.corde1.point_on_s().distance(self.corde2.point_on_s()));
        }

        issol
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:188-236`).
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
        let nplan = *self.corde1.n_plan();

        let (_, du1, dv1) = self.gen.surf1.d1(u1, v1);
        *norm_first = du1.crossed(&dv1);

        let (_, du2, dv2) = self.gen.surf2.d1(u2, v2);
        *norm_last = du2.crossed(&dv2);

        *tg_first = nplan.crossed(norm_first).normalized();
        *tg_last = nplan.crossed(norm_last).normalized();

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

impl MathFunctionSetWithDerivatives for Chamfer {
    /// `NbVariables()` = `4`, not overridden.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` = `4` (`GenChamfer.cxx:44-47`).
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:78-96`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f);
        true
    }

    /// `Derivatives(X, D)` (`cxx:100-131`).
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

impl BlendAppFunction for Chamfer {
    /// `Set(Param)` (`cxx:48-52`).
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

    /// `IsSolution(Sol, Tol)` (`cxx:56-75`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`GenChamfer.cxx:88-91`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Ok(self.gen.minimal_distance())
    }

    /// `Pnt1()` = `PointOnS1()` (`Blend_Function.cxx:25-28`).
    fn pnt1(&self) -> &GpPnt {
        self.corde1.point_on_s()
    }

    /// `Pnt2()` = `PointOnS2()` (`Blend_Function.cxx:30-33`).
    fn pnt2(&self) -> &GpPnt {
        self.corde2.point_on_s()
    }

    /// `IsRational()` (`GenChamfer.cxx:125-128`).
    fn is_rational(&self) -> bool {
        self.gen.is_rational()
    }

    /// `GetSectionSize()` throws
    /// `Standard_NotImplemented("BlendFunc_Chamfer::GetSectionSize()")`
    /// (`cxx:242-246`).
    fn get_section_size(&self) -> f64 {
        panic!("Standard_NotImplemented: BlendFunc_Chamfer::GetSectionSize()");
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
        self.gen.resolution(ic2d, tol, tol_u, tol_v);
    }
}

impl BlendFunction for Chamfer {
    /// `PointOnS1()` (`cxx:133-137`).
    fn point_on_s1(&self) -> &GpPnt {
        self.corde1.point_on_s()
    }

    /// `PointOnS2()` (`cxx:141-145`).
    fn point_on_s2(&self) -> &GpPnt {
        self.corde2.point_on_s()
    }

    /// `IsTangencyPoint()` (`cxx:149-152`): both `Corde` must be tangent.
    fn is_tangency_point(&self) -> bool {
        self.corde1.is_tangency_point() && self.corde2.is_tangency_point()
    }

    /// `TangentOnS1()` (`cxx:156-159`).
    fn tangent_on_s1(&self) -> &GpVec {
        self.corde1.tangent_on_s()
    }

    /// `Tangent2dOnS1()` (`cxx:170-173`).
    fn tangent2d_on_s1(&self) -> &GpVec2d {
        self.corde1.tangent2d_on_s()
    }

    /// `TangentOnS2()` (`cxx:163-166`).
    fn tangent_on_s2(&self) -> &GpVec {
        self.corde2.tangent_on_s()
    }

    /// `Tangent2dOnS2()` (`cxx:177-180`).
    fn tangent2d_on_s2(&self) -> &GpVec2d {
        self.corde2.tangent2d_on_s()
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:188-236`).
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
