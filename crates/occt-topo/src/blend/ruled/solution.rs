//! `IsSolution` / `Tangent` / `Resolution` of `BlendFunc_Ruled`, plus the
//! parametric tolerance and the parameter box.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_Ruled.cxx:65-93`
//! (`GetTolerance`, `GetBounds`), `:95-168` (`IsSolution`), `:396-414`
//! (`Tangent`), `:731-744` (`Resolution`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `IsSolution` solves `gradsol * secmember` by `math_Gauss` with the
//!   default `MinPivot` (`cxx:149`) and, unlike `ConstRad` / `EvolRad`, has no
//!   `math_SVD` fallback.
//! - The tangency guard raises `Standard_DomainError` in
//!   `TangentOnS1` / `TangentOnS2` (`cxx:432-449`) and in the 2d accessors
//!   (`cxx:451-466`); see the module note of [`super`].

use occt_core::gp::{GpVec, GpVec2d};
use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_math::{Gauss, MathMatrix, MathVector};

use super::{Ruled, GAUSS_MIN_PIVOT};

impl Ruled {
    /// `GetTolerance(Tolerance, Tol)` (`cxx:65-72`).
    pub(super) fn fill_parametric_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, approx_same_parameter::u_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(2, approx_same_parameter::v_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(3, approx_same_parameter::u_resolution(self.surf2.as_ref(), tol));
        tolerance.set_value(4, approx_same_parameter::v_resolution(self.surf2.as_ref(), tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:74-93`): the parameter boxes of
    /// the two surfaces, each widened by its own range when both ends are
    /// finite.
    pub(super) fn fill_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        let (u1min, u1max) = self.surf1.u_range();
        let (v1min, v1max) = self.surf1.v_range();
        let (u2min, u2max) = self.surf2.u_range();
        let (v2min, v2max) = self.surf2.v_range();

        inf_bound.set_value(1, u1min);
        inf_bound.set_value(2, v1min);
        inf_bound.set_value(3, u2min);
        inf_bound.set_value(4, v2min);
        sup_bound.set_value(1, u1max);
        sup_bound.set_value(2, v1max);
        sup_bound.set_value(3, u2max);
        sup_bound.set_value(4, v2max);

        for i in 1..=4 {
            if !Precision::is_infinite(inf_bound.value(i))
                && !Precision::is_infinite(sup_bound.value(i))
            {
                let range = sup_bound.value(i) - inf_bound.value(i);
                inf_bound.set_value(i, inf_bound.value(i) - range);
                sup_bound.set_value(i, sup_bound.value(i) + range);
            }
        }
    }

    /// `IsSolution(Sol, Tol)` (`cxx:95-168`): all four residuals are within
    /// `Tol`, then the tangents are solved for from the jacobian
    /// (`math_Gauss`, default pivot). Returns `true` even when the linear
    /// system is singular — `istangent` reports that case.
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 4);
        let mut gradsol = MathMatrix::new(1, 4, 1, 4);

        self.compute_values(sol, Some(&mut valsol), &mut gradsol);

        if !(valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= tol
            && valsol.value(4).abs() <= tol)
        {
            self.istangent = true;
            return false;
        }

        // `pts1` / `pts2` are already the ones of `Values` (`cxx:110-111`);
        // only the first derivatives are needed here.
        let (_pts1, d1u1, d1v1) = self.surf1.d1(sol.value(1), sol.value(2));
        let (_pts2, d1u2, d1v2) = self.surf2.d1(sol.value(3), sol.value(4));

        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let mut secmember = MathVector::new(1, 4);
        secmember.set_value(
            1,
            self.normtg - dnplan.dot(&GpVec::from_pnts(&self.ptgui, &self.pts1)),
        );
        secmember.set_value(
            2,
            self.normtg - dnplan.dot(&GpVec::from_pnts(&self.ptgui, &self.pts2)),
        );

        let p1p2 = GpVec::from_pnts(&self.pts1, &self.pts2);

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        // Derivative of nor1 with respect to the guide parameter
        // (`cxx:123-129`).
        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let temp = GpVec::zero().set_linear_form_3(
            (dnplan.dot(&ns) - grosterme * ndotns) / norm,
            &self.nplan,
            ndotns / norm,
            &dnplan,
            grosterme / norm,
            &ns,
        );

        secmember.set_value(3, -temp.dot(&p1p2));

        let ns = d1u2.crossed(&d1v2);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        // Derivative of nor2 with respect to the guide parameter
        // (`cxx:138-144`).
        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let temp = GpVec::zero().set_linear_form_3(
            (dnplan.dot(&ns) - grosterme * ndotns) / norm,
            &self.nplan,
            ndotns / norm,
            &dnplan,
            grosterme / norm,
            &ns,
        );

        secmember.set_value(4, -temp.dot(&p1p2));

        let resol = Gauss::new(&gradsol, GAUSS_MIN_PIVOT);
        if resol.is_done() {
            if let Some(x) = resol.solve(&secmember) {
                secmember = x;
            }

            self.tg1 = GpVec::zero().set_linear_form_2(
                secmember.value(1),
                &d1u1,
                secmember.value(2),
                &d1v1,
            );
            self.tg2 = GpVec::zero().set_linear_form_2(
                secmember.value(3),
                &d1u2,
                secmember.value(4),
                &d1v2,
            );
            self.tg12d = GpVec2d::new(secmember.value(1), secmember.value(2));
            self.tg22d = GpVec2d::new(secmember.value(3), secmember.value(4));
            self.istangent = false;
        } else {
            self.istangent = true;
        }
        true
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:396-414`): the
    /// tangents are both `Vec(pts1, pts2)`, the normals are the surface
    /// normals `d1u ^ d1v`.
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
        let (_bid, d1u, d1v) = self.surf2.d1(u2, v2);
        *norm_last = d1u.crossed(&d1v);

        let (_bid, d1u, d1v) = self.surf1.d1(u1, v1);
        *norm_first = d1u.crossed(&d1v);

        let tg = GpVec::from_pnts(&self.pts1, &self.pts2);
        *tg_first = tg;
        *tg_last = tg;
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:731-744`).
    pub(super) fn compute_resolution(
        &self,
        ic2d: i32,
        tol: f64,
        tol_u: &mut f64,
        tol_v: &mut f64,
    ) {
        if ic2d == 1 {
            *tol_u = approx_same_parameter::u_resolution(self.surf1.as_ref(), tol);
            *tol_v = approx_same_parameter::v_resolution(self.surf1.as_ref(), tol);
        } else {
            *tol_u = approx_same_parameter::u_resolution(self.surf2.as_ref(), tol);
            *tol_v = approx_same_parameter::v_resolution(self.surf2.as_ref(), tol);
        }
    }
}
