//! `IsSolution` / `Tangent` / `Resolution` of `BlendFunc_ConstRad`, plus the
//! parametric tolerance and the parameter box.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRad.cxx:797-828`
//! (`GetTolerance`, `GetBounds`), `:832-962` (`IsSolution`), `:1064-1114`
//! (`Tangent`), `:1942-1956` (`Resolution`) (OCCT 8.0.0).

use occt_core::gp::{GpPnt, GpVec, GpVec2d};
use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_math::{svd_padded, svd_solve_occt, Gauss, MathMatrix, MathVector};

use super::constants::{EPS, GAUSS_SOLUTION_EPS, SVD_EPS};
use super::ConstRad;

/// `A * x` for a `math_Matrix` and a `math_Vector`. OCCT spells this
/// `math_Vector operator*(const math_Matrix&, const math_Vector&)`
/// (`math_Matrix.hxx`); the port has no such operator, so it is written out.
pub(super) fn mat_vec(a: &MathMatrix, x: &MathVector) -> MathVector {
    let mut r = MathVector::new(1, a.row_count());
    for i in 1..=a.row_count() {
        let mut sum = 0.0;
        for j in 1..=a.col_count() {
            sum += a.value(i, j) * x.value(j);
        }
        r.set_value(i, sum);
    }
    r
}

impl ConstRad {
    /// `GetTolerance(Tolerance, Tol)` (`cxx:797-804`).
    pub(super) fn fill_parametric_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, approx_same_parameter::u_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(2, approx_same_parameter::v_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(3, approx_same_parameter::u_resolution(self.surf2.as_ref(), tol));
        tolerance.set_value(4, approx_same_parameter::v_resolution(self.surf2.as_ref(), tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:808-828`): the parameter boxes of
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

    /// `IsSolution(Sol, Tol)` (`cxx:832-962`): the four values are within
    /// tolerance, then the tangents are solved for from the jacobian
    /// (`math_Gauss` with `1.e-9`, falling back to `math_SVD` with `1.e-6`),
    /// and `maxang` / `minang` / `distmin` are updated from the section angle.
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        // `Ok = ComputeValues(Sol, 1, true, param)` (`cxx:838`); the result is
        // what the successful branch returns (`cxx:959`).
        let ok = self.compute_values(sol, 1, true, self.param);

        if !(self.e.value(1).abs() <= tol
            && self.e.value(2) * self.e.value(2)
                + self.e.value(3) * self.e.value(3)
                + self.e.value(4) * self.e.value(4)
                <= tol * tol)
        {
            self.istangent = true;
            return false;
        }

        // ns1, ns2 and np are copied locally to avoid crushing the fields
        // (`cxx:843-846`).
        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let np = self.nplan;

        let mut norm = np.crossed(&ns1).magnitude();
        if norm < EPS {
            // Unsatisfactory, but it is not necessary to stop (`cxx:851`).
            norm = 1.0;
        }
        ns1 = GpVec::zero().set_linear_form_2(np.dot(&ns1) / norm, &np, -1.0 / norm, &ns1);

        norm = np.crossed(&ns2).magnitude();
        if norm < EPS {
            norm = 1.0;
        }
        ns2 = GpVec::zero().set_linear_form_2(np.dot(&ns2) / norm, &np, -1.0 / norm, &ns2);

        let mut tolerances = MathVector::new(1, 4);
        self.fill_parametric_tolerance(&mut tolerances, tol);

        let mut solution = MathVector::new(1, 4);
        self.istangent = true;

        let resol = Gauss::new(&self.dedx, GAUSS_SOLUTION_EPS);
        if resol.is_done() {
            let minus_dedt = self.dedt.multiplied_scalar(-1.0);
            // `Resol.Solve(-DEDT, solution)` cannot fail once `IsDone()`
            // (`cxx:870`).
            if let Some(x) = resol.solve(&minus_dedt) {
                solution = x;
            }
            self.istangent = false;
            let controle = self.dedt.added(&mat_vec(&self.dedx, &solution));
            if controle.value(1).abs() > tolerances.value(1)
                || controle.value(2).abs() > tolerances.value(2)
                || controle.value(3).abs() > tolerances.value(3)
                || controle.value(4).abs() > tolerances.value(4)
            {
                self.istangent = true;
            }
        }

        if self.istangent {
            let sing_rs = svd_padded(&self.dedx);
            if sing_rs.is_done {
                let minus_dedt = self.dedt.multiplied_scalar(-1.0);
                svd_solve_occt(&sing_rs, &minus_dedt, &mut solution, SVD_EPS);
                self.istangent = false;
                let controle = self.dedt.added(&mat_vec(&self.dedx, &solution));
                if controle.value(1).abs() > tolerances.value(1)
                    || controle.value(2).abs() > tolerances.value(2)
                    || controle.value(3).abs() > tolerances.value(3)
                    || controle.value(4).abs() > tolerances.value(4)
                {
                    // `std::cout << "Cheminement : echec calcul des derivees"`
                    // (`cxx:891`).
                    self.istangent = true;
                }
            }
        }

        if !self.istangent {
            self.tg1 = GpVec::zero().set_linear_form_2(
                solution.value(1),
                &self.d1u1,
                solution.value(2),
                &self.d1v1,
            );
            self.tg2 = GpVec::zero().set_linear_form_2(
                solution.value(3),
                &self.d1u2,
                solution.value(4),
                &self.d1v2,
            );
            self.tg12d = GpVec2d::new(solution.value(1), solution.value(2));
            self.tg22d = GpVec2d::new(solution.value(3), solution.value(4));
        }

        // Update of maxang (`cxx:907-955`).
        if self.ray1 > 0.0 {
            ns1.reverse();
        }
        if self.ray2 > 0.0 {
            ns2.reverse();
        }
        let mut cosa = ns1.dot(&ns2);
        let mut sina = np.dot(&ns1.crossed(&ns2));
        if self.choix % 2 != 0 {
            // nplan is changed in -nplan (`cxx:921`).
            sina = -sina;
        }
        if cosa > 1.0 {
            cosa = 1.0;
            sina = 0.0;
        }
        let mut angle = cosa.acos();
        // Reframing on ]-pi/2, 3pi/2] (`cxx:931-943`).
        if sina < 0.0 {
            if cosa > 0.0 {
                angle = -angle;
            } else {
                angle = 2.0 * std::f64::consts::PI - angle;
            }
        }

        if angle.abs() > self.maxang {
            self.maxang = angle.abs();
        }
        if angle.abs() < self.minang {
            self.minang = angle.abs();
        }
        self.distmin = self.distmin.min(self.pts1.distance(&self.pts2));

        ok
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:1064-1114`): the
    /// tangents at the two extremities of the section and the normals of the
    /// two surfaces there.
    ///
    /// The surface normals are only re-evaluated when the parameters differ
    /// from the cached `xval` (`cxx:1080-1093`).
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
        let ns1;
        if u1 != self.x_val.value(1)
            || v1 != self.x_val.value(2)
            || u2 != self.x_val.value(3)
            || v2 != self.x_val.value(4)
        {
            let (_, d1u, d1v) = self.surf1.d1(u1, v1);
            ns1 = d1u.crossed(&d1v);
            *norm_first = ns1;
            let (_, d1u, d1v) = self.surf2.d1(u2, v2);
            *norm_last = d1u.crossed(&d1v);
        } else {
            ns1 = self.nsurf1;
            *norm_first = self.nsurf1;
            *norm_last = self.nsurf2;
        }

        let mut invnorm1 = self.nplan.crossed(&ns1).magnitude();
        if invnorm1 < EPS {
            invnorm1 = 1.0;
        } else {
            invnorm1 = 1.0 / invnorm1;
        }

        let ns1 = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns1) * invnorm1,
            &self.nplan,
            -invnorm1,
            &ns1,
        );
        let center = GpPnt::from_xyz(&self.pts1.coord.added(&ns1.coord.multiply(self.ray1)));

        // `gp_Vec(Center, pts1)` / `gp_Vec(Center, pts2)`.
        *tg_first = self.nplan.crossed(&GpVec::from_pnts(&center, &self.pts1));
        *tg_last = self.nplan.crossed(&GpVec::from_pnts(&center, &self.pts2));
        if self.choix % 2 == 1 {
            tg_first.reverse();
            tg_last.reverse();
        }
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1942-1956`).
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
