//! `IsSolution` / `Tangent` / `Resolution` of `BlendFunc_CSConstRad`, plus the
//! parametric tolerance and the parameter box.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSConstRad.cxx:101-133`
//! (`GetTolerance`, `GetBounds`), `:135-224` (`IsSolution`), `:447-463`
//! (`Tangent`), `:983-988` (`Resolution`) (OCCT 8.0.0).

use std::f64::consts::PI;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_math::{Gauss, MathVector};

use super::constants::GAUSS_MIN_PIVOT;
use super::CsConstRad;

impl CsConstRad {
    /// `GetTolerance(Tolerance, Tol)` (`cxx:101-106`).
    pub(super) fn fill_parametric_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, approx_same_parameter::u_resolution(self.surf.as_ref(), tol));
        tolerance.set_value(2, approx_same_parameter::v_resolution(self.surf.as_ref(), tol));
        tolerance.set_value(3, self.curv.resolution(tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:111-133`): the parameter box of the
    /// surface and the parameter range of the curve. Only the two surface
    /// parameters are widened by their own range (`cxx:119-132`).
    pub(super) fn fill_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        let (umin, umax) = self.surf.u_range();
        let (vmin, vmax) = self.surf.v_range();

        inf_bound.set_value(1, umin);
        inf_bound.set_value(2, vmin);
        inf_bound.set_value(3, self.curv.first_parameter());
        sup_bound.set_value(1, umax);
        sup_bound.set_value(2, vmax);
        sup_bound.set_value(3, self.curv.last_parameter());

        for i in 1..=2 {
            if !Precision::is_infinite(inf_bound.value(i))
                && !Precision::is_infinite(sup_bound.value(i))
            {
                let range = sup_bound.value(i) - inf_bound.value(i);
                inf_bound.set_value(i, inf_bound.value(i) - range);
                sup_bound.set_value(i, sup_bound.value(i) + range);
            }
        }
    }

    /// `IsSolution(Sol, Tol)` (`cxx:135-224`): the three values are within
    /// tolerance, then the tangents are solved for from the jacobian
    /// (`math_Gauss` with the default `1.e-20` pivot, no SVD fallback), and
    /// `maxang` / `minang` are updated from the section angle.
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        // `Values(Sol, valsol, gradsol)` (`cxx:143`).
        self.evaluate_values(sol);

        if self.e.value(1).abs() > tol
            || self.e.value(2).abs() > tol
            || self.e.value(3).abs() > tol * tol
        {
            self.istangent = true;
            return false;
        }

        // Calculation of the tangents (`cxx:149-160`).
        self.pt2d = GpPnt2d::new(sol.value(1), sol.value(2));
        self.prmc = sol.value(3);
        let (pts, d1u1, d1v1) = self.surf.d1(sol.value(1), sol.value(2));
        self.pts = pts;
        let (ptc, d1) = self.curv.d1(sol.value(3));
        self.ptc = ptc;

        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let mut secmember = MathVector::new(1, 3);
        let temp = GpVec::from_pnts(&self.ptgui, &self.pts);
        secmember.set_value(1, self.normtg - dnplan.dot(&temp));

        let temp = GpVec::from_pnts(&self.ptgui, &self.ptc);
        secmember.set_value(2, self.normtg - dnplan.dot(&temp));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let temp = GpVec::zero().set_linear_form_3(
            self.ray / norm * (dnplan.dot(&ns) - grosterme * ndotns),
            &self.nplan,
            self.ray * ndotns / norm,
            &dnplan,
            self.ray * grosterme / norm,
            &ns,
        );

        let ns = GpVec::zero().set_linear_form_2(
            ndotns / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );
        let resul = ns
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));
        secmember.set_value(3, -2.0 * temp.dot(&resul));

        self.istangent = true;
        let resol = Gauss::new(&self.dedx, GAUSS_MIN_PIVOT);
        if resol.is_done() {
            if let Some(x) = resol.solve(&secmember) {
                secmember = x;
            }
            self.tgs = GpVec::zero().set_linear_form_2(
                secmember.value(1),
                &d1u1,
                secmember.value(2),
                &d1v1,
            );
            self.tgc = d1.multiplied_scalar(secmember.value(3));
            self.tg2d = GpVec2d::new(secmember.value(1), secmember.value(2));
            self.istangent = false;
        }

        // Update of maxang (`cxx:194-218`).
        let ns2 = resul.normalized().reversed();
        let cosa = ns.dot(&ns2);
        let mut sina = self.nplan.dot(&ns.crossed(&ns2));
        if self.choix % 2 != 0 {
            // nplan is changed in -nplan (`cxx:200-203`).
            sina = -sina;
        }

        let mut angle = cosa.acos();
        if sina < 0.0 {
            angle = 2.0 * PI - angle;
        }

        if angle > self.maxang {
            self.maxang = angle;
        }
        if angle < self.minang {
            self.minang = angle;
        }

        true
    }

    /// `Tangent(U, V, TgS, NmS)` (`cxx:447-463`): the tangent at the section and
    /// the normal of the surface at `(U, V)`.
    pub(super) fn compute_tangent(&self, u: f64, v: f64, tg_s: &mut GpVec, norm_s: &mut GpVec) {
        let (bid, d1u, d1v) = self.surf.d1(u, v);

        let ns = d1u.crossed(&d1v);
        *norm_s = ns;

        let norm = self.nplan.crossed(&ns).magnitude();
        let ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );

        let center = GpPnt::from_xyz(&bid.coord.added(&ns.coord.multiply(self.ray)));
        *tg_s = self.nplan.crossed(&GpVec::from_pnts(&center, &bid));
        if self.choix % 2 == 1 {
            tg_s.reverse();
        }
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:983-988`).
    pub(super) fn compute_resolution(
        &self,
        _ic2d: i32,
        tol: f64,
        tol_u: &mut f64,
        tol_v: &mut f64,
    ) {
        *tol_u = approx_same_parameter::u_resolution(self.surf.as_ref(), tol);
        *tol_v = approx_same_parameter::v_resolution(self.surf.as_ref(), tol);
    }
}
