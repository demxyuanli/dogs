//! The second-order section of `BlendFunc_ConstRad`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRad.cxx:1527-1905`
//! (OCCT 8.0.0). This is the `Section(P, Poles, DPoles, D2Poles, Poles2d,
//! DPoles2d, D2Poles2d, Weights, DWeights, D2Weights)` overload used for the
//! first and last section of the approximation.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::{svd_padded, svd_solve_occt, Gauss, MathMatrix, MathVector};

use occt_geom::geom_fill::section_convert;

use crate::blend::{BlendFuncSectionShape, BlendPoint};

use super::constants::{EPS, GAUSS_SOLUTION_EPS, POINTS_EQUAL_TOL, SVD_EPS};
use super::solution::mat_vec;
use super::ConstRad;

impl ConstRad {
    /// `Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weights, DWeights, D2Weights)` (`cxx:1527-1905`): the first/last section
    /// with first and second derivatives. Returns whether the derivatives were
    /// computed.
    ///
    /// `ComputeValues(X, 2, true, prm)` (`cxx:1598`) is what fills `D2EDX2`,
    /// `D2EDXDT` and `D2EDT2`; the two linear solves then use the same
    /// `math_Gauss` (`1.e-9`) / `math_SVD` (`1.e-6`) fallback as the
    /// first-order section.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn section_with_second_derivatives(
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
        let prm = p.parameter();
        let low = 0usize;
        let upp = poles.len() - 1;

        let (u1, v1) = p
            .parameters_on_s1()
            .expect("BlendFunc_ConstRad::Section: point is not on the first surface");
        let (u2, v2) = p
            .parameters_on_s2()
            .expect("BlendFunc_ConstRad::Section: point is not on the second surface");
        let mut x = MathVector::new(1, 4);
        x.set_value(1, u1);
        x.set_value(2, v1);
        x.set_value(3, u2);
        x.set_value(4, v2);

        // Calculation of equations (`cxx:1597-1599`).
        self.compute_values(&x, 2, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(&self.pts2));

        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let mut np = self.nplan;
        let mut dnp = self.dnplan;
        let mut d2np = self.d2nplan;

        let mut sol = MathVector::new(1, 4);
        let mut secmember = MathVector::new(1, 4);
        let mut d2dx_dsdt = MathMatrix::new(1, 4, 1, 4);
        let mut istgt = true;

        // Calculation of derivatives (`cxx:1704-1721`).
        if self.pts1.distance(&self.pts2) > POINTS_EQUAL_TOL {
            // Precise tolerance !!!!! (`cxx:1707`).
            let resol = Gauss::new(&self.dedx, GAUSS_SOLUTION_EPS);
            if resol.is_done() {
                let minus_dedt = self.dedt.multiplied_scalar(-1.0);
                if let Some(x) = resol.solve(&minus_dedt) {
                    sol = x;
                }
                self.d2edx2.multiply(&sol, &mut d2dx_dsdt);
                secmember = self.second_member(&d2dx_dsdt, &sol);
                if let Some(x) = resol.solve(&secmember) {
                    secmember = x;
                }
                istgt = false;
            }
        }

        if istgt {
            let sing_rs = svd_padded(&self.dedx);
            if sing_rs.is_done {
                let minus_dedt = self.dedt.multiplied_scalar(-1.0);
                svd_solve_occt(&sing_rs, &minus_dedt, &mut sol, SVD_EPS);
                self.d2edx2.multiply(&sol, &mut d2dx_dsdt);
                let vbis = self.second_member(&d2dx_dsdt, &sol);
                svd_solve_occt(&sing_rs, &vbis, &mut secmember, SVD_EPS);
                istgt = false;
            }
        }

        let mut dnorm1w = GpVec::zero();
        let mut dnorm2w = GpVec::zero();
        let mut d2norm1w = GpVec::zero();
        let mut d2norm2w = GpVec::zero();
        let mut dtg1 = GpVec::zero();
        let mut dtg2 = GpVec::zero();

        if !istgt {
            let s1 = sol.value(1);
            let s2 = sol.value(2);
            let s3 = sol.value(3);
            let s4 = sol.value(4);

            self.tg1 = GpVec::zero().set_linear_form_2(s1, &self.d1u1, s2, &self.d1v1);
            self.tg2 = GpVec::zero().set_linear_form_2(s3, &self.d1u2, s4, &self.d1v2);

            dnorm1w = GpVec::zero().set_linear_form_3(
                s1,
                &self.dndu1,
                s2,
                &self.dndv1,
                1.0,
                &self.dn1w,
            );
            dnorm2w = GpVec::zero().set_linear_form_3(
                s3,
                &self.dndu2,
                s4,
                &self.dndv2,
                1.0,
                &self.dn2w,
            );

            let temp = GpVec::zero().set_linear_form_3(
                s1 * s1,
                &self.d2u1,
                2.0 * s1 * s2,
                &self.d2uv1,
                s2 * s2,
                &self.d2v1,
            );
            dtg1 = GpVec::zero()
                .set_linear_form_3(
                    secmember.value(1),
                    &self.d1u1,
                    secmember.value(2),
                    &self.d1v1,
                    1.0,
                    &temp,
                );

            let temp = GpVec::zero().set_linear_form_3(
                s3 * s3,
                &self.d2u2,
                2.0 * s3 * s4,
                &self.d2uv2,
                s4 * s4,
                &self.d2v2,
            );
            dtg2 = GpVec::zero().set_linear_form_3(
                secmember.value(3),
                &self.d1u2,
                secmember.value(4),
                &self.d1v2,
                1.0,
                &temp,
            );

            let temp = GpVec::zero().set_linear_form_3(
                s1 * s1,
                &self.d2ndu1,
                2.0 * s1 * s2,
                &self.d2nduv1,
                s2 * s2,
                &self.d2ndv1,
            );
            let tempbis = GpVec::zero().set_linear_form_3(
                2.0 * s1,
                &self.d2ndtu1,
                2.0 * s2,
                &self.d2ndtv1,
                1.0,
                &self.d2n1w,
            );
            let temp = temp.added(&tempbis);
            d2norm1w = GpVec::zero()
                .set_linear_form_2(secmember.value(1), &self.dndu1, secmember.value(2), &self.dndv1)
                .added(&temp);

            let temp = GpVec::zero().set_linear_form_3(
                s3 * s3,
                &self.d2ndu2,
                2.0 * s3 * s4,
                &self.d2nduv2,
                s4 * s4,
                &self.d2ndv2,
            );
            let tempbis = GpVec::zero().set_linear_form_3(
                2.0 * s3,
                &self.d2ndtu2,
                2.0 * s4,
                &self.d2ndtv2,
                1.0,
                &self.d2n2w,
            );
            let temp = temp.added(&tempbis);
            d2norm2w = GpVec::zero()
                .set_linear_form_2(secmember.value(3), &self.dndu2, secmember.value(4), &self.dndv2)
                .added(&temp);
        }

        // Tops 2d (`cxx:1770-1779`).
        poles2d[0] = GpPnt2d::new(x.value(1), x.value(2));
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(x.value(3), x.value(4));
        if !istgt {
            dpoles2d[0] = GpVec2d::new(sol.value(1), sol.value(2));
            dpoles2d[last2d] = GpVec2d::new(sol.value(3), sol.value(4));
            d2poles2d[0] = GpVec2d::new(secmember.value(1), secmember.value(2));
            d2poles2d[last2d] = GpVec2d::new(secmember.value(3), secmember.value(4));
        }

        // The linear case is processed (`cxx:1781-1800`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts1;
            poles[upp] = self.pts2;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            if !istgt {
                // OCCT writes `DPoles` twice here (`cxx:1790-1793`): the second
                // pair of assignments overwrites the tangents with `dtg1` /
                // `dtg2`. Kept verbatim.
                dpoles[low] = self.tg1;
                dpoles[upp] = self.tg2;
                dpoles[low] = dtg1;
                dpoles[upp] = dtg2;
                dweights[low] = 0.0;
                dweights[upp] = 0.0;
                d2weights[low] = 0.0;
                d2weights[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle (`cxx:1802-1905`).
        let mut norm1 = self.nplan.crossed(&ns1).magnitude();
        let mut norm2 = self.nplan.crossed(&ns2).magnitude();
        if norm1 < EPS {
            norm1 = 1.0;
        }
        if norm2 < EPS {
            norm2 = 1.0;
        }

        ns1 = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns1) / norm1,
            &self.nplan,
            -1.0 / norm1,
            &ns1,
        );
        ns2 = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns2) / norm2,
            &self.nplan,
            -1.0 / norm2,
            &ns2,
        );

        let center = GpPnt::from_xyz(&self.pts1.coord.added(&ns1.coord.multiply(self.ray1)));
        let (tgc, dtgc) = if !istgt {
            (
                GpVec::zero().set_linear_form_2(self.ray1, &dnorm1w, 1.0, &self.tg1),
                GpVec::zero().set_linear_form_2(self.ray1, &d2norm1w, 1.0, &dtg1),
            )
        } else {
            (GpVec::zero(), GpVec::zero())
        };

        // ns1 is oriented from the center to pts1 and ns2 from the center to
        // pts2; the trihedron ns1, ns2, nplan is made direct
        // (`cxx:1832-1855`).
        if self.ray1 > 0.0 {
            ns1.reverse();
            if !istgt {
                dnorm1w.reverse();
                d2norm1w.reverse();
            }
        }
        if self.ray2 > 0.0 {
            ns2.reverse();
            if !istgt {
                dnorm2w.reverse();
                d2norm2w.reverse();
            }
        }
        if self.choix % 2 != 0 {
            np.reverse();
            dnp.reverse();
            d2np.reverse();
        }

        if !istgt {
            section_convert::get_circle_d2(
                self.my_tconv,
                &ns1,
                &ns2,
                &dnorm1w,
                &dnorm2w,
                &d2norm1w,
                &d2norm2w,
                &np,
                &dnp,
                &d2np,
                &self.pts1,
                &self.pts2,
                &self.tg1,
                &self.tg2,
                &dtg1,
                &dtg2,
                self.ray1.abs(),
                0.0,
                0.0,
                &center,
                &tgc,
                &dtgc,
                poles,
                dpoles,
                d2poles,
                weights,
                dweights,
                d2weights,
            );
            true
        } else {
            section_convert::get_circle(
                self.my_tconv,
                &ns1,
                &ns2,
                &np,
                &self.pts1,
                &self.pts2,
                self.ray1.abs(),
                &center,
                poles,
                weights,
            );
            false
        }
    }

    /// `secmember = -(D2EDT2 + (2 * D2EDXDT + D2DXdSdt) * sol)`
    /// (`cxx:1714`, `:1728`).
    fn second_member(&self, d2dx_dsdt: &MathMatrix, sol: &MathVector) -> MathVector {
        let mut two_d2edxdt = self.d2edxdt.multiplied_scalar(2.0);
        two_d2edxdt.add(d2dx_dsdt);
        let mv = mat_vec(&two_d2edxdt, sol);
        self.d2edt2.added(&mv).multiplied_scalar(-1.0)
    }
}
