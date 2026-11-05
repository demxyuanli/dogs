//! Section generation of `BlendFunc_ConstRad` — the analytic circle, the
//! section without derivatives, the first-order section, and the
//! shape/tolerance/interval queries the approximation driver asks for.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRad.cxx:1142-1199`
//! (`Section(Param, ...)`), `:1204-1276` (`IsRational` ... `Mults`),
//! `:1284-1352` (`Section(P, Poles, ...)`), `:1356-1524`
//! (`Section(P, ..., DPoles, ...)`) (OCCT 8.0.0). The second-order section is
//! in [`super::section_order2`].

use std::f64::consts::PI;

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;
use occt_geom::geom_fill::section_convert;
use occt_math::{svd_padded, svd_solve_occt, Gauss, MathVector};

use crate::blend::{
    get_minimal_weights, get_shape as blend_get_shape, next_shape, BlendFuncSectionShape,
    BlendPoint,
};

use super::constants::{EPS, GAUSS_SOLUTION_EPS, POINTS_EQUAL_TOL, SVD_EPS};
use super::ConstRad;

/// `gp_Ax2(Center, np, ns)` from two vectors: OCCT raises when they are
/// parallel, and the port panics at the same place.
fn ax2_from_vecs(center: &GpPnt, z: &GpVec, x: &GpVec) -> GpAx2 {
    let z_dir = GpDir::from_vec(z).expect("BlendFunc_ConstRad::Section: degenerate axis");
    let x_dir = GpDir::from_vec(x).expect("BlendFunc_ConstRad::Section: degenerate axis");
    GpAx2::new(*center, z_dir, x_dir)
        .expect("BlendFunc_ConstRad::Section: gp_Ax2 construction error")
}

impl ConstRad {
    /// `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)` (`cxx:1142-1199`): the
    /// cross-section circle at the guide parameter `param`, for the solution
    /// `(U1, V1, U2, V2)`. Used for a quick approximate visualization.
    ///
    /// OCCT ignores the `ComputeValues` result here (`cxx:1160`), so a failed
    /// evaluation continues with whatever the fields hold.
    #[allow(clippy::too_many_arguments)]
    pub fn section_circle(
        &mut self,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> (f64, f64, GpCirc) {
        let mut x = MathVector::new(1, 4);
        x.set_value(1, u1);
        x.set_value(2, v1);
        x.set_value(3, u2);
        x.set_value(4, v2);
        let prm = param;
        self.compute_values(&x, 0, true, prm);

        let mut ns1 = self.nsurf1;
        let mut np = self.nplan;

        let mut norm1 = self.nplan.crossed(&ns1).magnitude();
        if norm1 < EPS {
            // Unsatisfactory, but it is not necessary to stop (`cxx:1169`).
            norm1 = 1.0;
        }
        ns1 = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns1) / norm1,
            &self.nplan,
            -1.0 / norm1,
            &ns1,
        );
        let center = GpPnt::from_xyz(&self.pts1.coord.added(&ns1.coord.multiply(self.ray1)));

        // ns1 is oriented from the center to pts1 (`cxx:1176`).
        if self.ray1 > 0.0 {
            ns1.reverse();
        }
        if self.choix % 2 != 0 {
            np.reverse();
        }
        let mut c = GpCirc::new(GpAx2::standard(), self.ray1.abs());
        c.set_radius(self.ray1.abs())
            .expect("BlendFunc_ConstRad::Section: invalid radius");
        c.set_position(&ax2_from_vecs(&center, &np, &ns1));
        let pdeb = 0.0;
        let mut pfin = clib::circle_parameter(&c.position(), &self.pts2);

        // Test negative and almost null angles : Singular case
        // (`cxx:1188-1197`).
        if pfin > 1.5 * PI {
            np.reverse();
            c.set_position(&ax2_from_vecs(&center, &np, &ns1));
            pfin = clib::circle_parameter(&c.position(), &self.pts2);
        }
        if pfin < PCONFUSION {
            pfin += PCONFUSION;
        }
        (pdeb, pfin, c)
    }

    /// `IsRational()` (`cxx:1204-1207`).
    pub(super) fn shape_is_rational(&self) -> bool {
        self.my_sshape == BlendFuncSectionShape::Rational
            || self.my_sshape == BlendFuncSectionShape::QuasiAngular
    }

    /// `GetSectionSize()` (`cxx:1211-1214`): the length of the maximum section.
    pub(super) fn section_size(&self) -> f64 {
        self.maxang * self.ray1.abs()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:1218-1221`). OCCT raises out of
    /// `GeomConvert` for an angle out of range (i.e. before any `IsSolution`);
    /// the port panics at the same place.
    pub(super) fn minimal_weights(&self, weights: &mut [f64]) {
        get_minimal_weights(
            self.my_sshape,
            self.my_tconv,
            self.minang,
            self.maxang,
            weights,
        )
        .expect("BlendFunc_ConstRad::GetMinimalWeight: section angle out of range");
    }

    /// `NbIntervals(S)` (`cxx:1225-1228`): `curv->NbIntervals(NextShape(S))`.
    pub(super) fn guide_nb_intervals(&self, s: Shape) -> i32 {
        self.curv.nb_intervals(next_shape(s) as u8)
    }

    /// `Intervals(T, S)` (`cxx:1232-1235`): `curv->Intervals(T, NextShape(S))`.
    pub(super) fn guide_intervals(&self, t: &mut [f64], s: Shape) {
        let src = self.curv.parameter_intervals(next_shape(s) as u8);
        let count = t.len().min(src.len());
        t[..count].copy_from_slice(&src[..count]);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:1239-1243`).
    pub(super) fn section_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        *nb_poles2d = 2;
        blend_get_shape(
            self.my_sshape,
            self.maxang,
            nb_poles,
            nb_knots,
            degree,
            &mut self.my_tconv,
        );
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)`
    /// (`cxx:1249-1262`). `GeomFill::GetTolerance` raises for an angle out of
    /// range in OCCT; the port panics at the same place.
    pub(super) fn fill_section_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut MathVector,
        tol1d: &mut MathVector,
    ) {
        let tol = section_convert::get_tolerance(
            self.my_tconv,
            self.minang,
            self.ray1.abs(),
            angle_tol,
            surf_tol,
        )
        .expect("BlendFunc_ConstRad::GetTolerance: section angle out of range");
        tol1d.init(surf_tol);
        tol3d.init(surf_tol);
        let low = tol3d.lower();
        let up = tol3d.upper();
        let inner = tol.min(surf_tol);
        tol3d.set_value(low + 1, inner);
        tol3d.set_value(up - 1, inner);
        let outer = tol.min(bound_tol);
        tol3d.set_value(low, outer);
        tol3d.set_value(up, outer);
    }

    /// `Knots(TKnots)` (`cxx:1266-1269`).
    pub(super) fn section_knots(&self, tknots: &mut [f64]) {
        section_convert::knots(self.my_tconv, tknots);
    }

    /// `Mults(TMults)` (`cxx:1273-1276`).
    pub(super) fn section_mults(&self, tmults: &mut [i32]) {
        section_convert::mults(self.my_tconv, tmults);
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:1284-1352`): the section
    /// without derivatives.
    pub(super) fn section_poles_impl(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        let prm = p.parameter();
        let low = 0usize;
        let upp = poles.len() - 1;

        let mut x = MathVector::new(1, 4);
        let (u1, v1) = p
            .parameters_on_s1()
            .expect("BlendFunc_ConstRad::Section: point is not on the first surface");
        let (u2, v2) = p
            .parameters_on_s2()
            .expect("BlendFunc_ConstRad::Section: point is not on the second surface");
        x.set_value(1, u1);
        x.set_value(2, v1);
        x.set_value(3, u2);
        x.set_value(4, v2);

        self.compute_values(&x, 0, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(&self.pts2));

        // ns1, ns2, np are copied locally to avoid crushing the fields
        // (`cxx:1301-1303`).
        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let mut np = self.nplan;

        poles2d[0] = GpPnt2d::new(x.value(1), x.value(2));
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(x.value(3), x.value(4));

        // Linear case (`cxx:1308-1316`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts1;
            poles[upp] = self.pts2;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            return;
        }

        let mut norm1 = self.nplan.crossed(&ns1).magnitude();
        let mut norm2 = self.nplan.crossed(&ns2).magnitude();
        if norm1 < EPS {
            // Unsatisfactory, but it is not necessary to stop (`cxx:1321`).
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

        // ns1 (resp. ns2) is oriented from the center to pts1 (resp. pts2),
        // and the trihedron ns1, ns2, nplan is made direct (`cxx:1335-1349`).
        if self.ray1 > 0.0 {
            ns1.reverse();
        }
        if self.ray2 > 0.0 {
            ns2.reverse();
        }
        if self.choix % 2 != 0 {
            np.reverse();
        }

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
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:1356-1524`): the first/last section with first derivatives.
    /// Returns whether the derivatives were computed.
    ///
    /// `nplan` / `dnplan` are copied locally: the `choix % 2 != 0` arm reverses
    /// the copies, not the fields (`cxx:1491-1496`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn section_with_derivatives(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        dpoles: &mut [GpVec],
        poles2d: &mut [GpPnt2d],
        dpoles2d: &mut [GpVec2d],
        weights: &mut [f64],
        dweights: &mut [f64],
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
        let mut sol = MathVector::new(1, 4);
        sol.set_value(1, u1);
        sol.set_value(2, v1);
        sol.set_value(3, u2);
        sol.set_value(4, v2);

        // Calculation of equations (`cxx:1377-1378`).
        self.compute_values(&sol, 1, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(&self.pts2));

        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let mut np = self.nplan;
        let mut dnp = self.dnplan;
        let mut dnorm1w = GpVec::zero();
        let mut dnorm2w = GpVec::zero();

        let mut secmember = MathVector::new(1, 4);
        let mut istgt = true;

        // Calculation of derivatives. Processing Normal (`cxx:1386-1408`).
        if self.pts1.distance(&self.pts2) > POINTS_EQUAL_TOL {
            let resol = Gauss::new(&self.dedx, GAUSS_SOLUTION_EPS);
            if resol.is_done() {
                let minus_dedt = self.dedt.multiplied_scalar(-1.0);
                if let Some(x) = resol.solve(&minus_dedt) {
                    secmember = x;
                }
                istgt = false;
            }
        }

        if istgt {
            let sing_rs = svd_padded(&self.dedx);
            if sing_rs.is_done {
                let minus_dedt = self.dedt.multiplied_scalar(-1.0);
                svd_solve_occt(&sing_rs, &minus_dedt, &mut secmember, SVD_EPS);
                istgt = false;
            }
        }

        if !istgt {
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

            dnorm1w = GpVec::zero().set_linear_form_3(
                secmember.value(1),
                &self.dndu1,
                secmember.value(2),
                &self.dndv1,
                1.0,
                &self.dn1w,
            );
            dnorm2w = GpVec::zero().set_linear_form_3(
                secmember.value(3),
                &self.dndu2,
                secmember.value(4),
                &self.dndv2,
                1.0,
                &self.dn2w,
            );
        }

        // Tops 2d (`cxx:1420-1428`).
        poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(sol.value(3), sol.value(4));
        if !istgt {
            dpoles2d[0] = GpVec2d::new(secmember.value(1), secmember.value(2));
            dpoles2d[last2d] = GpVec2d::new(secmember.value(3), secmember.value(4));
        }

        // The linear case is processed (`cxx:1429-1445`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts1;
            poles[upp] = self.pts2;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            if !istgt {
                dpoles[low] = self.tg1;
                dpoles[upp] = self.tg2;
                dweights[low] = 0.0;
                dweights[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle (`cxx:1447-1524`).
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
        let tgc = if !istgt {
            // `= tg1.Added(ray1*dn1w)` (`cxx:1469`).
            dnorm1w
                .multiplied_scalar(self.ray1)
                .added(&self.tg1)
        } else {
            GpVec::zero()
        };

        // ns1 is oriented from the center to pts1, ns2 from the center to
        // pts2, and the trihedron ns1, ns2, nplan is made direct
        // (`cxx:1472-1496`).
        if self.ray1 > 0.0 {
            ns1.reverse();
            if !istgt {
                dnorm1w.reverse();
            }
        }
        if self.ray2 > 0.0 {
            ns2.reverse();
            if !istgt {
                dnorm2w.reverse();
            }
        }
        if self.choix % 2 != 0 {
            np.reverse();
            dnp.reverse();
        }

        if !istgt {
            section_convert::get_circle_d1(
                self.my_tconv,
                &ns1,
                &ns2,
                &dnorm1w,
                &dnorm2w,
                &np,
                &dnp,
                &self.pts1,
                &self.pts2,
                &self.tg1,
                &self.tg2,
                self.ray1.abs(),
                0.0,
                &center,
                &tgc,
                poles,
                dpoles,
                weights,
                dweights,
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
}
