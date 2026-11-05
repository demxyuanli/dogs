//! Section generation of `BlendFunc_Ruled`: the section is the straight
//! segment between the two extremities, so the "poles" are just the two ends
//! and their tangents are the segment direction.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_Ruled.cxx:469-568`
//! (`GetSection`), `:572-637` (`IsRational` ... `Mults`), `:638-711` (the three
//! `Section` overloads), `:713-730` (`AxeRot`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `GetSection` (`cxx:469-568`) and `AxeRot` (`cxx:713-730`) have no caller
//!   in the OCCT tree; they are ported as inherent methods. `GetSection`
//!   raises `Standard_RangeError` for a too short or mismatched array pair
//!   (`cxx:479-482`) — a panic here.
//! - `P.ParametersOnS1(u, v)` returns `false` when the point carries no 2d
//!   parameters and leaves `u` / `v` untouched in OCCT, which reads them
//!   anyway (`cxx:666-670`, `:704-708`); the port panics at that point.
//! - `Knots` writes `TKnots(Lower())` and `TKnots(Upper())` only (`cxx:628-632`),
//!   `Mults` writes the first and last entry (`cxx:634-637`).

use occt_core::gp::{GpAx1, GpDir, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_math::{Gauss, MathMatrix, MathVector};

use crate::blend::{next_shape, BlendPoint};

use super::{Ruled, GAUSS_MIN_PIVOT};

impl Ruled {
    /// `NbIntervals(S)` (`cxx:593-596`): `curv->NbIntervals(NextShape(S))`.
    pub(super) fn guide_nb_intervals(&self, s: Shape) -> i32 {
        self.curv.nb_intervals(next_shape(s) as u8)
    }

    /// `Intervals(T, S)` (`cxx:600-603`): `curv->Intervals(T, NextShape(S))`.
    pub(super) fn guide_intervals(&self, t: &mut [f64], s: Shape) {
        let src = self.curv.parameter_intervals(next_shape(s) as u8);
        let count = t.len().min(src.len());
        t[..count].copy_from_slice(&src[..count]);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:606-612`): a
    /// degree 1 segment, `2` poles and `2` poles 2d.
    pub(super) fn section_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        *nb_poles = 2;
        *nb_knots = 2;
        *degree = 1;
        *nb_poles2d = 2;
    }

    /// `Knots(TKnots)` (`cxx:628-632`): `0.` and `1.`.
    pub(super) fn section_knots(&mut self, tknots: &mut [f64]) {
        tknots[0] = 0.0;
        let last = tknots.len() - 1;
        tknots[last] = 1.0;
    }

    /// `Mults(TMults)` (`cxx:634-637`): `2` at both ends.
    pub(super) fn section_mults(&mut self, tmults: &mut [i32]) {
        tmults[0] = 2;
        let last = tmults.len() - 1;
        tmults[last] = 2;
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:652-690`): the two extremities and their 2d parameters; the
    /// derivatives only when `P` is not a tangency point.
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
        poles[0] = *p.point_on_s1();
        poles[1] = *p.point_on_s2();

        let (u, v) = p
            .parameters_on_s1()
            .expect("BlendFunc_Ruled::Section: point is not on the first surface");
        poles2d[0] = GpPnt2d::new(u, v);
        let (u, v) = p
            .parameters_on_s2()
            .expect("BlendFunc_Ruled::Section: point is not on the second surface");
        poles2d[1] = GpPnt2d::new(u, v);

        weights[0] = 1.0;
        weights[1] = 1.0;

        if !p.is_tangency_point() {
            dpoles[0] = *p
                .tangent_on_s1()
                .expect("BlendFunc_Ruled::Section: no tangent on the first surface");
            dpoles[1] = *p
                .tangent_on_s2()
                .expect("BlendFunc_Ruled::Section: no tangent on the second surface");

            dpoles2d[0] = p
                .tangent2d_on_s1()
                .expect("BlendFunc_Ruled::Section: no 2d tangent on the first surface");
            dpoles2d[1] = p
                .tangent2d_on_s2()
                .expect("BlendFunc_Ruled::Section: no 2d tangent on the second surface");

            dweights[0] = 0.0;
            dweights[1] = 0.0;

            return true;
        }

        false
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:692-711`).
    pub(super) fn section_poles_impl(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        poles[0] = *p.point_on_s1();
        poles[1] = *p.point_on_s2();

        let (u, v) = p
            .parameters_on_s1()
            .expect("BlendFunc_Ruled::Section: point is not on the first surface");
        poles2d[0] = GpPnt2d::new(u, v);
        let (u, v) = p
            .parameters_on_s2()
            .expect("BlendFunc_Ruled::Section: point is not on the second surface");
        poles2d[1] = GpPnt2d::new(u, v);

        weights[0] = 1.0;
        weights[1] = 1.0;
    }

    /// `GetSection(Param, U1, V1, U2, V2, tabP, tabV)` (`cxx:469-568`): fills
    /// the section polyline for the solution `(U1, V1, U2, V2)` at the guide
    /// parameter `param`. `false` when the tangent system is singular.
    #[allow(clippy::too_many_arguments)]
    pub fn get_section(
        &mut self,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tab_p: &mut [GpPnt],
        tab_v: &mut [GpVec],
    ) -> bool {
        let nb_point = tab_p.len();
        if nb_point != tab_v.len() || nb_point < 2 {
            panic!("BlendFunc_Ruled::GetSection: Standard_RangeError");
        }

        let mut sol = MathVector::new(1, 4);
        let mut valsol = MathVector::new(1, 4);
        let mut gradsol = MathMatrix::new(1, 4, 1, 4);

        let (ptgui, d1gui, d2gui) = self.curv.d2(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        self.the_d = -self.nplan.dot(&GpVec::from_xyz(self.ptgui.xyz()));

        sol.set_value(1, u1);
        sol.set_value(2, v1);
        sol.set_value(3, u2);
        sol.set_value(4, v2);

        self.compute_values(&sol, Some(&mut valsol), &mut gradsol);

        let (_pts1, d1u1, d1v1) = self.surf1.d1(sol.value(1), sol.value(2));
        let (_pts2, d1u2, d1v2) = self.surf2.d1(sol.value(3), sol.value(4));
        // `pts1` / `pts2` come from `Values` above (`cxx:504-505`).

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
        // (`cxx:515-522`).
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
        // (`cxx:531-540`).
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

            tab_p[0] = self.pts1;
            tab_p[nb_point - 1] = self.pts2;

            tab_v[0] = self.tg1;
            tab_v[nb_point - 1] = self.tg2;

            for i in 2..=nb_point - 1 {
                let lambda = (i - 1) as f64 / (nb_point - 1) as f64;
                let xyz = self
                    .pts1
                    .xyz()
                    .multiplied(1.0 - lambda)
                    .added(&self.pts2.xyz().multiplied(lambda));
                tab_p[i - 2] = GpPnt::from_xyz(&xyz);
                tab_v[i - 2] = GpVec::zero().set_linear_form_2(
                    1.0 - lambda,
                    &self.tg1,
                    lambda,
                    &self.tg2,
                );
            }
            return true;
        }
        false
    }

    /// `AxeRot(Prm)` (`cxx:713-730`): the rotation axis of the section frame
    /// at the guide parameter `prm`.
    pub fn axe_rot(&mut self, prm: f64) -> GpAx1 {
        let (ptgui, d1gui, d2gui) = self.curv.d2(prm);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;

        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let dirax = self.nplan.crossed(&dnplan);
        let mut axrot = GpAx1::default();
        axrot.set_direction(
            GpDir::from_vec(&dirax)
                .expect("BlendFunc_Ruled::AxeRot: gp_Dir from a zero direction"),
        );
        let oriax = self
            .ptgui
            .xyz()
            .added(&dnplan.normalized().xyz().multiplied(self.normtg / dnplan.magnitude()));
        axrot.set_location(GpPnt::from_xyz(&oriax));
        axrot
    }
}
