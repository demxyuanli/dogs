//! Section generation of `BlendFunc_CSConstRad` — the analytic circle, the
//! section without derivatives, the first-order section, and the
//! shape/tolerance/interval queries the approximation driver asks for.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSConstRad.cxx:467-504`
//! (`Section(Param, ...)`), `:506-517` (the throwing second-order overload),
//! `:523-670` (`GetSection`), `:673-747` (`IsRational` ... `Mults`),
//! `:750-802` (`Section(P, Poles, ...)`), `:804-981`
//! (`Section(P, ..., DPoles, ...)`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Section(Param, ...)` sets `nplan` from `guide->D1`, while `GetSection` and
//!   both `Section(P, ...)` overloads use `guide->D2` because they need
//!   `dnplan`; each keeps its own `guide` call.
//! - `nplan` is reversed **in place** in `Section(P, ...)` (`cxx:797-800`,
//!   `:950-953`) but `nplan.Reversed()` — a copy — is used in
//!   `Section(Param, ...)` (`cxx:496-501`); the distinction is kept.
//! - The projected normal has no `Eps` guard in this class: `|nplan x ns|` is
//!   divided by as written.
//! - `GeomFill::GetCircle` answers `Result` for a pole array of the wrong
//!   length (OCCT's `Standard_OutOfRange`); the array always comes from
//!   `GetShape`, so the failure is discarded exactly like the other `blend`
//!   section routines.

use std::f64::consts::PI;

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::geom_fill::section_convert;
use occt_math::{Gauss, MathVector};

use crate::blend::{
    get_minimal_weights, get_shape as blend_get_shape, next_shape, BlendFuncSectionShape,
    BlendPoint,
};

use super::constants::GAUSS_MIN_PIVOT;
use super::CsConstRad;

/// `gp_Ax2(Center, np, ns)` from two vectors: OCCT raises when they are
/// parallel, and the port panics at the same place.
fn ax2_from_vecs(center: &GpPnt, z: &GpVec, x: &GpVec) -> GpAx2 {
    let z_dir = GpDir::from_vec(z).expect("BlendFunc_CSConstRad::Section: degenerate axis");
    let x_dir = GpDir::from_vec(x).expect("BlendFunc_CSConstRad::Section: degenerate axis");
    GpAx2::new(*center, z_dir, x_dir)
        .expect("BlendFunc_CSConstRad::Section: gp_Ax2 construction error")
}

impl CsConstRad {
    /// `Section(Param, U, V, W, Pdeb, Pfin, C)` (`cxx:467-504`): the
    /// cross-section circle at the guide parameter `param`, for the solution
    /// `(U, V, W)`. Used for a quick approximate visualization.
    #[allow(clippy::too_many_arguments)]
    pub fn section_circle(&mut self, param: f64, u: f64, v: f64, w: f64) -> (f64, f64, GpCirc) {
        let (ptgui, d1gui) = self.guide.d1(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.nplan = self.d1gui.normalized();

        let (pts, d1u1, d1v1) = self.surf.d1(u, v);
        self.pts = pts;
        self.ptc = self.curv.value(w);

        let ns = d1u1.crossed(&d1v1);
        let norm = self.nplan.crossed(&ns).magnitude();
        let ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );

        let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));
        let mut c = GpCirc::new(GpAx2::standard(), self.ray.abs());
        c.set_radius(self.ray.abs())
            .expect("BlendFunc_CSConstRad::Section: invalid radius");

        let np = if self.choix % 2 == 0 {
            self.nplan
        } else {
            self.nplan.reversed()
        };
        c.set_position(&ax2_from_vecs(&center, &np, &ns));

        let pdeb = 0.0;
        let pfin = clib::circle_parameter(&c.position(), &self.ptc);
        (pdeb, pfin, c)
    }

    /// `GetSection(Param, U, V, W, tabP, tabV)` (`cxx:523-670`): the poles and
    /// the tangents of the section, used by the approximation when the section
    /// is evaluated from a stored solution instead of a `Blend_Point`.
    pub fn get_section(
        &mut self,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
        tab_p: &mut [GpPnt],
        tab_v: &mut [GpVec],
    ) -> bool {
        let nb_point = tab_p.len();
        if nb_point != tab_v.len() || nb_point < 2 {
            panic!("Standard_RangeError");
        }

        let (ptgui, d1gui, d2gui) = self.guide.d2(param);
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

        let mut sol = MathVector::new(1, 3);
        sol.set_value(1, u);
        sol.set_value(2, v);
        sol.set_value(3, w);
        self.evaluate_values(&sol);

        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(u, v);
        self.pts = pts;
        let (ptc, d1) = self.curv.d1(w);
        self.ptc = ptc;

        let mut secmember = MathVector::new(1, 3);
        let temp = GpVec::from_pnts(&self.ptgui, &self.pts);
        secmember.set_value(1, self.normtg - dnplan.dot(&temp));

        let temp = GpVec::from_pnts(&self.ptgui, &self.ptc);
        secmember.set_value(2, self.normtg - dnplan.dot(&temp));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        // Derivative of n1 in relation to w (`cxx:574-582`).
        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let mut dnw = GpVec::zero().set_linear_form_3(
            (dnplan.dot(&ns) - grosterme * ndotns) / norm,
            &self.nplan,
            ndotns / norm,
            &dnplan,
            grosterme / norm,
            &ns,
        );

        let temp = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);
        let resul = temp
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));
        // `jag 950105 il manquait ray` (`cxx:585`).
        secmember.set_value(3, -2.0 * self.ray * dnw.dot(&resul));

        let resol = Gauss::new(&self.dedx, GAUSS_MIN_PIVOT);
        if !resol.is_done() {
            return false;
        }

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

        // Derivative of n1 in relation to u (`cxx:595-603`).
        let mut temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
        let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
        let resulu = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns - self.nplan.dot(&temp)) / norm,
            &self.nplan,
            grosterme / norm,
            &ns,
            -1.0 / norm,
            &temp,
        );

        // Derivative of n1 in relation to v (`cxx:605-614`).
        temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
        let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
        let resulv = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns - self.nplan.dot(&temp)) / norm,
            &self.nplan,
            grosterme / norm,
            &ns,
            -1.0 / norm,
            &temp,
        );

        dnw = GpVec::zero().set_linear_form_3(
            secmember.value(1),
            &resulu,
            secmember.value(2),
            &resulv,
            1.0,
            &dnw,
        );
        let ns = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);

        let dn2w = GpVec::zero().set_linear_form_2plus(self.ray, &dnw, -1.0, &self.tgc, &self.tgs);
        let dn2w = dn2w.divide(resul.magnitude());
        let ns2 = resul.normalized().reversed();
        let dn2w = GpVec::zero().set_linear_form_2(ns2.dot(&dn2w), &ns2, -1.0, &dn2w);

        let mut np = self.nplan;
        let mut dnp = dnplan;
        if self.choix % 2 != 0 {
            np.reverse();
            dnp.reverse();
        }

        tab_p[0] = self.pts;
        tab_p[nb_point - 1] = self.ptc;

        tab_v[0] = self.tgs;
        tab_v[nb_point - 1] = self.tgc;

        if nb_point > 2 {
            let cosa = ns.dot(&ns2);
            let sina = np.dot(&ns.crossed(&ns2));
            let mut angle = cosa.acos();
            if sina < 0.0 {
                angle = 2.0 * PI - angle;
            }
            let dangle = -(dnw.dot(&ns2) + ns.dot(&dn2w)) / sina;
            let ncrn = np.crossed(&ns);
            let dncrn = dnp.crossed(&ns).added(&np.crossed(&dnw));

            for i in 2..=nb_point - 1 {
                let lambda = (i - 1) as f64 / (nb_point - 1) as f64;
                let cosa = (lambda * angle).cos();
                let sina = (lambda * angle).sin();
                let offset = GpVec::zero().set_linear_form_2(cosa - 1.0, &ns, sina, &ncrn);
                tab_p[i - 1] = self
                    .pts
                    .translated_vec(&offset.multiplied_scalar(self.ray.abs()));

                let mut temp = GpVec::zero().set_linear_form_2(-sina, &ns, cosa, &ncrn);
                temp = temp.multiplied_scalar(lambda * dangle);
                temp = temp.added(&GpVec::zero().set_linear_form_2(
                    cosa - 1.0,
                    &dnw,
                    sina,
                    &dncrn,
                ));
                temp = temp.multiplied_scalar(self.ray.abs());
                tab_v[i - 1] = temp.added(&self.tgs);
            }
        }
        true
    }

    /// `IsRational()` (`cxx:673-676`).
    pub(super) fn shape_is_rational(&self) -> bool {
        self.my_sshape == BlendFuncSectionShape::Rational
            || self.my_sshape == BlendFuncSectionShape::QuasiAngular
    }

    /// `GetSectionSize()` (`cxx:680-683`): the length of the maximum section.
    pub(super) fn section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:687-691`). OCCT raises out of
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
        .expect("BlendFunc_CSConstRad::GetMinimalWeight: section angle out of range");
    }

    /// `NbIntervals(S)` (`cxx:695-698`): `curv->NbIntervals(NextShape(S))`.
    pub(super) fn curve_nb_intervals(&self, s: Shape) -> i32 {
        self.curv.nb_intervals(next_shape(s) as u8)
    }

    /// `Intervals(T, S)` (`cxx:702-705`): `curv->Intervals(T, NextShape(S))`.
    pub(super) fn curve_intervals(&self, t: &mut [f64], s: Shape) {
        let src = self.curv.parameter_intervals(next_shape(s) as u8);
        let count = t.len().min(src.len());
        t[..count].copy_from_slice(&src[..count]);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:710-714`).
    pub(super) fn section_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        *nb_poles2d = 1;
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
    /// (`cxx:718-733`). `GeomFill::GetTolerance` raises for an angle out of
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
            self.ray.abs(),
            angle_tol,
            surf_tol,
        )
        .expect("BlendFunc_CSConstRad::GetTolerance: section angle out of range");
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

    /// `Knots(TKnots)` (`cxx:736-739`).
    pub(super) fn section_knots(&self, tknots: &mut [f64]) {
        section_convert::knots(self.my_tconv, tknots);
    }

    /// `Mults(TMults)` (`cxx:742-745`).
    pub(super) fn section_mults(&self, tmults: &mut [i32]) {
        section_convert::mults(self.my_tconv, tmults);
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:750-802`): the section
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

        let (ptgui, d1gui) = self.guide.d1(prm);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.nplan = self.d1gui.normalized();

        let (u1, v1) = p
            .parameters_on_s()
            .expect("BlendFunc_CSConstRad::Section: point is not on the surface");
        let w = p
            .parameter_on_c()
            .expect("BlendFunc_CSConstRad::Section: point is not on the curve");

        let (pts, d1u1, d1v1) = self.surf.d1(u1, v1);
        self.pts = pts;
        self.ptc = self.curv.value(w);

        poles2d[0] = GpPnt2d::new(u1, v1);

        // Linear case (`cxx:777-785`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptc;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            return;
        }

        let ns = d1u1.crossed(&d1v1);
        let norm = self.nplan.crossed(&ns).magnitude();
        let ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );

        let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));
        let ns2 = GpVec::from_pnts(&center, &self.ptc).normalized();

        if self.choix % 2 != 0 {
            self.nplan.reverse();
        }

        section_convert::get_circle(
            self.my_tconv,
            &ns,
            &ns2,
            &self.nplan,
            &self.pts,
            &self.ptc,
            self.ray.abs(),
            &center,
            poles,
            weights,
        );
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:804-981`): the first/last section with first derivatives.
    /// Returns whether the derivatives were computed.
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

        let (ptgui, d1gui, d2gui) = self.guide.d2(prm);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        let mut dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let (u, v) = p
            .parameters_on_s()
            .expect("BlendFunc_CSConstRad::Section: point is not on the surface");
        let mut sol = MathVector::new(1, 3);
        sol.set_value(1, u);
        sol.set_value(2, v);
        sol.set_value(
            3,
            p.parameter_on_c()
                .expect("BlendFunc_CSConstRad::Section: point is not on the curve"),
        );

        // Calculation of equations (`cxx:840`).
        self.evaluate_values(&sol);

        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(sol.value(1), sol.value(2));
        self.pts = pts;
        let (ptc, d1) = self.curv.d1(sol.value(3));
        self.ptc = ptc;

        let mut secmember = MathVector::new(1, 3);
        let temp = GpVec::from_pnts(&self.ptgui, &self.pts);
        secmember.set_value(1, self.normtg - dnplan.dot(&temp));

        let temp = GpVec::from_pnts(&self.ptgui, &self.ptc);
        secmember.set_value(2, self.normtg - dnplan.dot(&temp));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        // Derivative of n1 in relation to w (`cxx:856-866`).
        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let mut dnw = GpVec::zero().set_linear_form_3(
            (dnplan.dot(&ns) - grosterme * ndotns) / norm,
            &self.nplan,
            ndotns / norm,
            &dnplan,
            grosterme / norm,
            &ns,
        );

        let temp = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);
        let resul = temp
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));
        // `jag 950105 il manquait ray` (`cxx:868`).
        secmember.set_value(3, -2.0 * self.ray * dnw.dot(&resul));

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

            // Derivative of n1 in relation to u (`cxx:880-889`).
            let mut temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
            let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
            let resulu = GpVec::zero().set_linear_form_3(
                -(grosterme * ndotns - self.nplan.dot(&temp)) / norm,
                &self.nplan,
                grosterme / norm,
                &ns,
                -1.0 / norm,
                &temp,
            );

            // Derivative of n1 in relation to v (`cxx:891-898`).
            temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
            let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
            let resulv = GpVec::zero().set_linear_form_3(
                -(grosterme * ndotns - self.nplan.dot(&temp)) / norm,
                &self.nplan,
                grosterme / norm,
                &ns,
                -1.0 / norm,
                &temp,
            );

            dnw = GpVec::zero().set_linear_form_3(
                secmember.value(1),
                &resulu,
                secmember.value(2),
                &resulv,
                1.0,
                &dnw,
            );
            let ns = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);

            let dn2w = GpVec::zero().set_linear_form_2plus(self.ray, &dnw, -1.0, &self.tgc, &self.tgs);
            let norm = resul.magnitude();
            let dn2w = dn2w.divide(norm);
            let ns2 = resul.normalized().reversed();
            let dn2w = GpVec::zero().set_linear_form_2(ns2.dot(&dn2w), &ns2, -1.0, &dn2w);

            // Poles 2d (`cxx:920-923`).
            poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));
            dpoles2d[0] = GpVec2d::new(secmember.value(1), secmember.value(2));

            // Linear case (`cxx:926-940`).
            if self.my_sshape == BlendFuncSectionShape::Linear {
                poles[low] = self.pts;
                poles[upp] = self.ptc;
                weights[low] = 1.0;
                weights[upp] = 1.0;
                dpoles[low] = self.tgs;
                dpoles[upp] = self.tgc;
                dweights[low] = 0.0;
                dweights[upp] = 0.0;
                return true;
            }

            // Case of the circle (`cxx:943-980`).
            let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));
            let tgct = self.tgs.added(&dnw.multiplied_scalar(self.ray));

            if self.choix % 2 != 0 {
                self.nplan.reverse();
                dnplan.reverse();
            }

            section_convert::get_circle_d1(
                self.my_tconv,
                &ns,
                &ns2,
                &dnw,
                &dn2w,
                &self.nplan,
                &dnplan,
                &self.pts,
                &self.ptc,
                &self.tgs,
                &self.tgc,
                self.ray.abs(),
                0.0,
                &center,
                &tgct,
                poles,
                dpoles,
                weights,
                dweights,
            );
            return true;
        }

        // The tangent system failed: the section is still built without
        // derivatives (`cxx:911-917`).
        let ns = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);
        let ns2 = resul.normalized().reversed();

        // Poles 2d (`cxx:920-923`).
        poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));

        // Linear case (`cxx:926-940`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptc;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            return false;
        }

        let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));

        if self.choix % 2 != 0 {
            self.nplan.reverse();
        }

        section_convert::get_circle(
            self.my_tconv,
            &ns,
            &ns2,
            &self.nplan,
            &self.pts,
            &self.ptc,
            self.ray.abs(),
            &center,
            poles,
            weights,
        );
        false
    }
}
