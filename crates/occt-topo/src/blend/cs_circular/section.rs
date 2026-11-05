//! Section generation of `BlendFunc_CSCircular` — the analytic circle, the
//! section without derivatives, the first-order section, and the
//! shape/tolerance/interval queries the approximation driver asks for.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSCircular.cxx:490-531`
//! (`Section(Param, ...)`), `:557-707` (`GetSection`), `:711-786`
//! (`IsRational` ... `Mults`), `:788-842` (`Section(P, Poles, ...)`),
//! `:845-1026` (`Section(P, ..., DPoles, ...)`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - Unlike `BlendFunc_CSConstRad`, the projected normal of `GetSection`,
//!   `Section(Param, ...)` and `Section(P, Poles, ...)` has **no** `Eps` guard
//!   (`cxx:511-513`, `:593-595`, `:822-824`); only the equation bodies and
//!   `IsSolution` clamp `norm` to `1.` below `Eps`.
//! - `Section(Param, ...)` calls `guide->D1` and does not refresh `normtg`;
//!   `GetSection` and both `Section(P, ...)` overloads call `Set(Param)`.
//! - The reversal of `ns` / `dnw` for a positive radius (`cxx:657-661`,
//!   `:981-987`, `:830-833`) is applied **after** the centre and the tangents
//!   are built, so the stored values keep the un-reversed normal.
//! - `GeomFill::GetCircle` answers `()` here (OCCT answers `true`); the section
//!   routines do not forward it.

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
use super::CsCircular;

/// `gp_Ax2(Center, np, ns)` from two vectors: OCCT raises when they are
/// parallel, and the port panics at the same place.
fn ax2_from_vecs(center: &GpPnt, z: &GpVec, x: &GpVec) -> GpAx2 {
    let z_dir = GpDir::from_vec(z).expect("BlendFunc_CSCircular::Section: degenerate axis");
    let x_dir = GpDir::from_vec(x).expect("BlendFunc_CSCircular::Section: degenerate axis");
    GpAx2::new(*center, z_dir, x_dir).expect("BlendFunc_CSCircular::Section: gp_Ax2 construction error")
}

impl CsCircular {
    /// `Section(Param, U, V, W, Pdeb, Pfin, C)` (`cxx:490-531`): the
    /// cross-section circle at the guide parameter `param`, for the solution
    /// `(U, V, W)`. Used for a quick approximate visualization.
    #[allow(clippy::too_many_arguments)]
    pub fn section_circle(&mut self, param: f64, u: f64, v: f64, w: f64) -> (f64, f64, GpCirc) {
        let (_ptgui, d1gui) = self.guide.d1(param);
        self.d1gui = d1gui;
        self.nplan = self.d1gui.normalized();

        let (pts, d1u1, d1v1) = self.surf.d1(u, v);
        self.pts = pts;
        self.ptc = self.curv.value(w);

        let ns = d1u1.crossed(&d1v1);
        let norm = self.nplan.crossed(&ns).magnitude();
        let mut ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );

        let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));
        let mut c = GpCirc::new(GpAx2::standard(), self.ray.abs());
        c.set_radius(self.ray.abs())
            .expect("BlendFunc_CSCircular::Section: invalid radius");

        if self.ray > 0.0 {
            ns.reverse();
        }

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

    /// `GetSection(Param, U, V, W, tabP, tabV)` (`cxx:557-707`): the poles and
    /// the tangents of the section, used by the approximation when the section
    /// is evaluated from a stored solution instead of a `Blend_Point`.
    pub fn get_section(
        &mut self,
        param: f64,
        u: f64,
        v: f64,
        _w: f64,
        tab_p: &mut [GpPnt],
        tab_v: &mut [GpVec],
    ) -> bool {
        let nb_point = tab_p.len();
        if nb_point != tab_v.len() || nb_point < 2 {
            panic!("Standard_RangeError");
        }

        self.set_param(param);
        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let (ptc, d1c, _d2c) = self.curv.d2(self.prmc);
        self.ptc = ptc;
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(u, v);
        self.pts = pts;

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        let mut secmember = MathVector::new(1, 2);
        let temp = GpVec::from_pnts(&self.ptc, &self.pts);
        secmember.set_value(
            1,
            self.dprmc * self.nplan.dot(&d1c) - dnplan.dot(&temp),
        );

        let ns2 = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);

        // Derivative of n1 in relation to w (`cxx:601-609`).
        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let mut dnw = GpVec::zero().set_linear_form_3(
            (dnplan.dot(&ns) - grosterme * ndotns) / norm,
            &self.nplan,
            ndotns / norm,
            &dnplan,
            grosterme / norm,
            &ns,
        );

        let temp = dnw
            .multiplied_scalar(self.ray)
            .subtracted(&d1c.multiplied_scalar(self.dprmc));
        let resul = ns2
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));
        secmember.set_value(2, -2.0 * resul.dot(&temp));

        let mut sol = MathVector::new(1, 2);
        sol.set_value(1, u);
        sol.set_value(2, v);
        self.evaluate_values(&sol);

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
        self.tgc = d1c.multiplied_scalar(self.dprmc);

        // Derivative of n1 in relation to u (`cxx:629-638`).
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

        // Derivative of n1 in relation to v (`cxx:639-648`).
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
        // `ns = ns2` (`cxx:650`) before the circle is built.
        let mut ns = ns2;

        let dn2w = GpVec::zero().set_linear_form_2plus(self.ray, &dnw, -1.0, &self.tgc, &self.tgs);
        let norm = resul.magnitude();
        let dn2w = dn2w.divide(norm);
        let ns2 = resul.normalized().reversed();
        let dn2w = GpVec::zero().set_linear_form_2(ns2.dot(&dn2w), &ns2, -1.0, &dn2w);

        if self.ray > 0.0 {
            ns.reverse();
            dnw.reverse();
        }

        if self.choix % 2 != 0 {
            self.nplan.reverse();
        }

        tab_p[0] = self.pts;
        tab_p[nb_point - 1] = self.ptc;

        tab_v[0] = self.tgs;
        tab_v[nb_point - 1] = self.tgc;

        if nb_point > 2 {
            let cosa = ns.dot(&ns2);
            let sina = self.nplan.dot(&ns.crossed(&ns2));
            let mut angle = cosa.acos();
            if sina < 0.0 {
                angle = 2.0 * PI - angle;
            }
            let dangle = -(dnw.dot(&ns2) + ns.dot(&dn2w)) / sina;
            let ncrn = self.nplan.crossed(&ns);
            let dncrn = dnplan.crossed(&ns).added(&self.nplan.crossed(&dnw));

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

    /// `IsRational()` (`cxx:713-716`).
    pub(super) fn shape_is_rational(&self) -> bool {
        self.my_sshape == BlendFuncSectionShape::Rational
            || self.my_sshape == BlendFuncSectionShape::QuasiAngular
    }

    /// `GetSectionSize()` (`cxx:719-722`): the length of the maximum section.
    pub(super) fn section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// `GetMinimalWeight(Weigths)` (`cxx:727-731`). OCCT raises out of
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
        .expect("BlendFunc_CSCircular::GetMinimalWeight: section angle out of range");
    }

    /// `NbIntervals(S)` (`cxx:736-739`): `curv->NbIntervals(NextShape(S))`.
    pub(super) fn curve_nb_intervals(&self, s: Shape) -> i32 {
        self.curv.nb_intervals(next_shape(s) as u8)
    }

    /// `Intervals(T, S)` (`cxx:743-746`): `curv->Intervals(T, NextShape(S))`.
    pub(super) fn curve_intervals(&self, t: &mut [f64], s: Shape) {
        let src = self.curv.parameter_intervals(next_shape(s) as u8);
        let count = t.len().min(src.len());
        t[..count].copy_from_slice(&src[..count]);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:749-753`).
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
    /// (`cxx:757-771`). `GeomFill::GetTolerance` raises for an angle out of
    /// range in OCCT; the port panics at the same place. The radius is passed
    /// **signed** here (`cxx:766`).
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
            self.ray,
            angle_tol,
            surf_tol,
        )
        .expect("BlendFunc_CSCircular::GetTolerance: section angle out of range");
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

    /// `Knots(TKnots)` (`cxx:774-777`).
    pub(super) fn section_knots(&self, tknots: &mut [f64]) {
        section_convert::knots(self.my_tconv, tknots);
    }

    /// `Mults(TMults)` (`cxx:782-785`).
    pub(super) fn section_mults(&self, tmults: &mut [i32]) {
        section_convert::mults(self.my_tconv, tmults);
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:788-842`): the section
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

        self.set_param(prm);

        let (u1, v1) = p
            .parameters_on_s()
            .expect("BlendFunc_CSCircular::Section: point is not on the surface");

        let (pts, d1u1, d1v1) = self.surf.d1(u1, v1);
        self.pts = pts;
        self.ptc = self.curv.value(self.prmc);

        poles2d[0] = GpPnt2d::new(u1, v1);

        // Linear case (`cxx:812-821`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptc;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            return;
        }

        let ns = d1u1.crossed(&d1v1);
        let norm = self.nplan.crossed(&ns).magnitude();
        let mut ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );

        let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));
        let ns2 = GpVec::from_pnts(&center, &self.ptc).normalized();

        if self.ray > 0.0 {
            ns.reverse();
        }

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
    /// (`cxx:845-1026`): the first/last section with first derivatives.
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

        self.set_param(prm);
        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let (ptc, d1, _d2) = self.curv.d2(self.prmc);
        self.ptc = ptc;

        let (u, v) = p
            .parameters_on_s()
            .expect("BlendFunc_CSCircular::Section: point is not on the surface");
        let mut sol = MathVector::new(1, 2);
        sol.set_value(1, u);
        sol.set_value(2, v);

        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(u, v);
        self.pts = pts;

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        let ns2 = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);
        let temp = GpVec::from_pnts(&self.ptc, &self.pts);
        let mut secmember = MathVector::new(1, 2);

        secmember.set_value(
            1,
            self.dprmc * self.nplan.dot(&d1) - dnplan.dot(&temp),
        );

        // Derivative of n1 in relation to w (`cxx:888-894`).
        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let mut dnw = GpVec::zero().set_linear_form_3(
            (dnplan.dot(&ns) - grosterme * ndotns) / norm,
            &self.nplan,
            ndotns / norm,
            &dnplan,
            grosterme / norm,
            &ns,
        );

        let temp = dnw
            .multiplied_scalar(self.ray)
            .subtracted(&d1.multiplied_scalar(self.dprmc));
        let resul = ns2
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));

        secmember.set_value(2, -2.0 * resul.dot(&temp));

        self.evaluate_values(&sol);
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
            self.tgc = d1.multiplied_scalar(self.dprmc);

            // Derivative of n1 in relation to u (`cxx:913-921`).
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

            // Derivative of n1 in relation to v (`cxx:923-931`).
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
            // `ns = ns2` (`cxx:934`).
            let mut ns = ns2;

            let dn2w = GpVec::zero().set_linear_form_2plus(self.ray, &dnw, -1.0, &self.tgc, &self.tgs);
            let norm = resul.magnitude();
            let dn2w = dn2w.divide(norm);
            let ns2 = resul.normalized().reversed();
            let dn2w = GpVec::zero().set_linear_form_2(ns2.dot(&dn2w), &ns2, -1.0, &dn2w);

            // Poles 2d (`cxx:952-956`).
            poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));
            dpoles2d[0] = GpVec2d::new(secmember.value(1), secmember.value(2));

            // Linear case (`cxx:959-972`).
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

            // Case of the circle (`cxx:974-1017`).
            let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));
            let tgct = self.tgs.added(&dnw.multiplied_scalar(self.ray));

            if self.ray > 0.0 {
                ns.reverse();
                dnw.reverse();
            }

            if self.choix % 2 != 0 {
                self.nplan.reverse();
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
        // derivatives (`cxx:943-947`).
        let ns = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);
        let ns2 = resul.normalized().reversed();

        // Poles 2d (`cxx:952-956`).
        poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));

        // Linear case (`cxx:959-972`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptc;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            return false;
        }

        let mut ns = ns;
        let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));

        if self.ray > 0.0 {
            ns.reverse();
        }

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
