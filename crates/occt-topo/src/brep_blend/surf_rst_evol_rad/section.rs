//! Section generation of `BRepBlend_SurfRstEvolRad`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstEvolRad.cxx:676-734`
//! (`Section(Param, ...)`), `:839-1070` (`Section(P, ..., DPoles, ...)`),
//! `:1076-1085` (the second-derivative `Section`, which returns `false`),
//! `:1087-1149` (`Section(P, ...)`), `:1151-1166` (`Resolution`) (OCCT 8.0.0).

use std::f64::consts::PI;

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_geom::geom_fill::section_convert;
use occt_math::{svd_padded, svd_solve_occt, Gauss, MathMatrix, MathVector};

use crate::blend::BlendPoint;

use super::equations::t3dto2d;
use super::SurfRstEvolRad;

/// `gp_Ax2(Center, np, ns)` from two vectors: OCCT raises when they are
/// parallel, and the port panics at the same place.
fn ax2_from_vecs(center: &GpPnt, z: &GpVec, x: &GpVec) -> GpAx2 {
    let z_dir = GpDir::from_vec(z).expect("BRepBlend_SurfRstEvolRad::Section: degenerate axis");
    let x_dir = GpDir::from_vec(x).expect("BRepBlend_SurfRstEvolRad::Section: degenerate axis");
    GpAx2::new(*center, z_dir, x_dir)
        .expect("BRepBlend_SurfRstEvolRad::Section: gp_Ax2 construction error")
}

impl SurfRstEvolRad {
    /// `Section(Param, U, V, W, Pdeb, Pfin, C)` (`cxx:676-734`): the
    /// cross-section circle at the guide parameter `Param` for the solution
    /// `(U, V, W)`. The radius is re-read from the law, `ray = sg1 *
    /// tevol->Value(Param)` (`cxx:687`).
    pub fn section_circle(&mut self, param: f64, u: f64, v: f64, w: f64) -> (f64, f64, GpCirc) {
        let (ptgui, d1gui) = self.tguide.d1(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        let mut np = self.d1gui.normalized();
        self.ray = self.sg1 * self.tevol.value(param);

        let (pts, d1u1, d1v1) = self.surf.d1(u, v);
        self.pts = pts;
        self.ptrst = self.cons.d0(w);

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
            .expect("BRepBlend_SurfRstEvolRad::Section: invalid radius");

        if self.ray > 0.0 {
            ns.reverse();
        }
        if self.choix % 2 != 0 {
            np.reverse();
        }
        c.set_position(&ax2_from_vecs(&center, &np, &ns));
        let pdeb = 0.0; // `ElCLib::Parameter(C, pts)` is commented out in OCCT.
        let mut pfin = clib::circle_parameter(&c.position(), &self.ptrst);

        // Test negative and almost null angles : Special case (`cxx:713-723`).
        if pfin > 1.5 * PI {
            np.reverse();
            c.set_position(&ax2_from_vecs(&center, &np, &ns));
            pfin = clib::circle_parameter(&c.position(), &self.ptrst);
        }
        if pfin < PCONFUSION {
            pfin += PCONFUSION;
        }
        (pdeb, pfin, c)
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:839-1070`): the first/last section with first derivatives.
    /// Returns whether the derivatives were computed.
    ///
    /// As in OCCT, `nplan` is the *member*: the `choix % 2 != 0` arm reverses
    /// it in place (`cxx:1028-1033`), so the reversal outlives the call.
    #[allow(clippy::too_many_arguments)]
    pub fn section_with_derivatives(
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

        let (ptgui, d1gui, d2gui) = self.tguide.d2(prm);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;

        let (ray, a_dray) = self.tevol.d1(prm);
        self.ray = self.sg1 * ray;
        let a_dray = self.sg1 * a_dray;

        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let (u1, v1) = p
            .parameters_on_s()
            .expect("BRepBlend_SurfRstEvolRad::Section: point is not on a surface");
        let mut sol = MathVector::new(1, 3);
        sol.set_value(1, u1);
        sol.set_value(2, v1);
        self.prmrst = p
            .parameter_on_c()
            .expect("BRepBlend_SurfRstEvolRad::Section: point is not on a curve");
        sol.set_value(3, self.prmrst);
        self.pt2drst = self.rst.d0(self.prmrst);

        let mut valsol = MathVector::new(1, 3);
        let mut gradsol = MathMatrix::new(1, 3, 1, 3);
        self.compute_values(&sol, &mut valsol, &mut gradsol);

        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(sol.value(1), sol.value(2));
        self.pts = pts;
        let (ptrst, d1) = self.cons.d1(sol.value(3));
        self.ptrst = ptrst;

        let mut secmember = MathVector::new(1, 3);
        let temp = GpVec::from_pnts(&self.ptgui, &self.pts);
        secmember.set_value(1, self.normtg - dnplan.dot(&temp));
        let temp = GpVec::from_pnts(&self.ptgui, &self.ptrst);
        secmember.set_value(2, self.normtg - dnplan.dot(&temp));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let mut norm = ncrossns.magnitude();
        if norm < super::EPS {
            // Not enough, but it is not necessary to stop (`cxx:893-899`).
            norm = 1.0;
        }

        // Derivative of n1 corresponding to w (`cxx:901-909`).
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
        let resul = GpVec::zero().set_linear_form_2(
            self.ray,
            &temp,
            1.0,
            &GpVec::from_pnts(&self.ptrst, &self.pts),
        );
        secmember.set_value(
            3,
            -2.0 * self.ray * dnw.dot(&resul) - 2.0 * a_dray * temp.dot(&resul)
                + 2.0 * self.ray * a_dray,
        );

        let resol = Gauss::new(&gradsol, super::GAUSS_DEFAULT_EPS);
        let istgt;
        if resol.is_done() {
            istgt = false;
            if let Some(x) = resol.solve(&secmember) {
                secmember = x;
            }
        } else {
            let sing_rs = svd_padded(&gradsol);
            if sing_rs.is_done {
                let dedt = secmember.clone();
                svd_solve_occt(&sing_rs, &dedt, &mut secmember, super::SVD_EPS);
                istgt = false;
            } else {
                istgt = true;
            }
        }

        let mut ns = ns;
        let ns2;
        let mut dn2w = GpVec::zero();
        if !istgt {
            self.tgs = GpVec::zero().set_linear_form_2(
                secmember.value(1),
                &d1u1,
                secmember.value(2),
                &d1v1,
            );
            self.tgrst = d1.multiplied_scalar(secmember.value(3));

            // Derivative of n1 corresponding to u1 (`cxx:944-950`).
            let temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
            let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
            let resulu = GpVec::zero().set_linear_form_3(
                -(grosterme * ndotns - self.nplan.dot(&temp)) / norm,
                &self.nplan,
                grosterme / norm,
                &ns,
                -1.0 / norm,
                &temp,
            );

            // Derivative of n1 corresponding to v1 (`cxx:952-958`).
            let temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
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
            ns = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);

            dn2w = GpVec::zero().set_linear_form_3(
                self.ray,
                &dnw,
                -1.0,
                &self.tgrst,
                1.0,
                &self.tgs,
            );
            dn2w = GpVec::zero().set_linear_form_2(a_dray, &ns, 1.0, &dn2w);
            let norm2 = resul.magnitude();
            dn2w = dn2w.divided(norm2);
            let ns2_local = resul.normalized().reversed();
            dn2w = GpVec::zero().set_linear_form_2(ns2_local.dot(&dn2w), &ns2_local, -1.0, &dn2w);
            ns2 = ns2_local;
        } else {
            ns = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);
            ns2 = resul.normalized().reversed();
        }

        // Tops 2D (`cxx:982-996`).
        poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(self.pt2drst.x(), self.pt2drst.y());
        if !istgt {
            dpoles2d[0] = GpVec2d::new(secmember.value(1), secmember.value(2));
            let (_, d1urst, d1vrst) = self.surfrst.d1(self.pt2drst.x(), self.pt2drst.y());
            let mut a = 0.0;
            let mut b = 0.0;
            t3dto2d(&mut a, &mut b, &self.tgrst, &d1urst, &d1vrst);
            dpoles2d[last2d] = GpVec2d::new(a, b);
        }

        // Linear Case (`cxx:998-1011`).
        if self.my_sshape == crate::blend::BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptrst;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            if !istgt {
                dpoles[low] = self.tgs;
                dpoles[upp] = self.tgrst;
                dweights[low] = 0.0;
                dweights[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle (`cxx:1013-1069`).
        let center = GpPnt::from_xyz(&self.pts.coord.added(&ns.coord.multiply(self.ray)));
        let tgct = if !istgt {
            GpVec::zero().set_linear_form_3(self.ray, &dnw, a_dray, &ns, 1.0, &self.tgs)
        } else {
            GpVec::zero()
        };

        if self.ray > 0.0 {
            ns.reverse();
            if !istgt {
                dnw.reverse();
            }
        }
        let mut dnplan = dnplan;
        if self.choix % 2 != 0 {
            self.nplan.reverse();
            dnplan.reverse();
        }
        if !istgt {
            // `rayprim` is `-aDray` when `ray < 0`, else `aDray`
            // (`cxx:1034-1042`).
            let rayprim = if self.ray < 0.0 { -a_dray } else { a_dray };
            section_convert::get_circle_d1(
                self.my_tconv,
                &ns,
                &ns2,
                &dnw,
                &dn2w,
                &self.nplan,
                &dnplan,
                &self.pts,
                &self.ptrst,
                &self.tgs,
                &self.tgrst,
                self.ray.abs(),
                rayprim,
                &center,
                &tgct,
                poles,
                dpoles,
                weights,
                dweights,
            );
            true
        } else {
            section_convert::get_circle(
                self.my_tconv,
                &ns,
                &ns2,
                &self.nplan,
                &self.pts,
                &self.ptrst,
                self.ray.abs(),
                &center,
                poles,
                weights,
            );
            false
        }
    }

    /// `Section(P, Poles, DPoles, D2Poles, ...)` (`cxx:1076-1085`): `false`, no
    /// second derivatives.
    pub(super) fn no_second_derivative_section(&mut self) -> bool {
        false
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:1087-1149`): the section
    /// without derivatives.
    pub fn section_poles_impl(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        let prm = p.parameter();
        let low = 0usize;
        let upp = poles.len() - 1;

        let (ptgui, d1gui) = self.tguide.d1(prm);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.ray = self.sg1 * self.tevol.value(prm);
        self.nplan = self.d1gui.normalized();

        let (u1, v1) = p
            .parameters_on_s()
            .expect("BRepBlend_SurfRstEvolRad::Section: point is not on a surface");
        // `jlr : point on curve not on surface` (`cxx:1110`).
        let w = p
            .parameter_on_c()
            .expect("BRepBlend_SurfRstEvolRad::Section: point is not on a curve");
        let pt2d = self.rst.d0(w);

        let (pts, d1u1, d1v1) = self.surf.d1(u1, v1);
        self.pts = pts;
        self.ptrst = self.cons.d0(w);
        self.distmin = self.distmin.min(self.pts.distance(&self.ptrst));

        poles2d[0] = GpPnt2d::new(u1, v1);
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(pt2d.x(), pt2d.y());

        // Linear case (`cxx:1121-1129`).
        if self.my_sshape == crate::blend::BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptrst;
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

        let ns2 = GpVec::from_pnts(&center, &self.ptrst).normalized();
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
            &self.ptrst,
            self.ray.abs(),
            &center,
            poles,
            weights,
        );
    }
}

impl SurfRstEvolRad {
    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1151-1166`).
    pub(super) fn compute_resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic2d == 1 {
            *tol_u = occt_geom::approx_same_parameter::u_resolution(self.surf.as_ref(), tol);
            *tol_v = occt_geom::approx_same_parameter::v_resolution(self.surf.as_ref(), tol);
        } else {
            *tol_u = occt_geom::approx_same_parameter::u_resolution(self.surfrst.as_ref(), tol);
            *tol_v = occt_geom::approx_same_parameter::v_resolution(self.surfrst.as_ref(), tol);
        }
    }
}
