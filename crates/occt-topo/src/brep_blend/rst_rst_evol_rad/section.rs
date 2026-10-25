//! Section generation of `BRepBlend_RstRstEvolRad`, plus its tolerances,
//! bounds, intervals and knots.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstEvolRad.cxx:613-655`
//! (`Section(Param, ...)`), `:773-836` (`Section(P, ...)`), `:839-1059`
//! (`Section(P, ..., DPoles, ...)`), `:1063-1074` (the second-derivative
//! `Section`, which returns `false`), `:1077-1092` (`Resolution`),
//! `:236-251` (`GetTolerance` / `GetBounds`), `:658-770` (`IsRational`,
//! `GetSectionSize`, `GetMinimalWeight`, `NbIntervals`, `Intervals`, `GetShape`,
//! `GetTolerance(BoundTol, ...)`, `Knots`, `Mults`) (OCCT 8.0.0).

use std::f64::consts::PI;

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_geom::geom_fill::section_convert;
use occt_math::{svd_padded, svd_solve_occt, Gauss, MathMatrix, MathVector};

use crate::blend::{
    get_minimal_weights, get_shape as blend_get_shape, next_shape, BlendFuncSectionShape,
    BlendPoint,
};

use super::equations::t3dto2d;
use super::SurfRstRstEvolRad;

/// `gp_Ax2(Center, np, ns)` from two vectors: OCCT raises when they are
/// parallel, and the port panics at the same place.
fn ax2_from_vecs(center: &GpPnt, z: &GpVec, x: &GpVec) -> GpAx2 {
    let z_dir = GpDir::from_vec(z).expect("BRepBlend_RstRstEvolRad::Section: degenerate axis");
    let x_dir = GpDir::from_vec(x).expect("BRepBlend_RstRstEvolRad::Section: degenerate axis");
    GpAx2::new(*center, z_dir, x_dir).expect("BRepBlend_RstRstEvolRad::Section: gp_Ax2 error")
}

impl SurfRstRstEvolRad {
    /// `Section(Param, U, V, Pdeb, Pfin, C)` (`cxx:613-655`): the cross-section
    /// circle at the guide parameter `Param`, through `cons1(U)` and
    /// `cons2(V)`. `U` is the parameter on `rst1`, `V` the one on `rst2`.
    pub fn section_circle(&mut self, param: f64, u: f64, v: f64) -> (f64, f64, GpCirc) {
        let (ptgui, d1gui) = self.tguide.d1(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.ray = self.tevol.value(param);
        let mut np = self.d1gui.normalized();
        self.ptrst1 = self.cons1.d0(u);
        self.ptrst2 = self.cons2.d0(v);

        // OCCT ignores the failure here (`cxx:629`).
        let (center, _not_used) = self.center_circle_rst1_rst2(&self.ptrst1, &self.ptrst2, &np);
        let center = center.unwrap_or_else(GpPnt::zero);

        let mut c = GpCirc::new(GpAx2::standard(), self.ray.abs());
        c.set_radius(self.ray.abs())
            .expect("BRepBlend_RstRstEvolRad::Section: invalid radius");
        let ns = GpVec::from_pnts(&center, &self.ptrst1).normalized();

        if self.choix % 2 != 0 {
            np.reverse();
        }
        c.set_position(&ax2_from_vecs(&center, &np, &ns));
        let pdeb = 0.0; // `ElCLib::Parameter(C, pts)` is commented out in OCCT.
        let mut pfin = clib::circle_parameter(&c.position(), &self.ptrst2);

        // Test of angles negative and almost null : Special Case
        // (`cxx:643-652`).
        if pfin > 1.5 * PI {
            np.reverse();
            c.set_position(&ax2_from_vecs(&center, &np, &ns));
            pfin = clib::circle_parameter(&c.position(), &self.ptrst2);
        }
        if pfin < PCONFUSION {
            pfin += PCONFUSION;
        }
        (pdeb, pfin, c)
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:839-1059`): the first/last section with first derivatives.
    /// Returns whether the derivatives were computed.
    ///
    /// As in OCCT, `nplan` is the *member*: the `choix % 2 != 0` arm reverses
    /// it in place (`cxx:1017-1021`), so the reversal outlives the call.
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
        let (ray, dray) = self.tevol.d1(prm);
        self.ray = ray;
        self.dray = dray;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        let mut dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let mut sol = MathVector::new(1, 2);
        self.prmrst1 = p
            .parameter_on_c1()
            .expect("BRepBlend_RstRstEvolRad::Section: point is not on a curve");
        sol.set_value(1, self.prmrst1);
        self.prmrst2 = p
            .parameter_on_c2()
            .expect("BRepBlend_RstRstEvolRad::Section: point is not on a curve");
        sol.set_value(2, self.prmrst2);
        self.pt2drst1 = self.rst1.d0(self.prmrst1);
        self.pt2drst2 = self.rst2.d0(self.prmrst2);

        let mut valsol = MathVector::new(1, 2);
        let mut gradsol = MathMatrix::new(1, 2, 1, 2);
        self.compute_values(&sol, &mut valsol, &mut gradsol);

        let (ptrst1, d11) = self.cons1.d1(sol.value(1));
        self.ptrst1 = ptrst1;
        let (ptrst2, d21) = self.cons2.d1(sol.value(2));
        self.ptrst2 = ptrst2;

        let mut secmember = MathVector::new(1, 2);
        let temp = GpVec::from_pnts(&self.ptgui, &self.ptrst1);
        secmember.set_value(1, self.normtg - dnplan.dot(&temp));
        let temp = GpVec::from_pnts(&self.ptgui, &self.ptrst2);
        secmember.set_value(2, self.normtg - dnplan.dot(&temp));

        let resol = Gauss::new(&gradsol, super::GAUSS_SECTION_EPS);
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

        let rst1rst2 = GpVec::from_pnts(&self.ptrst1, &self.ptrst2);
        let (center, med) = self.center_circle_rst1_rst2(&self.ptrst1, &self.ptrst2, &self.nplan);
        let center = match center {
            Some(c) => c,
            None => return false,
        };

        let normmed = med.magnitude();
        let med = med.normalized();
        // OCCT normalises these two only after the derivative block below
        // (`cxx:921-974`), so the unnormalised radii enter `d1n1` / `d1n2`.
        let n1 = GpVec::from_pnts(&center, &self.ptrst1);
        let n2 = GpVec::from_pnts(&center, &self.ptrst2);

        // OCCT leaves `d1n1` / `d1n2` unset in the `istgt` arm and never reads
        // them there; the port starts them at zero.
        let mut d1n1 = GpVec::zero();
        let mut d1n2 = GpVec::zero();
        if !istgt {
            // `secmember` contains the derivatives of the parameters on the
            // curves with respect to `t` (`cxx:927-929`).
            self.tgrst1 = d11.multiplied_scalar(secmember.value(1));
            self.tgrst2 = d21.multiplied_scalar(secmember.value(2));

            let norm2 = rst1rst2.square_magnitude();
            let mut d1rst1rst2 = self.tgrst2.subtracted(&self.tgrst1);
            let mut dist = self.ray * self.ray - 0.25 * norm2;
            let invdray = self.dray / self.ray;

            if dist > super::DIST_EPS {
                let d1p1p2_cros_np = d1rst1rst2
                    .crossed(&self.nplan)
                    .added(&rst1rst2.crossed(&dnplan));
                // Derivative of the perpendicular bisector (`cxx:941-943`).
                let mut dmed = d1p1p2_cros_np
                    .subtracted(&med.multiplied_scalar(med.dot(&d1p1p2_cros_np)));
                dmed = dmed.divided(normmed);
                dist = dist.sqrt();
                let d1_dist =
                    (self.ray * self.dray - 0.25 * rst1rst2.dot(&d1rst1rst2)) / dist;

                if self.choix > 2 {
                    dmed.reverse();
                }

                // The derivative of the coefficient `Dist` is located in `dmed`
                // (`cxx:945-947`).
                let dmed = GpVec::zero().set_linear_form_2(dist, &dmed, d1_dist, &med);
                d1rst1rst2 = d1rst1rst2.multiplied_scalar(0.5);
                // Derivative of the normal to the curve in `P1` / `P2`
                // (`cxx:949-960`).
                d1n1 = d1rst1rst2
                    .added(&dmed)
                    .added(&n1.multiplied_scalar(invdray))
                    .divided(-self.ray);
                d1n2 = d1rst1rst2
                    .subtracted(&dmed)
                    .subtracted(&n2.multiplied_scalar(invdray))
                    .divided(self.ray);
            } else {
                d1rst1rst2 = d1rst1rst2.multiplied_scalar(0.5);
                d1n1 = d1rst1rst2
                    .added(&n1.multiplied_scalar(invdray))
                    .divided(-self.ray);
                d1n2 = d1rst1rst2
                    .subtracted(&n2.multiplied_scalar(invdray))
                    .divided(self.ray);
            }
        }

        let n1 = n1.normalized();
        let n2 = n2.normalized();

        // Tops 2d (`cxx:977-990`).
        poles2d[0] = GpPnt2d::new(self.pt2drst1.x(), self.pt2drst1.y());
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(self.pt2drst2.x(), self.pt2drst2.y());
        if !istgt {
            let (_, d1urst, d1vrst) = self.surf1.d1(self.pt2drst1.x(), self.pt2drst1.y());
            let mut a = 0.0;
            let mut b = 0.0;
            t3dto2d(&mut a, &mut b, &self.tgrst1, &d1urst, &d1vrst);
            dpoles2d[0] = GpVec2d::new(a, b);

            let (_, d1urst, d1vrst) = self.surf2.d1(self.pt2drst2.x(), self.pt2drst2.y());
            t3dto2d(&mut a, &mut b, &self.tgrst2, &d1urst, &d1vrst);
            dpoles2d[last2d] = GpVec2d::new(a, b);
        }

        // Linear case (`cxx:993-1008`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.ptrst1;
            poles[upp] = self.ptrst2;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            if !istgt {
                dpoles[low] = self.tgrst1;
                dpoles[upp] = self.tgrst2;
                dweights[low] = 0.0;
                dweights[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle (`cxx:1011-1058`); the tangent at the centre of
        // the circle is built before the `choix` reversal.
        let tgct = if !istgt {
            GpVec::zero().set_linear_form_3(
                -self.ray,
                &d1n1,
                -self.dray,
                &n1,
                1.0,
                &self.tgrst1,
            )
        } else {
            GpVec::zero()
        };

        if self.choix % 2 != 0 {
            self.nplan.reverse();
            dnplan.reverse();
        }

        if !istgt {
            section_convert::get_circle_d1(
                self.my_tconv,
                &n1,
                &n2,
                &d1n1,
                &d1n2,
                &self.nplan,
                &dnplan,
                &self.ptrst1,
                &self.ptrst2,
                &self.tgrst1,
                &self.tgrst2,
                self.ray.abs(),
                self.dray,
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
                &n1,
                &n2,
                &self.nplan,
                &self.ptrst1,
                &self.ptrst2,
                self.ray.abs(),
                &center,
                poles,
                weights,
            );
            false
        }
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:773-836`): the section
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
        self.ray = self.tevol.value(prm);
        self.nplan = self.d1gui.normalized();

        let u = p
            .parameter_on_c1()
            .expect("BRepBlend_RstRstEvolRad::Section: point is not on a curve");
        let v = p
            .parameter_on_c2()
            .expect("BRepBlend_RstRstEvolRad::Section: point is not on a curve");

        let pt2d1 = self.rst1.d0(u);
        let pt2d2 = self.rst2.d0(v);

        self.ptrst1 = self.cons1.d0(u);
        self.ptrst2 = self.cons2.d0(v);
        self.distmin = self.distmin.min(self.ptrst1.distance(&self.ptrst2));

        poles2d[0] = GpPnt2d::new(pt2d1.x(), pt2d1.y());
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(pt2d2.x(), pt2d2.y());

        // Linear case (`cxx:805-812`).
        if self.my_sshape == BlendFuncSectionShape::Linear {
            poles[low] = self.ptrst1;
            poles[upp] = self.ptrst2;
            weights[low] = 1.0;
            weights[upp] = 1.0;
            return;
        }

        // Calculate the center of the circle (`cxx:815`).
        let (center, _not_used) =
            self.center_circle_rst1_rst2(&self.ptrst1, &self.ptrst2, &self.nplan);
        let center = center.unwrap_or_else(GpPnt::zero);

        // Normals to the section with points (`cxx:818-819`).
        let ns = GpVec::from_pnts(&center, &self.ptrst1).normalized();
        let ns2 = GpVec::from_pnts(&center, &self.ptrst2).normalized();

        if self.choix % 2 != 0 {
            self.nplan.reverse();
        }

        section_convert::get_circle(
            self.my_tconv,
            &ns,
            &ns2,
            &self.nplan,
            &self.ptrst1,
            &self.ptrst2,
            self.ray.abs(),
            &center,
            poles,
            weights,
        );
    }

    /// `Section(P, ..., D2Poles, ...)` (`cxx:1063-1074`): `false`, no second
    /// derivatives.
    pub(super) fn no_second_derivative_section(&mut self) -> bool {
        false
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:1077-1092`).
    pub(super) fn compute_resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic2d == 1 {
            *tol_u = occt_geom::approx_same_parameter::u_resolution(self.surf1.as_ref(), tol);
            *tol_v = occt_geom::approx_same_parameter::v_resolution(self.surf1.as_ref(), tol);
        } else {
            *tol_u = occt_geom::approx_same_parameter::u_resolution(self.surf2.as_ref(), tol);
            *tol_v = occt_geom::approx_same_parameter::v_resolution(self.surf2.as_ref(), tol);
        }
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:236-240`).
    pub(super) fn fill_parametric_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, self.cons1.resolution(tol));
        tolerance.set_value(2, self.cons2.resolution(tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:244-251`): no widening here, unlike
    /// the surface/restriction siblings.
    pub(super) fn fill_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        inf_bound.set_value(1, self.cons1.first_parameter());
        inf_bound.set_value(2, self.cons2.first_parameter());
        sup_bound.set_value(1, self.cons1.last_parameter());
        sup_bound.set_value(2, self.cons2.last_parameter());
    }

    /// `IsRational()` (`cxx:658-661`).
    pub(super) fn shape_is_rational(&self) -> bool {
        self.my_sshape == BlendFuncSectionShape::Rational
            || self.my_sshape == BlendFuncSectionShape::QuasiAngular
    }

    /// `GetSectionSize()` (`cxx:665-668`): the length of the maximum section,
    /// using the radius the law gave at the last `Set` / `Section`.
    pub(super) fn section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:672-676`). OCCT raises out of
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
        .expect("BRepBlend_RstRstEvolRad::GetMinimalWeight: section angle out of range");
    }

    /// `NbIntervals(S)` (`cxx:680-698`): `guide->NbIntervals(NextShape(S))`
    /// merged with `fevol->NbIntervals(S)`; when the law has a single interval
    /// the guide's count is returned unchanged.
    pub(super) fn nb_intervals_impl(&self, s: Shape) -> i32 {
        let nb_curve = self.guide.nb_intervals(next_shape(s) as u8);
        let nb_law = self
            .fevol
            .nb_intervals(s)
            .expect("BRepBlend_RstRstEvolRad::NbIntervals: law NbIntervals failed");

        if nb_law == 1 {
            return nb_curve;
        }

        let int_c = self.guide.parameter_intervals(next_shape(s) as u8);
        let int_l = self
            .fevol
            .intervals(s)
            .expect("BRepBlend_RstRstEvolRad::NbIntervals: law Intervals failed");
        let mut inter = Vec::new();
        super::queries::fusionne_intervalles(&int_c, &int_l, &mut inter);
        inter.len() as i32 - 1
    }

    /// `Intervals(T, S)` (`cxx:703-726`).
    pub(super) fn intervals_impl(&self, t: &mut [f64], s: Shape) {
        let nb_law = self
            .fevol
            .nb_intervals(s)
            .expect("BRepBlend_RstRstEvolRad::Intervals: law NbIntervals failed");

        if nb_law == 1 {
            let src = self.guide.parameter_intervals(next_shape(s) as u8);
            let count = t.len().min(src.len());
            t[..count].copy_from_slice(&src[..count]);
        } else {
            let int_c = self.guide.parameter_intervals(next_shape(s) as u8);
            let int_l = self
                .fevol
                .intervals(s)
                .expect("BRepBlend_RstRstEvolRad::Intervals: law Intervals failed");
            let mut inter = Vec::new();
            super::queries::fusionne_intervalles(&int_c, &int_l, &mut inter);
            let count = t.len().min(inter.len());
            t[..count].copy_from_slice(&inter[..count]);
        }
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:731-735`).
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

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)` (`cxx:742-756`).
    /// `GeomFill::GetTolerance` raises for an angle out of range in OCCT; the
    /// port panics at the same place.
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
        .expect("BRepBlend_RstRstEvolRad::GetTolerance: section angle out of range");
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

    /// `Knots(TKnots)` (`cxx:759-762`).
    pub(super) fn section_knots(&self, tknots: &mut [f64]) {
        section_convert::knots(self.my_tconv, tknots);
    }

    /// `Mults(TMults)` (`cxx:766-769`).
    pub(super) fn section_mults(&self, tmults: &mut [i32]) {
        section_convert::mults(self.my_tconv, tmults);
    }
}
