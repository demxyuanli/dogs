//! Section generation of `BlendFunc_EvolRad` — the analytic circle, the
//! section without derivatives, the first-order section, and the
//! shape/tolerance/interval queries the approximation driver asks for.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_EvolRad.cxx:39-98`
//! (`FusionneIntervalles`), `:1136-1158` (`Section(Param, ...)`),
//! `:1262-1279` / `:1332-1372` (`IsRational` ... `Mults`), `:1376-1450`
//! (`Section(P, Poles, ...)`), `:1452-1634` (`Section(P, ..., DPoles, ...)`)
//! (OCCT 8.0.0). The second-order section is in [`super::section_order2`].

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
use super::EvolRad;

/// `gp_Ax2(Center, np, ns)` from two vectors: OCCT raises when they are
/// parallel, and the port panics at the same place.
fn ax2_from_vecs(center: &GpPnt, z: &GpVec, x: &GpVec) -> GpAx2 {
    let z_dir = GpDir::from_vec(z).expect("BlendFunc_EvolRad::Section: degenerate axis");
    let x_dir = GpDir::from_vec(x).expect("BlendFunc_EvolRad::Section: degenerate axis");
    GpAx2::new(*center, z_dir, x_dir)
        .expect("BlendFunc_EvolRad::Section: gp_Ax2 construction error")
}

/// `FusionneIntervalles(I1, I2, Seq)` (`cxx:39-98`): the sorted union of two
/// parameter lists, collapsing the pairs closer than `PConfusion() * 0.99`
/// into their mean. OCCT appends into a `NCollection_Sequence`; the port
/// answers the sequence as a `Vec`.
fn fuse_intervals(i1: &[f64], i2: &[f64]) -> Vec<f64> {
    let epspar = PCONFUSION * 0.99;
    let mut seq = Vec::with_capacity(i1.len() + i2.len());
    let mut ind1 = 0usize;
    let mut ind2 = 0usize;

    while ind1 < i1.len() && ind2 < i2.len() {
        let v1 = i1[ind1];
        let v2 = i2[ind2];
        if (v1 - v2).abs() <= epspar {
            // Here elements of I1 and I2 are suitable (`cxx:58-63`).
            seq.push((v1 + v2) / 2.0);
            ind1 += 1;
            ind2 += 1;
        } else if v1 < v2 {
            // Here the element of I1 is suitable (`cxx:64-69`).
            seq.push(v1);
            ind1 += 1;
        } else {
            // Here the element of TABLE2 is suitable (`cxx:70-75`).
            seq.push(v2);
            ind2 += 1;
        }
    }

    if ind1 >= i1.len() {
        // Here I1 is empty, to be completed with the end of TABLE2
        // (`cxx:78-86`).
        while ind2 < i2.len() {
            seq.push(i2[ind2]);
            ind2 += 1;
        }
    }

    if ind2 >= i2.len() {
        // Here I2 is empty, to be completed with the end of I1 (`cxx:88-96`).
        while ind1 < i1.len() {
            seq.push(i1[ind1]);
            ind1 += 1;
        }
    }

    seq
}

impl EvolRad {
    /// `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)` (`cxx:1136-1193`): the
    /// cross-section circle at the guide parameter `param`, for the solution
    /// `(U1, V1, U2, V2)`. Used for a quick approximate visualization.
    ///
    /// OCCT ignores the `ComputeValues` result here (`cxx:1158`), so a failed
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
            // Unsatisfactory, but it is not necessary to stop (`cxx:1163`).
            norm1 = 1.0;
        }
        ns1 = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns1) / norm1,
            &self.nplan,
            -1.0 / norm1,
            &ns1,
        );
        let center = GpPnt::from_xyz(&self.pts1.coord.added(&ns1.coord.multiply(self.sg1 * self.ray)));

        // ns1 is oriented from the center to pts1 (`cxx:1170-1173`).
        if self.sg1 > 0.0 {
            ns1.reverse();
        }
        if self.choix % 2 != 0 {
            np.reverse();
        }
        let mut c = GpCirc::new(GpAx2::standard(), self.ray.abs());
        c.set_radius(self.ray.abs())
            .expect("BlendFunc_EvolRad::Section: invalid radius");
        c.set_position(&ax2_from_vecs(&center, &np, &ns1));
        let pdeb = 0.0;
        let mut pfin = clib::circle_parameter(&c.position(), &self.pts2);

        // Test of negative and almost null angles : Single Case (`cxx:1182-1191`).
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

    /// `IsRational()` (`cxx:1262-1266`).
    pub(super) fn shape_is_rational(&self) -> bool {
        self.my_sshape == BlendFuncSectionShape::Rational
            || self.my_sshape == BlendFuncSectionShape::QuasiAngular
    }

    /// `GetSectionSize()` (`cxx:1269-1272`): `lengthmax`, the length of the
    /// maximum section.
    pub(super) fn section_size(&self) -> f64 {
        self.lengthmax
    }

    /// `GetMinimalWeight(Weights)` (`cxx:1276-1279`). OCCT raises out of
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
        .expect("BlendFunc_EvolRad::GetMinimalWeight: section angle out of range");
    }

    /// `NbIntervals(S)` (`cxx:1283-1303`): the curve intervals refined by the
    /// law intervals when the law is not a single interval.
    pub(super) fn guide_nb_intervals(&self, s: Shape) -> i32 {
        let nb_int_courbe = self.curv.nb_intervals(next_shape(s) as u8);
        let nb_int_loi = self
            .fevol
            .nb_intervals(s)
            .expect("BlendFunc_EvolRad::NbIntervals: fevol->NbIntervals failed");
        if nb_int_loi == 1 {
            return nb_int_courbe;
        }

        let int_c = self.curv.parameter_intervals(next_shape(s) as u8);
        let int_l = self
            .fevol
            .intervals(s)
            .expect("BlendFunc_EvolRad::NbIntervals: fevol->Intervals failed");
        let inter = fuse_intervals(&int_c, &int_l);
        inter.len() as i32 - 1
    }

    /// `Intervals(T, S)` (`cxx:1307-1330`).
    pub(super) fn guide_intervals(&self, t: &mut [f64], s: Shape) {
        let nb_int_loi = self
            .fevol
            .nb_intervals(s)
            .expect("BlendFunc_EvolRad::Intervals: fevol->NbIntervals failed");
        if nb_int_loi == 1 {
            let src = self.curv.parameter_intervals(next_shape(s) as u8);
            let count = t.len().min(src.len());
            t[..count].copy_from_slice(&src[..count]);
            return;
        }

        let int_c = self.curv.parameter_intervals(next_shape(s) as u8);
        let int_l = self
            .fevol
            .intervals(s)
            .expect("BlendFunc_EvolRad::Intervals: fevol->Intervals failed");
        let inter = fuse_intervals(&int_c, &int_l);
        // `T(ii) = Inter(ii)` (`cxx:1325-1328`); a `T` that is too short
        // raises `Standard_OutOfRange` in OCCT.
        for (ii, val) in inter.iter().enumerate() {
            t[ii] = *val;
        }
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:1334-1338`).
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
    /// (`cxx:1346-1358`). `rayon = lengthmin / maxang` — "a radius is
    /// subtracted" (`cxx:1352`), not `|ray|`. `GeomFill::GetTolerance` raises
    /// for an angle out of range in OCCT; the port panics at the same place.
    pub(super) fn fill_section_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut MathVector,
        tol1d: &mut MathVector,
    ) {
        let rayon = self.lengthmin / self.maxang;
        let tol = section_convert::get_tolerance(
            self.my_tconv,
            self.maxang,
            rayon,
            angle_tol,
            surf_tol,
        )
        .expect("BlendFunc_EvolRad::GetTolerance: section angle out of range");
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

    /// `Knots(TKnots)` (`cxx:1362-1365`).
    pub(super) fn section_knots(&self, tknots: &mut [f64]) {
        section_convert::knots(self.my_tconv, tknots);
    }

    /// `Mults(TMults)` (`cxx:1369-1372`).
    pub(super) fn section_mults(&self, tmults: &mut [i32]) {
        section_convert::mults(self.my_tconv, tmults);
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:1376-1450`): the section
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
            .expect("BlendFunc_EvolRad::Section: point is not on the first surface");
        let (u2, v2) = p
            .parameters_on_s2()
            .expect("BlendFunc_EvolRad::Section: point is not on the second surface");
        x.set_value(1, u1);
        x.set_value(2, v1);
        x.set_value(3, u2);
        x.set_value(4, v2);

        // Calculation and storage of distmin (`cxx:1395-1396`).
        self.compute_values(&x, 0, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(&self.pts2));

        // ns1, ns2, np are copied locally to avoid crashing the fields
        // (`cxx:1398-1400`).
        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let mut np = self.nplan;

        poles2d[0] = GpPnt2d::new(x.value(1), x.value(2));
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(x.value(3), x.value(4));

        // Linear case (`cxx:1406-1413`).
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
            // Unsatisfactory, but it is not necessary to stop (`cxx:1418`).
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

        let center = GpPnt::from_xyz(&self.pts1.coord.added(&ns1.coord.multiply(self.sg1 * self.ray)));

        // ns1 (resp. ns2) is oriented from the center to pts1 (resp. pts2),
        // and the trihedron ns1, ns2, nplan is made direct (`cxx:1436-1450`).
        if self.sg1 > 0.0 {
            ns1.reverse();
        }
        if self.sg2 > 0.0 {
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
            self.ray.abs(),
            &center,
            poles,
            weights,
        );
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:1452-1634`): the first/last section with first derivatives.
    /// Returns whether the derivatives were computed.
    ///
    /// `nplan` / `dnplan` are copied locally: the `choix % 2 != 0` arm reverses
    /// the copies, not the fields (`cxx:1595-1600`).
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
            .expect("BlendFunc_EvolRad::Section: point is not on the first surface");
        let (u2, v2) = p
            .parameters_on_s2()
            .expect("BlendFunc_EvolRad::Section: point is not on the second surface");
        let mut sol = MathVector::new(1, 4);
        sol.set_value(1, u1);
        sol.set_value(2, v1);
        sol.set_value(3, u2);
        sol.set_value(4, v2);

        // Calculation of equations (`cxx:1480-1481`).
        self.compute_values(&sol, 1, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(&self.pts2));

        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let mut np = self.nplan;
        let mut dnp = self.dnplan;
        let mut rayprim = self.dray;
        let mut dnorm1w = GpVec::zero();
        let mut dnorm2w = GpVec::zero();

        let mut secmember = MathVector::new(1, 4);
        let mut istgt = true;

        // Calculation of derived Normal processing (`cxx:1491-1511`).
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

        // Tops 2D (`cxx:1523-1532`).
        poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));
        let last2d = poles2d.len() - 1;
        poles2d[last2d] = GpPnt2d::new(sol.value(3), sol.value(4));
        if !istgt {
            dpoles2d[0] = GpVec2d::new(secmember.value(1), secmember.value(2));
            dpoles2d[last2d] = GpVec2d::new(secmember.value(3), secmember.value(4));
        }

        // The linear case is processed (`cxx:1534-1547`).
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

        // Case of the circle (`cxx:1550-1634`).
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

        let center = GpPnt::from_xyz(&self.pts1.coord.added(&ns1.coord.multiply(self.sg1 * self.ray)));
        let tgc = if !istgt {
            // `tgc = sg1*ray*dnorm1w + sg1*dray*ns1 + tg1` (`cxx:1574`).
            GpVec::zero().set_linear_form_3(
                self.sg1 * self.ray,
                &dnorm1w,
                self.sg1 * self.dray,
                &ns1,
                1.0,
                &self.tg1,
            )
        } else {
            GpVec::zero()
        };

        // ns1 is oriented from the center to pts1, ns2 from the center to
        // pts2, and the trihedron ns1, ns2, nplan is made direct
        // (`cxx:1580-1600`).
        if self.sg1 > 0.0 {
            ns1.reverse();
            if !istgt {
                dnorm1w.reverse();
            }
        }
        if self.sg2 > 0.0 {
            ns2.reverse();
            if !istgt {
                dnorm2w.reverse();
            }
        }
        if self.choix % 2 != 0 {
            np.reverse();
            dnp.reverse();
        }

        // To avoid `std::abs(dray)` a few lines below (`cxx:1602-1605`).
        if self.ray < 0.0 {
            rayprim = -rayprim;
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
                self.ray.abs(),
                rayprim,
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
                self.ray.abs(),
                &center,
                poles,
                weights,
            );
            false
        }
    }
}
