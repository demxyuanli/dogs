//! `BlendFunc_ConstRad::ComputeValues` — the mandatory passage of every
//! equation, orders 0 and 1.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRad.cxx:136-432`
//! (OCCT 8.0.0). The `Order == 2` arm is in
//! [`super::compute_values_order2`].
//!
//! AUDIT (OCCT 8.0.1 `static` -> `thread_local` scratch): OCCT reused process-
//! static `gp_Vec` / `gp_Pnt` locals across calls; the port keeps those on
//! `ConstRad` (`ptgui`, `d1gui`, …) guarded by `my_x_order` / `my_t_order`.
//! No process-global scratch — no code change.

use occt_core::gp::{GpPnt2d, GpVec};
use occt_math::MathVector;

use crate::blend::blend_func;

use super::constants::EPS;
use super::{pnt_vec, ConstRad};

/// The order-0 quantities `ComputeValues` derives from `nplan`, `nsurf1` and
/// `nsurf2` once and reuses in the order-1 and order-2 arms
/// (`cxx:292-297`). OCCT keeps them in locals of the enclosing method; the
/// order-2 arm lives in its own file here, so they travel in one value.
pub(super) struct SectionPlaneTerms {
    /// `ncrossns1 = nplan.Crossed(nsurf1)` (`cxx:292`).
    pub ncrossns1: GpVec,
    /// `ncrossns2 = nplan.Crossed(nsurf2)` (`cxx:293`).
    pub ncrossns2: GpVec,
    /// `invnorm1` (`cxx:295-309`): `1 / |ncrossns1|`, `1` when singular.
    pub invnorm1: f64,
    /// `invnorm2` (`cxx:311-322`): `1 / |ncrossns2|`, `1` when singular.
    pub invnorm2: f64,
    /// `ndotns1` (`cxx:324`).
    pub ndotns1: f64,
    /// `ndotns2` (`cxx:325`).
    pub ndotns2: f64,
}

impl SectionPlaneTerms {
    /// `cxx:292-325`.
    pub(super) fn new(nplan: &GpVec, nsurf1: &GpVec, nsurf2: &GpVec) -> Self {
        let ncrossns1 = nplan.crossed(nsurf1);
        let ncrossns2 = nplan.crossed(nsurf2);

        let mut invnorm1 = ncrossns1.magnitude();
        if invnorm1 > EPS {
            invnorm1 = 1.0 / invnorm1;
        } else {
            // Unsatisfactory, but it is not necessary to crash (`cxx:305-308`).
            invnorm1 = 1.0;
        }
        let mut invnorm2 = ncrossns2.magnitude();
        if invnorm2 > EPS {
            invnorm2 = 1.0 / invnorm2;
        } else {
            invnorm2 = 1.0;
        }

        Self {
            ncrossns1,
            ncrossns2,
            invnorm1,
            invnorm2,
            ndotns1: nplan.dot(nsurf1),
            ndotns2: nplan.dot(nsurf2),
        }
    }
}

impl ConstRad {
    /// `ComputeValues(X, Order, ByParam, Param)` (`cxx:136-777`).
    ///
    /// `ByParam == false` evaluates the guide at the field `param` and caches
    /// the derivative order as `0` (`cxx:170-181`).
    pub(super) fn compute_values(
        &mut self,
        x: &MathVector,
        order: i32,
        by_param: bool,
        param: f64,
    ) -> bool {
        let mut t = param;
        // Case of implicite parameter (`cxx:149-153`).
        if !by_param {
            t = self.param;
        }

        // Is the work already done ? (`cxx:156-163`).
        let mut my_x_ok = order <= self.my_x_order;
        let mut ii = 1;
        while ii <= x.len() && my_x_ok {
            my_x_ok = x.value(ii) == self.x_val.value(ii);
            ii += 1;
        }

        let t_ok = t == self.tval && (order <= self.my_t_order || !by_param);

        if my_x_ok && t_ok {
            return true;
        }

        // Processing of t (`cxx:169-211`).
        if !t_ok {
            self.tval = t;
            self.my_t_order = if by_param { order } else { 0 };
            match self.my_t_order {
                0 => {
                    let (ptgui, d1gui) = self.tcurv.d1(t);
                    self.ptgui = ptgui;
                    self.d1gui = d1gui;
                    self.nplan = self.d1gui.normalized();
                }
                1 => {
                    let (ptgui, d1gui, d2gui) = self.tcurv.d2(t);
                    self.ptgui = ptgui;
                    self.d1gui = d1gui;
                    self.d2gui = d2gui;
                    self.nplan = self.d1gui.normalized();
                    self.invnormtg = 1.0 / self.d1gui.magnitude();
                    let nd = self.nplan.dot(&self.d2gui);
                    self.dnplan = GpVec::zero().set_linear_form_2(
                        self.invnormtg,
                        &self.d2gui,
                        -self.invnormtg * nd,
                        &self.nplan,
                    );
                }
                2 => {
                    let (ptgui, d1gui, d2gui, d3gui) = self.tcurv.d3(t);
                    self.ptgui = ptgui;
                    self.d1gui = d1gui;
                    self.d2gui = d2gui;
                    self.nplan = self.d1gui.normalized();
                    self.invnormtg = 1.0 / self.d1gui.magnitude();
                    let nd = self.nplan.dot(&self.d2gui);
                    self.dnplan = GpVec::zero().set_linear_form_2(
                        self.invnormtg,
                        &self.d2gui,
                        -self.invnormtg * nd,
                        &self.nplan,
                    );
                    let dinvnormtg = -nd * self.invnormtg * self.invnormtg;
                    self.d2nplan = GpVec::zero().set_linear_form_2(
                        self.invnormtg,
                        &d3gui,
                        dinvnormtg,
                        &self.d2gui,
                    );
                    let aux = dinvnormtg * nd
                        + self.invnormtg * (self.dnplan.dot(&self.d2gui) + self.nplan.dot(&d3gui));
                    self.d2nplan = GpVec::zero().set_linear_form_2(
                        -self.invnormtg * nd,
                        &self.dnplan,
                        -aux,
                        &self.nplan,
                    );
                }
                _ => return false,
            }
        }

        // Processing of X (`cxx:213-281`).
        if !my_x_ok {
            self.x_val.initialized(x);
            self.my_x_order = order;
            // -------------- Positioning on surfaces -----------------
            match self.my_x_order {
                0 => {
                    let (pts1, d1u1, d1v1) = self.surf1.d1(x.value(1), x.value(2));
                    self.pts1 = pts1;
                    self.d1u1 = d1u1;
                    self.d1v1 = d1v1;
                    self.nsurf1 = self.d1u1.crossed(&self.d1v1);
                    let (pts2, d1u2, d1v2) = self.surf2.d1(x.value(3), x.value(4));
                    self.pts2 = pts2;
                    self.d1u2 = d1u2;
                    self.d1v2 = d1v2;
                    self.nsurf2 = self.d1u2.crossed(&self.d1v2);
                }
                1 => {
                    let (pts1, d1u1, d1v1, d2u1, d2v1, d2uv1) =
                        self.surf1.d2(x.value(1), x.value(2));
                    self.pts1 = pts1;
                    self.d1u1 = d1u1;
                    self.d1v1 = d1v1;
                    self.d2u1 = d2u1;
                    self.d2v1 = d2v1;
                    self.d2uv1 = d2uv1;
                    self.nsurf1 = self.d1u1.crossed(&self.d1v1);
                    self.dns1u1 = self
                        .d2u1
                        .crossed(&self.d1v1)
                        .added(&self.d1u1.crossed(&self.d2uv1));
                    self.dns1v1 = self
                        .d2uv1
                        .crossed(&self.d1v1)
                        .added(&self.d1u1.crossed(&self.d2v1));
                    let (pts2, d1u2, d1v2, d2u2, d2v2, d2uv2) =
                        self.surf2.d2(x.value(3), x.value(4));
                    self.pts2 = pts2;
                    self.d1u2 = d1u2;
                    self.d1v2 = d1v2;
                    self.d2u2 = d2u2;
                    self.d2v2 = d2v2;
                    self.d2uv2 = d2uv2;
                    self.nsurf2 = self.d1u2.crossed(&self.d1v2);
                    self.dns1u2 = self
                        .d2u2
                        .crossed(&self.d1v2)
                        .added(&self.d1u2.crossed(&self.d2uv2));
                    self.dns1v2 = self
                        .d2uv2
                        .crossed(&self.d1v2)
                        .added(&self.d1u2.crossed(&self.d2v2));
                }
                2 => {
                    let sd1 = self
                        .surf1
                        .d3(x.value(1), x.value(2))
                        .expect("UNPORTED BlendFunc_ConstRad: surf1->D3 (cxx:240)");
                    self.pts1 = sd1.point;
                    self.d1u1 = sd1.d1u;
                    self.d1v1 = sd1.d1v;
                    self.d2u1 = sd1.d2u;
                    self.d2v1 = sd1.d2v;
                    self.d2uv1 = sd1.d2uv;
                    self.d3u1 = sd1.d3u;
                    self.d3v1 = sd1.d3v;
                    self.d3uuv1 = sd1.d3_uuv;
                    self.d3uvv1 = sd1.d3_uvv;
                    self.nsurf1 = self.d1u1.crossed(&self.d1v1);
                    self.dns1u1 = self
                        .d2u1
                        .crossed(&self.d1v1)
                        .added(&self.d1u1.crossed(&self.d2uv1));
                    self.dns1v1 = self
                        .d2uv1
                        .crossed(&self.d1v1)
                        .added(&self.d1u1.crossed(&self.d2v1));
                    let sd2 = self
                        .surf2
                        .d3(x.value(3), x.value(4))
                        .expect("UNPORTED BlendFunc_ConstRad: surf2->D3 (cxx:245)");
                    self.pts2 = sd2.point;
                    self.d1u2 = sd2.d1u;
                    self.d1v2 = sd2.d1v;
                    self.d2u2 = sd2.d2u;
                    self.d2v2 = sd2.d2v;
                    self.d2uv2 = sd2.d2uv;
                    self.d3u2 = sd2.d3u;
                    self.d3v2 = sd2.d3v;
                    self.d3uuv2 = sd2.d3_uuv;
                    self.d3uvv2 = sd2.d3_uvv;
                    self.nsurf2 = self.d1u2.crossed(&self.d1v2);
                    self.dns1u2 = self
                        .d2u2
                        .crossed(&self.d1v2)
                        .added(&self.d1u2.crossed(&self.d2uv2));
                    self.dns1v2 = self
                        .d2uv2
                        .crossed(&self.d1v2)
                        .added(&self.d1u2.crossed(&self.d2v2));
                }
                _ => return false,
            }

            // Case of degenerated surfaces (`cxx:254-280`).
            if self.nsurf1.magnitude() < EPS {
                let p2d = GpPnt2d::new(x.value(1), x.value(2));
                if order == 0 {
                    let normal = blend_func::compute_normal(self.surf1.as_ref(), &p2d)
                        .expect("BlendFunc_ConstRad: surf1->DN failed (cxx:261)");
                    if let Some(n) = normal {
                        self.nsurf1 = n;
                    }
                } else {
                    let answer = blend_func::compute_dnormal(self.surf1.as_ref(), &p2d)
                        .expect("BlendFunc_ConstRad: surf1->DN failed (cxx:265)");
                    if let Some((n, dnu, dnv)) = answer {
                        self.nsurf1 = n;
                        self.dns1u1 = dnu;
                        self.dns1v1 = dnv;
                    }
                }
            }
            if self.nsurf2.magnitude() < EPS {
                let p2d = GpPnt2d::new(x.value(3), x.value(4));
                if order == 0 {
                    let normal = blend_func::compute_normal(self.surf2.as_ref(), &p2d)
                        .expect("BlendFunc_ConstRad: surf2->DN failed (cxx:271)");
                    if let Some(n) = normal {
                        self.nsurf2 = n;
                    }
                } else {
                    let answer = blend_func::compute_dnormal(self.surf2.as_ref(), &p2d)
                        .expect("BlendFunc_ConstRad: surf2->DN failed (cxx:275)");
                    if let Some((n, dnu, dnv)) = answer {
                        self.nsurf2 = n;
                        self.dns1u2 = dnu;
                        self.dns1v2 = dnv;
                    }
                }
            }
        }

        // -------------------- Positioning of order 0 ---------------------
        // (`cxx:283-336`).
        let terms = SectionPlaneTerms::new(&self.nplan, &self.nsurf1, &self.nsurf2);
        let the_d = -self.nplan.dot(&pnt_vec(&self.ptgui));

        self.e.set_value(
            1,
            (self.nplan.x() * (self.pts1.x() + self.pts2.x())
                + self.nplan.y() * (self.pts1.y() + self.pts2.y())
                + self.nplan.z() * (self.pts1.z() + self.pts2.z()))
                / 2.0
                + the_d,
        );

        let temp = GpVec::zero().set_linear_form_2(terms.ndotns1, &self.nplan, -1.0, &self.nsurf1);
        let temp = temp.multiplied_scalar(terms.invnorm1);
        let mut resul = GpVec::zero().set_linear_form_2(
            self.ray1,
            &temp,
            1.0,
            &GpVec::from_pnts(&self.pts2, &self.pts1),
        );
        let temp2 = GpVec::zero()
            .set_linear_form_2(terms.ndotns2, &self.nplan, -1.0, &self.nsurf2)
            .multiplied_scalar(terms.invnorm2);
        resul = resul.subtracted(&temp2.multiplied_scalar(self.ray2));

        self.e.set_value(2, resul.x());
        self.e.set_value(3, resul.y());
        self.e.set_value(4, resul.z());

        // -------------------- Positioning of order 1 ---------------------
        if order >= 1 {
            self.compute_first_order(&terms, by_param);
        }

        // ------   Positioning of order 2  -----------------------------
        if order == 2 {
            self.compute_second_order(&terms, by_param);
        }

        true
    }

    /// `cxx:338-432`: `DEDX`, and `DEDT` when `ByParam`.
    fn compute_first_order(&mut self, terms: &SectionPlaneTerms, by_param: bool) {
        self.dedx
            .set_value(1, 1, self.nplan.dot(&self.d1u1) / 2.0);
        self.dedx
            .set_value(1, 2, self.nplan.dot(&self.d1v1) / 2.0);
        self.dedx
            .set_value(1, 3, self.nplan.dot(&self.d1u2) / 2.0);
        self.dedx
            .set_value(1, 4, self.nplan.dot(&self.d1v2) / 2.0);

        let mut cube = terms.invnorm1 * terms.invnorm1 * terms.invnorm1;
        // Derived in relation to u1 (`cxx:349-356`).
        let mut grosterme = -terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1u1))
            * cube;
        self.dndu1 = GpVec::zero().set_linear_form_3(
            grosterme * terms.ndotns1 + terms.invnorm1 * self.nplan.dot(&self.dns1u1),
            &self.nplan,
            -grosterme,
            &self.nsurf1,
            -terms.invnorm1,
            &self.dns1u1,
        );
        let mut resul = self.dndu1.multiplied_scalar(self.ray1).added(&self.d1u1);
        self.dedx.set_value(2, 1, resul.x());
        self.dedx.set_value(3, 1, resul.y());
        self.dedx.set_value(4, 1, resul.z());

        // Derived in relation to v1 (`cxx:363-372`).
        grosterme = -terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1v1))
            * cube;
        self.dndv1 = GpVec::zero().set_linear_form_3(
            grosterme * terms.ndotns1 + terms.invnorm1 * self.nplan.dot(&self.dns1v1),
            &self.nplan,
            -grosterme,
            &self.nsurf1,
            -terms.invnorm1,
            &self.dns1v1,
        );
        resul = self.dndv1.multiplied_scalar(self.ray1).added(&self.d1v1);
        self.dedx.set_value(2, 2, resul.x());
        self.dedx.set_value(3, 2, resul.y());
        self.dedx.set_value(4, 2, resul.z());

        cube = terms.invnorm2 * terms.invnorm2 * terms.invnorm2;
        // Derived in relation to u2 (`cxx:379-387`).
        grosterme = -terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1u2))
            * cube;
        self.dndu2 = GpVec::zero().set_linear_form_3(
            grosterme * terms.ndotns2 + terms.invnorm2 * self.nplan.dot(&self.dns1u2),
            &self.nplan,
            -grosterme,
            &self.nsurf2,
            -terms.invnorm2,
            &self.dns1u2,
        );
        resul = GpVec::zero().set_linear_form_2(-self.ray2, &self.dndu2, -1.0, &self.d1u2);
        self.dedx.set_value(2, 3, resul.x());
        self.dedx.set_value(3, 3, resul.y());
        self.dedx.set_value(4, 3, resul.z());

        // Derived in relation to v2 (`cxx:395-407`).
        grosterme = -terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1v2))
            * cube;
        self.dndv2 = GpVec::zero().set_linear_form_3(
            grosterme * terms.ndotns2 + terms.invnorm2 * self.nplan.dot(&self.dns1v2),
            &self.nplan,
            -grosterme,
            &self.nsurf2,
            -terms.invnorm2,
            &self.dns1v2,
        );
        resul = GpVec::zero().set_linear_form_2(-self.ray2, &self.dndv2, -1.0, &self.d1v2);
        self.dedx.set_value(2, 4, resul.x());
        self.dedx.set_value(3, 4, resul.y());
        self.dedx.set_value(4, 4, resul.z());

        if by_param {
            let temp = GpVec::from_xyz(
                &self
                    .pts1
                    .coord
                    .added(&self.pts2.coord)
                    .multiply(0.5)
                    .subtracted(&self.ptgui.coord),
            );
            // Derived from n1 in relation to w (`cxx:411-417`).
            grosterme = terms
                .ncrossns1
                .dot(&self.dnplan.crossed(&self.nsurf1))
                * terms.invnorm1
                * terms.invnorm1;
            self.dn1w = GpVec::zero().set_linear_form_3(
                (self.dnplan.dot(&self.nsurf1) - grosterme * terms.ndotns1) * terms.invnorm1,
                &self.nplan,
                terms.ndotns1 * terms.invnorm1,
                &self.dnplan,
                grosterme * terms.invnorm1,
                &self.nsurf1,
            );

            // Derivee from n2 in relation to w (`cxx:420-426`).
            grosterme = terms
                .ncrossns2
                .dot(&self.dnplan.crossed(&self.nsurf2))
                * terms.invnorm2
                * terms.invnorm2;
            self.dn2w = GpVec::zero().set_linear_form_3(
                (self.dnplan.dot(&self.nsurf2) - grosterme * terms.ndotns2) * terms.invnorm2,
                &self.nplan,
                terms.ndotns2 * terms.invnorm2,
                &self.dnplan,
                grosterme * terms.invnorm2,
                &self.nsurf2,
            );

            self.dedt
                .set_value(1, self.dnplan.dot(&temp) - 1.0 / self.invnormtg);
            self.dedt.set_value(2, self.ray1 * self.dn1w.x() - self.ray2 * self.dn2w.x());
            self.dedt.set_value(3, self.ray1 * self.dn1w.y() - self.ray2 * self.dn2w.y());
            self.dedt.set_value(4, self.ray1 * self.dn1w.z() - self.ray2 * self.dn2w.z());
        }
    }
}
