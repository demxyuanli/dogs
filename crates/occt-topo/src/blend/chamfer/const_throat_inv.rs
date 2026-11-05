//! `BlendFunc_ConstThroatInv` — the inverse of [`super::ConstThroat`]: a
//! solution on a restriction (curve on surface) of one of the two supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstThroatInv.hxx:26-76`,
//! `BlendFunc_ConstThroatInv.cxx:22-137` (OCCT 8.0.0).
//!
//! The four equations are the same as `ConstThroat`, with `X(2)` the guide-line
//! parameter and the partner point read from `X(3), X(4)` (or the curve on
//! surface, depending on `first`).
//!
//! Translation notes:
//! - `Value` (`cxx:89-134`) recomputes the guide-line frame and stores
//!   `pts1` / `pts2`; `Derivatives` (`cxx:139-230`) then builds `temp1` /
//!   `temp2` / `tempmid` from those stored points *before* its own `D1` calls
//!   overwrite them. Both orders are kept.
//! - `sign1` / `sign2` are written by `Set(Throat, dis1, Choix)` and never read
//!   by any body of this class; they stay as `#[allow(dead_code)]` members.
//! - `csurf` is an `Option`; a null handle makes `Value` / `Derivatives`
//!   answer `false` instead of faulting.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendFuncInv;

use super::gen_chamf_inv::GenChamfInv;
use super::shared::pnt_vec;

/// `class BlendFunc_ConstThroatInv : public BlendFunc_GenChamfInv`
/// (`hxx:26-76`).
pub struct ConstThroatInv {
    /// The `BlendFunc_GenChamfInv` base sub-object.
    pub(super) gen: GenChamfInv,
    /// `Throat` (`hxx:53`).
    pub(super) throat: f64,
    /// `param` (`hxx:55`).
    pub(super) param: f64,
    /// `sign1` (`hxx:56`): written by `Set`, never read.
    #[allow(dead_code)]
    pub(super) sign1: f64,
    /// `sign2` (`hxx:57`): written by `Set`, never read.
    #[allow(dead_code)]
    pub(super) sign2: f64,
    /// `ptgui` (`hxx:59`).
    pub(super) ptgui: GpPnt,
    /// `nplan` (`hxx:60`).
    pub(super) nplan: GpVec,
    /// `normtg` (`hxx:61`).
    pub(super) normtg: f64,
    /// `theD` (`hxx:62`).
    pub(super) the_d: f64,
    /// `d1gui` (`hxx:63`).
    pub(super) d1gui: GpVec,
    /// `d2gui` (`hxx:64`).
    pub(super) d2gui: GpVec,
    /// `pts1` (`hxx:66`).
    pub(super) pts1: GpPnt,
    /// `pts2` (`hxx:67`).
    pub(super) pts2: GpPnt,
    /// `d1u1` (`hxx:68`).
    pub(super) d1u1: GpVec,
    /// `d1v1` (`hxx:69`).
    pub(super) d1v1: GpVec,
    /// `d1u2` (`hxx:70`).
    pub(super) d1u2: GpVec,
    /// `d1v2` (`hxx:71`).
    pub(super) d1v2: GpVec,
}

impl ConstThroatInv {
    /// `BlendFunc_ConstThroatInv(S1, S2, C)` (`cxx:22-33`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            gen: GenChamfInv::new(surf1, surf2, curv),
            throat: 0.0,
            param: 0.0,
            sign1: 0.0,
            sign2: 0.0,
            ptgui: GpPnt::zero(),
            nplan: GpVec::zero(),
            normtg: 0.0,
            the_d: 0.0,
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            pts1: GpPnt::zero(),
            pts2: GpPnt::zero(),
            d1u1: GpVec::zero(),
            d1v1: GpVec::zero(),
            d1u2: GpVec::zero(),
            d1v2: GpVec::zero(),
        }
    }

    /// `Set(theThroat, dis1, Choix)` (`cxx:37-73`): the "quadrant" only selects
    /// the two signs, which no body of this class reads.
    pub fn set_throat(&mut self, the_throat: f64, _dis1: f64, choix: i32) {
        self.throat = the_throat;
        self.gen.choix = choix;
        match choix {
            1 | 2 => {
                self.sign1 = -1.0;
                self.sign2 = -1.0;
            }
            3 | 4 => {
                self.sign1 = 1.0;
                self.sign2 = -1.0;
            }
            5 | 6 => {
                self.sign1 = 1.0;
                self.sign2 = 1.0;
            }
            7 | 8 => {
                self.sign1 = -1.0;
                self.sign2 = 1.0;
            }
            _ => {
                self.sign1 = -1.0;
                self.sign2 = -1.0;
            }
        }
    }

    /// `Throat` (`hxx:53`).
    pub fn throat(&self) -> f64 {
        self.throat
    }

    /// `Set(OnFirst, COnSurf)` (`GenChamfInv.cxx:45-49`).
    pub fn set_first(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.gen.set_first(on_first, c_on_surf);
    }

    /// `Value(X, F)` (`cxx:89-134`).
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let Some(csurf) = self.gen.csurf.clone() else {
            return false;
        };
        let (p2d, _v2d) = csurf.d1(x.value(1));
        self.param = x.value(2);
        let (ptgui, d1gui, d2gui) = self.gen.curv.d2(self.param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        self.the_d = -self.nplan.dot(&pnt_vec(&self.ptgui));

        let xx = self.gen.partner_parameters(&p2d, x);
        self.pts1 = self.gen.surf1.value(xx.0, xx.1);
        self.pts2 = self.gen.surf2.value(xx.2, xx.3);

        f.set_value(1, self.nplan.dot(&pnt_vec(&self.pts1)) + self.the_d);
        f.set_value(2, self.nplan.dot(&pnt_vec(&self.pts2)) + self.the_d);

        let mid = super::shared::mid_pnt(&self.pts1, &self.pts2);
        let vmid = GpVec::from_pnts(&self.ptgui, &mid);
        f.set_value(3, vmid.square_magnitude() - self.throat * self.throat);

        let vref1 = GpVec::from_pnts(&self.ptgui, &self.pts1);
        let vref2 = GpVec::from_pnts(&self.ptgui, &self.pts2);
        f.set_value(4, vref1.square_magnitude() - vref2.square_magnitude());

        true
    }

    /// `Derivatives(X, D)` (`cxx:139-230`).
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let Some(csurf) = self.gen.csurf.clone() else {
            return false;
        };
        let (p2d, v2d) = csurf.d1(x.value(1));
        self.param = x.value(2);
        let (ptgui, d1gui, d2gui) = self.gen.curv.d2(self.param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        self.the_d = -self.nplan.dot(&pnt_vec(&self.ptgui));

        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let temp1 = GpVec::from_pnts(&self.ptgui, &self.pts1);
        let temp2 = GpVec::from_pnts(&self.ptgui, &self.pts2);
        let mid = super::shared::mid_pnt(&self.pts1, &self.pts2);
        let tempmid = GpVec::from_pnts(&self.ptgui, &mid);

        let xx = self.gen.partner_parameters(&p2d, x);
        let (p1, u1, v1) = self.gen.surf1.d1(xx.0, xx.1);
        self.pts1 = p1;
        self.d1u1 = u1;
        self.d1v1 = v1;
        let (p2, u2, v2) = self.gen.surf2.d1(xx.2, xx.3);
        self.pts2 = p2;
        self.d1u2 = u2;
        self.d1v2 = v2;

        let mid_after = super::shared::mid_pnt(&self.pts1, &self.pts2);

        if self.gen.first {
            let temp = GpVec::zero().set_linear_form_2(
                v2d.x(),
                &self.d1u1,
                v2d.y(),
                &self.d1v1,
            );

            d.set_value(1, 1, self.nplan.dot(&temp));
            d.set_value(2, 1, 0.0);
            d.set_value(3, 1, GpVec::from_pnts(&self.ptgui, &self.pts1).dot(&temp));
            d.set_value(4, 1, 2.0 * GpVec::from_pnts(&self.ptgui, &self.pts1).dot(&temp));

            d.set_value(1, 3, 0.0);
            d.set_value(1, 4, 0.0);
            d.set_value(2, 3, self.nplan.dot(&self.d1u2));
            d.set_value(2, 4, self.nplan.dot(&self.d1v2));
            d.set_value(3, 3, GpVec::from_pnts(&self.ptgui, &mid_after).dot(&self.d1u2));
            d.set_value(3, 4, GpVec::from_pnts(&self.ptgui, &mid_after).dot(&self.d1v2));
            d.set_value(4, 3, -2.0 * GpVec::from_pnts(&self.ptgui, &self.pts2).dot(&self.d1u2));
            d.set_value(4, 4, -2.0 * GpVec::from_pnts(&self.ptgui, &self.pts2).dot(&self.d1v2));
        } else {
            let temp = GpVec::zero().set_linear_form_2(
                v2d.x(),
                &self.d1u2,
                v2d.y(),
                &self.d1v2,
            );

            d.set_value(1, 1, 0.0);
            d.set_value(2, 1, self.nplan.dot(&temp));
            d.set_value(3, 1, GpVec::from_pnts(&self.ptgui, &self.pts2).dot(&temp));
            d.set_value(4, 1, -2.0 * GpVec::from_pnts(&self.ptgui, &self.pts2).dot(&temp));

            d.set_value(1, 3, self.nplan.dot(&self.d1u1));
            d.set_value(1, 4, self.nplan.dot(&self.d1v1));
            d.set_value(2, 3, 0.0);
            d.set_value(2, 4, 0.0);
            d.set_value(3, 3, GpVec::from_pnts(&self.ptgui, &mid_after).dot(&self.d1u1));
            d.set_value(3, 4, GpVec::from_pnts(&self.ptgui, &mid_after).dot(&self.d1v1));
            d.set_value(4, 3, 2.0 * GpVec::from_pnts(&self.ptgui, &self.pts1).dot(&self.d1u1));
            d.set_value(4, 4, 2.0 * GpVec::from_pnts(&self.ptgui, &self.pts1).dot(&self.d1v1));
        }

        d.set_value(1, 2, dnplan.dot(&temp1) - self.nplan.dot(&self.d1gui));
        d.set_value(2, 2, dnplan.dot(&temp2) - self.nplan.dot(&self.d1gui));
        d.set_value(3, 2, -2.0 * self.d1gui.dot(&tempmid));
        d.set_value(4, 2, 2.0 * self.d1gui.dot(&temp1) - 2.0 * self.d1gui.dot(&temp2));

        true
    }
}

impl MathFunctionSetWithDerivatives for ConstThroatInv {
    /// `NbVariables()` = `4` (`Blend_FuncInv.cxx:18-21`).
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`GenChamfInv.cxx:36-39`).
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:89-134`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:139-230`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d)
    }

    /// `Values(X, F, D)` (`GenChamfInv.cxx:119-122`): both are evaluated and
    /// the result is `true` regardless of `Value`'s answer.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.value(x, f);
        self.derivatives(x, d);
        true
    }
}

impl BlendFuncInv for ConstThroatInv {
    /// `Set(OnFirst, COnSurf)` (`GenChamfInv.cxx:45-49`).
    fn set(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.set_first(on_first, c_on_surf);
    }

    /// `GetTolerance(Tolerance, Tol)` (`GenChamfInv.cxx:52-65`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.gen.fill_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`GenChamfInv.cxx:69-114`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.gen.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:78-84`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 4);
        self.compute_value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= tol * tol
            && valsol.value(4).abs() <= tol * tol
    }
}
