//! `BlendFunc_ConstThroatWithPenetrationInv` — the inverse of
//! [`super::ConstThroatWithPenetration`]: a solution on a restriction (curve on
//! surface) of one of the two supports.
//!
//! Source:
//! `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstThroatWithPenetrationInv.hxx:24-49`,
//! `BlendFunc_ConstThroatWithPenetrationInv.cxx:22-199` (OCCT 8.0.0).
//!
//! Only `IsSolution`, `Value` and `Derivatives` are overridden; the frame
//! update, `Set`, the tolerances and the bounds are inherited from
//! [`super::ConstThroatInv`].
//!
//! Translation notes:
//! - `temp1` / `temp2` / `temp3` are built from the points stored by the last
//!   `Value` call, before `Derivatives`'s own `D1` calls overwrite them
//!   (`cxx:111-114`); that order is kept.
//! - A null `csurf` makes `Value` / `Derivatives` answer `false`.

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendFuncInv;

use super::const_throat_inv::ConstThroatInv;
use super::shared::pnt_vec;

/// `class BlendFunc_ConstThroatWithPenetrationInv : public
/// BlendFunc_ConstThroatInv` (`hxx:24-49`).
pub struct ConstThroatWithPenetrationInv {
    /// The `BlendFunc_ConstThroatInv` base sub-object.
    pub(super) parent: ConstThroatInv,
}

impl ConstThroatWithPenetrationInv {
    /// `BlendFunc_ConstThroatWithPenetrationInv(S1, S2, C)` (`cxx:22-27`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            parent: ConstThroatInv::new(surf1, surf2, curv),
        }
    }

    /// `Set(theThroat, dis1, Choix)` (`ConstThroatInv.cxx:37-73`).
    pub fn set_throat(&mut self, the_throat: f64, dis1: f64, choix: i32) {
        self.parent.set_throat(the_throat, dis1, choix);
    }

    /// `Throat` (`ConstThroatInv.hxx:53`).
    pub fn throat(&self) -> f64 {
        self.parent.throat()
    }

    /// `Set(OnFirst, COnSurf)` (`GenChamfInv.cxx:45-49`).
    pub fn set_first(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.parent.set_first(on_first, c_on_surf);
    }

    /// `Value(X, F)` (`cxx:44-87`).
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let Some(csurf) = self.parent.gen.csurf.clone() else {
            return false;
        };
        let (p2d, _v2d) = csurf.d1(x.value(1));
        self.parent.param = x.value(2);
        let (ptgui, d1gui, d2gui) = self.parent.gen.curv.d2(self.parent.param);
        self.parent.ptgui = ptgui;
        self.parent.d1gui = d1gui;
        self.parent.d2gui = d2gui;
        self.parent.normtg = self.parent.d1gui.magnitude();
        self.parent.nplan = self.parent.d1gui.normalized();
        self.parent.the_d = -self.parent.nplan.dot(&pnt_vec(&self.parent.ptgui));

        let xx = self.parent.gen.partner_parameters(&p2d, x);
        self.parent.pts1 = self.parent.gen.surf1.value(xx.0, xx.1);
        self.parent.pts2 = self.parent.gen.surf2.value(xx.2, xx.3);

        f.set_value(
            1,
            self.parent.nplan.dot(&pnt_vec(&self.parent.pts1)) + self.parent.the_d,
        );
        f.set_value(
            2,
            self.parent.nplan.dot(&pnt_vec(&self.parent.pts2)) + self.parent.the_d,
        );

        let vref = GpVec::from_pnts(&self.parent.ptgui, &self.parent.pts1);
        f.set_value(
            3,
            vref.square_magnitude() - self.parent.throat * self.parent.throat,
        );

        let vec12 = GpVec::from_pnts(&self.parent.pts1, &self.parent.pts2);
        f.set_value(4, vref.dot(&vec12));

        true
    }

    /// `Derivatives(X, D)` (`cxx:92-199`).
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let Some(csurf) = self.parent.gen.csurf.clone() else {
            return false;
        };
        let (p2d, v2d) = csurf.d1(x.value(1));
        self.parent.param = x.value(2);
        let (ptgui, d1gui, d2gui) = self.parent.gen.curv.d2(self.parent.param);
        self.parent.ptgui = ptgui;
        self.parent.d1gui = d1gui;
        self.parent.d2gui = d2gui;
        self.parent.normtg = self.parent.d1gui.magnitude();
        self.parent.nplan = self.parent.d1gui.normalized();
        self.parent.the_d = -self.parent.nplan.dot(&pnt_vec(&self.parent.ptgui));

        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.parent.normtg,
            &self.parent.d2gui,
            -1.0 / self.parent.normtg * self.parent.nplan.dot(&self.parent.d2gui),
            &self.parent.nplan,
        );

        let temp1 = GpVec::from_pnts(&self.parent.ptgui, &self.parent.pts1);
        let temp2 = GpVec::from_pnts(&self.parent.ptgui, &self.parent.pts2);
        let temp3 = GpVec::from_pnts(&self.parent.pts1, &self.parent.pts2);

        let xx = self.parent.gen.partner_parameters(&p2d, x);
        let (p1, u1, v1) = self.parent.gen.surf1.d1(xx.0, xx.1);
        self.parent.pts1 = p1;
        self.parent.d1u1 = u1;
        self.parent.d1v1 = v1;
        let (p2, u2, v2) = self.parent.gen.surf2.d1(xx.2, xx.3);
        self.parent.pts2 = p2;
        self.parent.d1u2 = u2;
        self.parent.d1v2 = v2;

        if self.parent.gen.first {
            let temp = GpVec::zero().set_linear_form_2(
                v2d.x(),
                &self.parent.d1u1,
                v2d.y(),
                &self.parent.d1v1,
            );

            d.set_value(1, 1, self.parent.nplan.dot(&temp));
            d.set_value(2, 1, 0.0);
            d.set_value(3, 1, 2.0 * temp1.dot(&temp));
            d.set_value(4, 1, temp.dot(&temp3) - temp.dot(&temp1));

            d.set_value(1, 3, 0.0);
            d.set_value(1, 4, 0.0);
            d.set_value(2, 3, self.parent.nplan.dot(&self.parent.d1u2));
            d.set_value(2, 4, self.parent.nplan.dot(&self.parent.d1v2));
            d.set_value(3, 3, 0.0);
            d.set_value(3, 4, 0.0);
            d.set_value(4, 3, temp1.dot(&self.parent.d1u2));
            d.set_value(4, 4, temp1.dot(&self.parent.d1v2));
        } else {
            let temp = GpVec::zero().set_linear_form_2(
                v2d.x(),
                &self.parent.d1u2,
                v2d.y(),
                &self.parent.d1v2,
            );

            d.set_value(1, 1, 0.0);
            d.set_value(2, 1, self.parent.nplan.dot(&temp));
            d.set_value(3, 1, 0.0);
            d.set_value(4, 1, temp1.dot(&temp));

            d.set_value(1, 3, self.parent.nplan.dot(&self.parent.d1u1));
            d.set_value(1, 4, self.parent.nplan.dot(&self.parent.d1v1));
            d.set_value(2, 3, 0.0);
            d.set_value(2, 4, 0.0);
            d.set_value(3, 3, 2.0 * temp1.dot(&self.parent.d1u1));
            d.set_value(3, 4, 2.0 * temp1.dot(&self.parent.d1v1));
            d.set_value(4, 3, self.parent.d1u1.dot(&temp3) - self.parent.d1u1.dot(&temp1));
            d.set_value(4, 4, self.parent.d1v1.dot(&temp3) - self.parent.d1v1.dot(&temp1));
        }

        d.set_value(1, 2, dnplan.dot(&temp1) - self.parent.nplan.dot(&self.parent.d1gui));
        d.set_value(2, 2, dnplan.dot(&temp2) - self.parent.nplan.dot(&self.parent.d1gui));
        d.set_value(3, 2, -2.0 * self.parent.d1gui.dot(&temp1));
        d.set_value(4, 2, -self.parent.d1gui.dot(&temp3));

        true
    }
}

impl MathFunctionSetWithDerivatives for ConstThroatWithPenetrationInv {
    /// `NbVariables()` = `4`, not overridden.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` = `4`, not overridden.
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:44-87`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:92-199`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_derivatives(x, d)
    }

    /// `Values(X, F, D)` (`GenChamfInv.cxx:119-122`), not overridden.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.value(x, f);
        self.derivatives(x, d);
        true
    }
}

impl BlendFuncInv for ConstThroatWithPenetrationInv {
    /// `Set(OnFirst, COnSurf)` (`GenChamfInv.cxx:45-49`), not overridden.
    fn set(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.parent.set_first(on_first, c_on_surf);
    }

    /// `GetTolerance(Tolerance, Tol)` (`GenChamfInv.cxx:52-65`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        self.parent.gen.fill_tolerance(tolerance, tol);
    }

    /// `GetBounds(InfBound, SupBound)` (`GenChamfInv.cxx:69-114`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        self.parent.gen.fill_bounds(inf_bound, sup_bound);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:32-39`): the last value is only compared
    /// with `Tol` (not `Tol^2`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 4);
        self.compute_value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= tol * tol
            && valsol.value(4).abs() <= tol
    }
}
