//! `BlendFunc_ChamfInv` — the inverse of [`super::Chamfer`]: a solution on a
//! restriction (curve on surface) of one of the two supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ChamfInv.hxx:29-79`,
//! `BlendFunc_ChamfInv.cxx:24-230` (OCCT 8.0.0).
//!
//! The vector `X` is `(t, w, U, V)`; the two [`super::Corde`] systems are the
//! ones of `Chamfer`, with `t` replacing the parameters of one support.
//!
//! Translation notes:
//! - `Derivatives` (`cxx:143-230`) takes the guide-line frame from `corde1`
//!   when `first`, from `corde2` otherwise, and calls `DerFguide` on both
//!   helpers; the `D1` call at the end is always on the *first* surface when
//!   `first` and on the second otherwise (`cxx:169`, `:179`).
//! - `IsSolution` short-circuits on the first `Corde` (`cxx:92-100`).
//! - `csurf` is an `Option`; a null handle makes `Value` / `Derivatives` /
//!   `IsSolution` answer `false` instead of faulting.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec2d};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendFuncInv;

use super::corde::Corde;
use super::gen_chamf_inv::GenChamfInv;

/// `class BlendFunc_ChamfInv : public BlendFunc_GenChamfInv` (`hxx:29-79`).
pub struct ChamfInv {
    /// The `BlendFunc_GenChamfInv` base sub-object.
    pub(super) gen: GenChamfInv,
    /// `corde1` (`hxx:77`).
    corde1: Corde,
    /// `corde2` (`hxx:78`).
    corde2: Corde,
}

impl ChamfInv {
    /// `BlendFunc_ChamfInv(S1, S2, C)` (`cxx:24-31`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        let corde1 = Corde::new(surf1.clone(), curv.clone());
        let corde2 = Corde::new(surf2.clone(), curv.clone());
        Self {
            gen: GenChamfInv::new(surf1, surf2, curv),
            corde1,
            corde2,
        }
    }

    /// `Set(Dist1, Dist2, Choix)` (`cxx:34-72`): the sign of each distance
    /// follows the "quadrant".
    pub fn set_distances(&mut self, dist1: f64, dist2: f64, choix: i32) {
        self.gen.choix = choix;
        let (dis1, dis2) = match choix {
            1 | 2 => (-dist1, -dist2),
            3 | 4 => (dist1, -dist2),
            5 | 6 => (dist1, dist2),
            7 | 8 => (-dist1, dist2),
            _ => (-dist1, -dist2),
        };
        self.corde1.set_dist(dis1);
        self.corde2.set_dist(dis2);
    }

    /// `Set(OnFirst, COnSurf)` (`GenChamfInv.cxx:45-49`).
    pub fn set_first(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.gen.set_first(on_first, c_on_surf);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:76-102`).
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let Some(csurf) = self.gen.csurf.clone() else {
            return false;
        };
        let (p2d, _v2d) = csurf.d1(sol.value(1));

        let mut sol1 = MathVector::new(1, 2);
        let mut sol2 = MathVector::new(1, 2);
        sol1.set_value(1, p2d.x());
        sol1.set_value(2, p2d.y());
        sol2.set_value(1, sol.value(3));
        sol2.set_value(2, sol.value(4));

        if self.gen.first {
            self.corde1.is_solution(&sol1, tol) && self.corde2.is_solution(&sol2, tol)
        } else {
            self.corde1.is_solution(&sol2, tol) && self.corde2.is_solution(&sol1, tol)
        }
    }

    /// `Value(X, F)` (`cxx:106-138`): always `true` when `csurf` is set.
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let Some(csurf) = self.gen.csurf.clone() else {
            return false;
        };
        let (p2d, _v2d) = csurf.d1(x.value(1));
        self.corde1.set_param(x.value(2));
        self.corde2.set_param(x.value(2));

        let mut x1 = MathVector::new(1, 2);
        let mut x2 = MathVector::new(1, 2);
        let mut f1 = MathVector::new(1, 2);
        let mut f2 = MathVector::new(1, 2);
        x1.set_value(1, p2d.x());
        x1.set_value(2, p2d.y());
        x2.set_value(1, x.value(3));
        x2.set_value(2, x.value(4));

        if self.gen.first {
            self.corde1.value(&x1, &mut f1);
            self.corde2.value(&x2, &mut f2);
        } else {
            self.corde1.value(&x2, &mut f1);
            self.corde2.value(&x1, &mut f2);
        }
        f.set_value(1, f1.value(1));
        f.set_value(2, f1.value(2));
        f.set_value(3, f2.value(1));
        f.set_value(4, f2.value(2));

        true
    }

    /// `Derivatives(X, D)` (`cxx:143-230`): always `true` when `csurf` is set.
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let Some(csurf) = self.gen.csurf.clone() else {
            return false;
        };
        let (p2d, v2d) = csurf.d1(x.value(1));
        self.corde1.set_param(x.value(2));
        self.corde2.set_param(x.value(2));

        let mut x1 = MathVector::new(1, 2);
        let mut x2 = MathVector::new(1, 2);
        x1.set_value(1, p2d.x());
        x1.set_value(2, p2d.y());
        x2.set_value(1, x.value(3));
        x2.set_value(2, x.value(4));

        let mut d1 = MathMatrix::new(1, 2, 1, 2);
        let mut d2 = MathMatrix::new(1, 2, 1, 2);
        let mut df1 = GpVec2d::new(0.0, 0.0);
        let mut df2 = GpVec2d::new(0.0, 0.0);

        let ptgui: GpPnt;
        let nplan: occt_core::gp::GpVec;
        let pts: GpPnt;
        let d1u: occt_core::gp::GpVec;
        let d1v: occt_core::gp::GpVec;

        if self.gen.first {
            ptgui = *self.corde1.point_on_guide();
            nplan = *self.corde1.n_plan();
            self.corde2.derivatives(&x2, &mut d2);
            self.corde1.der_fguide(&x1, &mut df1);
            self.corde2.der_fguide(&x2, &mut df2);
            let r = self.gen.surf1.d1(x1.value(1), x1.value(2));
            pts = r.0;
            d1u = r.1;
            d1v = r.2;
        } else {
            ptgui = *self.corde2.point_on_guide();
            nplan = *self.corde2.n_plan();
            self.corde1.derivatives(&x2, &mut d1);
            self.corde1.der_fguide(&x2, &mut df1);
            self.corde2.der_fguide(&x1, &mut df2);
            let r = self.gen.surf2.d1(x1.value(1), x1.value(2));
            pts = r.0;
            d1u = r.1;
            d1v = r.2;
        }

        let vref = occt_core::gp::GpVec::from_pnts(&ptgui, &pts);

        // Derivatives with respect to T (`cxx:184-198`).
        let temp = occt_core::gp::GpVec::zero().set_linear_form_2(v2d.x(), &d1u, v2d.y(), &d1v);
        if self.gen.first {
            d.set_value(1, 1, nplan.dot(&temp));
            d.set_value(2, 1, 2.0 * vref.dot(&temp));
            d.set_value(3, 1, 0.0);
            d.set_value(4, 1, 0.0);
        } else {
            d.set_value(1, 1, 0.0);
            d.set_value(2, 1, 0.0);
            d.set_value(3, 1, nplan.dot(&temp));
            d.set_value(4, 1, 2.0 * vref.dot(&temp));
        }

        // Derivatives with respect to W (`cxx:200-204`).
        d.set_value(1, 2, df1.x());
        d.set_value(2, 2, df1.y());
        d.set_value(3, 2, df2.x());
        d.set_value(4, 2, df2.y());

        // Derivatives with respect to U and V (`cxx:206-226`).
        if self.gen.first {
            for i in 1..=2 {
                for j in 3..=4 {
                    d.set_value(i, j, 0.0);
                    d.set_value(i + 2, j, d2.value(i, j - 2));
                }
            }
        } else {
            for i in 1..=2 {
                for j in 3..=4 {
                    d.set_value(i, j, d1.value(i, j - 2));
                    d.set_value(i + 2, j, 0.0);
                }
            }
        }

        true
    }
}

impl MathFunctionSetWithDerivatives for ChamfInv {
    /// `NbVariables()` = `4`, not overridden.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` = `4` (`GenChamfInv.cxx:36-39`).
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:106-138`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:143-230`).
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

impl BlendFuncInv for ChamfInv {
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

    /// `IsSolution(Sol, Tol)` (`cxx:76-102`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }
}
