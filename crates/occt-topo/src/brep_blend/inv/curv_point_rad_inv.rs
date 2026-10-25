//! `BRepBlend_CurvPointRadInv` — reframing between a point and a curve.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_CurvPointRadInv.hxx:31-80`,
//! `.cxx:17-130` (OCCT 8.0.0).
//!
//! Valid for constant and progressive radius (`hxx:40`): the radius does not
//! enter the equations, which only state that the given point and the point of
//! `curv2` lie in the section plane of `curv1` at `w` — the plane normal is the
//! unit tangent of `curv1` and its offset `theD = -(nplan . ptcur1)` puts the
//! origin on `ptcur1`.
//!
//! The vector `X` is `(w, U)`: `w` on the guide curve 1, `U` on the partner
//! curve 2 (`hxx:45-48`). `NbEquations` is `2` (`cxx:53-56`).
//!
//! Translation notes:
//! - `Set(Choix)` (`cxx:48-51`) only records the quadrant; no method of this
//!   class reads it, exactly as in OCCT.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendCurvPointFuncInv;

use super::pnt_vec;

/// `class BRepBlend_CurvPointRadInv : public Blend_CurvPointFuncInv`
/// (`hxx:31-80`).
pub struct CurvPointRadInv {
    /// `curv1` (`hxx:76`).
    curv1: Arc<dyn Curve>,
    /// `curv2` (`hxx:77`).
    curv2: Arc<dyn Curve>,
    /// `point` (`hxx:78`): default-initialised to the origin by the OCCT
    /// constructor (`cxx:21-27` does not list it).
    point: GpPnt,
    /// `choix` (`hxx:79`).
    choix: i32,
}

impl CurvPointRadInv {
    /// `BRepBlend_CurvPointRadInv(C1, C2)` (`cxx:21-27`): `choix = 0`.
    pub fn new(curv1: Arc<dyn Curve>, curv2: Arc<dyn Curve>) -> Self {
        Self {
            curv1,
            curv2,
            point: GpPnt::zero(),
            choix: 0,
        }
    }

    /// `Set(Choix)` (`cxx:48-51`).
    pub fn set_choix(&mut self, choix: i32) {
        self.choix = choix;
    }
}

impl MathFunctionSetWithDerivatives for CurvPointRadInv {
    /// `Blend_CurvPointFuncInv::NbVariables()` = `2`
    /// (`Blend_CurvPointFuncInv.cxx:21-24`).
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`cxx:53-56`): `2`.
    fn nb_equations(&self) -> usize {
        2
    }

    /// `Value(X, F)` (`cxx:58-73`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let (ptcur1, d1cur1) = self.curv1.d1(x.value(1));
        let nplan = d1cur1.normalized();
        // `theD = -(nplan.Dot(ptcur1.XYZ()))` (`cxx:63`).
        let the_d = -nplan.dot(&pnt_vec(&ptcur1));
        let (ptcur2, _d1cur2) = self.curv2.d1(x.value(2));
        f.set_value(1, nplan.dot(&pnt_vec(&self.point)) + the_d);
        f.set_value(2, nplan.dot(&pnt_vec(&ptcur2)) + the_d);
        true
    }

    /// `Derivatives(X, D)` (`cxx:77-99`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (ptcur1, d1cur1, d2cur1) = self.curv1.d2(x.value(1));
        let normd1cur1 = d1cur1.magnitude();
        let unsurnormd1cur1 = 1.0 / normd1cur1;
        let nplan = d1cur1.multiplied_scalar(unsurnormd1cur1);
        let mut dnplan =
            GpVec::zero().set_linear_form_2(-nplan.dot(&d2cur1), &nplan, 1.0, &d2cur1);
        dnplan = dnplan.multiplied_scalar(unsurnormd1cur1);
        let dthe_d = -nplan.dot(&d1cur1) - dnplan.dot(&pnt_vec(&ptcur1));
        d.set_value(1, 1, dnplan.dot(&pnt_vec(&self.point)) + dthe_d);
        d.set_value(1, 2, 0.0);
        let (ptcur2, d1cur2) = self.curv2.d1(x.value(2));
        d.set_value(2, 1, dnplan.dot(&pnt_vec(&ptcur2)) + dthe_d);
        d.set_value(2, 2, nplan.dot(&d1cur2));
        true
    }

    /// `Values(X, F, D)` (`cxx:103-110`): `Value` then `Derivatives`.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.value(x, f);
        self.derivatives(x, d);
        true
    }
}

impl BlendCurvPointFuncInv for CurvPointRadInv {
    /// `Set(P)` (`cxx:114-117`).
    fn set(&mut self, p: GpPnt) {
        self.point = p;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:119-124`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, self.curv1.resolution(tol));
        tolerance.set_value(2, self.curv2.resolution(tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:126-136`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        inf_bound.set_value(1, self.curv1.first_parameter());
        sup_bound.set_value(1, self.curv1.last_parameter());
        inf_bound.set_value(2, self.curv2.first_parameter());
        sup_bound.set_value(2, self.curv2.last_parameter());
    }

    /// `IsSolution(Sol, Tol)` (`cxx:138-146`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 2);
        self.value(sol, &mut valsol);
        valsol.value(1).abs() <= tol && valsol.value(2).abs() <= tol
    }
}
