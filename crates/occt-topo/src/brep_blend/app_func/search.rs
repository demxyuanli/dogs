//! `BRepBlend_AppFuncRoot::SearchPoint` / `SearchLocation`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_AppFuncRoot.cxx:284-432`
//! (OCCT 8.0.0).
//!
//! `SearchPoint` picks a starting solution from the already computed points of
//! the line, refines it with `math_FunctionSetRoot`, and stores the point back
//! into the line when the refinement took more than three iterations.
//!
//! Translation notes:
//! - `myLine->Point(i)` (`cxx:312`, `:329`, `:333`, `:337`, ...) is
//!   `Line::point(i)` returning `Option`; `line_point` reproduces OCCT's
//!   `Standard_OutOfRange` with a panic, since every index here is bounded by
//!   `NbPoints()`.
//! - `mySolver->Perform(Func, XInit, X1, X2)` (`cxx:351`) is
//!   `perform_with_bounds(..., false)`: OCCT's
//!   `math_FunctionSetRoot::Perform(F, X, Inf, Sup)` leaves `StopOnDivergent`
//!   at its default `false`.
//! - `Func` reaches the solver as `FuncSetRef`, because
//!   `MathFunctionSetRoot::perform_with_bounds` is generic over a sized `F`
//!   and `myFunc` is a `dyn BlendAppFunction` in the port.
//! - `#ifdef BREPBLEND_DEB` (`cxx:355-358`) and `#ifdef OCCT_DEBUG`
//!   (`cxx:367-371`, which also calls `mySolver->Dump(std::cout)`) are not
//!   translated.

use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendAppFunction;
use crate::brep_blend::Line;

use super::AppFuncRoot;

/// `myLine->Point(Index)` (`BRepBlend_Line.lxx:57-60`): OCCT raises
/// `Standard_OutOfRange` outside `1..=NbPoints()`.
pub fn line_point(line: &Line, index: usize) -> crate::blend::BlendPoint {
    line.point(index)
        .cloned()
        .expect("Standard_OutOfRange: BRepBlend_Line::Point")
}

/// Sized adapter that lets `math_FunctionSetRoot` drive the trait object
/// `myFunc` (`BRepBlend_AppFuncRoot.cxx:349-353`).
pub(super) struct FuncSetRef<'a>(pub &'a mut dyn BlendAppFunction);

impl MathFunctionSetWithDerivatives for FuncSetRef<'_> {
    fn nb_variables(&self) -> usize {
        self.0.nb_variables()
    }

    fn nb_equations(&self) -> usize {
        self.0.nb_equations()
    }

    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.0.value(x, f)
    }

    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.0.derivatives(x, d)
    }

    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.0.values(x, f, d)
    }

    fn state_number(&mut self) -> i32 {
        self.0.state_number()
    }
}

impl AppFuncRoot {
    /// `SearchPoint(Func, Param, Pnt)` (`cxx:300-376`).
    ///
    /// Algorithm (`cxx:290-296`): (1) an approximative solution is taken from
    /// the points already in the line, (2) it is refined by a Newton-type
    /// method, (3) the point is stored, (4) a slow evaluation is inserted into
    /// the line.
    pub(super) fn search_point(&self, param: f64) -> bool {
        let dim = self.dim;
        let nb_points = self.line.borrow().nb_points();
        if nb_points == 0 {
            return false;
        }

        // (1.a) `Param < Point(I1).Parameter()` / `Param > Point(I2).Parameter()`
        // (`cxx:310-318`).
        let t_first = line_point(&self.line.borrow(), 1).parameter();
        if param < t_first {
            return false;
        }
        let t_last = line_point(&self.line.borrow(), nb_points).parameter();
        if param > t_last {
            return false;
        }

        // (1.b) Find the interval (`cxx:320-321`).
        let (trouve, index) = self.search_location(param, 1, nb_points);

        let mut x_init = self.x_init.borrow_mut();
        if trouve {
            // (1.c) If the point is already calculated it is returned.
            let pnt = line_point(&self.line.borrow(), index);
            *self.pnt.borrow_mut() = pnt.clone();
            if self.contact.vec(&mut x_init, &pnt).is_err() {
                return false;
            }
        } else {
            // (1.d) Initialisation by linear interpolation (`cxx:330-347`).
            let pnt1 = line_point(&self.line.borrow(), index);
            *self.pnt.borrow_mut() = pnt1.clone();
            if self.contact.vec(&mut self.x1.borrow_mut(), &pnt1).is_err() {
                return false;
            }
            let t1 = pnt1.parameter();

            let pnt2 = line_point(&self.line.borrow(), index + 1);
            *self.pnt.borrow_mut() = pnt2.clone();
            if self.contact.vec(&mut self.x2.borrow_mut(), &pnt2).is_err() {
                return false;
            }
            let t2 = pnt2.parameter();

            let param_mt1 = (param - t1) / (t2 - t1);
            let t2_mparam = (t2 - param) / (t2 - t1);
            for i in 1..=dim {
                let x2 = self.x2.borrow().value(i);
                let x1 = self.x1.borrow().value(i);
                x_init.set_value(i, x2 * param_mt1 + x1 * t2_mparam);
            }
        }

        // (2) Calculation of the solution (`cxx:349-360`).
        let mut func = self.func.borrow_mut();
        func.set(param);
        func.get_bounds(&mut self.x1.borrow_mut(), &mut self.x2.borrow_mut());
        let sol = {
            let mut solver = self.solver.borrow_mut();
            let mut adapter = FuncSetRef(&mut *func);
            let inf = self.x1.borrow().clone();
            let sup = self.x2.borrow().clone();
            solver.perform_with_bounds(&mut adapter, &x_init, &inf, &sup, false);
            if !solver.is_done() {
                return false;
            }
            solver.root()
        };
        drop(x_init);

        // (3) Storage of the point (`cxx:362-363`).
        let pnt = match self.contact.point(&*func, param, &sol) {
            Ok(p) => p,
            Err(_) => return false,
        };
        *self.pnt.borrow_mut() = pnt.clone();
        *self.sol.borrow_mut() = sol;

        // (4) Insertion of the point if the calculation seems long
        // (`cxx:365-375`).
        let iterations = {
            let solver = self.solver.borrow();
            solver.nb_iterations()
        };
        if !trouve && iterations > 3 {
            self.line.borrow_mut().insert_before(index + 1, pnt);
        }
        true
    }

    /// `SearchLocation(Param, FirstIndex, LastIndex, ParamIndex)`
    /// (`cxx:387-432`): binary search of the parametric interval containing
    /// `Param`; answers `true` when `Param` is already a stored point.
    pub(super) fn search_location(
        &self,
        param: f64,
        first_index: usize,
        last_index: usize,
    ) -> (bool, usize) {
        let line = self.line.borrow();
        let mut ideb = first_index;
        let mut ifin = last_index;

        let valeur = line_point(&line, ideb).parameter();
        if param == valeur {
            return (true, ideb);
        }

        let valeur = line_point(&line, ifin).parameter();
        if param == valeur {
            return (true, ifin);
        }

        while ideb + 1 != ifin {
            let idemi = (ideb + ifin) / 2;
            let valeur = line_point(&line, idemi).parameter();
            if valeur < param {
                ideb = idemi;
            } else if valeur > param {
                ifin = idemi;
            } else {
                return (true, idemi);
            }
        }

        (false, ideb)
    }
}
