//! `math_NewtonFunctionSetRoot` — damped-free Newton iteration on a square
//! function set, with bound clamping and no line search.
//!
//! Source: `FoundationClasses/TKMath/math/math_NewtonFunctionSetRoot.{hxx,lxx}`
//! and `math_NewtonFunctionSetRoot.cxx:29-174` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `math_Vector` / `math_Matrix` are [`MathVector`] / [`MathMatrix`].
//! - `LU_Decompose(Jacobian, Indx, d, Scratch, 1.0e-30)` (`cxx:120-125`) is
//!   [`crate::lu::lu_decompose`], whose zero-pivot threshold is the same
//!   `1.0e-30`; `LU_Solve` (`cxx:127`) is [`crate::lu::lu_solve`].
//! - `IsSolutionReached` (`lxx:16-22`) is not virtual here: OCCT declares it
//!   `virtual` only so subclasses can weaken the test, and no subclass of this
//!   solver is ported.
//! - The tolerance-less constructor (`cxx:52-68`) leaves `TolX` at the
//!   `math_Vector(Lower, Upper)` default, which is **uninitialised memory**
//!   (`math_VectorBase.lxx:21-27` never calls `Init`). The port stores zeros
//!   there so `is_solution_reached` reduces to the `FValues` test alone; the
//!   uninitialised value itself is not representable.
//! - `Dump` (`cxx:156-174`) and the output-parameter accessors (`lxx:37-65`)
//!   are replaced by cloning accessors, as in the rest of this crate.

use crate::function_set_root::MathFunctionSetWithDerivatives;
use crate::lu::{lu_decompose, lu_solve};
use crate::{MathMatrix, MathVector};

/// `math_NewtonFunctionSetRoot` (`math_NewtonFunctionSetRoot.hxx:33-126`).
pub struct MathNewtonFunctionSetRoot {
    /// `TolX` (`hxx:74`).
    tol_x: MathVector,
    /// `TolF` (`hxx:75`).
    tol_f: f64,
    /// `Sol` (`hxx:83`).
    sol: MathVector,
    /// `DeltaX` (`hxx:85`).
    delta_x: MathVector,
    /// `FValues` (`hxx:86`).
    f_values: MathVector,
    /// `Jacobian` (`hxx:87`).
    jacobian: MathMatrix,
    /// `Done` (`hxx:88`).
    done: bool,
    /// `State` (`hxx:89`).
    state: i32,
    /// `Iter` (`hxx:90`).
    iter: i32,
    /// `Itermax` (`hxx:91`).
    iter_max: i32,
}

impl MathNewtonFunctionSetRoot {
    /// `math_NewtonFunctionSetRoot(F, XTolerance, FTolerance, NbIterations)`
    /// (`cxx:29-48`). `SetTolerance(theXTolerance)` (`cxx:76-82`) copies the
    /// vector component-wise, so a shorter `tol_x` leaves the rest untouched.
    pub fn new(
        function: &mut dyn MathFunctionSetWithDerivatives,
        tol_x: &MathVector,
        tol_f: f64,
        nb_iterations: i32,
    ) -> Self {
        let n = function.nb_variables();
        let mut this = Self::allocate(n, tol_f, nb_iterations);
        this.set_tolerance(tol_x);
        this
    }

    /// `math_NewtonFunctionSetRoot(F, FTolerance, NbIterations)`
    /// (`cxx:52-68`): `TolX` keeps its (uninitialised in OCCT) default; the
    /// port's default is zero.
    pub fn with_free_tolerance(
        function: &mut dyn MathFunctionSetWithDerivatives,
        tol_f: f64,
        nb_iterations: i32,
    ) -> Self {
        Self::allocate(function.nb_variables(), tol_f, nb_iterations)
    }

    fn allocate(n: usize, tol_f: f64, nb_iterations: i32) -> Self {
        Self {
            tol_x: MathVector::with_init(1, n, 0.0),
            tol_f,
            sol: MathVector::with_init(1, n, 0.0),
            delta_x: MathVector::with_init(1, n, 0.0),
            f_values: MathVector::with_init(1, n, 0.0),
            jacobian: MathMatrix::with_init(1, n, 1, n, 0.0),
            done: false,
            state: 0,
            iter: 0,
            iter_max: nb_iterations,
        }
    }

    /// `SetTolerance` (`cxx:76-82`): `for i = 1..TolX.Length(): TolX(i) =
    /// theXTolerance(i)`.
    pub fn set_tolerance(&mut self, tol_x: &MathVector) {
        for i in 1..=self.tol_x.len() {
            self.tol_x.set_value(i, tol_x.value(i));
        }
    }

    /// `Perform(F, StartingPoint)` (`cxx:87-94`): `[-RealFirst, RealLast]`
    /// bounds.
    pub fn perform(
        &mut self,
        function: &mut dyn MathFunctionSetWithDerivatives,
        starting_point: &MathVector,
    ) {
        let n = function.nb_variables();
        let inf = MathVector::with_init(1, n, -f64::MAX);
        let sup = MathVector::with_init(1, n, f64::MAX);
        self.perform_bounded(function, starting_point, &inf, &sup);
    }

    /// `Perform(F, StartingPoint, InfBound, SupBound)` (`cxx:96-153`).
    pub fn perform_bounded(
        &mut self,
        function: &mut dyn MathFunctionSetWithDerivatives,
        starting_point: &MathVector,
        inf_bound: &MathVector,
        sup_bound: &MathVector,
    ) {
        self.done = false;
        self.sol = starting_point.clone();
        self.iter = 0;

        let mut ok = function.values(&self.sol, &mut self.f_values, &mut self.jacobian);
        if !ok {
            return;
        }
        for iter in 1..=self.iter_max {
            self.iter = iter;
            for k in 1..=self.delta_x.len() {
                self.delta_x.set_value(k, -self.f_values.value(k));
            }
            // `LU_Decompose(Jacobian, Indx, d, Scratch, 1.0e-30)`: a non-zero
            // return means the Jacobian is singular and the solve is skipped.
            let lu = match lu_decompose(&self.jacobian) {
                Ok(lu) => lu,
                Err(_) => return,
            };
            let step = match lu_solve(&lu, &self.delta_x) {
                Ok(step) => step,
                Err(_) => return,
            };
            for i in 1..=self.sol.len() {
                let mut v = self.sol.value(i) + step.value(i);
                // Limitation de Sol dans les bornes [InfBound, SupBound].
                if v <= inf_bound.value(i) {
                    v = inf_bound.value(i);
                }
                if v >= sup_bound.value(i) {
                    v = sup_bound.value(i);
                }
                self.sol.set_value(i, v);
            }

            ok = function.values(&self.sol, &mut self.f_values, &mut self.jacobian);
            if !ok {
                return;
            }
            if self.is_solution_reached() {
                self.state = function.state_number();
                self.done = true;
                return;
            }
        }
    }

    /// `IsSolutionReached` (`lxx:16-22`).
    fn is_solution_reached(&self) -> bool {
        for i in 1..=self.delta_x.len() {
            if self.delta_x.value(i).abs() > self.tol_x.value(i)
                || self.f_values.value(i).abs() > self.tol_f
            {
                return false;
            }
        }
        true
    }

    /// `IsDone` (`lxx:24-27`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Root()` (`lxx:36-41`): `StdFail_NotDone` when the solve failed.
    pub fn root(&self) -> Result<MathVector, String> {
        if !self.done {
            return Err("StdFail_NotDone: math_NewtonFunctionSetRoot::Root".into());
        }
        Ok(self.sol.clone())
    }

    /// `Derivative()` (`lxx:49-54`).
    pub fn derivative(&self) -> Result<MathMatrix, String> {
        if !self.done {
            return Err("StdFail_NotDone: math_NewtonFunctionSetRoot::Derivative".into());
        }
        Ok(self.jacobian.clone())
    }

    /// `FunctionSetErrors()` (`lxx:63-68`).
    pub fn function_set_errors(&self) -> Result<MathVector, String> {
        if !self.done {
            return Err("StdFail_NotDone: math_NewtonFunctionSetRoot::FunctionSetErrors".into());
        }
        Ok(self.delta_x.clone())
    }

    /// `NbIterations()` (`lxx:77-82`).
    pub fn nb_iterations(&self) -> Result<i32, String> {
        if !self.done {
            return Err("StdFail_NotDone: math_NewtonFunctionSetRoot::NbIterations".into());
        }
        Ok(self.iter)
    }

    /// `GetStateNumber` of the last evaluated function set (`cxx:147`).
    pub fn state_number(&self) -> i32 {
        self.state
    }
}
