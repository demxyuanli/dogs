//! `math_FunctionRoot` (`math_FunctionRoot.cxx:73-136`, `.hxx:32-93`,
//! `.lxx:17-50`) and the `math_MyFunctionSetWithDerivatives` adapter.

use super::*;

/// `math_MyFunctionSetWithDerivatives` (`math_FunctionRoot.cxx:27-71`): wraps a
/// one-variable function so that `math_FunctionSetRoot` can solve it.
struct MyFunctionSetWithDerivatives<'a> {
    f: &'a mut dyn FunctionWithDerivative,
}

impl<'a> MathFunctionSetWithDerivatives for MyFunctionSetWithDerivatives<'a> {
    fn nb_variables(&self) -> usize {
        1
    }
    fn nb_equations(&self) -> usize {
        1
    }
    fn value(&mut self, x: &MathVector, out: &mut MathVector) -> bool {
        let mut v = 0.0;
        let ok = self.f.value(x.value(1), &mut v);
        if ok {
            out.set_value(1, v);
        }
        ok
    }
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let mut v = 0.0;
        let ok = self.f.derivative(x.value(1), &mut v);
        if ok {
            d.set_value(1, 1, v);
        }
        ok
    }
    fn values(&mut self, x: &MathVector, out: &mut MathVector, d: &mut MathMatrix) -> bool {
        let mut fv = 0.0;
        let mut dv = 0.0;
        let ok = self.f.values(x.value(1), &mut fv, &mut dv);
        if ok {
            out.set_value(1, fv);
            d.set_value(1, 1, dv);
        }
        ok
    }
}
// ---------------------------------------------------------------------------
// math_FunctionRoot (math_FunctionRoot.cxx:73-136, .hxx:32-93, .lxx:17-50)
// ---------------------------------------------------------------------------

/// `math_FunctionRoot` (`math_FunctionRoot.hxx:32-93`): root of a function of
/// one variable near an initial guess, computed through
/// `math_FunctionSetRoot` and `math_MyFunctionSetWithDerivatives`.
pub struct MathFunctionRoot {
    done: bool,
    root: f64,
    error: f64,
    derivative: f64,
    nb_iter: i32,
}

impl MathFunctionRoot {
    /// `math_FunctionRoot(F, Guess, Tolerance, NbIterations = 100)`
    /// (`math_FunctionRoot.cxx:73-93`). The unbounded `Perform` uses
    /// `InfBound = RealFirst()`, `SupBound = RealLast()`
    /// (`math_FunctionSetRoot.cxx:721-722`).
    pub fn new(
        f: &mut dyn FunctionWithDerivative,
        guess: f64,
        tolerance: f64,
        nb_iterations: i32,
    ) -> Self {
        Self::new_with_bounds(f, guess, tolerance, REAL_FIRST, REAL_LAST, nb_iterations)
    }

    /// `math_FunctionRoot(F, Guess, Tolerance, A, B, NbIterations = 100)`
    /// (`math_FunctionRoot.cxx:95-119`).
    pub fn new_with_bounds(
        f: &mut dyn FunctionWithDerivative,
        guess: f64,
        tolerance: f64,
        a: f64,
        b: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut v = MathVector::new(1, 1);
        v.set_value(1, guess);
        let mut tol = MathVector::new(1, 1);
        tol.set_value(1, tolerance);
        let mut aa = MathVector::new(1, 1);
        aa.set_value(1, a);
        let mut bb = MathVector::new(1, 1);
        bb.set_value(1, b);

        let (done, nb_iter, root, derivative) = {
            let mut ff = MyFunctionSetWithDerivatives { f: &mut *f };
            let mut sol = MathFunctionSetRoot::new(&ff, &tol, nb_iterations);
            sol.perform_with_bounds(&mut ff, &v, &aa, &bb, false);
            let done = sol.is_done();
            if done {
                (
                    true,
                    sol.nb_iterations(),
                    sol.root().value(1),
                    sol.derivative().value(1, 1),
                )
            } else {
                (false, 0, 0.0, 0.0)
            }
        };

        if done {
            f.state_number();
            let mut error = 0.0;
            f.value(root, &mut error);
            Self {
                done: true,
                root,
                error,
                derivative,
                nb_iter,
            }
        } else {
            Self {
                done: false,
                root: 0.0,
                error: 0.0,
                derivative: 0.0,
                nb_iter: 0,
            }
        }
    }

    /// `IsDone` (`math_FunctionRoot.lxx:17-20`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Root` (`math_FunctionRoot.lxx:28-32`).
    pub fn root(&self) -> f64 {
        assert!(self.done, "math_FunctionRoot: NotDone");
        self.root
    }

    /// `Derivative` (`math_FunctionRoot.lxx:34-38`).
    pub fn derivative(&self) -> f64 {
        assert!(self.done, "math_FunctionRoot: NotDone");
        self.derivative
    }

    /// `Value` (`math_FunctionRoot.lxx:40-44`).
    pub fn value(&self) -> f64 {
        assert!(self.done, "math_FunctionRoot: NotDone");
        self.error
    }

    /// `NbIterations` (`math_FunctionRoot.lxx:46-50`).
    pub fn nb_iterations(&self) -> i32 {
        assert!(self.done, "math_FunctionRoot: NotDone");
        self.nb_iter
    }
}
