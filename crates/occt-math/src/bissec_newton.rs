//! `math_BissecNewton` (`math_BissecNewton.cxx:21-169`, `.hxx:31-118`,
//! `.lxx:17-58`).
//!
//! A bracketed root finder that switches between bisection and Newton steps.
//! `IsSolutionReached` is `virtual` in OCCT (`.hxx:74-75`); the base body
//! (`.lxx:19-22`) is the only one implemented here, exposed as
//! [`MathBissecNewton::is_solution_reached`] so that a caller can reuse the same
//! predicate without subclassing.

use crate::status::MathStatus;
use crate::function_set_root::FunctionWithDerivative;

/// `math_BissecNewton` (`math_BissecNewton.hxx:31-118`).
#[derive(Debug, Clone)]
pub struct MathBissecNewton {
    /// `TheStatus` (`.hxx:96`).
    pub the_status: MathStatus,
    /// `XTol` (`.hxx:97`).
    pub x_tol: f64,
    /// `x` (`.hxx:98`).
    pub x: f64,
    /// `dx` (`.hxx:99`).
    pub dx: f64,
    /// `f` (`.hxx:100`).
    pub f: f64,
    /// `df` (`.hxx:101`).
    pub df: f64,
    /// `Done` (`.hxx:102`).
    pub done: bool,
}

impl MathBissecNewton {
    /// `math_BissecNewton(theXTolerance)` (`cxx:21-31`).
    pub fn new(the_x_tolerance: f64) -> Self {
        Self {
            the_status: MathStatus::NotBracketed,
            x_tol: the_x_tolerance,
            x: 0.0,
            dx: 0.0,
            f: 0.0,
            df: 0.0,
            done: false,
        }
    }

    /// `IsSolutionReached` base body (`.lxx:19-22`): `|dx| <= XTol`.
    pub fn is_solution_reached(&self) -> bool {
        self.dx.abs() <= self.x_tol
    }

    /// `Perform(F, Bound1, Bound2, NbIterations)` (`cxx:36-151`).
    pub fn perform<F: FunctionWithDerivative + ?Sized>(
        &mut self,
        f: &mut F,
        bound1: f64,
        bound2: f64,
        nb_iterations: i32,
    ) {
        // `RealEpsilon()` (`Standard_Real.hxx`): DBL_EPSILON.
        const A_F_TOL: f64 = f64::EPSILON;

        let mut fl = 0.0;
        let mut fh = 0.0;
        let mut good = f.values(bound1, &mut fl, &mut self.df);
        if !good {
            self.done = false;
            self.the_status = MathStatus::FunctionError;
            return;
        }
        good = f.values(bound2, &mut fh, &mut self.df);
        if !good {
            self.done = false;
            self.the_status = MathStatus::FunctionError;
            return;
        }

        if fl * fh > A_F_TOL * A_F_TOL {
            self.done = false;
            self.the_status = MathStatus::NotBracketed;
            return;
        }

        let (mut xl, mut xh);
        if fl < -A_F_TOL || (fl < A_F_TOL && fh < -A_F_TOL) {
            xl = bound1;
            xh = bound2;
        } else {
            xl = bound2;
            xh = bound1;
        }

        self.x = 0.5 * (bound1 + bound2);
        let mut dxold = (bound2 - bound1).abs();
        self.dx = dxold;
        good = f.values(self.x, &mut self.f, &mut self.df);
        if !good {
            self.done = false;
            self.the_status = MathStatus::FunctionError;
            return;
        }

        let mut j = 1;
        while j <= nb_iterations {
            let bisect = (((self.x - xh) * self.df - self.f) * ((self.x - xl) * self.df - self.f)
                >= 0.0)
                || ((2.0 * self.f).abs() > (dxold * self.df).abs());
            if bisect {
                dxold = self.dx;
                self.dx = 0.5 * (xh - xl);
                self.x = xl + self.dx;
                if self.dx.abs() < self.x_tol {
                    self.the_status = MathStatus::Ok;
                    self.done = true;
                    return;
                }
            } else {
                dxold = self.dx;
                self.dx = self.f / self.df;
                let temp = self.x;
                self.x -= self.dx;
                if temp == self.x {
                    self.the_status = MathStatus::Ok;
                    self.done = true;
                    return;
                }
            }

            if self.is_solution_reached() {
                self.the_status = MathStatus::Ok;
                self.done = true;
                return;
            }

            good = f.values(self.x, &mut self.f, &mut self.df);
            if !good {
                self.done = false;
                self.the_status = MathStatus::FunctionError;
                return;
            }
            if self.f < 0.0 {
                xl = self.x;
            } else if self.f > 0.0 {
                xh = self.x;
            } else {
                self.the_status = MathStatus::Ok;
                self.done = true;
                return;
            }
            j += 1;
        }

        self.the_status = MathStatus::TooManyIterations;
        self.done = false;
    }

    /// `IsDone` (`.lxx:34-37`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Root` (`.lxx:39-43`); panics where OCCT raises `StdFail_NotDone`.
    pub fn root(&self) -> f64 {
        assert!(self.done, "math_BissecNewton: NotDone");
        self.x
    }

    /// `Derivative` (`.lxx:45-49`); panics where OCCT raises `StdFail_NotDone`.
    pub fn derivative(&self) -> f64 {
        assert!(self.done, "math_BissecNewton: NotDone");
        self.df
    }

    /// `Value` (`.lxx:51-55`); panics where OCCT raises `StdFail_NotDone`.
    pub fn value(&self) -> f64 {
        assert!(self.done, "math_BissecNewton: NotDone");
        self.f
    }
}
