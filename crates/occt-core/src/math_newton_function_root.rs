//! `math_NewtonFunctionRoot` and `math_TrigonometricEquationFunction`.
//! Sources:
//! * `src/FoundationClasses/TKMath/math/math_NewtonFunctionRoot.cxx`,
//!   `math_NewtonFunctionRoot.hxx`, `math_NewtonFunctionRoot.lxx`.
//! * `src/FoundationClasses/TKMath/math/math_TrigonometricEquationFunction.hxx`.
//!
//! `math_NewtonFunctionRoot` is used by `math_TrigonometricFunctionRoots`
//! (`math_TrigonometricFunctionRoots.cxx:464-468`).

use crate::math_fn::{MathFunction, MathFunctionWithDerivative};

/// `RealLast()` / `RealFirst()` (`Standard_Real.hxx:128-130`).
const REAL_LAST: f64 = f64::MAX;
const REAL_FIRST: f64 = -f64::MAX;

/// `math_NewtonFunctionRoot` (`math_NewtonFunctionRoot.hxx:29-108`).
#[derive(Debug, Clone, Copy)]
pub struct NewtonFunctionRoot {
    done: bool,
    x: f64,
    fx: f64,
    dfx: f64,
    it: i32,
    epsilon_x: f64,
    epsilon_f: f64,
    itermax: i32,
    binf: f64,
    bsup: f64,
}

impl NewtonFunctionRoot {
    /// Ctor `(F, Guess, EpsX, EpsF, NbIterations = 100)` (`cxx:60-76`).
    pub fn new<F: MathFunctionWithDerivative>(
        f: &mut F,
        guess: f64,
        eps_x: f64,
        eps_f: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut r = Self {
            epsilon_x: eps_x,
            epsilon_f: eps_f,
            itermax: nb_iterations,
            binf: REAL_FIRST,
            bsup: REAL_LAST,
            done: false,
            x: REAL_LAST,
            dfx: 0.0,
            fx: REAL_LAST,
            it: 0,
        };
        r.perform(f, guess);
        r
    }

    /// Ctor `(F, Guess, EpsX, EpsF, A, B, NbIterations = 100)` (`cxx:19-37`).
    #[allow(clippy::too_many_arguments)]
    pub fn new_bounded<F: MathFunctionWithDerivative>(
        f: &mut F,
        guess: f64,
        eps_x: f64,
        eps_f: f64,
        a: f64,
        b: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut r = Self {
            epsilon_x: eps_x,
            epsilon_f: eps_f,
            binf: a,
            bsup: b,
            itermax: nb_iterations,
            done: false,
            x: REAL_LAST,
            dfx: 0.0,
            fx: REAL_LAST,
            it: 0,
        };
        r.perform(f, guess);
        r
    }

    /// Protected ctor `(A, B, EpsX, EpsF, NbIterations = 100)`
    /// (`cxx:40-57`): initialises the fields without performing.
    pub fn new_unperformed(a: f64, b: f64, eps_x: f64, eps_f: f64, nb_iterations: i32) -> Self {
        Self {
            binf: a,
            bsup: b,
            epsilon_x: eps_x,
            epsilon_f: eps_f,
            itermax: nb_iterations,
            done: false,
            x: REAL_LAST,
            dfx: 0.0,
            fx: REAL_LAST,
            it: 0,
        }
    }

    /// `Perform` (`cxx:78-155`).
    pub fn perform<F: MathFunctionWithDerivative>(&mut self, f: &mut F, guess: f64) {
        let mut dx;

        // lbr le 12 Nov 97: the best estimate is not saved and a worse
        // solution than Guess is returned.
        let mut best_x = self.x;
        let mut best_fx = REAL_LAST;

        let (aa, bb) = if self.binf < self.bsup {
            (self.binf, self.bsup)
        } else {
            (self.bsup, self.binf)
        };

        dx = REAL_LAST;
        self.fx = REAL_LAST;
        self.x = guess;
        self.it = 1;

        while self.it <= self.itermax && (dx.abs() > self.epsilon_x || self.fx.abs() > self.epsilon_f) {
            let mut fx = 0.0;
            let mut dfx = 0.0;
            let ok = f.values(self.x, &mut fx, &mut dfx);
            self.fx = fx;
            self.dfx = dfx;

            let abs_fx = self.fx.abs();
            if abs_fx < best_fx {
                best_fx = abs_fx;
                best_x = self.x;
            }

            if ok {
                if self.dfx == 0.0 {
                    self.done = false;
                    self.it = self.itermax + 1;
                } else {
                    dx = self.fx / self.dfx;
                    self.x -= dx;
                    // Limit the variations of X.
                    if self.x <= aa {
                        self.x = aa;
                    }
                    if self.x >= bb {
                        self.x = bb;
                    }
                    self.it += 1;
                }
            } else {
                self.done = false;
                self.it = self.itermax + 1;
            }
        }

        self.x = best_x;

        self.done = self.it <= self.itermax;
    }

    /// `IsDone` (`lxx:17-20`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Root` (`lxx:28-32`); `StdFail_NotDone` when not done.
    pub fn root(&self) -> f64 {
        assert!(self.done, "StdFail_NotDone in math_NewtonFunctionRoot::Root");
        self.x
    }

    /// `Derivative` (`lxx:34-38`); `StdFail_NotDone` when not done.
    pub fn derivative(&self) -> f64 {
        assert!(
            self.done,
            "StdFail_NotDone in math_NewtonFunctionRoot::Derivative"
        );
        self.dfx
    }

    /// `Value` (`lxx:40-44`); `StdFail_NotDone` when not done.
    pub fn value(&self) -> f64 {
        assert!(self.done, "StdFail_NotDone in math_NewtonFunctionRoot::Value");
        self.fx
    }

    /// `NbIterations` (`lxx:48-52`); `StdFail_NotDone` when not done.
    pub fn nb_iterations(&self) -> i32 {
        assert!(
            self.done,
            "StdFail_NotDone in math_NewtonFunctionRoot::NbIterations"
        );
        self.it
    }
}

/// `math_TrigonometricEquationFunction`
/// (`math_TrigonometricEquationFunction.hxx:27-79`):
/// `A*cos(x)^2 + 2*B*cos(x)*sin(x) + C*cos(x) + D*sin(x) + E = 0`.
#[derive(Debug, Clone, Copy)]
pub struct TrigonometricEquationFunction {
    aa: f64,
    bb: f64,
    cc: f64,
    dd: f64,
    ee: f64,
}

impl TrigonometricEquationFunction {
    /// `math_TrigonometricEquationFunction(A, B, C, D, E)` (`hxx:33-45`).
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        Self {
            aa: a,
            bb: b,
            cc: c,
            dd: d,
            ee: e,
        }
    }
}

impl MathFunction for TrigonometricEquationFunction {
    /// `Value` (`hxx:47-54`).
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        let cn = x.cos();
        let sn = x.sin();
        *f = cn * (self.aa * cn + (self.bb + self.bb) * sn + self.cc) + self.dd * sn + self.ee;
        true
    }
}

impl MathFunctionWithDerivative for TrigonometricEquationFunction {
    /// `Derivative` (`hxx:56-63`).
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let cn = x.cos();
        let sn = x.sin();
        let mut dd = -self.aa * cn * sn + self.bb * (cn * cn - sn * sn);
        dd += dd;
        dd += -self.cc * sn + self.dd * cn;
        *d = dd;
        true
    }

    /// `Values` (`hxx:65-78`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let cn = x.cos();
        let sn = x.sin();
        let aa_cn = self.aa * cn;
        let bb_sn = self.bb * sn;

        *f = aa_cn * cn + bb_sn * (cn + cn) + self.cc * cn + self.dd * sn + self.ee;
        let mut dd = -aa_cn * sn + self.bb * (cn * cn - sn * sn);
        dd += dd;
        dd += -self.cc * sn + self.dd * cn;
        *d = dd;
        true
    }
}
