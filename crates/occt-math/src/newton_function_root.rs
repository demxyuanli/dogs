//! `math_NewtonFunctionRoot` — one-dimensional Newton root search on a
//! function that also answers its derivative.
//!
//! Source: `FoundationClasses/TKMath/math/math_NewtonFunctionRoot.cxx:17-171`
//! and `.hxx` (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `math_NewtonFunctionRoot(F, Guess, EpsX, EpsF, A, B, NbIterations)` | `cxx:19-38` |
//! | `math_NewtonFunctionRoot(A, B, EpsX, EpsF, NbIterations)` | `cxx:41-57` |
//! | `math_NewtonFunctionRoot(F, Guess, EpsX, EpsF, NbIterations)` | `cxx:59-74` |
//! | `Perform(F, Guess)` | `cxx:76-153` |
//!
//! Translation notes:
//! - `math_FunctionWithDerivative` is the OCCT abstract base; the port already
//!   has it as [`crate::function_set_root::FunctionWithDerivative`]
//!   (`function_set_root/mod.rs:80-87`), whose `values` is the
//!   `F.Values(X, Fx, DFx)` of `Perform`,
//! - `RealFirst()` / `RealLast()` are `-DBL_MAX` / `DBL_MAX`
//!   (`Standard_Real.hxx:167-183`),
//! - `Dump` (`cxx:155-171`) is `UNPORTED` (no `Standard_OStream` in the port).

use crate::function_set_root::FunctionWithDerivative;

/// `math_NewtonFunctionRoot` (`math_NewtonFunctionRoot.hxx`).
#[derive(Debug, Clone)]
pub struct NewtonFunctionRoot {
    /// `EpsilonX` / `EpsilonF` (`cxx:25-26`, `:47-50`).
    epsilon_x: f64,
    epsilon_f: f64,
    /// `Binf` / `Bsup` (`cxx:28-30`, `:46-48`, `:68-69`).
    binf: f64,
    bsup: f64,
    /// `Itermax` (`cxx:31`).
    itermax: i32,
    /// `Done` (`cxx:32`).
    done: bool,
    /// `X` (`cxx:33`).
    x: f64,
    /// `DFx` / `Fx` (`cxx:34-36`).
    dfx: f64,
    fx: f64,
    /// `It` (`cxx:37`).
    it: i32,
}

impl NewtonFunctionRoot {
    /// A fresh, not-yet-performed instance with the bounds
    /// `math_NewtonFunctionRoot(A, B, EpsX, EpsF, NbIterations)` stores
    /// (`cxx:41-57`). The caller runs [`Self::perform`] itself.
    pub fn with_bounds(a: f64, b: f64, eps_x: f64, eps_f: f64, nb_iterations: i32) -> Self {
        Self {
            epsilon_x: eps_x,
            epsilon_f: eps_f,
            binf: a,
            bsup: b,
            itermax: nb_iterations,
            done: false,
            x: f64::MAX,
            dfx: 0.0,
            fx: f64::MAX,
            it: 0,
        }
    }

    /// `math_NewtonFunctionRoot(F, Guess, EpsX, EpsF, NbIterations)`
    /// (`cxx:59-74`): unbounded, performed on construction.
    pub fn new<F: FunctionWithDerivative>(
        f: &mut F,
        guess: f64,
        eps_x: f64,
        eps_f: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut r = Self::with_bounds(f64::MIN, f64::MAX, eps_x, eps_f, nb_iterations);
        r.perform(f, guess);
        r
    }

    /// `Perform(F, Guess)` (`cxx:76-153`).
    ///
    /// The best iterate seen so far (`BestX`, `cxx:89`) is what the OCCT body
    /// returns, not the last one: `X = BestX; Done = It <= Itermax;`
    /// (`cxx:151-152`).
    pub fn perform<F: FunctionWithDerivative>(&mut self, f: &mut F, guess: f64) {
        // `cxx:89`
        let mut best_x = self.x;
        let mut best_fx = f64::MAX;

        // `cxx:92-101`
        let (aa, bb) = if self.binf < self.bsup {
            (self.binf, self.bsup)
        } else {
            (self.bsup, self.binf)
        };

        let mut dx = f64::MAX;
        self.fx = f64::MAX;
        self.x = guess;
        self.it = 1;
        while self.it <= self.itermax && (dx.abs() > self.epsilon_x || self.fx.abs() > self.epsilon_f)
        {
            // `cxx:109`: `Ok = F.Values(X, Fx, DFx);`.
            let ok = f.values(self.x, &mut self.fx, &mut self.dfx);

            // `cxx:111-120`: keep the point with the smallest residual.
            let abs_fx = self.fx.abs();
            if abs_fx < best_fx {
                best_fx = abs_fx;
                best_x = self.x;
            }

            if ok {
                if self.dfx == 0.0 {
                    // `cxx:123-127`
                    self.done = false;
                    self.it = self.itermax + 1;
                } else {
                    // `cxx:129-142`
                    dx = self.fx / self.dfx;
                    self.x -= dx;
                    if self.x <= aa {
                        self.x = aa;
                    }
                    if self.x >= bb {
                        self.x = bb;
                    }
                    self.it += 1;
                }
            } else {
                // `cxx:144-148`
                self.done = false;
                self.it = self.itermax + 1;
            }
        }

        // `cxx:151-152`
        self.x = best_x;
        self.done = self.it <= self.itermax;
    }

    /// `IsDone()`.
    pub fn is_done(&self) -> bool { self.done }

    /// `Root()` (`= X`).
    pub fn root(&self) -> f64 { self.x }

    /// `Value()` (`= Fx`).
    pub fn value(&self) -> f64 { self.fx }

    /// `Derivative()` (`= DFx`).
    pub fn derivative(&self) -> f64 { self.dfx }

    /// `NbIterations()` (`= It`).
    pub fn iterations(&self) -> i32 { self.it }
}
