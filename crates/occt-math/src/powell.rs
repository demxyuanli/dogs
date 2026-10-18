//! `math_Powell` (`math_Powell.cxx:107-229`, `math_Powell.hxx:31-101`,
//! `math_Powell.lxx:18-22`).
//!
//! Powell's direction-set method: no gradient is needed, each iteration walks
//! the current direction set with a one-dimensional `math_BracketMinimum` +
//! `math_BrentMinimum` line search and then replaces the direction of largest
//! decrease by the iteration displacement.
//!
//! Faithful points worth naming:
//!
//! * `DirFunctionBis` (`cxx:35-78`) is the one-dimensional restriction
//!   `F(P0 + x * Dir)`; its `Value` propagates the objective's `false`
//!   (`cxx:76-77`), which is what makes `MinimizeDirection` fail.
//! * `MinimizeDirection` (`cxx:80-105`) brackets on `[0.0, 1.0]` and refines
//!   with `math_BrentMinimum(1.0e-10)`; on success it scales `Dir` in place and
//!   adds it to `P` (`cxx:98-101`).
//! * `IsSolutionReached` (`lxx:18-22`) is
//!   `2 |PreviousMinimum - TheMinimum| <= XTol (|PreviousMinimum| + |TheMinimum| + EPSZ)`.
//! * The objective's `Value` result is discarded at `cxx:157`/`:205` exactly as
//!   OCCT discards it (a failed evaluation leaves the target untouched).

use crate::{
    BracketMinimum, BrentMinimum, MathMatrix, MathMultipleVarFunction, MathStatus, MathVector,
};

/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;

fn sqr(a: f64) -> f64 {
    a * a
}

/// `F.Value(X, out)` with OCCT's "leave the target untouched on failure"
/// semantics (`math_MultipleVarFunction::Value` returns false and does not
/// write `F`).
fn eval_into(f: &dyn MathMultipleVarFunction, x: &MathVector, out: &mut f64) {
    let mut v = 0.0;
    if f.value(x, &mut v) {
        *out = v;
    }
}

/// `DirFunctionBis` (`math_Powell.cxx:35-78`): the one-dimensional function
/// `F(P0 + x * Dir)` handed to `math_BracketMinimum` / `math_BrentMinimum`.
struct DirFunctionBis<'a> {
    p0: MathVector,
    dir: MathVector,
    f: &'a dyn MathMultipleVarFunction,
}

impl<'a> DirFunctionBis<'a> {
    fn new(nb_variables: usize, f: &'a dyn MathMultipleVarFunction) -> Self {
        Self {
            p0: MathVector::new(1, nb_variables),
            dir: MathVector::new(1, nb_variables),
            f,
        }
    }

    /// `DirFunctionBis::Initialize(p0, dir)` (`cxx:64-68`).
    fn initialize(&mut self, p0: &MathVector, dir: &MathVector) {
        self.p0 = p0.clone();
        self.dir = dir.clone();
    }

    /// `DirFunctionBis::Value(x, fval)` (`cxx:70-78`): `P = P0 + x * Dir`, then
    /// the objective at `P`. `None` is the `false` of `F->Value`.
    fn value(&self, x: f64) -> Option<f64> {
        let p = self.p0.added(&self.dir.multiplied_scalar(x));
        let mut fval = 0.0;
        if self.f.value(&p, &mut fval) {
            Some(fval)
        } else {
            None
        }
    }
}

/// `MinimizeDirection(P, Dir, Result, F)` (`math_Powell.cxx:80-105`).
///
/// On success `Dir` becomes `Dir * Scale` and `P` becomes `P + Dir` — the
/// in-place updates of `cxx:98-100`.
fn minimize_direction(
    p: &mut MathVector,
    dir: &mut MathVector,
    result: &mut f64,
    f_dir: &mut DirFunctionBis<'_>,
) -> bool {
    f_dir.initialize(p, dir);

    // math_BracketMinimum Bracket(F, 0.0, 1.0);
    let mut bracket = BracketMinimum::new(0.0, 1.0);
    let mut objective = |x: f64| f_dir.value(x);
    bracket.perform_fallible(&mut objective);
    if bracket.is_done() {
        if let Some((ax, xx, bx)) = bracket.values() {
            // math_BrentMinimum Sol(1.0e-10);
            let mut sol = BrentMinimum::with_tolerance(1.0e-10);
            sol.perform_fallible(&mut objective, ax, xx, bx);
            if sol.is_done() {
                let scale = sol.location().unwrap_or(0.0);
                *result = sol.minimum().unwrap_or(0.0);
                dir.multiply_scalar(scale);
                p.add(dir);
                return true;
            }
        }
    }
    false
}

/// `math_Powell` (`math_Powell.hxx:31`).
pub struct Powell {
    location: MathVector,
    minimum: f64,
    location_error: f64,
    previous_minimum: f64,
    x_tol: f64,
    eps_z: f64,
    done: bool,
    iter: usize,
    status: MathStatus,
    directions: MathMatrix,
    state: i32,
    iter_max: usize,
}

impl Powell {
    /// `math_Powell(theFunction, theTolerance, theNbIterations = 200,
    /// theZEPS = 1.0e-12)` (`math_Powell.cxx:107-127`, `math_Powell.hxx:38-41`).
    pub fn new(nb_variables: usize, tolerance: f64, nb_iterations: usize, zeps: f64) -> Self {
        Self {
            location: MathVector::new(1, nb_variables),
            minimum: REAL_LAST,
            location_error: REAL_LAST,
            previous_minimum: REAL_LAST,
            x_tol: tolerance,
            eps_z: zeps,
            done: false,
            iter: 0,
            status: MathStatus::NotBracketed,
            directions: MathMatrix::new(1, nb_variables, 1, nb_variables),
            state: 0,
            iter_max: nb_iterations,
        }
    }

    /// The `math_Powell(*myFunc, 1e-10)` form of `math_GlobOptMin.cxx:322`.
    pub fn with_tolerance(nb_variables: usize, tolerance: f64) -> Self {
        Self::new(nb_variables, tolerance, 200, 1.0e-12)
    }

    /// `IsSolutionReached(F)` (`math_Powell.lxx:18-22`).
    pub fn is_solution_reached(&self) -> bool {
        2.0 * (self.previous_minimum - self.minimum).abs()
            <= self.x_tol * (self.previous_minimum.abs() + self.minimum.abs() + self.eps_z)
    }

    /// `IsDone()` (`math_Powell.lxx:24-27`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Location()` (`math_Powell.lxx:29-33`). OCCT raises `StdFail_NotDone`
    /// when `!IsDone()`; the location is returned as-is here.
    pub fn location(&self) -> MathVector {
        self.location.clone()
    }

    /// `Minimum()` (`math_Powell.lxx:35-39`).
    pub fn minimum(&self) -> f64 {
        self.minimum
    }

    /// `NbIterations()` (`math_Powell.lxx:49-53`).
    pub fn nb_iterations(&self) -> usize {
        self.iter
    }

    /// `GetStatus()`.
    pub fn status(&self) -> MathStatus {
        self.status
    }

    /// `State` (`math_Powell.hxx:99`): the function's `GetStateNumber()` at the
    /// terminating iteration (`cxx:185`).
    pub fn state(&self) -> i32 {
        self.state
    }

    /// `LocationError` (`math_Powell.hxx:89`) — written by no arm of the `.cxx`;
    /// exposed because OCCT carries the member.
    pub fn location_error(&self) -> f64 {
        self.location_error
    }

    /// `Perform(F, StartingPoint, StartingDirections)`
    /// (`math_Powell.cxx:136-229`).
    pub fn perform(
        &mut self,
        f: &dyn MathMultipleVarFunction,
        starting_point: &MathVector,
        starting_directions: &MathMatrix,
    ) {
        self.done = false;
        let n = self.location.len();

        let mut pt = MathVector::new(1, n);
        let mut ptt = MathVector::new(1, n);
        let mut xit = MathVector::new(1, n);
        let mut f_dir = DirFunctionBis::new(n, f);

        self.location = starting_point.clone();
        self.directions = starting_directions.clone();
        pt = self.location.clone(); // sauvegarde du point initial

        for iter in 1..=self.iter_max {
            self.iter = iter;

            eval_into(f, &self.location, &mut self.previous_minimum);
            let mut ibig = 0usize;
            let mut del = 0.0;

            for i in 1..=n {
                for j in 1..=n {
                    xit.set_value(j, self.directions.value(j, i));
                }
                let mut fptt = 0.0;
                eval_into(f, &self.location, &mut fptt);
                let is_good = minimize_direction(&mut self.location, &mut xit, &mut self.minimum, &mut f_dir);

                if !is_good {
                    self.done = false;
                    self.status = MathStatus::DirectionSearchError;
                    return;
                }

                if (fptt - self.minimum).abs() > del {
                    del = (fptt - self.minimum).abs();
                    ibig = i;
                }
            }

            if self.is_solution_reached() {
                // Termination criterion (`cxx:182-189`).
                self.state = f.state_number();
                self.done = true;
                self.status = MathStatus::Ok;
                return;
            }

            if iter == self.iter_max {
                self.done = false;
                self.status = MathStatus::TooManyIterations;
                return;
            }

            // ptt = 2.0 * TheLocation - pt; xit = TheLocation - pt; pt = TheLocation;
            ptt = self.location.multiplied_scalar(2.0).subtracted(&pt);
            xit = self.location.subtracted(&pt);
            pt = self.location.clone();

            // Value of the function at the extrapolated point (`cxx:203-205`).
            let mut fptt = 0.0;
            eval_into(f, &ptt, &mut fptt);

            if fptt < self.previous_minimum {
                let t = 2.0
                    * (self.previous_minimum - 2.0 * self.minimum + fptt)
                    * sqr(self.previous_minimum - self.minimum - del)
                    - del * sqr(self.previous_minimum - fptt);
                if t < 0.0 {
                    // Minimisation along the direction (`cxx:213-226`).
                    let is_good = minimize_direction(
                        &mut self.location,
                        &mut xit,
                        &mut self.minimum,
                        &mut f_dir,
                    );
                    if !is_good {
                        self.done = false;
                        self.status = MathStatus::FunctionError;
                        return;
                    }

                    for j in 1..=n {
                        self.directions.set_value(j, ibig, xit.value(j));
                    }
                }
            }
        }
    }
}

impl Default for Powell {
    fn default() -> Self {
        Self::new(2, 1.0e-7, 200, 1.0e-12)
    }
}
