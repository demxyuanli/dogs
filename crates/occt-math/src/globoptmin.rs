//! Deterministic global optimization (Evtushenko's non-uniform-mesh algorithm).
//! Source: `math_GlobOptMin.cxx` / `math_GlobOptMin.hxx`.
//!
//! The search box is swept by a deterministic nested grid whose step adapts to
//! the Lipschitz estimate; promising cells are refined by a local (gradient)
//! descent. Finds all local minima inside the box by default; a single-solution
//! mode is available via `find_single_solution`.

use crate::{BFGS, MathVector};

/// `Precision::Infinite()` in OCCT.
const INF: f64 = 2.0e100;

/// Evtushenko global optimizer over a hyper-rectangle.
///
/// Defaults follow the OCCT constructor: Lipschitz constant 9,
/// discretization tolerance 1e-2, functional indifference tolerance 1e-7.
pub struct GlobOptMin {
    /// Parameter-space discretization tolerance (grid fineness).
    pub discretization_tol: f64,
    /// Functional-value indifference tolerance.
    pub same_tol: f64,
    /// Lipschitz constant; auto-estimated during `perform` unless locked.
    pub lip_const: f64,
    /// If true, stop as soon as the functional minimal value is reached.
    pub find_single_solution: bool,
    /// If true, the Lipschitz constant is not re-estimated by `perform`.
    pub lip_const_locked: bool,
    /// Lower bound on the objective considered "good enough"
    /// (`myFunctionalMinimalValue`, `math_GlobOptMin.hxx:238`).
    pub functional_minimal_value: f64,
    /// `myCont` (`math_GlobOptMin.hxx:264`), default 2 (`cxx:69`).
    pub cont: i32,
    /// Stand-in for OCCT's
    /// `dynamic_cast<math_MultipleVarFunctionWithGradient*>(myFunc)`
    /// (`math_GlobOptMin.cxx:294`); see [`GlobOptMin::compute_local_extremum`].
    /// The port's objective is a plain value closure, but any such closure can be
    /// numerically differentiated, so this defaults to `true` — matching the
    /// runtime type OCCT actually sees, since `Extrema_GlobOptFuncCCC2` is a
    /// `math_MultipleVarFunctionWithGradient` (`Extrema_GlobOptFuncCC.hxx`).
    pub has_gradient: bool,
    /// Stand-in for OCCT's
    /// `dynamic_cast<math_MultipleVarFunctionWithHessian*>(myFunc)`
    /// (`math_GlobOptMin.cxx:273`). Defaults to `false`: the port has no
    /// Hessian-input local engine, so the Newton arm is UNPORTED (see
    /// [`GlobOptMin::compute_local_extremum`]).
    pub has_hessian: bool,
    /// Pending `SetLocalParams` sub-box (`myLocalA`/`myLocalB` are not stored in
    /// OCCT — `SetLocalParams` writes `myA`/`myB` directly; the port keeps the
    /// request until `perform` reaches the same point of the sequence).
    local_params: Option<(MathVector, MathVector)>,

    // State.
    done: bool,
    n: usize,
    /// Current best objective value.
    f: f64,
    /// Flat storage of found minima: `n` coordinates per point.
    y: Vec<f64>,
    sol_count: usize,
    /// Initial Lipschitz constant supplied by the user.
    init_lip_const: f64,

    // Algorithm buffers.
    a: MathVector,
    b: MathVector,
    glob_a: MathVector,
    glob_b: MathVector,
    x: MathVector,
    tmp: MathVector,
    v: MathVector,
    max_v: MathVector,
    z: f64,
    e1: f64,
    e2: f64,
    e3: f64,
    last_step: f64,
}

impl Default for GlobOptMin {
    fn default() -> Self {
        Self::new()
    }
}

impl GlobOptMin {
    /// Create with OCCT default parameters.
    pub fn new() -> Self {
        Self {
            discretization_tol: 1.0e-2,
            same_tol: 1.0e-7,
            lip_const: 9.0,
            find_single_solution: false,
            lip_const_locked: false,
            functional_minimal_value: -INF,
            cont: 2,
            has_gradient: true,
            has_hessian: false,
            local_params: None,
            done: false,
            n: 0,
            f: INF,
            y: Vec::new(),
            sol_count: 0,
            init_lip_const: 9.0,
            a: MathVector::new(1, 1),
            b: MathVector::new(1, 1),
            glob_a: MathVector::new(1, 1),
            glob_b: MathVector::new(1, 1),
            x: MathVector::new(1, 1),
            tmp: MathVector::new(1, 1),
            v: MathVector::new(1, 1),
            max_v: MathVector::new(1, 1),
            z: -1.0,
            e1: 0.0,
            e2: 0.0,
            e3: 0.0,
            last_step: 0.0,
        }
    }

    /// Set the two tolerances (discretization, indifference).
    pub fn set_tolerance(&mut self, discretization_tol: f64, same_tol: f64) {
        self.discretization_tol = discretization_tol;
        self.same_tol = same_tol;
    }

    /// Set (and unlock) the Lipschitz constant.
    pub fn set_lip_const(&mut self, c: f64) {
        self.lip_const = c;
        self.init_lip_const = c;
    }

    /// `math_GlobOptMin::SetLipConstState(theFlag)` (`math_GlobOptMin.cxx:116-123`):
    /// when locked, `Perform` does **not** re-estimate the Lipschitz constant
    /// (`cxx:223-227`).
    pub fn set_lip_const_state(&mut self, flag: bool) {
        self.lip_const_locked = flag;
    }

    /// `math_GlobOptMin::GetLipConstState()` (`math_GlobOptMin.hxx:118`).
    pub fn lip_const_state(&self) -> bool {
        self.lip_const_locked
    }

    /// `math_GlobOptMin::SetContinuity(theCont)` (`math_GlobOptMin.hxx:102`).
    ///
    /// OCCT uses it to choose the local-refinement engine inside
    /// `computeLocalExtremum` (`cxx:266-339`: `myCont >= 2` → `math_NewtonMinimum`
    /// on a `math_MultipleVarFunctionWithHessian`, `myCont >= 1` → `math_BFGS`
    /// on a `WithGradient`, else `math_Powell`). The port has no function-class
    /// hierarchy, so the two `dynamic_cast` tests are the explicit
    /// [`GlobOptMin::has_gradient`] / [`GlobOptMin::has_hessian`] flags and the
    /// engine selection lives in [`GlobOptMin::compute_local_extremum`].
    pub fn set_continuity(&mut self, the_cont: i32) {
        self.cont = the_cont;
    }

    /// `math_GlobOptMin::GetContinuity()` (`math_GlobOptMin.hxx:104`).
    pub fn continuity(&self) -> i32 {
        self.cont
    }

    /// `math_GlobOptMin::SetFunctionalMinimalValue(theMinimalValue)`
    /// (`math_GlobOptMin.hxx:107-110`): the value the stop criterion
    /// `CheckFunctionalStopCriteria` (`cxx:600-604`) compares against.
    pub fn set_functional_minimal_value(&mut self, the_minimal_value: f64) {
        self.functional_minimal_value = the_minimal_value;
    }

    /// `math_GlobOptMin::GetFunctionalMinimalValue()` (`math_GlobOptMin.hxx:112`).
    pub fn functional_minimal_value(&self) -> f64 {
        self.functional_minimal_value
    }

    /// `math_GlobOptMin::SetLocalParams(theLocalA, theLocalB)`
    /// (`math_GlobOptMin.cxx:154-171`): restrict the search box to a sub-box
    /// **without** re-running `ComputeInitSol` — OCCT relies on that: the
    /// constructor's `SetGlobalParams` (`cxx:112-148`) already ran
    /// `initCellSize()` + `ComputeInitSol()` on the *global* box, and
    /// `Extrema_GGenExtCC::Perform` then calls `SetLocalParams` per interval
    /// pair before every `Perform` (`Extrema_GGenExtCC.hxx:648-649`).
    ///
    /// The port's [`GlobOptMin::perform`] applies this pending sub-box at the
    /// same point of the sequence (after the global `compute_init_sol`, before
    /// the `e1/e2/e3` setup of `Perform`).
    pub fn set_local_params(&mut self, the_local_a: &MathVector, the_local_b: &MathVector) {
        self.local_params = Some((the_local_a.clone(), the_local_b.clone()));
    }

    /// The pending `SetLocalParams` box, if any.
    pub fn local_params(&self) -> Option<(&MathVector, &MathVector)> {
        self.local_params.as_ref().map(|(a, b)| (a, b))
    }

    /// 2-D convenience wrapper: minimize `f(x, y)` over
    /// `[xmin, xmax] × [ymin, ymax]`. Returns `(x, y, f(x, y))` of the global
    /// minimum. For arbitrary dimensions use [`GlobOptMin::perform`].
    pub fn perform2<F>(&mut self, f: F, xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> Result<(f64, f64, f64), String>
    where
        F: Fn(f64, f64) -> f64,
    {
        let lower = MathVector::from_slice(&[xmin, ymin]);
        let upper = MathVector::from_slice(&[xmax, ymax]);
        let wrap = |x: &MathVector| f(x.value(1), x.value(2));
        self.perform(&wrap, &lower, &upper)?;
        if self.sol_count == 0 {
            return Err("no extrema found".into());
        }
        let p = self.point(0);
        Ok((p.value(1), p.value(2), self.f))
    }

    /// Search for global minima of `f` inside `[lower, upper]`^n.
    ///
    /// This is the port's convenience combination of OCCT's
    /// `SetGlobalParams` + `Perform` (the constructor's form). Callers that need
    /// `math_GlobOptMin`'s incremental drive — one `SetLocalParams` + `Perform`
    /// per parameter sub-box, with the best value and the solution list carried
    /// across the calls, as `Extrema_GGenExtCC::Perform` does
    /// (`Extrema_GGenExtCC.hxx:639-649`) — use [`GlobOptMin::set_global_params`]
    /// followed by [`GlobOptMin::perform_local`].
    pub fn perform<F>(&mut self, f: F, lower: &MathVector, upper: &MathVector) -> Result<(), String>
    where
        F: Fn(&MathVector) -> f64,
    {
        let c = self.lip_const;
        self.set_global_params(&f, lower, upper, c)?;
        self.perform_local(&f)
    }

    /// `math_GlobOptMin::SetGlobalParams` (`cxx:112-148`): install the objective,
    /// the Lipschitz estimate `the_c`, the global box and the working box,
    /// `myMaxV = (b-a)/3`, then `initCellSize()` + `ComputeInitSol()` and
    /// `myDone = false`.
    pub fn set_global_params<F>(
        &mut self,
        f: &F,
        lower: &MathVector,
        upper: &MathVector,
        the_c: f64,
    ) -> Result<(), String>
    where
        F: Fn(&MathVector) -> f64,
    {
        if lower.len() == 0 || lower.len() != upper.len() {
            return Err("dimension mismatch".into());
        }
        if self.discretization_tol <= 0.0 {
            return Err("discretization tolerance must be positive".into());
        }
        let n = lower.len();
        self.n = n;
        self.a = MathVector::new(1, n);
        self.b = MathVector::new(1, n);
        self.glob_a = MathVector::new(1, n);
        self.glob_b = MathVector::new(1, n);
        self.x = MathVector::new(1, n);
        self.tmp = MathVector::new(1, n);
        self.v = MathVector::new(1, n);
        self.max_v = MathVector::new(1, n);
        for i in 1..=n {
            let lo = lower.value(i);
            let hi = upper.value(i);
            self.glob_a.set_value(i, lo);
            self.glob_b.set_value(i, hi);
            self.a.set_value(i, lo);
            self.b.set_value(i, hi);
            self.max_v.set_value(i, (hi - lo) / 3.0);
        }
        self.lip_const = the_c;
        self.init_lip_const = the_c;
        self.z = -1.0;
        self.f = INF;
        self.y.clear();
        self.sol_count = 0;
        self.done = false;

        // `ComputeInitSol()` runs here, i.e. before a pending `SetLocalParams`
        // narrows the box — OCCT's sequence.
        self.compute_init_sol(f)
    }

    /// `math_GlobOptMin::Perform(isFindSingleSolution)` (`cxx:192-262`) driven by
    /// the state left by [`GlobOptMin::set_global_params`] /
    /// [`GlobOptMin::set_local_params`]. **Does not reset** the best value or the
    /// solution list, so successive calls accumulate — which is exactly how
    /// `Extrema_GGenExtCC` uses it.
    pub fn perform_local<F>(&mut self, f: &F) -> Result<(), String>
    where
        F: Fn(&MathVector) -> f64,
    {
        let n = self.n;

        // `math_GlobOptMin::SetLocalParams` (`cxx:154-171`): override the
        // working box (and `myMaxV`), reset `myZ`, keep `myDone = false`.
        if let Some((la, lb)) = self.local_params.clone() {
            for i in 1..=n {
                self.a.set_value(i, la.value(i));
                self.b.set_value(i, lb.value(i));
                self.max_v
                    .set_value(i, (self.b.value(i) - self.a.value(i)) / 3.0);
            }
            self.z = -1.0;
            self.done = false;
        }

        let mut min_length = f64::MAX;
        let mut max_length = f64::MIN;
        for i in 1..=n {
            let l = self.b.value(i) - self.a.value(i);
            min_length = min_length.min(l);
            max_length = max_length.max(l);
            self.v.set_value(i, 0.0);
        }
        if min_length < 1.0e-9 {
            // Degenerate parameter space (`Precision::PConfusion()`, `cxx:214-221`;
            // OCCT returns with `myDone = false`, the port reports it).
            return Err("degenerated parameter space".into());
        }
        if !self.lip_const_locked {
            self.compute_initial_values(f);
        }
        self.e1 = min_length * self.discretization_tol;
        self.e2 = max_length * self.discretization_tol;
        if self.find_single_solution {
            self.e3 = 0.0;
        } else if self.lip_const > 1.0 {
            self.e3 = -max_length * self.discretization_tol / 4.0;
        } else {
            self.e3 = -max_length * self.discretization_tol * self.lip_const / 4.0;
        }

        if self.check_functional_stop_criteria() {
            self.done = true;
            return Ok(());
        }
        self.last_step = 0.0;
        self.compute_global_extremum(f, n);
        self.done = true;
        Ok(())
    }

    /// True once [`GlobOptMin::perform`] has completed.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Number of found global extrema.
    pub fn nb_extrema(&self) -> usize {
        self.sol_count
    }

    /// Best functional value found.
    pub fn minimal_value(&self) -> f64 {
        self.f
    }

    /// X coordinate of the best point (2-D convenience; 0 if absent).
    pub fn x(&self) -> f64 {
        let p = self.point(0);
        if self.n >= 1 {
            p.value(1)
        } else {
            0.0
        }
    }

    /// Y coordinate of the best point (2-D convenience; 0 if absent).
    pub fn y(&self) -> f64 {
        let p = self.point(0);
        if self.n >= 2 {
            p.value(2)
        } else {
            0.0
        }
    }

    /// The `index`-th solution point (0-based); zero vector if out of range.
    pub fn point(&self, index: usize) -> MathVector {
        let mut p = MathVector::new(1, self.n);
        if index < self.sol_count {
            for j in 1..=self.n {
                p.set_value(j, self.y[index * self.n + (j - 1)]);
            }
        }
        p
    }

    /// All solution points found.
    pub fn points(&self) -> Vec<MathVector> {
        (0..self.sol_count).map(|i| self.point(i)).collect()
    }

    // ------------------------------------------------------------------
    // Internal algorithm (faithful port of `math_GlobOptMin.cxx`).
    // ------------------------------------------------------------------

    fn is_inside(&self, pnt: &MathVector) -> bool {
        for i in 1..=self.n {
            if pnt.value(i) < self.glob_a.value(i) || pnt.value(i) > self.glob_b.value(i) {
                return false;
            }
        }
        true
    }

    fn check_functional_stop_criteria(&self) -> bool {
        self.find_single_solution
            && (self.f - self.functional_minimal_value).abs() < self.same_tol * 0.01
    }

    fn compute_initial_values<F: Fn(&MathVector) -> f64>(&mut self, f: &F) {
        const A_MIN_LC: f64 = 0.01;
        const A_MAX_LC: f64 = 1000.0;
        const A_MIN_EPS: f64 = 0.1;
        const A_MAX_EPS: f64 = 100.0;
        const A_PNT_NB: i32 = 13;

        let mut curr_pnt = MathVector::new(1, self.n);
        let mut param_step = MathVector::new(1, self.n);
        let mut prev_val_diag = f(&self.a);
        let mut prev_val_proj = prev_val_diag;
        let mut lip_const: f64 = 0.0;
        let a_step = self.b.subtracted(&self.a).norm() / A_PNT_NB as f64;
        param_step = self.b.subtracted(&self.a).divided_scalar(A_PNT_NB as f64);
        for i in 1..=A_PNT_NB {
            curr_pnt = self.a.added(&param_step.multiplied_scalar(i as f64));
            let val = f(&curr_pnt);
            lip_const = lip_const.max((val - prev_val_diag).abs());
            prev_val_diag = val;

            curr_pnt.set_value(1, self.a.value(1));
            let val = f(&curr_pnt);
            lip_const = lip_const.max((val - prev_val_proj).abs());
            prev_val_proj = val;
        }
        self.lip_const = self.init_lip_const;
        lip_const *= (self.n as f64).sqrt() / a_step;
        if lip_const < self.lip_const * A_MIN_EPS {
            self.lip_const = (lip_const * A_MIN_EPS).max(A_MIN_LC);
        } else if lip_const > self.lip_const * A_MAX_EPS {
            self.lip_const = (self.lip_const * A_MAX_EPS).min(A_MAX_LC);
        }
    }

    fn compute_init_sol<F: Fn(&MathVector) -> f64>(&mut self, f: &F) -> Result<(), String> {
        // Midpoint; protects against local optimizers failing on all inputs.
        let mut pnt = self.glob_a.added(&self.glob_b).multiplied_scalar(0.5);
        let val = f(&pnt);
        if !val.is_finite() {
            return Err("function returned non-finite value".into());
        }
        self.check_add_candidate(&pnt, val);

        // Local optimization from lower corner, midpoint and upper corner.
        for i in 1..=3 {
            let t = (i as f64 - 1.0) / 2.0;
            let start = self.a.added(&self.b.subtracted(&self.a).multiplied_scalar(t));
            if let Some((out_pnt, out_val)) =
                self.compute_local_extremum(f, &start, self.has_gradient, self.has_hessian)
            {
                self.check_add_candidate(&out_pnt, out_val);
            }
        }
        Ok(())
    }

    fn compute_global_extremum<F: Fn(&MathVector) -> f64>(&mut self, f: &F, j: usize) {
        let mut d = f64::MAX;
        let mut a_prev_val = f64::MAX;
        let mut val = f64::MAX;
        let mut a_step_best_value = f64::MAX;
        let mut a_step_best_point = MathVector::new(1, self.n);
        let mut is_inside = false;
        let mut is_reached = false;

        self.x.set_value(j, self.a.value(j) + self.e1);
        while !is_reached {
            if self.x.value(j) > self.b.value(j) {
                self.x.set_value(j, self.b.value(j));
                is_reached = true;
            }
            if self.check_functional_stop_criteria() {
                return;
            }
            if j == 1 {
                is_inside = false;
                a_prev_val = d;
                d = f(&self.x);
                // Evtushenko estimate.
                let r1 = (d + self.z * self.lip_const * self.last_step - self.f) * self.z;
                // Shubert / Piyavsky estimate.
                let r2 = ((d + a_prev_val - self.lip_const * self.last_step) * 0.5 - self.f) * self.z;
                let r = r1.min(r2);
                if r > self.e3 {
                    let save_param = self.x.value(1);
                    // Piyavsky midpoint estimate.
                    let mut a_param = (2.0 * self.x.value(1) - self.v.value(1)) * 0.5
                        + (a_prev_val - d) * 0.5 / self.lip_const;
                    if a_prev_val.is_infinite() {
                        a_param = self.x.value(1) - self.v.value(1) * 0.5;
                    }
                    self.x.set_value(1, a_param);
                    let a_val = f(&self.x);
                    self.x.set_value(1, save_param);

                    if (a_val < d && a_val < a_prev_val)
                        || distance_to_border(&self.x, &self.a, &self.b) < self.e1
                    {
                        if let Some((out_pnt, out_val)) =
                            self.compute_local_extremum(f, &self.x, self.has_gradient, self.has_hessian)
                        {
                            is_inside = true;
                            val = out_val;
                            self.tmp = out_pnt;
                        }
                    }
                }
                a_step_best_value = if is_inside && val < d { val } else { d };
                a_step_best_point = if is_inside && val < d {
                    self.tmp.clone()
                } else {
                    self.x.clone()
                };
                self.check_add_candidate(&a_step_best_point, a_step_best_value);
                if self.check_functional_stop_criteria() {
                    return;
                }
                let step = (self.e2 + (self.f - d).abs() / self.lip_const).min(self.max_v.value(1));
                self.v.set_value(1, step);
                self.last_step = step;
            } else {
                self.v.set_value(j, f64::MAX / 2.0);
                self.compute_global_extremum(f, j - 1);
                // Nullify steps on lower dimensions.
                for i in 1..j {
                    self.v.set_value(i, 0.0);
                }
            }
            if j < self.n {
                let upper_dim_step = self.v.value(j).max(self.e2);
                if self.v.value(j + 1) > upper_dim_step {
                    if upper_dim_step > self.max_v.value(j + 1) {
                        self.v.set_value(j + 1, self.max_v.value(j + 1));
                    } else {
                        self.v.set_value(j + 1, upper_dim_step);
                    }
                }
            }
            if is_reached {
                break;
            }
            self.x.set_value(j, self.x.value(j) + self.v.value(j));
        }
    }

    /// `math_GlobOptMin::computeLocalExtremum` (`math_GlobOptMin.cxx:266-339`):
    /// run a local descent from `pnt` and return the minimizing point and its
    /// value, or `None` when every available engine failed or left the global
    /// box.
    ///
    /// OCCT picks the engine from `myCont` and the runtime type of `myFunc`
    /// (the `dynamic_cast`s at `cxx:273` and `cxx:294`). The port has no
    /// function-class hierarchy, so the two casts are the explicit
    /// `has_hessian` / `has_gradient` arguments:
    ///
    /// * Newton (`cxx:272-291`): `myCont >= 2` + `WithHessian` →
    ///   `math_NewtonMinimum` with `SetBoundary(myGlobA, myGlobB)` — **UNPORTED**,
    ///   see below.
    /// * BFGS (`cxx:293-312`): `myCont >= 1` + `WithGradient` → `math_BFGS` with
    ///   `SetBoundary(myGlobA, myGlobB)` (`cxx:299`). Landed here; the gradient
    ///   is a central difference of the value closure.
    /// * Powell (`cxx:314-336`): the base-class cast always succeeds in OCCT, so
    ///   this is the unconditional fallback — **UNPORTED**, see below.
    ///
    /// Unlike OCCT, a BFGS run that is not `IsDone` does **not** fall through to
    /// the Newton/Powell arms (`cxx:281-291` / `cxx:326-335`): those arms do not
    /// exist in the port, so the function returns `None` exactly as if their
    /// `dynamic_cast` had failed.
    fn compute_local_extremum<F: Fn(&MathVector) -> f64>(
        &self,
        f: &F,
        pnt: &MathVector,
        has_gradient: bool,
        has_hessian: bool,
    ) -> Option<(MathVector, f64)> {
        // Newton method — cxx:272-291.
        //
        // UNPORTED: `math_NewtonMinimum::Perform` (`math_NewtonMinimum.cxx:77-263`)
        // needs `math_MultipleVarFunctionWithHessian::Values(pnt, val, grad, hess)`
        // (`cxx:100`), a Jacobi eigen-decomposition for the convexity treatment
        // (`cxx:115-142`) and a Gauss solve of the Hessian (`cxx:146-153`), plus
        // `SetBoundary` projection (`math_NewtonMinimum.cxx:155-205`). The port's
        // `newton::NewtonMinimum` (`newton.rs:94-169`) is a BFGS-update descent
        // without a Hessian input, so grafting that boundary control flow onto it
        // would invent an algorithm. The arm is skipped, as a failed cast would be.
        let _ = has_hessian;

        // BFGS method used — cxx:293-312.
        if self.cont >= 1 && has_gradient {
            let mut bfgs = BFGS::new();
            bfgs.set_boundary(&self.glob_a, &self.glob_b); // cxx:299
            let grad = |x: &MathVector| numeric_gradient(f, x);
            if let Ok(x) = bfgs.minimize(f, grad, pnt) {
                if self.is_inside(&x) {
                    let val = f(&x);
                    if val.is_finite() {
                        return Some((x, val));
                    }
                }
            }
        }

        // Powell method used — cxx:314-336.
        //
        // UNPORTED: OCCT builds an identity direction matrix and runs
        // `math_Powell(*myFunc, 1e-10); powell.Perform(*myFunc, thePnt, m)`
        // (`cxx:317-324`). The port's `Powell` (`powell.rs:8-39`) neither accepts a
        // direction matrix nor follows `math_Powell::Perform`, and OCCT's Powell
        // has no `SetBoundary`, so it cannot stand in for the missing arm.
        None
    }

    fn is_stored(&self, pnt: &MathVector) -> bool {
        let mut a_tol = MathVector::new(1, self.n);
        for i in 1..=self.n {
            a_tol.set_value(i, (self.b.value(i) - self.a.value(i)) * self.same_tol);
        }
        for i in 0..self.sol_count {
            let mut is_same = true;
            for j in 1..=self.n {
                if (pnt.value(j) - self.y[i * self.n + (j - 1)]).abs() > a_tol.value(j) {
                    is_same = false;
                    break;
                }
            }
            if is_same {
                return true;
            }
        }
        false
    }

    fn check_add_candidate(&mut self, pnt: &MathVector, value: f64) {
        // Points within the indifference tolerance of the current best are
        // collected as separate extrema (unless single-solution mode).
        if (value - self.f).abs() < self.same_tol * 0.01 && !self.find_single_solution {
            if !self.is_stored(pnt) {
                if (value - self.f) * self.z > 0.0 {
                    self.f = value;
                }
                for j in 1..=self.n {
                    self.y.push(pnt.value(j));
                }
                self.sol_count += 1;
            }
        }

        // A new best solution.
        let delta = (value - self.f) * self.z;
        if delta > self.same_tol * 0.01 || (delta > 0.0 && self.find_single_solution) {
            self.f = value;
            self.y.clear();
            for j in 1..=self.n {
                self.y.push(pnt.value(j));
            }
            self.sol_count = 1;
        }
    }
}

fn distance_to_border(x: &MathVector, min: &MathVector, max: &MathVector) -> f64 {
    let mut dist = f64::MAX;
    for idx in 1..=x.len() {
        let d1 = (x.value(idx) - min.value(idx)).abs();
        let d2 = (x.value(idx) - max.value(idx)).abs();
        dist = dist.min(d1.min(d2));
    }
    dist
}

/// Central-difference numeric gradient of `f` at `x`.
fn numeric_gradient<F: Fn(&MathVector) -> f64>(f: &F, x: &MathVector) -> MathVector {
    let n = x.len();
    let h = 1e-6;
    let mut g = MathVector::new(1, n);
    for i in 1..=n {
        let mut xp = x.clone();
        let mut xm = x.clone();
        xp.set_value(i, x.value(i) + h);
        xm.set_value(i, x.value(i) - h);
        g.set_value(i, (f(&xp) - f(&xm)) / (2.0 * h));
    }
    g
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn globoptmin_rosenbrock() {
        let mut gom = GlobOptMin::new();
        let rosenbrock = |x: f64, y: f64| {
            (1.0 - x).powi(2) + 100.0 * (y - x * x).powi(2)
        };
        let (x, y, fval) = gom.perform2(rosenbrock, -2.0, 2.0, -1.0, 3.0).expect("converges");
        assert!(fval < 1e-4, "f={}", fval);
        assert!((x - 1.0).abs() < 1e-2, "x={}", x);
        assert!((y - 1.0).abs() < 1e-2, "y={}", y);
        assert!(gom.nb_extrema() >= 1);
        assert!(gom.is_done());
    }

    #[test]
    fn globoptmin_paraboloid() {
        let mut gom = GlobOptMin::new();
        let f = |x: f64, y: f64| x * x + y * y;
        let (x, y, fval) = gom.perform2(f, -1.0, 1.0, -1.0, 1.0).unwrap();
        assert!(fval < 1e-6, "fval={}", fval);
        assert!((x - 0.0).abs() < 1e-2, "x={}", x);
        assert!((y - 0.0).abs() < 1e-2, "y={}", y);
    }

    #[test]
    fn globoptmin_points_and_queries() {
        let mut gom = GlobOptMin::new();
        let f = |x: &MathVector| (x.value(1) - 1.0).powi(2) + (x.value(2) - 2.0).powi(2);
        gom.perform(&f, &MathVector::from_slice(&[0.0, 0.0]), &MathVector::from_slice(&[2.0, 3.0]))
            .unwrap();
        assert_eq!(gom.minimal_value(), gom.minimal_value()); // f is finite
        assert!(gom.minimal_value() < 1e-6);
        assert!(gom.nb_extrema() >= 1);
        let pts = gom.points();
        assert!(!pts.is_empty());
        assert!((pts[0].value(1) - 1.0).abs() < 1e-2);
        assert!((pts[0].value(2) - 2.0).abs() < 1e-2);
    }

    #[test]
    fn globoptmin_degenerate() {
        let mut gom = GlobOptMin::new();
        let f = |x: f64, _y: f64| x * x;
        assert!(gom.perform2(f, 1.0, 1.0, 0.0, 1.0).is_err());
    }

    #[test]
    fn globoptmin_dimension_mismatch() {
        let mut gom = GlobOptMin::new();
        let f = |x: &MathVector| x.value(1) * x.value(1);
        let lo = MathVector::from_slice(&[0.0, 0.0]);
        let hi = MathVector::from_slice(&[1.0]);
        assert!(gom.perform(&f, &lo, &hi).is_err());
    }

    #[test]
    fn globoptmin_nan() {
        let mut gom = GlobOptMin::new();
        let f = |_: f64, _: f64| f64::NAN;
        assert!(gom.perform2(f, -1.0, 1.0, -1.0, 1.0).is_err());
    }

    #[test]
    fn globoptmin_zero_tolerance() {
        let mut gom = GlobOptMin::new();
        gom.set_tolerance(0.0, 1e-7);
        let f = |x: f64, _y: f64| x * x;
        assert!(gom.perform2(f, -1.0, 1.0, -1.0, 1.0).is_err());
    }
}
