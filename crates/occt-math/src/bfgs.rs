//! `math_BFGS` — BFGS quasi-Newton minimization, including OCCT's conditional
//! (hyper-parallelepiped) optimization.
//!
//! Source: `math_BFGS.cxx` (`Perform` `:327-443`, `IsSolutionReached` `:449-454`,
//! constructor `:458-476`, `SetBoundary` `:505-510`, `MinimizeDirection`
//! `:203-321`, `ComputeInitScale` `:115-135`, `ComputeMinMaxScale` `:144-196`),
//! `math_BFGS.hxx:37-129` and `math_BFGS.lxx`.
//!
//! `SetBoundary` (`math_BFGS.cxx:505-510`) makes the algorithm refuse to
//! evaluate the objective outside `[myLeft, myRight]`: the bounds enter
//! `MinimizeDirection` through `ComputeMinMaxScale` (`cxx:222-255`), which
//! turns the box into a feasible interval `[aMinLambda, aMaxLambda]` for the
//! 1-D step, and through `math_BracketMinimum::SetLimits` (`cxx:265-268`) and
//! the "bracket failed on the border" fallback (`cxx:292-319`). `Perform`
//! installs `StartingPoint` verbatim (`cxx:344`) — there is no pre-clamp; a
//! start outside the box makes `ComputeMinMaxScale` return `false` (`cxx:181`).
//!
//! The port keeps OCCT's field names for the boundary state (`myIsBoundsDefined`,
//! `myLeft`, `myRight`, `math_BFGS.hxx:122-124`). `minimize` takes `&self` (the
//! existing port API), so the remaining mutable members of `Perform`
//! (`TheLocation`, `TheGradient`, `PreviousMinimum`, `TheMinimum`, `nbiter`) are
//! locals here; the control flow is unchanged.

use crate::{BracketMinimum, BrentMinimum, MathMatrix, MathStatus, MathVector};

/// `RealSmall()` = `DBL_MIN` (`Standard_Real.hxx:132-135`).
const REAL_SMALL: f64 = f64::MIN_POSITIVE;
/// `Precision::PConfusion()` (`Precision.hxx:334`).
const P_CONFUSION: f64 = 1.0e-9;
/// `Precision::Infinite()` (`Precision.hxx:371`).
const INFINITE: f64 = 2.0e100;

/// BFGS optimizer for (optionally bounded) minimization.
///
/// Defaults follow the OCCT constructor (`math_BFGS.hxx:49-52`): tolerance
/// `1e-8`, at most 200 iterations, ZEPS `1e-12`.
pub struct BFGS {
    /// `Itermax` (`math_BFGS.hxx:128`).
    pub max_iter: usize,
    /// `XTol` (`math_BFGS.hxx:119`).
    pub tolerance: f64,
    /// `EPSZ` (`math_BFGS.hxx:120`).
    pub zeps: f64,
    /// `myIsBoundsDefined` (`math_BFGS.hxx:122`).
    bounds_defined: bool,
    /// `myLeft` (`math_BFGS.hxx:123`).
    left: MathVector,
    /// `myRight` (`math_BFGS.hxx:124`).
    right: MathVector,
}

impl BFGS {
    /// `math_BFGS(NbVariables, Tolerance = 1.0e-8, NbIterations = 200, ZEPS = 1.0e-12)`
    /// (`math_BFGS.cxx:458-476`). The dimension is taken from the starting
    /// point in `minimize`, so the vectors start empty.
    pub fn new() -> Self {
        Self {
            max_iter: 200,
            tolerance: 1.0e-8,
            zeps: 1.0e-12,
            bounds_defined: false,
            left: MathVector::new(1, 1),
            right: MathVector::new(1, 1),
        }
    }

    /// `math_BFGS::SetBoundary(theLeftBorder, theRightBorder)`
    /// (`math_BFGS.cxx:505-510`).
    pub fn set_boundary(&mut self, the_left_border: &MathVector, the_right_border: &MathVector) {
        self.left = the_left_border.clone();
        self.right = the_right_border.clone();
        self.bounds_defined = true;
    }

    /// `myIsBoundsDefined` (`math_BFGS.hxx:122`): true once `SetBoundary` ran.
    pub fn is_boundary(&self) -> bool {
        self.bounds_defined
    }

    /// `math_BFGS::Perform(F, StartingPoint)` (`math_BFGS.cxx:327-443`).
    ///
    /// `f` supplies `F.Value` and `grad` supplies `F.Values`' gradient; the port
    /// has no `math_MultipleVarFunctionWithGradient` hierarchy, so the two
    /// closures replace the single `F.Values(pnt, val, grad)` call. The return
    /// value maps the `math_Status` values of the `.cxx`:
    /// `math_DirectionSearchError` (`cxx:381`) and `math_TooManyIterations`
    /// (`cxx:442`) become `Err`; `math_FunctionError` cannot occur because the
    /// port's closures are infallible.
    pub fn minimize<F, G>(&self, f: F, grad: G, x0: &MathVector) -> Result<MathVector, MathStatus>
    where
        F: Fn(&MathVector) -> f64,
        G: Fn(&MathVector) -> MathVector,
    {
        let n = x0.len();

        let mut location = x0.clone();
        let mut gradient = grad(&location);
        // Good = F.Values(TheLocation, PreviousMinimum, TheGradient) — cxx:344-351
        let mut previous_minimum = f(&location);

        let mut hessin = MathMatrix::new(1, n, 1, n);
        hessin.init(0.0);
        let mut xi = MathVector::new(1, n);
        for i in 1..=n {
            hessin.set_value(i, i, 1.0);
            xi.set_value(i, -gradient.value(i));
        }

        // for (nbiter = 1; nbiter <= Itermax; nbiter++) — cxx:358
        for _nbiter in 1..=self.max_iter {
            // TheMinimum = PreviousMinimum — cxx:360
            let mut minimum = previous_minimum;
            let is_good = self.minimize_direction(
                &mut location,
                previous_minimum,
                &gradient,
                &mut xi,
                &mut minimum,
                &f,
            );

            // if (IsSolutionReached(F)) { Done = true; return; } — cxx:371-376
            if 2.0 * (minimum - previous_minimum).abs()
                <= self.tolerance * (minimum.abs() + previous_minimum.abs() + self.zeps)
            {
                return Ok(location);
            }
            if !is_good {
                // cxx:378-383
                return Err(MathStatus::DirectionSearchError);
            }
            previous_minimum = minimum;

            // dg = TheGradient — cxx:386
            let mut dg = gradient.clone();

            // Good = F.Values(TheLocation, TheMinimum, TheGradient) — cxx:388-394.
            // The value written to TheMinimum there is dead: cxx:360 overwrites
            // it before the next read, so only the fresh gradient is kept.
            let _ = f(&location);
            gradient = grad(&location);

            for i in 1..=n {
                dg.set_value(i, gradient.value(i) - dg.value(i));
            }

            // hdg = hessin * dg — cxx:401-408
            let mut hdg = MathVector::new(1, n);
            for i in 1..=n {
                let mut s = 0.0;
                for j in 1..=n {
                    s += hessin.value(i, j) * dg.value(j);
                }
                hdg.set_value(i, s);
            }

            let mut fac = 0.0;
            let mut fae = 0.0;
            for i in 1..=n {
                fac += dg.value(i) * xi.value(i);
                fae += dg.value(i) * hdg.value(i);
            }
            fac = 1.0 / fac;
            let fad = 1.0 / fae;

            for i in 1..=n {
                dg.set_value(i, fac * xi.value(i) - fad * hdg.value(i));
            }

            for i in 1..=n {
                for j in 1..=n {
                    let a_new = hessin.value(i, j) + fac * xi.value(i) * xi.value(j)
                        - fad * hdg.value(i) * hdg.value(j)
                        + fae * dg.value(i) * dg.value(j);
                    hessin.set_value(i, j, a_new);
                }
            }

            for i in 1..=n {
                let mut s = 0.0;
                for j in 1..=n {
                    s += hessin.value(i, j) * gradient.value(j);
                }
                xi.set_value(i, -s);
            }
        }
        // Done = false; TheStatus = math_TooManyIterations; — cxx:441-442
        Err(MathStatus::TooManyIterations)
    }

    /// `MinimizeDirection` (`math_BFGS.cxx:203-321`): solves the 1-D problem
    /// along `dir`, then advances `p` to the bracket minimum.
    fn minimize_direction<F>(
        &self,
        p: &mut MathVector,
        f0: f64,
        gr: &MathVector,
        dir: &mut MathVector,
        result: &mut f64,
        f: &F,
    ) -> bool
    where
        F: Fn(&MathVector) -> f64,
    {
        let n = p.len();

        // cxx:213-217
        let mut lambda = 0.0;
        if !compute_init_scale(f0, dir, gr, &mut lambda) {
            return false;
        }

        // by default the scaling range is unlimited — cxx:219-221
        let mut a_min_lambda = -INFINITE;
        let mut a_max_lambda = INFINITE;
        if self.bounds_defined {
            // limit the scaling range taking into account the bounds — cxx:222-228
            if !compute_min_max_scale(
                p,
                dir,
                &self.left,
                &self.right,
                &mut a_min_lambda,
                &mut a_max_lambda,
            ) {
                return false;
            }

            if a_min_lambda > -P_CONFUSION && a_max_lambda < P_CONFUSION {
                // Point is on the border and the direction shows outside.
                // Make direction to go along the border. — cxx:230-241
                for an_idx in 1..=n {
                    let on_right = (p.value(an_idx) - self.right.value(an_idx)).abs() < P_CONFUSION
                        && dir.value(an_idx) > 0.0;
                    let on_left = (p.value(an_idx) - self.left.value(an_idx)).abs() < P_CONFUSION
                        && dir.value(an_idx) < 0.0;
                    if on_right || on_left {
                        dir.set_value(an_idx, 0.0);
                    }
                }

                // re-compute scale values with new direction — cxx:243-251
                if !compute_init_scale(f0, dir, gr, &mut lambda) {
                    return false;
                }
                if !compute_min_max_scale(
                    p,
                    dir,
                    &self.left,
                    &self.right,
                    &mut a_min_lambda,
                    &mut a_max_lambda,
                ) {
                    return false;
                }
            }
            lambda = lambda.min(a_max_lambda);
            lambda = lambda.max(a_min_lambda);
        }

        // F.Initialize(P, Dir); F.Value(lambda, F1) — cxx:257-262
        let base = p.clone();
        let direction = dir.clone();
        let eval = |t: f64| -> f64 {
            let mut q = direction.multiplied_scalar(t);
            q.add(&base);
            f(&q)
        };
        let f1 = eval(lambda);

        // math_BracketMinimum Bracket(0.0, lambda); — cxx:264
        let mut bracket = BracketMinimum::new(0.0, lambda);
        if self.bounds_defined {
            // Bracket.SetLimits(aMinLambda, aMaxLambda); — cxx:265-268
            bracket.set_limits(a_min_lambda, a_max_lambda);
        }
        bracket.set_fa(f0);
        bracket.set_fb(f1);
        bracket.perform(&eval);
        if bracket.is_done() {
            // find minimum inside the bracket — cxx:272-291
            let (ax, xx, bx) = bracket.values().unwrap_or((0.0, 0.0, 0.0));
            let (_fax, fxx, _fbx) = bracket.function_values().unwrap_or((0.0, 0.0, 0.0));

            let mut sol = BrentMinimum::with_fbx(1.0e-3, fxx, 100, 1.0e-08);
            sol.perform(&eval, ax, xx, bx);
            if sol.is_done() {
                let scale = sol.location().unwrap_or(0.0);
                *result = sol.minimum().unwrap_or(fxx);
                dir.multiply_scalar(scale);
                p.add(dir);
                return true;
            }
        } else if self.bounds_defined {
            // Bracket definition is failure. If the bounds are defined then
            // set current point to intersection with bounds. — cxx:292-319
            let a_f_min = eval(a_min_lambda);
            let a_f_max = eval(a_max_lambda);
            let a_best_lambda;
            if a_f_min < a_f_max {
                a_best_lambda = a_min_lambda;
                *result = a_f_min;
            } else {
                a_best_lambda = a_max_lambda;
                *result = a_f_max;
            }
            dir.multiply_scalar(a_best_lambda);
            p.add(dir);
            return true;
        }
        false
    }
}

impl Default for BFGS {
    fn default() -> Self {
        Self::new()
    }
}

/// `ComputeInitScale` (`math_BFGS.cxx:115-135`): initial value of the scale
/// factor to apply to the direction.
fn compute_init_scale(f0: f64, dir: &MathVector, gr: &MathVector, scale: &mut f64) -> bool {
    let dy1 = gr.dot(dir);
    if dy1.abs() < REAL_SMALL {
        return false;
    }

    let a_hnr1 = dir.norm2();
    let alfa = 0.7 * (-f0) / dy1;
    *scale = 0.015 / a_hnr1.sqrt();
    if *scale > alfa {
        *scale = alfa;
    }
    true
}

/// `ComputeMinMaxScale` (`math_BFGS.cxx:144-196`): for a given point and
/// direction, and bounding box, find min and max scale factors with which the
/// point reaches borders. Returns `false` if the point is out of bounds.
fn compute_min_max_scale(
    point: &MathVector,
    dir: &MathVector,
    left: &MathVector,
    right: &MathVector,
    min_scale: &mut f64,
    max_scale: &mut f64,
) -> bool {
    for an_idx in 1..=left.upper() {
        let a_left = left.value(an_idx) - point.value(an_idx);
        let a_right = right.value(an_idx) - point.value(an_idx);
        if dir.value(an_idx).abs() > REAL_SMALL {
            // Use PConfusion to get off a little from the bounds to prevent
            // possible refuse in Value function.
            let a_l_scale = (a_left + P_CONFUSION) / dir.value(an_idx);
            let a_r_scale = (a_right - P_CONFUSION) / dir.value(an_idx);
            if a_left.abs() < P_CONFUSION {
                // Point is on the left border. — cxx:161-166
                *max_scale = (*max_scale).min(0.0f64.max(a_r_scale));
                *min_scale = (*min_scale).max(0.0f64.min(a_r_scale));
            } else if a_right.abs() < P_CONFUSION {
                // Point is on the right border. — cxx:167-172
                *max_scale = (*max_scale).min(0.0f64.max(a_l_scale));
                *min_scale = (*min_scale).max(0.0f64.min(a_l_scale));
            } else if a_left * a_right < 0.0 {
                // Point is inside allowed range. — cxx:173-178
                *max_scale = (*max_scale).min(a_l_scale.max(a_r_scale));
                *min_scale = (*min_scale).max(a_l_scale.min(a_r_scale));
            } else {
                // point is out of bounds — cxx:179-183
                return false;
            }
        } else {
            // Direction is parallel to the border.
            // Check that the point is not out of bounds — cxx:185-193
            if a_left > P_CONFUSION || a_right < -P_CONFUSION {
                return false;
            }
        }
    }
    true
}
