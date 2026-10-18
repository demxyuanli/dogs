//! `math_NewtonMinimum` (`math_NewtonMinimum.cxx:34-276`,
//! `math_NewtonMinimum.hxx:34-129`, `math_NewtonMinimum.lxx:19-70`) and the
//! `math_Gauss` linear solver it relies on (`math_Gauss.cxx:23-56`) over
//! `LU_Decompose`/`LU_Solve` (`math_Recipes.cxx:202-296`, `:309-345`).
//!
//! Rust note: the objective is passed to [`NewtonMinimum::perform`] as
//! `&dyn MathMultipleVarFunctionWithHessian`, and `NbVariables` is passed to
//! [`NewtonMinimum::new`] because OCCT reads it from the function in the
//! constructor.

use crate::{
    MathMatrix, MathMultipleVarFunctionWithHessian, MathStatus, MathVector, Jacobi,
};
use occt_core::precision::{PCONFUSION, RESOLUTION};

/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;

/// `math_Gauss` (`math_Gauss.cxx:23-56`): LU decomposition with scaled partial
/// pivoting, exposing `IsDone()`/`Solve(B, X)`.
pub struct Gauss {
    lu: MathMatrix,
    index: Vec<usize>,
    d: f64,
    done: bool,
}

impl Gauss {
    /// `math_Gauss(A, MinPivot)` (`math_Gauss.cxx:23-35`). `MinPivot` is the
    /// `TINY` argument of `LU_Decompose` (`math_Recipes.cxx:202-296`): a row
    /// whose largest element is `<= MinPivot`, or a pivot `<= MinPivot`, makes
    /// the decomposition singular and `IsDone()` false.
    ///
    /// Note: OCCT returns `math_Status_SingularMatrix` from `LU_Decompose` for
    /// the singular case; this port's [`MathStatus`] has no `SingularMatrix`
    /// variant, so the failure is reported as [`MathStatus::FunctionError`].
    /// Only `IsDone()` is consumed, so the distinction is not observable here.
    pub fn new(a: &MathMatrix, min_pivot: f64) -> Self {
        let n = a.row_count();
        let mut lu = a.clone();
        let mut index = vec![0usize; n + 1];
        let mut vv = MathVector::new(1, n);
        let mut d = 1.0;
        let err = lu_decompose(&mut lu, &mut index, &mut d, &mut vv, min_pivot);
        Self {
            lu,
            index,
            d,
            done: err == MathStatus::Ok,
        }
    }

    /// `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `D` (`math_Gauss.cxx:31`), the row-swap sign product.
    pub fn d(&self) -> f64 {
        self.d
    }

    /// `Solve(B, X)` (`math_Gauss.cxx:37-43`): `X = B` then in-place `LU_Solve`.
    /// Returns `None` where OCCT raises `StdFail_NotDone` when `!IsDone()`.
    pub fn solve(&self, b: &MathVector) -> Option<MathVector> {
        if !self.done {
            return None;
        }
        let mut x = b.clone();
        lu_solve(&self.lu, &self.index, &mut x);
        Some(x)
    }
}

/// `LU_Decompose` (`math_Recipes.cxx:202-296`) with the row-scaling vector `vv`
/// (`TINY` is OCCT's `MinPivot`).
///
/// The pivot-choice test is deliberately written as OCCT has it
/// (`:257-262`, "comparison is made so as to have imax updated even if argument
/// is NAN, Inf or IND"): `imax` advances whenever `dum >= big`.
fn lu_decompose(
    a: &mut MathMatrix,
    indx: &mut [usize],
    d: &mut f64,
    vv: &mut MathVector,
    tiny: f64,
) -> MathStatus {
    let n = a.row_count();
    *d = 1.0;

    for i in 1..=n {
        let mut big = 0.0;
        for j in 1..=n {
            let temp = a.value(i, j).abs();
            if temp > big {
                big = temp;
            }
        }
        if big <= tiny {
            return MathStatus::FunctionError;
        }
        vv.set_value(i, 1.0 / big);
    }

    for j in 1..=n {
        for i in 1..j {
            let mut sum = a.value(i, j);
            for k in 1..i {
                sum -= a.value(i, k) * a.value(k, j);
            }
            a.set_value(i, j, sum);
        }

        let mut big = 0.0;
        let mut imax = 0usize;
        for i in j..=n {
            let mut sum = a.value(i, j);
            for k in 1..j {
                sum -= a.value(i, k) * a.value(k, j);
            }
            a.set_value(i, j, sum);

            let dum = vv.value(i) * sum.abs();
            if dum < big {
                continue;
            }
            big = dum;
            imax = i;
        }

        if j != imax {
            for k in 1..=n {
                let dum = a.value(imax, k);
                a.set_value(imax, k, a.value(j, k));
                a.set_value(j, k, dum);
            }
            *d = -*d;
            let t = vv.value(imax);
            vv.set_value(imax, vv.value(j));
            vv.set_value(j, t);
        }

        indx[j] = imax;

        if a.value(j, j).abs() <= tiny {
            return MathStatus::FunctionError;
        }

        if j != n {
            let dum = 1.0 / a.value(j, j);
            for i in (j + 1)..=n {
                let v = a.value(i, j) * dum;
                a.set_value(i, j, v);
            }
        }
    }

    MathStatus::Ok
}

/// `LU_Solve` (`math_Recipes.cxx:309-345`) for 1-based vectors (`nblow = 0`).
fn lu_solve(a: &MathMatrix, indx: &[usize], b: &mut MathVector) {
    let n = a.row_count();
    let mut ii = 0usize;

    for i in 1..=n {
        let ip = indx[i];
        let mut sum = b.value(ip);
        b.set_value(ip, b.value(i));
        if ii != 0 {
            for j in ii..i {
                sum -= a.value(i, j) * b.value(j);
            }
        } else if sum != 0.0 {
            ii = i;
        }
        b.set_value(i, sum);
    }

    for i in (1..=n).rev() {
        let mut sum = b.value(i);
        for j in (i + 1)..=n {
            sum -= a.value(i, j) * b.value(j);
        }
        b.set_value(i, sum / a.value(i, i));
    }
}

/// `math_NewtonMinimum` (`math_NewtonMinimum.hxx:34-129`).
pub struct NewtonMinimum {
    status: MathStatus,
    location: MathVector,
    gradient: MathVector,
    step: MathVector,
    hessian: MathMatrix,
    previous_minimum: f64,
    minimum: f64,
    min_eigen_value: f64,
    x_tol: f64,
    c_tol: f64,
    nbiter: usize,
    no_convex_treatement: bool,
    convex: bool,
    is_bounds_defined: bool,
    left: MathVector,
    right: MathVector,
    done: bool,
    iter_max: usize,
}

impl NewtonMinimum {
    /// `math_NewtonMinimum(theFunction, theTolerance = Precision::Confusion(),
    /// theNbIterations = 40, theConvexity = 1.0e-6, theWithSingularity = true)`
    /// (`math_NewtonMinimum.cxx:34-58`, `math_NewtonMinimum.hxx:41-45`).
    pub fn new(
        nb_variables: usize,
        tolerance: f64,
        nb_iterations: usize,
        convexity: f64,
        with_singularity: bool,
    ) -> Self {
        Self {
            status: MathStatus::NotBracketed,
            location: MathVector::new(1, nb_variables),
            gradient: MathVector::new(1, nb_variables),
            step: MathVector::with_init(1, nb_variables, 10.0 * tolerance),
            hessian: MathMatrix::new(1, nb_variables, 1, nb_variables),
            previous_minimum: 0.0,
            minimum: 0.0,
            min_eigen_value: 0.0,
            x_tol: tolerance,
            c_tol: convexity,
            nbiter: 0,
            no_convex_treatement: with_singularity,
            convex: true,
            is_bounds_defined: false,
            left: MathVector::with_init(1, nb_variables, 0.0),
            right: MathVector::with_init(1, nb_variables, 0.0),
            done: false,
            iter_max: nb_iterations,
        }
    }

    /// OCCT's constructor defaults: `Precision::Confusion()`, 40 iterations,
    /// convexity `1.0e-6`, singularity treatment enabled.
    pub fn with_defaults(nb_variables: usize) -> Self {
        Self::new(nb_variables, occt_core::precision::CONFUSION, 40, 1.0e-6, true)
    }

    /// `SetBoundary(theLeftBorder, theRightBorder)`
    /// (`math_NewtonMinimum.cxx:67-73`).
    pub fn set_boundary(&mut self, left_border: MathVector, right_border: MathVector) {
        self.left = left_border;
        self.right = right_border;
        self.is_bounds_defined = true;
    }

    /// `IsConverged()` (`math_NewtonMinimum.lxx:19-23`).
    pub fn is_converged(&self) -> bool {
        self.step.norm() <= self.x_tol
            || (self.minimum - self.previous_minimum).abs()
                <= self.x_tol * self.previous_minimum.abs()
    }

    /// `IsDone()` (`math_NewtonMinimum.lxx:25-28`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsConvex()` (`math_NewtonMinimum.hxx:63`).
    pub fn is_convex(&self) -> bool {
        self.convex
    }

    /// `Location()` (`math_NewtonMinimum.lxx:30-34`). OCCT raises
    /// `StdFail_NotDone` when `!IsDone()`; the location is returned as-is here.
    pub fn location(&self) -> &MathVector {
        &self.location
    }

    /// `Minimum()` (`math_NewtonMinimum.lxx:42-46`).
    pub fn minimum(&self) -> f64 {
        self.minimum
    }

    /// `NbIterations()` (`math_NewtonMinimum.lxx:60-64`).
    pub fn nb_iterations(&self) -> usize {
        self.nbiter
    }

    /// `GetStatus()` (`math_NewtonMinimum.lxx:66-69`).
    pub fn status(&self) -> MathStatus {
        self.status
    }

    /// `Perform(theFunction, theStartingPoint)`
    /// (`math_NewtonMinimum.cxx:77-263`).
    pub fn perform(
        &mut self,
        f: &dyn MathMultipleVarFunctionWithHessian,
        starting_point: &MathVector,
    ) {
        let n = self.location.upper() - self.location.lower() + 1;
        let mut point1 = starting_point.clone();
        let mut point2 = MathVector::new(1, n);

        let mut ok = true;
        let mut nb_conv = 0usize;
        self.done = true;
        self.status = MathStatus::Ok;
        self.nbiter = 0;

        while ok && nb_conv < 2 {
            self.nbiter += 1;

            // Positionnement (`math_NewtonMinimum.cxx:100-111`).
            let mut v_precedent = 0.0;
            ok = f.values_hessian(
                &point1,
                &mut v_precedent,
                &mut self.gradient,
                &mut self.hessian,
            );
            if !ok {
                self.done = false;
                self.status = MathStatus::FunctionError;
                return;
            }
            if self.nbiter == 1 {
                self.previous_minimum = v_precedent;
                self.minimum = v_precedent;
            }

            // Traitement de la non convexite (`math_NewtonMinimum.cxx:115-142`).
            let calcul_vp = Jacobi::new(&self.hessian);
            if !calcul_vp.is_done {
                self.done = false;
                self.status = MathStatus::FunctionError;
                return;
            }
            self.min_eigen_value = calcul_vp.values.value(calcul_vp.values.min_index());
            if self.min_eigen_value < self.c_tol {
                self.convex = false;
                if self.no_convex_treatement && self.min_eigen_value.abs() > self.c_tol {
                    let delta = self.c_tol + 0.1 * self.min_eigen_value.abs() - self.min_eigen_value;
                    for ii in 1..=n {
                        let v = self.hessian.value(ii, ii) + delta;
                        self.hessian.set_value(ii, ii, v);
                    }
                } else {
                    self.done = false;
                    self.status = MathStatus::FunctionError;
                    return;
                }
            }

            // Schemas de Newton (`math_NewtonMinimum.cxx:144-153`).
            let lu = Gauss::new(&self.hessian, self.c_tol / 100.0);
            if !lu.is_done() {
                self.done = false;
                self.status = MathStatus::DirectionSearchError;
                return;
            }
            match lu.solve(&self.gradient) {
                Some(solved) => self.step = solved,
                None => {
                    self.done = false;
                    self.status = MathStatus::DirectionSearchError;
                    return;
                }
            }

            if self.is_bounds_defined {
                // Project point on bounds or nullify TheStep coords if point
                // lies on boundary (`math_NewtonMinimum.cxx:155-206`).
                for an_idx in 1..=n {
                    point2.set_value(an_idx, point1.value(an_idx) - self.step.value(an_idx));
                }
                let mut mult = REAL_LAST;
                for an_idx in 1..=n {
                    let an_abs_step = self.step.value(an_idx).abs();
                    if an_abs_step < RESOLUTION {
                        continue;
                    }
                    if point2.value(an_idx) < self.left.value(an_idx) {
                        let value_val = (point1.value(an_idx) - self.left.value(an_idx)).abs()
                            / an_abs_step;
                        mult = mult.min(value_val);
                    }
                    if point2.value(an_idx) > self.right.value(an_idx) {
                        let value_val = (point1.value(an_idx) - self.right.value(an_idx)).abs()
                            / an_abs_step;
                        mult = mult.min(value_val);
                    }
                }

                if mult != REAL_LAST {
                    if mult > PCONFUSION {
                        // Project point into param space.
                        self.step = self.step.multiplied_scalar(mult);
                    } else {
                        // Old point on border and new point out of border:
                        // nullify corresponding TheStep indexes.
                        for an_idx in 1..=n {
                            let on_right = (point1.value(an_idx) - self.right.value(an_idx)).abs()
                                < PCONFUSION;
                            let on_left = (point1.value(an_idx) - self.left.value(an_idx)).abs()
                                < PCONFUSION;
                            if (on_right && self.step.value(an_idx) < 0.0)
                                || (on_left && self.step.value(an_idx) > 0.0)
                            {
                                self.step.set_value(an_idx, 0.0);
                            }
                        }
                    }
                }
            }

            let mut has_problem;
            loop {
                for an_idx in 1..=n {
                    point2.set_value(an_idx, point1.value(an_idx) - self.step.value(an_idx));
                }

                has_problem = !f.value(&point2, &mut self.minimum);
                if has_problem {
                    self.step = self.step.divided_scalar(2.0);
                }
                if !has_problem {
                    break;
                }
            }

            if self.is_converged() {
                nb_conv += 1;
            } else {
                nb_conv = 0;
            }

            // Controle et corrections (`math_NewtonMinimum.cxx:230-261`).
            let mut v_itere = self.minimum;
            self.minimum = self.previous_minimum;
            let mut n_reduction = 0usize;
            while v_itere > v_precedent && n_reduction < 10 {
                self.step = self.step.multiplied_scalar(0.4);
                for an_idx in 1..=n {
                    point2.set_value(an_idx, point1.value(an_idx) - self.step.value(an_idx));
                }
                f.value(&point2, &mut v_itere);
                n_reduction += 1;
            }

            if v_itere <= v_precedent {
                std::mem::swap(&mut point1, &mut point2);
                self.previous_minimum = v_precedent;
                self.minimum = v_itere;
                ok = self.nbiter < self.iter_max;
                if !ok && nb_conv < 2 {
                    self.status = MathStatus::TooManyIterations;
                }
            } else {
                ok = false;
                self.status = MathStatus::DirectionSearchError;
            }
        }

        self.location = point1;
    }
}
