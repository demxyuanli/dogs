//! `MyDirFunction`, `SearchDirection` and `Bounds` of
//! `math_FunctionSetRoot.cxx` (`:70-705`).

use super::linear::*;
use super::*;

// ---------------------------------------------------------------------------
// MyDirFunction (math_FunctionSetRoot.cxx:70-195)
// ---------------------------------------------------------------------------

/// `MyDirFunction` (`math_FunctionSetRoot.cxx:70-195`): a one-dimensional
/// restriction of the function set along a direction, plus the
/// value/derivative pack used by `Perform`.
pub(super) struct DirFunction {
    p0: MathVector,
    dir: MathVector,
    p: MathVector,
    fv: MathVector,
}

impl DirFunction {
    pub(super) fn new(ninc: usize, neq: usize) -> Self {
        Self {
            p0: MathVector::new(1, ninc),
            dir: MathVector::new(1, ninc),
            p: MathVector::new(1, ninc),
            fv: MathVector::new(1, neq),
        }
    }

    /// `MyDirFunction::Initialize` (`math_FunctionSetRoot.cxx:114-118`).
    pub(super) fn initialize(&mut self, p0: &MathVector, dir: &MathVector) {
        self.p0 = p0.clone();
        self.dir = dir.clone();
    }

    /// `MyDirFunction::Value(x, fval)` (`math_FunctionSetRoot.cxx:120-150`).
    /// `None` mirrors OCCT's `false` return.
    pub(super) fn value_scalar<F: MathFunctionSetWithDerivatives>(
        &mut self,
        f: &mut F,
        x: f64,
    ) -> Option<f64> {
        for i in self.p.lower()..=self.p.upper() {
            let p = self.dir.value(i);
            self.p.set_value(i, p * x + self.p0.value(i));
        }
        if f.value(&self.p, &mut self.fv) {
            for i in self.fv.lower()..=self.fv.upper() {
                let a = self.fv.value(i);
                if a <= -1.0e100 || a >= 1.0e100 {
                    return None;
                }
            }
            Some(0.5 * self.fv.norm2())
        } else {
            None
        }
    }

    /// `MyDirFunction::Value(Sol, FF, DF, GH, F2, Gnr1)`
    /// (`math_FunctionSetRoot.cxx:152-195`). `None` mirrors OCCT's `false`.
    pub(super) fn value_full<F: MathFunctionSetWithDerivatives>(
        &mut self,
        f: &mut F,
        sol: &MathVector,
        ff: &mut MathVector,
        df: &mut MathMatrix,
        gh: &mut MathVector,
    ) -> Option<(f64, f64)> {
        if f.values(sol, ff, df) {
            for i in ff.lower()..=ff.upper() {
                let a = ff.value(i);
                if a < 0.0 {
                    if a <= -1.0e100 {
                        return None;
                    }
                } else if a >= 1.0e100 {
                    return None;
                }
            }
            let f2 = 0.5 * ff.norm2();
            // GH.TMultiply(DF, FF)
            gh.multiply_vec_mat(ff, df);
            for i in gh.lower()..=gh.upper() {
                if is_infinite(gh.value(i)) {
                    return None;
                }
            }
            let gnr1 = gh.norm2();
            Some((f2, gnr1))
        } else {
            None
        }
    }
}

// ---------------------------------------------------------------------------
// SearchDirection (math_FunctionSetRoot.cxx:439-620)
// ---------------------------------------------------------------------------

/// `SearchDirection(DF, GH, FF, ChangeDirection, InvLengthMax, Direction, Dy)`
/// (`math_FunctionSetRoot.cxx:439-531`). Returns the updated
/// `ChangeDirection`.
pub(super) fn search_direction(
    df: &MathMatrix,
    gh: &MathVector,
    ff: &MathVector,
    change_direction: bool,
    inv_length_max: &MathVector,
    direction: &mut MathVector,
    dy: &mut f64,
) -> bool {
    let ninc = df.col_count();
    let neq = df.row_count();
    let eps = 1.0e-32;
    let mut change = change_direction;
    if !change {
        if ninc == neq {
            for i in ff.lower()..=ff.upper() {
                direction.set_value(i, -ff.value(i));
            }
            let g = lu_decompose(df, 1.0e-9); // math_Gauss(DF, 1.e-9)
            if g.done {
                lu_solve(&g.lu, &g.indx, direction);
            } else {
                let svd = svd_padded(df);
                if svd.is_done {
                    let neg = ff.opposite();
                    svd_solve_occt(&svd, &neg, direction);
                } else {
                    change = true;
                }
            }
        } else if ninc > neq {
            let svd = svd_padded(df);
            if svd.is_done {
                let neg = ff.opposite();
                svd_solve_occt(&svd, &neg, direction);
            } else {
                change = true;
            }
        } else {
            // Ninc < Neq: GaussLeastSquare.
            let gls = GaussLeastSquare::new(df, 1.0e-20);
            if gls.done {
                let neg = ff.opposite();
                gls.solve(&neg, direction);
            } else {
                change = true;
            }
        }
    }
    // Cap the step length (math_FunctionSetRoot.cxx:502-516).
    let mut ratio = (direction.value(direction.lower()) * inv_length_max.value(direction.lower()))
        .abs();
    for i in (direction.lower() + 1)..=direction.upper() {
        ratio = ratio.max((direction.value(i) * inv_length_max.value(i)).abs());
    }
    if ratio > 1.0 {
        direction.divide_scalar(ratio);
    }

    *dy = direction.dot(gh);
    if *dy >= -eps {
        change = true;
    }
    if change {
        for i in direction.lower()..=direction.upper() {
            direction.set_value(i, -gh.value(i));
        }
        *dy = -gh.norm2();
    }
    change
}

/// `SearchDirection(DF, GH, FF, Constraints, X, ChangeDirection, InvLengthMax,
/// Direction, Dy)` (`math_FunctionSetRoot.cxx:534-620`). The unused `X`
/// argument of OCCT is kept as `_x`.
pub(super) fn search_direction_constrained(
    df: &MathMatrix,
    gh: &MathVector,
    ff: &MathVector,
    constraints: &[i32],
    _x: &MathVector,
    change_direction: bool,
    inv_length_max: &MathVector,
    direction: &mut MathVector,
    dy: &mut f64,
) {
    let ninc = df.col_count();
    let neq = df.row_count();
    let mut cons = 0;
    for i in 1..=ninc {
        if constraints[i] != 0 {
            cons += 1;
        }
    }

    if cons == 0 {
        // math_FunctionSetRoot.cxx:564 — ChangeDirection is by value, so the
        // recursion's internal switch is not propagated back.
        search_direction(df, gh, ff, change_direction, inv_length_max, direction, dy);
    } else if cons == ninc {
        for i in direction.lower()..=direction.upper() {
            direction.set_value(i, 0.0);
        }
        *dy = 0.0;
    } else {
        let mut df2 = MathMatrix::new(1, neq, 1, ninc - cons);
        let mut gh_val = MathVector::new(1, ninc - cons);
        let mut direction_val = MathVector::new(1, ninc - cons);
        let mut inv_length_max_val = MathVector::new(1, ninc);
        let mut k = 1;
        for i in 1..=ninc {
            if constraints[i] == 0 {
                gh_val.set_value(k, gh.value(i));
                inv_length_max_val.set_value(k, inv_length_max.value(i));
                direction_val.set_value(k, direction.value(i));
                for j in 1..=neq {
                    df2.set_value(j, k, df.value(j, i));
                }
                k += 1;
            }
        }
        // math_FunctionSetRoot.cxx:596 — the sub-problem's ChangeDirection is by
        // value; the reconstruction below reads this function's own parameter.
        search_direction(
            &df2,
            &gh_val,
            ff,
            change_direction,
            &inv_length_max_val,
            &mut direction_val,
            dy,
        );
        let mut k = 1;
        for i in 1..=ninc {
            if constraints[i] == 0 {
                if !change_direction {
                    direction.set_value(i, direction_val.value(k));
                } else {
                    direction.set_value(i, -gh.value(i));
                }
                k += 1;
            } else {
                direction.set_value(i, 0.0);
            }
        }
    }
}

/// `Bounds(InfBound, SupBound, Tol, Sol, SolSave, Constraints, Delta, IsNewSol)`
/// (`math_FunctionSetRoot.cxx:623-705`). Returns `(Out, IsNewSol)`.
pub(super) fn bounds(
    inf_bound: &MathVector,
    sup_bound: &MathVector,
    tol: &MathVector,
    sol: &mut MathVector,
    sol_save: &MathVector,
    constraints: &mut [i32],
    delta: &mut MathVector,
) -> (bool, bool) {
    let mut out = false;
    let ninc = sol.len();
    let mut monratio: f64 = 1.0;
    let mut is_new_sol = true;

    for i in 1..=ninc {
        constraints[i] = 0;
        delta.set_value(i, sol.value(i) - sol_save.value(i));
        if inf_bound.value(i) == sup_bound.value(i) {
            constraints[i] = 1;
            out = true;
        } else if sol.value(i) < inf_bound.value(i) {
            constraints[i] = 1;
            out = true;
            if -delta.value(i) > tol.value(i) {
                monratio = monratio.min(
                    (inf_bound.value(i) - sol_save.value(i)) / delta.value(i),
                );
            }
        } else if sol.value(i) > sup_bound.value(i) {
            constraints[i] = 1;
            out = true;
            if delta.value(i) > tol.value(i) {
                monratio = monratio.min(
                    (sup_bound.value(i) - sol_save.value(i)) / delta.value(i),
                );
            }
        }
    }

    if out {
        if monratio == 0.0 {
            is_new_sol = false;
            *sol = sol_save.clone();
            delta.init(0.0);
        } else {
            delta.multiply_scalar(monratio);
            *sol = sol_save.added(delta);
            for i in 1..=ninc {
                if sol.value(i) < inf_bound.value(i) {
                    sol.set_value(i, inf_bound.value(i));
                    delta.set_value(i, sol.value(i) - sol_save.value(i));
                } else if sol.value(i) > sup_bound.value(i) {
                    sol.set_value(i, sup_bound.value(i));
                    delta.set_value(i, sol.value(i) - sol_save.value(i));
                }
            }
        }
    }
    (out, is_new_sol)
}
