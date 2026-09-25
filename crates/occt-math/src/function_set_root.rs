//! `math_FunctionSetRoot` / `math_FunctionRoot` — root of a set of N functions
//! of M variables (Newton direction, gradient fallback, boundary handling).
//!
//! Ported from:
//! * `math_FunctionSetRoot.cxx` / `.hxx` — the local `MyDirFunction`
//!   (`:70-195`), `MinimizeDirection` (`:198-264`, `:267-436`),
//!   `SearchDirection` (`:439-531`, `:534-620`), `Bounds` (`:623-705`),
//!   the constructors (`:709-769`) and `Perform` (`:796-1417`).
//! * `math_FunctionRoot.cxx` / `.hxx` / `.lxx` — `math_MyFunctionSetWithDerivatives`
//!   (`:27-71`) and the two constructor wrappers (`:73-119`).
//! * `math_Recipes.cxx` — `LU_Decompose` (`:202-296`), `LU_Solve` (`:309-344`),
//!   `SVD_Solve` (`:712-748`).
//! * `math_Gauss.cxx:23-56`, `math_GaussLeastSquare.cxx:28-49`,
//!   `math_SVD.cxx:29-68` for the linear solvers used by `SearchDirection`.
//!
//! UNPORTED / not representable:
//! * `math_FunctionSetRoot::IsSolutionReached` is `virtual`
//!   (`math_FunctionSetRoot.hxx:69-79`); only the base implementation is ported
//!   (`MathFunctionSetRoot::is_solution_reached`). Rust has no virtual dispatch
//!   on a non-generic struct, so subclass overrides cannot be expressed.
//! * `Dump` (`math_FunctionSetRoot.cxx:1421-1434`) and the output-parameter
//!   `Root`/`FunctionSetErrors` overloads (`:1438-1452`) are replaced by
//!   cloning accessors.
//! * `math_FunctionSetWithDerivatives::Derivatives` is exposed for interface
//!   fidelity but is never called by `math_FunctionSetRoot.cxx` (only `Values`
//!   and `Value` are).

use crate::{BrentMinimum, MathMatrix, MathVector, SVD};

/// `Precision::Infinite()` (`Precision.hxx:371`).
const PRECISION_INFINITE: f64 = 2.0e100;
/// `RealFirst()` (`Standard_Real.hxx:167`).
const REAL_FIRST: f64 = -f64::MAX;
/// `RealLast()` (`Standard_Real.hxx:179`).
const REAL_LAST: f64 = f64::MAX;
/// `Precision::SquareConfusion()` (`Precision.hxx:169`).
const SQUARE_CONFUSION: f64 = 1.0e-14;

/// `Precision::IsInfinite(R)` (`Precision.hxx:350-353`) = `|R| >= 0.5*Infinite()`.
fn is_infinite(r: f64) -> bool {
    r.abs() >= 0.5 * PRECISION_INFINITE
}

/// `Epsilon(theValue)` (`Standard_Real.hxx:242-246`) — magnitude of the ULP of
/// `theValue` in the direction of the larger magnitude.
fn epsilon_of(v: f64) -> f64 {
    if v == 0.0 {
        f64::from_bits(1)
    } else if v >= 0.0 {
        f64::from_bits(v.to_bits() + 1) - v
    } else {
        v - f64::from_bits(v.to_bits() - 1)
    }
}

/// `math_FunctionSetWithDerivatives` (`math_FunctionSetWithDerivatives.hxx:30-58`)
/// plus `math_FunctionSet::GetStateNumber` (`math_FunctionSet.hxx:59`).
pub trait MathFunctionSetWithDerivatives {
    fn nb_variables(&self) -> usize;
    fn nb_equations(&self) -> usize;
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool;
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool;
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool;
    fn get_state_number(&mut self) -> i32 {
        0
    }
}

/// `math_FunctionWithDerivative` (`math_FunctionWithDerivative.hxx:31-52`) plus
/// `math_Function::GetStateNumber` (`math_Function.hxx:57`).
pub trait MathFunctionWithDerivative {
    fn value(&mut self, x: f64, f: &mut f64) -> bool;
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool;
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool;
    fn get_state_number(&mut self) -> i32 {
        0
    }
}

/// `math_MyFunctionSetWithDerivatives` (`math_FunctionRoot.cxx:27-71`): wraps a
/// one-variable function so that `math_FunctionSetRoot` can solve it.
struct MyFunctionSetWithDerivatives<'a> {
    f: &'a mut dyn MathFunctionWithDerivative,
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
// Linear solvers used by SearchDirection
// ---------------------------------------------------------------------------

/// Result of `LU_Decompose` (`math_Recipes.cxx:202-296`).
struct LuDecomp {
    lu: MathMatrix,
    indx: Vec<usize>,
    done: bool,
}

/// `LU_Decompose(a, indx, d, TINY)` (`math_Recipes.cxx:202-296`), the implicit
/// row-scaling variant used by `math_Gauss` and `math_GaussLeastSquare`.
fn lu_decompose(a: &MathMatrix, tiny: f64) -> LuDecomp {
    let n = a.row_count();
    let mut lu = a.clone();
    let mut indx = vec![0usize; n + 1];
    let mut vv = vec![0.0f64; n + 1];

    for i in 1..=n {
        let mut big = 0.0;
        for j in 1..=n {
            let t = lu.value(i, j).abs();
            if t > big {
                big = t;
            }
        }
        if big <= tiny {
            return LuDecomp {
                lu,
                indx,
                done: false,
            };
        }
        vv[i] = 1.0 / big;
    }

    for j in 1..=n {
        for i in 1..j {
            let mut sum = lu.value(i, j);
            for k in 1..i {
                sum -= lu.value(i, k) * lu.value(k, j);
            }
            lu.set_value(i, j, sum);
        }
        let mut big = 0.0;
        let mut imax = j;
        for i in j..=n {
            let mut sum = lu.value(i, j);
            for k in 1..j {
                sum -= lu.value(i, k) * lu.value(k, j);
            }
            lu.set_value(i, j, sum);
            let dum = vv[i] * sum.abs();
            if dum < big {
                continue;
            }
            big = dum;
            imax = i;
        }
        if j != imax {
            for k in 1..=n {
                let t = lu.value(imax, k);
                lu.set_value(imax, k, lu.value(j, k));
                lu.set_value(j, k, t);
            }
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if lu.value(j, j).abs() <= tiny {
            return LuDecomp {
                lu,
                indx,
                done: false,
            };
        }
        if j != n {
            let dum = 1.0 / lu.value(j, j);
            for i in (j + 1)..=n {
                lu.set_value(i, j, lu.value(i, j) * dum);
            }
        }
    }
    LuDecomp {
        lu,
        indx,
        done: true,
    }
}

/// `LU_Solve(a, indx, b)` (`math_Recipes.cxx:309-344`).
fn lu_solve(a: &MathMatrix, indx: &[usize], b: &mut MathVector) {
    let n = a.row_count();
    let nblow = b.lower() - 1;
    let mut ii = 0usize;
    for i in 1..=n {
        let ip = indx[i];
        let mut sum = b.value(ip + nblow);
        let t = b.value(i + nblow);
        b.set_value(ip + nblow, t);
        if ii != 0 {
            for j in ii..i {
                sum -= a.value(i, j) * b.value(j + nblow);
            }
        } else if sum != 0.0 {
            ii = i;
        }
        b.set_value(i + nblow, sum);
    }
    for i in (1..=n).rev() {
        let mut sum = b.value(i + nblow);
        for j in (i + 1)..=n {
            sum -= a.value(i, j) * b.value(j + nblow);
        }
        b.set_value(i + nblow, sum / a.value(i, i));
    }
}

/// `math_SVD` constructor (`math_SVD.cxx:29-39`): pads A with zero rows to
/// `max(RowNumber, ColNumber)` before calling `SVD_Decompose`.
fn svd_padded(a: &MathMatrix) -> SVD {
    let m = a.row_count();
    let n = a.col_count();
    let size = m.max(n);
    let mut u = MathMatrix::new(1, size, 1, n);
    for i in 1..=m {
        for j in 1..=n {
            u.set_value(i, j, a.value(i, j));
        }
    }
    SVD::new(&u)
}

/// `math_SVD::Solve` (`math_SVD.cxx:41-68`) and `SVD_Solve`
/// (`math_Recipes.cxx:712-748`), with the default `Eps = 1.0e-6`
/// (`math_SVD.hxx:53`).
fn svd_solve_occt(svd: &SVD, b: &MathVector, x: &mut MathVector) {
    let m = svd.u.row_count();
    let n = svd.v.row_count();
    let mut wmax = 0.0f64;
    for i in 1..=svd.w.len() {
        wmax = wmax.max(svd.w.value(i));
    }
    let wmin = 1.0e-6 * wmax;
    let mut w = svd.w.clone();
    for i in 1..=w.len() {
        if w.value(i) < wmin {
            w.set_value(i, 0.0);
        }
    }
    // BB(1..U.RowNumber()) = B padded with zeros.
    let mut bb = MathVector::new(1, m);
    for i in 1..=b.len() {
        bb.set_value(i, b.value(i));
    }
    let mut tmp = MathVector::new(1, n);
    for j in 1..=n {
        let mut s = 0.0;
        if w.value(j) != 0.0 {
            for i in 1..=m {
                s += svd.u.value(i, j) * bb.value(i);
            }
            s /= w.value(j);
        }
        tmp.set_value(j, s);
    }
    for j in 1..=n {
        let mut s = 0.0;
        for jj in 1..=n {
            s += svd.v.value(j, jj) * tmp.value(jj);
        }
        x.set_value(j, s);
    }
}

/// `math_GaussLeastSquare` (`math_GaussLeastSquare.cxx:28-49`): LU of
/// `A^T A` plus the `A^T B` right-hand side.
struct GaussLeastSquare {
    lu: MathMatrix,
    indx: Vec<usize>,
    a2: MathMatrix,
    done: bool,
}

impl GaussLeastSquare {
    fn new(a: &MathMatrix, min_pivot: f64) -> Self {
        let neq = a.row_count();
        let ninc = a.col_count();
        // A2 = A.Transposed() (ninc x neq).
        let mut a2 = MathMatrix::new(1, ninc, 1, neq);
        for i in 1..=neq {
            for j in 1..=ninc {
                a2.set_value(j, i, a.value(i, j));
            }
        }
        let mut lu = MathMatrix::new(1, ninc, 1, ninc);
        lu.multiply_mat(&a2, a); // LU = A2 * A
        let d = lu_decompose(&lu, min_pivot);
        Self {
            lu: d.lu,
            indx: d.indx,
            a2,
            done: d.done,
        }
    }

    /// `Solve(B, X)` (`math_GaussLeastSquare.cxx:40-49`).
    fn solve(&self, b: &MathVector, x: &mut MathVector) {
        x.multiply_mat_vec(&self.a2, b); // X = A2 * B
        lu_solve(&self.lu, &self.indx, x);
    }
}

// ---------------------------------------------------------------------------
// MyDirFunction (math_FunctionSetRoot.cxx:70-195)
// ---------------------------------------------------------------------------

/// `MyDirFunction` (`math_FunctionSetRoot.cxx:70-195`): a one-dimensional
/// restriction of the function set along a direction, plus the
/// value/derivative pack used by `Perform`.
struct DirFunction {
    p0: MathVector,
    dir: MathVector,
    p: MathVector,
    fv: MathVector,
}

impl DirFunction {
    fn new(ninc: usize, neq: usize) -> Self {
        Self {
            p0: MathVector::new(1, ninc),
            dir: MathVector::new(1, ninc),
            p: MathVector::new(1, ninc),
            fv: MathVector::new(1, neq),
        }
    }

    /// `MyDirFunction::Initialize` (`math_FunctionSetRoot.cxx:114-118`).
    fn initialize(&mut self, p0: &MathVector, dir: &MathVector) {
        self.p0 = p0.clone();
        self.dir = dir.clone();
    }

    /// `MyDirFunction::Value(x, fval)` (`math_FunctionSetRoot.cxx:120-150`).
    /// `None` mirrors OCCT's `false` return.
    fn value_scalar<F: MathFunctionSetWithDerivatives>(
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
    fn value_full<F: MathFunctionSetWithDerivatives>(
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
fn search_direction(
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
fn search_direction_constrained(
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
        let mut my_gh = MathVector::new(1, ninc - cons);
        let mut my_direction = MathVector::new(1, ninc - cons);
        let mut my_inv_length_max = MathVector::new(1, ninc);
        let mut k = 1;
        for i in 1..=ninc {
            if constraints[i] == 0 {
                my_gh.set_value(k, gh.value(i));
                my_inv_length_max.set_value(k, inv_length_max.value(i));
                my_direction.set_value(k, direction.value(i));
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
            &my_gh,
            ff,
            change_direction,
            &my_inv_length_max,
            &mut my_direction,
            dy,
        );
        let mut k = 1;
        for i in 1..=ninc {
            if constraints[i] == 0 {
                if !change_direction {
                    direction.set_value(i, my_direction.value(k));
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
fn bounds(
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

// ---------------------------------------------------------------------------
// MinimizeDirection (math_FunctionSetRoot.cxx:198-436)
// ---------------------------------------------------------------------------

/// `MinimizeDirection(P0, P1, P2, F1, Delta, Tol, F)`
/// (`math_FunctionSetRoot.cxx:198-264`): minimisation from three points.
fn minimize_direction_3pt<F: MathFunctionSetWithDerivatives>(
    p0: &MathVector,
    p1: &MathVector,
    p2: &MathVector,
    f1: f64,
    delta: &mut MathVector,
    tol: &MathVector,
    dir: &mut DirFunction,
    f: &mut F,
) -> bool {
    // (1) 1D parametric tolerance.
    let mut tol1d: f64 = 2.1;
    let eps = 1.0e-16;
    for ii in 1..=tol.len() {
        let invnorme = delta.value(ii).abs();
        if invnorme > eps {
            tol1d = tol1d.min(tol.value(ii) / invnorme);
        }
    }
    if tol1d > 1.9 {
        return false; // Pas la peine de se fatiguer
    }
    tol1d /= 3.0;

    *delta = p1.subtracted(p0);
    let mut invnorme = delta.norm();
    if invnorme <= eps {
        return false;
    }
    invnorme = 1.0 / invnorme;

    dir.initialize(p1, delta);

    // (2) On minimise.
    let ax = -1.0;
    let bx = 0.0;
    let cx = p2.subtracted(p1).norm() * invnorme;
    if cx < 1.0e-2 {
        return false;
    }

    let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
    let mut closure = |x: f64| dir.value_scalar(f, x);
    sol.perform_fallible(&mut closure, ax, bx, cx);

    if sol.is_done() {
        let tsol = sol.location().unwrap();
        if sol.minimum().unwrap() < f1 {
            delta.multiply_scalar(tsol);
            return true;
        }
    }
    false
}

/// `MinimizeDirection(P, Dir, PValue, PDirValue, Gradient, DGradient, Tol, F)`
/// (`math_FunctionSetRoot.cxx:267-436`): minimisation from two points and a
/// derivative.
fn minimize_direction_2pt<F: MathFunctionSetWithDerivatives>(
    p: &MathVector,
    dir_vec: &mut MathVector,
    p_value: f64,
    p_dir_value: f64,
    gradient: &MathVector,
    d_gradient: &MathVector,
    tol: &MathVector,
    dir: &mut DirFunction,
    f: &mut F,
) -> bool {
    if is_infinite(p_value) || is_infinite(p_dir_value) {
        return false;
    }
    // (0) 1D parametric tolerance.
    let mut good = false;
    let eps = 1.0e-20;
    let mut tol1d: f64 = 1.1;
    let mut result = p_value;
    for ii in 1..=tol.len() {
        let absdir = dir_vec.value(ii).abs();
        if absdir > eps {
            tol1d = tol1d.min(tol.value(ii) / absdir);
        }
    }
    if tol1d > 0.9 {
        return false;
    }

    // (1) Premiere interpolation quadratique.
    let df1 = gradient.dot(dir_vec);
    let df2 = d_gradient.dot(dir_vec);
    let mut tsol;
    if df1 < -eps && df2 > eps {
        // cuvette
        tsol = -df1 / (df2 - df1);
    } else {
        let cx = p_value;
        let bx = df1;
        let ax = p_dir_value - (bx + cx);
        if ax.abs() <= eps {
            // cas lineaire
            if bx.abs() >= eps {
                tsol = -cx / bx;
            } else {
                tsol = 0.0;
            }
        } else {
            // cas quadratique
            let mut disc = bx * bx - 4.0 * ax * cx;
            if disc > 1.0e-9 {
                disc = disc.sqrt();
                tsol = -(bx + disc);
                let tsolbis = disc - bx;
                if tsolbis.abs() < tsol.abs() {
                    tsol = tsolbis;
                }
                tsol /= 2.0 * ax;
            } else {
                // pas ou peu de racine : on "extremise"
                tsol = -(0.5 * bx) / ax;
            }
        }
    }

    if tsol.abs() >= 1.0 {
        return false; // resultat sans interet
    }

    dir.initialize(p, dir_vec);
    // math_FunctionSetRoot.cxx:358 ignores the bool returned by
    // MyDirFunction::Value, leaving fsol indeterminate on failure. We model the
    // failure as "no value": the progress test below is skipped and the Brent
    // search is entered.
    let fsol = dir.value_scalar(f, tsol);
    if let Some(v) = fsol {
        if v < p_value {
            good = true;
            result = v;
        }
    }

    // (2) Recherche en bonne et due forme si l'on n'a pas assez progresse.
    let need_brent = fsol.map_or(true, |v| v > 0.2 * p_value) && tol1d < 0.5;
    if need_brent {
        let (ax, bx, cx);
        if tsol < 0.0 {
            ax = tsol;
            bx = 0.0;
            cx = 1.0;
        } else {
            ax = 0.0;
            bx = tsol;
            cx = 1.0;
        }
        let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
        let mut closure = |x: f64| dir.value_scalar(f, x);
        sol.perform_fallible(&mut closure, ax, bx, cx);
        if sol.is_done() {
            if sol.minimum().unwrap() <= result {
                tsol = sol.location().unwrap();
                good = true;
                result = sol.minimum().unwrap();

                // Objective function changes too fast -> extra computations.
                if gradient.norm2() > 1.0 / SQUARE_CONFUSION && tsol > ax && tsol < cx {
                    sol.perform_fallible(&mut closure, ax, (ax + tsol) / 2.0, tsol);
                    if sol.is_done() && sol.minimum().unwrap() <= result {
                        tsol = sol.location().unwrap();
                        good = true;
                        result = sol.minimum().unwrap();
                    }
                    sol.perform_fallible(&mut closure, tsol, (cx + tsol) / 2.0, cx);
                    if sol.is_done() && sol.minimum().unwrap() <= result {
                        tsol = sol.location().unwrap();
                        good = true;
                        result = sol.minimum().unwrap();
                    }
                }
            }
        }
    }

    if good {
        // mise a jour du Delta
        dir_vec.multiply_scalar(tsol);
    }
    good
}

// ---------------------------------------------------------------------------
// math_FunctionSetRoot (math_FunctionSetRoot.hxx:36-202, .cxx:709-1452)
// ---------------------------------------------------------------------------

/// `math_FunctionSetRoot` (`math_FunctionSetRoot.hxx:36-202`).
///
/// The protected OCCT fields `Delta`/`Sol`/`DF`/`Tol` are public here; the
/// private fields keep OCCT's names.
pub struct MathFunctionSetRoot {
    pub delta: MathVector,
    pub sol: MathVector,
    pub df: MathMatrix,
    pub tol: MathVector,
    done: bool,
    kount: i32,
    state: i32,
    itermax: i32,
    inf_bound: MathVector,
    sup_bound: MathVector,
    sol_save: MathVector,
    gh: MathVector,
    dh: MathVector,
    dhsave: MathVector,
    ff: MathVector,
    previous_solution: MathVector,
    /// `Save(0, theNbIterations)` (`math_FunctionSetRoot.cxx:729`); stored
    /// 0-based, index `k` is OCCT's `Save(k)`.
    save: Vec<f64>,
    my_is_divergent: bool,
}

impl MathFunctionSetRoot {
    /// `math_FunctionSetRoot(F, Tolerance, NbIterations = 100)`
    /// (`math_FunctionSetRoot.cxx:709-738`).
    pub fn new<F: MathFunctionSetWithDerivatives>(
        f: &F,
        tolerance: &MathVector,
        nb_iterations: i32,
    ) -> Self {
        let mut s = Self::with_iterations(f, nb_iterations);
        s.set_tolerance(tolerance);
        s
    }

    /// `math_FunctionSetRoot(F, NbIterations = 100)`
    /// (`math_FunctionSetRoot.cxx:742-769`): the tolerance must be set
    /// separately with `SetTolerance`.
    pub fn with_iterations<F: MathFunctionSetWithDerivatives>(
        f: &F,
        nb_iterations: i32,
    ) -> Self {
        let ninc = f.nb_variables();
        let neq = f.nb_equations();
        let save_len = if nb_iterations >= 0 {
            nb_iterations as usize + 1
        } else {
            1
        };
        Self {
            delta: MathVector::new(1, ninc),
            sol: MathVector::new(1, ninc),
            df: MathMatrix::new(1, neq, 1, ninc),
            tol: MathVector::new(1, ninc),
            done: false,
            kount: 0,
            state: 0,
            itermax: nb_iterations,
            inf_bound: MathVector::with_init(1, ninc, REAL_FIRST),
            sup_bound: MathVector::with_init(1, ninc, REAL_LAST),
            sol_save: MathVector::new(1, ninc),
            gh: MathVector::new(1, ninc),
            dh: MathVector::new(1, ninc),
            dhsave: MathVector::new(1, ninc),
            ff: MathVector::new(1, neq),
            previous_solution: MathVector::new(1, ninc),
            save: vec![0.0; save_len],
            my_is_divergent: false,
        }
    }

    /// `SetTolerance` (`math_FunctionSetRoot.cxx:777-783`).
    pub fn set_tolerance(&mut self, tolerance: &MathVector) {
        for i in 1..=self.tol.len() {
            self.tol.set_value(i, tolerance.value(i));
        }
    }

    /// `IsDone` (`math_FunctionSetRoot.hxx:100`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbIterations` (`math_FunctionSetRoot.hxx:105-109`); panics where OCCT
    /// raises `StdFail_NotDone`.
    pub fn nb_iterations(&self) -> i32 {
        assert!(self.done, "math_FunctionSetRoot: NotDone");
        self.kount
    }

    /// `StateNumber` (`math_FunctionSetRoot.hxx:113-117`).
    pub fn state_number(&self) -> i32 {
        assert!(self.done, "math_FunctionSetRoot: NotDone");
        self.state
    }

    /// `Root()` (`math_FunctionSetRoot.hxx:121-125`).
    pub fn root(&self) -> MathVector {
        assert!(self.done, "math_FunctionSetRoot: NotDone");
        self.sol.clone()
    }

    /// `Derivative()` (`math_FunctionSetRoot.hxx:135-139`).
    pub fn derivative(&self) -> MathMatrix {
        assert!(self.done, "math_FunctionSetRoot: NotDone");
        self.df.clone()
    }

    /// `FunctionSetErrors()` (`math_FunctionSetRoot.hxx:156-160`).
    pub fn function_set_errors(&self) -> MathVector {
        assert!(self.done, "math_FunctionSetRoot: NotDone");
        self.delta.clone()
    }

    /// `IsDivergent` (`math_FunctionSetRoot.hxx:174`).
    pub fn is_divergent(&self) -> bool {
        self.my_is_divergent
    }

    /// `IsSolutionReached` base implementation
    /// (`math_FunctionSetRoot.hxx:69-79`).
    pub fn is_solution_reached(&self) -> bool {
        for i in 1..=self.sol.len() {
            if self.delta.value(i).abs() > self.tol.value(i) {
                return false;
            }
        }
        true
    }

    /// `Perform(F, StartingPoint, StopOnDivergent = false)`
    /// (`math_FunctionSetRoot.cxx:787-792`).
    pub fn perform<F: MathFunctionSetWithDerivatives>(
        &mut self,
        f: &mut F,
        starting_point: &MathVector,
        stop_on_divergent: bool,
    ) {
        let inf = self.inf_bound.clone();
        let sup = self.sup_bound.clone();
        self.perform_with_bounds(f, starting_point, &inf, &sup, stop_on_divergent);
    }

    /// `Perform(F, StartingPoint, InfBound, SupBound, StopOnDivergent = false)`
    /// (`math_FunctionSetRoot.cxx:796-1417`).
    pub fn perform_with_bounds<F: MathFunctionSetWithDerivatives>(
        &mut self,
        f: &mut F,
        starting_point: &MathVector,
        the_inf_bound: &MathVector,
        the_sup_bound: &MathVector,
        the_stop_on_divergent: bool,
    ) {
        let ninc = f.nb_variables();
        let neq = f.nb_equations();

        if neq == 0
            || starting_point.len() != ninc
            || the_inf_bound.len() != ninc
            || the_sup_bound.len() != ninc
        {
            // math_FunctionSetRoot.cxx:804-808 throws Standard_DimensionError.
            panic!("math_FunctionSetRoot::Perform: dimension error");
        }

        let epsk = 1.0e-16; // EpsSqrt
        let eps = 1.0e-32; // Eps
        let eps2 = 1.0e-64; // Eps2
        let progres = 0.005; // Progres

        let mut change_direction = false;
        let mut sort = false;
        let mut is_new_sol = false;

        // InvLengthMax / aConstraints (math_FunctionSetRoot.cxx:817-824).
        let mut inv_length_max = MathVector::new(1, ninc);
        for i in 1..=ninc {
            let a_sup_bound = the_sup_bound.value(i).min(PRECISION_INFINITE);
            let an_inf_bound = the_inf_bound.value(i).max(-PRECISION_INFINITE);
            inv_length_max.set_value(
                i,
                1.0 / ((a_sup_bound - an_inf_bound) / 4.0).max(1.0e-9),
            );
        }
        let mut a_constraints = vec![0i32; ninc + 1];

        let mut dir = DirFunction::new(ninc, neq);

        self.done = false;
        self.sol = starting_point.clone();
        self.kount = 0;

        self.my_is_divergent = false;
        for i in 1..=ninc {
            self.my_is_divergent = self.my_is_divergent
                || self.sol.value(i) < the_inf_bound.value(i)
                || self.sol.value(i) > the_sup_bound.value(i);
        }
        if the_stop_on_divergent && self.my_is_divergent {
            return;
        }

        // Recentrage sur les bornes (math_FunctionSetRoot.cxx:846-856).
        for i in 1..=ninc {
            if self.sol.value(i) <= the_inf_bound.value(i) {
                self.sol.set_value(i, the_inf_bound.value(i));
            } else if self.sol.value(i) > the_sup_bound.value(i) {
                self.sol.set_value(i, the_sup_bound.value(i));
            }
        }

        // Premiere valeur de F et de son gradient.
        let (mut f2, mut gnr1) = match dir.value_full(
            f,
            &self.sol,
            &mut self.ff,
            &mut self.df,
            &mut self.gh,
        ) {
            Some(v) => v,
            None => {
                self.done = false;
                if !the_stop_on_divergent || !self.my_is_divergent {
                    self.state = f.get_state_number();
                }
                return;
            }
        };
        let mut ambda2 = gnr1;
        self.save[0] = f2.max(epsk);
        let a_tol_func = epsilon_of(f2);

        if f2 <= eps || gnr1 <= eps2 {
            self.done = false;
            if !the_stop_on_divergent || !self.my_is_divergent {
                self.done = true;
                self.state = f.get_state_number();
            }
            return;
        }

        let mut dy = 0.0f64;
        self.kount = 1;
        while self.kount <= self.itermax {
            let previous_minimum = f2;
            let oldgr = gnr1;
            self.previous_solution = self.sol.clone();
            self.sol_save = self.sol.clone();

            // math_FunctionSetRoot.cxx:895 — ChangeDirection is passed by value
            // to SearchDirection and is not updated by it.
            search_direction(
                &self.df,
                &self.gh,
                &self.ff,
                change_direction,
                &inv_length_max,
                &mut self.dh,
                &mut dy,
            );
            if dy.abs() <= eps {
                self.done = false;
                if !the_stop_on_divergent || !self.my_is_divergent {
                    self.done = true;
                    // modified by jgv, 31.08.2011: update F before GetStateNumber
                    f.value(&self.sol, &mut self.ff);
                    self.state = f.get_state_number();
                }
                return;
            }
            let mut ambda;
            if change_direction {
                ambda = ambda2 / dy.abs().sqrt();
                if ambda > 1.0 {
                    ambda = 1.0;
                }
            } else {
                ambda = 1.0;
                ambda2 = 0.5 * ambda / self.dh.norm();
            }

            for i in 1..=ninc {
                self.sol
                    .set_value(i, self.sol.value(i) + ambda * self.dh.value(i));
            }
            for i in 1..=ninc {
                self.my_is_divergent = self.my_is_divergent
                    || self.sol.value(i) < the_inf_bound.value(i)
                    || self.sol.value(i) > the_sup_bound.value(i);
            }
            if the_stop_on_divergent && self.my_is_divergent {
                return;
            }

            let br = bounds(
                the_inf_bound,
                the_sup_bound,
                &self.tol,
                &mut self.sol,
                &self.sol_save,
                &mut a_constraints,
                &mut self.delta,
            );
            sort = br.0;
            is_new_sol = br.1;

            self.dhsave = self.gh.clone();
            if is_new_sol {
                match dir.value_full(
                    f,
                    &self.sol,
                    &mut self.ff,
                    &mut self.df,
                    &mut self.gh,
                ) {
                    Some(v) => {
                        f2 = v.0;
                        gnr1 = v.1;
                    }
                    None => {
                        self.done = false;
                        if !the_stop_on_divergent || !self.my_is_divergent {
                            self.state = f.get_state_number();
                        }
                        return;
                    }
                }
            }

            if f2 <= eps || gnr1 <= eps2 {
                self.done = false;
                if !the_stop_on_divergent || !self.my_is_divergent {
                    self.done = true;
                    f.value(&self.sol, &mut self.ff);
                    self.state = f.get_state_number();
                }
                return;
            }

            if sort || (f2 / previous_minimum > progres) {
                dy = self.gh.dot(&self.dh);
                let mut old_f = previous_minimum;
                let mut stop = false;
                let mut good = false;
                let mut descente_iter = 0;
                let mut sort_bis;

                // Standard processing without boundary handling.
                if !sort {
                    while (f2 / previous_minimum > progres) && !stop {
                        if f2 < old_f && dy < 0.0 {
                            descente_iter += 1;
                            self.sol_save = self.sol.clone();
                            old_f = f2;
                            for i in 1..=ninc {
                                self.sol.set_value(
                                    i,
                                    self.sol.value(i) + ambda * self.dh.value(i),
                                );
                            }
                            for i in 1..=ninc {
                                self.my_is_divergent = self.my_is_divergent
                                    || self.sol.value(i) < the_inf_bound.value(i)
                                    || self.sol.value(i) > the_sup_bound.value(i);
                            }
                            if the_stop_on_divergent && self.my_is_divergent {
                                return;
                            }
                            let br = bounds(
                                the_inf_bound,
                                the_sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut a_constraints,
                                &mut self.delta,
                            );
                            stop = br.0;
                            is_new_sol = br.1;
                            ambda *= 1.7;
                        } else {
                            if f2 >= old_f || f2 >= previous_minimum {
                                good = false;
                                if descente_iter == 0 {
                                    // C'est le premier pas qui flanche.
                                    descente_iter += 1;
                                    good = minimize_direction_2pt(
                                        &self.sol_save,
                                        &mut self.delta,
                                        old_f,
                                        f2,
                                        &self.dhsave,
                                        &self.gh,
                                        &self.tol,
                                        &mut dir,
                                        f,
                                    );
                                } else if change_direction
                                    || descente_iter > 1
                                    || old_f > previous_minimum
                                {
                                    descente_iter += 1;
                                    good = minimize_direction_3pt(
                                        &self.previous_solution,
                                        &self.sol_save,
                                        &self.sol,
                                        old_f,
                                        &mut self.delta,
                                        &self.tol,
                                        &mut dir,
                                        f,
                                    );
                                }
                                if !good {
                                    self.sol = self.sol_save.clone();
                                    f2 = old_f;
                                } else {
                                    self.sol = self.sol_save.added(&self.delta);
                                    for i in 1..=ninc {
                                        self.my_is_divergent = self.my_is_divergent
                                            || self.sol.value(i) < the_inf_bound.value(i)
                                            || self.sol.value(i) > the_sup_bound.value(i);
                                    }
                                    if the_stop_on_divergent && self.my_is_divergent {
                                        return;
                                    }
                                    let br = bounds(
                                        the_inf_bound,
                                        the_sup_bound,
                                        &self.tol,
                                        &mut self.sol,
                                        &self.sol_save,
                                        &mut a_constraints,
                                        &mut self.delta,
                                    );
                                    sort = br.0;
                                    is_new_sol = br.1;
                                }
                                sort = false; // On a rejete le point sur la frontiere
                            }
                            stop = true;
                        }
                        self.dhsave = self.gh.clone();
                        if is_new_sol {
                            match dir.value_full(
                                f,
                                &self.sol,
                                &mut self.ff,
                                &mut self.df,
                                &mut self.gh,
                            ) {
                                Some(v) => {
                                    f2 = v.0;
                                    gnr1 = v.1;
                                }
                                None => {
                                    self.done = false;
                                    if !the_stop_on_divergent || !self.my_is_divergent {
                                        self.state = f.get_state_number();
                                    }
                                    return;
                                }
                            }
                        }
                        dy = self.gh.dot(&self.dh);
                        if dy.abs() <= eps {
                            if f2 > old_f {
                                self.sol = self.sol_save.clone();
                            }
                            self.done = false;
                            if !the_stop_on_divergent || !self.my_is_divergent {
                                self.done = true;
                                f.value(&self.sol, &mut self.ff);
                                self.state = f.get_state_number();
                            }
                            return;
                        }
                        if descente_iter >= 100 {
                            stop = true;
                        }
                    }
                }
                // ------------------------------------
                //  on passe au traitement des bords
                // ------------------------------------
                if sort {
                    stop = f2 > 1.001 * old_f; // Pour ne pas progresser sur le bord
                    sort_bis = sort;
                    descente_iter = 0;
                    while sort_bis && (f2 < old_f || descente_iter == 0) && !stop {
                        descente_iter += 1;
                        // On essaye de progresser sur le bord.
                        self.sol_save = self.sol.clone();
                        old_f = f2;
                        search_direction_constrained(
                            &self.df,
                            &self.gh,
                            &self.ff,
                            &a_constraints,
                            &self.sol,
                            change_direction,
                            &inv_length_max,
                            &mut self.dh,
                            &mut dy,
                        );
                        if dy < -eps {
                            // Pour eviter des calculs inutiles et des /0...
                            if change_direction {
                                ambda = ambda2 / (-dy).sqrt();
                                if ambda > 1.0 {
                                    ambda = 1.0;
                                }
                            } else {
                                ambda = 1.0;
                                ambda2 = 0.5 * ambda / self.dh.norm();
                            }
                            for i in 1..=ninc {
                                self.sol.set_value(
                                    i,
                                    self.sol.value(i) + ambda * self.dh.value(i),
                                );
                            }
                            for i in 1..=ninc {
                                self.my_is_divergent = self.my_is_divergent
                                    || self.sol.value(i) < the_inf_bound.value(i)
                                    || self.sol.value(i) > the_sup_bound.value(i);
                            }
                            if the_stop_on_divergent && self.my_is_divergent {
                                return;
                            }
                            let br = bounds(
                                the_inf_bound,
                                the_sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut a_constraints,
                                &mut self.delta,
                            );
                            sort_bis = br.0;
                            is_new_sol = br.1;

                            self.dhsave = self.gh.clone();
                            if is_new_sol {
                                match dir.value_full(
                                    f,
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                ) {
                                    Some(v) => {
                                        f2 = v.0;
                                        gnr1 = v.1;
                                    }
                                    None => {
                                        self.done = false;
                                        if !the_stop_on_divergent || !self.my_is_divergent {
                                            self.state = f.get_state_number();
                                        }
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                        } else {
                            stop = true;
                        }

                        while (f2 / previous_minimum > progres) && (f2 < old_f) && !stop {
                            descente_iter += 1;
                            if f2 < old_f && dy < 0.0 {
                                // On essaye de progresser dans cette direction.
                                self.sol_save = self.sol.clone();
                                old_f = f2;
                                for i in 1..=ninc {
                                    self.sol.set_value(
                                        i,
                                        self.sol.value(i) + ambda * self.dh.value(i),
                                    );
                                }
                                for i in 1..=ninc {
                                    self.my_is_divergent = self.my_is_divergent
                                        || self.sol.value(i) < the_inf_bound.value(i)
                                        || self.sol.value(i) > the_sup_bound.value(i);
                                }
                                if the_stop_on_divergent && self.my_is_divergent {
                                    return;
                                }
                                let br = bounds(
                                    the_inf_bound,
                                    the_sup_bound,
                                    &self.tol,
                                    &mut self.sol,
                                    &self.sol_save,
                                    &mut a_constraints,
                                    &mut self.delta,
                                );
                                sort_bis = br.0;
                                is_new_sol = br.1;
                            }
                            self.dhsave = self.gh.clone();
                            if is_new_sol {
                                match dir.value_full(
                                    f,
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                ) {
                                    Some(v) => {
                                        f2 = v.0;
                                        gnr1 = v.1;
                                    }
                                    None => {
                                        self.done = false;
                                        if !the_stop_on_divergent || !self.my_is_divergent {
                                            self.state = f.get_state_number();
                                        }
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                            dy = self.gh.dot(&self.dh);
                            stop = dy >= 0.0 || descente_iter >= 10 || sort_bis;
                        }
                        stop = dy >= 0.0 || descente_iter >= 10;
                    }
                    if ((f2 / previous_minimum > progres) && (f2 >= old_f))
                        || (f2 >= previous_minimum)
                    {
                        // On minimise par Brent.
                        descente_iter += 1;
                        good = minimize_direction_2pt(
                            &self.sol_save,
                            &mut self.delta,
                            old_f,
                            f2,
                            &self.dhsave,
                            &self.gh,
                            &self.tol,
                            &mut dir,
                            f,
                        );
                        if !good {
                            self.sol = self.sol_save.clone();
                            sort = false;
                        } else {
                            self.sol = self.sol_save.added(&self.delta);
                            for i in 1..=ninc {
                                self.my_is_divergent = self.my_is_divergent
                                    || self.sol.value(i) < the_inf_bound.value(i)
                                    || self.sol.value(i) > the_sup_bound.value(i);
                            }
                            if the_stop_on_divergent && self.my_is_divergent {
                                return;
                            }
                            let br = bounds(
                                the_inf_bound,
                                the_sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut a_constraints,
                                &mut self.delta,
                            );
                            sort = br.0;
                            is_new_sol = br.1;
                            if is_new_sol {
                                match dir.value_full(
                                    f,
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                ) {
                                    Some(v) => {
                                        f2 = v.0;
                                        gnr1 = v.1;
                                    }
                                    None => {
                                        self.done = false;
                                        if !the_stop_on_divergent || !self.my_is_divergent {
                                            self.state = f.get_state_number();
                                        }
                                        return;
                                    }
                                }
                            }
                        }
                        dy = self.gh.dot(&self.dh);
                    }
                }

            }
            // ---------------------------------------------
            //  on passe aux tests d'ARRET
            // ---------------------------------------------
            self.save[self.kount as usize] = f2;
            let verif: bool;
            if change_direction {
                verif = true;
                // Gradient : Il faut eviter de boucler.
            } else {
                if self.kount > 1 {
                    // Pour accelerer les cas quasi-quadratique.
                    verif = self.save[(self.kount - 1) as usize]
                        < 1.0e-4 * self.save[(self.kount - 2) as usize];
                } else {
                    // Pour les cas dejas solutions.
                    verif = f2 < 1.0e-6 * self.save[0];
                }
            }
            if verif {
                for i in 1..=ninc {
                    self.delta.set_value(
                        i,
                        self.previous_solution.value(i) - self.sol.value(i),
                    );
                }
                if self.is_solution_reached() {
                    if previous_minimum < f2 {
                        self.sol = self.sol_save.clone();
                    }
                    self.done = false;
                    if !the_stop_on_divergent || !self.my_is_divergent {
                        self.done = true;
                        f.value(&self.sol, &mut self.ff);
                        self.state = f.get_state_number();
                    }
                    return;
                }
            }

            // Analyse de la progression...
            if (f2 - previous_minimum) <= a_tol_func {
                if self.kount > 5 {
                    // L'historique est il bon ?
                    if f2 >= 0.95 * self.save[(self.kount - 5) as usize] {
                        if !change_direction {
                            change_direction = true;
                        } else {
                            self.done = false;
                            if !the_stop_on_divergent || !self.my_is_divergent {
                                self.done = true;
                                self.state = f.get_state_number();
                            }
                            // si un gain inf a 5% on sort
                            return;
                        }
                    } else {
                        change_direction = false; // If yes we restart
                    }
                } else {
                    change_direction = false; // No history, we continue
                }
                // If the gradient does not decrease sufficiently with Newton, we
                // try the gradient method unless f decreases.
                if (gnr1 > 0.9 * oldgr) && (f2 > 0.5 * previous_minimum) {
                    change_direction = true;
                }

                if !change_direction && !verif {
                    for i in 1..=ninc {
                        self.delta.set_value(
                            i,
                            self.previous_solution.value(i) - self.sol.value(i),
                        );
                    }
                    if self.is_solution_reached() {
                        self.done = false;
                        if !the_stop_on_divergent || !self.my_is_divergent {
                            self.done = true;
                            f.value(&self.sol, &mut self.ff);
                            self.state = f.get_state_number();
                        }
                        return;
                    }
                }
            } else {
                // Cas de regression.
                if !change_direction {
                    // On passe au gradient.
                    change_direction = true;
                    self.sol = self.previous_solution.clone();
                    match dir.value_full(
                        f,
                        &self.sol,
                        &mut self.ff,
                        &mut self.df,
                        &mut self.gh,
                    ) {
                        Some(v) => {
                            f2 = v.0;
                            gnr1 = v.1;
                        }
                        None => {
                            self.done = false;
                            if !the_stop_on_divergent || !self.my_is_divergent {
                                self.state = f.get_state_number();
                            }
                            return;
                        }
                    }
                } else {
                    if !the_stop_on_divergent || !self.my_is_divergent {
                        self.state = f.get_state_number();
                    }
                    // y a plus d'issues
                    return;
                }
            }

            self.kount += 1;
        }
        if !the_stop_on_divergent || !self.my_is_divergent {
            self.state = f.get_state_number();
        }
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
    the_root: f64,
    the_error: f64,
    the_derivative: f64,
    nb_iter: i32,
}

impl MathFunctionRoot {
    /// `math_FunctionRoot(F, Guess, Tolerance, NbIterations = 100)`
    /// (`math_FunctionRoot.cxx:73-93`). The unbounded `Perform` uses
    /// `InfBound = RealFirst()`, `SupBound = RealLast()`
    /// (`math_FunctionSetRoot.cxx:721-722`).
    pub fn new(
        f: &mut dyn MathFunctionWithDerivative,
        guess: f64,
        tolerance: f64,
        nb_iterations: i32,
    ) -> Self {
        Self::new_with_bounds(f, guess, tolerance, REAL_FIRST, REAL_LAST, nb_iterations)
    }

    /// `math_FunctionRoot(F, Guess, Tolerance, A, B, NbIterations = 100)`
    /// (`math_FunctionRoot.cxx:95-119`).
    pub fn new_with_bounds(
        f: &mut dyn MathFunctionWithDerivative,
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

        let (done, nb_iter, the_root, the_derivative) = {
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
            f.get_state_number();
            let mut the_error = 0.0;
            f.value(the_root, &mut the_error);
            Self {
                done: true,
                the_root,
                the_error,
                the_derivative,
                nb_iter,
            }
        } else {
            Self {
                done: false,
                the_root: 0.0,
                the_error: 0.0,
                the_derivative: 0.0,
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
        self.the_root
    }

    /// `Derivative` (`math_FunctionRoot.lxx:34-38`).
    pub fn derivative(&self) -> f64 {
        assert!(self.done, "math_FunctionRoot: NotDone");
        self.the_derivative
    }

    /// `Value` (`math_FunctionRoot.lxx:40-44`).
    pub fn value(&self) -> f64 {
        assert!(self.done, "math_FunctionRoot: NotDone");
        self.the_error
    }

    /// `NbIterations` (`math_FunctionRoot.lxx:46-50`).
    pub fn nb_iterations(&self) -> i32 {
        assert!(self.done, "math_FunctionRoot: NotDone");
        self.nb_iter
    }
}







