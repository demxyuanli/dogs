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

use crate::{MathMatrix, MathVector};

mod direction;
mod function_root;
mod linear;
mod minimize;
mod perform;

pub use self::function_root::MathFunctionRoot;

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
    fn state_number(&mut self) -> i32 {
        0
    }
}

/// `math_FunctionWithDerivative` (`math_FunctionWithDerivative.hxx:31-52`) plus
/// `math_Function::GetStateNumber` (`math_Function.hxx:57`).
pub trait FunctionWithDerivative {
    fn value(&mut self, x: f64, f: &mut f64) -> bool;
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool;
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool;
    fn state_number(&mut self) -> i32 {
        0
    }
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
    is_divergent: bool,
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
            is_divergent: false,
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
        self.is_divergent
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
}
