//! `ProjLib_PrjResolve` - solve the normal-projection function set in a window.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_PrjResolve.cxx:27-172`
//! plus `ProjLib_PrjResolve.hxx:20-58` (OCCT 8.0.0). [`PrjFunc`] carries the
//! equations.
//!
//! Translation notes:
//! - the two solvers are used in OCCT's order: `math_NewtonFunctionSetRoot`
//!   first, and only when it fails `math_FunctionSetRoot` (`cxx:96-107`). The
//!   Newton run keeps OCCT's default `NbIterations = 100` (the
//!   `math_NewtonFunctionSetRoot(F, XTolerance, FTolerance, NbIterations)`
//!   overload `cxx:96` selects).
//! - `Perform`'s last parameter is unnamed in the header and unused in the
//!   body (`cxx:51`), so [`PrjResolve::perform`] keeps the name
//!   `_strict_inside` instead of inventing a branch.
//! - `Solution()` raises `StdFail_NotDone` when `!IsDone()` (`cxx:161-169`);
//!   the port answers `None`.
//! - the clamping block (`cxx:115-138`) and the `ExtraU` / `ExtraV` window test
//!   (`cxx:139`) are transcribed with the same `2 * Tol2d` margins; the
//!   residual test (`cxx:140-155`) runs only when `FuncTol > 0` **and** the
//!   Newton run failed, exactly as in OCCT.

use occt_core::gp::GpVec2d;
use occt_math::{
    MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathNewtonFunctionSetRoot, MathVector,
};

use crate::{Curve, Surface};

use super::prj_func::PrjFunc;

/// `ProjLib_PrjResolve` (`ProjLib_PrjResolve.hxx:20-58`).
pub struct PrjResolve<'a> {
    /// `myDone` (`hxx:48`).
    done: bool,
    /// `myFix` (`hxx:50`).
    fix: i32,
    /// `mySolution` (`hxx:47`).
    solution: GpVec2d,
    /// `myCurve` (`hxx:52`).
    curve: &'a dyn Curve,
    /// `mySurface` (`hxx:53`).
    surface: &'a dyn Surface,
}

impl<'a> PrjResolve<'a> {
    /// `ProjLib_PrjResolve(C, S, Fix)` (`cxx:27-39`).
    ///
    /// `Fix` outside `1..=3` raises `Standard_ConstructionError` (`cxx:32-35`).
    pub fn new(curve: &'a dyn Curve, surface: &'a dyn Surface, fix: i32) -> Result<Self, String> {
        if !(1..=3).contains(&fix) {
            return Err("ProjLib_PrjResolve: myFix must be in 1..=3 (cxx:32-35)".into());
        }
        Ok(Self {
            done: false,
            fix,
            solution: GpVec2d::new(0.0, 0.0),
            curve,
            surface,
        })
    }

    /// `Perform(t, U, V, Tol2d, Inf, Sup, FuncTol, StrictInside)`
    /// (`cxx:45-155`).
    pub fn perform(
        &mut self,
        t: f64,
        u: f64,
        v: f64,
        tol2d: (f64, f64),
        inf: (f64, f64),
        sup: (f64, f64),
        func_tol: f64,
        _strict_inside: bool,
    ) {
        self.done = false;
        let fix_val;
        let ext_u = 10.0 * tol2d.0;
        let ext_v = 10.0 * tol2d.1;
        // `ExtInf` / `ExtSup` (`cxx:60-62`).
        let ext_inf = (inf.0 - ext_u, inf.1 - ext_v);
        let ext_sup = (sup.0 + ext_u, sup.1 + ext_v);
        let b_inf = MathVector::from_slice(&[ext_inf.0, ext_inf.1]);
        let b_sup = MathVector::from_slice(&[ext_sup.0, ext_sup.1]);
        let tol = MathVector::from_slice(&[tol2d.0, tol2d.1]);

        let (start_u, start_v) = match self.fix {
            1 => {
                fix_val = t;
                (u, v)
            }
            2 => {
                fix_val = u;
                (t, v)
            }
            _ => {
                fix_val = v;
                (t, u)
            }
        };
        let start = MathVector::from_slice(&[start_u, start_v]);

        let mut f = match PrjFunc::new(self.curve, fix_val, self.surface, self.fix) {
            Ok(f) => f,
            Err(_) => return,
        };

        // `math_NewtonFunctionSetRoot SR(F, Tol, FuncTol)` (`cxx:96-97`).
        let mut sr = MathNewtonFunctionSetRoot::new(&mut f, &tol, func_tol, 100);
        sr.perform_bounded(&mut f, &start, &b_inf, &b_sup);
        if !sr.is_done() {
            // `math_FunctionSetRoot S1(F, Tol)` (`cxx:100-101`).
            let mut s1 = MathFunctionSetRoot::new(&f, &tol, 100);
            s1.perform_with_bounds(&mut f, &start, &b_inf, &b_sup, false);
            if !s1.is_done() {
                return;
            }
        }

        let (sx, sy) = f.solution();
        self.solution = GpVec2d::new(sx, sy);
        self.done = true;

        // `ExtraU` / `ExtraV` (`cxx:115-118`).
        let extra_u = 2.0 * tol2d.0;
        let extra_v = 2.0 * tol2d.1;
        let mut x = self.solution.x();
        let mut y = self.solution.y();
        if x > inf.0 - tol2d.0 && x < inf.0 {
            x = inf.0;
        }
        if x > sup.0 && x < sup.0 + tol2d.0 {
            x = sup.0;
        }
        if y > inf.1 - tol2d.1 && y < inf.1 {
            y = inf.1;
        }
        if y > sup.1 && y < sup.1 + tol2d.1 {
            y = sup.1;
        }
        self.solution = GpVec2d::new(x, y);

        if x < inf.0 - extra_u
            || x > sup.0 + extra_u
            || y < inf.1 - extra_v
            || y > sup.1 + extra_v
        {
            self.done = false;
        } else if func_tol > 0.0 {
            let xx = MathVector::from_slice(&[self.solution.x(), self.solution.y()]);
            let mut fval = MathVector::with_init(1, 2, 0.0);
            f.value(&xx, &mut fval);
            if !sr.is_done() {
                if fval.value(1) * fval.value(1) + fval.value(2) * fval.value(2) > func_tol {
                    self.done = false;
                }
            }
        }
    }

    /// `IsDone()` (`cxx:158-161`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Solution()` (`cxx:163-171`): `None` is `StdFail_NotDone`.
    pub fn solution(&self) -> Option<(f64, f64)> {
        if !self.done {
            return None;
        }
        Some((self.solution.x(), self.solution.y()))
    }
}
