//! `IntPatch_CurvIntSurf` - the curve-on-surface / surface intersection point
//! computed from the `IntPatch_CSFunction` residual.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_CurvIntSurf.hxx`,
//! `IntPatch_CurvIntSurf_0.cxx` (which only specialises the template) and the
//! template body `TKGeomAlgo/IntImp/IntImp_IntCS.gxx`.
//!
//! Ported from `IntImp_IntCS.gxx`:
//! * both constructors (`:27-92`), including the `MarginCoef` extension of the
//!   UV bounds (`:58-77`).
//! * `Perform(U, V, W, Rsnld, u0, u1, v0, v1, w0, w1)` (`:94-154`): the three
//!   Newton attempts that restart from `W`, `w0` and `w1`.
//! * `IsDone` / `IsEmpty` / `Point` / `ParameterOnSurface` / `ParameterOnCurve`
//!   / `Function` (`:156-200`).
//!
//! OCCT copies the `IntPatch_CSFunction` by value (`gxx:35`); the port borrows
//! the two surfaces and the curve for the lifetime of the solver, which is the
//! same data (`IntPatchCSFunction` owns only the last residual and midpoint).
//!
//! The `StdFail_NotDone` / `Standard_DomainError` raises of the accessors are
//! `Err` returns, per the crate convention.

use occt_core::gp::GpPnt;
use occt_core::precision::{Precision, SQUARE_CONFUSION};
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_math::function_set_root::MathFunctionSetRoot;
use occt_math::MathVector;

use crate::intpatch_cs_function::IntPatchCSFunction;
use crate::intpatch_h_curve2d_tool as curve2d_tool;

/// `IntPatch_CurvIntSurf` (`IntPatch_CurvIntSurf.hxx:26-102`).
pub struct IntPatchCurvIntSurf<'a> {
    /// `done` (`hxx:96`); false while the last `Perform` found no root.
    done: bool,
    /// `empty` (`hxx:97`); true until a root within tolerance is accepted.
    empty: bool,
    /// `myFunction` (`hxx:98`).
    my_function: IntPatchCSFunction<'a>,
    /// `w` (`hxx:99`): the parameter on the 2d curve.
    w: f64,
    /// `u` (`hxx:100`).
    u: f64,
    /// `v` (`hxx:101`).
    v: f64,
    /// `tol` (`hxx:102`): `TolTangency * TolTangency`, floored at
    /// `Precision::SquareConfusion()` (`gxx:38-43`, `:85-89`).
    tol: f64,
}

impl<'a> IntPatchCurvIntSurf<'a> {
    /// `IntPatch_CurvIntSurf(U, V, W, F, TolTangency, MarginCoef = 0.)`
    /// (`gxx:27-79`).
    pub fn new(
        u: f64,
        v: f64,
        w: f64,
        f: &IntPatchCSFunction<'a>,
        tol_tangency: f64,
        margin_coef: f64,
    ) -> Self {
        let mut this = Self::with_tolerance(f, tol_tangency);

        let s = this.my_function.auxillar_surface();
        let c = this.my_function.auxillar_curve();

        let w0 = curve2d_tool::first_parameter(c);
        let w1 = curve2d_tool::last_parameter(c);

        let mut u0 = s.u_range().0;
        let mut v0 = s.v_range().0;
        let mut u1 = s.u_range().1;
        let mut v1 = s.v_range().1;

        if margin_coef > 0.0 {
            if !Precision::is_infinite(u0) && !Precision::is_infinite(u1) {
                let mut marg = (u1 - u0) * margin_coef;
                if u0 > u1 {
                    marg = -marg;
                }
                u0 -= marg;
                u1 += marg;
            }
            if !Precision::is_infinite(v0) && !Precision::is_infinite(v1) {
                let mut marg = (v1 - v0) * margin_coef;
                if v0 > v1 {
                    marg = -marg;
                }
                v0 -= marg;
                v1 += marg;
            }
        }

        let mut rsnld = MathFunctionSetRoot::with_iterations(&this.my_function, 100);
        this.perform(u, v, w, &mut rsnld, u0, u1, v0, v1, w0, w1);
        this
    }

    /// `IntPatch_CurvIntSurf(F, TolTangency)` (`gxx:81-92`): no bounds, no
    /// `Perform`; `u` / `v` / `w` keep their OCCT uninitialised slots, which the
    /// port reads as `0.` until a `Perform` overwrites them.
    pub fn with_tolerance(f: &IntPatchCSFunction<'a>, tol_tangency: f64) -> Self {
        let mut tol = tol_tangency * tol_tangency;
        if tol < SQUARE_CONFUSION {
            tol = SQUARE_CONFUSION;
        }
        Self {
            done: true,
            empty: true,
            my_function: *f,
            w: 0.0,
            u: 0.0,
            v: 0.0,
            tol,
        }
    }

    /// `Perform(U, V, W, Rsnld, u0, u1, v0, v1, w0, w1)` (`gxx:94-154`).
    pub fn perform(
        &mut self,
        u: f64,
        v: f64,
        w: f64,
        rsnld: &mut MathFunctionSetRoot,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
        w0: f64,
        w1: f64,
    ) {
        self.done = true;
        let mut uvap = MathVector::new(1, 3);
        uvap.set_value(1, u);
        uvap.set_value(2, v);
        uvap.set_value(3, w);

        let mut born_inf = MathVector::new(1, 3);
        let mut born_sup = MathVector::new(1, 3);
        let mut tolerance = MathVector::new(1, 3);

        let s = self.my_function.auxillar_surface();
        let c = self.my_function.auxillar_curve();

        born_inf.set_value(1, u0);
        born_inf.set_value(2, v0);
        born_sup.set_value(1, u1);
        born_sup.set_value(2, v1);

        born_inf.set_value(3, w0);
        born_sup.set_value(3, w1);

        tolerance.set_value(1, u_resolution(s, Precision::CONFUSION));
        tolerance.set_value(2, v_resolution(s, Precision::CONFUSION));
        tolerance.set_value(3, curve2d_tool::resolution(c, Precision::CONFUSION));
        rsnld.set_tolerance(&tolerance);
        let mut autretentative = 0;
        self.done = false;
        loop {
            if autretentative == 1 {
                uvap.set_value(3, w0);
            } else if autretentative == 2 {
                uvap.set_value(3, w1);
            }
            autretentative += 1;
            rsnld.perform_with_bounds(&mut self.my_function, &uvap, &born_inf, &born_sup, false);
            if rsnld.is_done() {
                let abs_root = self.my_function.root().abs();
                if abs_root <= self.tol {
                    uvap = rsnld.root();
                    self.u = uvap.value(1);
                    self.v = uvap.value(2);
                    self.w = uvap.value(3);
                    self.empty = false;
                    self.done = true;
                }
            }
            if self.done || autretentative >= 3 {
                break;
            }
        }
    }

    /// `IsDone()` (`gxx:156-159`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsEmpty()` (`gxx:161-166`): `Err` is the `StdFail_NotDone` raise.
    pub fn is_empty(&self) -> Result<bool, String> {
        if !self.done {
            return Err("StdFail_NotDone".into());
        }
        Ok(self.empty)
    }

    /// `Point()` (`gxx:168-175`): `Err` is `StdFail_NotDone` /
    /// `Standard_DomainError`.
    pub fn point(&self) -> Result<&GpPnt, String> {
        if !self.done {
            return Err("StdFail_NotDone".into());
        }
        if self.empty {
            return Err("Standard_DomainError".into());
        }
        Ok(self.my_function.point())
    }

    /// `ParameterOnSurface(U, V)` (`gxx:177-186`).
    pub fn parameter_on_surface(&self) -> Result<(f64, f64), String> {
        if !self.done {
            return Err("StdFail_NotDone".into());
        }
        if self.empty {
            return Err("Standard_DomainError".into());
        }
        Ok((self.u, self.v))
    }

    /// `ParameterOnCurve()` (`gxx:188-195`).
    pub fn parameter_on_curve(&self) -> Result<f64, String> {
        if !self.done {
            return Err("StdFail_NotDone".into());
        }
        if self.empty {
            return Err("Standard_DomainError".into());
        }
        Ok(self.w)
    }

    /// `Function()` (`gxx:197-200`).
    pub fn function(&mut self) -> &mut IntPatchCSFunction<'a> {
        &mut self.my_function
    }
}
