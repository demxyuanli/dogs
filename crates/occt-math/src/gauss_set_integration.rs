//! `math_GaussSetIntegration` (`math_GaussSetIntegration.cxx:29-99`) and the
//! `math_FunctionSet` interface it integrates
//! (`math_FunctionSet.hxx:30-62`).
//!
//! Translation notes:
//! - `math::GaussPoints(Ordre, GaussP)` / `math::GaussWeights(Ordre, GaussW)`
//!   (`math.cxx:25-31`, `:1300-1360`) hand out the tabulated OCCT nodes in
//!   decreasing order, positive half first; the port uses
//!   [`gauss_legendre`] on `[-1, 1]`, whose nodes ascend. The two sample sets
//!   of `cxx:79-95` (`Xm + Xr*P(i)` and `Xm - Xr*P(i)`) are the same, and the
//!   weights are symmetric, so the sums coincide.
//! - `Standard_NotImplemented_Raise_if` (`cxx:43-45`) is an assertion here: the
//!   single-variable / tabulated-order contract is a caller obligation.

use crate::gauss::gauss_legendre;
use crate::MathVector;

/// `math::GaussPointsMax()` (`math.cxx:25-28`): the largest tabulated order.
pub const GAUSS_POINTS_MAX: usize = 61;

/// `math_FunctionSet` (`math_FunctionSet.hxx:30-62`): the vector-valued
/// function `math_GaussSetIntegration` integrates.
pub trait MathFunctionSet {
    /// `NbVariables()` (`math_FunctionSet.hxx:34`).
    fn nb_variables(&self) -> usize;

    /// `NbEquations()` (`math_FunctionSet.hxx:38`).
    fn nb_equations(&self) -> usize;

    /// `Value(X, F)` (`math_FunctionSet.hxx:44`): `false` aborts the
    /// integration exactly as OCCT's early `return`s do (`cxx:68-94`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool;

    /// `GetStateNumber()` (`math_FunctionSet.hxx:59`).
    fn state_number(&mut self) -> i32 {
        0
    }
}

/// `math_GaussSetIntegration`
/// (`math_GaussSetIntegration.hxx:22-60`, `.cxx:29-99`).
pub struct GaussSetIntegration {
    val: MathVector,
    done: bool,
}

impl GaussSetIntegration {
    /// `math_GaussSetIntegration(F, Lower, Upper, Order)`
    /// (`math_GaussSetIntegration.cxx:29-99`).
    pub fn new<F: MathFunctionSet + ?Sized>(
        f: &mut F,
        lower: &MathVector,
        upper: &MathVector,
        order: usize,
    ) -> Self {
        let nb_equa = f.nb_equations();
        let nb_var = f.nb_variables();
        let mut val = MathVector::new(1, nb_equa);
        let mut f_val1 = MathVector::new(1, nb_equa);
        let mut f_val2 = MathVector::new(1, nb_equa);
        let mut t_val = MathVector::new(1, nb_var);

        // `Standard_NotImplemented_Raise_if` (`cxx:43-45`).
        assert!(
            nb_var == 1 && order <= GAUSS_POINTS_MAX,
            "GaussSetIntegration"
        );

        let mut done = false;

        let x_deb = lower.value(lower.lower());
        let x_fin = upper.value(upper.lower());
        let ordre = order;
        let xm = 0.5 * (x_deb + x_fin);
        let xr = 0.5 * (x_fin - x_deb);

        // `math::GaussPoints` / `math::GaussWeights` (`cxx:56-57`).
        let (gauss_p, gauss_w) = gauss_legendre(-1.0, 1.0, ordre);

        let ind = ordre / 2;
        let ind1 = (ordre + 1) / 2;
        if ind1 > ind {
            // odd case: the middle node is exactly `Xm` (`cxx:65-74`).
            t_val.set_value(1, xm);
            if !f.value(&t_val, &mut val) {
                return Self { val, done };
            }
            val.multiply_scalar(gauss_w[ind1 - 1]);
        } else {
            val.init(0.0);
        }

        for i in 1..=ind {
            t_val.set_value(1, xm + xr * gauss_p[i - 1]);
            if !f.value(&t_val, &mut f_val1) {
                return Self { val, done };
            }
            t_val.set_value(1, xm - xr * gauss_p[i - 1]);
            if !f.value(&t_val, &mut f_val2) {
                return Self { val, done };
            }
            f_val1.add(&f_val2);
            f_val1.multiply_scalar(gauss_w[i - 1]);
            val.add(&f_val1);
        }
        val.multiply_scalar(xr);

        done = true;
        Self { val, done }
    }

    /// `IsDone()` (`math_GaussSetIntegration.hxx:33`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Value()` (`math_GaussSetIntegration.hxx:38`). OCCT raises
    /// `StdFail_NotDone` when `!IsDone()`; the vector is returned as-is here.
    pub fn value(&self) -> &MathVector {
        &self.val
    }
}
