//! `CSLib_NormalPolyDef` — the trigonometric polynomial whose roots decide
//! whether a singular point has a single normal.
//!
//! Source: `FoundationClasses/TKMath/CSLib/CSLib_NormalPolyDef.cxx:24-111`
//! (OCCT 8.0.0), declared `CSLib_NormalPolyDef.hxx:27-46`.
//!
//! The polynomial is
//! `F(X) = sum_{i=0}^{k0} Bin(k0, i) cos^i(X) sin^(k0-i)(X) li(i)` (`cxx:49-54`)
//! and its derivative is evaluated with the closed form of `cxx:75-81`.
//!
//! Translation notes:
//! - `theF` / `theD` are out-parameters in OCCT and become `&mut f64` here,
//! - the three members always answer `true`, as in the source; the negative
//!   exponents of `pow` at `i == k0` (resp. `i == 0`) are reproduced as-is,
//!   guarded by the same `RealSmall()` test (`cxx:44-47`, `:66-70`, `:94-98`).

use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::precision::REAL_SMALL;

/// `CSLib_NormalPolyDef` (`CSLib_NormalPolyDef.hxx:27-46`).
#[derive(Debug, Clone)]
pub struct NormalPolyDef {
    /// `myK0`
    k0: i32,
    /// `myTABli` — `NCollection_Array1<double>` indexed `0..=myK0`.
    tabli: Vec<f64>,
}

impl NormalPolyDef {
    /// `CSLib_NormalPolyDef(theK0, theLi)` (`cxx:24-32`).
    pub fn new(k0: i32, li: &[f64]) -> Self {
        let mut tabli = vec![0.0f64; (k0 + 1) as usize];
        for i in 0..=k0 {
            tabli[i as usize] = li[i as usize];
        }
        Self { k0, tabli }
    }
}

impl MathFunction for NormalPolyDef {
    /// `CSLib_NormalPolyDef::Value(theX, theF)` (`cxx:36-56`).
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = 0.0;
        let a_cos = x.cos();
        let a_sin = x.sin();

        if a_cos.abs() <= REAL_SMALL || a_sin.abs() <= REAL_SMALL {
            return true;
        }

        for i in 0..=self.k0 {
            *f += occt_core::bspl::plib::bin(self.k0, i)
                * a_cos.powi(i)
                * a_sin.powi(self.k0 - i)
                * self.tabli[i as usize];
        }
        true
    }
}

impl MathFunctionWithDerivative for NormalPolyDef {
    /// `CSLib_NormalPolyDef::Derivative(theX, theD)` (`cxx:60-82`).
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        *d = 0.0;
        let a_cos = x.cos();
        let a_sin = x.sin();

        if a_cos.abs() <= REAL_SMALL || a_sin.abs() <= REAL_SMALL {
            return true;
        }

        for i in 0..=self.k0 {
            *d += occt_core::bspl::plib::bin(self.k0, i)
                * a_cos.powi(i - 1)
                * a_sin.powi(self.k0 - i - 1)
                * (f64::from(self.k0) * a_cos * a_cos - f64::from(i))
                * self.tabli[i as usize];
        }
        true
    }

    /// `CSLib_NormalPolyDef::Values(theX, theF, theD)` (`cxx:86-111`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        *f = 0.0;
        *d = 0.0;
        let a_cos = x.cos();
        let a_sin = x.sin();

        if a_cos.abs() <= REAL_SMALL || a_sin.abs() <= REAL_SMALL {
            return true;
        }

        for i in 0..=self.k0 {
            let a_bin_coeff = occt_core::bspl::plib::bin(self.k0, i);
            let a_li_coeff = self.tabli[i as usize];

            *f += a_bin_coeff * a_cos.powi(i) * a_sin.powi(self.k0 - i) * a_li_coeff;
            *d += a_bin_coeff
                * a_cos.powi(i - 1)
                * a_sin.powi(self.k0 - i - 1)
                * (f64::from(self.k0) * a_cos * a_cos - f64::from(i))
                * a_li_coeff;
        }
        true
    }
}
