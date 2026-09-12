//! `CSLib_NormalPolyDef`.
//! Source: `CSLib_NormalPolyDef.cxx:24-112`.

use crate::bspl::plib;
use crate::math_fn::{MathFunction, MathFunctionWithDerivative};
use crate::precision::REAL_SMALL;

/// `CSLib_NormalPolyDef` (`CSLib_NormalPolyDef.hxx:39-82`).
pub struct NormalPolyDef {
    k0: i32,
    li: Vec<f64>,
}

impl NormalPolyDef {
    /// `CSLib_NormalPolyDef::CSLib_NormalPolyDef` (`cxx:24-32`).
    pub fn new(k0: i32, li: &[f64]) -> Self {
        let mut tab = vec![0.0; (k0 as usize) + 1];
        for i in 0..=k0 {
            tab[i as usize] = li[i as usize];
        }
        Self { k0, li: tab }
    }
}

impl MathFunction for NormalPolyDef {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = 0.0;
        let a_cos = x.cos();
        let a_sin = x.sin();
        if a_cos.abs() <= REAL_SMALL || a_sin.abs() <= REAL_SMALL {
            return true;
        }
        for i in 0..=self.k0 {
            *f += plib::bin(self.k0, i)
                * a_cos.powi(i)
                * a_sin.powi(self.k0 - i)
                * self.li[i as usize];
        }
        true
    }
}

impl MathFunctionWithDerivative for NormalPolyDef {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        *d = 0.0;
        let a_cos = x.cos();
        let a_sin = x.sin();
        if a_cos.abs() <= REAL_SMALL || a_sin.abs() <= REAL_SMALL {
            return true;
        }
        for i in 0..=self.k0 {
            *d += plib::bin(self.k0, i)
                * a_cos.powi(i - 1)
                * a_sin.powi(self.k0 - i - 1)
                * (self.k0 as f64 * a_cos * a_cos - i as f64)
                * self.li[i as usize];
        }
        true
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        *f = 0.0;
        *d = 0.0;
        let a_cos = x.cos();
        let a_sin = x.sin();
        if a_cos.abs() <= REAL_SMALL || a_sin.abs() <= REAL_SMALL {
            return true;
        }
        for i in 0..=self.k0 {
            let bin = plib::bin(self.k0, i);
            let li = self.li[i as usize];
            *f += bin * a_cos.powi(i) * a_sin.powi(self.k0 - i) * li;
            *d += bin
                * a_cos.powi(i - 1)
                * a_sin.powi(self.k0 - i - 1)
                * (self.k0 as f64 * a_cos * a_cos - i as f64)
                * li;
        }
        true
    }
}
