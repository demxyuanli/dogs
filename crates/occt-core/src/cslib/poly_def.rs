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
        let cos_val = x.cos();
        let sin_val = x.sin();
        if cos_val.abs() <= REAL_SMALL || sin_val.abs() <= REAL_SMALL {
            return true;
        }
        for i in 0..=self.k0 {
            *f += plib::bin(self.k0, i)
                * cos_val.powi(i)
                * sin_val.powi(self.k0 - i)
                * self.li[i as usize];
        }
        true
    }
}

impl MathFunctionWithDerivative for NormalPolyDef {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        *d = 0.0;
        let cos_val = x.cos();
        let sin_val = x.sin();
        if cos_val.abs() <= REAL_SMALL || sin_val.abs() <= REAL_SMALL {
            return true;
        }
        for i in 0..=self.k0 {
            *d += plib::bin(self.k0, i)
                * cos_val.powi(i - 1)
                * sin_val.powi(self.k0 - i - 1)
                * (self.k0 as f64 * cos_val * cos_val - i as f64)
                * self.li[i as usize];
        }
        true
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        *f = 0.0;
        *d = 0.0;
        let cos_val = x.cos();
        let sin_val = x.sin();
        if cos_val.abs() <= REAL_SMALL || sin_val.abs() <= REAL_SMALL {
            return true;
        }
        for i in 0..=self.k0 {
            let bin = plib::bin(self.k0, i);
            let li = self.li[i as usize];
            *f += bin * cos_val.powi(i) * sin_val.powi(self.k0 - i) * li;
            *d += bin
                * cos_val.powi(i - 1)
                * sin_val.powi(self.k0 - i - 1)
                * (self.k0 as f64 * cos_val * cos_val - i as f64)
                * li;
        }
        true
    }
}
