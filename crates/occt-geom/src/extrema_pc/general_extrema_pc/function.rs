//! `ExtPcSolution` and `ExtPcf` with its `MathFunction` implementations.

use super::*;

/// One `Extrema_ExtPC` solution (`SquareDistance` / `IsMin` / `Point`).
#[derive(Debug, Clone, Copy)]
pub struct ExtPcSolution {
    pub u: f64,
    pub point: GpPnt,
    pub sq_dist: f64,
    pub is_min: bool,
}

/// `Extrema_GFuncExtPC` for 3D point–curve (`Extrema_PCFOfEPCOfExtPC`).
pub(crate) struct ExtPcf<'a> {
    pub(crate) curve: &'a dyn Curve,
    pub(crate) point: GpPnt,
    pub(crate) u: f64,
    pub(crate) pc: GpPnt,
    pub(crate) d1f: f64,
    pub(crate) u_inf: f64,
    pub(crate) u_sup: f64,
    /// `myMaxDerivOrder` (`hxx:67` / Initialize GetType switch).
    pub(crate) max_deriv_order: i32,
    /// `myTol` — SearchOfTolerance or MinTol (`hxx:92-98`).
    pub(crate) tol: f64,
    pub(crate) sols: Vec<ExtPcSolution>,
}

impl MathFunction for ExtPcf<'_> {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        self.value_at(x, f)
    }

    fn state_number(&mut self) -> i32 {
        // `GFuncExtPC::GetStateNumber` (`hxx:357-378`).
        let sq = self.pc.square_distance(&self.point);
        let mut ff = 0.0;
        let mut dd = 0.0;
        let _ = self.values_at(self.u, &mut ff, &mut dd);
        let is_min = self.d1f > 0.0;
        self.sols.push(ExtPcSolution {
            u: self.u,
            point: self.pc,
            sq_dist: sq,
            is_min,
        });
        0
    }
}

impl MathFunctionWithDerivative for ExtPcf<'_> {
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let mut f = 0.0;
        self.values_at(x, &mut f, d)
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        self.values_at(x, f, d)
    }
}

