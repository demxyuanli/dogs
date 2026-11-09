//! Convert piecewise polynomials to B-spline poles.
//! Source: `Convert_CompPolynomialToPoles.cxx`.

use crate::bspl::banded_interp::{interpolate, knot_sequence, schoenberg_points};
use crate::bspl::plib_eval::eval_poly0;

/// Result of `Convert_CompPolynomialToPoles`.
pub struct CompPolynomialToPoles {
    pub poles: Vec<Vec<f64>>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
    pub done: bool,
}

impl CompPolynomialToPoles {
    /// Continuity-per-span constructor (`cxx:89-138`).
    pub fn new(
        num_curves: i32,
        dimension: i32,
        max_degree: i32,
        continuity: &[i32],
        num_coeff: &[i32],
        coefficients: &[f64],
        poly_intervals: &[[f64; 2]],
        true_intervals: &[f64],
    ) -> Result<Self, ()> {
        if num_curves <= 0 || max_degree <= 0 || dimension <= 0 {
            return Err(());
        }
        let mut degree = 0i32;
        for i in 0..num_curves as usize {
            degree = degree.max(num_coeff[i] - 1);
        }
        let n = num_curves as usize;
        let mut knots = vec![0.0; n + 1];
        for i in 0..=n {
            knots[i] = true_intervals[i];
        }
        let mut mults = vec![0i32; n + 1];
        for i in 1..n {
            if continuity[i] > degree && n > 1 {
                return Err(());
            }
            mults[i] = degree - continuity[i];
        }
        mults[0] = degree + 1;
        mults[n] = degree + 1;
        let mut s = Self {
            poles: Vec::new(),
            knots,
            mults,
            degree,
            done: false,
        };
        s.perform(
            num_curves,
            max_degree,
            dimension,
            num_coeff,
            coefficients,
            poly_intervals,
            true_intervals,
        )?;
        Ok(s)
    }

    fn perform(
        &mut self,
        num_curves: i32,
        max_degree: i32,
        dimension: i32,
        num_coeff: &[i32],
        coefficients: &[f64],
        poly_intervals: &[[f64; 2]],
        true_intervals: &[f64],
    ) -> Result<(), ()> {
        let dim = dimension as usize;
        let mut num_flat = 2 * self.degree + 2;
        for i in 1..self.mults.len() - 1 {
            num_flat += self.mults[i];
        }
        let num_poles = num_flat - self.degree - 1;
        let flat = knot_sequence(&self.knots, &self.mults, self.degree);
        let params = schoenberg_points(self.degree as usize, &flat, num_poles as usize);
        let mut poles = vec![0.0; (num_poles as usize) * dim];
        let mut index = 2i32;
        let mut tindex = 1usize;
        let mut pindex = 0usize;
        for ii in 0..num_poles as usize {
            while params[ii] >= true_intervals[tindex] && index <= num_curves {
                index += 1;
                tindex += 1;
                pindex += 1;
            }
            let mut nv = params[ii] - true_intervals[tindex - 1];
            nv /= true_intervals[tindex] - true_intervals[tindex - 1];
            nv = (1.0 - nv) * poly_intervals[pindex][0] + nv * poly_intervals[pindex][1];
            let coeff_index =
                ((index - 2) as usize) * dim * (max_degree.max(self.degree) as usize + 1);
            let deg = (num_coeff[(index - 2) as usize] - 1) as i32;
            let mut out = vec![0.0; dim];
            eval_poly0(&coefficients[coeff_index..], deg, dim, nv, &mut out);
            for d in 0..dim {
                poles[ii * dim + d] = out[d];
            }
        }
        interpolate(
            self.degree as usize,
            &flat,
            &params,
            &mut poles,
            dim,
        )
        .map_err(|_| ())?;
        self.poles = (0..num_poles as usize)
            .map(|i| poles[i * dim..(i + 1) * dim].to_vec())
            .collect();
        self.done = true;
        Ok(())
    }
}
