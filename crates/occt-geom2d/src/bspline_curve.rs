//! 2D B-spline curve. Source: `Geom2d_BSplineCurve.hxx`

use crate::curve::Curve2d;
use occt_core::gp::{GpPnt2d, GpVec2d, GpTrsf2d};
use occt_core::bspl::knots;

/// 2D B-spline curve (non-rational), stored as separate x/y pole arrays.
#[derive(Clone)]
pub struct Geom2dBSplineCurve {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub knots: Vec<f64>,
    pub degree: usize,
}

impl Geom2dBSplineCurve {
    /// Build a 2D B-spline. xs/ys lengths must match and knot count must be
    /// poles + degree + 1.
    pub fn new(xs: Vec<f64>, ys: Vec<f64>, knots: Vec<f64>, degree: usize) -> Result<Self, &'static str> {
        if xs.len() != ys.len() {
            return Err("Geom2dBSplineCurve: xs/ys length mismatch");
        }
        knots::check_degree(xs.len(), degree, knots.len())?;
        Ok(Self { xs, ys, knots, degree })
    }

    pub fn nb_poles(&self) -> usize { self.xs.len() }
    pub fn degree(&self) -> usize { self.degree }
    pub fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    pub fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }

    /// De Boor triangular evaluation on (x, y) pole pairs.
    fn de_boor(&self, u: f64) -> GpPnt2d {
        let n = self.xs.len();
        if n == 0 { return GpPnt2d::zero(); }
        let idx = knots::hunt(&self.knots, u).max(self.degree).min(n - 1);
        let mut x = vec![0.0f64; self.degree + 1];
        let mut y = vec![0.0f64; self.degree + 1];
        for k in 0..=self.degree {
            let pi = (idx - self.degree + k).min(n - 1);
            x[k] = self.xs[pi];
            y[k] = self.ys[pi];
        }
        for r in 1..=self.degree {
            for i in (r..=self.degree).rev() {
                let k0 = idx + i - self.degree;
                let k1 = k0 + self.degree + 1 - r;
                let alpha = (u - self.knots[k0]) / (self.knots[k1] - self.knots[k0]);
                if alpha.is_finite() {
                    x[i] = (1.0 - alpha) * x[i - 1] + alpha * x[i];
                    y[i] = (1.0 - alpha) * y[i - 1] + alpha * y[i];
                }
            }
        }
        GpPnt2d::new(x[self.degree], y[self.degree])
    }
}

impl Curve2d for Geom2dBSplineCurve {
    fn d0(&self, u: f64) -> GpPnt2d { self.de_boor(u) }

    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let h = 1e-6;
        let p = self.d0(u);
        let p1 = self.d0(u + h);
        let p2 = self.d0(u - h);
        (p, GpVec2d::new((p1.x() - p2.x()) / (2.0 * h), (p1.y() - p2.y()) / (2.0 * h)))
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1) = self.d1(u);
        let h = 1e-6;
        let p1 = self.d0(u + h);
        let p0 = self.d0(u);
        let p2 = self.d0(u - h);
        let d2 = GpVec2d::new(
            (p1.x() - 2.0 * p0.x() + p2.x()) / (h * h),
            (p1.y() - 2.0 * p0.y() + p2.y()) / (h * h),
        );
        (p, d1, d2)
    }

    fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }
    fn continuity(&self) -> u8 { if self.degree >= 2 { 3 } else { 1 } }

    fn transform(&mut self, t: &GpTrsf2d) {
        for i in 0..self.xs.len() {
            let mut p = GpPnt2d::new(self.xs[i], self.ys[i]);
            p.transform(t);
            self.xs[i] = p.x();
            self.ys[i] = p.y();
        }
    }

    fn reverse(&mut self) {
        self.xs.reverse();
        self.ys.reverse();
        let n = self.knots.len();
        for i in 0..n / 2 { self.knots.swap(i, n - 1 - i); }
        for k in self.knots.iter_mut() { *k = 1.0 - *k; }
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_midpoint() {
        let c = Geom2dBSplineCurve::new(
            vec![0.0, 1.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0, 1.0, 1.0],
            1,
        ).unwrap();
        let p = c.d0(0.5);
        assert!((p.x() - 0.5).abs() < 1e-12);
        assert!((p.y() - 0.0).abs() < 1e-12);
    }
}
