//! 2D B-spline curve. Source: `Geom2d_BSplineCurve.hxx`

use crate::curve::Curve2d;
use occt_core::bspl::{eval, knots};
use occt_core::gp::{GpPnt, GpPnt2d, GpTrsf2d, GpVec2d};

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

    /// Lift `(x, y)` poles to `z = 0` for `BSplCLib` evaluators.
    fn poles_3d(&self) -> Vec<GpPnt> {
        self.xs
            .iter()
            .zip(self.ys.iter())
            .map(|(&x, &y)| GpPnt::new(x, y, 0.0))
            .collect()
    }

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
        // `Geom2d_BSplineCurve::D1` / `BSplCLib::D1`.
        let poles = self.poles_3d();
        let (p, d) = eval::eval_curve_d1(&poles, &self.knots, self.degree, u);
        (GpPnt2d::new(p.x(), p.y()), GpVec2d::new(d.x(), d.y()))
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        // `Geom2d_BSplineCurve::D2` / `BSplCLib::D2`.
        let poles = self.poles_3d();
        let (p, d1, d2) = eval::eval_curve_d2(&poles, &self.knots, self.degree, u);
        (
            GpPnt2d::new(p.x(), p.y()),
            GpVec2d::new(d1.x(), d1.y()),
            GpVec2d::new(d2.x(), d2.y()),
        )
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
        // `Geom2d_BSplineCurve::Reverse` / `BSplCLib::Reverse` on a flat
        // knot sequence: `k' = umax - k` after reversing the array.
        // `1 - k` is only valid when the last knot is 1 (Shape-2 pcurves
        // use V knots on `[0, 150]`).
        self.xs.reverse();
        self.ys.reverse();
        let n = self.knots.len();
        if n == 0 {
            return;
        }
        let umax = self.knots[n - 1];
        for i in 0..n / 2 {
            self.knots.swap(i, n - 1 - i);
        }
        for k in self.knots.iter_mut() {
            *k = umax - *k;
        }
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
