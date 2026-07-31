//! 2D Bezier curve.

use crate::curve::Curve2d;
use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};

/// Rational-free polynomial Bezier curve in the plane.
#[derive(Clone)]
pub struct Geom2dBezierCurve {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
}

fn de_casteljau(xs: &[f64], ys: &[f64], u: f64) -> (f64, f64) {
    let mut x = xs.to_vec();
    let mut y = ys.to_vec();
    while x.len() > 1 {
        for i in 0..(x.len() - 1) {
            x[i] = x[i] + (x[i + 1] - x[i]) * u;
            y[i] = y[i] + (y[i + 1] - y[i]) * u;
        }
        x.pop();
        y.pop();
    }
    (x[0], y[0])
}

impl Geom2dBezierCurve {
    pub fn new(xs: Vec<f64>, ys: Vec<f64>) -> Result<Self, &'static str> {
        if xs.len() != ys.len() {
            return Err("Geom2dBezierCurve: xs and ys must have equal lengths");
        }
        if xs.len() < 2 {
            return Err("Geom2dBezierCurve requires at least 2 poles");
        }
        Ok(Self { xs, ys })
    }

    pub fn nb_poles(&self) -> usize {
        self.xs.len()
    }

    pub fn degree(&self) -> usize {
        self.xs.len() - 1
    }

    pub fn pole(&self, i: usize) -> GpPnt2d {
        GpPnt2d::new(self.xs[i], self.ys[i])
    }
}

impl Curve2d for Geom2dBezierCurve {
    fn d0(&self, u: f64) -> GpPnt2d {
        let (x, y) = de_casteljau(&self.xs, &self.ys, u);
        GpPnt2d::new(x, y)
    }

    fn d1(&self, u: f64) -> GpVec2d {
        let deg = self.degree();
        if deg == 0 {
            return GpVec2d::new(0.0, 0.0);
        }
        let d = deg as f64;
        let mut dx = Vec::with_capacity(deg);
        let mut dy = Vec::with_capacity(deg);
        for i in 0..deg {
            dx.push((self.xs[i + 1] - self.xs[i]) * d);
            dy.push((self.ys[i + 1] - self.ys[i]) * d);
        }
        let (x, y) = de_casteljau(&dx, &dy, u);
        GpVec2d::new(x, y)
    }

    fn d2(&self, u: f64) -> GpVec2d {
        let deg = self.degree();
        if deg < 2 {
            return GpVec2d::new(0.0, 0.0);
        }
        let d = deg as f64;
        let mut dx = Vec::with_capacity(deg);
        let mut dy = Vec::with_capacity(deg);
        for i in 0..deg {
            dx.push((self.xs[i + 1] - self.xs[i]) * d);
            dy.push((self.ys[i + 1] - self.ys[i]) * d);
        }
        let mut ddx = Vec::with_capacity(deg - 1);
        let mut ddy = Vec::with_capacity(deg - 1);
        for i in 0..(deg - 1) {
            ddx.push((dx[i + 1] - dx[i]) * (d - 1.0));
            ddy.push((dy[i + 1] - dy[i]) * (d - 1.0));
        }
        let (x, y) = de_casteljau(&ddx, &ddy, u);
        GpVec2d::new(x, y)
    }

    fn first_parameter(&self) -> f64 {
        0.0
    }

    fn last_parameter(&self) -> f64 {
        1.0
    }

    fn continuity(&self) -> usize {
        3
    }

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
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_bezier() {
        let c = Geom2dBezierCurve::new(vec![0.0, 1.0, 2.0], vec![0.0, 2.0, 0.0]).unwrap();
        let p = c.d0(0.5);
        assert!((p.x() - 1.0).abs() < 1e-12);
        assert!((p.y() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_mismatched_lengths() {
        assert!(Geom2dBezierCurve::new(vec![0.0, 1.0], vec![0.0]).is_err());
    }
}
