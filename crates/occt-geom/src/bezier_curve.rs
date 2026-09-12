//! 3D Bezier curve.

use crate::curve::Curve;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

/// Rational-free polynomial Bezier curve in 3D.
#[derive(Clone)]
pub struct GeomBezierCurve {
    pub poles: Vec<GpPnt>,
}

fn lerp(a: &GpPnt, b: &GpPnt, t: f64) -> GpPnt {
    GpPnt::new(
        a.x() + (b.x() - a.x()) * t,
        a.y() + (b.y() - a.y()) * t,
        a.z() + (b.z() - a.z()) * t,
    )
}

fn lerp_vec(a: &GpVec, b: &GpVec, t: f64) -> GpVec {
    GpVec::new(
        a.x() + (b.x() - a.x()) * t,
        a.y() + (b.y() - a.y()) * t,
        a.z() + (b.z() - a.z()) * t,
    )
}

fn de_casteljau(poles: &[GpPnt], u: f64) -> GpPnt {
    let mut pts: Vec<GpPnt> = poles.to_vec();
    while pts.len() > 1 {
        for i in 0..(pts.len() - 1) {
            pts[i] = lerp(&pts[i], &pts[i + 1], u);
        }
        pts.pop();
    }
    pts[0]
}

fn de_casteljau_vec(vs: &[GpVec], u: f64) -> GpVec {
    let mut pts: Vec<GpVec> = vs.to_vec();
    while pts.len() > 1 {
        for i in 0..(pts.len() - 1) {
            pts[i] = lerp_vec(&pts[i], &pts[i + 1], u);
        }
        pts.pop();
    }
    pts[0]
}

impl GeomBezierCurve {
    pub fn new(poles: Vec<GpPnt>) -> Result<Self, &'static str> {
        if poles.len() < 2 {
            return Err("GeomBezierCurve requires at least 2 poles");
        }
        Ok(Self { poles })
    }

    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }

    pub fn degree(&self) -> usize {
        self.poles.len() - 1
    }

    pub fn pole(&self, i: usize) -> &GpPnt {
        &self.poles[i]
    }

    pub fn set_pole(&mut self, i: usize, p: GpPnt) {
        self.poles[i] = p;
    }

    fn tangent(&self, u: f64) -> GpVec {
        let deg = self.degree();
        if deg == 0 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let d = deg as f64;
        let mut dp: Vec<GpVec> = Vec::with_capacity(deg);
        for i in 0..deg {
            let p0 = &self.poles[i];
            let p1 = &self.poles[i + 1];
            dp.push(GpVec::new(
                (p1.x() - p0.x()) * d,
                (p1.y() - p0.y()) * d,
                (p1.z() - p0.z()) * d,
            ));
        }
        de_casteljau_vec(&dp, u)
    }

    fn second_deriv(&self, u: f64) -> GpVec {
        let deg = self.degree();
        if deg < 2 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let d = deg as f64;
        let mut dp: Vec<GpVec> = Vec::with_capacity(deg);
        for i in 0..deg {
            let p0 = &self.poles[i];
            let p1 = &self.poles[i + 1];
            dp.push(GpVec::new(
                (p1.x() - p0.x()) * d,
                (p1.y() - p0.y()) * d,
                (p1.z() - p0.z()) * d,
            ));
        }
        let mut ddp: Vec<GpVec> = Vec::with_capacity(deg - 1);
        for i in 0..(deg - 1) {
            ddp.push(GpVec::new(
                (dp[i + 1].x() - dp[i].x()) * (d - 1.0),
                (dp[i + 1].y() - dp[i].y()) * (d - 1.0),
                (dp[i + 1].z() - dp[i].z()) * (d - 1.0),
            ));
        }
        de_casteljau_vec(&ddp, u)
    }
}

impl Curve for GeomBezierCurve {
    fn d0(&self, u: f64) -> GpPnt {
        de_casteljau(&self.poles, u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        (self.d0(u), self.tangent(u))
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        (self.d0(u), self.tangent(u), self.second_deriv(u))
    }

    fn first_parameter(&self) -> f64 {
        0.0
    }

    fn last_parameter(&self) -> f64 {
        1.0
    }

    fn continuity(&self) -> u8 {
        3
    }

    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.poles {
            *p = p.transformed(t);
        }
    }

    fn reverse(&mut self) {
        self.poles.reverse();
    }

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }

    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        Some(&self.poles)
    }
    fn nurbs_degree(&self) -> Option<usize> {
        Some(self.poles.len().saturating_sub(1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: &GpPnt, b: &GpPnt) -> bool {
        (a.x() - b.x()).abs() < 1e-12
            && (a.y() - b.y()).abs() < 1e-12
            && (a.z() - b.z()).abs() < 1e-12
    }

    #[test]
    fn quadratic_bezier() {
        let c = GeomBezierCurve::new(vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 2.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
        ])
        .unwrap();

        assert!(approx(&c.d0(0.5), &GpPnt::new(1.0, 1.0, 0.0)));
        assert!(approx(&c.d0(0.0), &GpPnt::new(0.0, 0.0, 0.0)));
        assert!(approx(&c.d0(1.0), &GpPnt::new(2.0, 0.0, 0.0)));
    }

    #[test]
    fn requires_two_poles() {
        assert!(GeomBezierCurve::new(vec![GpPnt::new(0.0, 0.0, 0.0)]).is_err());
    }
}
