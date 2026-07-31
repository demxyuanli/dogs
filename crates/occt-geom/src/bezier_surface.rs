//! 3D Bezier surface.

use crate::surface::Surface;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

/// Tensor-product (non-rational) Bezier surface in 3D.
///
/// Poles are stored row-major: `poles[i * n_v + j]` is the control point for
/// the `i`-th u-index and `j`-th v-index.
#[derive(Clone)]
pub struct GeomBezierSurface {
    pub poles: Vec<GpPnt>,
    pub n_u: usize,
    pub n_v: usize,
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

impl GeomBezierSurface {
    pub fn new(poles: Vec<GpPnt>, n_u: usize, n_v: usize) -> Result<Self, &'static str> {
        if poles.len() != n_u * n_v {
            return Err("GeomBezierSurface: poles.len() must equal n_u * n_v");
        }
        if poles.len() < 4 {
            return Err("GeomBezierSurface requires at least 4 poles");
        }
        Ok(Self { poles, n_u, n_v })
    }

    pub fn nb_u_poles(&self) -> usize {
        self.n_u
    }

    pub fn nb_v_poles(&self) -> usize {
        self.n_v
    }

    pub fn u_degree(&self) -> usize {
        self.n_u - 1
    }

    pub fn v_degree(&self) -> usize {
        self.n_v - 1
    }

    fn du(&self, u: f64, v: f64) -> GpVec {
        if self.n_u < 2 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let nu = self.n_u as f64;
        // Differentiate the u-direction Bernstein basis: net is (n_u-1) x n_v.
        let mut net: Vec<GpVec> = Vec::with_capacity((self.n_u - 1) * self.n_v);
        for i in 0..(self.n_u - 1) {
            for j in 0..self.n_v {
                let p0 = &self.poles[i * self.n_v + j];
                let p1 = &self.poles[(i + 1) * self.n_v + j];
                net.push(GpVec::new(
                    (p1.x() - p0.x()) * nu,
                    (p1.y() - p0.y()) * nu,
                    (p1.z() - p0.z()) * nu,
                ));
            }
        }
        // Evaluate each row (along v) then along u.
        let mut pts = Vec::with_capacity(self.n_u - 1);
        for i in 0..(self.n_u - 1) {
            let start = i * self.n_v;
            pts.push(de_casteljau_vec(&net[start..start + self.n_v], v));
        }
        de_casteljau_vec(&pts, u)
    }

    fn dv(&self, u: f64, v: f64) -> GpVec {
        if self.n_v < 2 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let nv = self.n_v as f64;
        // Differentiate the v-direction Bernstein basis: net is n_u x (n_v-1).
        let mut net: Vec<GpVec> = Vec::with_capacity(self.n_u * (self.n_v - 1));
        for i in 0..self.n_u {
            for j in 0..(self.n_v - 1) {
                let p0 = &self.poles[i * self.n_v + j];
                let p1 = &self.poles[i * self.n_v + j + 1];
                net.push(GpVec::new(
                    (p1.x() - p0.x()) * nv,
                    (p1.y() - p0.y()) * nv,
                    (p1.z() - p0.z()) * nv,
                ));
            }
        }
        // Evaluate each row (along v) then along u.
        let mut pts = Vec::with_capacity(self.n_u);
        for i in 0..self.n_u {
            let start = i * (self.n_v - 1);
            pts.push(de_casteljau_vec(&net[start..start + self.n_v - 1], v));
        }
        de_casteljau_vec(&pts, u)
    }
}

impl Surface for GeomBezierSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        // Evaluate each u-row along v, then the resulting points along u.
        let mut rows = Vec::with_capacity(self.n_u);
        for i in 0..self.n_u {
            let start = i * self.n_v;
            rows.push(de_casteljau(&self.poles[start..start + self.n_v], v));
        }
        de_casteljau(&rows, u)
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        (self.d0(u, v), self.du(u, v), self.dv(u, v))
    }

    fn u_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    fn v_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    fn continuity(&self) -> usize {
        3
    }

    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.poles {
            *p = p.transformed(t);
        }
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
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
    fn bilinear_surface() {
        let s = GeomBezierSurface::new(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
                GpPnt::new(1.0, 1.0, 0.0),
            ],
            2,
            2,
        )
        .unwrap();

        assert!(approx(&s.d0(0.5, 0.5), &GpPnt::new(0.5, 0.5, 0.0)));
        assert!(approx(&s.d0(0.0, 0.0), &GpPnt::new(0.0, 0.0, 0.0)));
        assert!(approx(&s.d0(1.0, 1.0), &GpPnt::new(1.0, 1.0, 0.0)));
    }

    #[test]
    fn rejects_mismatched_pole_count() {
        assert!(GeomBezierSurface::new(vec![GpPnt::new(0.0, 0.0, 0.0)], 2, 2).is_err());
    }
}
