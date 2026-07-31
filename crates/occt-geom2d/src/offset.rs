//! Offset 2D curve. Source: `Geom2d_OffsetCurve.hxx`
use std::sync::Arc;
use occt_core::gp::{GpPnt2d, GpVec2d, GpTrsf2d};
use crate::curve::Curve2d;

#[derive(Debug, Clone)]
pub struct Geom2dOffsetCurve {
    basis: Arc<dyn Curve2d>,
    offset: f64,
}

impl Geom2dOffsetCurve {
    pub fn new(curve: Arc<dyn Curve2d>, offset: f64) -> Self { Self { basis: curve, offset } }
    pub fn basis_curve(&self) -> &Arc<dyn Curve2d> { &self.basis }
    pub fn offset_value(&self) -> f64 { self.offset }
}

impl Curve2d for Geom2dOffsetCurve {
    fn d0(&self, u: f64) -> GpPnt2d {
        let p = self.basis.d0(u);
        let (_, d) = self.basis.d1(u);
        let len = (d.x()*d.x() + d.y()*d.y()).sqrt();
        if len < 1e-30 { return p; }
        let nx = -d.y() / len; let ny = d.x() / len;
        GpPnt2d::new(p.x() + self.offset * nx, p.y() + self.offset * ny)
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (p, d) = self.basis.d1(u);
        let len = (d.x()*d.x() + d.y()*d.y()).sqrt();
        if len < 1e-30 { return (p, d); }
        let nx = -d.y() / len; let ny = d.x() / len;
        let pt = GpPnt2d::new(p.x() + self.offset * nx, p.y() + self.offset * ny);
        (pt, d) // first derivative same as basis (offset doesn't change tangent)
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1, d2) = self.basis.d2(u);
        let len = (d1.x()*d1.x() + d1.y()*d1.y()).sqrt();
        if len < 1e-30 { return (p, d1, d2); }
        let nx = -d1.y() / len; let ny = d1.x() / len;
        let pt = GpPnt2d::new(p.x() + self.offset * nx, p.y() + self.offset * ny);
        (pt, d1, d2)
    }
    fn first_parameter(&self) -> f64 { self.basis.first_parameter() }
    fn last_parameter(&self) -> f64 { self.basis.last_parameter() }
    fn is_periodic(&self) -> bool { self.basis.is_periodic() }
    fn period(&self) -> f64 { self.basis.period() }
    fn continuity(&self) -> u8 { self.basis.continuity().min(1) } // offset reduces to C1 at best
    fn transform(&mut self, t: &GpTrsf2d) { self.offset *= t.scale_factor().abs(); }
    fn reverse(&mut self) { self.offset = -self.offset; }
    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
}
