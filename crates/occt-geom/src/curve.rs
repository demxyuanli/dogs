//! Abstract 3D parametric curve. Source: `Geom_Curve.hxx`
use occt_core::gp::{GpPnt, GpVec, GpTrsf};

pub trait Curve: Send + Sync {
    fn d0(&self, u: f64) -> GpPnt;
    fn d1(&self, u: f64) -> (GpPnt, GpVec);
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec);
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (p, d1, d2) = self.d2(u);
        (p, d1, d2, GpVec::zero())
    }
    fn value(&self, u: f64) -> GpPnt { self.d0(u) }
    fn first_parameter(&self) -> f64;
    fn last_parameter(&self) -> f64;
    fn is_periodic(&self) -> bool { false }
    fn period(&self) -> f64 { 0.0 }
    fn continuity(&self) -> u8;
    fn transform(&mut self, t: &GpTrsf);
    fn reverse(&mut self);
    fn clone_dyn(&self) -> Box<dyn Curve>;
    fn transformed(&self, t: &GpTrsf) -> Box<dyn Curve> { let mut c = self.clone_dyn(); c.transform(t); c }
    fn reversed(&self) -> Box<dyn Curve> { let mut c = self.clone_dyn(); c.reverse(); c }
}
