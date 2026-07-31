//! Abstract 3D parametric surface. Source: `Geom_Surface.hxx`
use occt_core::gp::{GpPnt, GpVec, GpTrsf};

pub trait Surface: Send + Sync {
    fn d0(&self, u: f64, v: f64) -> GpPnt;
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec);
    fn value(&self, u: f64, v: f64) -> GpPnt { self.d0(u, v) }
    fn u_range(&self) -> (f64, f64);
    fn v_range(&self) -> (f64, f64);
    fn is_u_periodic(&self) -> bool { false }
    fn is_v_periodic(&self) -> bool { false }
    fn continuity(&self) -> u8;
    fn transform(&mut self, t: &GpTrsf);
    fn clone_dyn(&self) -> Box<dyn Surface>;
    fn transformed(&self, t: &GpTrsf) -> Box<dyn Surface> { let mut s = self.clone_dyn(); s.transform(t); s }
}
