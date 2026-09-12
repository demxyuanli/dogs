use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpTorus, GpPnt, GpVec, GpTrsf};
use crate::surface::Surface;

#[derive(Debug, Clone)]
pub struct GeomTorus { pos: GpTorus }

impl GeomTorus {
    pub fn new(t: GpTorus) -> Self { Self { pos: t } }
    /// `Geom_ToroidalSurface::Torus`.
    pub fn torus(&self) -> &GpTorus { &self.pos }
}

impl Surface for GeomTorus {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::torus_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        surface_eval::torus_d1(&self.pos, u, v)
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn is_u_periodic(&self) -> bool { true }
    fn is_v_periodic(&self) -> bool { true }
    fn gp_torus(&self) -> Option<GpTorus> { Some(self.pos) }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}
