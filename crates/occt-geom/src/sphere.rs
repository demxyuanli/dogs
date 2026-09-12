use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpSphere, GpPnt, GpVec, GpTrsf};
use crate::surface::Surface;

#[derive(Debug, Clone)]
pub struct GeomSphere { pos: GpSphere }

impl GeomSphere {
    pub fn new(s: GpSphere) -> Self { Self { pos: s } }
    /// `Geom_SphericalSurface::Sphere`.
    pub fn sphere(&self) -> &GpSphere { &self.pos }
}

impl Surface for GeomSphere {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::sphere_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        surface_eval::sphere_d1(&self.pos, u, v)
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2) }
    fn is_u_periodic(&self) -> bool { true }
    fn gp_sphere(&self) -> Option<GpSphere> { Some(self.pos) }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}
