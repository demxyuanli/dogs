use std::sync::Arc;
use occt_core::gp::{GpCone, GpPnt, GpVec, GpTrsf};
use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::line::GeomLine;
use crate::surface::Surface;
use occt_core::elib::{slib, surface_eval};

#[derive(Debug, Clone)]
pub struct GeomCone { pos: GpCone }

impl GeomCone {
    pub fn new(c: GpCone) -> Self { Self { pos: c } }
    /// `Geom_ConicalSurface::Cone`.
    pub fn cone(&self) -> &GpCone { &self.pos }
}

impl Surface for GeomCone {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::cone_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        surface_eval::cone_d1(&self.pos, u, v)
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn is_u_periodic(&self) -> bool { true }
    fn cone_ref(&self) -> Option<(f64, f64)> {
        Some((self.pos.radius, self.pos.semi_angle))
    }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
    fn gp_cone(&self) -> Option<GpCone> { Some(self.pos.clone()) }

    /// `Geom_ConicalSurface::UIso` (`Geom_ConicalSurface.cxx:337-341`):
    /// `Geom_Line(ElSLib::ConeUIso(pos, radius, semiAngle, U))`.
    /// `ElSLib::ConeUIso` (`ElSLib.cxx:1727-1733`) is the generatrix through
    /// `ConeValue(U, 0)` directed along `ConeD1(U, 0).DV`.
    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        Some(Arc::new(GeomLine::new(slib::cone_u_iso(
            &self.pos.pos,
            self.pos.radius,
            self.pos.semi_angle,
            u,
        ))))
    }

    /// `Geom_ConicalSurface::VIso` (`Geom_ConicalSurface.cxx:345-349`):
    /// `Geom_Circle(ElSLib::ConeVIso(pos, radius, semiAngle, V))`.
    /// `ElSLib::ConeVIso` (`ElSLib.cxx:1793-1810`) is the circle of radius
    /// `radius + V*sin(semiAngle)` centered `V*cos(semiAngle)` along the axis;
    /// a negative radius reverses both X and Y directions of the placement.
    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        Some(Arc::new(GeomCircle::new(slib::cone_v_iso(
            &self.pos.pos,
            self.pos.radius,
            self.pos.semi_angle,
            v,
        ))))
    }
}
