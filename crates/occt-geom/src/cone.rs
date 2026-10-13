use std::sync::Arc;
use occt_core::gp::{GpCone, GpPnt, GpVec, GpTrsf};
use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::line::GeomLine;
use crate::surface::Surface;
use occt_core::elib::{slib, surface_eval, surface_eval_d3};

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
    /// `Geom_ConicalSurface::D2` → `ElSLib::ConeD2` (`ElSLib.cxx:867-...`):
    /// `Vuu = -R·Vxy`, `Vvv = 0`, `Vuv = sinA·DVxy`.
    /// Replaces the trait's central-difference default (audit A15).
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        surface_eval::cone_d2(&self.pos, u, v)
    }
    /// `Geom_ConicalSurface::EvalDN` (`Geom_ConicalSurface.cxx:319-333`).
    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        if nu + nv < 1 || nu < 0 || nv < 0 {
            return Err("Geom_UndefinedDerivative".into());
        }
        if nv > 1 {
            return Ok(GpVec::zero());
        }
        let pos = self.pos.position();
        Ok(surface_eval::cone_dn(
            u,
            v,
            &pos,
            self.pos.radius(),
            self.pos.semi_angle(),
            nu,
            nv,
        ))
    }
    /// `Geom_ConicalSurface::EvalD3` (`Geom_ConicalSurface.cxx:296-317`).
    fn try_eval_d3_uv(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        let d3 = surface_eval_d3::cone_d3(&self.pos, u, v);
        match (nu, nv) {
            (2, 1) => Ok(d3.d3_uuv),
            (1, 2) => Ok(d3.d3_uvv),
            _ => Ok(GpVec::zero()),
        }
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn is_u_periodic(&self) -> bool { true }
    fn cone_ref(&self) -> Option<(f64, f64)> {
        Some((self.pos.radius, self.pos.semi_angle))
    }
    fn continuity(&self) -> u8 { 6 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
    fn gp_cone(&self) -> Option<GpCone> { Some(self.pos.clone()) }

    /// `Geom_ConicalSurface::UReversed` (`Geom_ConicalSurface.cxx`): a copy with
    /// `UReverse()` applied (`gp_Cone.hxx:115` -> `pos.YReverse()`).
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        let mut pos = self.pos.clone();
        pos.u_reverse();
        Some(Arc::new(GeomCone::new(pos)))
    }

    /// `Geom_ConicalSurface::VReversed` (`Geom_Surface.cxx:42-48`) = `Copy()` +
    /// `VReverse` (`gp_Cone.hxx:119` -> `ZReverse`).
    fn v_reversed(&self) -> Option<Arc<dyn Surface>> {
        let mut pos = self.pos.clone();
        pos.v_reverse();
        Some(Arc::new(GeomCone::new(pos)))
    }

    /// `Geom_ConicalSurface::UReversedParameter`
    /// (`Geom_ConicalSurface.cxx:90-93`): `2*pi - U`.
    fn u_reversed_parameter(&self, u: f64) -> f64 { 2.0 * std::f64::consts::PI - u }

    /// `Geom_ConicalSurface::VReversedParameter`
    /// (`Geom_ConicalSurface.cxx:97-100`): `-V`.
    fn v_reversed_parameter(&self, v: f64) -> f64 { -v }

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
