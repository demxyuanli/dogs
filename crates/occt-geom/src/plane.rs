use std::sync::Arc;
use occt_core::gp::{GpLin, GpPln, GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;
use crate::line::GeomLine;
use crate::surface::Surface;
use occt_core::elib::{slib, surface_eval, surface_eval_d3};

#[derive(Debug, Clone)]
pub struct GeomPlane { pos: GpPln }

impl GeomPlane {
    pub fn new(pl: GpPln) -> Self { Self { pos: pl } }
    pub fn pln(&self) -> &GpPln { &self.pos }
}

impl Surface for GeomPlane {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::plane_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) { slib::plane_d1(&self.pos, u, v) }
    /// `Geom_Plane::D2` → `ElSLib::PlaneD2` (all second derivatives vanish).
    /// Replaces the trait's central-difference default (audit A15).
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        surface_eval::plane_d2(&self.pos, u, v)
    }
    /// `Geom_Plane::EvalDN` (`Geom_Plane.cxx:242-257`).
    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        if nu + nv < 1 || nu < 0 || nv < 0 {
            return Err("Geom_UndefinedDerivative".into());
        }
        let pos = self.pos.position();
        Ok(surface_eval::plane_dn(u, v, &pos, nu, nv))
    }
    /// `Geom_Plane::EvalD3` (`Geom_Plane.cxx:226-238`).
    /// ApproxSurface order 3 stores `D3UUV` / `D3UVV` (`GeomConvert_ApproxSurface.cxx:193-205`).
    fn try_eval_d3_uv(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        let d3 = surface_eval_d3::plane_d3(&self.pos, u, v);
        match (nu, nv) {
            (2, 1) => Ok(d3.d3_uuv),
            (1, 2) => Ok(d3.d3_uvv),
            _ => Ok(GpVec::zero()),
        }
    }
    fn u_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn v_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn continuity(&self) -> u8 { 6 }
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
    fn gp_pln(&self) -> Option<GpPln> { Some(self.pos.clone()) }

    /// `Geom_Plane::UReversed` (`Geom_Plane.cxx`): a copy with `UReverse()`
    /// applied (`gp_Pln.hxx:101` -> `myPosition.XReverse()`), which flips the
    /// plane's coefficients (`gp_Pln::Coefficients` negates them on an indirect
    /// placement).
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        let mut pos = self.pos.clone();
        pos.u_reverse();
        Some(Arc::new(GeomPlane::new(pos)))
    }

    /// `Geom_Plane::UIso` (`Geom_Plane.cxx:261-265`):
    /// `Geom_Line(ElSLib::PlaneUIso(pos, U))`.
    /// `ElSLib::PlaneUIso` (`ElSLib.cxx:1705-1712`) is a line through the
    /// placement point translated by `U * XDirection`, directed along
    /// `YDirection`.
    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        let mut l = GpLin::from_pnt_dir(self.pos.location(), *self.pos.y_axis().direction());
        let ve = GpVec::from_xyz(self.pos.x_axis().direction().xyz()).multiplied_scalar(u);
        l.translate_vec(&ve);
        Some(Arc::new(GeomLine::new(l)))
    }

    /// `Geom_Plane::VIso` (`Geom_Plane.cxx:269-273`):
    /// `Geom_Line(ElSLib::PlaneVIso(pos, V))`.
    /// `ElSLib::PlaneVIso` (`ElSLib.cxx:1770-1777`) is a line through the
    /// placement point translated by `V * YDirection`, directed along
    /// `XDirection`.
    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        let mut l = GpLin::from_pnt_dir(self.pos.location(), *self.pos.x_axis().direction());
        let ve = GpVec::from_xyz(self.pos.y_axis().direction().xyz()).multiplied_scalar(v);
        l.translate_vec(&ve);
        Some(Arc::new(GeomLine::new(l)))
    }
}
