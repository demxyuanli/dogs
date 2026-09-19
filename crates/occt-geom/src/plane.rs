use std::sync::Arc;
use occt_core::gp::{GpLin, GpPln, GpPnt, GpVec, GpTrsf};
use crate::curve::Curve;
use crate::line::GeomLine;
use crate::surface::Surface;
use occt_core::elib::slib;

#[derive(Debug, Clone)]
pub struct GeomPlane { pos: GpPln }

impl GeomPlane {
    pub fn new(pl: GpPln) -> Self { Self { pos: pl } }
    pub fn pln(&self) -> &GpPln { &self.pos }
}

impl Surface for GeomPlane {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::plane_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) { slib::plane_d1(&self.pos, u, v) }
    fn u_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn v_range(&self) -> (f64, f64) { (f64::NEG_INFINITY, f64::INFINITY) }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
    fn gp_pln(&self) -> Option<GpPln> { Some(self.pos.clone()) }

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
