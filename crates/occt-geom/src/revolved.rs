//! Surface of revolution. Source: `Geom_SurfaceOfRevolution.hxx`
use std::sync::Arc;
use occt_core::gp::{GpPnt, GpVec, GpTrsf, GpAx1};
use crate::surface::Surface;
use crate::curve::Curve;

/// Surface obtained by revolving a curve around an axis.
#[derive(Clone)]
pub struct GeomRevolvedSurface {
    basis: Arc<dyn Curve>,
    axis: GpAx1,
}

impl GeomRevolvedSurface {
    pub fn new(curve: Arc<dyn Curve>, axis: GpAx1) -> Self { Self { basis: curve, axis } }
    pub fn basis_curve(&self) -> &Arc<dyn Curve> { &self.basis }
    pub fn axis(&self) -> &GpAx1 { &self.axis }

    fn transform_point(&self, p: &GpPnt, angle: f64) -> GpPnt {
        // Project p onto axis, rotate around axis
        let ax_loc = self.axis.location();
        let ax_dir = self.axis.direction();
        let v = p.coord.subtracted(&ax_loc.coord);
        let proj_len = v.dot(ax_dir.xyz());
        let proj = ax_dir.xyz().multiplied(proj_len);
        let rad = v.subtracted(&proj);
        let cos_a = angle.cos(); let sin_a = angle.sin();
        let rot = rad.multiplied(cos_a).added(&ax_dir.xyz().crossed(&rad).multiplied(sin_a));
        GpPnt::from_xyz(&ax_loc.coord.added(&proj).added(&rot))
    }
}

impl Surface for GeomRevolvedSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        let pt = self.basis.d0(u);
        self.transform_point(&pt, v)
    }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let pt = self.basis.d0(u);
        let (_, d1_u) = self.basis.d1(u);
        let p = self.transform_point(&pt, v);
        let du = GpVec::new(d1_u.x(), d1_u.y(), d1_u.z());
        let dv = GpVec::new(pt.x(), pt.y(), pt.z()); // simplified derivative
        (p, du, dv)
    }
    fn u_range(&self) -> (f64, f64) { (self.basis.first_parameter(), self.basis.last_parameter()) }
    fn v_range(&self) -> (f64, f64) { (0.0, 2.0 * std::f64::consts::PI) }
    fn is_v_periodic(&self) -> bool { true }
    fn is_surface_of_revolution(&self) -> bool { true }
    fn revolution_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        Some(self.basis.clone())
    }
    fn continuity(&self) -> u8 { self.basis.continuity() }
    fn transform(&mut self, _t: &GpTrsf) {}
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}
