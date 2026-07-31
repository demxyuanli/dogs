//! Offset 3D curve. Source: `Geom_OffsetCurve.hxx`
use std::sync::Arc;
use occt_core::gp::{GpPnt, GpVec, GpDir, GpTrsf};
use crate::curve::Curve;

#[derive(Clone)]
pub struct GeomOffsetCurve { basis: Arc<dyn Curve>, offset: f64, direction: GpDir }

impl GeomOffsetCurve {
    pub fn new(curve: Arc<dyn Curve>, offset: f64, dir: GpDir) -> Self { Self { basis: curve, offset, direction: dir } }
}

impl Curve for GeomOffsetCurve {
    fn d0(&self, u: f64) -> GpPnt {
        let p = self.basis.d0(u);
        GpPnt::from_xyz(&p.coord.added(&self.direction.xyz().multiplied(self.offset)))
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) { let (p, d) = self.basis.d1(u); (self.d0(u), d) }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) { let (p, d1, d2) = self.basis.d2(u); (self.d0(u), d1, d2) }
    fn first_parameter(&self) -> f64 { self.basis.first_parameter() }
    fn last_parameter(&self) -> f64 { self.basis.last_parameter() }
    fn continuity(&self) -> u8 { self.basis.continuity() }
    fn transform(&mut self, t: &GpTrsf) { self.offset *= t.scale_factor().abs(); }
    fn reverse(&mut self) { self.offset = -self.offset; }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}
