//! 2D parabola curve. Source: `Geom2d_Parabola.hxx`
use occt_core::gp::{GpParab2d, GpPnt2d, GpVec2d, GpTrsf2d, GpAx22d};
use crate::curve::Curve2d;

#[derive(Debug, Clone)]
pub struct Geom2dParabola { pos: GpParab2d }

impl Geom2dParabola {
    pub fn new(p: GpParab2d) -> Self { Self { pos: p } }
    pub fn from_axes(pos: GpAx22d, focal: f64) -> Self { Self { pos: GpParab2d::new(pos, focal) } }
}

impl Curve2d for Geom2dParabola {
    fn d0(&self, u: f64) -> GpPnt2d {
        let f = self.pos.focal;
        let cx = self.pos.pos.point.x(); let cy = self.pos.pos.point.y();
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        let x = u*u/(4.0*f);
        GpPnt2d::new(cx + x*xd.x + u*yd.x, cy + x*xd.y + u*yd.y)
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let f = self.pos.focal;
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        let p = self.d0(u);
        (p, GpVec2d::new(u/(2.0*f)*xd.x + yd.x, u/(2.0*f)*xd.y + yd.y))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let f = self.pos.focal;
        let xd = self.pos.pos.vxdir;
        let p = self.d0(u);
        let d1 = GpVec2d::new(u/(2.0*f)*xd.x + self.pos.pos.vydir.x, u/(2.0*f)*xd.y + self.pos.pos.vydir.y);
        let d2 = GpVec2d::new(xd.x/(2.0*f), xd.y/(2.0*f));
        (p, d1, d2)
    }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 }
    fn clone_dyn(&self) -> Box<dyn crate::curve::Curve2d> { Box::new(self.clone()) }
    fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.focal = -self.pos.focal; }
}
