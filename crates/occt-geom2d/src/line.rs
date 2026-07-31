//! 2D line curve. Source: `Geom2d_Line.hxx`
use occt_core::gp::{GpLin2d, GpPnt2d, GpVec2d, GpDir2d, GpTrsf2d, GpAx2d};
use crate::curve::Curve2d;

#[derive(Debug, Clone)]
pub struct Geom2dLine { pos: GpLin2d }

impl Geom2dLine {
    pub fn new(a: GpAx2d) -> Self { Self { pos: GpLin2d::new(a) } }
    pub fn from_pnt_dir(p: GpPnt2d, d: GpDir2d) -> Self { Self { pos: GpLin2d::from_pnt_dir(p, d) } }
    pub fn lin(&self) -> &GpLin2d { &self.pos }
}

impl Curve2d for Geom2dLine {
    fn d0(&self, u: f64) -> GpPnt2d {
        let lx = self.pos.pos.loc.x(); let ly = self.pos.pos.loc.y();
        let dx = self.pos.pos.vdir.x; let dy = self.pos.pos.vdir.y;
        GpPnt2d::new(lx + u * dx, ly + u * dy)
    }
    fn d1(&self, _u: f64) -> (GpPnt2d, GpVec2d) {
        let p = self.d0(0.0);
        (p, GpVec2d::new(self.pos.pos.vdir.x, self.pos.pos.vdir.y))
    }
    fn d2(&self, _u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1) = self.d1(0.0);
        (p, d1, GpVec2d::zero())
    }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 }
    fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.pos.vdir.reverse(); }
    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
}
