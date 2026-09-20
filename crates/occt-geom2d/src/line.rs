//! 2D line curve. Source: `Geom2d_Line.hxx`
use occt_core::gp::{GpLin2d, GpPnt2d, GpVec2d, GpDir2d, GpTrsf2d, GpAx2d};
use occt_core::elib::clib2d;
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
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        // `Geom2d_Line::D1`: `ElCLib::LineValue(U)` + constant direction.
        let p = self.d0(u);
        (p, GpVec2d::new(self.pos.pos.vdir.x, self.pos.pos.vdir.y))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, d1) = self.d1(u);
        (p, d1, GpVec2d::zero())
    }
    /// `Geom2d_Line::EvalDN` (`Geom2d_Line.cxx:221-233`) → `ElCLib::LineDN` 2d
    /// (`clib2d.rs:329-335`, `ElCLib.cxx:1049-1055`): `N == 1` is the direction,
    /// any other order the null vector.
    fn eval_dn(&self, _u: f64, n: i32) -> GpVec2d {
        if n < 1 { return GpVec2d::zero(); }
        clib2d::line_dn_ax2d(_u, &self.pos.pos, n)
    }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    fn reverse(&mut self) { self.pos.pos.vdir.reverse(); }
    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
    fn is_line(&self) -> bool { true }
    fn gp_lin2d(&self) -> Option<GpLin2d> { Some(self.pos) }
}
