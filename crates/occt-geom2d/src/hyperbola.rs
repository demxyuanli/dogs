//! 2D hyperbola curve. Source: `Geom2d_Hyperbola.hxx`
use occt_core::elib::{clib, clib2d};
use occt_core::gp::{GpHypr2d, GpPnt2d, GpVec2d, GpTrsf2d, GpAx22d};
use crate::curve::Curve2d;

#[derive(Debug, Clone)]
pub struct Geom2dHyperbola { pos: GpHypr2d }

impl Geom2dHyperbola {
    pub fn new(h: GpHypr2d) -> Self { Self { pos: h } }
    pub fn from_axes(pos: GpAx22d, major: f64, minor: f64) -> Self { Self { pos: GpHypr2d::new(pos, major, minor) } }
}

impl Curve2d for Geom2dHyperbola {
    fn d0(&self, u: f64) -> GpPnt2d {
        let a = self.pos.major_radius; let b = self.pos.minor_radius;
        let cx = self.pos.pos.point.x(); let cy = self.pos.pos.point.y();
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        GpPnt2d::new(cx + a*u.cosh()*xd.x + b*u.sinh()*yd.x, cy + a*u.cosh()*xd.y + b*u.sinh()*yd.y)
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let a = self.pos.major_radius; let b = self.pos.minor_radius;
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        let p = self.d0(u);
        (p, GpVec2d::new(a*u.sinh()*xd.x + b*u.cosh()*yd.x, a*u.sinh()*xd.y + b*u.cosh()*yd.y))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let a = self.pos.major_radius; let b = self.pos.minor_radius;
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        let p = self.d0(u);
        let d1 = GpVec2d::new(a*u.sinh()*xd.x + b*u.cosh()*yd.x, a*u.sinh()*xd.y + b*u.cosh()*yd.y);
        let d2 = GpVec2d::new(a*u.cosh()*xd.x + b*u.sinh()*yd.x, a*u.cosh()*xd.y + b*u.sinh()*yd.y);
        (p, d1, d2)
    }
    /// `Geom2d_Hyperbola::EvalD3` → `ElCLib::HyperbolaD3` (`ElCLib.cxx:878-905`),
    /// where OCCT sets `V3 = V1`. Needed by `Geom2d_OffsetCurve`'s `CalculateD2`.
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        clib::hyperbola2d_d3(&self.pos, u)
    }
    /// `Geom2d_Hyperbola::EvalDN` (`Geom2d_Hyperbola.cxx:329-336`) →
    /// `ElCLib::HyperbolaDN` 2d (`clib2d.rs:394-413`, `ElCLib.cxx:1137-1158`).
    /// `N < 1` returns a zero vector instead of the OCCT throw.
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        if n < 1 { return GpVec2d::zero(); }
        clib2d::hyperbola_dn_ax22d(u, &self.pos.pos, self.pos.major_radius, self.pos.minor_radius, n)
    }
    fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    fn last_parameter(&self) -> f64 { f64::INFINITY }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn clone_dyn(&self) -> Box<dyn crate::curve::Curve2d> { Box::new(self.clone()) }
    fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    fn reverse(&mut self) { /* swap branches */ }
}
