//! 2D ellipse curve. Source: `Geom2d_Ellipse.hxx`
use occt_core::elib::{clib, clib2d};
use occt_core::gp::{GpElips2d, GpPnt2d, GpVec2d, GpTrsf2d, GpAx22d};
use crate::curve::Curve2d;

#[derive(Debug, Clone)]
pub struct Geom2dEllipse { pos: GpElips2d }

impl Geom2dEllipse {
    pub fn new(e: GpElips2d) -> Self { Self { pos: e } }
    pub fn from_axes(pos: GpAx22d, major: f64, minor: f64) -> Self { Self { pos: GpElips2d::new(pos, major, minor) } }
}

impl Curve2d for Geom2dEllipse {
    /// `ElCLib::EllipseValue(U, gp_Ax22d, Major, Minor)` (`cxx:543-555`):
    /// `P = Loc + Major*cos(U)*Xd + Minor*sin(U)*Yd`.
    fn d0(&self, u: f64) -> GpPnt2d {
        let a = self.pos.major_radius;
        let b = self.pos.minor_radius;
        let cx = self.pos.pos.point.x(); let cy = self.pos.pos.point.y();
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        GpPnt2d::new(cx + a*u.cos()*xd.x + b*u.sin()*yd.x, cy + a*u.cos()*xd.y + b*u.sin()*yd.y)
    }
    /// `ElCLib::EllipseD1(U, gp_Ax22d, Major, Minor)` (`cxx:623-650`):
    /// `V1 = -Major*sin(U)*Xd + Minor*cos(U)*Yd`.
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let a = self.pos.major_radius; let b = self.pos.minor_radius;
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        let p = self.d0(u);
        (p, GpVec2d::new(-a*u.sin()*xd.x + b*u.cos()*yd.x, -a*u.sin()*xd.y + b*u.cos()*yd.y))
    }
    /// `ElCLib::EllipseD2(U, gp_Ax22d, Major, Minor)`:
    /// `V2 = -Major*cos(U)*Xd - Minor*sin(U)*Yd`.
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let a = self.pos.major_radius; let b = self.pos.minor_radius;
        let xd = self.pos.pos.vxdir; let yd = self.pos.pos.vydir;
        let p = self.d0(u);
        let d1 = GpVec2d::new(-a*u.sin()*xd.x + b*u.cos()*yd.x, -a*u.sin()*xd.y + b*u.cos()*yd.y);
        let d2 = GpVec2d::new(-a*u.cos()*xd.x - b*u.sin()*yd.x, -a*u.cos()*xd.y - b*u.sin()*yd.y);
        (p, d1, d2)
    }
    /// `Geom2d_Ellipse::EvalD3` → `ElCLib::EllipseD3` (`ElCLib.cxx:843-875`).
    /// Needed by `Geom2d_OffsetCurve`'s `CalculateD2` (`D2Ndir`).
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        clib::ellipse2d_d3(&self.pos, u)
    }
    /// `Geom2d_Ellipse::EvalDN` (`Geom2d_Ellipse.cxx:297-304`) →
    /// `ElCLib::EllipseDN` 2d (`clib2d.rs:363-391`, `ElCLib.cxx:1095-1133`).
    /// `N < 1` returns a zero vector instead of the OCCT throw.
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        if n < 1 { return GpVec2d::zero(); }
        clib2d::ellipse_dn_ax22d(u, &self.pos.pos, self.pos.major_radius, self.pos.minor_radius, n)
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn is_periodic(&self) -> bool { true }
    fn period(&self) -> f64 { 2.0 * std::f64::consts::PI }
    fn continuity(&self) -> u8 { 6 }
    fn clone_dyn(&self) -> Box<dyn crate::curve::Curve2d> { Box::new(self.clone()) }
    fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    fn reverse(&mut self) { std::mem::swap(&mut self.pos.major_radius, &mut self.pos.minor_radius); }
}
