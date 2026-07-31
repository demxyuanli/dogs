//! Elementary Curves Library. Source: `ElCLib.hxx`
use crate::gp::{GpLin,GpCirc,GpElips,GpHypr,GpParab,GpPnt,GpVec,GpPnt2d,GpVec2d,GpLin2d,GpCirc2d,GpElips2d};

// Helper: GpVec from scaled & added GpXyz refs
fn vec_add(a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64) -> GpVec {
    let r = a.multiplied(sa).added(&b.multiplied(sb));
    GpVec::from_xyz(&r)
}

fn pt_add(org: &crate::gp::GpXyz, a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64) -> GpPnt {
    let r = org.added(&a.multiplied(sa)).added(&b.multiplied(sb));
    GpPnt::from_xyz(&r)
}
fn pt_add3(org: &crate::gp::GpXyz, a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64, c: &crate::gp::GpXyz, sc: f64) -> GpPnt {
    let r = org.added(&a.multiplied(sa)).added(&b.multiplied(sb)).added(&c.multiplied(sc));
    GpPnt::from_xyz(&r)
}

pub fn line_value(l: &GpLin, u: f64) -> GpPnt {
    GpPnt::from_xyz(&l.pos.loc.coord.added(&l.pos.vdir.xyz().multiplied(u)))
}
pub fn line_d1(l: &GpLin, _u: f64) -> (GpPnt, GpVec) {
    (line_value(l, 0.0), GpVec::from_xyz(l.pos.vdir.xyz()))
}
pub fn line_d2(l: &GpLin, _u: f64) -> (GpPnt, GpVec, GpVec) {
    let (p, d1) = line_d1(l, 0.0); (p, d1, GpVec::zero())
}

pub fn circle_value(c: &GpCirc, u: f64) -> GpPnt {
    let r = c.radius;
    pt_add(&c.pos.location().coord, c.pos.x_direction().xyz(), r*u.cos(), c.pos.y_direction().xyz(), r*u.sin())
}
pub fn circle_d1(c: &GpCirc, u: f64) -> (GpPnt, GpVec) {
    let r = c.radius;
    let p = circle_value(c, u);
    (p, vec_add(c.pos.x_direction().xyz(), -r*u.sin(), c.pos.y_direction().xyz(), r*u.cos()))
}
pub fn circle_d2(c: &GpCirc, u: f64) -> (GpPnt, GpVec, GpVec) {
    let r = c.radius;
    let p = circle_value(c, u);
    let d1 = vec_add(c.pos.x_direction().xyz(), -r*u.sin(), c.pos.y_direction().xyz(), r*u.cos());
    let d2 = vec_add(c.pos.x_direction().xyz(), -r*u.cos(), c.pos.y_direction().xyz(), -r*u.sin());
    (p, d1, d2)
}

pub fn ellipse_value(e: &GpElips, u: f64) -> GpPnt {
    let a = e.major_radius; let b = e.minor_radius;
    pt_add(&e.location().coord, e.pos.x_direction().xyz(), a*u.cos(), e.pos.y_direction().xyz(), -b*u.sin())
}
pub fn ellipse_d1(e: &GpElips, u: f64) -> (GpPnt, GpVec) {
    let a = e.major_radius; let b = e.minor_radius;
    let p = ellipse_value(e, u);
    (p, vec_add(e.pos.x_direction().xyz(), -a*u.sin(), e.pos.y_direction().xyz(), -b*u.cos()))
}
pub fn ellipse_d2(e: &GpElips, u: f64) -> (GpPnt, GpVec, GpVec) {
    let a = e.major_radius; let b = e.minor_radius;
    let p = ellipse_value(e, u);
    let d1 = vec_add(e.pos.x_direction().xyz(), -a*u.sin(), e.pos.y_direction().xyz(), -b*u.cos());
    let d2 = vec_add(e.pos.x_direction().xyz(), -a*u.cos(), e.pos.y_direction().xyz(), b*u.sin());
    (p, d1, d2)
}

pub fn hyperbola_value(h: &GpHypr, u: f64) -> GpPnt {
    let a = h.major_radius; let b = h.minor_radius;
    pt_add(&h.location().coord, h.pos.x_direction().xyz(), a*u.cosh(), h.pos.y_direction().xyz(), b*u.sinh())
}
pub fn hyperbola_d1(h: &GpHypr, u: f64) -> (GpPnt, GpVec) {
    let a = h.major_radius; let b = h.minor_radius;
    let p = hyperbola_value(h, u);
    (p, vec_add(h.pos.x_direction().xyz(), a*u.sinh(), h.pos.y_direction().xyz(), b*u.cosh()))
}

pub fn parabola_value(p: &GpParab, u: f64) -> GpPnt {
    let f = p.focal;
    pt_add(&p.location().coord, p.pos.x_direction().xyz(), u*u/(4.0*f), p.pos.y_direction().xyz(), u)
}
pub fn parabola_d1(p: &GpParab, u: f64) -> (GpPnt, GpVec) {
    let f = p.focal;
    let pt = parabola_value(p, u);
    (pt, vec_add(p.pos.x_direction().xyz(), u/(2.0*f), p.pos.y_direction().xyz(), 1.0))
}

// 2D
pub fn line2d_value(l: &GpLin2d, u: f64) -> GpPnt2d { GpPnt2d::new(l.pos.loc.x()+u*l.pos.vdir.x, l.pos.loc.y()+u*l.pos.vdir.y) }
pub fn circle2d_value(c: &GpCirc2d, u: f64) -> GpPnt2d {
    let r=c.radius; let cx=c.pos.point.x(); let cy=c.pos.point.y(); let xd=c.pos.vxdir; let yd=c.pos.vydir;
    GpPnt2d::new(cx+r*u.cos()*xd.x-r*u.sin()*yd.x, cy+r*u.cos()*xd.y-r*u.sin()*yd.y)
}
pub fn ellipse2d_value(e: &GpElips2d, u: f64) -> GpPnt2d {
    let a=e.major_radius; let b=e.minor_radius; let cx=e.pos.point.x(); let cy=e.pos.point.y(); let xd=e.pos.vxdir; let yd=e.pos.vydir;
    GpPnt2d::new(cx+a*u.cos()*xd.x-b*u.sin()*yd.x, cy+a*u.cos()*xd.y-b*u.sin()*yd.y)
}
