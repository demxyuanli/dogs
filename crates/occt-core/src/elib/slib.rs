//! Elementary Surfaces Library. Source: `ElSLib.hxx`
use crate::gp::{GpPln,GpCylinder,GpCone,GpSphere,GpTorus,GpPnt,GpVec};

fn pt_add2(o: &crate::gp::GpXyz, a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64) -> GpPnt {
    GpPnt::from_xyz(&o.added(&a.multiplied(sa)).added(&b.multiplied(sb)))
}
fn pt_add3(o: &crate::gp::GpXyz, a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64, c: &crate::gp::GpXyz, sc: f64) -> GpPnt {
    GpPnt::from_xyz(&o.added(&a.multiplied(sa)).added(&b.multiplied(sb)).added(&c.multiplied(sc)))
}

pub fn plane_value(pl: &GpPln, u: f64, v: f64) -> GpPnt {
    pt_add2(&pl.location().coord, pl.pos.x_direction().xyz(), u, pl.pos.y_direction().xyz(), v)
}
pub fn plane_d1(pl: &GpPln, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    (plane_value(pl, u, v), GpVec::from_xyz(pl.pos.x_direction().xyz()), GpVec::from_xyz(pl.pos.y_direction().xyz()))
}

pub fn cylinder_value(cy: &GpCylinder, u: f64, v: f64) -> GpPnt {
    let r = cy.radius;
    pt_add3(&cy.location().coord, cy.pos.x_direction().xyz(), r*u.cos(), cy.pos.y_direction().xyz(), r*u.sin(), cy.pos.direction().xyz(), v)
}

pub fn cone_value(co: &GpCone, u: f64, v: f64) -> GpPnt {
    let r = co.radius; let a = co.semi_angle;
    let r0 = r + v * a.sin(); let z0 = v * a.cos();
    pt_add3(&co.apex().coord, co.pos.x_direction().xyz(), r0*u.cos(), co.pos.y_direction().xyz(), r0*u.sin(), co.pos.direction().xyz(), z0)
}

pub fn sphere_value(s: &GpSphere, u: f64, v: f64) -> GpPnt {
    let r = s.radius;
    let x = r * v.cos() * u.cos(); let y = r * v.cos() * u.sin(); let z = r * v.sin();
    pt_add3(&s.location().coord, s.pos.x_direction().xyz(), x, s.pos.y_direction().xyz(), y, s.pos.direction().xyz(), z)
}

pub fn torus_value(t: &GpTorus, u: f64, v: f64) -> GpPnt {
    let maj_r = t.major_radius; let min_r = t.minor_radius;
    let r0 = maj_r + min_r * v.cos();
    pt_add3(&t.location().coord, t.pos.x_direction().xyz(), r0*u.cos(), t.pos.y_direction().xyz(), r0*u.sin(), t.pos.direction().xyz(), min_r*v.sin())
}
