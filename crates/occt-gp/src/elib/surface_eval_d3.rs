//! `ElSLib::*D3` (`ElSLib.cxx:1147-1543`).
//!
//! Order matches the Fortran-style outputs of ConeD3 / CylinderD3 / SphereD3 /
//! TorusD3: point, first, second, then `Vuuu`, `Vvvv`, `Vuuv`, `Vuvv`.

use crate::gp::{GpAx3, GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTorus, GpVec};

/// Third-order surface jet (`Geom_Surface::ResD3` field order).
#[derive(Clone, Debug)]
pub struct SurfD3 {
    pub point: GpPnt,
    pub d1u: GpVec,
    pub d1v: GpVec,
    pub d2u: GpVec,
    pub d2v: GpVec,
    pub d2uv: GpVec,
    pub d3u: GpVec,
    pub d3v: GpVec,
    pub d3_uuv: GpVec,
    pub d3_uvv: GpVec,
}

fn vec3(x: f64, y: f64, z: f64) -> GpVec {
    GpVec::new(x, y, z)
}

fn pnt3(x: f64, y: f64, z: f64) -> GpPnt {
    GpPnt::new(x, y, z)
}

/// `Geom_Plane::EvalD3` (`Geom_Plane.cxx:226-238`): `PlaneD1` plus zero higher
/// derivatives. There is no `ElSLib::PlaneD3`.
pub fn plane_d3(pl: &GpPln, u: f64, v: f64) -> SurfD3 {
    let (point, d1u, d1v) = super::surface_eval::plane_d1(pl, u, v);
    SurfD3 {
        point,
        d1u,
        d1v,
        d2u: GpVec::zero(),
        d2v: GpVec::zero(),
        d2uv: GpVec::zero(),
        d3u: GpVec::zero(),
        d3v: GpVec::zero(),
        d3_uuv: GpVec::zero(),
        d3_uvv: GpVec::zero(),
    }
}

/// `ElSLib::ConeD3` (`ElSLib.cxx:1147-1244`).
pub fn cone_d3(co: &GpCone, u: f64, v: f64) -> SurfD3 {
    cone_d3_ax3(&co.position(), co.radius(), co.semi_angle(), u, v)
}

fn cone_d3_ax3(pos: &GpAx3, radius: f64, s_angle: f64, u: f64, v: f64) -> SurfD3 {
    let xdir = pos.x_direction().xyz();
    let ydir = pos.y_direction().xyz();
    let zdir_owned = pos.direction();
    let zdir = zdir_owned.xyz();
    let ploc = pos.location();
    let (px, py, pz) = (ploc.x(), ploc.y(), ploc.z());
    let (cos_u, sin_u) = (u.cos(), u.sin());
    let (cos_a, sin_a) = (s_angle.cos(), s_angle.sin());
    let r = radius + v * sin_a;
    let a3 = v * cos_a;
    let a1 = r * cos_u;
    let a2 = r * sin_u;
    let r1 = sin_a * cos_u;
    let r2 = sin_a * sin_u;
    let som1_x = a1 * xdir.x() + a2 * ydir.x();
    let som1_y = a1 * xdir.y() + a2 * ydir.y();
    let som1_z = a1 * xdir.z() + a2 * ydir.z();
    let som2_x = r1 * xdir.x() + r2 * ydir.x();
    let som2_y = r1 * xdir.y() + r2 * ydir.y();
    let som2_z = r1 * xdir.z() + r2 * ydir.z();
    let dif1_x = a2 * xdir.x() - a1 * ydir.x();
    let dif1_y = a2 * xdir.y() - a1 * ydir.y();
    let dif1_z = a2 * xdir.z() - a1 * ydir.z();
    SurfD3 {
        point: pnt3(
            som1_x + a3 * zdir.x() + px,
            som1_y + a3 * zdir.y() + py,
            som1_z + a3 * zdir.z() + pz,
        ),
        d1u: vec3(-dif1_x, -dif1_y, -dif1_z),
        d1v: vec3(
            som2_x + cos_a * zdir.x(),
            som2_y + cos_a * zdir.y(),
            som2_z + cos_a * zdir.z(),
        ),
        d2u: vec3(-som1_x, -som1_y, -som1_z),
        d2v: GpVec::zero(),
        d2uv: vec3(
            -r2 * xdir.x() + r1 * ydir.x(),
            -r2 * xdir.y() + r1 * ydir.y(),
            -r2 * xdir.z() + r1 * ydir.z(),
        ),
        d3u: vec3(dif1_x, dif1_y, dif1_z),
        d3v: GpVec::zero(),
        d3_uuv: vec3(-som2_x, -som2_y, -som2_z),
        d3_uvv: GpVec::zero(),
    }
}

/// `ElSLib::CylinderD3` (`ElSLib.cxx:1246-1303`).
pub fn cylinder_d3(cy: &GpCylinder, u: f64, v: f64) -> SurfD3 {
    cylinder_d3_ax3(&cy.position(), cy.radius(), u, v)
}

fn cylinder_d3_ax3(pos: &GpAx3, radius: f64, u: f64, v: f64) -> SurfD3 {
    let xdir = pos.x_direction().xyz();
    let ydir = pos.y_direction().xyz();
    let zdir_owned = pos.direction();
    let zdir = zdir_owned.xyz();
    let ploc = pos.location();
    let (px, py, pz) = (ploc.x(), ploc.y(), ploc.z());
    let a1 = radius * u.cos();
    let a2 = radius * u.sin();
    let som1_x = a1 * xdir.x() + a2 * ydir.x();
    let som1_y = a1 * xdir.y() + a2 * ydir.y();
    let som1_z = a1 * xdir.z() + a2 * ydir.z();
    let dif1_x = a2 * xdir.x() - a1 * ydir.x();
    let dif1_y = a2 * xdir.y() - a1 * ydir.y();
    let dif1_z = a2 * xdir.z() - a1 * ydir.z();
    SurfD3 {
        point: pnt3(
            som1_x + v * zdir.x() + px,
            som1_y + v * zdir.y() + py,
            som1_z + v * zdir.z() + pz,
        ),
        d1u: vec3(-dif1_x, -dif1_y, -dif1_z),
        d1v: vec3(zdir.x(), zdir.y(), zdir.z()),
        d2u: vec3(-som1_x, -som1_y, -som1_z),
        d2v: GpVec::zero(),
        d2uv: GpVec::zero(),
        d3u: vec3(dif1_x, dif1_y, dif1_z),
        d3v: GpVec::zero(),
        d3_uuv: GpVec::zero(),
        d3_uvv: GpVec::zero(),
    }
}

/// `ElSLib::SphereD3` (`ElSLib.cxx:1305-1401`).
pub fn sphere_d3(s: &GpSphere, u: f64, v: f64) -> SurfD3 {
    sphere_d3_ax3(s.position(), s.radius(), u, v)
}

fn sphere_d3_ax3(pos: &GpAx3, radius: f64, u: f64, v: f64) -> SurfD3 {
    let xdir = pos.x_direction().xyz();
    let ydir = pos.y_direction().xyz();
    let zdir_owned = pos.direction();
    let zdir = zdir_owned.xyz();
    let ploc = pos.location();
    let (px, py, pz) = (ploc.x(), ploc.y(), ploc.z());
    let (cos_u, sin_u) = (u.cos(), u.sin());
    let r1 = radius * v.cos();
    let r2 = radius * v.sin();
    let a1 = r1 * cos_u;
    let a2 = r1 * sin_u;
    let a3 = r2 * cos_u;
    let a4 = r2 * sin_u;
    let som1_x = a1 * xdir.x() + a2 * ydir.x();
    let som1_y = a1 * xdir.y() + a2 * ydir.y();
    let som1_z = a1 * xdir.z() + a2 * ydir.z();
    let som3_x = a3 * xdir.x() + a4 * ydir.x();
    let som3_y = a3 * xdir.y() + a4 * ydir.y();
    let som3_z = a3 * xdir.z() + a4 * ydir.z();
    let dif1_x = a2 * xdir.x() - a1 * ydir.x();
    let dif1_y = a2 * xdir.y() - a1 * ydir.y();
    let dif1_z = a2 * xdir.z() - a1 * ydir.z();
    let r1z_x = r1 * zdir.x();
    let r1z_y = r1 * zdir.y();
    let r1z_z = r1 * zdir.z();
    let r2z_x = r2 * zdir.x();
    let r2z_y = r2 * zdir.y();
    let r2z_z = r2 * zdir.z();
    SurfD3 {
        point: pnt3(
            som1_x + r2z_x + px,
            som1_y + r2z_y + py,
            som1_z + r2z_z + pz,
        ),
        d1u: vec3(-dif1_x, -dif1_y, -dif1_z),
        d1v: vec3(-som3_x + r1z_x, -som3_y + r1z_y, -som3_z + r1z_z),
        d2u: vec3(-som1_x, -som1_y, -som1_z),
        d2v: vec3(-som1_x - r2z_x, -som1_y - r2z_y, -som1_z - r2z_z),
        d2uv: vec3(
            a4 * xdir.x() - a3 * ydir.x(),
            a4 * xdir.y() - a3 * ydir.y(),
            a4 * xdir.z() - a3 * ydir.z(),
        ),
        d3u: vec3(dif1_x, dif1_y, dif1_z),
        d3v: vec3(som3_x - r1z_x, som3_y - r1z_y, som3_z - r1z_z),
        d3_uuv: vec3(som3_x, som3_y, som3_z),
        d3_uvv: vec3(dif1_x, dif1_y, dif1_z),
    }
}

/// `ElSLib::TorusD3` (`ElSLib.cxx:1403-1543`).
pub fn torus_d3(t: &GpTorus, u: f64, v: f64) -> SurfD3 {
    torus_d3_ax3(
        t.position(),
        t.major_radius(),
        t.minor_radius(),
        u,
        v,
    )
}

fn torus_d3_ax3(pos: &GpAx3, major_radius: f64, minor_radius: f64, u: f64, v: f64) -> SurfD3 {
    let xdir = pos.x_direction().xyz();
    let ydir = pos.y_direction().xyz();
    let zdir_owned = pos.direction();
    let zdir = zdir_owned.xyz();
    let ploc = pos.location();
    let (px, py, pz) = (ploc.x(), ploc.y(), ploc.z());
    let (cos_u, sin_u) = (u.cos(), u.sin());
    let r1 = minor_radius * v.cos();
    let r2 = minor_radius * v.sin();
    let r = major_radius + r1;
    let mut a1 = r * cos_u;
    let mut a2 = r * sin_u;
    let mut a3 = r2 * cos_u;
    let mut a4 = r2 * sin_u;
    let mut a5 = r1 * cos_u;
    let mut a6 = r1 * sin_u;
    // OCC620 snap (`ElSLib.cxx:1466-1496`).
    let eps = 10.0 * (minor_radius + major_radius) * f64::EPSILON;
    if a1.abs() <= eps {
        a1 = 0.0;
    }
    if a2.abs() <= eps {
        a2 = 0.0;
    }
    if a3.abs() <= eps {
        a3 = 0.0;
    }
    if a4.abs() <= eps {
        a4 = 0.0;
    }
    if a5.abs() <= eps {
        a5 = 0.0;
    }
    if a6.abs() <= eps {
        a6 = 0.0;
    }
    let som1_x = a1 * xdir.x() + a2 * ydir.x();
    let som1_y = a1 * xdir.y() + a2 * ydir.y();
    let som1_z = a1 * xdir.z() + a2 * ydir.z();
    let som3_x = a3 * xdir.x() + a4 * ydir.x();
    let som3_y = a3 * xdir.y() + a4 * ydir.y();
    let som3_z = a3 * xdir.z() + a4 * ydir.z();
    let dif1_x = a2 * xdir.x() - a1 * ydir.x();
    let dif1_y = a2 * xdir.y() - a1 * ydir.y();
    let dif1_z = a2 * xdir.z() - a1 * ydir.z();
    let r1z_x = r1 * zdir.x();
    let r1z_y = r1 * zdir.y();
    let r1z_z = r1 * zdir.z();
    let r2z_x = r2 * zdir.x();
    let r2z_y = r2 * zdir.y();
    let r2z_z = r2 * zdir.z();
    SurfD3 {
        point: pnt3(
            som1_x + r2z_x + px,
            som1_y + r2z_y + py,
            som1_z + r2z_z + pz,
        ),
        d1u: vec3(-dif1_x, -dif1_y, -dif1_z),
        d1v: vec3(-som3_x + r1z_x, -som3_y + r1z_y, -som3_z + r1z_z),
        d2u: vec3(-som1_x, -som1_y, -som1_z),
        d2v: vec3(
            -a5 * xdir.x() - a6 * ydir.x() - r2z_x,
            -a5 * xdir.y() - a6 * ydir.y() - r2z_y,
            -a5 * xdir.z() - a6 * ydir.z() - r2z_z,
        ),
        d2uv: vec3(
            a4 * xdir.x() - a3 * ydir.x(),
            a4 * xdir.y() - a3 * ydir.y(),
            a4 * xdir.z() - a3 * ydir.z(),
        ),
        d3u: vec3(dif1_x, dif1_y, dif1_z),
        d3v: vec3(som3_x - r1z_x, som3_y - r1z_y, som3_z - r1z_z),
        d3_uuv: vec3(som3_x, som3_y, som3_z),
        d3_uvv: vec3(
            a6 * xdir.x() - a5 * ydir.x(),
            a6 * xdir.y() - a5 * ydir.y(),
            a6 * xdir.z() - a5 * ydir.z(),
        ),
    }
}
