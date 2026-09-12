//! Elementary Surfaces Library. Source: `ElSLib.hxx`
use std::f64::consts::PI;

use crate::gp::{
    GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpTrsf,
    GpVec,
};
use crate::precision::RESOLUTION;

fn pt_add2(o: &crate::gp::GpXyz, a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64) -> GpPnt {
    GpPnt::from_xyz(&o.added(&a.multiplied(sa)).added(&b.multiplied(sb)))
}
fn pt_add3(
    o: &crate::gp::GpXyz,
    a: &crate::gp::GpXyz,
    sa: f64,
    b: &crate::gp::GpXyz,
    sb: f64,
    c: &crate::gp::GpXyz,
    sc: f64,
) -> GpPnt {
    GpPnt::from_xyz(&o.added(&a.multiplied(sa)).added(&b.multiplied(sb)).added(&c.multiplied(sc)))
}

fn normalize_angle(mut u: f64) -> f64 {
    let two_pi = 2.0 * PI;
    while u < -RESOLUTION {
        u += two_pi;
    }
    while u > two_pi * (1.0 + RESOLUTION) {
        u -= two_pi;
    }
    if u < 0.0 {
        0.0
    } else {
        u
    }
}

fn to_local(pos: &GpAx3, p: &GpPnt) -> GpPnt {
    let mut t = GpTrsf::identity();
    t.set_transformation(pos);
    p.transformed(&t)
}

pub fn plane_value(pl: &GpPln, u: f64, v: f64) -> GpPnt {
    pt_add2(
        &pl.location().coord,
        pl.pos.x_direction().xyz(),
        u,
        pl.pos.y_direction().xyz(),
        v,
    )
}
pub fn plane_d1(pl: &GpPln, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    (
        plane_value(pl, u, v),
        GpVec::from_xyz(pl.pos.x_direction().xyz()),
        GpVec::from_xyz(pl.pos.y_direction().xyz()),
    )
}

/// `ElSLib::PlaneParameters`.
pub fn plane_parameters(pos: &GpAx3, p: &GpPnt) -> (f64, f64) {
    let loc = to_local(pos, p);
    (loc.x(), loc.y())
}

pub fn cylinder_value(cy: &GpCylinder, u: f64, v: f64) -> GpPnt {
    let r = cy.radius;
    pt_add3(
        &cy.location().coord,
        cy.pos.x_direction().xyz(),
        r * u.cos(),
        cy.pos.y_direction().xyz(),
        r * u.sin(),
        cy.pos.direction().xyz(),
        v,
    )
}

/// `ElSLib::CylinderParameters`.
pub fn cylinder_parameters(pos: &GpAx3, p: &GpPnt) -> (f64, f64) {
    let loc = to_local(pos, p);
    (normalize_angle(loc.y().atan2(loc.x())), loc.z())
}

pub fn cone_value(co: &GpCone, u: f64, v: f64) -> GpPnt {
    let r = co.radius;
    let a = co.semi_angle;
    let r0 = r + v * a.sin();
    let z0 = v * a.cos();
    // STEP/OCCT convention (ElSLib::ConeValue): v=0 is the placement plane
    // (radius = RefRadius), the apex is RefRadius/tan(semi_angle) below it
    // along the axis. Base on the placement, not `apex()` (which in this port
    // returns the placement + (r/tan α)·axis, a non-OCCT convention).
    pt_add3(
        &co.location().coord,
        co.pos.x_direction().xyz(),
        r0 * u.cos(),
        co.pos.y_direction().xyz(),
        r0 * u.sin(),
        co.pos.direction().xyz(),
        z0,
    )
}

/// `ElSLib::ConeParameters`.
pub fn cone_parameters(pos: &GpAx3, radius: f64, s_angle: f64, p: &GpPnt) -> (f64, f64) {
    let loc = to_local(pos, p);
    let mut u = if loc.x().abs() < RESOLUTION && loc.y().abs() < RESOLUTION {
        0.0
    } else if -radius > loc.z() * s_angle.tan() {
        (-loc.y()).atan2(-loc.x())
    } else {
        loc.y().atan2(loc.x())
    };
    u = normalize_angle(u);
    let v = s_angle.sin() * (loc.x() * u.cos() + loc.y() * u.sin() - radius)
        + s_angle.cos() * loc.z();
    (u, v)
}

pub fn sphere_value(s: &GpSphere, u: f64, v: f64) -> GpPnt {
    let r = s.radius;
    let x = r * v.cos() * u.cos();
    let y = r * v.cos() * u.sin();
    let z = r * v.sin();
    pt_add3(
        &s.location().coord,
        s.pos.x_direction().xyz(),
        x,
        s.pos.y_direction().xyz(),
        y,
        s.pos.direction().xyz(),
        z,
    )
}

/// `ElSLib::SphereParameters`.
pub fn sphere_parameters(pos: &GpAx3, p: &GpPnt) -> (f64, f64) {
    let loc = to_local(pos, p);
    let l = (loc.x() * loc.x() + loc.y() * loc.y()).sqrt();
    if l < RESOLUTION {
        let v = if loc.z() > 0.0 {
            std::f64::consts::FRAC_PI_2
        } else {
            -std::f64::consts::FRAC_PI_2
        };
        (0.0, v)
    } else {
        (normalize_angle(loc.y().atan2(loc.x())), (loc.z() / l).atan())
    }
}

pub fn torus_value(t: &GpTorus, u: f64, v: f64) -> GpPnt {
    let maj_r = t.major_radius;
    let min_r = t.minor_radius;
    let r0 = maj_r + min_r * v.cos();
    pt_add3(
        &t.location().coord,
        t.pos.x_direction().xyz(),
        r0 * u.cos(),
        t.pos.y_direction().xyz(),
        r0 * u.sin(),
        t.pos.direction().xyz(),
        min_r * v.sin(),
    )
}

/// `ElSLib::TorusParameters`.
pub fn torus_parameters(pos: &GpAx3, major: f64, minor: f64, p: &GpPnt) -> (f64, f64) {
    let loc = to_local(pos, p);
    let mut u = loc.y().atan2(loc.x());
    if major < minor {
        let cosu = u.cos();
        let sinu = u.sin();
        let z2 = loc.z() * loc.z();
        let min_r2 = minor * minor;
        let xm = loc.x() - major * cosu;
        let ym = loc.y() - major * sinu;
        let xp = loc.x() + major * cosu;
        let yp = loc.y() + major * sinu;
        let d1 = (xm * xm + ym * ym + z2 - min_r2).abs();
        let d2 = (xp * xp + yp * yp + z2 - min_r2).abs();
        if d2 < d1 {
            u += PI;
        }
    }
    u = normalize_angle(u);
    let radial = (loc.x() - major * u.cos()) * u.cos() + (loc.y() - major * u.sin()) * u.sin();
    (u, loc.z().atan2(radial))
}

/// `ElSLib::SphereUIso` (`ElSLib.cxx:1738`).
pub fn sphere_u_iso(pos: &GpAx3, radius: f64, u: f64) -> GpCirc {
    let dx = GpVec::from_xyz(pos.x_direction().xyz());
    let dy = GpVec::from_xyz(pos.y_direction().xyz());
    let (su, cu) = u.sin_cos();
    let cx_vec = dx.multiplied_scalar(cu).added(&dy.multiplied_scalar(su));
    let cx = GpDir::from_vec(&cx_vec).unwrap_or_else(|_| *pos.x_direction());
    let n = cx.crossed(&pos.direction()).unwrap_or_else(|_| pos.direction());
    let axes = GpAx2::new(pos.location(), n, cx).unwrap_or_else(|_| pos.ax2());
    GpCirc::new(axes, radius)
}

/// `ElSLib::SphereVIso` (`ElSLib.cxx:1815`).
pub fn sphere_v_iso(pos: &GpAx3, radius: f64, v: f64) -> GpCirc {
    let mut axes = pos.ax2();
    let ve = GpVec::from_xyz(pos.direction().xyz()).multiplied_scalar(radius * v.sin());
    axes.set_location(axes.location().translated_vec(&ve));
    let mut r = radius * v.cos();
    if r < 0.0 {
        axes.set_direction(axes.direction().reversed());
        r = -r;
    }
    GpCirc::new(axes, r)
}

/// `ElSLib::TorusUIso` (`ElSLib.cxx:1751-1765`).
pub fn torus_u_iso(pos: &GpAx3, major: f64, minor: f64, u: f64) -> GpCirc {
    let dx = GpVec::from_xyz(pos.x_direction().xyz());
    let dy = GpVec::from_xyz(pos.y_direction().xyz());
    let (su, cu) = u.sin_cos();
    let cx_vec = dx.multiplied_scalar(cu).added(&dy.multiplied_scalar(su));
    let cx = GpDir::from_vec(&cx_vec).unwrap_or_else(|_| *pos.x_direction());
    let n = cx.crossed(&pos.direction()).unwrap_or_else(|_| pos.direction());
    let mut axes = GpAx2::new(pos.location(), n, cx).unwrap_or_else(|_| pos.ax2());
    let ve = GpVec::from_xyz(cx.xyz()).multiplied_scalar(major);
    axes.set_location(axes.location().translated_vec(&ve));
    GpCirc::new(axes, minor)
}

/// `ElSLib::TorusVIso` (`ElSLib.cxx:1836-1853`).
pub fn torus_v_iso(pos: &GpAx3, major: f64, minor: f64, v: f64) -> GpCirc {
    let mut axes = pos.ax2();
    let ve = GpVec::from_xyz(pos.direction().xyz()).multiplied_scalar(minor * v.sin());
    axes.set_location(axes.location().translated_vec(&ve));
    let mut r = major + minor * v.cos();
    if r < 0.0 {
        let x = axes.x_direction().reversed();
        axes = GpAx2::new(axes.location(), axes.direction(), x).unwrap_or(axes);
        r = -r;
    }
    GpCirc::new(axes, r)
}

/// `ElSLib::CylinderVIso` (`ElSLib.cxx:1781-1788`).
pub fn cylinder_v_iso(pos: &GpAx3, radius: f64, v: f64) -> GpCirc {
    let mut axes = pos.ax2();
    let ve = GpVec::from_xyz(pos.direction().xyz()).multiplied_scalar(v);
    axes.set_location(axes.location().translated_vec(&ve));
    GpCirc::new(axes, radius)
}

/// `ElSLib::CylinderUIso` (`ElSLib.cxx:1716-1723`).
pub fn cylinder_u_iso(pos: &GpAx3, radius: f64, u: f64) -> GpLin {
    let cy = GpCylinder {
        pos: *pos,
        radius,
    };
    let p = cylinder_value(&cy, u, 0.0);
    GpLin::from_pnt_dir(p, pos.direction())
}

/// `ElSLib::ConeUIso` (`ElSLib.cxx:1727-1733`).
pub fn cone_u_iso(pos: &GpAx3, radius: f64, s_angle: f64, u: f64) -> GpLin {
    let co = GpCone {
        pos: *pos,
        radius,
        semi_angle: s_angle,
    };
    let p = cone_value(&co, u, 0.0);
    let (su, cu) = u.sin_cos();
    let sa = s_angle.sin();
    let ca = s_angle.cos();
    let dv = GpVec::from_xyz(
        &pos.x_direction()
            .xyz()
            .multiplied(sa * cu)
            .added(&pos.y_direction().xyz().multiplied(sa * su))
            .added(&pos.direction().xyz().multiplied(ca)),
    );
    let dir = GpDir::from_vec(&dv).unwrap_or(pos.direction());
    GpLin::from_pnt_dir(p, dir)
}

/// `ElSLib::ConeVIso` (`ElSLib.cxx:1793-1810`).
pub fn cone_v_iso(pos: &GpAx3, radius: f64, s_angle: f64, v: f64) -> GpCirc {
    let mut axes = *pos;
    let ve = GpVec::from_xyz(pos.direction().xyz()).multiplied_scalar(v * s_angle.cos());
    axes.set_location(axes.location().translated_vec(&ve));
    let mut r = radius + v * s_angle.sin();
    if r < 0.0 {
        // OCCT `XReverse` then `YReverse` (`gp_Ax3.hxx:125-128`) reverse vx
        // and vy independently. `GpAx3::x_reverse`/`y_reverse` currently reverse
        // both axes and would cancel here.
        axes.vxdir.reverse();
        axes.vydir.reverse();
        r = -r;
    }
    GpCirc::new(axes.ax2(), r)
}
