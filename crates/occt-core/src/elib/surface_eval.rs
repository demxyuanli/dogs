//! Elementary surface evaluation with full derivatives and normals.
//! Source: `ElSLib.cxx` D1/D2/Norm functions.
use crate::gp::{GpAx3, GpPln, GpCylinder, GpCone, GpSphere, GpTorus, GpPnt, GpVec};
use crate::elib::slib;

/// Evaluate plane D1: point + du + dv.
pub fn plane_d1(pl: &GpPln, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    slib::plane_d1(pl, u, v)
}

/// Plane normal = Xdir × Ydir (constant).
pub fn plane_normal(pl: &GpPln) -> GpVec {
    GpVec::from_xyz(&pl.pos.x_direction().xyz().crossed(pl.pos.y_direction().xyz()))
}

/// Cylinder D1: point + du (tangential) + dv (axial).
pub fn cylinder_d1(cy: &GpCylinder, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    let r = cy.radius;
    let du = GpVec::from_xyz(&cy.pos.x_direction().xyz().multiplied(-r*u.sin())
        .added(&cy.pos.y_direction().xyz().multiplied(r*u.cos())));
    let dv = GpVec::from_xyz(cy.pos.direction().xyz());
    (slib::cylinder_value(cy, u, v), du, dv)
}

/// Cylinder outward normal (radial direction).
pub fn cylinder_normal(cy: &GpCylinder, u: f64) -> GpVec {
    GpVec::from_xyz(&cy.pos.x_direction().xyz().multiplied(u.cos())
        .added(&cy.pos.y_direction().xyz().multiplied(u.sin())))
}

/// Sphere D1: point + du + dv.
pub fn sphere_d1(s: &GpSphere, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    let r = s.radius;
    let cv = v.cos(); let sv = v.sin();
    let du = GpVec::from_xyz(&s.pos.x_direction().xyz().multiplied(-r*cv*u.sin())
        .added(&s.pos.y_direction().xyz().multiplied(r*cv*u.cos())));
    let dv = GpVec::from_xyz(&s.pos.x_direction().xyz().multiplied(-r*sv*u.cos())
        .added(&s.pos.y_direction().xyz().multiplied(-r*sv*u.sin()))
        .added(&s.pos.direction().xyz().multiplied(r*cv)));
    (slib::sphere_value(s, u, v), du, dv)
}

/// Sphere outward normal (radial from center to point).
pub fn sphere_normal(s: &GpSphere, u: f64, v: f64) -> GpVec {
    let p = slib::sphere_value(s, u, v);
    let n = p.coord.subtracted(&s.location().coord);
    let m = n.modulus();
    if m > 1e-30 { GpVec::from_xyz(&n.divided(m)) } else { GpVec::zero() }
}

/// Torus D1.
pub fn torus_d1(t: &GpTorus, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    let r0 = t.major_radius + t.minor_radius * v.cos();
    let du = GpVec::from_xyz(&t.pos.x_direction().xyz().multiplied(-r0*u.sin())
        .added(&t.pos.y_direction().xyz().multiplied(r0*u.cos())));
    let dv = GpVec::from_xyz(&t.pos.x_direction().xyz().multiplied(-t.minor_radius*v.sin()*u.cos())
        .added(&t.pos.y_direction().xyz().multiplied(-t.minor_radius*v.sin()*u.sin()))
        .added(&t.pos.direction().xyz().multiplied(t.minor_radius*v.cos())));
    (slib::torus_value(t, u, v), du, dv)
}

/// General surface normal from D1 (du × dv normalized).
pub fn normal_from_d1(du: &GpVec, dv: &GpVec) -> GpVec {
    let n = du.xyz().crossed(dv.xyz());
    let m = n.modulus();
    if m > 1e-30 { GpVec::from_xyz(&n.divided(m)) } else { GpVec::zero() }
}

/// Cone D1. r0 = radius + v*sin(semi_angle).
pub fn cone_d1(co: &GpCone, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    let a = co.semi_angle;
    let r0 = co.radius + v * a.sin();
    let du = GpVec::from_xyz(&co.pos.x_direction().xyz().multiplied(-r0*u.sin())
        .added(&co.pos.y_direction().xyz().multiplied(r0*u.cos())));
    let dv = GpVec::from_xyz(&co.pos.x_direction().xyz().multiplied(a.sin()*u.cos())
        .added(&co.pos.y_direction().xyz().multiplied(a.sin()*u.sin()))
        .added(&co.pos.direction().xyz().multiplied(a.cos())));
    (slib::cone_value(co, u, v), du, dv)
}

/// `ElSLib::PlaneD2`: the plane is affine, so all second derivatives vanish.
pub fn plane_d2(pl: &GpPln, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
    let (p, du, dv) = plane_d1(pl, u, v);
    (p, du, dv, GpVec::zero(), GpVec::zero(), GpVec::zero())
}

/// `ElSLib::CylinderD2`: `Vuu = -R·Vxy`, `Vvv = Vuv = 0`.
pub fn cylinder_d2(cy: &GpCylinder, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
    let r = cy.radius;
    let duu = GpVec::from_xyz(&cy.pos.x_direction().xyz().multiplied(-r*u.cos())
        .added(&cy.pos.y_direction().xyz().multiplied(-r*u.sin())));
    let (p, du, dv) = cylinder_d1(cy, u, v);
    (p, du, dv, duu, GpVec::zero(), GpVec::zero())
}

/// `ElSLib::ConeD2` (`ElSLib.cxx:867-...`). With `R = Radius + V·sinA`,
/// `Vxy = cosU·X + sinU·Y`, `DVxy = -sinU·X + cosU·Y`:
/// `Vuu = -R·Vxy`, `Vvv = 0`, `Vuv = sinA·DVxy`.
pub fn cone_d2(co: &GpCone, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
    let a = co.semi_angle;
    let big_r = co.radius + v * a.sin();
    let (cu, su) = (u.cos(), u.sin());
    let (p, du, dv) = cone_d1(co, u, v);
    let vxy = co.pos.x_direction().xyz().multiplied(cu)
        .added(&co.pos.y_direction().xyz().multiplied(su));
    let dvxy = co.pos.x_direction().xyz().multiplied(-su)
        .added(&co.pos.y_direction().xyz().multiplied(cu));
    let duu = GpVec::from_xyz(&vxy.multiplied(-big_r));
    let duv = GpVec::from_xyz(&dvxy.multiplied(a.sin()));
    (p, du, dv, duu, GpVec::zero(), duv)
}

/// `ElSLib::SphereD2` (`ElSLib.cxx:975-1037`): `Vuu = -R·cosV·Vxy`,
/// `Vvv = -R·cosV·Vxy - R·sinV·Z`, `Vuv = -R·sinV·DVxy`.
pub fn sphere_d2(s: &GpSphere, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
    let r = s.radius;
    let (cv, sv) = (v.cos(), v.sin());
    let (cu, su) = (u.cos(), u.sin());
    let (p, du, dv) = sphere_d1(s, u, v);
    let vxy = s.pos.x_direction().xyz().multiplied(cu)
        .added(&s.pos.y_direction().xyz().multiplied(su));
    let dvxy = s.pos.x_direction().xyz().multiplied(-su)
        .added(&s.pos.y_direction().xyz().multiplied(cu));
    let duu = GpVec::from_xyz(&vxy.multiplied(-r * cv));
    let dvv = GpVec::from_xyz(&vxy.multiplied(-r * cv)
        .added(&s.pos.direction().xyz().multiplied(-r * sv)));
    let duv = GpVec::from_xyz(&dvxy.multiplied(-r * sv));
    (p, du, dv, duu, dvv, duv)
}

/// `ElSLib::TorusD2` (`ElSLib.cxx:1039-1100`): `Vuu = -R·Vxy`,
/// `Vvv = -r·cosV·Vxy - r·sinV·Z`, `Vuv = r·sinV·(sinU·X - cosU·Y)`.
pub fn torus_d2(t: &GpTorus, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
    let r = t.minor_radius;
    let (cv, sv) = (v.cos(), v.sin());
    let (cu, su) = (u.cos(), u.sin());
    let big_r = t.major_radius + r * cv;
    let (p, du, dv) = torus_d1(t, u, v);
    let vxy = t.pos.x_direction().xyz().multiplied(cu)
        .added(&t.pos.y_direction().xyz().multiplied(su));
    // `Vuv = r·sinV·(sinU·X − cosU·Y) = −r·sinV·DVxy`.
    let dvxy = t.pos.x_direction().xyz().multiplied(-su)
        .added(&t.pos.y_direction().xyz().multiplied(cu));
    let duu = GpVec::from_xyz(&vxy.multiplied(-big_r));
    let dvv = GpVec::from_xyz(&vxy.multiplied(-r * cv)
        .added(&t.pos.direction().xyz().multiplied(-r * sv)));
    let duv = GpVec::from_xyz(&dvxy.multiplied(-r * sv));
    (p, du, dv, duu, dvv, duv)
}

/// Evaluate any elementary surface D2 (point + du + dv + duu + dvv + duv).
///
/// Delegates to the per-type `ElSLib::*D2` transcriptions above. The previous
/// body (before audit A10) returned **zero** `duu/dvv/duv` for sphere and torus,
/// silently dropping curvature; `SurfaceRef` has no `Cone` variant, so a cone is
/// evaluated through [`cone_d2`] directly.
pub fn surface_d2(
    s: &SurfaceRef, u: f64, v: f64,
) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
    match s {
        SurfaceRef::Plane(pl) => plane_d2(pl, u, v),
        SurfaceRef::Cylinder(cy) => cylinder_d2(cy, u, v),
        SurfaceRef::Sphere(s) => sphere_d2(s, u, v),
        SurfaceRef::Torus(t) => torus_d2(t, u, v),
    }
}

/// Reference to an elementary surface for polymorphic evaluation.
#[derive(Debug, Clone, Copy)]
pub enum SurfaceRef<'a> {
    Plane(&'a GpPln),
    Cylinder(&'a GpCylinder),
    Sphere(&'a GpSphere),
    Torus(&'a GpTorus),
}

impl<'a> SurfaceRef<'a> {
    /// Evaluate point.
    pub fn value(&self, u: f64, v: f64) -> GpPnt {
        match self {
            SurfaceRef::Plane(p) => slib::plane_value(p, u, v),
            SurfaceRef::Cylinder(c) => slib::cylinder_value(c, u, v),
            SurfaceRef::Sphere(s) => slib::sphere_value(s, u, v),
            SurfaceRef::Torus(t) => slib::torus_value(t, u, v),
        }
    }
    /// Outward normal.
    pub fn normal(&self, u: f64, v: f64) -> GpVec {
        match self {
            SurfaceRef::Plane(p) => plane_normal(p),
            SurfaceRef::Cylinder(c) => cylinder_normal(c, u),
            SurfaceRef::Sphere(s) => sphere_normal(s, u, v),
            SurfaceRef::Torus(t) => {
                let (_, du, dv) = torus_d1(t, u, v);
                normal_from_d1(&du, &dv)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ElSLib `*DN`: n-th derivative of the elementary surfaces
// Source: `ElSLib.cxx:169-557`.
// ---------------------------------------------------------------------------

/// `IsOdd` / `IsEven` (`Standard_Integer.hxx:33-43`).
fn is_even(v: i32) -> bool {
    v % 2 == 0
}
fn is_odd(v: i32) -> bool {
    v % 2 == 1
}

/// `ElSLib::PlaneDN(U, V, Pos, Nu, Nv)` (`ElSLib.cxx:169-180`).
pub fn plane_dn(_u: f64, _v: f64, pos: &GpAx3, nu: i32, nv: i32) -> GpVec {
    if nu == 0 && nv == 1 {
        GpVec::from_xyz(&pos.y_direction().xyz())
    } else if nu == 1 && nv == 0 {
        GpVec::from_xyz(&pos.x_direction().xyz())
    } else {
        GpVec::new(0.0, 0.0, 0.0)
    }
}

/// `ElSLib::ConeDN(U, V, Pos, Radius, SAngle, Nu, Nv)` (`ElSLib.cxx:182-215`).
pub fn cone_dn(u: f64, v: f64, pos: &GpAx3, radius: f64, s_angle: f64, nu: i32, nv: i32) -> GpVec {
    let xdir = pos.x_direction().xyz();
    let ydir = pos.y_direction().xyz();
    let um = u + nu as f64 * std::f64::consts::FRAC_PI_2;
    let mut xyz = xdir
        .multiply_scalar(um.cos())
        .added(&ydir.multiply_scalar(um.sin()));
    if nv == 0 {
        xyz = xyz.multiply_scalar(radius + v * s_angle.sin());
        if nu == 0 {
            xyz = xyz.added(&pos.location().coord);
        }
        return GpVec::from_xyz(&xyz);
    } else if nv == 1 {
        xyz = xyz.multiply_scalar(s_angle.sin());
        if nu == 0 {
            let zdir = pos.direction();
            xyz = xyz.added(&zdir.xyz().multiply_scalar(s_angle.cos()));
        }
        return GpVec::from_xyz(&xyz);
    }
    GpVec::new(0.0, 0.0, 0.0)
}

/// `ElSLib::CylinderDN(U, V, Pos, Radius, Nu, Nv)` (`ElSLib.cxx:217-265`).
pub fn cylinder_dn(u: f64, _v: f64, pos: &GpAx3, radius: f64, nu: i32, nv: i32) -> GpVec {
    if nu + nv < 1 || nu < 0 || nv < 0 {
        return GpVec::new(0.0, 0.0, 0.0);
    }
    if nv == 0 {
        let r_cos_u = radius * u.cos();
        let r_sin_u = radius * u.sin();
        let mut xdir = pos.x_direction().xyz().clone();
        let mut ydir = pos.y_direction().xyz().clone();
        if (nu + 6) % 4 == 0 {
            xdir = xdir.multiply_scalar(-r_cos_u);
            ydir = ydir.multiply_scalar(-r_sin_u);
        } else if (nu + 5) % 4 == 0 {
            xdir = xdir.multiply_scalar(r_sin_u);
            ydir = ydir.multiply_scalar(-r_cos_u);
        } else if (nu + 3) % 4 == 0 {
            xdir = xdir.multiply_scalar(-r_sin_u);
            ydir = ydir.multiply_scalar(r_cos_u);
        } else if nu % 4 == 0 {
            xdir = xdir.multiply_scalar(r_cos_u);
            ydir = ydir.multiply_scalar(r_sin_u);
        }
        GpVec::from_xyz(&xdir.added(&ydir))
    } else if nv == 1 && nu == 0 {
        GpVec::from_xyz(&pos.direction().xyz())
    } else {
        GpVec::new(0.0, 0.0, 0.0)
    }
}

/// `ElSLib::SphereDN(U, V, Pos, Radius, Nu, Nv)` (`ElSLib.cxx:267-365`).
pub fn sphere_dn(u: f64, v: f64, pos: &GpAx3, radius: f64, nu: i32, nv: i32) -> GpVec {
    if nu + nv < 1 || nu < 0 || nv < 0 {
        return GpVec::new(0.0, 0.0, 0.0);
    }
    let cos_u = u.cos();
    let sin_u = u.sin();
    let r_cos_v = radius * v.cos();
    let xdir = pos.x_direction().xyz();
    let ydir = pos.y_direction().xyz();
    let zdir_ = pos.direction();
    let zdir = zdir_.xyz();
    let mut x;
    let mut y;
    let mut z;
    if nu == 0 {
        let r_sin_v = radius * v.sin();
        let (a1, a2, a3) = if is_odd(nv) {
            (-r_sin_v * cos_u, -r_sin_v * sin_u, r_cos_v)
        } else {
            (-r_cos_v * cos_u, -r_cos_v * sin_u, -r_sin_v)
        };
        x = a1 * xdir.x() + a2 * ydir.x() + a3 * zdir.x();
        y = a1 * xdir.y() + a2 * ydir.y() + a3 * zdir.y();
        z = a1 * xdir.z() + a2 * ydir.z() + a3 * zdir.z();
        if (nv + 2) % 4 != 0 && (nv + 3) % 4 != 0 {
            x = -x;
            y = -y;
            z = -z;
        }
    } else if nv == 0 {
        let (a1, a2) = if is_odd(nu) {
            (-r_cos_v * sin_u, r_cos_v * cos_u)
        } else {
            (r_cos_v * cos_u, r_cos_v * sin_u)
        };
        x = a1 * xdir.x() + a2 * ydir.x();
        y = a1 * xdir.y() + a2 * ydir.y();
        z = a1 * xdir.z() + a2 * ydir.z();
        if (nu + 2) % 4 == 0 || (nu + 1) % 4 == 0 {
            x = -x;
            y = -y;
            z = -z;
        }
    } else {
        let r_sin_v = radius * v.sin();
        let (a1, a2) = if is_odd(nu) {
            (-sin_u, cos_u)
        } else {
            (-cos_u, -sin_u)
        };
        let a3 = if is_odd(nv) { -r_sin_v } else { -r_cos_v };
        // OCCT multiplies the in-plane combination by A3 (`ElSLib.cxx:355-357`);
        // reproduced literally.
        x = (a1 * xdir.x() + a2 * ydir.x()) * a3;
        y = (a1 * xdir.y() + a2 * ydir.y()) * a3;
        z = (a1 * xdir.z() + a2 * ydir.z()) * a3;
        if ((nu + 2) % 4 != 0
            && (nu + 3) % 4 != 0
            && ((nv + 2) % 4 == 0 || (nv + 3) % 4 == 0))
            || (((nu + 2) % 4 == 0 || (nu + 3) % 4 == 0)
                && (nv + 2) % 4 != 0
                && (nv + 3) % 4 != 0)
        {
            x = -x;
            y = -y;
            z = -z;
        }
    }
    GpVec::new(x, y, z)
}

/// `ElSLib::TorusDN(U, V, Pos, MajorRadius, MinorRadius, Nu, Nv)`
/// (`ElSLib.cxx:367-557`).
pub fn torus_dn(
    u: f64,
    v: f64,
    pos: &GpAx3,
    major_radius: f64,
    minor_radius: f64,
    nu: i32,
    nv: i32,
) -> GpVec {
    if nu + nv < 1 || nu < 0 || nv < 0 {
        return GpVec::new(0.0, 0.0, 0.0);
    }
    let cos_u = u.cos();
    let sin_u = u.sin();
    let xdir = pos.x_direction().xyz();
    let ydir = pos.y_direction().xyz();
    let zdir_ = pos.direction();
    let zdir = zdir_.xyz();
    // `eps = 10. * (MinorRadius + MajorRadius) * RealEpsilon()` (`cxx:386`).
    let eps = 10.0 * (minor_radius + major_radius) * f64::EPSILON;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;

    if nv == 0 {
        let r = major_radius + minor_radius * v.cos();
        let (mut a1, mut a2) = if is_odd(nu) {
            (-r * sin_u, r * cos_u)
        } else {
            (-r * cos_u, -r * sin_u)
        };
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        x = a1 * xdir.x() + a2 * ydir.x();
        y = a1 * xdir.y() + a2 * ydir.y();
        z = a1 * xdir.z() + a2 * ydir.z();
        if (nu + 2) % 4 != 0 && (nu + 3) % 4 != 0 {
            x = -x;
            y = -y;
            z = -z;
        }
    } else if nu == 0 {
        let r_cos_v = minor_radius * v.cos();
        let r_sin_v = minor_radius * v.sin();
        let (mut a1, mut a2, mut a3) = if is_odd(nv) {
            (-r_sin_v * cos_u, -r_sin_v * sin_u, r_cos_v)
        } else {
            (-r_cos_v * cos_u, -r_cos_v * sin_u, -r_sin_v)
        };
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        if a3.abs() <= eps {
            a3 = 0.0;
        }
        x = a1 * xdir.x() + a2 * ydir.x() + a3 * zdir.x();
        y = a1 * xdir.y() + a2 * ydir.y() + a3 * zdir.y();
        z = a1 * xdir.z() + a2 * ydir.z() + a3 * zdir.z();
        if (nv + 2) % 4 != 0 && (nv + 3) % 4 != 0 {
            x = -x;
            y = -y;
            z = -z;
        }
    } else if is_odd(nu) && is_odd(nv) {
        let r_sin_v = minor_radius * v.sin();
        let (mut a1, mut a2) = (r_sin_v * sin_u, -r_sin_v * cos_u);
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        x = a1 * xdir.x() + a2 * ydir.x();
        y = a1 * xdir.y() + a2 * ydir.y();
        z = a1 * xdir.z() + a2 * ydir.z();
    } else if is_even(nu) && is_even(nv) {
        let r_cos_v = minor_radius * v.cos();
        let (mut a1, mut a2) = (r_cos_v * cos_u, r_cos_v * sin_u);
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        x = a1 * xdir.x() + a2 * ydir.x();
        y = a1 * xdir.y() + a2 * ydir.y();
        z = a1 * xdir.z() + a2 * ydir.z();
    } else if is_even(nv) && is_odd(nu) {
        let r_cos_v = minor_radius * v.cos();
        let (mut a1, mut a2) = (r_cos_v * sin_u, -r_cos_v * cos_u);
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        x = a1 * xdir.x() + a2 * ydir.x();
        y = a1 * xdir.y() + a2 * ydir.y();
        z = a1 * xdir.z() + a2 * ydir.z();
        if (nv + nu + 3) % 4 == 0 {
            x = -x;
            y = -y;
            z = -z;
        }
    } else if is_odd(nv) && is_even(nu) {
        let r_sin_v = minor_radius * v.sin();
        let (mut a1, mut a2) = (r_sin_v * cos_u, r_sin_v * sin_u);
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        x = a1 * xdir.x() + a2 * ydir.x();
        y = a1 * xdir.y() + a2 * ydir.y();
        z = a1 * xdir.z() + a2 * ydir.z();
        if (nu + nv + 3) % 4 == 0 {
            x = -x;
            y = -y;
            z = -z;
        }
    }
    GpVec::new(x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{GpAx3, GpDir};

    fn ax3_z() -> GpAx3 {
        GpAx3::new(GpPnt::new(0.,0.,0.), GpDir::from_axis(crate::gp::dir::DirAxis::Z), &GpDir::from_axis(crate::gp::dir::DirAxis::X)).unwrap()
    }

    #[test]
    fn plane_normal_constant() {
        let pl = GpPln::new(ax3_z());
        let n = plane_normal(&pl);
        assert!((n.z() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn cylinder_normal_at_angle() {
        let cy = GpCylinder::new(ax3_z(), 2.0).unwrap();
        // At u=0, normal = +X
        let n = cylinder_normal(&cy, 0.0);
        assert!((n.x() - 1.0).abs() < 1e-14);
        // At u=π/2, normal = +Y
        let n2 = cylinder_normal(&cy, std::f64::consts::FRAC_PI_2);
        assert!((n2.y() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn sphere_normal_outward() {
        let s = GpSphere::new(ax3_z(), 3.0).unwrap();
        // At u=0, v=0 → point (3,0,0), normal = +X
        let n = sphere_normal(&s, 0.0, 0.0);
        assert!((n.x() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn surface_ref_value() {
        let s = SurfaceRef::Plane(&GpPln::new(ax3_z()));
        let p = s.value(1.0, 2.0);
        assert!((p.x() - 1.0).abs() < 1e-14);
        assert!((p.y() - 2.0).abs() < 1e-14);
    }
}
