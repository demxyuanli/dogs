//! Elementary surface evaluation with full derivatives and normals.
//! Source: `ElSLib.cxx` D1/D2/Norm functions.
use crate::gp::{GpPln, GpCylinder, GpCone, GpSphere, GpTorus, GpPnt, GpVec, GpXyz};
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
