//! Quadric recovery and `IntSurf_Quadric::Normale` for ImpImp.
//! Source: `IntPatch_ImpImpIntersection.cxx` SetQuad / AdjustToSeam,
//! `IntSurf_Quadric.cxx` Normale.

use std::f64::consts::PI;

use occt_core::gp::{
    GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec,
};
use occt_core::precision::ANGULAR;
use occt_geom::Surface;

use crate::pcurve_full::{cone_params, torus_params};

use super::super::plane_from_surface;
use super::super::sphere_params;

/// Plane=1, Cylinder=2, Cone=3, Sphere=4, Torus=5 (`SetQuad` iRet).
#[derive(Debug, Clone)]
pub(crate) enum ImplicitQuad {
    Plane(GpPln),
    Cylinder(GpCylinder),
    Cone(GpCone),
    Sphere(GpSphere),
    Torus(GpTorus),
}

impl ImplicitQuad {
    pub(crate) fn code(&self) -> i32 {
        match self {
            ImplicitQuad::Plane(_) => 1,
            ImplicitQuad::Cylinder(_) => 2,
            ImplicitQuad::Cone(_) => 3,
            ImplicitQuad::Sphere(_) => 4,
            ImplicitQuad::Torus(_) => 5,
        }
    }

    pub(crate) fn as_plane(&self) -> Option<&GpPln> {
        match self {
            ImplicitQuad::Plane(p) => Some(p),
            _ => None,
        }
    }

    pub(crate) fn as_cylinder(&self) -> Option<&GpCylinder> {
        match self {
            ImplicitQuad::Cylinder(c) => Some(c),
            _ => None,
        }
    }

    pub(crate) fn as_cone(&self) -> Option<&GpCone> {
        match self {
            ImplicitQuad::Cone(c) => Some(c),
            _ => None,
        }
    }

    pub(crate) fn as_sphere(&self) -> Option<&GpSphere> {
        match self {
            ImplicitQuad::Sphere(s) => Some(s),
            _ => None,
        }
    }

    pub(crate) fn as_torus(&self) -> Option<&GpTorus> {
        match self {
            ImplicitQuad::Torus(t) => Some(t),
            _ => None,
        }
    }

    /// `IntSurf_Quadric::Parameters` / `ElSLib::Parameters`.
    pub(crate) fn parameters(&self, p: &GpPnt) -> (f64, f64) {
        use occt_core::elib::slib;
        match self {
            ImplicitQuad::Plane(pl) => slib::plane_parameters(&pl.position(), p),
            ImplicitQuad::Cylinder(cy) => slib::cylinder_parameters(&cy.position(), p),
            ImplicitQuad::Cone(co) => {
                slib::cone_parameters(&co.position(), co.radius(), co.semi_angle(), p)
            }
            ImplicitQuad::Sphere(sp) => slib::sphere_parameters(sp.position(), p),
            ImplicitQuad::Torus(to) => {
                slib::torus_parameters(to.position(), to.major_radius(), to.minor_radius(), p)
            }
        }
    }
}

/// `SetQuad`. Returns None when the surface is not an implicit quadric.
pub(crate) fn set_quad(s: &dyn Surface) -> Option<ImplicitQuad> {
    if let Some(p) = plane_from_surface(s) {
        return Some(ImplicitQuad::Plane(p));
    }
    if let Some(c) = cylinder_from_surface(s) {
        return Some(ImplicitQuad::Cylinder(c));
    }
    if let Some(c) = cone_from_surface(s) {
        return Some(ImplicitQuad::Cone(c));
    }
    if let Some(sp) = sphere_from_surface(s) {
        return Some(ImplicitQuad::Sphere(sp));
    }
    if let Some(t) = torus_from_surface(s) {
        return Some(ImplicitQuad::Torus(t));
    }
    None
}

fn perp_x_dir(z: &GpDir) -> GpDir {
    let base = if z.x().abs() < 0.9 {
        GpDir::new(1.0, 0.0, 0.0).expect("x")
    } else {
        GpDir::new(0.0, 1.0, 0.0).expect("y")
    };
    base.cross(z).unwrap_or(base)
}

fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let axis = GpVec::from_pnts(&p0, &p1);
    let m = axis.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axis.divided(m);
    let q0 = s.d0(0.0, 0.0);
    let q1 = s.d0(PI, 0.0);
    let r = q0.distance(&q1) / 2.0;
    if r <= 1e-9 {
        return None;
    }
    let center = GpPnt::new(
        0.5 * (q0.x() + q1.x()),
        0.5 * (q0.y() + q1.y()),
        0.5 * (q0.z() + q1.z()),
    );
    let (nu, nv) = (8, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = (j as f64 - nv as f64 / 2.0) / 2.0;
            let p = s.d0(u, v);
            if (dist_to_axis(&p, &center, &ax) - r).abs() > 1e-3 * r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, r))
}

fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

fn cylinder_from_surface(s: &dyn Surface) -> Option<GpCylinder> {
    let (center, ax, r) = cylinder_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(center, z, &x).ok()?;
    GpCylinder::new(ax3, r).ok()
}

fn cone_from_surface(s: &dyn Surface) -> Option<GpCone> {
    let (apex, ax, alpha) = cone_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(apex, z, &x).ok()?;
    GpCone::new(ax3, 0.0, alpha).ok()
}

fn sphere_from_surface(s: &dyn Surface) -> Option<GpSphere> {
    let (c, r) = sphere_params(s)?;
    let z = GpDir::new(0.0, 0.0, 1.0).ok()?;
    let x = GpDir::new(1.0, 0.0, 0.0).ok()?;
    let ax3 = GpAx3::new(c, z, &x).ok()?;
    GpSphere::new(ax3, r).ok()
}

fn torus_from_surface(s: &dyn Surface) -> Option<GpTorus> {
    let (center, ax, major, minor) = torus_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(center, z, &x).ok()?;
    GpTorus::new(ax3, major, minor).ok()
}

/// `IntSurf_QuadricTool::Tolerance`.
pub(crate) fn quadric_tolerance(q: &ImplicitQuad) -> f64 {
    match q {
        ImplicitQuad::Sphere(sp) => 2.0e-6 * sp.radius(),
        ImplicitQuad::Cylinder(cy) => 2.0e-6 * cy.radius(),
        _ => 1.0e-6,
    }
}

/// `IntSurf_Quadric::Gradient(P)`.
pub(crate) fn gradient(q: &ImplicitQuad, p: &GpPnt) -> GpVec {
    match q {
        ImplicitQuad::Plane(pl) => GpVec::from_xyz(pl.axis().direction().xyz()),
        ImplicitQuad::Cylinder(cy) => {
            let loc = cy.location();
            let az = GpVec::from_xyz(cy.axis().direction().xyz());
            let w = GpVec::from_pnts(&loc, p);
            let pp = loc.translated_vec(&az.multiplied_scalar(w.dot(&az)));
            let g = GpVec::from_pnts(&pp, p);
            let n = g.magnitude();
            if n > 1e-14 {
                g.divided(n)
            } else {
                GpVec::zero()
            }
        }
        ImplicitQuad::Sphere(sp) => {
            let g = GpVec::from_pnts(&sp.location(), p);
            let n = g.magnitude();
            if n > 1e-14 {
                g.divided(n)
            } else {
                GpVec::zero()
            }
        }
        ImplicitQuad::Cone(_) | ImplicitQuad::Torus(_) => {
            let g = normale(q, p);
            let n = g.magnitude();
            if n > 1e-14 {
                g.divided(n)
            } else {
                GpVec::zero()
            }
        }
    }
}

/// `IntSurf_Quadric::Normale(P)`.
pub(crate) fn normale(q: &ImplicitQuad, p: &GpPnt) -> GpVec {
    match q {
        ImplicitQuad::Plane(pl) => GpVec::from_xyz(pl.axis().direction().xyz()),
        ImplicitQuad::Cylinder(cy) => {
            let lin = GpLin::from_pnt_dir(cy.location(), *cy.axis().direction());
            GpVec::from_xyz(lin.normal(p).direction().xyz())
        }
        ImplicitQuad::Sphere(sp) => GpVec::from_pnts(&sp.location(), p),
        ImplicitQuad::Cone(co) => {
            let apex = co.apex();
            let gen = GpVec::from_pnts(&apex, p);
            if gen.magnitude() < 1e-14 {
                return GpVec::zero();
            }
            let az = GpVec::from_xyz(co.axis().direction().xyz());
            let circ = az.crossed(&gen);
            if circ.magnitude() < 1e-14 {
                return GpVec::zero();
            }
            circ.crossed(&gen)
        }
        ImplicitQuad::Torus(to) => {
            let o = to.location();
            let oz = GpVec::from_xyz(to.position().direction().xyz());
            let op = GpVec::from_pnts(&o, p);
            let pp = p.translated_vec(&oz.multiplied_scalar(-op.dot(&oz)));
            let dopp = if o.square_distance(&pp) < 1e-14 {
                GpVec::from_xyz(to.position().x_direction().xyz())
            } else {
                GpVec::from_pnts(&o, &pp)
            };
            let Ok(dopp) = GpDir::from_vec(&dopp) else {
                return GpVec::zero();
            };
            let pt = o.translated_vec(&GpVec::from_xyz(dopp.xyz()).multiplied_scalar(to.major_radius()));
            if pt.square_distance(p) < 1e-14 {
                return oz;
            }
            GpVec::from_pnts(&pt, p)
        }
    }
}

/// `SeamPosition` + `AdjustToSeam` for a circle on a quadric of revolution.
pub(crate) fn adjust_circ_to_seam(circ: &mut GpCirc, pos: &GpAx3) {
    if let Ok(ax2) = GpAx2::new(circ.location(), pos.direction(), *pos.x_direction()) {
        circ.set_position(&ax2);
    }
}

pub(crate) fn adjust_sphere_circ(circ: &mut GpCirc, sph: &GpSphere) {
    let qdir = sph.position().direction();
    if circ.axis().direction().is_parallel(&qdir) {
        adjust_circ_to_seam(circ, sph.position());
    }
}

/// Axis-to-axis distance when directions are parallel; otherwise the common
/// perpendicular length.
pub(crate) fn axis_distance(a: &occt_core::gp::GpAx1, b: &occt_core::gp::GpAx1) -> f64 {
    let d1 = GpVec::from_xyz(a.direction().xyz());
    let d2 = GpVec::from_xyz(b.direction().xyz());
    let n = d1.crossed(&d2);
    let nm = n.magnitude();
    let w = GpVec::from_pnts(a.location(), b.location());
    if nm < ANGULAR {
        w.subtracted(&d1.multiplied_scalar(w.dot(&d1))).magnitude()
    } else {
        w.dot(&n).abs() / nm
    }
}

pub(crate) fn dist_point_axis(p: &GpPnt, ax: &occt_core::gp::GpAx1) -> f64 {
    GpLin::from_pnt_dir(*ax.location(), *ax.direction()).distance(p)
}

/// `IntSurf_Quadric::Distance(P)` — signed distance used by `IntPatch_ArcFunction`.
pub(crate) fn distance(q: &ImplicitQuad, p: &GpPnt) -> f64 {
    match q {
        ImplicitQuad::Plane(pl) => {
            let n = GpVec::from_xyz(pl.axis().direction().xyz());
            let loc = pl.location();
            n.x() * (p.x() - loc.x()) + n.y() * (p.y() - loc.y()) + n.z() * (p.z() - loc.z())
        }
        ImplicitQuad::Cylinder(cy) => {
            let lin = GpLin::from_pnt_dir(cy.location(), *cy.axis().direction());
            lin.distance(p) - cy.radius()
        }
        ImplicitQuad::Sphere(sp) => p.distance(&sp.location()) - sp.radius(),
        ImplicitQuad::Cone(co) => {
            let lin = GpLin::from_pnt_dir(*co.axis().location(), *co.axis().direction());
            let dist = lin.distance(p);
            let az = GpVec::from_xyz(co.axis().direction().xyz());
            let vpar = GpVec::from_pnts(&co.apex(), p).dot(&az);
            let r_at = (co.radius() + vpar * co.semi_angle().tan()).abs();
            let cosa = co.semi_angle().cos();
            if cosa.abs() < 1e-14 {
                dist - r_at
            } else {
                (dist - r_at) / cosa
            }
        }
        ImplicitQuad::Torus(to) => {
            let o = to.location();
            let oz = GpVec::from_xyz(to.position().direction().xyz());
            let op = GpVec::from_pnts(&o, p);
            let pp = p.translated_vec(&oz.multiplied_scalar(-op.dot(&oz)));
            let dopp = if o.square_distance(&pp) < 1e-14 {
                GpVec::from_xyz(to.position().x_direction().xyz())
            } else {
                GpVec::from_pnts(&o, &pp)
            };
            let Ok(dopp) = GpDir::from_vec(&dopp) else {
                return 0.0;
            };
            let pt = o.translated_vec(
                &GpVec::from_xyz(dopp.xyz()).multiplied_scalar(to.major_radius()),
            );
            p.distance(&pt) - to.minor_radius()
        }
    }
}
