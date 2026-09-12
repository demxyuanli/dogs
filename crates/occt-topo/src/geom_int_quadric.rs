//! `IntSurf_Quadric::Parameters` / `ElSLib::Parameters` used by LineConstructor.

use std::f64::consts::PI;

use occt_core::gp::{GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};

/// UV of `pt` on `s`. Analytic kinds use ElSLib-style frames; others project.
pub fn surface_parameters(s: &dyn Surface, pt: &GpPnt) -> Option<(f64, f64)> {
    match classify_surface(s) {
        SurfaceKind::Plane => plane_from_samples(s, pt),
        SurfaceKind::Cylinder => cylinder_from_samples(s, pt),
        SurfaceKind::Sphere => sphere_from_samples(s, pt),
        SurfaceKind::Cone | SurfaceKind::Torus | SurfaceKind::Other => {
            project_point_on_surface(s, pt, 0.0).map(|p| (p.u, p.v))
        }
    }
}

fn plane_from_samples(s: &dyn Surface, pt: &GpPnt) -> Option<(f64, f64)> {
    let o = s.d0(0.0, 0.0);
    let ( _, du, dv) = s.d1(0.0, 0.0);
    let nx = du.magnitude();
    let ny = dv.magnitude();
    if nx < 1e-30 || ny < 1e-30 {
        return project_point_on_surface(s, pt, 0.0).map(|p| (p.u, p.v));
    }
    let xu = du.divided(nx);
    let yv = dv.divided(ny);
    let w = GpVec::from_pnts(&o, pt);
    Some((w.dot(&xu), w.dot(&yv)))
}

fn cylinder_from_samples(s: &dyn Surface, pt: &GpPnt) -> Option<(f64, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let axis = GpVec::from_pnts(&p0, &p1);
    let m = axis.magnitude();
    if m < 1e-30 {
        return project_point_on_surface(s, pt, 0.0).map(|p| (p.u, p.v));
    }
    let az = axis.divided(m);
    let q0 = s.d0(0.0, 0.0);
    let q1 = s.d0(PI, 0.0);
    let center = occt_core::gp::GpPnt::new(
        0.5 * (q0.x() + q1.x()),
        0.5 * (q0.y() + q1.y()),
        0.5 * (q0.z() + q1.z()),
    );
    let px = s.d0(0.0, 0.0);
    let rel0 = GpVec::from_pnts(&center, &px);
    let xdir = rel0.subtracted(&az.multiplied_scalar(rel0.dot(&az)));
    let xn = xdir.magnitude();
    if xn < 1e-30 {
        return project_point_on_surface(s, pt, 0.0).map(|p| (p.u, p.v));
    }
    let x = xdir.divided(xn);
    let y = az.crossed(&x);
    let w = GpVec::from_pnts(&center, pt);
    let v = w.dot(&az);
    let rx = w.dot(&x);
    let ry = w.dot(&y);
    let mut u = ry.atan2(rx);
    if u < 0.0 {
        u += 2.0 * PI;
    }
    Some((u, v))
}

fn sphere_from_samples(s: &dyn Surface, pt: &GpPnt) -> Option<(f64, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(PI, 0.0);
    let center = occt_core::gp::GpPnt::new(
        0.5 * (p0.x() + p1.x()),
        0.5 * (p0.y() + p1.y()),
        0.5 * (p0.z() + p1.z()),
    );
    let r = center.distance(&p0);
    if r < 1e-30 {
        return project_point_on_surface(s, pt, 0.0).map(|p| (p.u, p.v));
    }
    let (_, du, dv) = s.d1(0.0, 0.0);
    let x = if du.magnitude() > 1e-30 {
        du.normalized()
    } else {
        GpVec::new(1.0, 0.0, 0.0)
    };
    let z = if dv.magnitude() > 1e-30 {
        dv.normalized()
    } else {
        GpVec::new(0.0, 0.0, 1.0)
    };
    let y = z.crossed(&x);
    let w = GpVec::from_pnts(&center, pt);
    let wx = w.dot(&x);
    let wy = w.dot(&y);
    let wz = w.dot(&z);
    let mut u = wy.atan2(wx);
    if u < 0.0 {
        u += 2.0 * PI;
    }
    let v = (wz / r).clamp(-1.0, 1.0).asin();
    Some((u, v))
}

/// `GeomInt_LineConstructor` AdjustPeriodic on both surfaces.
pub fn adjust_periodic_uv(
    s1: &dyn Surface,
    s2: &dyn Surface,
    u1: &mut f64,
    v1: &mut f64,
    u2: &mut f64,
    v2: &mut f64,
) {
    adjust_one(s1, u1, v1);
    adjust_one(s2, u2, v2);
}

fn adjust_one(s: &dyn Surface, u: &mut f64, v: &mut f64) {
    let kind = classify_surface(s);
    let (u_per, v_per) = match kind {
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere => (true, false),
        SurfaceKind::Torus => (true, true),
        _ => (false, false),
    };
    let two_pi = PI + PI;
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    if u_per {
        let (nu, _) = crate::geom_int::adjust_periodic(*u, uf, ul, two_pi, 0.0);
        *u = nu;
    }
    if v_per {
        let (nv, _) = crate::geom_int::adjust_periodic(*v, vf, vl, two_pi, 0.0);
        *v = nv;
    }
}

/// Keep gp constructors available for GLine evaluation helpers.
#[allow(dead_code)]
pub fn unused_gp_markers(_p: &GpPln, _c: &GpCylinder, _k: &GpCone, _s: &GpSphere, _t: &GpTorus) {}
