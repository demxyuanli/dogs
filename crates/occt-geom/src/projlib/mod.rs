//! `ProjLib` — project 3D points and elementary curves onto analytic surfaces.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib` (`ProjLib.cxx`). Point
//! projection is `ElSLib::Parameters`. [`projected_samples`] is the sample
//! half of `ProjLib_ProjectedCurve` (periodic-seam unwrap); the 2D B-spline
//! fit lives with pcurve construction in `occt-topo`.
//!
//! The `ProjLib_Projector` family (the vtable `ProjLib_ProjectedCurve` talks
//! to) lives in the sibling modules of this file, one OCCT class each:
//! [`projector`] (`ProjLib_Projector.cxx:34-251`), [`plane`]
//! (`ProjLib_Plane.cxx:24-169`), [`cylinder`] (`ProjLib_Cylinder.cxx:28-175`),
//! [`cone`] (`ProjLib_Cone.cxx:30-180`), [`sphere`]
//! (`ProjLib_Sphere.cxx:33-249`) and [`torus`] (`ProjLib_Torus.cxx:28-193`).

mod cone;
mod comp_projected_curve;
mod compute_approx;
mod cylinder;
mod plane;
mod prj_func;
mod prj_resolve;
mod projected_curve;
mod projected_curve_general;
mod projector;
mod sphere;
mod torus;

pub use cone::Cone;
pub use comp_projected_curve::CompProjectedCurve;
pub use comp_projected_curve::ProjectedCurve2d;
pub use compute_approx::{ComputeApprox, ProjLibFunction, UvBounds};
pub use cylinder::Cylinder;
pub use plane::Plane;
pub use prj_func::PrjFunc;
pub use prj_resolve::PrjResolve;
pub use projected_curve::{
    compute_tol_u, compute_tol_v, curve_type, iso_is_deg, project, projector_value, surface_type,
    trim_c3d, Constraint, ProjectedCurve,
};
pub use projector::{Projector, ProjectorState};
pub use sphere::Sphere;
pub use torus::Torus;

use occt_core::elib::slib;
use occt_core::gp::{
    GpCirc, GpCone, GpCylinder, GpDir2d, GpLin, GpLin2d, GpPln, GpPnt, GpPnt2d, GpSphere, GpTorus,
};

use crate::curve::Curve;
use crate::surface::Surface;

/// `ProjLib::IsAnaSurf(Adaptor3d_Surface&)` (`ProjLib.cxx`): `true` for the
/// five elementary surfaces only — a Bezier / BSpline / swept / offset surface
/// answers `false` even though it is "analytical" in the loose sense.
pub fn is_ana_surf(s: &dyn Surface) -> bool {
    matches!(
        projected_curve::surface_type(s),
        occt_core::kernel::geomabs::SurfaceType::Plane
            | occt_core::kernel::geomabs::SurfaceType::Cylinder
            | occt_core::kernel::geomabs::SurfaceType::Cone
            | occt_core::kernel::geomabs::SurfaceType::Sphere
            | occt_core::kernel::geomabs::SurfaceType::Torus
    )
}

/// `ProjLib::Project(gp_Pln, gp_Pnt)`.
pub fn project_pln_pnt(pl: &GpPln, p: &GpPnt) -> GpPnt2d {
    let (u, v) = slib::plane_parameters(&pl.pos, p);
    GpPnt2d::new(u, v)
}

/// `ProjLib_ProjectOnPlane::ProjectPnt` (`cxx:486-497`) with `Dir = plane Z`.
pub fn project_pnt_on_plane(pl: &GpPln, p: &GpPnt) -> GpPnt {
    let loc = pl.location();
    let z = pl.pos.direction();
    let alpha = (loc.x() - p.x()) * z.x() + (loc.y() - p.y()) * z.y() + (loc.z() - p.z()) * z.z();
    GpPnt::new(p.x() + alpha * z.x(), p.y() + alpha * z.y(), p.z() + alpha * z.z())
}

/// `ProjLib_Plane::EvalPnt2d` (`cxx:87-92`).
pub fn eval_pln_pnt2d(pl: &GpPln, p: &GpPnt) -> GpPnt2d {
    let loc = pl.location();
    let x = pl.pos.x_direction();
    let y = pl.pos.y_direction();
    let ox = p.x() - loc.x();
    let oy = p.y() - loc.y();
    let oz = p.z() - loc.z();
    GpPnt2d::new(
        ox * x.x() + oy * x.y() + oz * x.z(),
        ox * y.x() + oy * y.y() + oz * y.z(),
    )
}

/// `ProjLib_Plane::EvalDir2d` (`cxx:94-97`).
pub fn eval_pln_dir2d(pl: &GpPln, d: &occt_core::gp::GpDir) -> (f64, f64) {
    let x = pl.pos.x_direction();
    let y = pl.pos.y_direction();
    (
        d.x() * x.x() + d.y() * x.y() + d.z() * x.z(),
        d.x() * y.x() + d.y() * y.y() + d.z() * y.z(),
    )
}

/// `ProjLib::Project(gp_Cylinder, gp_Pnt)`.
pub fn project_cyl_pnt(cy: &GpCylinder, p: &GpPnt) -> GpPnt2d {
    let (u, v) = slib::cylinder_parameters(&cy.pos, p);
    GpPnt2d::new(u, v)
}

/// `ProjLib::Project(gp_Cone, gp_Pnt)`.
pub fn project_cone_pnt(co: &GpCone, p: &GpPnt) -> GpPnt2d {
    let (u, v) = slib::cone_parameters(&co.pos, co.radius, co.semi_angle, p);
    GpPnt2d::new(u, v)
}

/// `ProjLib::Project(gp_Sphere, gp_Pnt)`.
pub fn project_sphere_pnt(sp: &GpSphere, p: &GpPnt) -> GpPnt2d {
    let (u, v) = slib::sphere_parameters(&sp.pos, p);
    GpPnt2d::new(u, v)
}

/// `ProjLib::Project(gp_Torus, gp_Pnt)`.
pub fn project_torus_pnt(to: &GpTorus, p: &GpPnt) -> GpPnt2d {
    let (u, v) = slib::torus_parameters(&to.pos, to.major_radius, to.minor_radius, p);
    GpPnt2d::new(u, v)
}

/// Project the location and a unit step of a 3D line onto a plane
/// (`ProjLib_Plane` line case).
pub fn project_pln_lin(pl: &GpPln, lin: &GpLin) -> (GpPnt2d, GpPnt2d) {
    let p0 = lin.location();
    let dir = lin.direction();
    let d = dir.xyz();
    let p1 = GpPnt::new(p0.x() + d.x, p0.y() + d.y, p0.z() + d.z);
    (project_pln_pnt(pl, &p0), project_pln_pnt(pl, &p1))
}

/// Project a circle's center onto a plane (`ProjLib_Plane` circle case).
pub fn project_pln_circ(pl: &GpPln, c: &GpCirc) -> GpPnt2d {
    project_pln_pnt(pl, &c.location())
}

/// `ProjLib_Plane::Project(const gp_Lin&)` (`ProjLib_Plane.cxx:102-107`, via
/// `ProjLib::Project(const gp_Pln&, const gp_Lin&)`, `ProjLib.cxx:59-63`): the
/// line's location and direction carried into the plane's local 2D frame.
/// `None` where OCCT's `gp_Dir2d(X, Y)` raises `Standard_ConstructionError`,
/// i.e. when the direction lies along the plane normal.
pub fn project_pln_lin2d(pl: &GpPln, lin: &GpLin) -> Option<GpLin2d> {
    let dir3 = lin.direction();
    let (dx, dy) = eval_pln_dir2d(pl, &dir3);
    let dir = GpDir2d::new(dx, dy).ok()?;
    Some(GpLin2d::from_pnt_dir(eval_pln_pnt2d(pl, &lin.location()), dir))
}

/// Sample `curve` on `[first, last]`, project each point, unwrap periodic
/// seams (`ProjLib_ProjectedCurve` sampling).
pub fn projected_samples(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    project: impl Fn(&GpPnt) -> GpPnt2d,
    u_periodic: bool,
    v_periodic: bool,
) -> Vec<GpPnt2d> {
    let n = 32usize;
    if !(last > first) {
        return Vec::new();
    }
    let mut pts: Vec<GpPnt2d> = Vec::with_capacity(n + 1);
    let mut prev_u = f64::NAN;
    let mut prev_v = f64::NAN;
    for i in 0..=n {
        let t = first + (last - first) * (i as f64) / (n as f64);
        let q = project(&curve.d0(t));
        let mut u = q.x();
        let mut v = q.y();
        if u_periodic && prev_u.is_finite() {
            while u - prev_u > std::f64::consts::PI {
                u -= 2.0 * std::f64::consts::PI;
            }
            while u - prev_u < -std::f64::consts::PI {
                u += 2.0 * std::f64::consts::PI;
            }
        }
        if v_periodic && prev_v.is_finite() {
            while v - prev_v > std::f64::consts::PI {
                v -= 2.0 * std::f64::consts::PI;
            }
            while v - prev_v < -std::f64::consts::PI {
                v += 2.0 * std::f64::consts::PI;
            }
        }
        prev_u = u;
        prev_v = v;
        pts.push(GpPnt2d::new(u, v));
    }
    pts
}
