//! `HLRBRep_SurfaceTool` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_SurfaceTool.hxx`,
//! `.cxx:22-155`).
//!
//! The `SurfaceTool` of `HLRBRep_InterCSurf`'s `gp_Lin` instantiation
//! (`HLRBRep_InterCSurf.cxx:47-48`), i.e. the surface half of the
//! `IntCurveSurface_Inter` template. Every member is a call on
//! [`HlrBRepSurface`], which is a `BRepAdaptor_Surface` plus the cached
//! `GetType()`, so the sampling counts and the UV window follow the **face**,
//! not the underlying `Geom_Surface`.
//!
//! Two numbers decide the HLR surface grid and neither is
//! `Adaptor3d_HSurfaceTool`'s:
//!
//! * `NbSamplesU(S, U1, U2)` / `NbSamplesV(S, V1, V2)` (`cxx:113-155`) divide by
//!   `(uf - ul)` (first minus last, i.e. a **negative** number) and truncate the
//!   ratio to an integer *before* multiplying (`n *= (int)(...)`), where
//!   `Adaptor3d_HSurfaceTool` (`...cxx:72-93`, `:95-116`) divides by the
//!   positive `(ul - uf)`, has no truncation step and clamps to `50`. For any
//!   sub-window the HLR ratio is a negative integer, so `n` goes negative and
//!   the `n < 5` floor answers `5`; the port reproduces that literally rather
//!   than reusing [`occt_topo::intcurvesurface::surface_nb_samples_u_range`].
//! * the per-type base counts (`cxx:22-108`) switch on
//!   `HLRBRep_Surface::GetType()` (the cached type), so a Bezier patch of degree
//!   1 in U and V is a **plane** here (hence `2` samples) while the adaptor
//!   reports `GeomAbs_BezierSurface` (`3 + NbUPoles`).
//!
//! `UTrim` (`IntCurveSurface_Inter.pxx:489-490`) is the only member of the
//! `IntCurveSurface` `SurfaceTool` protocol with no `HLRBRep_SurfaceTool`
//! counterpart: `HLRBRep_SurfaceTool` declares no trim, and the one place
//! `Inter.pxx` calls it (the B-spline arm of `InternalPerformPolygonBounds`) has
//! no HLR equivalent, so nothing is ported for it.

use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::hlrbrep_surface::HlrBRepSurface;

/// `GetType()` (`HLRBRep_Surface.lxx:191-194`): the cached classified type.
pub fn get_type(s: &HlrBRepSurface) -> GeomAbsSurfaceType {
    s.get_type()
}

/// `FirstUParameter()` (`HLRBRep_Surface.lxx:37-40`): the face's UMin.
pub fn first_u_parameter(s: &HlrBRepSurface) -> f64 {
    s.first_u_parameter()
}

/// `LastUParameter()` (`HLRBRep_Surface.lxx:44-47`): the face's UMax.
pub fn last_u_parameter(s: &HlrBRepSurface) -> f64 {
    s.last_u_parameter()
}

/// `FirstVParameter()` (`HLRBRep_Surface.lxx:51-54`): the face's VMin.
pub fn first_v_parameter(s: &HlrBRepSurface) -> f64 {
    s.first_v_parameter()
}

/// `LastVParameter()` (`HLRBRep_Surface.lxx:58-61`): the face's VMax.
pub fn last_v_parameter(s: &HlrBRepSurface) -> f64 {
    s.last_v_parameter()
}

/// `NbSamplesU(S)` (`cxx:22-65`): the per-type base count.
pub fn nb_samples_u(s: &HlrBRepSurface) -> i64 {
    match get_type(s) {
        GeomAbsSurfaceType::Plane => 2,
        GeomAbsSurfaceType::BezierSurface => 3 + s.nb_u_poles().unwrap_or(0) as i64,
        GeomAbsSurfaceType::BSplineSurface => {
            let nbs = s.nb_u_knots().unwrap_or(0) as i64 * s.u_degree().unwrap_or(0) as i64;
            if nbs < 2 {
                2
            } else {
                nbs
            }
        }
        GeomAbsSurfaceType::Torus => 20,
        // `cxx:49-57`: Cylinder / Cone / Sphere / SurfaceOfRevolution /
        // SurfaceOfExtrusion.
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::SurfaceOfRevolution
        | GeomAbsSurfaceType::SurfaceOfExtrusion => 10,
        // `cxx:58-61`.
        _ => 10,
    }
}

/// `NbSamplesV(S)` (`cxx:68-108`): the per-type base count.
pub fn nb_samples_v(s: &HlrBRepSurface) -> i64 {
    match get_type(s) {
        GeomAbsSurfaceType::Plane => 2,
        GeomAbsSurfaceType::BezierSurface => 3 + s.nb_v_poles().unwrap_or(0) as i64,
        GeomAbsSurfaceType::BSplineSurface => {
            let nbs = s.nb_v_knots().unwrap_or(0) as i64 * s.v_degree().unwrap_or(0) as i64;
            if nbs < 2 {
                2
            } else {
                nbs
            }
        }
        // `cxx:91-99`.
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::Torus
        | GeomAbsSurfaceType::SurfaceOfRevolution
        | GeomAbsSurfaceType::SurfaceOfExtrusion => 15,
        // `cxx:101-104`.
        _ => 10,
    }
}

/// `NbSamplesU(S, U1, U2)` (`cxx:113-133`): `n *= (int)((U2-U1)/(Uf-Ul))`, then
/// the `n > nbs` clamp and the `n < 5` floor.
pub fn nb_samples_u_range(s: &HlrBRepSurface, u1: f64, u2: f64) -> i64 {
    let nbs = nb_samples_u(s);
    let mut n = nbs;
    if nbs > 10 {
        let uf = first_u_parameter(s);
        let ul = last_u_parameter(s);
        // `(int)` truncation of the ratio, then the multiplication.
        let ratio = ((u2 - u1) / (uf - ul)).trunc();
        n *= ratio as i64;
        if n > nbs {
            n = nbs;
        }
        if n < 5 {
            n = 5;
        }
    }
    n
}

/// `NbSamplesV(S, V1, V2)` (`cxx:136-155`).
pub fn nb_samples_v_range(s: &HlrBRepSurface, v1: f64, v2: f64) -> i64 {
    let nbs = nb_samples_v(s);
    let mut n = nbs;
    if nbs > 10 {
        let vf = first_v_parameter(s);
        let vl = last_v_parameter(s);
        let ratio = ((v2 - v1) / (vf - vl)).trunc();
        n *= ratio as i64;
        if n > nbs {
            n = nbs;
        }
        if n < 5 {
            n = 5;
        }
    }
    n
}

/// `IsUPeriodic()` (`lxx:104-107`).
pub fn is_u_periodic(s: &HlrBRepSurface) -> bool {
    s.is_u_periodic()
}

/// `IsVPeriodic()` (`lxx:118-121`).
pub fn is_v_periodic(s: &HlrBRepSurface) -> bool {
    s.is_v_periodic()
}

/// `UPeriod()` (`lxx:111-114`).
pub fn u_period(s: &HlrBRepSurface) -> f64 {
    s.u_period()
}

/// `VPeriod()` (`lxx:125-128`).
pub fn v_period(s: &HlrBRepSurface) -> f64 {
    s.v_period()
}

/// `Plane()` (`HLRBRep_Surface.cxx:319-335`): the tangent plane at `(0.5, 0.5)`
/// for a Bezier patch, `GeomAdaptor_Surface::Plane` otherwise.
pub fn plane(s: &HlrBRepSurface) -> Option<occt_core::gp::GpPln> {
    s.plane()
}

/// `Cylinder()` (`lxx:198-201`).
pub fn cylinder(s: &HlrBRepSurface) -> Option<occt_core::gp::GpCylinder> {
    s.cylinder()
}

/// `Cone()` (`lxx:205-208`).
pub fn cone(s: &HlrBRepSurface) -> Option<occt_core::gp::GpCone> {
    s.cone()
}

/// `Sphere()` (`lxx:212-215`).
pub fn sphere(s: &HlrBRepSurface) -> Option<occt_core::gp::GpSphere> {
    s.sphere()
}

/// `Torus()` (`lxx:219-222`).
pub fn torus(s: &HlrBRepSurface) -> Option<occt_core::gp::GpTorus> {
    s.torus()
}
