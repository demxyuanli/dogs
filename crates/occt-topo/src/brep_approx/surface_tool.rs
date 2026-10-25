//! `BRepApprox_SurfaceTool` — the `ThePSurfaceTool` of the `BRepApprox`
//! instantiation: static functions over `BRepAdaptor_Surface`.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepApprox/BRepApprox_SurfaceTool.hxx:38-168`,
//! `.lxx:22-251`, `.cxx:19-40` (OCCT 8.0.0). Every `.lxx` body is a one-line
//! forward to the adaptor, so the functions here forward to
//! [`SurfaceAdaptor`](crate::topopebrep::tool::adaptor::SurfaceAdaptor) the same
//! way.
//!
//! Ported: `FirstUParameter` / `FirstVParameter` / `LastUParameter` /
//! `LastVParameter` (`lxx:22-44`), `Value` (`lxx:143-149`), `D0` (`lxx:151-157`),
//! `D1` (`lxx:159-168`), `GetType` (`lxx:222-226`), the four `NbSamples*`
//! (`cxx:19-40`), `IsUPeriodic` / `UPeriod` / `IsVPeriodic` / `VPeriod`
//! (`lxx:129-141`), `UResolution` / `VResolution` (`lxx:206-220`) and
//! `Plane` / `Cylinder` / `Sphere` / `Cone` / `Torus` (`lxx:228-251`).
//!
//! `GeomAdaptor_Surface::load` replaces `mySurface` with the basis of a
//! `Geom_RectangularTrimmedSurface` (`GeomAdaptor_Surface.cxx:423-425`), so
//! every query that OCCT forwards to `mySurface` sees the peeled basis; the
//! helpers below peel the same way ([`SurfaceTool::basis`]).
//!
//! UNPORTED (no `BRepApprox` caller in the port yet, and the underlying
//! `SurfaceAdaptor` slice does not expose them):
//! - `NbUIntervals` / `NbVIntervals` / `UIntervals` / `VIntervals`
//!   (`lxx:46-87`, `BRepAdaptor_Surface::NbUIntervals` family);
//! - `UTrim` / `VTrim` (`lxx:89-113`, they return an `Adaptor3d_Surface` handle);
//! - `IsUClosed` / `IsVClosed` (`lxx:115-127`);
//! - `D2` / `D3` / `DN` (`lxx:170-204`);
//! - `Bezier` / `BSpline` / `AxeOfRevolution` / `Direction` / `BasisCurve`
//!   (`lxx:241-251`).
//!
//! The two `ThePSurfaceTool` calls of `ApproxInt_PrmPrmSvSurfaces.gxx`
//! (`:180`, `:207`) are `D1`, and those of `ApproxInt_ImpPrmSvSurfaces.gxx` are
//! `D1` / `GetType` / `Plane` / `Cylinder` / `Cone` / `Sphere` / `Torus` /
//! `Value` / `First*Parameter` / `Last*Parameter` / `Is*Periodic` / `*Period` /
//! `*Resolution`, which is exactly the slice below.

use std::sync::Arc;

use occt_core::gp::{GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom::surface::Surface;

use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `BRepApprox_SurfaceTool` (`hxx:38`).
pub struct SurfaceTool;

impl SurfaceTool {
    /// `FirstUParameter(S)` (`lxx:22-25`).
    pub fn first_u_parameter(s: &SurfaceAdaptor) -> f64 {
        s.first_u_parameter()
    }

    /// `FirstVParameter(S)` (`lxx:27-30`).
    pub fn first_v_parameter(s: &SurfaceAdaptor) -> f64 {
        s.first_v_parameter()
    }

    /// `LastUParameter(S)` (`lxx:32-35`).
    pub fn last_u_parameter(s: &SurfaceAdaptor) -> f64 {
        s.last_u_parameter()
    }

    /// `LastVParameter(S)` (`lxx:37-40`).
    pub fn last_v_parameter(s: &SurfaceAdaptor) -> f64 {
        s.last_v_parameter()
    }

    /// `Value(S, U, V)` (`lxx:143-149`).
    pub fn value(s: &SurfaceAdaptor, u: f64, v: f64) -> GpPnt {
        let (p, _, _) = s.d1(u, v).unwrap_or((GpPnt::zero(), GpVec::new(0.0, 0.0, 0.0), GpVec::new(0.0, 0.0, 0.0)));
        p
    }

    /// `D0(S, U, V, P)` (`lxx:151-157`).
    pub fn d0(s: &SurfaceAdaptor, u: f64, v: f64, p: &mut GpPnt) {
        *p = Self::value(s, u, v);
    }

    /// `D1(S, U, V, P, D1U, D1V)` (`lxx:159-168`).
    pub fn d1(s: &SurfaceAdaptor, u: f64, v: f64, p: &mut GpPnt, d1u: &mut GpVec, d1v: &mut GpVec) {
        if let Some((pt, du, dv)) = s.d1(u, v) {
            *p = pt;
            *d1u = du;
            *d1v = dv;
        }
    }

    /// `GetType(S)` (`lxx:222-226`).
    pub fn get_type(s: &SurfaceAdaptor) -> GeomAbsSurfaceType {
        s.get_type()
    }

    /// `NbSamplesU(S)` (`cxx:19-22`).
    pub fn nb_samples_u(_s: &SurfaceAdaptor) -> i32 {
        10
    }

    /// `NbSamplesV(S)` (`cxx:24-27`).
    pub fn nb_samples_v(_s: &SurfaceAdaptor) -> i32 {
        10
    }

    /// `NbSamplesU(S, U1, U2)` (`cxx:29-33`).
    pub fn nb_samples_u_range(_s: &SurfaceAdaptor, _u1: f64, _u2: f64) -> i32 {
        10
    }

    /// `NbSamplesV(S, V1, V2)` (`cxx:35-40`).
    pub fn nb_samples_v_range(_s: &SurfaceAdaptor, _v1: f64, _v2: f64) -> i32 {
        10
    }

    /// The surface `GeomAdaptor_Surface` keeps after peeling a
    /// `Geom_RectangularTrimmedSurface` (`GeomAdaptor_Surface.cxx:423-425`).
    /// `None` stands for a face without geometry.
    pub fn basis(s: &SurfaceAdaptor) -> Option<Arc<dyn Surface>> {
        let mut cur = s.surface()?.clone();
        while let Some(b) = cur.rectangular_trimmed_basis() {
            cur = b;
        }
        Some(cur)
    }

    /// `IsUPeriodic(S)` (`lxx:129-132`): `GeomAdaptor_Surface::IsUPeriodic`
    /// (`GeomAdaptor_Surface.cxx:946-949`).
    pub fn is_u_periodic(s: &SurfaceAdaptor) -> bool {
        Self::basis(s).map(|b| b.is_u_periodic()).unwrap_or(false)
    }

    /// `UPeriod(S)` (`lxx:134-137`): `GeomAdaptor_Surface::UPeriod`
    /// (`GeomAdaptor_Surface.cxx:953-957`) raises for a non-periodic surface;
    /// the port returns `0.0` because its only caller is guarded by
    /// [`SurfaceTool::is_u_periodic`].
    pub fn u_period(s: &SurfaceAdaptor) -> f64 {
        Self::basis(s).map(|b| b.u_period()).unwrap_or(0.0)
    }

    /// `IsVPeriodic(S)` (`lxx:139-141`).
    pub fn is_v_periodic(s: &SurfaceAdaptor) -> bool {
        Self::basis(s).map(|b| b.is_v_periodic()).unwrap_or(false)
    }

    /// `VPeriod(S)` (`lxx:139-141`).
    pub fn v_period(s: &SurfaceAdaptor) -> f64 {
        Self::basis(s).map(|b| b.v_period()).unwrap_or(0.0)
    }

    /// `UResolution(S, R3d)` (`lxx:206-209`):
    /// `GeomAdaptor_Surface::UResolution` (`GeomAdaptor_Surface.cxx:1818-1896`).
    pub fn u_resolution(s: &SurfaceAdaptor, r3d: f64) -> f64 {
        match Self::basis(s) {
            Some(b) => u_resolution(b.as_ref(), r3d),
            None => 0.0,
        }
    }

    /// `VResolution(S, R3d)` (`lxx:211-214`):
    /// `GeomAdaptor_Surface::VResolution` (`GeomAdaptor_Surface.cxx:1900-1959`).
    pub fn v_resolution(s: &SurfaceAdaptor, r3d: f64) -> f64 {
        match Self::basis(s) {
            Some(b) => v_resolution(b.as_ref(), r3d),
            None => 0.0,
        }
    }

    /// `Plane(S)` (`lxx:228-232`): `GeomAdaptor_Surface::Plane`
    /// (`GeomAdaptor_Surface.cxx:1963-1971`) raises `Standard_NoSuchObject`
    /// unless `GetType() == GeomAbs_Plane`; the port panics on the same
    /// condition.
    pub fn plane(s: &SurfaceAdaptor) -> GpPln {
        Self::basis(s)
            .and_then(|b| b.gp_pln())
            .expect("BRepApprox_SurfaceTool::Plane: GeomAdaptor_Surface.cxx:1965")
    }

    /// `Cylinder(S)` (`lxx:234-238`).
    pub fn cylinder(s: &SurfaceAdaptor) -> GpCylinder {
        Self::basis(s)
            .and_then(|b| b.gp_cylinder())
            .expect("BRepApprox_SurfaceTool::Cylinder: GeomAdaptor_Surface.cxx:1980")
    }

    /// `Cone(S)` (`lxx:240-244`).
    pub fn cone(s: &SurfaceAdaptor) -> GpCone {
        Self::basis(s)
            .and_then(|b| b.gp_cone())
            .expect("BRepApprox_SurfaceTool::Cone: GeomAdaptor_Surface.cxx:1994")
    }

    /// `Sphere(S)` (`lxx:246-250`).
    pub fn sphere(s: &SurfaceAdaptor) -> GpSphere {
        Self::basis(s)
            .and_then(|b| b.gp_sphere())
            .expect("BRepApprox_SurfaceTool::Sphere: GeomAdaptor_Surface.cxx:2008")
    }

    /// `Torus(S)` (`lxx:246-250`).
    pub fn torus(s: &SurfaceAdaptor) -> GpTorus {
        Self::basis(s)
            .and_then(|b| b.gp_torus())
            .expect("BRepApprox_SurfaceTool::Torus: GeomAdaptor_Surface.cxx:2022")
    }
}
