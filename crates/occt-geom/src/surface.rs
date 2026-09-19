//! Abstract 3D parametric surface. Source: `Geom_Surface.hxx`
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpCone, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus, GpVec, GpTrsf};

use crate::curve::Curve;

pub trait Surface: Send + Sync {
    fn d0(&self, u: f64, v: f64) -> GpPnt;
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec);
    /// Default D2 by central difference of D1. Analytic surfaces override.
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        let h = 1e-6;
        let (p, du, dv) = self.d1(u, v);
        let (_, du_u, dv_u) = self.d1(u + h, v);
        let (_, du_v, dv_v) = self.d1(u, v + h);
        let d2u = du_u.subtracted(&du).divided(h);
        let d2v = dv_v.subtracted(&dv).divided(h);
        let d2uv = dv_u.subtracted(&dv).divided(h);
        let _ = du_v;
        (p, du, dv, d2u, d2v, d2uv)
    }
    fn value(&self, u: f64, v: f64) -> GpPnt { self.d0(u, v) }
    fn u_range(&self) -> (f64, f64);
    fn v_range(&self) -> (f64, f64);
    fn is_u_periodic(&self) -> bool { false }
    fn is_v_periodic(&self) -> bool { false }
    fn continuity(&self) -> u8;
    fn transform(&mut self, t: &GpTrsf);
    fn clone_dyn(&self) -> Box<dyn Surface>;
    fn transformed(&self, t: &GpTrsf) -> Box<dyn Surface> { let mut s = self.clone_dyn(); s.transform(t); s }
    /// `STANDARD_TYPE(Geom_BSplineSurface)` — extra `BRepTools::AddUVBounds` probe.
    fn is_bspline_surface(&self) -> bool { false }
    /// Clone as `Geom_BSplineSurface` when this is one (or a Bezier converted
    /// the way `Geom_OsculatingSurface::Init` does).
    fn osculating_bspline(&self) -> Option<crate::bspline_surface::GeomBSplineSurface> {
        None
    }
    /// `Geom_Surface::IsUClosed`.
    fn is_u_closed(&self) -> bool { self.is_u_periodic() }
    /// `Geom_Surface::IsVClosed`.
    fn is_v_closed(&self) -> bool { self.is_v_periodic() }
    /// `(RefRadius, SemiAngle)` when this is a `Geom_ConicalSurface`.
    /// Source: `Adaptor3d_Surface::GetType() == GeomAbs_Cone` then `Cone()`.
    fn cone_ref(&self) -> Option<(f64, f64)> { None }
    /// `Geom_SphericalSurface::Sphere`.
    fn gp_sphere(&self) -> Option<GpSphere> { None }
    /// `STANDARD_TYPE(Geom_ToroidalSurface)` then `Torus()`.
    fn gp_torus(&self) -> Option<GpTorus> { None }
    /// `Adaptor3d_Surface::GetType() == GeomAbs_Cylinder` then `Cylinder()`.
    fn gp_cylinder(&self) -> Option<GpCylinder> { None }
    /// `Adaptor3d_Surface::GetType() == GeomAbs_Plane` then `Plane()`.
    fn gp_pln(&self) -> Option<GpPln> { None }
    /// `Adaptor3d_Surface::GetType() == GeomAbs_Cone` then `Cone()`.
    fn gp_cone(&self) -> Option<GpCone> { None }
    /// `STANDARD_TYPE(Geom_SurfaceOfRevolution)`.
    fn is_surface_of_revolution(&self) -> bool { false }
    /// `STANDARD_TYPE(Geom_OffsetSurface)`.
    fn is_offset_surface(&self) -> bool { false }
    /// `Geom_OffsetSurface::BasisSurface`.
    fn offset_basis_surface(&self) -> Option<Arc<dyn Surface>> { None }
    /// `Geom_OffsetSurface::Offset`.
    fn offset_distance(&self) -> Option<f64> { None }
    /// `STANDARD_TYPE(Geom_RectangularTrimmedSurface)` then `BasisSurface()`.
    fn rectangular_trimmed_basis(&self) -> Option<Arc<dyn Surface>> { None }
    /// `STANDARD_TYPE(Geom_SurfaceOfLinearExtrusion)`.
    fn is_surface_of_linear_extrusion(&self) -> bool { false }
    /// `Geom_Surface::UIso`: for an analytic surface the partner is
    /// `GeomAdaptor_Surface::UIso` (`GeomAdaptor_Surface.cxx`) over
    /// `Geom_Plane::UIso` / `Geom_CylindricalSurface::UIso` /
    /// `Geom_ConicalSurface::UIso` / `Geom_SphericalSurface::UIso` /
    /// `Geom_ToroidalSurface::UIso`; for a parametric patch it is
    /// `Geom_BSplineSurface::UIso` / `Geom_BezierSurface::UIso`.
    fn u_iso_curve(&self, _u: f64) -> Option<Arc<dyn Curve>> { None }
    /// `Geom_Surface::VIso`; the mirror of [`Surface::u_iso_curve`].
    fn v_iso_curve(&self, _v: f64) -> Option<Arc<dyn Curve>> { None }
    /// `Geom_Surface::UPeriod`.
    fn u_period(&self) -> f64 {
        let (a, b) = self.u_range();
        if a.is_finite() && b.is_finite() { (b - a).abs() } else { 0.0 }
    }
    /// `Geom_Surface::VPeriod`.
    fn v_period(&self) -> f64 {
        let (a, b) = self.v_range();
        if a.is_finite() && b.is_finite() { (b - a).abs() } else { 0.0 }
    }
    /// `Adaptor3d_Surface::NbUPoles` (Bezier / BSpline). Default 2 so
    /// `NbPoles - 1` is a single interval.
    fn nb_u_poles(&self) -> i32 { 2 }
    /// `Adaptor3d_Surface::NbVPoles`.
    fn nb_v_poles(&self) -> i32 { 2 }
    /// `Adaptor3d_Surface::NbUIntervals`. Default one span.
    fn nb_u_intervals(&self, _continuity: u8) -> i32 { 1 }
    /// `Adaptor3d_Surface::NbVIntervals`.
    fn nb_v_intervals(&self, _continuity: u8) -> i32 { 1 }
    /// `Adaptor3d_Surface::UIntervals` including the range ends.
    fn u_intervals(&self, _continuity: u8) -> Vec<f64> {
        let (a, b) = self.u_range();
        vec![a, b]
    }
    /// `Adaptor3d_Surface::VIntervals` including the range ends.
    fn v_intervals(&self, _continuity: u8) -> Vec<f64> {
        let (a, b) = self.v_range();
        vec![a, b]
    }
    /// `Adaptor3d_Surface::UDegree` (Bezier / BSpline). Default 1.
    fn u_degree(&self) -> i32 { 1 }
    /// `Adaptor3d_Surface::VDegree` (Bezier / BSpline). Default 1.
    fn v_degree(&self) -> i32 { 1 }
    /// `Adaptor3d_Surface::BasisCurve` on `Geom_SurfaceOfRevolution`.
    fn revolution_basis_curve(&self) -> Option<Arc<dyn Curve>> { None }
    /// `Adaptor3d_Surface::AxeOfRevolution` on `Geom_SurfaceOfRevolution`.
    fn revolution_axis(&self) -> Option<GpAx1> { None }
    /// `Adaptor3d_Surface::BasisCurve` on `Geom_SurfaceOfLinearExtrusion`.
    fn extrusion_basis_curve(&self) -> Option<Arc<dyn Curve>> { None }
    /// `Geom_SweptSurface::Direction` of a `Geom_SurfaceOfLinearExtrusion`.
    fn extrusion_direction(&self) -> Option<GpDir> { None }
    /// `Geom_BSplineSurface::Resolution` / `GeomAdaptor_Surface::UResolution`
    /// + `VResolution`. `None` falls back to a finite-difference estimate.
    fn uv_resolution(&self, _r3d: f64) -> Option<(f64, f64)> { None }
    /// `Geom_Surface::EvalDN` (`Geom_Surface.hxx:260`, pure virtual).
    /// Default is the unported arm (zero). `GeomBSplineSurface` overrides.
    fn eval_dn(&self, _u: f64, _v: f64, _nu: i32, _nv: i32) -> GpVec {
        GpVec::new(0.0, 0.0, 0.0)
    }
}
