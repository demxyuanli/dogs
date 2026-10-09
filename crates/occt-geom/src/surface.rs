//! Abstract 3D parametric surface. Source: `Geom_Surface.hxx`
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpCone, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus, GpVec, GpTrsf};

use crate::curve::Curve;

pub trait Surface: Send + Sync {
    fn d0(&self, u: f64, v: f64) -> GpPnt;
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec);
    /// Default D2 by central difference of D1.
    ///
    /// **UNPORTED (audit A15)**: OCCT's `Geom_Surface::D2` is pure virtual — every
    /// concrete class computes it analytically (`ElSLib::*D2` for the elementary
    /// surfaces, `BSplSLib::D2` for B-spline/Bezier, `Geom_OffsetSurfaceUtils`
    /// for offsets) and no finite-difference fallback exists anywhere in OCCT.
    /// The five elementary surfaces (plane / cylinder / cone / sphere / torus),
    /// `GeomBSplineSurface`, `GeomRectangularTrimmedSurface`,
    /// `GeomSurfaceOfRevolution` and `GeomSurfaceOfLinearExtrusion` override this
    /// (`Geom_RevolutionUtils::CalculateD2` / `Geom_ExtrusionUtils::CalculateD2`);
    /// the remaining implementors are `GeomBezierSurface`, `GeomOffsetSurface` and
    /// the port-only `surface_fit::GridSurface` (audit A15; the abandoned
    /// `surface_to_grid` module was deleted in batch 65 as a dead non-OCCT type).
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
    /// `Adaptor3d_Surface::GetType() == GeomAbs_BezierSurface`
    /// (`Geom_BezierSurface`; `GeomAdaptor_Surface.cxx:480-487`).
    fn is_bezier_surface(&self) -> bool { false }
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
    /// `Geom_Surface::UReversed` — a copy of the same dynamic type whose U
    /// parametrisation is reversed.
    ///
    /// **UNPORTED** for the surface classes without a branch here: OCCT defines
    /// it on every `Geom_Surface` (the five elementary surfaces return a copy
    /// with `UReverse()` applied — `gp_Pln.hxx:101` uses `XReverse`, every other
    /// elementary surface `YReverse`, `gp_Cylinder.hxx:86`;
    /// `Geom_BSplineSurface::UReverse` reverses the poles and knots in U
    /// (`Geom_BSplineSurface_1.cxx:1572-1592`), ported for a non-periodic U only;
    /// `Geom_SurfaceOfRevolution::UReverse` reverses the axis direction;
    /// `Geom_SurfaceOfLinearExtrusion::UReverse` reverses the basis curve;
    /// `Geom_OffsetSurface::UReverse` U-reverses the basis and negates the
    /// offset. `None` means the
    /// caller cannot take `BRepToIGES_BRShell::TransferFace`'s REVERSED branch
    /// (`BRepToIGES_BRShell.cxx:135-141`) for this surface.
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> { None }
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

    /// `Geom_BSplineSurface::Poles()` — the pole grid, U-major (`[u][v]`).
    /// Needed by the IGES/STEP writers, which emit the poles as data
    /// (`GeomToIGES_GeomSurface::TransferBSplineSurface`).
    fn bspline_surface_poles(&self) -> Option<&[Vec<GpPnt>]> { None }

    /// `Geom_BSplineSurface::UKnots()` (the flattened U knot vector).
    fn bspline_surface_uknots(&self) -> Option<&[f64]> { None }

    /// `Geom_BSplineSurface::VKnots()` (the flattened V knot vector).
    fn bspline_surface_vknots(&self) -> Option<&[f64]> { None }

    /// `Geom_BSplineSurface::Weights()`; `None` for a non-rational surface
    /// (`Geom_BSplineSurface::IsRational() == false`), whose weights are all 1.
    fn bspline_surface_weights(&self) -> Option<&[Vec<f64>]> { None }
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
