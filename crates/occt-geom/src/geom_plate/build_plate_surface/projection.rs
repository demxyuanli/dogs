//! `ProjectCurve` / `ProjectedCurve` / `ProjectPoint`
//! (`GeomPlate_BuildPlateSurface.cxx:253-376`).
//!
//! `ProjLib_HCompProjectedCurve` is `projlib::CompProjectedCurve`; the
//! `Adaptor2d_Curve2d` the `Value(U)` call sites see is
//! `projlib::ProjectedCurve2d`, the same wrapper `ProjLib_CompProjectedCurve`
//! supplies in OCCT because it *is* an `Adaptor2d_Curve2d` (`hxx:36`).
//!
//! `myProj` (`Extrema_ExtPS`) is re-initialised by `Perform` on the current
//! `mySurfInit` with `(Uinf, Usup, Vinf, Vsup, TolU, TolV)` (`cxx:516-517`,
//! `:542-543`, `:565-566`, `:1799-1801`, `:1876-1878`); the port builds the
//! engine per call from the same surface and the stored `tol_u` / `tol_v`,
//! because [`crate::extrema_surf::ExtPs`] borrows the surface.
//!
//! **UNPORTED / deviation**: `ProjectedCurve`'s
//! `HProjector->Trim(First2, Last2, Precision::PConfusion())` (`cxx:334-336`)
//! only trims the adaptor to the range `Bounds(1, ...)` already returned, so the
//! port returns the untrimmed wrapper. The `NbCurves() != 1` and the
//! `max(myTolU, myTolV)` tests are unchanged.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::PCONFUSION;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::curve::Curve2d;

use crate::approx::curve_on_surface::ApproxCurveOnSurface;
use crate::curve::Curve;
use crate::extrema_surf::ExtPs;
use crate::projlib::{CompProjectedCurve, ProjectedCurve2d};

use super::BuildPlateSurface;

impl BuildPlateSurface {
    /// `ProjectCurve(Curv)` (`cxx:253-301`): projects a curve on the initial
    /// (planar) surface and returns the 2D approximation, or a degenerate 2D
    /// Bezier when the projection collapses to a point.
    pub(super) fn project_curve(
        &self,
        curv: &Arc<dyn Curve>,
    ) -> Result<Option<Arc<dyn Curve2d>>, String> {
        let surf = self
            .surf_init
            .clone()
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial surface".to_string())?;

        // `cxx:258-262`.
        let mut projector = CompProjectedCurve::with_tolerances(
            surf.clone(),
            curv.clone(),
            self.tol3d / 10.0,
            self.tol3d / 10.0,
        );

        // `cxx:265-268`.
        let udeb_check = curv.first_parameter();
        let ufin_check = curv.last_parameter();
        let (proj_udeb, proj_ufin) = projector
            .bounds(1)
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no projected part".to_string())?;

        // `cxx:269-287`.
        if projector.nb_curves() != 1
            || (udeb_check - proj_udeb).abs() > PCONFUSION
            || (ufin_check - proj_ufin).abs() > PCONFUSION
        {
            let (p2d, is_single) = projector
                .is_single_pnt(1)
                .ok_or_else(|| "GeomPlate_BuildPlateSurface: no projected part".to_string())?;
            if is_single {
                // `cxx:276-279`: `NCollection_Array1<gp_Pnt2d> poles(1, 2)`
                // initialised with the same point.
                let poles = vec![p2d, p2d];
                let curve = Geom2dBezierCurve::new(poles)
                    .map_err(|e| format!("Geom2d_BezierCurve: {e}"))?;
                Ok(Some(Arc::new(curve)))
            } else {
                // `cxx:282`: no continuous solution.
                Ok(None)
            }
        } else {
            // `cxx:290-299`.
            let max_seg = 20 + projector.nb_intervals(5).unwrap_or(0);
            let u2d: Arc<dyn Curve2d> = Arc::new(ProjectedCurve2d::new(projector));
            let mut appr =
                ApproxCurveOnSurface::new(u2d, surf, proj_udeb, proj_ufin, self.tol3d);
            // `Perform(MaxSeg, 10, GeomAbs_C1, false, true)`: `Shape::C1` is
            // the `continuity` encoding of `Approx_CurveOnSurface::Perform`
            // (`approx/curve_on_surface.rs:355-361`).
            appr.perform(max_seg, 10, 2, false, true)
                .map_err(|e| e.to_string())?;
            // `cxx:299`: `Curve2d = appr.Curve2d()`; a failed `Perform` leaves
            // the handle null.
            Ok(appr
                .curve2d()
                .cloned()
                .map(|c| Arc::new(c) as Arc<dyn Curve2d>))
        }
    }

    /// `ProjectedCurve(Curv)` (`cxx:306-346`): the projector itself when it is
    /// one continuous part spanning the whole curve.
    pub(super) fn projected_curve(
        &self,
        curv: &Arc<dyn Curve>,
    ) -> Result<Option<Arc<dyn Curve2d>>, String> {
        let surf = self
            .surf_init
            .clone()
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial surface".to_string())?;

        // `cxx:313-315`.
        let projector = CompProjectedCurve::with_tolerances(
            surf,
            curv.clone(),
            self.tol_u / 10.0,
            self.tol_v / 10.0,
        );

        // `cxx:316-319`.
        if projector.nb_curves() != 1 {
            return Ok(None);
        }

        // `cxx:325-341`.
        let first1 = curv.first_parameter();
        let last1 = curv.last_parameter();
        let Some((first2, last2)) = projector.bounds(1) else {
            return Ok(None);
        };
        let tol = self.tol_u.max(self.tol_v);
        if (first1 - first2).abs() <= tol && (last1 - last2).abs() <= tol {
            // `cxx:334-336` trims to `[first2, last2]`, which is the range
            // `bounds(1)` just reported (see the module note).
            Ok(Some(Arc::new(ProjectedCurve2d::new(projector))))
        } else {
            Ok(None)
        }
    }

    /// `ProjectPoint(p3d)` (`cxx:353-376`): the nearest parameter of the point
    /// on the initial surface, in the `(TolU, TolV)` window `Perform` set up.
    pub(super) fn project_point(&self, p3d: &GpPnt) -> Result<GpPnt2d, String> {
        let surf = self
            .surf_init
            .clone()
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial surface".to_string())?;

        // `myProj.Perform(p3d)` (`cxx:358`).
        let mut ext = ExtPs::with_surface(p3d, surf.as_ref(), self.tol_u, self.tol_v);
        if !ext.is_done() {
            return Err("GeomPlate_BuildPlateSurface: Extrema_ExtPS not done".to_string());
        }
        // `cxx:359-368`: the smallest square distance of the solution set.
        let mut nearest = 1usize;
        let nb = ext.nb_ext();
        if nb > 1 {
            let mut dist2mini = ext.square_distance(1);
            for i in 2..=nb {
                let d = ext.square_distance(i);
                if d < dist2mini {
                    dist2mini = d;
                    nearest = i;
                }
            }
        }
        // `cxx:369-372`: `Extrema_POnSurf::Parameter(u, v)`.
        let (u, v, _) = ext.point(nearest);
        Ok(GpPnt2d::new(u, v))
    }
}
