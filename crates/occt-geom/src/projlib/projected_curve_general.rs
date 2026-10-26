//! `ProjLib_ProjectedCurve::Perform`, `default:` arm — the approximation of a
//! curve on a surface outside the five elementary ones and the two patch
//! types.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ProjectedCurve.cxx:539-715`
//! (OCCT 8.0.0).
//!
//! The arm has three early returns of its own: an unavailable
//! `ComputeTolU`/`ComputeTolV`, an empty `ProjLib_HCompProjectedCurve` result
//! list, and a failed `Approx_CurveOnSurface::Perform`. [`general_approx`]
//! reports them as `None`, which the caller turns into the same bare return
//! (the shared `cxx:717-886` tail is skipped, as in OCCT).
//!
//! UNPORTED, with the OCCT lines they occupy:
//! - `cxx:553-610`: the singularity handling of a surface of revolution —
//!   `AxeOfRevolution` / `BasisCurve` on the adaptor, `Extrema_ExtPC` over the
//!   basis curve for the singular parameter, and `Adaptor3d_Curve::Trim`.
//!   `IsTrimmed` therefore stays `false`.
//! - `cxx:658-712`: the `ExtendC2d` / `GeomLib::SameRange` / `RemoveKnot` tail,
//!   which only runs when `IsTrimmed` is set.
//!
//! A curve whose end lies on the revolution axis is thus approximated without
//! the end extension.

use std::sync::Arc;

use occt_core::precision;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;

use crate::approx::curve_on_surface::ApproxCurveOnSurface;
use crate::curve::Curve;
use crate::projlib::comp_projected_curve::{CompProjectedCurve, ProjectedCurve2d};
use crate::projlib::projected_curve::{compute_tol_u, compute_tol_v, Constraint};
use crate::surface::Surface;

/// `ProjLib_ProjectedCurve::Perform`'s `default:` arm (`cxx:539-715`): the
/// 2D B-spline of the normal projection plus the tolerance it implies.
pub(crate) fn general_approx(
    surface: &Arc<dyn Surface>,
    curve: &Arc<dyn Curve>,
    tolerance: f64,
    max_dist: f64,
    deg_max: i32,
    max_segments: i32,
    bnd_pnt: Constraint,
) -> Option<(Arc<Geom2dBSplineCurve>, f64)> {
    // `cxx:616-624`: the 2D tolerances and their diagonal.
    let tol_u = compute_tol_u(surface.as_ref(), tolerance)?;
    let tol_v = compute_tol_v(surface.as_ref(), tolerance)?;
    let a_tol_u = tol_u.max(precision::CONFUSION);
    let a_tol_v = tol_v.max(precision::CONFUSION);
    let a_tol2d = (a_tol_u * a_tol_u + a_tol_v * a_tol_v).sqrt();

    // `cxx:617-627`: `ProjLib_HCompProjectedCurve(mySurface, myCurve, aTolU,
    // aTolV, aMaxDist)`.
    let mut a_max_dist = 100.0 * tolerance;
    if max_dist > 0.0 {
        a_max_dist = max_dist;
    }
    let hprojector = CompProjectedCurve::with_tolerances_and_max_dist(
        surface.clone(),
        curve.clone(),
        a_tol_u,
        a_tol_v,
        a_max_dist,
    );

    // `cxx:630-640`: within `ProjLib` the projection must be one single
    // curve, so part 1 is what gets approximated.
    if hprojector.nb_curves() <= 0 {
        return None;
    }
    let (udeb, ufin) = hprojector.bounds(1)?;

    // `cxx:643-664`: `Only3d = false`, `Only2d = true`, `C1` (`C0` for
    // `AppParCurves_PassPoint`), degree 14 and 16 segments unless overridden.
    let continuity = if bnd_pnt == Constraint::PassPoint { 0 } else { 2 };
    let max_degree = if deg_max > 0 { deg_max } else { 14 };
    let max_seg = if max_segments > 0 { max_segments } else { 16 };

    // `cxx:666-667`: `Approx_CurveOnSurface(HProjector, mySurface, Udeb, Ufin,
    // myTolerance)` — note the *untrimmed* surface, the widened window lives
    // inside the projection adaptor.
    let hpc = Arc::new(ProjectedCurve2d::new(hprojector));
    let mut appr = ApproxCurveOnSurface::new(hpc, surface.clone(), udeb, ufin, tolerance);
    if appr
        .perform(max_seg, max_degree, continuity, false, true)
        .is_err()
    {
        // `Approx_CurveOnSurface::Perform` raises
        // `Standard_ConstructionError` for `Only3d && Only2d` (`cxx:367-370`)
        // and lets the `AdvApprox` exceptions through; OCCT then leaves
        // `myResult` undone.
        return None;
    }

    // `cxx:668-676`: `myTolerance *= sqrt(errU^2 + errV^2) / aTol2d`.
    let a_res = appr.curve2d()?;
    let err_u = appr.max_error2du();
    let err_v = appr.max_error2dv();
    let a_new_tol2d = (err_u * err_u + err_v * err_v).sqrt();
    Some((Arc::new(a_res.clone()), tolerance * a_new_tol2d / a_tol2d))
}
