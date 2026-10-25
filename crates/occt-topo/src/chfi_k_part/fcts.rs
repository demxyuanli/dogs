//! `ChFiKPart_ComputeData_Fcts` — the file statics every `ChFiKPart` surface
//! builder shares.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_Fcts.cxx:34-155`
//! (OCCT 8.0.0).
//!
//! The ported functions are `ChFiKPart_InPeriod`, `ChFiKPart_PCurve`,
//! `ChFiKPart_IndexCurveInDS` and `ChFiKPart_IndexSurfaceInDS`.
//! `ChFiKPart_ProjPC` is UNPORTED: it reads the projection result back through
//! `ProjLib_ProjectedCurve::Line()` / `::Bezier()` / `::BSpline()`
//! (`cxx:83-135`), and the port's `occt_geom::projlib::ProjectedCurve` keeps
//! that result in a `ProjectorState` with no such accessors yet. Its two call
//! sites are the spherical corners (`ChFiKPart_ComputeData_Sphere.cxx:139-147`,
//! `ChFiKPart_ComputeData_Rotule.cxx:143-151`), so the gap is closed with that
//! batch.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;

use crate::topopebrep::ds::{Curve, DataStructure, Surface};

/// `ChFiKPart_InPeriod(U, UFirst, ULast, Eps)` (`cxx:36-52`): shift `U` by
/// whole periods until it lies in `[UFirst, ULast]`, then clamp it up to
/// `UFirst` when it fell below.
pub fn ch_fi_kpart_in_period(u: f64, u_first: f64, u_last: f64, eps: f64) -> f64 {
    let period = u_last - u_first;
    let mut u = u;
    while eps < u_first - u {
        u += period;
    }
    while eps > u_last - u {
        u -= period;
    }
    if u < u_first {
        u = u_first;
    }
    u
}

/// `ChFiKPart_PCurve(UV1, UV2, Pardeb, Parfin)` (`cxx:60-75`): the straight
/// line between two `(u, v)` points as a degree-1 B-spline with knots
/// `{Pardeb, Pardeb, Parfin, Parfin}` — OCCT builds it as
/// `Geom2d_BSplineCurve(p, k, m, 1)` with `k = {Pardeb, Parfin}` and
/// `m = {2, 2}` (`cxx:66-74`), which `updateKnots` expands to exactly that
/// flat vector.
pub fn ch_fi_kpart_p_curve(
    uv1: GpPnt2d,
    uv2: GpPnt2d,
    pardeb: f64,
    parfin: f64,
) -> Result<Geom2dBSplineCurve, String> {
    Geom2dBSplineCurve::new(
        vec![uv1.x(), uv2.x()],
        vec![uv1.y(), uv2.y()],
        vec![pardeb, pardeb, parfin, parfin],
        1,
    )
    .map_err(|e| e.to_string())
}

/// `ChFiKPart_IndexCurveInDS(C, DStr)` (`cxx:142-145`):
/// `DStr.AddCurve(TopOpeBRepDS_Curve(C, 0.))`. `C` is a handle, so OCCT also
/// accepts a null one — `ChFiKPart_ComputeData_FilPlnCyl.cxx:492` /
/// `_FilPlnCon.cxx:240` pass the still-null `GCircPln` of their sphere arm.
pub fn ch_fi_kpart_index_curve_in_ds(
    c: Option<Arc<dyn GeomCurve>>,
    ds: &mut DataStructure,
) -> i32 {
    let mut curve = Curve::with_curve(c, 0.0, false);
    ds.add_curve(&mut curve)
}

/// `ChFiKPart_IndexSurfaceInDS(S, DStr)` (`cxx:152-155`):
/// `DStr.AddSurface(TopOpeBRepDS_Surface(S, 0.))`.
pub fn ch_fi_kpart_index_surface_in_ds(s: Arc<dyn GeomSurface>, ds: &mut DataStructure) -> i32 {
    ds.add_surface(&Surface::with_tolerance(s, 0.0))
}
