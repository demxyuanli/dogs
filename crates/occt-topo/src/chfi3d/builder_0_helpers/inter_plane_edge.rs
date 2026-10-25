//! `ChFi3d_InterPlaneEdge`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:718-767`
//! (OCCT 8.0.0).
//!
//! The function intersects a curve with a plane and answers the parameter of the
//! first intersection inside the curve's own range (plus `tolc`):
//!
//! ```text
//! IntCurveSurface_HInter Intersection;
//! Intersection.Perform(C, Plan);
//! for (iip = 1; iip <= nbp; iip++) {
//!   CW = Point(iip).W();
//!   if (C->IsPeriodic()) CW = ElCLib::InPeriod(CW, uf - tolc, uf - tolc + C->Period());
//!   if (uf - tolc <= CW && ul + tolc >= CW) {            // keep it
//!     if (isol == 0) { isol = iip; W = CW; }
//!     else if (Sens && CW < W)  { W = CW; isol = iip; }   // first one
//!     else if (!Sens && CW > W) { W = CW; isol = iip; }   // last one
//!   }
//! }
//! ```
//!
//! `false` is OCCT's `isol == 0`; `w` is left untouched there.
//!
//! Divergences that are not heuristics:
//! - `IntCurveSurface_HInter::Perform(C, Plan)` becomes
//!   [`intcurvesurface::perform_curve_surface`] (the port of the same class).
//!   That port needs the surface's `(u0, v0, u1, v1)` window while OCCT's
//!   `GeomAdaptor_Surface(Geom_Plane)` has infinite bounds; the caller passes
//!   the window, and [`UNBOUNDED_PLANE_WINDOW`] is the representable stand-in
//!   for "no bound" (it only feeds the `ComputeAppendPoint` window test, which
//!   an infinite plane never fails),
//! - the curve range is the adaptor's own `FirstParameter` / `LastParameter`,
//!   exactly what `Intersection.Perform(C, Plan)` reads,
//! - `Intersection.Point(iip).W()` is the curve parameter, i.e. the port's
//!   `IntersectionPoint::param`.

use occt_core::elib::clib::in_period;
use occt_geom::{Curve, Surface};

use crate::intcurvesurface::perform_curve_surface;

/// Stand-in for the infinite `U`/`V` bounds of `GeomAdaptor_Surface` on a
/// `Geom_Plane` (`GeomAdaptor_Surface.cxx`, unbounded plane arm). The port's
/// `IntCurveSurface_HInter` requires a finite window; `1e100` keeps every
/// comparison in `ComputeAppendPoint` inside the window while leaving the
/// Euclidean arithmetic untouched.
pub const UNBOUNDED_PLANE_WINDOW: f64 = 1.0e100;

/// The plane window used for an unbounded plane surface.
pub fn unbounded_plane_window() -> (f64, f64, f64, f64) {
    let r = UNBOUNDED_PLANE_WINDOW;
    (-r, r, -r, r)
}

/// `ChFi3d_InterPlaneEdge(Plan, C, W, Sens, tolc)` (`cxx:718-767`).
///
/// `sens` selects the smallest (`true`) or the largest (`false`) intersection
/// parameter inside `[uf - tolc, ul + tolc]`. `false` is OCCT's `isol == 0`.
pub fn inter_plane_edge(
    plan: &dyn Surface,
    c: &dyn Curve,
    w: &mut f64,
    sens: bool,
    tolc: f64,
) -> bool {
    let uf = c.first_parameter();
    let ul = c.last_parameter();
    let dir = (uf, ul);
    let Some(intersection) = perform_curve_surface_plane(plan, c, dir) else {
        return false;
    };

    let mut isol = 0usize;
    for (iip, &param) in intersection.iter().enumerate() {
        let mut cw = param;
        if c.is_periodic() {
            cw = in_period(cw, uf - tolc, uf - tolc + c.period());
        }
        if uf - tolc <= cw && ul + tolc >= cw {
            if isol == 0 {
                isol = iip + 1;
                *w = cw;
            } else if sens && cw < *w {
                *w = cw;
                isol = iip + 1;
            } else if !sens && cw > *w {
                *w = cw;
                isol = iip + 1;
            }
        }
    }
    isol != 0
}

/// `Intersection.Perform(C, Plan)` + the `NbPoints` / `Point(iip).W()` reads,
/// answering the point parameters in the port's sorted order.
fn perform_curve_surface_plane(
    plan: &dyn Surface,
    c: &dyn Curve,
    dir: (f64, f64),
) -> Option<Vec<f64>> {
    let res = perform_curve_surface(c, plan, dir, unbounded_plane_window()).ok()?;
    Some(res.points().iter().map(|p| p.param).collect())
}
