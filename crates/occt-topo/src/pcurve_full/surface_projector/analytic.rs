//! Analytic projected-curve helpers: the `generateCurvePoints` sample count (`ShapeConstruct_ProjectCurveOnSurface.cxx:531-556`) and the basis-surface unwrap behind the `ProjLib` projection route.

use super::*;

/// Continuous projection of a 3D curve onto an analytic surface — the OCCT
/// `ProjLib` `projected curve` route. Each sample is projected and then
/// unwrapped across the `u`/`v` seam so the polyline stays on the *same sheet*
/// of a self-intersecting surface (a spindle torus), which a point-wise
/// projection cannot guarantee. Restricted to the torus: for the other analytic
/// surfaces the point-wise grid search is already exact (their `d0` inverse has
/// no sheet ambiguity), so the caller keeps that path. Returns `None` otherwise.
/// `ShapeConstruct_ProjectCurveOnSurface.cxx:71`.
pub(crate) const THE_NCONTROL: usize = 23;

/// `generateCurvePoints` count (`cxx:531-556`).
pub(crate) fn generate_curve_point_count(curve: &dyn Curve, first: f64, last: f64) -> usize {
    let mut n = THE_NCONTROL;
    if curve.bspline_poles().is_some() {
        let knots = curve.parameter_intervals(0);
        let mut used = 0usize;
        for w in knots.windows(2) {
            if w[1] > first && w[0] < last {
                used += 1;
            }
        }
        let deg = curve.nurbs_degree().unwrap_or(1);
        let min_pnt = used * (deg + 1);
        while n < min_pnt {
            n += THE_NCONTROL - 1;
        }
    }
    n.max(2)
}

/// `ShapeConstruct_ProjectCurveOnSurface::projectAnalytic`
/// (`cxx:848-897`): the plane behind the surface itself, or behind one
/// `Geom_RectangularTrimmedSurface` / `Geom_OffsetSurface` wrapper. The wrapper
/// is resolved with `else if`, exactly as the C++ does, so only one level is
/// peeled.
///
/// The projected curve returned by `ProjLib_ProjectedCurve` on the wrapper is
/// the same `(u, v)` as the plane's because both wrappers keep the basis
/// parameterization: `Geom_RectangularTrimmedSurface::Value` clamps the
/// arguments but returns the basis point, and
/// `Geom_OffsetSurface::Value = basis + Offset*Normal`, so only the origin of
/// the frame moves (along the normal, which the `(u, v)` projection ignores).
pub(crate) fn plane_behind_wrapper(surf: &dyn Surface) -> Option<GpPln> {
    if let Some(pln) = surf.gp_pln() {
        return Some(pln);
    }
    if let Some(basis) = surf.rectangular_trimmed_basis() {
        if let Some(pln) = basis.gp_pln() {
            return Some(pln);
        }
    } else if let Some(basis) = surf.offset_basis_surface() {
        if let Some(pln) = basis.gp_pln() {
            return Some(pln);
        }
    }
    None
}

/// `ShapeConstruct_ProjectCurveOnSurface::projectAnalytic(aCurve3DTrim)`
/// (`cxx:848-897`): a `None` result means "no analytic plane projection", and
/// the caller falls back to the sampled approximation of `Perform`.
pub(crate) fn project_analytic(curve: &dyn Curve, surf: &dyn Surface) -> Option<Arc<dyn Curve2d>> {
    let pln = plane_behind_wrapper(surf)?;
    crate::pcurve::project_curve_on_plane(curve, &pln)
}

/// `GeomAdaptor_Surface::load` (`GeomAdaptor_Surface.cxx:423-430`) recurses
/// through every `Geom_RectangularTrimmedSurface` and stores the *basis*
/// geometry, so `SurfAdapt.Cylinder()` / `Cone()` / `Sphere()` / `Torus()`
/// answer for a trimmed face. The Rust `Surface` trait accessors do not peel
/// that wrapper, so the analytic arms of [`value_of_uv`] must look through it.
pub(crate) fn analytic_basis(s: &dyn Surface) -> Option<Arc<dyn Surface>> {
    let mut deepest: Option<Arc<dyn Surface>> = None;
    let mut cur = s.rectangular_trimmed_basis();
    while let Some(b) = cur {
        cur = b.rectangular_trimmed_basis();
        deepest = Some(b);
    }
    deepest
}

