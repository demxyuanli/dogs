//! The generatrix / stop-surface intersection of `GeomFill_LocationDraft`.
//!
//! Source: `GeomFill_LocationDraft.cxx:168-169` (`IntCurveSurface_HInter
//! ::Perform(GeomAdaptor_Curve(Geom_Line), mySurf)`), `cxx:173-187` (the
//! minimum-`W` selection over `Int.NbPoints()`).
//!
//! Only the `GeomAbs_Plane` arm is ported: it is the exact path
//! `IntCurveSurface_Inter::PerformConicSurf` takes for a line and a plane,
//! whose `IntSurf_Quadric::Distance` is `n.(P - O)`
//! (`IntSurf_Quadric.cxx:173-179`). The generatrix is a `Geom_Line`, so that
//! distance is affine in the line parameter, its single root is
//! `-n.(P - O) / (n.d)`, and the surface parameters are the `gp_Pln` frame
//! coordinates `(u, v) = ((Q - O).XDir, (Q - O).YDir)`, i.e. exactly the
//! `U` / `V` that `IntCurveSurface_IntersectionPoint::SetValues` receives from
//! the plane adaptor.
//!
//! **UNPORTED**: every other surface type. `IntCurveSurface_HInter` itself is
//! ported in `occt_topo::intcurvesurface::perform_curve_surface`, which is the
//! higher crate; `occt-geom` cannot call it, and the conic / quadric arms are
//! not replicated here as approximate scans.

use occt_core::gp::{GpPnt, GpVec};

use crate::surface::Surface;

/// The `UNPORTED` message for a non-planar stop surface.
pub(crate) const UNPORTED_NON_PLANE: &str =
    "UNPORTED: GeomFill_LocationDraft.cxx:168 IntCurveSurface_HInter for a non-plane stop surface \
     lives in occt-topo::intcurvesurface and is unreachable from occt-geom";

/// `IntCurveSurface_HInter::Perform` for the lines this class builds: the
/// `(U, V, W)` triples of every intersection, `W` being the line parameter.
pub(crate) fn line_surface_intersections(
    p: &GpPnt,
    d: &GpVec,
    surf: &dyn Surface,
) -> Result<Vec<(f64, f64, f64)>, String> {
    let Some(pln) = surf.gp_pln() else {
        return Err(UNPORTED_NON_PLANE.to_string());
    };
    let pos = pln.position();
    let o = pos.location();
    let n = GpVec::from_xyz(pos.direction().xyz());
    let denom = n.dot(d);
    // `IntSurf_Quadric::Distance` is constant along the line: no root.
    if denom == 0.0 {
        return Ok(Vec::new());
    }
    let t = -GpVec::from_pnts(&o, p).dot(&n) / denom;
    let q = GpPnt::new(p.x() + t * d.coord.x, p.y() + t * d.coord.y, p.z() + t * d.coord.z);
    let rel = GpVec::from_pnts(&o, &q);
    let u = rel.dot(&GpVec::from_xyz(pos.x_direction().xyz()));
    let v = rel.dot(&GpVec::from_xyz(pos.y_direction().xyz()));
    Ok(vec![(u, v, t)])
}

/// The `W`-minimum point of `cxx:171-188`, or `None` when `Int.NbPoints() == 0`.
pub(crate) fn closest_intersection(points: &[(f64, f64, f64)]) -> Option<(f64, f64, f64)> {
    let mut best = *points.first()?;
    for q in &points[1..] {
        if best.2 > q.2 {
            best = *q;
        }
    }
    Some(best)
}
