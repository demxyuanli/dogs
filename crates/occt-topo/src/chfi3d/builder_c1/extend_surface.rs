//! `ChFi3d_ExtendSurface`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:421-449`
//! (OCCT 8.0.0).
//!
//! ```text
//! if (prol) return;
//! prol = S->IsKind(Geom_BSplineSurface) ? 1 : S->IsKind(Geom_BezierSurface) ? 2 : 0;
//! if (!prol) return;
//! S->Bounds(umin, umax, vmin, vmax);
//! S->D0(umin, vmin, P1); S->D0(umax, vmax, P2);
//! length = P1.Distance(P2);
//! aBS = down_cast<Geom_BoundedSurface>(S);
//! GeomLib::ExtendSurfByLength(aBS, length, 1, false, true);
//! GeomLib::ExtendSurfByLength(aBS, length, 1, true,  true);
//! GeomLib::ExtendSurfByLength(aBS, length, 1, false, false);
//! GeomLib::ExtendSurfByLength(aBS, length, 1, true,  false);
//! S = aBS;
//! ```
//!
//! `prol` is the in/out `int` that records which arm ran (`0` = not a bounded
//! surface, `1` = B-spline, `2` = Bezier); the four extensions grow the surface
//! by its own diagonal length at `C1` on every side.
//!
//! Divergences that are not heuristics:
//! - OCCT mutates the `Geom_Surface` handle in place; the port takes the surface
//!   by `Arc` and returns the (possibly new) surface, chaining each
//!   [`occt_geom::geom_lib_extend_surf::extend_surf_by_length`] result into the
//!   next call - the same four extensions of the same object,
//! - `S->IsKind(STANDARD_TYPE(Geom_BSplineSurface))` is
//!   [`occt_geom::surface::Surface::is_bspline_surface`] and
//!   `IsKind(STANDARD_TYPE(Geom_BezierSurface))` is
//!   [`occt_geom::surface::Surface::is_bezier_surface`]; both are false for a
//!   `Geom_RectangularTrimmedSurface`, as `IsKind` is,
//! - an arm `extend_surf_by_length` reports as `None` (its own
//!   `default: break`) leaves the surface untouched and the remaining three
//!   calls still run, exactly like a no-op `ExtendSurfByLength`.

use std::sync::Arc;

use occt_geom::geom_lib_extend_surf::extend_surf_by_length;
use occt_geom::Surface;

/// `ChFi3d_ExtendSurface(S, prol)` (`cxx:421-449`): extends a B-spline /
/// Bezier surface by its diagonal length on all four sides and answers the
/// surface to keep using. `prol` is updated in place.
pub fn ch_fi3d_extend_surface(s: &Arc<dyn Surface>, prol: &mut i32) -> Arc<dyn Surface> {
    if *prol != 0 {
        return Arc::clone(s);
    }

    *prol = if s.is_bspline_surface() {
        1
    } else if s.is_bezier_surface() {
        2
    } else {
        0
    };
    if *prol == 0 {
        return Arc::clone(s);
    }

    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    let p1 = s.d0(umin, vmin);
    let p2 = s.d0(umax, vmax);
    let length = p1.distance(&p2);

    let mut current = Arc::clone(s);
    // `GeomLib::ExtendSurfByLength(aBS, length, 1, InU, After)` x4.
    for (in_u, after) in [(false, true), (true, true), (false, false), (true, false)] {
        if let Some(extended) = extend_surf_by_length(&current, length, 1, in_u, after) {
            current = extended;
        }
    }
    current
}
