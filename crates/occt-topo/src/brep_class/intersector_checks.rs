//! Static helpers of `BRepClass_Intersector.cxx`: `CheckOn`, `CheckSkip`,
//! `RefineTolerance`, `GetTangentAsChord`, `IsInter`, `MaxTol2DCurEdge` and
//! the `ElCLib` line evaluation / parameter pair they use.

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dTransition,
};
use occt_core::precision::{Precision, CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::extrema2d::point_curve_extrema2d_all;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::int_tools_wline::{u_resolution, v_resolution};
use crate::shape::{Face, Vertex};

use super::edge::{vertices_cum_ori, BRepClassEdge};

/// `Precision::SquarePConfusion()`.
const SQUARE_PCONFUSION: f64 = PCONFUSION * PCONFUSION;

/// `ElCLib::Value(U, L)` for `gp_Lin2d` (`ElCLib.cxx`).
pub(super) fn lin_value(l: &GpLin2d, u: f64) -> GpPnt2d {
    let loc = l.location();
    let dir = l.direction();
    GpPnt2d::new(loc.x() + u * dir.x(), loc.y() + u * dir.y())
}

/// `ElCLib::Parameter(L, P)` for `gp_Lin2d` (`ElCLib.cxx`): the projection
/// of `P - Location` on the line direction.
pub(super) fn lin_parameter(l: &GpLin2d, p: &GpPnt2d) -> f64 {
    let loc = l.location();
    let dir = l.direction();
    (p.x() - loc.x()) * dir.x() + (p.y() - loc.y()) * dir.y()
}

/// `Geom2dAdaptor_Curve(C, UFirst, ULast)`: a view of `curve` restricted to
/// `[first, last]`. The trimmed-curve constructor unwraps an existing
/// trimmed basis, as `Geom2dAdaptor_Curve::Load` does, and does not check the
/// range against the basis (`adjust_periodic = false`).
pub(super) fn adaptor_curve(curve: &Arc<dyn Curve2d>, first: f64, last: f64) -> Arc<dyn Curve2d> {
    Arc::new(Geom2dTrimmedCurve::new_sense(
        curve.clone(),
        first,
        last,
        true,
        false,
    ))
}

/// `MaxTol2DCurEdge` (`BRepClass_Intersector.cxx:72-117`).
pub(super) fn max_tol_2d_cur_edge(
    v1: Option<&Vertex>,
    v2: Option<&Vertex>,
    face: &Face,
    tol: f64,
) -> f64 {
    let tol_v1 = v1.map_or(0.0, BRepTool::vertex_tolerance);
    let tol_v2 = v2.map_or(0.0, BRepTool::vertex_tolerance);
    let tol_v3d = tol_v1.max(tol_v2);
    let Some(surf) = BRepTool::face_surface(face) else {
        return tol;
    };
    // `BRepAdaptor_Surface aS(theF, false)`: UResolution / VResolution.
    let u_res = u_resolution(surf.as_ref(), tol_v3d);
    let v_res = v_resolution(surf.as_ref(), tol_v3d);
    u_res.max(v_res).max(tol)
}

/// `IsInter` (`BRepClass_Intersector.cxx:124-137`).
pub(super) fn is_inter(bond: &BndBox2d, l: &GpLin2d, p: f64) -> bool {
    if Precision::is_infinite(p) {
        !bond.is_out_lin(l)
    } else {
        let pnt_f = l.location();
        let pnt_l = lin_value(l, p);
        !bond.is_out_segment(&pnt_f, &pnt_l)
    }
}

/// `GetTangentAsChord` (`BRepClass_Intersector.cxx:520-550`). `None` when the
/// chord is shorter than `SquarePConfusion` (OCCT leaves the output
/// unchanged).
pub(super) fn get_tangent_as_chord(
    pcurve: &dyn Curve2d,
    param: f64,
    first: f64,
    last: f64,
) -> Option<GpDir2d> {
    let mut offset = 0.1 * (last - first);
    if last - param < PCONFUSION {
        // `theParam == theLast`
        offset *= -1.0;
    } else if param + offset > last {
        // `theParam` is close to `theLast`
        offset = 0.5 * (last - param);
    }
    let p0 = pcurve.d0(param);
    let p1 = pcurve.d0(param + offset);
    let mut chord = GpVec2d::new(p1.x() - p0.x(), p1.y() - p0.y());
    if offset < 0.0 {
        chord.reverse();
    }
    if chord.square_magnitude() > SQUARE_PCONFUSION {
        GpDir2d::from_vec2d(&chord).ok()
    } else {
        None
    }
}

/// `RefineTolerance` (`BRepClass_Intersector.cxx:478-515`): for cylindrical
/// faces the tolerance is tightened by the parametric resolution along the
/// curve tangent at `t`.
pub(super) fn refine_tolerance(face: &Face, cur: &dyn Curve2d, t: f64, tol_z: &mut f64) {
    let Some(surf) = BRepTool::face_surface(face) else {
        return;
    };
    if classify_surface(surf.as_ref()) != SurfaceKind::Cylinder {
        return;
    }
    let u_res = u_resolution(surf.as_ref(), *tol_z);
    let v_res = v_resolution(surf.as_ref(), *tol_z);
    let (_p, v2d) = cur.d1(t);
    let dir = GpDir2d::from_vec2d(&v2d).expect("RefineTolerance: null tangent");
    let mut tol_x = u_res * dir.y() + v_res * dir.x();
    if tol_x < 0.0 {
        tol_x = -tol_x;
    }
    if tol_x < CONFUSION {
        tol_x = CONFUSION;
    }
    if tol_x < *tol_z {
        *tol_z = tol_x;
    }
}

/// `CheckOn` (`BRepClass_Intersector.cxx:140-198`): the nearest point of the
/// curve to the line is a hit when it lies within `tol_z`. The distance
/// search is `Extrema_ExtPC2d`, here `point_curve_extrema2d_all`.
pub(super) fn check_on(
    face: &Face,
    l: &GpLin2d,
    cur: &dyn Curve2d,
    tol_z: &mut f64,
    fin: f64,
    deb: f64,
) -> Option<IntRes2dIntersectionPoint> {
    let extrema = point_curve_extrema2d_all(cur, &l.location());
    // `RealLast()` start; strict `<` keeps the first minimum as in OCCT.
    let mut min_dist = f64::MAX;
    let mut min_ind: Option<usize> = None;
    for (i, ext) in extrema.iter().enumerate() {
        if ext.distance < min_dist {
            min_dist = ext.distance;
            min_ind = Some(i);
        }
    }
    let idx = min_ind?;
    if min_dist > *tol_z {
        return None;
    }
    let ext = &extrema[idx];
    // `Point(i).Value()` is the curve point; `Parameter()` its curve
    // parameter (`pair2d_c` stores it in `u2`).
    let pnt_exact = ext.p2;
    let par = ext.u2;
    refine_tolerance(face, cur, par, tol_z);
    if min_dist > *tol_z {
        return None;
    }
    let tr_on_lin = IntRes2dTransition::undecided(IntRes2dPosition::Head);
    let pos_on_curve = if (par - deb).abs() <= CONFUSION || par < deb {
        IntRes2dPosition::Head
    } else if (par - fin).abs() <= CONFUSION || par > fin {
        IntRes2dPosition::End
    } else {
        IntRes2dPosition::Middle
    };
    let tr_on_curve = IntRes2dTransition::undecided(pos_on_curve);
    Some(IntRes2dIntersectionPoint::with_transitions(
        &pnt_exact,
        0.0,
        par,
        &tr_on_lin,
        &tr_on_curve,
        false,
    ))
}

/// `CheckSkip` (`BRepClass_Intersector.cxx:201-325`). When the line misses
/// the edge because it passes a vertex of high tolerance, the gap between
/// this edge's last point and the next edge's first point is intersected
/// with the line instead. `deb` / `fin` are the edge range from `Perform`.
/// Returns the replacement intersection, or `None` when nothing applies.
#[allow(clippy::too_many_arguments)]
pub(super) fn check_skip(
    l: &GpLin2d,
    gl: &dyn Curve2d,
    e: &BRepClassEdge,
    c2d: &dyn Curve2d,
    dl: &IntRes2dDomain,
    deb: f64,
    fin: f64,
    max_tol: f64,
) -> Option<Geom2dIntGInter> {
    let (Some(edge), Some(face)) = (e.edge(), e.face()) else {
        return None;
    };
    // `TopExp::LastVertex(theE.Edge(), true)`.
    let vl = vertices_cum_ori(edge).1?;
    if !(BRepTool::vertex_tolerance(&vl) > max_tol) {
        return None;
    }
    let next = e.next_edge()?;
    let (lc2d, ldeb, lfin) = curve_on_surface_range(next, face)?;

    let (a, b, c) = l.coefficients();

    let mut at1 = fin;
    if edge.orientation() != Orientation::Forward {
        at1 = deb;
    }
    let mut at2 = ldeb;
    if next.orientation() != Orientation::Forward {
        at2 = lfin;
    }
    let p1 = c2d.d0(at1);
    let p2 = lc2d.d0(at2);

    // Both ends must lie strictly inside the DL domain.
    let par1 = lin_parameter(l, &p1);
    let par2 = lin_parameter(l, &p2);
    if par1 <= dl.first_parameter()
        || par1 >= dl.last_parameter()
        || par2 <= dl.first_parameter()
        || par2 >= dl.last_parameter()
    {
        return None;
    }

    let fv = a * p1.x() + b * p1.y() + c;
    let sv = a * p2.x() + b * p2.y() + c;
    // Same sign: the gap does not cross the line.
    if fv * sv >= 0.0 {
        return None;
    }

    // `GC_MakeSegment2d(P1, P2)`: fails when the points coincide within
    // `gp::Resolution()` (`RealSmall`).
    let dist = p1.distance(&p2);
    if dist <= f64::MIN_POSITIVE {
        return None;
    }
    let dir = GpDir2d::from_vec2d(&GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())).ok()?;
    let seg_line = Geom2dLine::from_pnt_dir(p1, dir);
    let skip: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(Arc::new(seg_line), 0.0, dist));

    // `theCur.Load(aSkipC2D)`: the range is the segment's own [0, dist].
    let skip_first = skip.first_parameter();
    let skip_last = skip.last_parameter();
    let pdeb = skip.d0(skip_first);
    let pfin = skip.d0(skip_last);
    let de = IntRes2dDomain::bounded(&pdeb, skip_first, 1.0e-5, &pfin, skip_last, 1.0e-5);

    let mut inter = Geom2dIntGInter::new();
    inter.perform(
        gl,
        dl,
        skip.as_ref(),
        &de,
        PCONFUSION,
        occt_core::precision::PINTERSECTION,
    );
    Some(inter)
}
