use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_geom::Curve;

use crate::int_tools_curve_box::add_curve_to_box;
use crate::inttools_data::IntRange;

// ---------------------------------------------------------------------------
// Range / box helpers
// ---------------------------------------------------------------------------

/// `SplitRangeOnSegments` (`IntTools_EdgeEdge.cxx:1366-1406`).
pub(in crate::edge_edge) fn split_range_on_segments(
    t1: f64,
    t2: f64,
    resolution: f64,
    nb_seg: i32,
) -> Vec<IntRange> {
    let diff = t2 - t1;
    if diff < resolution || nb_seg == 1 {
        return vec![IntRange::new_unchecked(t1, t2)];
    }
    let mut nb_segments = nb_seg;
    let mut dt = diff / nb_segments as f64;
    if dt < resolution {
        let seg = diff / resolution;
        nb_segments = seg as i32 + 1;
        dt = diff / nb_segments as f64;
    }
    let mut out = Vec::new();
    let mut t1x = t1;
    for _ in 1..nb_segments {
        let t2x = t1x + dt;
        out.push(IntRange::new_unchecked(t1x, t2x));
        t1x = t2x;
    }
    out.push(IntRange::new_unchecked(t1x, t2));
    out
}

/// `BndBuildBox` (`IntTools_EdgeEdge.cxx:1410-1419`) — the OCCT call is
/// `BndLib_Add3dCurve::Add(theBAC, aT1, aT2, theTol, aB)`.
pub(in crate::edge_edge) fn bnd_build_box(c: &dyn Curve, t1: f64, t2: f64, tol: f64) -> BndBox {
    let mut b = BndBox::new();
    add_curve_to_box(c, t1, t2, tol, &mut b);
    b
}

/// `PointBoxDistance` (`IntTools_EdgeEdge.cxx:1423-1452`).
///
/// `Bnd_Box::Get` already adds the box gap, exactly as `BndBox::distance`
/// does, so this is the same formula. A void box is `+inf` here (OCCT's `Get`
/// on a void box yields `RealLast` coordinates, i.e. a huge finite distance).
pub(in crate::edge_edge) fn point_box_distance(b: &BndBox, p: &GpPnt) -> f64 {
    b.distance(p)
}
