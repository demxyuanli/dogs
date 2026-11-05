//! The four analytical fillet cases of `ChFi2d_AnaFilletAlgo.cxx` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `SegmentFilletSegment` (`cxx:584-641`) | 58 | [`segment_fillet_segment`] |
//! | `SegmentFilletArc` (`cxx:643-788`) | 146 | [`segment_fillet_arc`] |
//! | `ArcFilletSegment` (`cxx:790-936`) | 147 | [`arc_fillet_segment`] |
//! | `ArcFilletArc` (`cxx:938-1048`) | 111 | [`arc_fillet_arc`] |
//!
//! Each function fills the `FilletDef` the C++ code passes as a reference
//! parameter list and answers the C++ `bool`. Everything is 2D: the caller
//! projected the two edges onto the working plane.
//!
//! Notes on the differences between the four bodies, all preserved from
//! OCCT rather than unified:
//! - `SegmentFilletArc` answers the point of the second circle with
//!   `fabs(d2 - Rf2) < Precision::Confusion()` (`cxx:778`), while
//!   `ArcFilletSegment` uses `Precision::SquareConfusion()` (`cxx:926`),
//! - `ArcFilletArc` never normalizes its bisectrix (`cxx:982`), and its
//!   `checkl` (`cxx:988`) is built and never read,
//! - the shrink of the second neighbour in `Perform` tests `e1`'s orientation
//!   (`cxx:538`), not `e2`'s; that quirk lives in `mod.rs`.

use std::f64::consts::PI;

use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::intana2d::IntAna2dAnaIntersection;
use occt_core::precision::{ANGULAR, CONFUSION, REAL_SMALL, SQUARE_CONFUSION};

use super::fillet_support::{is_equal_2d};
use super::{FilletDef, Neighbour};

/// Column vector of the X axis, `gp::DX2d()`.
fn dx2d() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("gp::DX2d")
}

/// `gp_Vec2d(from, to)` (`gp_Vec2d.hxx`): `to - from`.
fn vec2d(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `gp_Lin2d(P, V)` (`gp_Lin2d.hxx`): origin plus direction, normalizing `V`.
fn lin2d(pt: &GpPnt2d, dir: &GpVec2d) -> Option<GpLin2d> {
    GpDir2d::from_vec2d(dir).ok().map(|d| GpLin2d::from_pnt_dir(*pt, d))
}

/// `bool ChFi2d_AnaFilletAlgo::SegmentFilletSegment(radius, xc, yc, cw, start,
/// end)` (`cxx:588-641`).
pub(super) fn segment_fillet_segment(
    left: &Neighbour,
    right: &Neighbour,
    radius: f64,
    def: &mut FilletDef,
) -> bool {
    // `cxx:594-596`.
    let p11 = GpPnt2d::new(left.x1, left.y1);
    let p12 = GpPnt2d::new(left.x2, left.y2);
    let p22 = GpPnt2d::new(right.x2, right.y2);

    // `cxx:600-603`.
    if is_equal_2d(&p12, &p11) || is_equal_2d(&p12, &p22) {
        return false;
    }

    // `cxx:606-609`: `gp_Vec2d v1(p12, p11); gp_Vec2d v2(p12, p22);` then both
    // normalized.
    let mut v1 = vec2d(&p12, &p11);
    let mut v2 = vec2d(&p12, &p22);
    let _ = v1.normalize();
    let _ = v2.normalize();

    // `cxx:612`: `gp_Vec2d bisec = 0.5 * (v1 + v2);`.
    let mut bisec = v1.added(&v2).multiplied_scalar(0.5);

    // `cxx:615-618`.
    if bisec.square_magnitude() < SQUARE_CONFUSION {
        return false;
    }

    // `cxx:621`.
    let _ = bisec.normalize();

    // `cxx:624`: the angle at the bisectrix.
    let beta = v1.angle(&bisec);

    // `cxx:627`: length along the bisectrix till the centre of the fillet.
    let l = radius / beta.abs().sin();

    // `cxx:630-631`: `gp_Pnt2d pc = p12.Translated(L * bisec); pc.Coord(xc, yc);`.
    let pc = p12.translated_vec(&bisec.multiplied_scalar(l));
    def.xc = pc.x();
    def.yc = pc.y();

    // `cxx:634-635`.
    def.start = (l * l - radius * radius).sqrt();
    def.end = def.start;

    // `cxx:638`.
    def.cw = beta > 0.0;
    true
}

/// `bool ChFi2d_AnaFilletAlgo::SegmentFilletArc(radius, xc, yc, cw, start, end,
/// xend, yend)` (`cxx:646-788`).
pub(super) fn segment_fillet_arc(
    left: &Neighbour,
    right: &Neighbour,
    radius: f64,
    def: &mut FilletDef,
) -> bool {
    // `cxx:657-659`.
    let p12 = GpPnt2d::new(left.x2, left.y2);
    let p11 = GpPnt2d::new(left.x1, left.y1);
    let pc2 = GpPnt2d::new(right.xc, right.yc);

    // `cxx:662-665`.
    if p11.square_distance(&p12) < REAL_SMALL {
        return false;
    }

    // `cxx:668-669`.
    let mut v1 = vec2d(&p12, &p11);
    let mut v2 = vec2d(&p12, &pc2);

    // `cxx:672-678`: rotate the arc vector to become tangential at p21.
    if right.cw {
        v2.rotate(PI / 2.0);
    } else {
        v2.rotate(-PI / 2.0);
    }

    // `cxx:682-686`.
    let mut angle = v1.angle(&v2);
    if angle.abs() < ANGULAR {
        return false;
    }

    // `cxx:689-691`.
    let _ = v2.normalize();
    let _ = v1.normalize();
    let mut bisec = v1.added(&v2).multiplied_scalar(0.5);

    // `cxx:695-698`.
    if bisec.square_magnitude() < REAL_SMALL {
        return false;
    }

    // `cxx:701-702`.
    let _ = bisec.normalize();
    let nearp = p12.translated_vec(&bisec.multiplied_scalar(radius));
    let Some(nearl) = lin2d(&p12, &bisec) else {
        return false;
    };

    // `cxx:706-711`: a line parallel to the segment through the "near" point.
    let mut d1 = v1;
    let Some(mut line) = lin2d(&p11, &d1.reversed()) else {
        return false;
    };
    d1.rotate(PI / 2.0);
    line.translate_vec(&d1.multiplied_scalar(radius));
    if line.distance(&nearp) > radius {
        line.translate_vec(&d1.multiplied_scalar(-2.0 * radius));
    }

    // `cxx:716-722`: a circle of radius `radius2 +/- radius`.
    let mut circ = GpCirc2d::new(GpAx22d::from_xdir(pc2, dx2d()), right.radius + radius);
    if right.radius > radius && circ.distance(&nearp) > radius {
        circ.set_radius(right.radius - radius);
    }

    // `cxx:725-729`.
    let mut intersector = IntAna2dAnaIntersection::new();
    intersector.perform_lin_circ(&line, &circ);
    if !intersector.is_done() || intersector.nb_points() == 0 {
        return false;
    }

    // `cxx:732-745`: the intersection point closest to the "near" line.
    let mut min_dist = f64::MAX;
    let mut xc = 0.0;
    let mut yc = 0.0;
    for i in 1..=intersector.nb_points() {
        let p = *intersector.point(i).value();
        let d = nearl.distance(&p);
        if d < min_dist {
            min_dist = d;
            xc = p.x();
            yc = p.y();
        }
    }

    // `cxx:748-751`: shrink of the segment.
    let pc = GpPnt2d::new(xc, yc);
    def.xc = xc;
    def.yc = yc;
    let l2 = pc.square_distance(&p12);
    let rf2 = radius * radius;
    def.start = (l2 - rf2).sqrt();

    // `cxx:754-755`: shrink of the arc.
    let pcc = vec2d(&pc2, &pc);
    def.end = vec2d(&pc2, &p12).angle(&pcc).abs();

    // `cxx:759-766`: the point on the arc that coincides with the fillet end.
    line.pos.set_location(pc2);
    match GpDir2d::from_vec2d(&pcc) {
        Ok(d) => line.pos.set_direction(d),
        Err(_) => return false,
    }
    circ.set_location(pc2);
    circ.set_radius(right.radius);
    intersector.perform_lin_circ(&line, &circ);
    if !intersector.is_done() || intersector.nb_points() == 0 {
        return false;
    }

    // `cxx:769-782`.
    def.xend = f64::MAX;
    def.yend = f64::MAX;
    for i in 1..=intersector.nb_points() {
        let p = *intersector.point(i).value();
        let d2 = p.square_distance(&pc);
        if (d2 - rf2).abs() < CONFUSION {
            def.xend = p.x();
            def.yend = p.y();
            break;
        }
    }

    // `cxx:785-786`.
    angle = v1.angle(&v2);
    def.cw = angle > 0.0;
    true
}

/// `bool ChFi2d_AnaFilletAlgo::ArcFilletSegment(radius, xc, yc, cw, start, end,
/// xstart, ystart)` (`cxx:793-936`).
pub(super) fn arc_fillet_segment(
    left: &Neighbour,
    right: &Neighbour,
    radius: f64,
    def: &mut FilletDef,
) -> bool {
    // `cxx:803-805`.
    let p12 = GpPnt2d::new(left.x2, left.y2);
    let p22 = GpPnt2d::new(right.x2, right.y2);
    let pc1 = GpPnt2d::new(left.xc, left.yc);

    // `cxx:809-812`.
    if p12.square_distance(&p22) < REAL_SMALL {
        return false;
    }

    // `cxx:815-816`.
    let mut v1 = vec2d(&p12, &pc1);
    let mut v2 = vec2d(&p12, &p22);

    // `cxx:819-826`: rotate the arc vector to become tangential at p21.
    if left.cw {
        v1.rotate(-PI / 2.0);
    } else {
        v1.rotate(PI / 2.0);
    }

    // `cxx:830-834`.
    let mut angle = v1.angle(&v2);
    if angle.abs() < ANGULAR {
        return false;
    }

    // `cxx:837-839`.
    let _ = v1.normalize();
    let _ = v2.normalize();
    let mut bisec = v1.added(&v2).multiplied_scalar(0.5);

    // `cxx:843-846`.
    if bisec.square_magnitude() < REAL_SMALL {
        return false;
    }

    // `cxx:849-850`.
    let _ = bisec.normalize();
    let near_point = p12.translated_vec(&bisec.multiplied_scalar(radius));
    let Some(near_line) = lin2d(&p12, &bisec) else {
        return false;
    };

    // `cxx:854-861`: a line parallel to the segment through the "near" point.
    let mut a_d2_vec = v2;
    let Some(mut line) = lin2d(&p22, &a_d2_vec.reversed()) else {
        return false;
    };
    a_d2_vec.rotate(PI / 2.0);
    line.translate_vec(&a_d2_vec.multiplied_scalar(radius));
    if line.distance(&near_point) > radius {
        line.translate_vec(&a_d2_vec.multiplied_scalar(-2.0 * radius));
    }

    // `cxx:864-871`: a circle of radius `radius1 +/- radius`.
    let mut circ = GpCirc2d::new(GpAx22d::from_xdir(pc1, dx2d()), left.radius + radius);
    if left.radius > radius && circ.distance(&near_point) > radius {
        circ.set_radius(left.radius - radius);
    }

    // `cxx:874-878`.
    let mut intersector = IntAna2dAnaIntersection::new();
    intersector.perform_lin_circ(&line, &circ);
    if !intersector.is_done() || intersector.nb_points() == 0 {
        return false;
    }

    // `cxx:881-893`: the intersection point closest to the "near" line.
    let mut min_dist = f64::MAX;
    let mut xc = 0.0;
    let mut yc = 0.0;
    for i in 1..=intersector.nb_points() {
        let p = *intersector.point(i).value();
        let d = near_line.distance(&p);
        if d < min_dist {
            min_dist = d;
            xc = p.x();
            yc = p.y();
        }
    }

    // `cxx:896-899`: shrink of the segment.
    let pc = GpPnt2d::new(xc, yc);
    def.xc = xc;
    def.yc = yc;
    let l2 = pc.square_distance(&p12);
    let rf2 = radius * radius;
    def.end = (l2 - rf2).sqrt();

    // `cxx:902-903`: shrink of the arc.
    let pcc = vec2d(&pc1, &pc);
    def.start = vec2d(&pc1, &p12).angle(&pcc).abs();

    // `cxx:907-914`: the point on the arc that coincides with the fillet start.
    line.pos.set_location(pc1);
    match GpDir2d::from_vec2d(&pcc) {
        Ok(d) => line.pos.set_direction(d),
        Err(_) => return false,
    }
    circ.set_location(pc1);
    circ.set_radius(left.radius);
    intersector.perform_lin_circ(&line, &circ);
    if !intersector.is_done() || intersector.nb_points() == 0 {
        return false;
    }

    // `cxx:917-930`.
    def.xstart = f64::MAX;
    def.ystart = f64::MAX;
    for i in 1..=intersector.nb_points() {
        let p = *intersector.point(i).value();
        let d2 = p.square_distance(&pc);
        if (d2 - rf2).abs() < SQUARE_CONFUSION {
            def.xstart = p.x();
            def.ystart = p.y();
            break;
        }
    }

    // `cxx:933-934`.
    angle = v2.angle(&v1);
    def.cw = angle < 0.0;
    true
}

/// `bool ChFi2d_AnaFilletAlgo::ArcFilletArc(radius, xc, yc, cw, start, end)`
/// (`cxx:944-1048`).
pub(super) fn arc_fillet_arc(
    left: &Neighbour,
    right: &Neighbour,
    radius: f64,
    def: &mut FilletDef,
) -> bool {
    // `cxx:951-953`.
    let pc1 = GpPnt2d::new(left.xc, left.yc);
    let pc2 = GpPnt2d::new(right.xc, right.yc);
    let p12 = GpPnt2d::new(left.x2, left.y2);

    // `cxx:956-957`.
    let mut v1 = vec2d(&pc1, &p12);
    let mut v2 = vec2d(&pc2, &p12);

    // `cxx:960-974`: rotate the vectors so that they are tangent at p12.
    if left.cw {
        v1.rotate(PI / 2.0);
    } else {
        v1.rotate(-PI / 2.0);
    }
    if right.cw {
        v2.rotate(-PI / 2.0);
    } else {
        v2.rotate(PI / 2.0);
    }

    // `cxx:978-984`. The bisectrix stays unnormalized here (see the module
    // notes), so `checkp` is only a distance reference.
    let _ = v1.normalize();
    let _ = v2.normalize();
    let bisec = v1.added(&v2).multiplied_scalar(0.5);
    if bisec.square_magnitude() < REAL_SMALL {
        return false;
    }

    // `cxx:988` builds `checkl(p12, bisec)` and never reads it.
    let checkp = p12.translated_vec(&bisec.multiplied_scalar(radius));

    // `cxx:992-997`: arc 1.
    let mut c1 = GpCirc2d::new(GpAx22d::from_xdir(pc1, dx2d()), left.radius + radius);
    if left.radius > radius && c1.distance(&checkp) > radius {
        c1.set_radius(left.radius - radius);
    }
    // `cxx:998-1005`: arc 2.
    let mut c2 = GpCirc2d::new(GpAx22d::from_xdir(pc2, dx2d()), right.radius + radius);
    if right.radius > radius && c2.distance(&checkp) > radius {
        c2.set_radius(right.radius - radius);
    }

    // `cxx:1009-1012`.
    let mut intersector = IntAna2dAnaIntersection::new();
    intersector.perform_circ_circ(&c1, &c2);
    if !intersector.is_done() || intersector.nb_points() == 0 {
        return false;
    }

    // `cxx:1015-1030`: the intersection point closest to the "check" point.
    let mut pc = GpPnt2d::zero();
    let mut min_dist = f64::MAX;
    for i in 1..=intersector.nb_points() {
        let p = *intersector.point(i).value();
        let d = checkp.square_distance(&p);
        if d < min_dist {
            min_dist = d;
            pc = p;
        }
    }
    def.xc = pc.x();
    def.yc = pc.y();

    // `cxx:1033-1043`: orientation of the fillet.
    let angle = v1.angle(&v2);
    if angle.abs() < ANGULAR {
        let a = vec2d(&pc, &pc1).angle(&vec2d(&pc, &pc2));
        def.cw = a < 0.0;
    } else {
        def.cw = angle > 0.0;
    }

    // `cxx:1046-1047`.
    def.start = vec2d(&pc1, &p12).angle(&vec2d(&pc1, &pc)).abs();
    def.end = vec2d(&pc2, &p12).angle(&vec2d(&pc2, &pc)).abs();
    true
}
