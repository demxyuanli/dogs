//! `AppParCurves_LeastSquare::Affect` and the two `CheckTangents` helpers
//! (`AppParCurves_LeastSquare.gxx:52-135`, `:1056-1206`).
//!
//! `Affect` reads the tangency (and, for the curvature constraint, the
//! curvature) of one end column out of the multi-line and may downgrade the
//! constraint when the tool cannot supply it.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::MathVector;

use super::*;

/// `CheckTangents(Pt1, Pt2, Tg3d, Tg2d)` (`gxx:64-107`): reverses every 3d
/// and 2d tangent when any 3d tangent points against `Pt2 - Pt1`.
fn check_tangents_3d(pts1: &[GpPnt], pts2: &[GpPnt], tg3d: &mut [GpVec], tg2d: &mut [GpVec2d]) {
    if pts1.len() != pts2.len() {
        return;
    }
    if tg3d.len() != pts1.len() {
        return;
    }

    let mut is_to_change_dir = false;
    for i in 0..pts1.len() {
        let v1 = GpVec::from_pnts(&pts1[i], &pts2[i]);
        if v1.dot(&tg3d[i]) < 0.0 {
            is_to_change_dir = true;
            break;
        }
    }

    if !is_to_change_dir {
        return;
    }

    for v in tg3d.iter_mut() {
        v.reverse();
    }
    for v in tg2d.iter_mut() {
        v.reverse();
    }
}

/// `CheckTangents(Pt1, Pt2, Tg2d)` (`gxx:115-136`): reverses the individual
/// 2d tangents that point against `Pt2 - Pt1`.
fn check_tangents_2d(pts1: &[GpPnt2d], pts2: &[GpPnt2d], tg2d: &mut [GpVec2d]) {
    if pts1.len() != pts2.len() {
        return;
    }

    for i in 0..pts1.len() {
        let v1 = GpVec2d::new(pts2[i].x() - pts1[i].x(), pts2[i].y() - pts1[i].y());
        if v1.dot(&tg2d[i]) < 0.0 {
            tg2d[i].reverse();
        }
    }
}

/// `Affect(SSP, Index, Cons, Vt, Vc)` (`gxx:1056-1206`).
///
/// `cons` is an in/out parameter exactly as in OCCT: `CurvaturePoint` may
/// fall back to `TangencyPoint`, and `TangencyPoint` to `PassPoint`, when the
/// tool has no such derivative at this column. `nb_p` / `nb_p2d` are the
/// `ToolLine::NbP3d` / `NbP2d` of the multi-line.
pub(super) fn affect<L: LineTool>(
    ml: &L,
    index: i32,
    cons: &mut Constraint,
    vt: &mut MathVector,
    vc: &mut MathVector,
    nb_p: i32,
    nb_p2d: i32,
) {
    // Vt: vector of tangent, Vc: vector of curvature.
    if !at_least_tangency(*cons) {
        return;
    }

    let mynb_p2d = if nb_p2d == 0 { 1 } else { nb_p2d };
    let mynb_p = if nb_p == 0 { 1 } else { nb_p };

    if *cons == Constraint::CurvaturePoint {
        let mut tab_v = vec![GpVec::zero(); mynb_p.max(0) as usize];
        let mut tab_v2d = vec![GpVec2d::new(0.0, 0.0); mynb_p2d.max(0) as usize];
        let found;
        if nb_p != 0 && nb_p2d != 0 {
            match ml.curvature(index) {
                Some((v, v2d)) => {
                    tab_v = v;
                    tab_v2d = v2d;
                    found = true;
                }
                None => found = false,
            }
        } else if nb_p2d != 0 {
            match ml.curvature_2d(index) {
                Some(v2d) => {
                    tab_v2d = v2d;
                    found = true;
                }
                None => found = false,
            }
        } else {
            match ml.curvature_3d(index) {
                Some(v) => {
                    tab_v = v;
                    found = true;
                }
                None => found = false,
            }
        }
        if !found {
            *cons = Constraint::TangencyPoint;
        }
        if found {
            let mut i2 = 1;
            for i in 1..=nb_p {
                let v = tab_v[(i - 1) as usize];
                vc.set_value(i2 as usize, v.x());
                vc.set_value((i2 + 1) as usize, v.y());
                vc.set_value((i2 + 2) as usize, v.z());
                i2 += 3;
            }
            for i in 1..=nb_p2d {
                let v = tab_v2d[(i - 1) as usize];
                vc.set_value(i2 as usize, v.x());
                vc.set_value((i2 + 1) as usize, v.y());
                i2 += 2;
            }
        }
    }

    if !at_least_tangency(*cons) {
        return;
    }

    let mut tab_v = vec![GpVec::zero(); mynb_p.max(0) as usize];
    let mut tab_v2d = vec![GpVec2d::new(0.0, 0.0); mynb_p2d.max(0) as usize];
    let found;
    if nb_p != 0 && nb_p2d != 0 {
        match ml.tangency(index) {
            Some((v, v2d)) => {
                tab_v = v;
                tab_v2d = v2d;
                found = true;
            }
            None => found = false,
        }
    } else if nb_p2d != 0 {
        match ml.tangency_2d(index) {
            Some(v2d) => {
                tab_v2d = v2d;
                found = true;
            }
            None => found = false,
        }
    } else {
        match ml.tangency_3d(index) {
            Some(v) => {
                tab_v = v;
                found = true;
            }
            None => found = false,
        }
    }

    if !found {
        *cons = Constraint::PassPoint;
    }

    if found {
        if nb_p != 0 {
            let (pts1, pts2) = if index < ml.last_point() {
                (
                    ml.value_3d(index),
                    ml.value_3d(index + 1),
                )
            } else {
                // (Index == ToolLine::LastPoint(theML))
                (
                    ml.value_3d(index - 1),
                    ml.value_3d(index),
                )
            };
            check_tangents_3d(&pts1, &pts2, &mut tab_v, &mut tab_v2d);
        } else if nb_p2d != 0 {
            let (pts1, pts2): (Vec<GpPnt2d>, Vec<GpPnt2d>) =
                if index < ml.last_point() {
                    let (_, a) = ml.value(index);
                    let (_, b) = ml.value(index + 1);
                    (a, b)
                } else {
                    // (Index == ToolLine::LastPoint(theML))
                    let (_, a) = ml.value(index - 1);
                    let (_, b) = ml.value(index);
                    (a, b)
                };
            check_tangents_2d(&pts1, &pts2, &mut tab_v2d);
        }

        let mut i2 = 1;
        for i in 1..=nb_p {
            let v = tab_v[(i - 1) as usize];
            vt.set_value(i2 as usize, v.x());
            vt.set_value((i2 + 1) as usize, v.y());
            vt.set_value((i2 + 2) as usize, v.z());
            i2 += 3;
        }
        for i in 1..=nb_p2d {
            let v = tab_v2d[(i - 1) as usize];
            vt.set_value(i2 as usize, v.x());
            vt.set_value((i2 + 1) as usize, v.y());
            i2 += 2;
        }
    }
}




