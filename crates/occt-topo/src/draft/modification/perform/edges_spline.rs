//! `Draft_Modification::Perform`, the `Geom_BSplineCurve` arm of the
//! surface / surface edge reconstruction (`Draft_Modification_1.cxx:1159-1325`).
//!
//! When `GeomInt_IntSS` answers with B-spline lines there is no pcurve to reuse
//! and no single line to pick: the reference starts from the piece that begins
//! at the old edge's first vertex, walks the remaining pieces end to end and
//! concatenates them with `GeomConvert_CompCurveToBSplineCurve`.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | seam filter on `S1` (cylinder / cone) | `_1.cxx:1161-1208` | [`glue_bspline_lines`] |
//! | first piece, nearest to `pfv` | `_1.cxx:1211-1230` | [`glue_bspline_lines`] |
//! | end-to-end walk (`ToGlue`) | `_1.cxx:1232-1273` | [`glue_bspline_lines`] |
//! | `CompCurveToBSplineCurve` concatenation | `_1.cxx:1280-1290` | [`glue_bspline_lines`] |
//! | `Extrema_ExtPC` re-parameterisation, orientation | `_1.cxx:1292-1313` | [`glue_bspline_lines`] |
//!
//! # Translation notes
//!
//! * `GeomAPI_ProjectPointOnSurf(Pnt, S1, Precision::Confusion())` +
//!   `LowerDistanceParameters` (`_1.cxx:1170-1185`) is
//!   `geom_api::project_point_on_surface`; the reference never tests
//!   `IsDone()`, and a failed read there is undefined, so the port answers
//!   "not on the seam" instead of inventing a branch.
//! * `aCurve->Reverse()` (`_1.cxx:1187`, `:1265`, `:1313`) mutates the shared
//!   `Geom_Curve` handle. The port's curves are immutable, so the reversed copy
//!   replaces the entry in its own list, which is what the later lists see.
//! * `FirstCurve != i2s.Line(i)` (`_1.cxx:1237`) compares handles; with
//!   copies the port compares the intersection line's index instead.
//! * The `FirstCurve.IsNull()` test (`_1.cxx:1275-1278`) cannot fire: the
//!   candidate list is non-empty and the first candidate always beats
//!   `RealLast()` in the distance loop. The port keeps the emptiness test only.
//! * `myExtPC.IsMin(i)` (`_1.cxx:1297`) is `ExtPcSolution::is_min`.
//! * The pcurve assignment after the gluing is commented out in the reference
//!   (`_1.cxx:1314-1323`), so the port sets no pcurve here.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::convert_bspl::CompCurveToBSplineCurve;
use occt_geom::curve::Curve;
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::geom_int::IntSS;

/// `Draft_Modification_1.cxx:1159-1325`: glue the B-spline intersection lines
/// into one curve starting at the old edge's first vertex.
pub(super) fn glue_bspline_lines(
    i2s: &IntSS,
    pfv: &GpPnt,
    d1fv: &GpVec,
    s1: &Arc<dyn Surface>,
) -> Option<Arc<dyn Curve>> {
    // `cxx:1162`: the intersection lines, kept with their 1-based indices so
    // the "everything but the first piece" list of `cxx:1237` can be rebuilt.
    let mut lines: Vec<(i32, Arc<dyn Curve>)> = (1..=i2s.nb_lines())
        .filter_map(|i| i2s.line(i).map(|l| (i, l.curve.clone())))
        .collect();
    if lines.is_empty() {
        return None;
    }

    // `cxx:1161-1208`: on a cylinder or cone only the pieces that reach the
    // `U = 0` / `U = 2 pi` seam can start the chain.
    let mut candidates: Vec<(i32, Arc<dyn Curve>)> = Vec::new();
    if s1.gp_cylinder().is_some() || s1.gp_cone().is_some() {
        for (idx, curve) in lines.iter_mut() {
            if end_on_seam(s1.as_ref(), &curve.d0(curve.first_parameter())) {
                candidates.push((*idx, curve.clone()));
            } else if end_on_seam(s1.as_ref(), &curve.d0(curve.last_parameter())) {
                *curve = Arc::from(curve.reversed());
                candidates.push((*idx, curve.clone()));
            }
        }
        if candidates.is_empty() {
            // `cxx:1192-1200`: nothing starts on the seam, keep the order the
            // intersection produced.
            candidates = lines.clone();
        }
    } else {
        // `cxx:1202-1208`
        candidates = lines.clone();
    }

    // `cxx:1211-1230`
    let (first_idx, first_curve) = if candidates.len() > 1 {
        let mut dist_min = Precision::INFINITE;
        let mut best: Option<(i32, Arc<dyn Curve>)> = None;
        for (idx, c) in &candidates {
            let pnt = c.d0(c.first_parameter());
            let dist = pnt.distance(pfv);
            if dist - dist_min < -CONFUSION {
                dist_min = dist;
                best = Some((*idx, c.clone()));
            }
        }
        best?
    } else {
        candidates.into_iter().next()?
    };

    // `cxx:1232-1240`
    let mut curves: Vec<(i32, Arc<dyn Curve>)> =
        lines.drain(..).filter(|(idx, _)| *idx != first_idx).collect();

    // `cxx:1242-1273`: walk the pieces end to end from the first one's end.
    let mut to_glue: Vec<Arc<dyn Curve>> = Vec::new();
    let mut end_point = first_curve.d0(first_curve.last_parameter());
    let mut added = true;
    while added {
        added = false;
        let mut i = 0;
        while i < curves.len() {
            let curve = curves[i].1.clone();
            let pfirst = curve.d0(curve.first_parameter());
            let plast = curve.d0(curve.last_parameter());
            if pfirst.distance(&end_point) <= CONFUSION {
                to_glue.push(curve);
                end_point = plast;
                curves.remove(i);
                added = true;
                break;
            }
            if plast.distance(&end_point) <= CONFUSION {
                to_glue.push(Arc::from(curve.reversed()));
                end_point = pfirst;
                curves.remove(i);
                added = true;
                break;
            }
            i += 1;
        }
    }

    // `cxx:1280-1290`: `GeomConvert_CompCurveToBSplineCurve Concat(FirstCurve)`
    // then `Add(ToGlue(i), Precision::Confusion(), true)`.
    let mut concat =
        CompCurveToBSplineCurve::from_curve(first_curve.as_ref(), ParameterisationType::TgtThetaOver2)
            .ok()?;
    for g in &to_glue {
        // OCCT keeps `myCurve` untouched when `Add` answers `false`; only a
        // conversion failure throws.
        concat
            .add(g.as_ref(), CONFUSION, true, true, 0)
            .expect("GeomConvert_CompCurveToBSplineCurve::Add");
    }
    let new_c: Arc<dyn Curve> = Arc::new(concat.bspline_curve()?.clone());

    // `cxx:1292-1308`: `Extrema_ExtPC myExtPC(pfv, TheCurve)` over the glued
    // curve, smallest `IsMin` solution.
    let (sols, is_done) =
        extrema_ext_pc_range_done(new_c.as_ref(), pfv, new_c.first_parameter(), new_c.last_parameter());
    if !is_done {
        panic!("Standard_Failure: Draft_Modification_1::Perform: ExtremaPC not done.");
    }
    let mut pmin = 0.0f64;
    let mut dist2_min = f64::MAX;
    for sol in &sols {
        if !sol.is_min {
            continue;
        }
        if sol.sq_dist < dist2_min {
            dist2_min = sol.sq_dist;
            pmin = sol.u;
        }
    }

    // `cxx:1308-1313`
    let (_, new_d1) = new_c.d1(pmin);
    if d1fv.dot(&new_d1) < 0.0 {
        return Some(Arc::from(new_c.reversed()));
    }
    Some(new_c)
}

/// The `GeomAPI_ProjectPointOnSurf` + `LowerDistanceParameters` seam test of
/// `_1.cxx:1169-1186`.
fn end_on_seam(s1: &dyn Surface, pnt: &GpPnt) -> bool {
    // The reference never tests `IsDone()`; a missing solution is undefined
    // there, so a failure keeps both comparisons false here.
    let Some(on) = project_point_on_surface(s1, pnt, CONFUSION) else {
        return false;
    };
    let u = on.u;
    u.abs() <= CONFUSION || (u - 2.0 * std::f64::consts::PI).abs() <= CONFUSION
}
