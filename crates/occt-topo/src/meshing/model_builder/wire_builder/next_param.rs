use super::pcurve::{pcurve_and_range, sample_pcurve};
use super::wire_edges::wire_edges_sewd;

use super::super::prelude::*;

/// `GetNextParamOnPC` (`BRepTools_WireExplorer.cxx:798-873`): the parameter of
/// the first sample of `pc` more than `(tol_u, tol_v)` away from `p_ref`,
/// stepping from the `reverse` end of `[f, l]` in `|l - f| / 1000` increments.
pub(super) fn next_param_on_pc(
    pc: &dyn Curve2d,
    p_ref: &GpPnt2d,
    f: f64,
    l: f64,
    tol_u: f64,
    tol_v: f64,
    reverse: bool,
) -> f64 {
    let d_p = (l - f).abs() / 1000.0;
    let on_ref = |p: &GpPnt2d| {
        (p_ref.x() - p.x()).abs() < tol_u && (p_ref.y() - p.y()).abs() < tol_v
    };
    if reverse {
        let mut result = f;
        let mut start_par = f;
        let mut found = false;
        while !found && start_par < l {
            start_par += d_p;
            if on_ref(&pc.d0(start_par)) {
                continue;
            }
            result = start_par;
            found = true;
        }
        if !found {
            result = l;
        }
        if result > l {
            result = l;
        }
        result
    } else {
        let mut result = l;
        let mut start_par = l;
        let mut found = false;
        while !found && start_par > f {
            start_par -= d_p;
            if on_ref(&pc.d0(start_par)) {
                continue;
            }
            result = start_par;
            found = true;
        }
        if !found {
            result = f;
        }
        if result < f {
            result = f;
        }
        result
    }
}

/// `ShapeAnalysis::TotCross2D` — signed 2D area of a wire's pcurves (trapezoid
/// rule over sampled pcurve points, sequence reversed per REVERSED edge).
pub(super) fn wire_area_2d(wire: &Wire, face: &Face) -> f64 {
    let mut totcross = 0.0;
    let mut uv0: Option<GpPnt2d> = None;
    let mut fuv = GpPnt2d::new(0.0, 0.0);
    let mut nbc = 0usize;
    for e in wire_edges_sewd(wire) {
        // `BRep_Tool::CurveOnSurface(edge, face)` picks PCurve2 on a reversed
        // seam (`BRep_Tool.cxx:354-357`) and returns the COS range; then
        // TotCross2D reverses the samples for a REVERSED edge.
        let Some((pc, a, b)) = pcurve_and_range(&e, face) else {
            continue;
        };
        let mut pts = sample_pcurve(pc.as_ref(), a, b);
        if e.0.orientation().is_reversed() {
            pts.reverse();
        }
        nbc += 1;
        if nbc == 1 {
            fuv = pts[0];
            uv0 = Some(pts[0]);
        }
        for p in &pts {
            totcross += (fuv.x() - p.x()) * (fuv.y() + p.y()) / 2.0;
            fuv = *p;
        }
    }
    if let Some(u0) = uv0 {
        totcross += (fuv.x() - u0.x()) * (fuv.y() + u0.y()) / 2.0;
    }
    totcross
}
