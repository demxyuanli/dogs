//! `IntPatch_PrmPrmIntersection::Perform(Surf1, D1, Surf2, D2, U1Depart,
//! V1Depart, U2Depart, V2Depart, TolTangency, Epsilon, Deflection,
//! Increment)` (`IntPatch_PrmPrmIntersection.cxx:2176-2325`).
//!
//! A single user-supplied start point: walk from it, keep the walking line
//! when it has more than two points and survives `PutVertexOnLine`, append it
//! directly (`cxx:2320`). There is no duplicate check against the lines
//! already in `SLin` and no `AddWLine` call.
//!
//! When `PutVertexOnLine` cuts the line below two points, the whole `Perform`
//! returns (`cxx:2291`, `:2298`) - the empty `returns`, not `continue`.

use occt_core::gp::GpVec;
use occt_geom::Surface;

use crate::brep_top_adaptor_topol_tool::TopolTool as BrepTopolTool;
use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::{PatchPoint, WLine, WLineWay};
use crate::intpatch_line::IntPatchLine;
use crate::intpatch_rst_int::put_vertex_on_line;
use crate::intsurf_transition::IntSurfTypeTrans;
use crate::intwalk::PWalking;

/// `Perform(...)` driven by one start point (`cxx:2176-2325`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_start_point(
    slin: &mut Vec<GeomIntLine>,
    s1: &dyn Surface,
    d1: &mut BrepTopolTool,
    s2: &dyn Surface,
    d2: &mut BrepTopolTool,
    start: [f64; 4],
    tol_tangency: f64,
    epsilon: f64,
    deflection: f64,
    increment: f64,
) {
    let mut pw = PWalking::new(s1, s2, tol_tangency, epsilon, deflection, increment);

    if pw.perform_first_point(start).is_none() {
        return;
    }
    pw.perform(start);
    if !pw.is_done() || pw.nb_points() <= 2 {
        return;
    }

    let p_debut = pw.value(1).value();
    let p_fin = pw.value(pw.nb_points()).value();

    // `cxx:2261-2282`: transitions from the normals at the tangent point.
    let (tg_dir, indextg) = pw.tangent_at_line();
    let tgline = GpVec::from_xyz(tg_dir.xyz());
    let (locu, locv) = pw.value(indextg).parameters_on_s1();
    let (_, d1u, d1v) = s1.d1(locu, locv);
    let norm1 = d1u.crossed(&d1v);
    let (locu, locv) = pw.value(indextg).parameters_on_s2();
    let (_, d2u, d2v) = s2.d1(locu, locv);
    let norm2 = d2u.crossed(&d2v);
    let (trans1, trans2) = if tgline.dot_cross(&norm2, &norm1) > 0.0 {
        (IntSurfTypeTrans::Out, IntSurfTypeTrans::In)
    } else {
        (IntSurfTypeTrans::In, IntSurfTypeTrans::Out)
    };

    let mut wl = WLine::new();
    wl.base = IntPatchLine::with_transitions(false, trans1, trans2);
    wl.set_creating_way(WLineWay::PrmPrm);
    for p in pw.line() {
        wl.add(*p);
    }

    // `PutVertexOnLine(wline, Surf1, D1, Surf2, true, TolTang)` (`cxx:2289`).
    let mut line = GeomIntLine::Walking(wl);
    if put_vertex_on_line(&mut line, s1, d1, s2, true, tol_tangency).is_err() {
        return;
    }
    if line.as_wline().map_or(true, |w| w.nb_pnts() < 2) {
        return;
    }
    // `PutVertexOnLine(wline, Surf2, D2, Surf1, false, TolTang)` (`cxx:2293`).
    if put_vertex_on_line(&mut line, s2, d2, s1, false, tol_tangency).is_err() {
        return;
    }
    let mut wl = match line {
        GeomIntLine::Walking(w) => w,
        _ => return,
    };
    if wl.nb_pnts() < 2 {
        return;
    }

    // `cxx:2301-2317`: fallback vertices when none were attached.
    if wl.nb_vertex() == 0 {
        let params = pw.value(1).parameters();
        let mut vtx = PatchPoint::new(p_debut, 1.0, params.0, params.1, params.2, params.3);
        vtx.set_tolerance(tol_tangency);
        wl.vertices.push(vtx);

        let n = wl.nb_pnts();
        let params = pw.value(n).parameters();
        let mut vtx = PatchPoint::new(p_fin, n as f64, params.0, params.1, params.2, params.3);
        vtx.set_tolerance(tol_tangency);
        wl.vertices.push(vtx);
    }

    // `SLin.Append(wline); empt = false;` (`cxx:2320-2321`).
    slin.push(GeomIntLine::Walking(wl));
}
