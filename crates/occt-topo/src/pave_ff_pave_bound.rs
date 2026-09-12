//! getBoundPaves, PutBoundPaveOnCurve, PutClosingPaveOnCurve.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`getBoundPaves` at 2255,
//! `PutBoundPaveOnCurve` at 2308, `PutClosingPaveOnCurve` at 3500).
//!
//! Closing paves require a non-empty `BRepLib::FindValidRange` after the
//! vertex grows to `max(tolV, dist + DTolerance)`. Bound paves skip a closed
//! curve that already has a vertex on either end.

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::algo_tools::{AlgoTools, D_TOLERANCE};
use crate::bopds::{BopdsDS, BopdsPave};
use crate::bopds_ff::BopdsCurve;
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::inttools_range;
use crate::pave_filler::PaveFiller;
use crate::pave_update_sd;
use crate::shape::{Face, Vertex};

/// `getBoundPaves` (`_6.cxx:2255`).
///
/// Extreme ext-paves compared with the curve ends through `ComputeVV`.
/// A miss (`compute_vv == 0` in this port, OCCT `iFlag != 0`) clears that end.
pub fn get_bound_paves(ds: &BopdsDS, nc: &BopdsCurve) -> [Option<usize>; 2] {
    let Some(pb) = nc.pave_blocks().first() else {
        return [None, None];
    };
    if pb.ext_paves().is_empty() {
        return [None, None];
    }
    let mut tmin = f64::MAX;
    let mut tmax = f64::MIN;
    let mut nv = [None, None];
    for p in pb.ext_paves() {
        if p.param < tmin {
            nv[0] = Some(p.index);
            tmin = p.param;
        }
        if p.param > tmax {
            nv[1] = Some(p.index);
            tmax = p.param;
        }
    }
    let (t0, t1) = nc.range();
    let Some(curve) = nc.curve() else {
        return nv;
    };
    let pts = [curve.d0(t0), curve.d0(t1)];
    let tol = nc.tolerance().max(nc.tangential_tolerance()) + CONFUSION;
    for j in 0..2 {
        let Some(n) = nv[j] else {
            continue;
        };
        let Some(v) = ds.shape(n) else {
            nv[j] = None;
            continue;
        };
        if AlgoTools::compute_vv(v, &pts[j], tol) == 0 {
            nv[j] = None;
        }
    }
    nv
}

/// `BOPAlgo_PaveFiller::PutBoundPaveOnCurve` (`_6.cxx:2308`).
pub fn put_bound_pave_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    face1: &Face,
    face2: &Face,
    ctx: &mut IntToolsContext,
) -> Result<Vec<usize>, String> {
    let (curve, range, tol_r3d) = {
        let nc = &f.ds().interf_ff()[i].curves()[j];
        let Some(c) = nc.curve().cloned() else {
            return Ok(Vec::new());
        };
        (
            c,
            nc.range(),
            nc.tolerance().max(nc.tangential_tolerance()),
        )
    };
    let (t0, t1) = range;
    let p0 = curve.d0(t0);
    let p1 = curve.d0(t1);
    let bnd = get_bound_paves(f.ds(), &f.ds().interf_ff()[i].curves()[j]);
    let is_closed = p1.distance(&p0) <= CONFUSION;
    if is_closed && (bnd[0].is_some() || bnd[1].is_some()) {
        return Ok(Vec::new());
    }
    let mut lbv = Vec::new();
    let ts = [t0, t1];
    let ps = [p0, p1];
    for j_end in 0..2 {
        if bnd[j_end].is_some() {
            continue;
        }
        if j_end == 1 && is_closed {
            continue;
        }
        let ok_f = ctx
            .is_point_in_on_face(face1, &ps[j_end], None, tol_r3d)
            .unwrap_or(false);
        let ok_g = ctx
            .is_point_in_on_face(face2, &ps[j_end], None, tol_r3d)
            .unwrap_or(false);
        if !ok_f || !ok_g {
            continue;
        }
        let mut vn = AlgoTools::make_new_vertex(&ps[j_end], tol_r3d)?;
        let dist = BRepTool::vertex_point(&Vertex(vn.clone())).distance(&ps[j_end]);
        let need = dist + D_TOLERANCE;
        if BRepTool::vertex_tolerance(&Vertex(vn.clone())) < need {
            Vertex(vn.clone()).set_tolerance(need);
        }
        let n_vn = f.ds_mut().append(vn)?;
        f.ds_mut().refresh_vertex_box(n_vn, need);
        let pb = f.ds_mut().interf_ff_mut()[i].change_curves()[j].change_pave_block1();
        pb.append_ext_pave(BopdsPave::new(n_vn, ts[j_end]));
        lbv.push(n_vn);
    }
    Ok(lbv)
}

/// `BOPAlgo_PaveFiller::PutClosingPaveOnCurve` (`_6.cxx:3500`).
pub fn put_closing_pave_on_curve(f: &mut PaveFiller, i: usize, j: usize) {
    let (curve, range, tol_r3d) = {
        let nc = &f.ds().interf_ff()[i].curves()[j];
        let Some(c) = nc.curve().cloned() else {
            return;
        };
        (c, nc.range(), nc.tolerance().max(nc.tangential_tolerance()))
    };
    let (t0, t1) = range;
    if !t0.is_finite() || !t1.is_finite() {
        return;
    }
    let p0 = curve.d0(t0);
    let p1 = curve.d0(t1);
    let ts = [t0, t1];
    let ps = [p0, p1];
    let ext: Vec<BopdsPave> = f.ds().interf_ff()[i].curves()[j]
        .pave_blocks()
        .first()
        .map(|pb| pb.ext_paves().to_vec())
        .unwrap_or_default();
    let mut found: Option<(usize, f64, GpPnt)> = None;
    for pave in &ext {
        for k in 0..2 {
            if (pave.param - ts[k]).abs() < PCONFUSION {
                found = Some((pave.index, ts[1 - k], ps[1 - k]));
                break;
            }
        }
        if found.is_some() {
            break;
        }
    }
    let Some((mut n_v, t_op, p_op)) = found else {
        return;
    };
    let Some(v_shape) = f.ds().shape(n_v).cloned() else {
        return;
    };
    let vtx = Vertex(v_shape);
    let a_tol_v = BRepTool::vertex_tolerance(&vtx);
    let pv = BRepTool::vertex_point(&vtx);
    let a_tol_p = tol_r3d + CONFUSION;
    let dist = pv.distance(&p_op);
    if dist > a_tol_v + a_tol_p {
        return;
    }
    let a_new_tol_v = a_tol_v.max(dist + D_TOLERANCE);
    if inttools_range::find_valid_range(
        curve.as_ref(),
        t0,
        t1,
        &ps[0],
        a_new_tol_v,
        &ps[1],
        a_new_tol_v,
    )
    .is_none()
    {
        return;
    }
    if a_new_tol_v > a_tol_v {
        if let Ok(n_vn) = pave_update_sd::update_vertex(f, n_v, a_new_tol_v) {
            n_v = n_vn;
        }
    }
    let pb = f.ds_mut().interf_ff_mut()[i].change_curves()[j].change_pave_block1();
    pb.append_ext_pave1(BopdsPave::new(n_v, t_op));
}
