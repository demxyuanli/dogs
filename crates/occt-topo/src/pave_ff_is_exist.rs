//! `IsExistingVertex` and both `IsExistingPaveBlock` overloads.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx`:
//! - `IsExistingVertex` at 1950
//! - shared-edge `IsExistingPaveBlock` at 1988
//! - ON/IN-tree `IsExistingPaveBlock` at 2047 (the linear AABB stand-in is
//!   already in `pave_ff_exist`; this module documents the OCCT gates that
//!   module implements and exposes the shared-edge overload used by
//!   `MakeBlocks`).
//!
//! Shared-edge overload: mid-point of the section block, box vs each shared
//! edge, `ComputePE(Pm, max(tolE, tolV) + Fuzzy, E)`. A hit returns the edge
//! index and the projection distance as the new tolerance.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_geom::Curve;

use crate::algo_tools::AlgoTools;
use crate::bopds::{BopdsDS, BopdsPaveBlock};
use crate::boptools_2d::intermediate_point;
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Vertex};

/// `BOPAlgo_PaveFiller::IsExistingVertex` (`_6.cxx:1950`).
///
/// OCCT `ComputeVV` returns 0 on coincidence; this port's `compute_vv`
/// returns 1, so a hit is `!= 0`.
pub fn is_existing_vertex(
    ds: &BopdsDS,
    p: &GpPnt,
    tol_r3d: f64,
    mv_on_in: &std::collections::HashSet<usize>,
    fuzzy: f64,
) -> bool {
    let tol_check = tol_r3d + fuzzy;
    let mut box_p = BndBox::new();
    box_p.add_point(p);
    box_p.enlarge(tol_r3d);
    for &n_v in mv_on_in {
        if let Some(bv) = ds.box_of(n_v) {
            if box_p.is_out_box(bv) {
                continue;
            }
        }
        let Some(v) = ds.shape(n_v) else {
            continue;
        };
        if AlgoTools::compute_vv(v, p, tol_check) != 0 {
            return true;
        }
    }
    false
}

/// Shared-edge `IsExistingPaveBlock` (`_6.cxx:1988`).
pub fn is_existing_pave_block_on_shared(
    ds: &BopdsDS,
    ctx: &IntToolsContext,
    pb: &BopdsPaveBlock,
    curve: &dyn Curve,
    lse: &[usize],
    fuzzy: f64,
) -> Option<(usize, f64)> {
    if lse.is_empty() {
        return None;
    }
    let (t1, t2) = pb.range();
    let (n_v1, n_v2) = pb.indices();
    let mut tol = 0.0f64;
    if let Some(v) = ds.shape(n_v1) {
        tol = tol.max(BRepTool::vertex_tolerance(&Vertex(v.clone())));
    }
    if let Some(v) = ds.shape(n_v2) {
        tol = tol.max(BRepTool::vertex_tolerance(&Vertex(v.clone())));
    }
    let tm = intermediate_point(t1, t2);
    let pm = curve.d0(tm);
    let mut box_pm = BndBox::new();
    box_pm.add_point(&pm);
    box_pm.enlarge(tol);
    for &n_e in lse {
        if n_e == 0 {
            continue;
        }
        if let Some(be) = ds.box_of(n_e) {
            if be.is_out_box(&box_pm) {
                continue;
            }
        }
        let Some(e_shape) = ds.shape(n_e) else {
            continue;
        };
        let e = Edge(e_shape.clone());
        let tol_e = BRepTool::edge_tolerance(&e);
        let tol_check = tol_e.max(tol) + fuzzy;
        let (i_flag, _tx, dist) = ctx.compute_pe_pnt(&pm, tol_check, &e);
        if i_flag == 0 {
            return Some((n_e, dist));
        }
    }
    None
}

/// Limit values used by the ON/IN `IsExistingPaveBlock` when both ends
/// classify as `iflag == 2` and the candidate has no common block with a
/// face (`_6.cxx:2108-2112`).
pub fn thin_face_max_tol_add(tol_check: f64) -> f64 {
    const A_MAX_TOL_ADD: f64 = 0.001;
    const A_COEFF_TOL_ADD: f64 = 10.0;
    A_MAX_TOL_ADD.min(A_COEFF_TOL_ADD * tol_check)
}

/// Whether two 1D ranges overlap (used by ProcessExisting distance lookup).
pub fn range_contains_param(t1: f64, t2: f64, f: f64, l: f64) -> bool {
    (t1 <= f && f <= t2) || (t1 <= l && l <= t2) || (f <= t1 && t1 <= l) || (f <= t2 && t2 <= l)
}

/// `ComputePE` hit: OCCT `iFlag == 0`.
pub fn compute_pe_hit(ctx: &IntToolsContext, p: &GpPnt, tol_check: f64, edge: &Edge) -> Option<(f64, f64)> {
    let (st, t, dist) = ctx.compute_pe_pnt(p, tol_check, edge);
    if st == 0 {
        Some((t, dist))
    } else {
        None
    }
}

/// Vertex-edge `ComputeVE` as used by the bound ProcessExisting overload.
/// OCCT returns 0 on a hit; this port's `compute_pe` matches that code.
pub fn compute_ve_hit(ctx: &IntToolsContext, vertex: &Vertex, edge: &Edge, fuzzy: f64) -> bool {
    ctx.compute_pe(vertex, edge, fuzzy) == 0
}
