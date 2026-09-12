//! Remainder of `BOPAlgo_PaveFiller_3.cxx`: `ForceInterfVE`, `GetPBBox`,
//! `UpdateVerticesOfCB`.
//!
//! Source: `BOPAlgo_PaveFiller_3.cxx` (`ForceInterfVE` at 828, `GetPBBox` at
//! 914, `UpdateVerticesOfCB` at 959).

use std::collections::HashMap;

use occt_core::bnd::BndBox;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::bbox_from_geometry::shape_bbox;
use crate::bopds::{BopdsPave, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::int_tools_curve_box::add_curve_to_box;
use crate::pave_common::update_vertex_sd;
use crate::pave_ef::has_interf_shape_sub_shapes;
use crate::pave_ff_exist::pb_key;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Vertex};

/// Bounding box of `edge` over `[first, last]`, enlarged by edge tolerance
/// plus `Precision::Confusion`. Port of `BndLib_Add3dCurve::Add` in
/// `GetPBBox` (`_3.cxx:948-951`).
fn pb_range_box(edge: &Edge, first: f64, last: f64) -> BndBox {
    let a_tol = BRepTool::edge_tolerance(edge) + CONFUSION;
    if let Some(c) = BRepTool::edge_curve_world(edge) {
        let mut box_ = BndBox::new();
        add_curve_to_box(c.as_ref(), first, last, a_tol, &mut box_);
        box_
    } else {
        let mut box_ = shape_bbox(&edge.0);
        box_.enlarge(a_tol);
        box_
    }
}

/// `BOPAlgo_PaveFiller::GetPBBox`.
pub fn get_pb_box(
    the_e: &Edge,
    the_pb: &BopdsPaveBlock,
    pb_box: &mut HashMap<(usize, u64, u64), BndBox>,
) -> Option<(f64, f64, f64, f64, BndBox)> {
    let (first, last) = the_pb.range();
    if last - first <= PCONFUSION {
        return None;
    }
    if the_pb.has_shrunk_data() {
        let (sf, sl, _) = the_pb.shrunk_data();
        return Some((first, last, sf, sl, pb_range_box(the_e, sf, sl)));
    }
    let key = pb_key(the_pb);
    if let Some(b) = pb_box.get(&key) {
        return Some((first, last, first, last, b.clone()));
    }
    let box_ = pb_range_box(the_e, first, last);
    pb_box.insert(key, box_.clone());
    Some((first, last, first, last, box_))
}

/// `BOPAlgo_PaveFiller::ForceInterfVE`.
pub fn force_interf_ve(
    f: &mut PaveFiller,
    n_v: usize,
    a_pb: &mut BopdsPaveBlock,
    the_m_edges: &mut Vec<usize>,
) -> Result<bool, String> {
    let n_e = a_pb.original_edge();
    let Some(si_e) = f.ds().shape_info(n_e).cloned() else {
        return Ok(false);
    };
    if si_e.has_subshape(n_v) {
        return Ok(true);
    }
    if f.ds().has_interf_pair(n_v, n_e) {
        return Ok(true);
    }
    if has_interf_shape_sub_shapes(f.ds(), n_v, n_e) {
        return Ok(true);
    }
    let (pv1, pv2) = a_pb.indices();
    if pv1 == n_v || pv2 == n_v {
        return Ok(true);
    }
    let n_vx0 = f.ds().get_same_domain_index(n_v);
    let Some(v_shape) = f.ds().shape(n_vx0).cloned() else {
        return Ok(false);
    };
    let Some(e_shape) = f.ds().shape(n_e).cloned() else {
        return Ok(false);
    };
    let v = Vertex(v_shape);
    let e = Edge(e_shape);
    let fuzzy = f.fuzzy_value();
    let (flag, t, dist) = f.context().compute_pe_pnt(
        &BRepTool::vertex_point(&v),
        BRepTool::vertex_tolerance(&v) + fuzzy,
        &e,
    );
    if flag != 0 && flag != -4 {
        return Ok(false);
    }
    if flag == -4 {
        return Ok(false);
    }
    f.ds_mut().add_interf_ve(n_v, n_e, None);
    let tol_new = BRepTool::vertex_tolerance(&v).max(dist);
    let n_vx = update_vertex_sd(f, n_v, tol_new)?;
    if f.ds().is_new_shape(n_vx) {
        for it in f.ds_mut().interf_ve_mut() {
            if it.contains(n_v) && it.contains(n_e) {
                it.set_index_new(n_vx);
                it.set_common_range(t, t);
                break;
            }
        }
    }
    a_pb.append_ext_pave(BopdsPave::new(n_vx, t));
    if !the_m_edges.contains(&n_e) {
        the_m_edges.push(n_e);
    }
    let i_rv = f.ds().rank(n_v);
    let i_re = f.ds().rank(n_e);
    if i_rv == i_re {
        f.add_warning(format!(
            "acquired self-interference: vertex {n_v} and edge {n_e} of one argument"
        ));
    }
    Ok(true)
}

/// Write a pave block that was cloned out of the DS back to its pool slot.
/// OCCT mutates the handle in place; the Rust clone must be stored.
pub fn write_pave_block(f: &mut PaveFiller, pb: &BopdsPaveBlock) {
    let orig = pb.original_edge();
    if orig >= f.ds().nb_shapes() || !f.ds().has_pave_blocks(orig) {
        return;
    }
    let key = pb_key(pb);
    let blocks = f.ds_mut().change_pave_blocks_mut(orig);
    for slot in blocks.iter_mut() {
        if pb_key(slot) == key {
            *slot = pb.clone();
            return;
        }
    }
}

/// `BOPAlgo_PaveFiller::UpdateVerticesOfCB`.
pub fn update_vertices_of_cb(f: &mut PaveFiller) -> Result<(), String> {
    let mut fence: Vec<(usize, u64, u64)> = Vec::new();
    let n = f.ds().nb_shapes();
    let mut jobs: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        if !f.ds().has_pave_blocks(i) {
            continue;
        }
        let blocks = f.ds().pave_blocks(i).to_vec();
        for pb in blocks {
            let Some(cb) = f.ds().common_block(&pb).cloned() else {
                continue;
            };
            let Some(pbr) = cb.pave_block1() else {
                continue;
            };
            let k = pb_key(pbr);
            if fence.contains(&k) {
                continue;
            }
            fence.push(k);
            let tol = cb.tolerance();
            if tol > 0.0 {
                let (n1, n2) = pbr.indices();
                jobs.push((n1, n2, tol));
            }
        }
    }
    for (n1, n2, tol) in jobs {
        update_vertex_sd(f, n1, tol)?;
        update_vertex_sd(f, n2, tol)?;
    }
    Ok(())
}
