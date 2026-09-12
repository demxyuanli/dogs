//! Remainder of `BOPAlgo_PaveFiller_5.cxx`: `CheckFacePaves`, `ForceInterfVF`,
//! `ReduceIntersectionRange`, and the `myDistances` write of `PerformEF`.
//!
//! Source: `BOPAlgo_PaveFiller_5.cxx` (`CheckFacePaves` at 596 and 605,
//! `ForceInterfVF` at 631, `ReduceIntersectionRange` at 685, `PerformEF`
//! distance branch at 349).

use std::collections::HashMap;

use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::bopds::{BopdsDS, BopdsPaveBlock};
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::int_tools_full::IntToolsContext;
use crate::pave_common::update_vertex_sd;
use crate::pave_filler::{EdgeRangeDistance, PaveFiller};
use crate::shape::{Edge, Face, Vertex};

/// `BOPAlgo_PaveFiller::CheckFacePaves` (vertex index vs On/In maps).
pub fn check_face_paves_index(n_vx: usize, on: &[usize], inn: &[usize]) -> bool {
    on.contains(&n_vx) || inn.contains(&n_vx)
}

/// `BOPAlgo_PaveFiller::CheckFacePaves` (new vertex vs On map).
pub fn check_face_paves_vertex(ds: &BopdsDS, v_new: &Vertex, on: &[usize]) -> bool {
    let p = BRepTool::vertex_point(v_new);
    let tol = BRepTool::vertex_tolerance(v_new);
    for &n_v in on {
        let Some(s) = ds.shape(n_v) else {
            continue;
        };
        if AlgoTools::compute_vv(s, &p, tol) == 1 {
            return true;
        }
    }
    false
}

/// `IntTools_Context::ComputeVF` returning `(flag, u, v, tol_v_new)`.
pub fn compute_vf_uv(
    ctx: &mut IntToolsContext,
    vertex: &Vertex,
    face: &Face,
    fuzzy: f64,
) -> (i32, f64, f64, f64) {
    let p = BRepTool::vertex_point(vertex);
    let Some(surf) = BRepTool::face_surface(face) else {
        return (-1, 0.0, 0.0, 0.0);
    };
    let Ok((u, v)) = ctx.project_point_on_face(face, &p) else {
        return (-1, 0.0, 0.0, 0.0);
    };
    let dist = surf.d0(u, v).distance(&p);
    let tol_sum = BRepTool::vertex_tolerance(vertex)
        + BRepTool::face_tolerance(face)
        + fuzzy.max(CONFUSION);
    if dist > tol_sum {
        return (-2, u, v, dist);
    }
    let state = match ctx.state_point_face(face, (u, v), fuzzy) {
        Ok(s) => s,
        Err(_) => return (-1, u, v, dist),
    };
    if state == FaceState::Out {
        return (-3, u, v, dist);
    }
    let tol_new = BRepTool::vertex_tolerance(vertex).max(dist);
    (0, u, v, tol_new)
}

/// `BOPAlgo_PaveFiller::ForceInterfVF`.
pub fn force_interf_vf(f: &mut PaveFiller, n_v: usize, n_f: usize) -> Result<bool, String> {
    let Some(v_shape) = f.ds().shape(n_v).cloned() else {
        return Ok(false);
    };
    let Some(f_shape) = f.ds().shape(n_f).cloned() else {
        return Ok(false);
    };
    let v = Vertex(v_shape);
    let face = Face(f_shape);
    let fuzzy = f.fuzzy_value();
    let (flag, u, vv, tol_new) = {
        let ctx = f.context_mut();
        compute_vf_uv(ctx, &v, &face, fuzzy)
    };
    if flag != 0 && flag != -2 {
        return Ok(false);
    }
    f.ds_mut().add_interf_vf(n_v, n_f, None);
    let n_vx = update_vertex_sd(f, n_v, tol_new)?;
    if f.ds().is_new_shape(n_vx) {
        let _ = f.ds_mut().bind_vf_new_vertex(n_v, n_f, n_vx);
    }
    if let Some(fi) = f.ds_mut().face_info_mut(n_f) {
        fi.add_vert_in(n_vx);
        fi.add_vert(n_vx, u, vv);
    }
    let i_rv = f.ds().rank(n_v);
    let i_rf = f.ds().rank(n_f);
    if i_rv == i_rf {
        f.add_warning(format!(
            "acquired self-interference: vertex {n_v} and face {n_f} of one argument"
        ));
    }
    Ok(true)
}

/// `BOPDS_DS::HasInterfShapeSubShapes`.
pub fn has_interf_shape_sub_shapes(ds: &BopdsDS, i1: usize, i2: usize) -> bool {
    let Some(si) = ds.shape_info(i2) else {
        return false;
    };
    for &sub in si.sub_shapes() {
        if ds.has_interf_pair(i1, sub) {
            return true;
        }
    }
    false
}

/// `BOPAlgo_PaveFiller::ReduceIntersectionRange`.
pub fn reduce_intersection_range(
    ds: &BopdsDS,
    the_v1: usize,
    the_v2: usize,
    the_e: usize,
    the_f: usize,
    ts1: &mut f64,
    ts2: &mut f64,
) {
    if !ds.is_new_shape(the_v1) && !ds.is_new_shape(the_v2) {
        return;
    }
    if !has_interf_shape_sub_shapes(ds, the_e, the_f) {
        return;
    }
    if ds.interf_ee().is_empty() {
        return;
    }
    let mut mfe: Vec<usize> = Vec::new();
    if let Some(si) = ds.shape_info(the_f) {
        for &n in si.sub_shapes() {
            if ds.shape_info(n).map(|s| s.shape_type()) == Some(ShapeType::Edge) {
                mfe.push(n);
            }
        }
    }
    for ee in ds.interf_ee() {
        let Some(n_v) = ee.get_index_new() else {
            continue;
        };
        if n_v != the_v1 && n_v != the_v2 {
            continue;
        }
        let (n_e1, n_e2) = ee.indices();
        if (the_e != n_e1 && the_e != n_e2) || (!mfe.contains(&n_e1) && !mfe.contains(&n_e2)) {
            continue;
        }
        let (tr1, tr2) = ee.common_range();
        if n_v == the_v1 {
            if *ts1 < tr2 {
                *ts1 = tr2;
            }
        } else if *ts2 > tr1 {
            *ts2 = tr1;
        }
    }
}

/// Write `myDistances` for an EF pair with no common part (`_5.cxx:349`).
pub fn record_ef_distance(
    distances: &mut HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    n_e: usize,
    n_f: usize,
    pb: &BopdsPaveBlock,
    dist: f64,
    tol_e: f64,
    tol_f: f64,
) {
    if dist >= f64::MAX || dist <= tol_e + tol_f {
        return;
    }
    let (t1, t2) = pb.range();
    distances
        .entry((n_e, n_f))
        .or_default()
        .push(EdgeRangeDistance::new(t1, t2, dist));
}

/// `IntTools_Tools::IsInRange`.
pub fn is_in_range(r1_first: f64, r1_last: f64, r2_first: f64, r2_last: f64, tol: f64) -> bool {
    let a = (r1_first - tol).min(r1_last - tol);
    let b = (r1_first + tol).max(r1_last + tol);
    let c = r2_first.min(r2_last);
    let d = r2_first.max(r2_last);
    c <= b && d >= a
}

/// Collect On/In vertex indices of a face from `FaceInfo`.
pub fn face_paves_on_in(ds: &BopdsDS, n_f: usize) -> (Vec<usize>, Vec<usize>) {
    let Some(fi) = ds.face_info(n_f) else {
        return (Vec::new(), Vec::new());
    };
    (fi.verts_on().to_vec(), fi.verts_in().to_vec())
}
