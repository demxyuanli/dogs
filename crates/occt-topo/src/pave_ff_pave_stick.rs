//! PutStickPavesOnCurve, PutEFPavesOnCurve, EstimatePaveOnCurve.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx` (`PutEFPavesOnCurve` at 2692,
//! `PutStickPavesOnCurve` at 2748, `EstimatePaveOnCurve` at 4056).
//!
//! Stick paves: unused stick vertices, at least one free curve end, both
//! pcurves present, vertex within `aDT2 = 2e-7` of a free end, crease
//! `1 - |n1·n2| <= 5e-9`. EF paves: Bezier/BSpline only, `RemoveUsedVertices`,
//! `ProjPT` then `PutPaveOnCurve` with the projection distance.

use std::collections::{HashMap, HashSet};

use crate::algo_tools3d;
use crate::bopds::BopdsDS;
use crate::bopds_ff::BopdsCurve;
use crate::brep_tool::BRepTool;
use crate::int_tools_vertex_line;
use crate::pave_ff_misc::remove_used_vertices;
use crate::pave_ff_pave_put::put_pave_on_curve;
use crate::pave_filler::PaveFiller;
use crate::shape::{Face, Vertex};

/// Stick-pave "rich" criteria — square distance to a bound (`_6.cxx:2793`).
const STICK_DT2: f64 = 2.0e-7;

/// Stick-pave crease criteria — `1 - |n1·n2|` (`_6.cxx:2794`).
const STICK_DSCPR: f64 = 5.0e-9;

/// `BOPAlgo_PaveFiller::PutStickPavesOnCurve` (`_6.cxx:2748`).
pub fn put_stick_paves_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    face1: &Face,
    face2: &Face,
    mi: &HashSet<usize>,
    mv_stick: &HashSet<usize>,
    mv_tol: &mut HashMap<usize, f64>,
    dmvlv: &mut HashMap<usize, Vec<usize>>,
) {
    let bnd = crate::pave_ff_pave_bound::get_bound_paves(f.ds(), &f.ds().interf_ff()[i].curves()[j]);
    if bnd[0].is_some() && bnd[1].is_some() {
        return;
    }
    let mut leftover: HashSet<usize> = mv_stick.clone();
    remove_used_vertices(f.ds().interf_ff()[i].curves(), &mut leftover);
    if leftover.is_empty() {
        return;
    }
    let (curve, range, pc1, pc2) = {
        let nc = &f.ds().interf_ff()[i].curves()[j];
        (
            nc.curve().cloned(),
            nc.range(),
            nc.pcurve1().cloned(),
            nc.pcurve2().cloned(),
        )
    };
    let Some(curve) = curve else {
        return;
    };
    let (Some(c2d0), Some(c2d1)) = (pc1, pc2) else {
        return;
    };
    let Some(s1) = BRepTool::face_surface(face1) else {
        return;
    };
    let Some(s2) = BRepTool::face_surface(face2) else {
        return;
    };
    let (t0, t1) = range;
    let ts = [t0, t1];
    let ps = [curve.d0(t0), curve.d0(t1)];
    let verts: Vec<usize> = leftover.into_iter().collect();
    for n_v in verts {
        let Some(v_shape) = f.ds().shape(n_v).cloned() else {
            continue;
        };
        let a_pv = BRepTool::vertex_point(&Vertex(v_shape));
        for m in 0..2 {
            if bnd[m].is_some() {
                continue;
            }
            let a_d2 = ps[m].square_distance(&a_pv);
            if a_d2 > STICK_DT2 {
                continue;
            }
            let mut nrm = [None, None];
            for n in 0..2 {
                let c2d = if n == 0 { c2d0.as_ref() } else { c2d1.as_ref() };
                let p2d = c2d.d0(ts[m]);
                let s = if n == 0 { s1.as_ref() } else { s2.as_ref() };
                nrm[n] = algo_tools3d::get_normal_to_surface(s, p2d.x(), p2d.y());
            }
            let (Some(n0), Some(n1)) = (nrm[0], nrm[1]) else {
                continue;
            };
            let mut sc = n0.dot(&n1);
            if sc < 0.0 {
                sc = -sc;
            }
            sc = 1.0 - sc;
            if sc > STICK_DSCPR {
                continue;
            }
            let a_d = a_d2.sqrt();
            put_pave_on_curve(f, i, j, n_v, a_d, mi, mv_tol, dmvlv, 0);
        }
    }
}

/// `BOPAlgo_PaveFiller::PutEFPavesOnCurve` (`_6.cxx:2692`).
pub fn put_ef_paves_on_curve(
    f: &mut PaveFiller,
    i: usize,
    j: usize,
    mi: &HashSet<usize>,
    mv_ef: &HashSet<usize>,
    mv_tol: &mut HashMap<usize, f64>,
    dmvlv: &mut HashMap<usize, Vec<usize>>,
) {
    if mv_ef.is_empty() {
        return;
    }
    let Some(curve) = f.ds().interf_ff()[i].curves()[j].curve().cloned() else {
        return;
    };
    if !int_tools_vertex_line::curve_is_bezier_or_bspline(curve.as_ref()) {
        return;
    }
    let mut leftover: HashSet<usize> = mv_ef.clone();
    remove_used_vertices(f.ds().interf_ff()[i].curves(), &mut leftover);
    if leftover.is_empty() {
        return;
    }
    for n_v in leftover {
        let Some(v_shape) = f.ds().shape(n_v).cloned() else {
            continue;
        };
        let a_pv = BRepTool::vertex_point(&Vertex(v_shape));
        let Some((_, a_dist)) =
            int_tools_vertex_line::project_point_on_curve_dist(curve.as_ref(), &a_pv)
        else {
            continue;
        };
        put_pave_on_curve(f, i, j, n_v, a_dist, mi, mv_tol, dmvlv, 2);
    }
}

/// `BOPAlgo_PaveFiller::EstimatePaveOnCurve` (`_6.cxx:4056`).
pub fn estimate_pave_on_curve(ds: &BopdsDS, n_v: usize, nc: &BopdsCurve, tol_r3d: f64) -> bool {
    let Some(v_s) = ds.shape(n_v) else {
        return false;
    };
    let Some(c) = nc.curve() else {
        return false;
    };
    int_tools_vertex_line::is_vertex_on_line(&Vertex(v_s.clone()), c.as_ref(), tol_r3d).is_some()
}
