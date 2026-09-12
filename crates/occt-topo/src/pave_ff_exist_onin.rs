//! ON/IN `IsExistingPaveBlock` (`_6.cxx:2047-2251`).
//!
//! AABB tree (`BOPTools_BoxTree`) on the first-point box. OCCT flags:
//! - `iFlag1/2 == 2` when the section vertex is shared with the candidate;
//! - `iFlag2 == 0` when the last-point box misses the candidate edge box;
//! - a common-block candidate inflates `aRealTol`, and a common-block-with-face
//!   doubles it again;
//! - both-ends-shared (`iFlag1 == iFlag2 == 2`) may run a tangent probe
//!   (`cos >= 0.9063`) and raise `aCoeff` to 2;
//! - both ends must `ComputePE` onto the candidate; the closest `aCoeff *
//!   aDistToSp` wins.

use occt_core::bnd::BndBox;
use occt_core::precision::RESOLUTION;
use occt_geom::Curve;

use crate::bopds::{BopdsDS, BopdsPaveBlock};
use crate::boptools_2d::intermediate_point;
use crate::brep_tool::BRepTool;
use crate::int_tools_full::IntToolsContext;
use crate::pave_ff_is_exist::thin_face_max_tol_add;
use crate::shape::{Edge, Vertex};

fn vertex_tol(ds: &BopdsDS, n_v: usize) -> f64 {
    ds.shape(n_v)
        .map(|s| BRepTool::vertex_tolerance(&Vertex(s.clone())))
        .unwrap_or(0.0)
}

fn tuples_eq(a: (usize, f64, f64), b: (usize, f64, f64)) -> bool {
    a.0 == b.0 && (a.1 - b.1).abs() <= 1e-7 && (a.2 - b.2).abs() <= 1e-7
}

fn pb_in_common(common: &[(usize, f64, f64)], pb: &BopdsPaveBlock) -> bool {
    let t = (pb.edge(), pb.first, pb.last);
    common.iter().copied().any(|x| tuples_eq(x, t))
}

fn is_geom_line(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let samples: [f64; 3] = if a.is_finite() && b.is_finite() && b > a {
        [a, 0.5 * (a + b), b]
    } else {
        [-1.0, 0.0, 1.0]
    };
    let t0 = c.d1(samples[0]).1;
    if t0.square_magnitude() < RESOLUTION * RESOLUTION {
        return false;
    }
    for &u in &samples[1..] {
        let t = c.d1(u).1;
        if t0.cross_magnitude(&t) > 1e-7 * t0.magnitude() * t.magnitude() {
            return false;
        }
    }
    true
}

fn edge_box_out(ds: &BopdsDS, n_e: usize, other: &BndBox) -> bool {
    ds.box_of(n_e).map(|b| b.is_out_box(other)).unwrap_or(false)
}

fn first_point_candidates<'a>(
    ds: &BopdsDS,
    on_in: &'a [BopdsPaveBlock],
    box_p1: &BndBox,
) -> Vec<&'a BopdsPaveBlock> {
    let mut tree = crate::bop_aabb_faces::AabbTree::new();
    tree.set_size(on_in.len());
    for (i, a_pb) in on_in.iter().enumerate() {
        if a_pb.edge() == 0 {
            continue;
        }
        let Some(b) = ds.box_of(a_pb.edge()) else {
            continue;
        };
        if b.is_void() {
            continue;
        }
        tree.add(i, crate::bnd_tools::bnd2bvh3d(&b));
    }
    tree.build();
    let mut candidates = Vec::new();
    for i in tree.select(box_p1) {
        if let Some(a_pb) = on_in.get(i) {
            if a_pb.edge() != 0 {
                candidates.push(a_pb);
            }
        }
    }
    candidates
}

fn flags_for_candidate(
    ds: &BopdsDS,
    n_v11: usize,
    n_v12: usize,
    n_v21: usize,
    n_v22: usize,
    n_e: usize,
    box_p2: &BndBox,
) -> (i32, i32) {
    let i_flag1: i32 = if n_v11 == n_v21 || n_v11 == n_v22 { 2 } else { 1 };
    let i_flag2: i32 = if n_v12 == n_v21 || n_v12 == n_v22 {
        2
    } else if !edge_box_out(ds, n_e, box_p2) {
        1
    } else {
        0
    };
    (i_flag1, i_flag2)
}

fn tangent_probe(
    ctx: &IntToolsContext,
    curve: &dyn Curve,
    a_sp: &Edge,
    pm: &occt_core::gp::GpPnt,
    v_tgt1: &occt_core::gp::GpVec,
    is_vtgt1_valid: bool,
    n_v11: usize,
    n_v12: usize,
    n_v21: usize,
    n_v22: usize,
    a_real_tol: f64,
    a_tol_v1: f64,
    a_tol_v2: f64,
    a_max_tol_add: f64,
) -> (i32, f64, f64, f64) {
    let skip = (n_v11 == n_v12 && n_v21 != n_v22) || (n_v11 != n_v12 && n_v21 == n_v22);
    if skip || !is_vtgt1_valid {
        return (1, 0.0, a_real_tol, 1.0);
    }
    let edge_is_line = BRepTool::edge_curve(a_sp)
        .map(|c| is_geom_line(c.as_ref()))
        .unwrap_or(false);
    if is_geom_line(curve) && edge_is_line {
        return (1, 0.0, a_real_tol, 1.0);
    }
    let a_tol_add = 2.0 * a_max_tol_add.min(a_real_tol.max(a_tol_v1.max(a_tol_v2)));
    let (st, tldp, dist) = ctx.compute_pe_pnt(pm, a_tol_add, a_sp);
    if st != 0 {
        return (st, dist, a_real_tol, 1.0);
    }
    if let Some(ec) = BRepTool::edge_curve(a_sp) {
        let (_, v_tgt2) = ec.d1(tldp);
        if v_tgt2.square_magnitude() > RESOLUTION {
            let a_cos = v_tgt1.dot(&v_tgt2.normalized());
            if a_cos.abs() >= 0.9063 {
                return (st, dist, a_tol_add, 2.0);
            }
        }
    }
    (st, dist, a_real_tol, 1.0)
}

/// ON/IN `IsExistingPaveBlock` (`_6.cxx:2047`).
pub fn is_existing_pave_block_on_in(
    ds: &BopdsDS,
    ctx: &IntToolsContext,
    the_pb: &BopdsPaveBlock,
    curve: &dyn Curve,
    the_tol_r3d: f64,
    the_mpb_on_in: &[BopdsPaveBlock],
    the_mpb_common: &[(usize, f64, f64)],
    fuzzy: f64,
) -> Option<(BopdsPaveBlock, f64)> {
    let (a_t1, a_t2) = the_pb.range();
    let (n_v11, n_v12) = the_pb.indices();
    let a_p1 = curve.d0(a_t1);
    let mut a_box_p1 = BndBox::new();
    a_box_p1.add_point(&a_p1);
    let a_tol_v11 = vertex_tol(ds, n_v11);
    a_box_p1.enlarge(a_tol_v11);

    let candidates = first_point_candidates(ds, the_mpb_on_in, &a_box_p1);
    if candidates.is_empty() {
        return None;
    }

    let a_tm = intermediate_point(a_t1, a_t2);
    let (a_pm, mut a_vtgt1) = curve.d1(a_tm);
    let mut a_box_pm = BndBox::new();
    a_box_pm.add_point(&a_pm);
    let is_vtgt1_valid = a_vtgt1.square_magnitude() > RESOLUTION;
    if is_vtgt1_valid {
        a_vtgt1.normalize();
    }

    let a_p2 = curve.d0(a_t2);
    let mut a_box_p2 = BndBox::new();
    a_box_p2.add_point(&a_p2);
    let a_tol_v12 = vertex_tol(ds, n_v12);
    a_box_p2.enlarge(a_tol_v12);

    let a_tol_v1 = a_tol_v11.max(a_tol_v12) + fuzzy;
    let a_tol_check = the_tol_r3d + fuzzy;
    let a_max_tol_add = thin_face_max_tol_add(a_tol_check);

    let mut found: Option<(BopdsPaveBlock, f64)> = None;
    let mut the_tol_new = f64::MAX;

    for a_pb in candidates {
        let (n_v21, n_v22) = a_pb.indices();
        let a_tol_v21 = vertex_tol(ds, n_v21);
        let a_tol_v22 = vertex_tol(ds, n_v22);
        let a_tol_v2 = a_tol_v21.max(a_tol_v22) + fuzzy;
        let n_e = a_pb.edge();
        let Some(sp_shape) = ds.shape(n_e) else {
            continue;
        };
        let a_sp = Edge(sp_shape.clone());
        let (mut i_flag1, mut i_flag2) =
            flags_for_candidate(ds, n_v11, n_v12, n_v21, n_v22, n_e, &a_box_p2);
        if i_flag2 == 0 {
            continue;
        }

        let mut a_coeff = 1.0;
        let mut a_dist_m1m2 = 0.0;
        let mut a_pe_status = 1i32;
        let mut a_real_tol = a_tol_check;
        if ds.is_common_block(a_pb) {
            a_real_tol = a_real_tol.max(a_tol_v1.max(a_tol_v2));
            if pb_in_common(the_mpb_common, a_pb) {
                a_real_tol *= 2.0;
            }
        } else if i_flag1 == 2 && i_flag2 == 2 {
            let (st, dist, real, coeff) = tangent_probe(
                ctx,
                curve,
                &a_sp,
                &a_pm,
                &a_vtgt1,
                is_vtgt1_valid,
                n_v11,
                n_v12,
                n_v21,
                n_v22,
                a_real_tol,
                a_tol_v1,
                a_tol_v2,
                a_max_tol_add,
            );
            a_pe_status = st;
            a_dist_m1m2 = dist;
            a_real_tol = real;
            a_coeff = coeff;
        }

        let mut a_box_tmp = a_box_pm;
        a_box_tmp.enlarge(a_real_tol);
        let mut a_dist_to_sp = 0.0;
        if edge_box_out(ds, n_e, &a_box_tmp) || a_pe_status < 0 {
            continue;
        } else if a_pe_status == 0 {
            a_dist_to_sp = a_dist_m1m2;
        } else if a_pe_status == 1 {
            let (st, _, dist) = ctx.compute_pe_pnt(&a_pm, a_real_tol, &a_sp);
            a_pe_status = st;
            if a_pe_status < 0 {
                continue;
            }
            a_dist_to_sp = dist;
        }
        if i_flag1 == 1 {
            let (st, _, dist) = ctx.compute_pe_pnt(&a_p1, a_real_tol, &a_sp);
            i_flag1 = if st == 0 { 1 } else { 0 };
            if i_flag1 != 0 && a_dist_to_sp < dist {
                a_dist_to_sp = dist;
            }
        }
        if i_flag2 == 1 {
            let (st, _, dist) = ctx.compute_pe_pnt(&a_p2, a_real_tol, &a_sp);
            i_flag2 = if st == 0 { 1 } else { 0 };
            if i_flag2 != 0 && a_dist_to_sp < dist {
                a_dist_to_sp = dist;
            }
        }
        if i_flag1 != 0 && i_flag2 != 0 && a_dist_to_sp < the_tol_new {
            the_tol_new = a_coeff * a_dist_to_sp;
            found = Some((a_pb.clone(), the_tol_new));
        }
    }
    found
}
