//! `BOPAlgo_Tools::EdgesToWires` / `WiresToFaces` and plane/tangent helpers.
//!
//! Source: `BOPAlgo_Tools.cxx` (`EdgesToWires` at 360, `WiresToFaces` at 665,
//! `MakeWires` at 805, `FindEdgeTangent` at 852 and 875, `FindPlane` at 910
//! and 976). Shared-edge glue (`BOPAlgo_Builder`) is left to the caller via
//! `the_shared`; when `the_shared` is false the edges are packed into a
//! compound without a nested boolean.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::{Curve, GeomPlane, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::builder_area::AreaBuilder;
use crate::builder_face::FaceBuilder;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of, vertices_of};

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

fn dirs_parallel(a: &GpDir, b: &GpDir, ang: f64) -> bool {
    let a0 = a.angle(b);
    a0 <= ang || (std::f64::consts::PI - a0) <= ang
}

fn dirs_normal(a: &GpDir, b: &GpDir, ang: f64) -> bool {
    (a.angle(b) - std::f64::consts::FRAC_PI_2).abs() <= ang
}

fn pln_from_pnt_dir(p: &GpPnt, d: GpDir) -> GpPln {
    GpPln::new(GpAx3::from_ax1(&GpAx1::new(*p, d)))
}

fn max_edge_tolerance(wire: &TopoShape) -> f64 {
    edges_of(wire)
        .iter()
        .map(|e| BRepTool::edge_tolerance(e))
        .fold(0.0, f64::max)
}

fn is_geometric_edge(e: &Edge) -> bool {
    !BRepTool::is_degenerated(e) && BRepTool::edge_curve(e).is_some()
}

/// `FindEdgeTangent` (curve overload, `BOPAlgo_Tools.cxx:875`).
pub fn find_edge_tangent_curve(the_curve: &dyn Curve) -> Option<GpVec> {
    let a_t1 = the_curve.first_parameter();
    let a_t2 = the_curve.last_parameter();
    if !a_t1.is_finite() || !a_t2.is_finite() {
        return None;
    }
    const A_NB_P: i32 = 11;
    let a_dt = (a_t2 - a_t1) / A_NB_P as f64;
    let mut a_t = a_t1 + a_dt;
    while a_t <= a_t2 {
        let (_p, tgt) = the_curve.d1(a_t);
        if tgt.magnitude() > CONFUSION {
            return Some(tgt);
        }
        a_t += a_dt;
    }
    None
}

/// `FindEdgeTangent` (edge + cache, `BOPAlgo_Tools.cxx:852`).
pub fn find_edge_tangent(
    the_edge: &Edge,
    dm_edge_tgt: &mut HashMap<usize, GpDir>,
) -> Option<GpDir> {
    let key = shape_key(&the_edge.0);
    if let Some(d) = dm_edge_tgt.get(&key) {
        return Some(*d);
    }
    let curve = BRepTool::edge_curve(the_edge)?;
    let a_vte = find_edge_tangent_curve(curve.as_ref())?;
    let d = GpDir::new(a_vte.x(), a_vte.y(), a_vte.z()).ok()?;
    dm_edge_tgt.insert(key, d);
    Some(d)
}

/// `FindPlane` (curve overload, `BOPAlgo_Tools.cxx:910`).
pub fn find_plane_curve(the_curve: &dyn Curve) -> Option<GpPln> {
    let a_t1 = the_curve.first_parameter();
    let a_t2 = the_curve.last_parameter();
    if !a_t1.is_finite() || !a_t2.is_finite() {
        return None;
    }
    const A_NB_P: i32 = 11;
    let a_dt = (a_t2 - a_t1) / A_NB_P as f64;
    let (p1, a_v1) = the_curve.d1(a_t1);
    let mut a_t = a_t1 + a_dt;
    while a_t <= a_t2 {
        let (_p2, a_v2) = the_curve.d1(a_t);
        let a_vn = a_v1.crossed(&a_v2);
        if a_vn.magnitude() > CONFUSION {
            let d = GpDir::new(a_vn.x(), a_vn.y(), a_vn.z()).ok()?;
            return Some(pln_from_pnt_dir(&p1, d));
        }
        a_t += a_dt;
    }
    let _ = p1;
    None
}

/// `FindPlane` (wire overload, `BOPAlgo_Tools.cxx:976`).
pub fn find_plane_wire(
    the_wire: &TopoShape,
    dm_edge_tgt: &mut HashMap<usize, GpDir>,
    m_edges_no_unique_plane: &mut HashSet<usize>,
) -> Option<GpPln> {
    let edges = edges_of(the_wire);
    if edges.is_empty() {
        return None;
    }
    for a_e1 in &edges {
        let Some(a_dte1) = find_edge_tangent(a_e1, dm_edge_tgt) else {
            continue;
        };
        for a_e2 in &edges {
            if a_e2.0.same_tshape(&a_e1.0) {
                continue;
            }
            let Some(a_dte2) = find_edge_tangent(a_e2, dm_edge_tgt) else {
                continue;
            };
            if dirs_parallel(&a_dte1, &a_dte2, ANGULAR) {
                continue;
            }
            let Ok(a_dn) = a_dte1.crossed(&a_dte2) else {
                continue;
            };
            let verts = vertices_of(&a_e1.0);
            let p = verts
                .first()
                .map(|v| BRepTool::vertex_point(v))
                .unwrap_or_default();
            return Some(pln_from_pnt_dir(&p, a_dn));
        }
    }
    for a_e in &edges {
        let k = shape_key(&a_e.0);
        if m_edges_no_unique_plane.contains(&k) {
            continue;
        }
        if let Some(c) = BRepTool::edge_curve(a_e) {
            if let Some(pln) = find_plane_curve(c.as_ref()) {
                return Some(pln);
            }
        }
        m_edges_no_unique_plane.insert(k);
    }
    None
}

fn make_connexity_edge_blocks(edges: &[TopoShape]) -> Vec<Vec<TopoShape>> {
    let mut v2e: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        for v in vertices_of(e) {
            v2e.entry(shape_key(&v.0)).or_default().push(i);
        }
    }
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for idxs in v2e.values() {
        for &a in idxs {
            for &b in idxs {
                if a != b {
                    adj.entry(a).or_default().push(b);
                }
            }
        }
    }
    let mut seen = HashSet::new();
    let mut blocks = Vec::new();
    for i in 0..edges.len() {
        if !seen.insert(i) {
            continue;
        }
        let mut stack = vec![i];
        let mut block = Vec::new();
        while let Some(n) = stack.pop() {
            block.push(edges[n].clone());
            if let Some(nbrs) = adj.get(&n) {
                for &q in nbrs {
                    if seen.insert(q) {
                        stack.push(q);
                    }
                }
            }
        }
        blocks.push(block);
    }
    blocks
}

/// `MakeWires` (`BOPAlgo_Tools.cxx:805`).
pub fn make_wires(
    the_edges: &[TopoShape],
    the_wires: &mut Vec<TopoShape>,
    the_check_unique_plane: bool,
    dm_edge_tgt: &mut HashMap<usize, GpDir>,
    m_edges_no_unique_plane: &mut HashSet<usize>,
) {
    if the_edges.is_empty() {
        return;
    }
    let blocks = make_connexity_edge_blocks(the_edges);
    let b = TopoBuilder::new();
    for a_cbe in blocks {
        if the_check_unique_plane {
            let dummy = TopoBuilder::new().make_compound_of(&a_cbe);
            if find_plane_wire(&dummy.0, dm_edge_tgt, m_edges_no_unique_plane).is_none() {
                continue;
            }
        }
        let edges: Vec<Edge> = a_cbe.iter().map(|s| Edge(s.clone())).collect();
        if edges.is_empty() {
            continue;
        }
        let w = b.make_wire(&edges);
        the_wires.push(w.0);
    }
}

/// `BOPAlgo_Tools::EdgesToWires` (`BOPAlgo_Tools.cxx:360`).
///
/// Returns `(error_code, wires_compound)`. `error_code == 1` means no edges;
/// `error_code == 2` is unused here (the nested boolean of `the_shared == false`
/// is the caller's responsibility).
pub fn edges_to_wires(
    the_edges: &TopoShape,
    the_shared: bool,
    the_ang_tol: f64,
) -> (i32, TopoShape) {
    let mut a_le: Vec<TopoShape> = Vec::new();
    for e in edges_of(the_edges) {
        if is_geometric_edge(&e) {
            a_le.push(e.0);
        }
    }
    let b = TopoBuilder::new();
    if a_le.is_empty() {
        return (1, b.make_compound_of(&[]).0);
    }
    if a_le.len() == 1 {
        let w = b.make_wire(&[Edge(a_le[0].clone())]);
        return (0, b.make_compound_of(&[w.0]).0);
    }
    let a_s_edges = if the_shared {
        b.make_compound_of(&a_le).0
    } else {
        b.make_compound_of(&a_le).0
    };
    let mut dm_edge_pln: HashMap<usize, GpPln> = HashMap::new();
    let mut dm_edge_tgt: HashMap<usize, GpDir> = HashMap::new();
    let mut a_m_no_plane: HashSet<usize> = HashSet::new();
    let mut a_l_edges: Vec<TopoShape> = Vec::new();
    for e in edges_of(&a_s_edges) {
        let Some(c) = BRepTool::edge_curve(&e) else {
            continue;
        };
        if let Some(pln) = find_plane_curve(c.as_ref()) {
            dm_edge_pln.insert(shape_key(&e.0), pln);
        } else if let Some(d) = find_edge_tangent(&e, &mut dm_edge_tgt) {
            dm_edge_tgt.insert(shape_key(&e.0), d);
            a_m_no_plane.insert(shape_key(&e.0));
            a_l_edges.push(e.0);
        }
    }
    let mut a_r_wires: Vec<TopoShape> = Vec::new();
    let mut a_lp_fence: Vec<GpDir> = Vec::new();
    let mut a_me_fence: HashSet<usize> = HashSet::new();
    let pln_keys: Vec<usize> = dm_edge_pln.keys().copied().collect();
    for (ii, &k_i) in pln_keys.iter().enumerate() {
        if !a_me_fence.insert(k_i) {
            continue;
        }
        let Some(a_ei) = a_s_edges_edge_by_key(&a_s_edges, k_i) else {
            continue;
        };
        let a_pln_i = dm_edge_pln[&k_i].clone();
        let a_di = a_pln_i.position().direction();
        a_lp_fence.push(a_di);
        let mut a_me_pln: Vec<TopoShape> = vec![a_ei.clone()];
        let mut a_mv: HashSet<usize> = vertices_of(&a_ei)
            .iter()
            .map(|v| shape_key(&v.0))
            .collect();
        for &k_j in pln_keys.iter().skip(ii + 1) {
            let a_dj = dm_edge_pln[&k_j].position().direction();
            if dirs_parallel(&a_di, &a_dj, the_ang_tol) {
                if let Some(a_ej) = a_s_edges_edge_by_key(&a_s_edges, k_j) {
                    a_me_pln.push(a_ej.clone());
                    a_me_fence.insert(k_j);
                    for v in vertices_of(&a_ej) {
                        a_mv.insert(shape_key(&v.0));
                    }
                }
            }
        }
        let mut a_ce_pln: Vec<TopoShape> = Vec::new();
        for (k_j, a_dj) in &dm_edge_tgt {
            if dirs_normal(&a_di, a_dj, the_ang_tol) {
                if let Some(e) = a_s_edges_edge_by_key(&a_s_edges, *k_j) {
                    a_ce_pln.push(e);
                }
            }
        }
        let blocks = make_connexity_edge_blocks(&a_ce_pln);
        for a_cbe in blocks {
            let mut b_add = a_cbe.iter().any(|s| {
                vertices_of(s)
                    .iter()
                    .any(|v| a_mv.contains(&shape_key(&v.0)))
            });
            if !b_add {
                let dummy = TopoBuilder::new().make_compound_of(&a_cbe);
                b_add = find_plane_wire(&dummy.0, &mut dm_edge_tgt, &mut a_m_no_plane).is_some();
            }
            if b_add {
                a_me_pln.extend(a_cbe);
            }
        }
        make_wires(
            &a_me_pln,
            &mut a_r_wires,
            false,
            &mut dm_edge_tgt,
            &mut a_m_no_plane,
        );
    }
    let mut a_dmve: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for e in &a_l_edges {
        for v in vertices_of(e) {
            a_dmve.entry(shape_key(&v.0)).or_default().push(e.clone());
        }
    }
    for a_lei in a_dmve.values() {
        if a_lei.len() < 2 {
            continue;
        }
        for a_ei1 in a_lei {
            let Some(a_di1) = dm_edge_tgt.get(&shape_key(a_ei1)).copied() else {
                continue;
            };
            for a_ei2 in a_lei {
                if a_ei2.same_tshape(a_ei1) {
                    continue;
                }
                let Some(a_di2) = dm_edge_tgt.get(&shape_key(a_ei2)).copied() else {
                    continue;
                };
                if dirs_parallel(&a_di1, &a_di2, the_ang_tol) {
                    continue;
                }
                let Ok(a_dni) = a_di1.crossed(&a_di2) else {
                    continue;
                };
                if a_lp_fence.iter().any(|d| dirs_parallel(&a_dni, d, the_ang_tol)) {
                    continue;
                }
                a_lp_fence.push(a_dni);
                let mut a_me_pln = vec![a_ei1.clone(), a_ei2.clone()];
                for (k_j, a_dj) in &dm_edge_tgt {
                    if dirs_normal(&a_dni, a_dj, the_ang_tol) {
                        if let Some(e) = a_s_edges_edge_by_key(&a_s_edges, *k_j) {
                            a_me_pln.push(e);
                        }
                    }
                }
                make_wires(
                    &a_me_pln,
                    &mut a_r_wires,
                    true,
                    &mut dm_edge_tgt,
                    &mut a_m_no_plane,
                );
            }
        }
    }
    let mut a_me_used: HashSet<usize> = HashSet::new();
    for w in &a_r_wires {
        for e in edges_of(w) {
            a_me_used.insert(shape_key(&e.0));
        }
    }
    let mut a_me_alone: Vec<TopoShape> = Vec::new();
    for (k, _) in &dm_edge_tgt {
        if !a_me_used.contains(k) {
            if let Some(e) = a_s_edges_edge_by_key(&a_s_edges, *k) {
                a_me_alone.push(e);
            }
        }
    }
    make_wires(
        &a_me_alone,
        &mut a_r_wires,
        false,
        &mut dm_edge_tgt,
        &mut a_m_no_plane,
    );
    (0, b.make_compound_of(&a_r_wires).0)
}

fn a_s_edges_edge_by_key(compound: &TopoShape, key: usize) -> Option<TopoShape> {
    edges_of(compound)
        .into_iter()
        .find(|e| shape_key(&e.0) == key)
        .map(|e| e.0)
}

/// `BOPAlgo_Tools::WiresToFaces` (`BOPAlgo_Tools.cxx:665`).
pub fn wires_to_faces(the_wires: &TopoShape, the_ang_tol: f64) -> (bool, TopoShape) {
    let mut a_m_fence: HashSet<usize> = HashSet::new();
    let mut a_r_faces: Vec<TopoShape> = Vec::new();
    let mut dm_edge_tgt: HashMap<usize, GpDir> = HashMap::new();
    let mut dm_wire_pln: Vec<(TopoShape, GpPln)> = Vec::new();
    let mut dm_wire_tol: HashMap<usize, f64> = HashMap::new();
    let mut a_m_no_plane: HashSet<usize> = HashSet::new();
    for w in ShapeIterator::of_shape(the_wires) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        if let Some(pln) = find_plane_wire(&w, &mut dm_edge_tgt, &mut a_m_no_plane) {
            dm_wire_tol.insert(shape_key(&w), max_edge_tolerance(&w));
            dm_wire_pln.push((w, pln));
        }
    }
    let a_nb = dm_wire_pln.len();
    for i in 0..a_nb {
        let a_wire_i = dm_wire_pln[i].0.clone();
        let k_i = shape_key(&a_wire_i);
        if a_m_fence.contains(&k_i) {
            continue;
        }
        let a_pln_i = dm_wire_pln[i].1.clone();
        let mut a_lw = vec![a_wire_i.clone()];
        a_m_fence.insert(k_i);
        let a_tol_i = dm_wire_tol.get(&k_i).copied().unwrap_or(0.0);
        for j in (i + 1)..a_nb {
            let a_wire_j = dm_wire_pln[j].0.clone();
            let k_j = shape_key(&a_wire_j);
            if a_m_fence.contains(&k_j) {
                continue;
            }
            let a_pln_j = &dm_wire_pln[j].1;
            if !dirs_parallel(
                &a_pln_i.position().direction(),
                &a_pln_j.position().direction(),
                the_ang_tol,
            ) {
                continue;
            }
            let a_dist = a_pln_i.location().distance(&a_pln_j.location());
            let a_tol_j = dm_wire_tol.get(&k_j).copied().unwrap_or(0.0);
            if a_dist > (a_tol_i + a_tol_j) {
                continue;
            }
            a_lw.push(a_wire_j);
            a_m_fence.insert(k_j);
        }
        let mut a_le: Vec<TopoShape> = Vec::new();
        for w in &a_lw {
            for e in edges_of(w) {
                let mut ef = e.0.clone();
                ef.set_orientation(Orientation::Forward);
                a_le.push(ef.clone());
                ef.set_orientation(Orientation::Reversed);
                a_le.push(ef);
            }
        }
        const A_MAX: f64 = 1.0e8;
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(a_pln_i.clone()));
        let mut a_ff = TopoBuilder::new().make_face(surf, &[]);
        a_ff.0.set_orientation(Orientation::Forward);
        let _ = A_MAX;
        let mut a_bf = FaceBuilder::new();
        a_bf.set_face(&a_ff);
        a_bf.set_shapes(&a_le);
        if a_bf.perform().is_err() {
            continue;
        }
        for a_fsp in a_bf.areas() {
            a_r_faces.push(a_fsp.clone());
        }
    }
    let out = TopoBuilder::new().make_compound_of(&a_r_faces).0;
    (!a_r_faces.is_empty(), out)
}

/// Bounding box of `shape` (used by classification helpers).
pub fn shape_box_of(s: &TopoShape) -> occt_core::bnd::BndBox {
    shape_bbox(s)
}
