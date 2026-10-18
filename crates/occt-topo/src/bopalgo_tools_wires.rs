//! `BOPAlgo_Tools::EdgesToWires` / `WiresToFaces` and plane/tangent helpers.
//!
//! Source: `BOPAlgo_Tools.cxx` (`EdgesToWires` at 360, `WiresToFaces` at 665,
//! `MakeWires` at 805, `FindEdgeTangent` at 852 and 875, `FindPlane` at 910
//! and 976). Shared-edge glue (`BOPAlgo_Builder`) is left to the caller via
//! `shared`; when `shared` is false the edges are packed into a
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
use crate::topo_tools::{edges_of, is_same, vertices_of};

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
pub fn find_edge_tangent_curve(curve: &dyn Curve) -> Option<GpVec> {
    let t1 = curve.first_parameter();
    let t2 = curve.last_parameter();
    if !t1.is_finite() || !t2.is_finite() {
        return None;
    }
    const A_NB_P: i32 = 11;
    let dt = (t2 - t1) / A_NB_P as f64;
    let mut t = t1 + dt;
    while t <= t2 {
        let (_p, tgt) = curve.d1(t);
        if tgt.magnitude() > CONFUSION {
            return Some(tgt);
        }
        t += dt;
    }
    None
}

/// `FindEdgeTangent` (edge + cache, `BOPAlgo_Tools.cxx:852`).
pub fn find_edge_tangent(
    edge: &Edge,
    dm_edge_tgt: &mut HashMap<usize, GpDir>,
) -> Option<GpDir> {
    let key = shape_key(&edge.0);
    if let Some(d) = dm_edge_tgt.get(&key) {
        return Some(*d);
    }
    let curve = BRepTool::edge_curve(edge)?;
    let vte = find_edge_tangent_curve(curve.as_ref())?;
    let d = GpDir::new(vte.x(), vte.y(), vte.z()).ok()?;
    dm_edge_tgt.insert(key, d);
    Some(d)
}

/// `FindPlane` (curve overload, `BOPAlgo_Tools.cxx:910`).
pub fn find_plane_curve(curve: &dyn Curve) -> Option<GpPln> {
    let t1 = curve.first_parameter();
    let t2 = curve.last_parameter();
    if !t1.is_finite() || !t2.is_finite() {
        return None;
    }
    const A_NB_P: i32 = 11;
    let dt = (t2 - t1) / A_NB_P as f64;
    let (p1, v1) = curve.d1(t1);
    let mut t = t1 + dt;
    while t <= t2 {
        let (_p2, v2) = curve.d1(t);
        let vn = v1.crossed(&v2);
        if vn.magnitude() > CONFUSION {
            let d = GpDir::new(vn.x(), vn.y(), vn.z()).ok()?;
            return Some(pln_from_pnt_dir(&p1, d));
        }
        t += dt;
    }
    let _ = p1;
    None
}

/// `FindPlane` (wire overload, `BOPAlgo_Tools.cxx:976`).
pub fn find_plane_wire(
    wire: &TopoShape,
    dm_edge_tgt: &mut HashMap<usize, GpDir>,
    m_edges_no_unique_plane: &mut HashSet<usize>,
) -> Option<GpPln> {
    let edges = edges_of(wire);
    if edges.is_empty() {
        return None;
    }
    for e1 in &edges {
        let Some(dte1) = find_edge_tangent(e1, dm_edge_tgt) else {
            continue;
        };
        for e2 in &edges {
            if is_same(&e2.0, &e1.0) {
                continue;
            }
            let Some(dte2) = find_edge_tangent(e2, dm_edge_tgt) else {
                continue;
            };
            if dirs_parallel(&dte1, &dte2, ANGULAR) {
                continue;
            }
            let Ok(dn) = dte1.crossed(&dte2) else {
                continue;
            };
            let verts = vertices_of(&e1.0);
            let p = verts
                .first()
                .map(|v| BRepTool::vertex_point(v))
                .unwrap_or_default();
            return Some(pln_from_pnt_dir(&p, dn));
        }
    }
    for e in &edges {
        let k = shape_key(&e.0);
        if m_edges_no_unique_plane.contains(&k) {
            continue;
        }
        if let Some(c) = BRepTool::edge_curve(e) {
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
    edges_par: &[TopoShape],
    wires: &mut Vec<TopoShape>,
    check_unique_plane: bool,
    dm_edge_tgt: &mut HashMap<usize, GpDir>,
    m_edges_no_unique_plane: &mut HashSet<usize>,
) {
    if edges_par.is_empty() {
        return;
    }
    let blocks = make_connexity_edge_blocks(edges_par);
    let b = TopoBuilder::new();
    for cbe in blocks {
        if check_unique_plane {
            let dummy = TopoBuilder::new().make_compound_of(&cbe);
            if find_plane_wire(&dummy.0, dm_edge_tgt, m_edges_no_unique_plane).is_none() {
                continue;
            }
        }
        let edges: Vec<Edge> = cbe.iter().map(|s| Edge(s.clone())).collect();
        if edges.is_empty() {
            continue;
        }
        let w = b.make_wire(&edges);
        wires.push(w.0);
    }
}

/// `BOPAlgo_Tools::EdgesToWires` (`BOPAlgo_Tools.cxx:360`).
///
/// Returns `(error_code, wires_compound)`. `error_code == 1` means no edges;
/// `error_code == 2` is unused here (the nested boolean of `shared == false`
/// is the caller's responsibility).
pub fn edges_to_wires(
    edges: &TopoShape,
    shared: bool,
    ang_tol: f64,
) -> (i32, TopoShape) {
    let mut le: Vec<TopoShape> = Vec::new();
    for e in edges_of(edges) {
        if is_geometric_edge(&e) {
            le.push(e.0);
        }
    }
    let b = TopoBuilder::new();
    if le.is_empty() {
        return (1, b.make_compound_of(&[]).0);
    }
    if le.len() == 1 {
        let w = b.make_wire(&[Edge(le[0].clone())]);
        return (0, b.make_compound_of(&[w.0]).0);
    }
    let s_edges = if shared {
        b.make_compound_of(&le).0
    } else {
        b.make_compound_of(&le).0
    };
    let mut dm_edge_pln: HashMap<usize, GpPln> = HashMap::new();
    let mut dm_edge_tgt: HashMap<usize, GpDir> = HashMap::new();
    let mut m_no_plane: HashSet<usize> = HashSet::new();
    let mut l_edges: Vec<TopoShape> = Vec::new();
    for e in edges_of(&s_edges) {
        let Some(c) = BRepTool::edge_curve(&e) else {
            continue;
        };
        if let Some(pln) = find_plane_curve(c.as_ref()) {
            dm_edge_pln.insert(shape_key(&e.0), pln);
        } else if let Some(d) = find_edge_tangent(&e, &mut dm_edge_tgt) {
            dm_edge_tgt.insert(shape_key(&e.0), d);
            m_no_plane.insert(shape_key(&e.0));
            l_edges.push(e.0);
        }
    }
    let mut r_wires: Vec<TopoShape> = Vec::new();
    let mut lp_fence: Vec<GpDir> = Vec::new();
    let mut me_fence: HashSet<usize> = HashSet::new();
    let pln_keys: Vec<usize> = dm_edge_pln.keys().copied().collect();
    for (ii, &k_i) in pln_keys.iter().enumerate() {
        if !me_fence.insert(k_i) {
            continue;
        }
        let Some(ei) = s_edges_edge_by_key(&s_edges, k_i) else {
            continue;
        };
        let pln_i = dm_edge_pln[&k_i].clone();
        let di = pln_i.position().direction();
        lp_fence.push(di);
        let mut me_pln: Vec<TopoShape> = vec![ei.clone()];
        let mut mv: HashSet<usize> = vertices_of(&ei)
            .iter()
            .map(|v| shape_key(&v.0))
            .collect();
        for &k_j in pln_keys.iter().skip(ii + 1) {
            let dj = dm_edge_pln[&k_j].position().direction();
            if dirs_parallel(&di, &dj, ang_tol) {
                if let Some(ej) = s_edges_edge_by_key(&s_edges, k_j) {
                    me_pln.push(ej.clone());
                    me_fence.insert(k_j);
                    for v in vertices_of(&ej) {
                        mv.insert(shape_key(&v.0));
                    }
                }
            }
        }
        let mut ce_pln: Vec<TopoShape> = Vec::new();
        for (k_j, dj) in &dm_edge_tgt {
            if dirs_normal(&di, dj, ang_tol) {
                if let Some(e) = s_edges_edge_by_key(&s_edges, *k_j) {
                    ce_pln.push(e);
                }
            }
        }
        let blocks = make_connexity_edge_blocks(&ce_pln);
        for cbe in blocks {
            let mut b_add = cbe.iter().any(|s| {
                vertices_of(s)
                    .iter()
                    .any(|v| mv.contains(&shape_key(&v.0)))
            });
            if !b_add {
                let dummy = TopoBuilder::new().make_compound_of(&cbe);
                b_add = find_plane_wire(&dummy.0, &mut dm_edge_tgt, &mut m_no_plane).is_some();
            }
            if b_add {
                me_pln.extend(cbe);
            }
        }
        make_wires(
            &me_pln,
            &mut r_wires,
            false,
            &mut dm_edge_tgt,
            &mut m_no_plane,
        );
    }
    let mut dmve: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for e in &l_edges {
        for v in vertices_of(e) {
            dmve.entry(shape_key(&v.0)).or_default().push(e.clone());
        }
    }
    for lei in dmve.values() {
        if lei.len() < 2 {
            continue;
        }
        for ei1 in lei {
            let Some(di1) = dm_edge_tgt.get(&shape_key(ei1)).copied() else {
                continue;
            };
            for ei2 in lei {
                if is_same(ei2, ei1) {
                    continue;
                }
                let Some(di2) = dm_edge_tgt.get(&shape_key(ei2)).copied() else {
                    continue;
                };
                if dirs_parallel(&di1, &di2, ang_tol) {
                    continue;
                }
                let Ok(dni) = di1.crossed(&di2) else {
                    continue;
                };
                if lp_fence.iter().any(|d| dirs_parallel(&dni, d, ang_tol)) {
                    continue;
                }
                lp_fence.push(dni);
                let mut me_pln = vec![ei1.clone(), ei2.clone()];
                for (k_j, dj) in &dm_edge_tgt {
                    if dirs_normal(&dni, dj, ang_tol) {
                        if let Some(e) = s_edges_edge_by_key(&s_edges, *k_j) {
                            me_pln.push(e);
                        }
                    }
                }
                make_wires(
                    &me_pln,
                    &mut r_wires,
                    true,
                    &mut dm_edge_tgt,
                    &mut m_no_plane,
                );
            }
        }
    }
    let mut me_used: HashSet<usize> = HashSet::new();
    for w in &r_wires {
        for e in edges_of(w) {
            me_used.insert(shape_key(&e.0));
        }
    }
    let mut me_alone: Vec<TopoShape> = Vec::new();
    for (k, _) in &dm_edge_tgt {
        if !me_used.contains(k) {
            if let Some(e) = s_edges_edge_by_key(&s_edges, *k) {
                me_alone.push(e);
            }
        }
    }
    make_wires(
        &me_alone,
        &mut r_wires,
        false,
        &mut dm_edge_tgt,
        &mut m_no_plane,
    );
    (0, b.make_compound_of(&r_wires).0)
}

fn s_edges_edge_by_key(compound: &TopoShape, key: usize) -> Option<TopoShape> {
    edges_of(compound)
        .into_iter()
        .find(|e| shape_key(&e.0) == key)
        .map(|e| e.0)
}

/// `BOPAlgo_Tools::WiresToFaces` (`BOPAlgo_Tools.cxx:665`).
pub fn wires_to_faces(wires: &TopoShape, ang_tol: f64) -> (bool, TopoShape) {
    let mut m_fence: HashSet<usize> = HashSet::new();
    let mut r_faces: Vec<TopoShape> = Vec::new();
    let mut dm_edge_tgt: HashMap<usize, GpDir> = HashMap::new();
    let mut dm_wire_pln: Vec<(TopoShape, GpPln)> = Vec::new();
    let mut dm_wire_tol: HashMap<usize, f64> = HashMap::new();
    let mut m_no_plane: HashSet<usize> = HashSet::new();
    for w in ShapeIterator::of_shape(wires) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        if let Some(pln) = find_plane_wire(&w, &mut dm_edge_tgt, &mut m_no_plane) {
            dm_wire_tol.insert(shape_key(&w), max_edge_tolerance(&w));
            dm_wire_pln.push((w, pln));
        }
    }
    let nb = dm_wire_pln.len();
    for i in 0..nb {
        let wire_i = dm_wire_pln[i].0.clone();
        let k_i = shape_key(&wire_i);
        if m_fence.contains(&k_i) {
            continue;
        }
        let pln_i = dm_wire_pln[i].1.clone();
        let mut lw = vec![wire_i.clone()];
        m_fence.insert(k_i);
        let tol_i = dm_wire_tol.get(&k_i).copied().unwrap_or(0.0);
        for j in (i + 1)..nb {
            let wire_j = dm_wire_pln[j].0.clone();
            let k_j = shape_key(&wire_j);
            if m_fence.contains(&k_j) {
                continue;
            }
            let pln_j = &dm_wire_pln[j].1;
            if !dirs_parallel(
                &pln_i.position().direction(),
                &pln_j.position().direction(),
                ang_tol,
            ) {
                continue;
            }
            let dist = pln_i.location().distance(&pln_j.location());
            let tol_j = dm_wire_tol.get(&k_j).copied().unwrap_or(0.0);
            if dist > (tol_i + tol_j) {
                continue;
            }
            lw.push(wire_j);
            m_fence.insert(k_j);
        }
        let mut le: Vec<TopoShape> = Vec::new();
        for w in &lw {
            for e in edges_of(w) {
                let mut ef = e.0.clone();
                ef.set_orientation(Orientation::Forward);
                le.push(ef.clone());
                ef.set_orientation(Orientation::Reversed);
                le.push(ef);
            }
        }
        const A_MAX: f64 = 1.0e8;
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln_i.clone()));
        let mut ff = TopoBuilder::new().make_face(surf, &[]);
        ff.0.set_orientation(Orientation::Forward);
        let _ = A_MAX;
        let mut bf = FaceBuilder::new();
        bf.set_face(&ff);
        bf.set_shapes(&le);
        if bf.perform().is_err() {
            continue;
        }
        for fsp in bf.areas() {
            r_faces.push(fsp.clone());
        }
    }
    let out = TopoBuilder::new().make_compound_of(&r_faces).0;
    (!r_faces.is_empty(), out)
}

/// Bounding box of `shape` (used by classification helpers).
pub fn wires_shape_box(s: &TopoShape) -> occt_core::bnd::BndBox {
    shape_bbox(s)
}
