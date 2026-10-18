//! `CheckSelfInterference` — acquired self-intersection warnings.
//!
//! Source: `BOPAlgo_PaveFiller_11.cxx`. In self-interference mode
//! (`arguments.len() == 1`) the check is skipped. For multi-argument runs
//! the algorithm walks each argument range, records connections of new
//! vertices / IN and section edges onto argument faces, and warns when
//! several faces or edges of the same argument share the same extra vertex
//! or edge, or when a common block contains several original edges of one
//! argument.

use std::collections::{HashMap, HashSet};

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::pave_filler::PaveFiller;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::is_same;

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

fn add_connection(
    mcsi: &mut HashMap<usize, (TopoShape, Vec<TopoShape>)>,
    connector: &TopoShape,
    owner: &TopoShape,
) {
    let key = shape_key(connector);
    let entry = mcsi
        .entry(key)
        .or_insert_with(|| (connector.clone(), Vec::new()));
    if !entry.1.iter().any(|s| is_same(s, owner)) {
        entry.1.push(owner.clone());
    }
}

fn warn_compound(f: &mut PaveFiller, shapes: &[TopoShape]) {
    if shapes.len() < 2 {
        return;
    }
    let _wc = TopoBuilder::new().make_compound_of(shapes);
    f.add_warning(format!(
        "acquired self-interference: {} sub-shapes of one argument share a new vertex or edge",
        shapes.len()
    ));
}

/// `BOPAlgo_PaveFiller::CheckSelfInterference` (`BOPAlgo_PaveFiller_11.cxx`).
pub fn check_self_interference(f: &mut PaveFiller) -> Result<bool, String> {
    if f.arguments().len() == 1 {
        return Ok(false);
    }
    let mut detected = false;
    let nb_r = f.ds().nb_ranges();
    for i in 0..nb_r {
        let Some(r) = f.ds().range(i) else {
            continue;
        };
        let mut mcsi: HashMap<usize, (TopoShape, Vec<TopoShape>)> = HashMap::new();
        let mut mcb_fence: HashSet<Vec<usize>> = HashSet::new();
        let r_last = r.last();
        let mut j = r.first();
        while j <= r_last {
            let Some(si) = f.ds().shape_info(j) else {
                j += 1;
                continue;
            };
            if !f.ds().has_pave_blocks(j) && f.ds().face_info(j).is_none() {
                j += 1;
                continue;
            }
            let Some(s) = f.ds().shape(j).cloned() else {
                j += 1;
                continue;
            };
            match si.shape_type() {
                ShapeType::Edge => {
                    if si.has_flag() {
                        j += 1;
                        continue;
                    }
                    let mut m_sub: HashSet<usize> = HashSet::new();
                    for &n_v in si.sub_shapes() {
                        let n_v = f.ds().has_shape_sd(n_v).unwrap_or(n_v);
                        m_sub.insert(n_v);
                    }
                    let lpb = f.ds().pave_blocks(j).to_vec();
                    let b_analyze_v = lpb.len() > 1;
                    for pb in &lpb {
                        if b_analyze_v {
                            let (n0, n1) = pb.indices();
                            for n_vk in [n0, n1] {
                                if !r.contains(n_vk) && !m_sub.contains(&n_vk) {
                                    if let Some(v) = f.ds().shape(n_vk).cloned() {
                                        add_connection(&mut mcsi, &v, &s);
                                    }
                                }
                            }
                        }
                        if f.ds().is_common_block(pb) {
                            if let Some(cb) = f.ds().common_block(pb) {
                                let pbs = if cb.pave_blocks().is_empty() {
                                    vec![pb.clone()]
                                } else {
                                    cb.pave_blocks().to_vec()
                                };
                                let mut fence: Vec<usize> =
                                    pbs.iter().map(|p| p.original_edge()).collect();
                                fence.sort_unstable();
                                fence.dedup();
                                if mcb_fence.insert(fence) {
                                    let mut le: Vec<usize> = Vec::new();
                                    for pbcb in &pbs {
                                        let n_e_or = pbcb.original_edge();
                                        if r.contains(n_e_or) {
                                            le.push(n_e_or);
                                        }
                                    }
                                    if le.len() > 1 {
                                        let mut wc: Vec<TopoShape> = Vec::new();
                                        for n_e1 in le {
                                            if let Some(e1) = f.ds().shape(n_e1).cloned() {
                                                wc.push(e1);
                                            }
                                        }
                                        detected = true;
                                        warn_compound(f, &wc);
                                    }
                                }
                            }
                        }
                    }
                }
                ShapeType::Face => {
                    if let Some(fi) = f.ds().face_info(j) {
                        for &n_v in fi.verts_in.iter().chain(fi.verts_sc.iter()) {
                            if let Some(v) = f.ds().shape(n_v).cloned() {
                                add_connection(&mut mcsi, &v, &s);
                            }
                        }
                        let pbs: Vec<(usize, f64, f64)> = fi
                            .paves_in()
                            .iter()
                            .chain(fi.paves().iter())
                            .copied()
                            .collect();
                        for (n_e, _, _) in pbs {
                            if n_e == 0 {
                                continue;
                            }
                            if let Some(e) = f.ds().shape(n_e).cloned() {
                                add_connection(&mut mcsi, &e, &s);
                            }
                        }
                    }
                }
                _ => {}
            }
            j += 1;
        }
        for (_k, (_conn, owners)) in mcsi {
            if owners.len() > 1 {
                detected = true;
                warn_compound(f, &owners);
            }
        }
    }
    Ok(detected)
}
