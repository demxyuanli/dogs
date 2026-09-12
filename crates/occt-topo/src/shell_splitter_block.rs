//! `BOPAlgo_ShellSplitter::SplitBlock` (`BOPAlgo_ShellSplitter.cxx:153`).
//!
//! Irregular connexity blocks (a TShape present FORWARD and REVERSED, or an
//! edge with other than two faces) are walked with `GetEdgeOff` / `GetFaceOff`.
//! Closed shells are kept; `RefineShell` splits a walk that still has a
//! multi-connected or same-orientation edge.

use std::collections::{HashMap, HashSet};

use crate::abs::Orientation;
use crate::algo_tools_face::{get_edge_off, get_face_off, CoupleOfShape};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::connexity_block::ConnexityBlock;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, Shell, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, wires_of_face};

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

fn ori_byte(o: Orientation) -> u8 {
    match o {
        Orientation::Forward => 0,
        Orientation::Reversed => 1,
        Orientation::Internal => 2,
        Orientation::External => 3,
    }
}

fn ori_key(s: &TopoShape) -> (usize, u8) {
    (shape_key(s), ori_byte(s.orientation()))
}

fn face_boundary_edges(face: &TopoShape) -> Vec<Edge> {
    let mut out = Vec::new();
    for w in wires_of_face(&Face(face.clone())) {
        out.extend(edges_of_wire(&w));
    }
    out
}

fn make_shell(faces: &[TopoShape]) -> TopoShape {
    let bld = TopoBuilder::new();
    let mut shell = Shell::new();
    for f in faces {
        bld.add(&mut shell.0, f);
    }
    shell.0
}

fn map_edges_and_faces(faces: &[TopoShape]) -> HashMap<usize, Vec<TopoShape>> {
    let mut mef: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for f in faces {
        for e in face_boundary_edges(f) {
            let k = shape_key(&e.0);
            let ent = mef.entry(k).or_default();
            if !ent.iter().any(|x| ori_key(x) == ori_key(f)) {
                ent.push(f.clone());
            }
        }
    }
    mef
}

fn merge_face_into_mef(mef: &mut HashMap<usize, Vec<TopoShape>>, face: &TopoShape) {
    for e in face_boundary_edges(face) {
        let k = shape_key(&e.0);
        let ent = mef.entry(k).or_default();
        if !ent.iter().any(|x| ori_key(x) == ori_key(face)) {
            ent.push(face.clone());
        }
    }
}

fn shell_is_closed(shell: &TopoShape) -> bool {
    !crate::bop_build_solids::geometrically_open(shell)
}

/// `RefineShell` (`ShellSplitter.cxx:443`). Empty stop-edge set returns the
/// input shell. Otherwise the walk does not cross stop edges.
fn refine_shell(
    shell: &TopoShape,
    mef: &HashMap<usize, Vec<TopoShape>>,
) -> Vec<TopoShape> {
    let faces: Vec<TopoShape> = crate::bop_occt_util::iter_children(shell)
        .into_iter()
        .filter(|c| c.is_face())
        .collect();
    if faces.is_empty() {
        return Vec::new();
    }
    let mut stop: HashSet<usize> = HashSet::new();
    for (ek, lf) in mef {
        if lf.len() > 2 {
            stop.insert(*ek);
            continue;
        }
        if lf.len() == 2 {
            let e1 = face_boundary_edges(&lf[0])
                .into_iter()
                .find(|e| shape_key(&e.0) == *ek)
                .map(|e| e.0);
            let e2 = face_boundary_edges(&lf[1])
                .into_iter()
                .find(|e| shape_key(&e.0) == *ek)
                .map(|e| e.0);
            if let (Some(a), Some(b)) = (e1, e2) {
                if a.orientation() == b.orientation() {
                    stop.insert(*ek);
                    continue;
                }
            }
        }
        let mut nb = 0usize;
        for f in lf {
            nb += 1;
            if face_boundary_edges(f).iter().any(|e| {
                shape_key(&e.0) == *ek && e.0.orientation() == Orientation::Internal
            }) {
                nb += 1;
            }
        }
        if nb > 2 {
            stop.insert(*ek);
        }
    }
    if stop.is_empty() {
        return vec![shell.clone()];
    }
    let mut processed: HashSet<(usize, u8)> = HashSet::new();
    let mut out: Vec<TopoShape> = Vec::new();
    for f1 in &faces {
        if !processed.insert(ori_key(f1)) {
            continue;
        }
        let mut block: Vec<TopoShape> = vec![f1.clone()];
        let mut in_block: HashSet<(usize, u8)> = HashSet::new();
        in_block.insert(ori_key(f1));
        let mut wave = vec![f1.clone()];
        loop {
            let mut next: Vec<TopoShape> = Vec::new();
            for fp in &wave {
                for e in face_boundary_edges(fp) {
                    if stop.contains(&shape_key(&e.0)) {
                        continue;
                    }
                    if e.0.orientation() == Orientation::Internal {
                        continue;
                    }
                    if BRepTool::is_degenerated(&e) {
                        continue;
                    }
                    let Some(lf) = mef.get(&shape_key(&e.0)) else {
                        continue;
                    };
                    for fp1 in lf {
                        if fp1.same_tshape(fp) {
                            continue;
                        }
                        if in_block.contains(&ori_key(fp1)) {
                            continue;
                        }
                        if processed.insert(ori_key(fp1)) {
                            in_block.insert(ori_key(fp1));
                            block.push(fp1.clone());
                            next.push(fp1.clone());
                        }
                    }
                }
            }
            if next.is_empty() {
                break;
            }
            wave = next;
        }
        if !block.is_empty() {
            out.push(make_shell(&block));
        }
    }
    out
}

/// `BOPAlgo_ShellSplitter::SplitBlock`.
pub fn split_block(block: &mut ConnexityBlock) {
    block.change_loops_mut().clear();
    let my_shapes = block.shapes().to_vec();
    let mut a_m_faces: HashMap<(usize, u8), TopoShape> = HashMap::new();
    for f in &my_shapes {
        a_m_faces.insert(ori_key(f), f.clone());
    }
    loop {
        let live: Vec<TopoShape> = a_m_faces.values().cloned().collect();
        let a_ef = map_edges_and_faces(&live);
        let begin = a_m_faces.len();
        for (ek, lf) in &a_ef {
            let Some(e_rep) = live.iter().find_map(|f| {
                face_boundary_edges(f)
                    .into_iter()
                    .find(|e| shape_key(&e.0) == *ek)
            }) else {
                continue;
            };
            if BRepTool::is_degenerated(&e_rep) || e_rep.0.orientation() == Orientation::Internal {
                continue;
            }
            if lf.len() == 1 {
                a_m_faces.remove(&ori_key(&lf[0]));
            }
        }
        let end = a_m_faces.len();
        if end == begin || end == 0 {
            break;
        }
    }
    if a_m_faces.is_empty() {
        return;
    }
    let mut a_lf_connected: Vec<TopoShape> = Vec::new();
    let mut a_boundary: HashSet<usize> = HashSet::new();
    for f in &my_shapes {
        if a_m_faces.contains_key(&ori_key(f)) {
            a_lf_connected.push(f.clone());
            if !a_boundary.insert(shape_key(f)) {
                a_boundary.remove(&shape_key(f));
            }
        }
    }
    let a_nb_shapes = a_lf_connected.len();
    let live: Vec<TopoShape> = a_m_faces.values().cloned().collect();
    let a_ef_map = map_edges_and_faces(&live);
    let mut added: HashSet<(usize, u8)> = HashSet::new();
    let mut ctx = IntToolsContext::new();
    let mut all_taken = false;
    for a_ff in &a_lf_connected {
        if all_taken {
            break;
        }
        if !added.insert(ori_key(a_ff)) {
            continue;
        }
        let mut shell_faces: Vec<TopoShape> = vec![a_ff.clone()];
        let mut a_mefp: HashMap<usize, Vec<TopoShape>> = HashMap::new();
        merge_face_into_mef(&mut a_mefp, a_ff);
        let mut i = 0usize;
        while i < shell_faces.len() {
            let a_f = shell_faces[i].clone();
            let is_boundary = a_boundary.contains(&shape_key(&a_f));
            for a_e in face_boundary_edges(&a_f) {
                if a_mefp
                    .get(&shape_key(&a_e.0))
                    .map(|l| l.len() > 1)
                    .unwrap_or(false)
                {
                    continue;
                }
                if a_e.0.orientation() == Orientation::Internal {
                    continue;
                }
                if BRepTool::is_degenerated(&a_e) {
                    continue;
                }
                let Some(a_lf) = a_ef_map.get(&shape_key(&a_e.0)) else {
                    continue;
                };
                if a_lf.is_empty() {
                    continue;
                }
                let mut a_lcs_off: Vec<CoupleOfShape> = Vec::new();
                let mut a_nb_ways_inside = 0i32;
                let mut a_sel: Option<TopoShape> = None;
                for a_fl in a_lf {
                    if a_f.same_tshape(a_fl) || added.contains(&ori_key(a_fl)) {
                        continue;
                    }
                    let Some(a_el) = get_edge_off(&a_e, &Face(a_fl.clone())) else {
                        continue;
                    };
                    if is_boundary && !a_boundary.contains(&shape_key(a_fl)) {
                        a_nb_ways_inside += 1;
                        a_sel = Some(a_fl.clone());
                    }
                    a_lcs_off.push(CoupleOfShape::new(a_el.0, a_fl.clone()));
                }
                if a_lcs_off.is_empty() {
                    continue;
                }
                if !is_boundary || a_nb_ways_inside != 1 {
                    if a_lcs_off.len() == 1 {
                        a_sel = Some(a_lcs_off[0].shape2.clone());
                    } else if let Some((f, _)) =
                        get_face_off(&a_e, &Face(a_f.clone()), &a_lcs_off, &ctx)
                    {
                        a_sel = Some(f.0);
                    }
                }
                if let Some(sel) = a_sel {
                    if added.insert(ori_key(&sel)) {
                        shell_faces.push(sel.clone());
                        merge_face_into_mef(&mut a_mefp, &sel);
                    }
                }
            }
            i += 1;
        }
        let walked = make_shell(&shell_faces);
        let mut a_l_sh_nc: Vec<TopoShape> = Vec::new();
        let refined = refine_shell(&walked, &a_mefp);
        let n_sp = refined.len();
        for sh in refined {
            if shell_is_closed(&sh) {
                let mut c = sh;
                c.set_closed(true);
                block.change_loops_mut().push(c);
            } else {
                a_l_sh_nc.push(sh);
            }
        }
        all_taken = added.len() == a_nb_shapes;
        if all_taken {
            break;
        }
        if n_sp == 1 {
            continue;
        }
        for sh in a_l_sh_nc {
            for f in crate::bop_occt_util::iter_children(&sh) {
                if f.is_face() {
                    added.remove(&ori_key(&f));
                }
            }
        }
    }
}
