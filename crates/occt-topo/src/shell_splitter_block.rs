//! `BOPAlgo_ShellSplitter::SplitBlock` (`BOPAlgo_ShellSplitter.cxx:153`).
//!
//! Irregular connexity blocks (a TShape present FORWARD and REVERSED, or an
//! edge with other than two faces) are walked with `GetEdgeOff` / `GetFaceOff`.
//! Closed shells are kept; `RefineShell` splits a walk that still has a
//! multi-connected or same-orientation edge.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools_face::{face_off, CoupleOfShape};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::connexity_block::ConnexityBlock;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, Shell, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, is_same, wires_of_face};

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

/// `BOPTools_AlgoTools::GetEdgeOff` (`BOPTools_AlgoTools.cxx:1099-1126`).
///
/// OCCT matches `theE1` inside `theF2` with `IsSame` (the same `TopoDS_Edge`
/// TShape) and requires the opposite orientation. Its BOP guarantees that two
/// faces meeting along a boundary carry that shared edge, so `IsSame` is the
/// geometric edge. This port does not run that share step, hence the same
/// geometric identity ([`edge_ek`]) is used instead of `IsSame`, exactly as in
/// the ancestor maps above. Two geometric views of one segment may carry
/// opposite intrinsic directions, so the orientation test is expressed as the
/// traversal test it stands for in OCCT ([`crate::shell_splitter::edge_traversal`]).
fn edge_off_geo(e1: &Edge, f2: &TopoShape) -> Option<Edge> {
    let k = edge_ek(e1);
    for w in crate::bop_occt_util::iter_children(f2) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        for e in crate::bop_occt_util::iter_children(&w) {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            let e2 = Edge(e.clone());
            if edge_ek(&e2) != k {
                continue;
            }
            // `GetEdgeOff` (`BOPTools_AlgoTools.cxx:1115-1123`) requires the
            // occurrence of `theE1` inside `theF2` to carry the opposite
            // orientation; the geometric views used here are compared by
            // traversal instead ([`crate::shell_splitter::same_traversal`]).
            let off = match crate::shell_splitter::same_traversal(e1, &e2) {
                Some(same) => !same,
                None => e2.0.orientation() == e1.0.orientation().reversed(),
            };
            if off {
                return Some(e2);
            }
        }
    }
    None
}

/// Identity of an edge for the ancestor maps of `SplitBlock`.
///
/// OCCT keys `aEFMap` / `aMEFP` with `TopExp::MapShapesAndAncestors`, i.e. the
/// `TopoDS_Edge` itself (`BOPAlgo_ShellSplitter.cxx:192`, `:310`). In OCCT the
/// BOP pipeline already guarantees that two faces meeting along a boundary
/// share that `TopoDS_Edge`, so the TShape *is* the geometric edge. This port
/// does not run that share step, so the equivalent identity is the geometric
/// one the rest of the splitter uses
/// ([`crate::shell_splitter::edge_key`], the module's documented convention and
/// the key [`crate::shell_splitter::ShellSplitter::perform`] builds the blocks
/// with). Keying `SplitBlock` by TShape instead makes every face whose
/// neighbours were built by a different producer look "free-edged" and drops
/// the whole block.
fn edge_ek(e: &Edge) -> crate::shell_splitter::EKey {
    crate::shell_splitter::edge_key(e)
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

/// edge → ancestor faces, appending the face once per **edge occurrence**
/// (`TopExp::MapShapesAndAncestors`, `TopExp.cxx:80-120`: the inner explorer
/// walks `TopAbs_EDGE` with the default `CumOri = true` and calls
/// `M(index).Append(anc)` for every occurrence). Deduplicating the face per edge
/// key makes a seam stored twice in one wire (Forward and Reversed) look like a
/// *free* edge (`aLF.Extent() == 1`), so `SplitBlock`'s free-edge pass
/// (`BOPAlgo_ShellSplitter.cxx:202-214`) removes the whole lateral face and the
/// removal cascades through its neighbours.
fn map_edges_and_faces(faces: &[TopoShape]) -> HashMap<crate::shell_splitter::EKey, Vec<TopoShape>> {
    let mut mef: HashMap<crate::shell_splitter::EKey, Vec<TopoShape>> = HashMap::new();
    for f in faces {
        for e in face_boundary_edges(f) {
            let k = edge_ek(&e);
            mef.entry(k).or_default().push(f.clone());
        }
    }
    mef
}

/// `TopExp::MapShapesAndAncestors(aShell, EDGE, FACE, aMEFP)`
/// (`BOPAlgo_ShellSplitter.cxx:266`): occurrence-counting append, as above.
fn merge_face_into_mef(
    mef: &mut HashMap<crate::shell_splitter::EKey, Vec<TopoShape>>,
    face: &TopoShape,
) {
    for e in face_boundary_edges(face) {
        let k = edge_ek(&e);
        mef.entry(k).or_default().push(face.clone());
    }
}

fn shell_is_closed(shell: &TopoShape) -> bool {
    !crate::bop_build_solids::geometrically_open(shell)
}

/// `RefineShell` (`ShellSplitter.cxx:443`). Empty stop-edge set returns the
/// input shell. Otherwise the walk does not cross stop edges.
fn refine_shell(
    shell: &TopoShape,
    mef: &HashMap<crate::shell_splitter::EKey, Vec<TopoShape>>,
) -> Vec<TopoShape> {
    let faces: Vec<TopoShape> = crate::bop_occt_util::iter_children(shell)
        .into_iter()
        .filter(|c| c.is_face())
        .collect();
    if faces.is_empty() {
        return Vec::new();
    }
    let mut stop: HashSet<crate::shell_splitter::EKey> = HashSet::new();
    for (ek, lf) in mef {
        if lf.len() > 2 {
            stop.insert(*ek);
            continue;
        }
        if lf.len() == 2 {
            let e1 = face_boundary_edges(&lf[0])
                .into_iter()
                .find(|e| edge_ek(e) == *ek);
            let e2 = face_boundary_edges(&lf[1])
                .into_iter()
                .find(|e| edge_ek(e) == *ek);
            // `RefineShell` (`BOPAlgo_ShellSplitter.cxx:470-481`) stops on the
            // edges whose two faces traverse them the same way. It compares the
            // orientations of the one shared `TopoDS_Edge`; the geometric views
            // used here are compared by traversal instead
            // ([`crate::shell_splitter::same_traversal`], which also covers the
            // closed ring whose endpoints carry no traversal).
            if let (Some(a), Some(b)) = (e1, e2) {
                let same_dir = match crate::shell_splitter::same_traversal(&a, &b) {
                    Some(same) => same,
                    None => a.0.orientation() == b.0.orientation(),
                };
                if same_dir {
                    stop.insert(*ek);
                    continue;
                }
            }
        }
        let mut nb = 0usize;
        for f in lf {
            nb += 1;
            if face_boundary_edges(f).iter().any(|e| {
                edge_ek(e) == *ek && e.0.orientation() == Orientation::Internal
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
                    if stop.contains(&edge_ek(&e)) {
                        continue;
                    }
                    if e.0.orientation() == Orientation::Internal {
                        continue;
                    }
                    if BRepTool::is_degenerated(&e) {
                        continue;
                    }
                    let Some(lf) = mef.get(&edge_ek(&e)) else {
                        continue;
                    };
                    for fp1 in lf {
                        if is_same(fp1, fp) {
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
pub fn split_shell_block(block: &mut ConnexityBlock) {
    block.change_loops_mut().clear();
    let shapes_val = block.shapes().to_vec();
    let mut m_faces: HashMap<(usize, u8), TopoShape> = HashMap::new();
    for f in &shapes_val {
        m_faces.insert(ori_key(f), f.clone());
    }
    loop {
        let live: Vec<TopoShape> = m_faces.values().cloned().collect();
        let ef = map_edges_and_faces(&live);
            let begin = m_faces.len();
        for (ek, lf) in &ef {
            let Some(e_rep) = live.iter().find_map(|f| {
                face_boundary_edges(f)
                    .into_iter()
                    .find(|e| edge_ek(e) == *ek)
            }) else {
                continue;
            };
            if BRepTool::is_degenerated(&e_rep) || e_rep.0.orientation() == Orientation::Internal {
                continue;
            }
            if lf.len() == 1 {
                m_faces.remove(&ori_key(&lf[0]));
            }
        }
        let end = m_faces.len();
        if end == begin || end == 0 {
            break;
        }
    }
    if m_faces.is_empty() {
        return;
    }

    let mut lf_connected: Vec<TopoShape> = Vec::new();
    let mut boundary: HashSet<usize> = HashSet::new();
    for f in &shapes_val {
        if m_faces.contains_key(&ori_key(f)) {
            lf_connected.push(f.clone());
            if !boundary.insert(shape_key(f)) {
                boundary.remove(&shape_key(f));
            }
        }
    }
    let nb_shapes = lf_connected.len();
    let live: Vec<TopoShape> = m_faces.values().cloned().collect();
    let ef_map = map_edges_and_faces(&live);
    let mut added: HashSet<(usize, u8)> = HashSet::new();
    let ctx = IntToolsContext::new();
    let mut all_taken = false;
    for ff in &lf_connected {
        if all_taken {
            break;
        }
        if !added.insert(ori_key(ff)) {
            continue;
        }
        let mut shell_faces: Vec<TopoShape> = vec![ff.clone()];
        let mut mefp: HashMap<crate::shell_splitter::EKey, Vec<TopoShape>> = HashMap::new();
        merge_face_into_mef(&mut mefp, ff);
        let mut i = 0usize;
        while i < shell_faces.len() {
            let f_val = shell_faces[i].clone();
            let is_boundary = boundary.contains(&shape_key(&f_val));
            for e_val in face_boundary_edges(&f_val) {
                if mefp
                    .get(&edge_ek(&e_val))
                    .map(|l| l.len() > 1)
                    .unwrap_or(false)
                {
                    continue;
                }
                if e_val.0.orientation() == Orientation::Internal {
                    continue;
                }
                if BRepTool::is_degenerated(&e_val) {
                    continue;
                }
                let Some(lf_val) = ef_map.get(&edge_ek(&e_val)) else {
                    continue;
                };
                if lf_val.is_empty() {
                    continue;
                }
                let mut lcs_off: Vec<CoupleOfShape> = Vec::new();
                let mut nb_ways_inside = 0i32;
                let mut sel_val: Option<TopoShape> = None;
                for fl in lf_val {
                    if is_same(&f_val, fl) || added.contains(&ori_key(fl)) {
                        continue;
                    }
                    let Some(el) = edge_off_geo(&e_val, fl) else {
                        continue;
                    };
                    if is_boundary && !boundary.contains(&shape_key(fl)) {
                        nb_ways_inside += 1;
                        sel_val = Some(fl.clone());
                    }
                    lcs_off.push(CoupleOfShape::new(el.0, fl.clone()));
                }
                if lcs_off.is_empty() {
                    continue;
                }
                if !is_boundary || nb_ways_inside != 1 {
                    if lcs_off.len() == 1 {
                        sel_val = Some(lcs_off[0].shape2.clone());
                    } else if let Some((f, _)) =
                        face_off(&e_val, &Face(f_val.clone()), &lcs_off, &ctx)
                    {
                        sel_val = Some(f.0);
                    }
                }
                if let Some(sel) = sel_val {
                    if added.insert(ori_key(&sel)) {
                        shell_faces.push(sel.clone());
                        merge_face_into_mef(&mut mefp, &sel);
                    }
                }
            }
            i += 1;
        }
        let walked = make_shell(&shell_faces);
        let mut l_sh_nc: Vec<TopoShape> = Vec::new();
        let refined = refine_shell(&walked, &mefp);
        let n_sp = refined.len();

        for sh in refined {
            if shell_is_closed(&sh) {
                let c = sh;
                c.set_closed(true);
                block.change_loops_mut().push(c);
            } else {
                l_sh_nc.push(sh);
            }
        }
        all_taken = added.len() == nb_shapes;
        if all_taken {
            break;
        }
        if n_sp == 1 {
            continue;
        }
        for sh in l_sh_nc {
            for f in crate::bop_occt_util::iter_children(&sh) {
                if f.is_face() {
                    added.remove(&ori_key(&f));
                }
            }
        }
    }
}
