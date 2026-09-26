//! `BOPAlgo_ShellSplitter::SplitBlock` (`BOPAlgo_ShellSplitter.cxx:153`).
//!
//! Irregular connexity blocks (a TShape present FORWARD and REVERSED, or an
//! edge with other than two faces) are walked with `GetEdgeOff` / `GetFaceOff`.
//! Closed shells are kept; `RefineShell` splits a walk that still has a
//! multi-connected or same-orientation edge.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools_face::{get_face_off, CoupleOfShape};
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
fn get_edge_off_geo(the_e1: &Edge, the_f2: &TopoShape) -> Option<Edge> {
    let k = edge_ek(the_e1);
    let tr1 = crate::shell_splitter::edge_traversal(the_e1);
    for w in crate::bop_occt_util::iter_children(the_f2) {
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
            let off = match (tr1, crate::shell_splitter::edge_traversal(&e2)) {
                (Some((a1, b1)), Some((a2, b2))) => a2 == b1 && b2 == a1,
                _ => e2.0.orientation() == the_e1.0.orientation().reversed(),
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
            // used here are compared by traversal instead.
            if let (Some(a), Some(b)) = (e1, e2) {
                let same_dir = match (
                    crate::shell_splitter::edge_traversal(&a),
                    crate::shell_splitter::edge_traversal(&b),
                ) {
                    // A closed ring (both endpoints the same vertex) carries no
                    // traversal direction, so its two views would always compare
                    // equal and every such edge would become a stop edge — the
                    // port's box∪cylinder shell was split on the section circle
                    // for exactly that reason. OCCT's `RefineShell` compares only
                    // the cumulated orientations of the two occurrence views
                    // (`BOPAlgo_ShellSplitter.cxx:475-482`, `FindShape`), which is
                    // the one signal a closed ring still carries.
                    (Some(ta), Some(tb)) if ta.0 != ta.1 => ta == tb,
                    _ => a.0.orientation() == b.0.orientation(),
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
                    .find(|e| edge_ek(e) == *ek)
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
    let ctx = IntToolsContext::new();
    let mut all_taken = false;
    for a_ff in &a_lf_connected {
        if all_taken {
            break;
        }
        if !added.insert(ori_key(a_ff)) {
            continue;
        }
        let mut shell_faces: Vec<TopoShape> = vec![a_ff.clone()];
        let mut a_mefp: HashMap<crate::shell_splitter::EKey, Vec<TopoShape>> = HashMap::new();
        merge_face_into_mef(&mut a_mefp, a_ff);
        let mut i = 0usize;
        while i < shell_faces.len() {
            let a_f = shell_faces[i].clone();
            let is_boundary = a_boundary.contains(&shape_key(&a_f));
            for a_e in face_boundary_edges(&a_f) {
                if a_mefp
                    .get(&edge_ek(&a_e))
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
                let Some(a_lf) = a_ef_map.get(&edge_ek(&a_e)) else {
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
                    let Some(a_el) = get_edge_off_geo(&a_e, a_fl) else {
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
                let c = sh;
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
