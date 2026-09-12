//! Face-edge connexity used by `BOPAlgo_FillIn3DParts`.
//!
//! Source: `BOPAlgo_Tools.cxx`
//! * `MapEdgesAndFaces` at 1523
//! * `MakeConnexityBlock` at 1555
//! * `FillInternals` shell assembly at 1867-1907 (`MakeConnexityBlocks`
//!   then INTERNAL shells)
//!
//! ClassifyFaces groups candidate faces into connexity blocks that do not
//! cross the solid's own edges (`theMEAvoid`). A whole block is classified
//! by one representative face (`theFaceToClassify`): a face that touches a
//! solid-boundary or degenerated edge, when one exists, otherwise the start
//! face. This is the OCCT reduction that avoids classifying every face of a
//! sheet that sits inside a solid.

use std::collections::{HashMap, HashSet, VecDeque};

use crate::abs::{Orientation, ShapeType};
use crate::bop_occt_util::{iter_children, shape_key};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, TopoShape};
use crate::topo_tools_full::edges_of;

/// Edge TShape key → (a representative edge, faces incident on that edge).
/// Mirrors `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`.
pub type EdgeFaceMap = HashMap<usize, (TopoShape, Vec<TopoShape>)>;

/// `BOPAlgo_FillIn3DParts::MapEdgesAndFaces` (`_cxx:1523`).
///
/// Walks wires of `the_f` and, for each edge, appends `the_f` to the edge's
/// face list. Degenerated edges are still mapped (MakeConnexityBlock later
/// refuses to walk through them).
pub fn map_edges_and_faces(the_f: &TopoShape, the_ef_map: &mut EdgeFaceMap) {
    for a_w in iter_children(the_f) {
        if a_w.shape_type() != ShapeType::Wire {
            continue;
        }
        for a_e in iter_children(&a_w) {
            if a_e.shape_type() != ShapeType::Edge {
                continue;
            }
            let k = shape_key(&a_e);
            let entry = the_ef_map
                .entry(k)
                .or_insert_with(|| (a_e.clone(), Vec::new()));
            if !entry.1.iter().any(|f| f.same_tshape(the_f)) {
                entry.1.push(the_f.clone());
            }
        }
    }
}

/// Map every face in `faces` into one shared EF map (the ClassifyFaces
/// preparation when `aNbFP > 1`, `_cxx:1422-1428`).
pub fn map_edges_and_faces_all(faces: &[TopoShape]) -> EdgeFaceMap {
    let mut m = EdgeFaceMap::new();
    for f in faces {
        map_edges_and_faces(f, &mut m);
    }
    m
}

/// Result of `MakeConnexityBlock`.
#[derive(Debug, Clone)]
pub struct ConnexityFaceBlock {
    /// Faces of the block (`theLCB`).
    pub faces: Vec<TopoShape>,
    /// Preferred classification face (`theFaceToClassify`); `None` means
    /// the start face (OCCT leaves it null and the caller substitutes `aFP`).
    pub face_to_classify: Option<TopoShape>,
}

/// `BOPAlgo_FillIn3DParts::MakeConnexityBlock` (`_cxx:1555`).
///
/// Seeds `the_lcb` with `the_f_start` and grows through edges that are
/// *not* in `the_me_avoid` and *not* degenerated. Growing faces are marked
/// in `the_mf_done` so the outer classification loop skips them. A face that
/// touches an avoided or degenerated edge becomes `theFaceToClassify` (first
/// such face wins).
pub fn make_connexity_block(
    the_f_start: &TopoShape,
    the_me_avoid: &HashSet<usize>,
    the_ef_map: &EdgeFaceMap,
    the_mf_done: &mut HashSet<usize>,
) -> ConnexityFaceBlock {
    let mut the_lcb: Vec<TopoShape> = vec![the_f_start.clone()];
    let mut face_to_classify: Option<TopoShape> = None;
    if the_ef_map.is_empty() {
        return ConnexityFaceBlock {
            faces: the_lcb,
            face_to_classify,
        };
    }

    let mut i = 0;
    while i < the_lcb.len() {
        let a_f = the_lcb[i].clone();
        i += 1;
        for a_w in iter_children(&a_f) {
            if a_w.shape_type() != ShapeType::Wire {
                continue;
            }
            for a_e_s in iter_children(&a_w) {
                if a_e_s.shape_type() != ShapeType::Edge {
                    continue;
                }
                let a_e = Edge(a_e_s);
                if the_me_avoid.contains(&shape_key(&a_e.0)) || BRepTool::is_degenerated(&a_e) {
                    if face_to_classify.is_none() {
                        face_to_classify = Some(a_f.clone());
                    }
                    continue;
                }
                let Some((_, p_lf)) = the_ef_map.get(&shape_key(&a_e.0)) else {
                    continue;
                };
                for a_f_to_add in p_lf {
                    if the_mf_done.insert(shape_key(a_f_to_add)) {
                        the_lcb.push(a_f_to_add.clone());
                    }
                }
            }
        }
    }
    ConnexityFaceBlock {
        faces: the_lcb,
        face_to_classify,
    }
}

/// Representative used for `IsInternalFace`: OCCT substitutes `aFP` when
/// `aFaceToClassify` is null (`_cxx:1493-1496`).
pub fn classification_face<'a>(
    block: &'a ConnexityFaceBlock,
    start: &'a TopoShape,
) -> &'a TopoShape {
    block.face_to_classify.as_ref().unwrap_or(start)
}

/// Build INTERNAL shells from IN faces (`FillInternals` `_cxx:1867-1907`).
///
/// Each connexity block of IN faces becomes one INTERNAL shell added to
/// `the_solid`.
pub fn add_internal_face_shells(the_solid: &mut TopoShape, in_faces: &[TopoShape]) {
    if in_faces.is_empty() {
        return;
    }
    let blocks = connexity_blocks_by_edge(in_faces);
    let bld = TopoBuilder::new();
    for block in blocks {
        let mut shell = Shell::new();
        for mut a_f in block {
            a_f.set_orientation(Orientation::Internal);
            if let Some(fc) = Face::wrap(a_f) {
                bld.add_face(&mut shell, &fc);
            }
        }
        bld.add(the_solid, &shell.0);
    }
}

/// `BOPTools_AlgoTools::MakeConnexityBlocks(compound, EDGE, FACE)`.
///
/// Two faces are connected when they share a non-degenerated edge. Blocks
/// are the connected components of that graph, ordered by the first face.
pub fn connexity_blocks_by_edge(faces: &[TopoShape]) -> Vec<Vec<TopoShape>> {
    let n = faces.len();
    if n == 0 {
        return Vec::new();
    }
    let mut e2f: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, f) in faces.iter().enumerate() {
        for e in edges_of(f) {
            if BRepTool::is_degenerated(&e) {
                continue;
            }
            e2f.entry(shape_key(&e.0)).or_default().push(i);
        }
    }
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for idxs in e2f.values() {
        for &a in idxs {
            for &b in idxs {
                if a != b && !adj[a].contains(&b) {
                    adj[a].push(b);
                }
            }
        }
    }
    let mut seen = vec![false; n];
    let mut blocks = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut block = Vec::new();
        let mut q = VecDeque::new();
        q.push_back(start);
        seen[start] = true;
        while let Some(i) = q.pop_front() {
            block.push(faces[i].clone());
            for &j in &adj[i] {
                if !seen[j] {
                    seen[j] = true;
                    q.push_back(j);
                }
            }
        }
        blocks.push(block);
    }
    blocks
}

/// Edge keys of a solid used as `theMEAvoid` (the solid's own edges, which
/// connexity must not cross).
pub fn solid_edge_avoid_set(solid: &TopoShape) -> HashSet<usize> {
    edges_of(solid).into_iter().map(|e| shape_key(&e.0)).collect()
}

/// Face keys of a solid plus its own INTERNAL faces (`aMSF` in Perform).
pub fn solid_face_avoid_set(solid: &TopoShape, own_if: &[TopoShape]) -> HashSet<usize> {
    let mut s: HashSet<usize> = crate::topo_tools_full::faces_of(solid)
        .into_iter()
        .map(|f| shape_key(&f.0))
        .collect();
    for f in own_if {
        s.insert(shape_key(f));
    }
    s
}

/// Filter candidate indices whose shapes are not already faces of the solid
/// (`aIVec` construction, `_cxx:1380-1393`).
pub fn faces_to_process(indices: &[usize], v_shapes: &[TopoShape], solid_faces: &HashSet<usize>) -> Vec<usize> {
    let mut out = Vec::new();
    for &i in indices {
        if i >= v_shapes.len() {
            continue;
        }
        if !solid_faces.contains(&shape_key(&v_shapes[i])) {
            out.push(i);
        }
    }
    out.sort_unstable();
    out
}

/// Grow every unprocessed face of `candidates` into a connexity block,
/// skipping faces already in `done`.
pub fn all_connexity_blocks(
    candidates: &[TopoShape],
    me_avoid: &HashSet<usize>,
    ef_map: &EdgeFaceMap,
) -> Vec<ConnexityFaceBlock> {
    let mut done: HashSet<usize> = HashSet::new();
    let mut blocks = Vec::new();
    for start in candidates {
        if !done.insert(shape_key(start)) {
            continue;
        }
        let block = make_connexity_block(start, me_avoid, ef_map, &mut done);
        for f in &block.faces {
            done.insert(shape_key(f));
        }
        blocks.push(block);
    }
    blocks
}

/// Whether every vertex box of `block` interferes with `solid_box`
/// (`_cxx:1467-1490`). A vertex whose box is out of the solid box makes the
/// whole block OUT without running `IsInternalFace`.
pub fn block_vertices_out_of_solid(
    block: &[TopoShape],
    solid_box: &occt_core::bnd::BndBox,
) -> bool {
    if solid_box.is_whole() {
        return false;
    }
    for s in block {
        for v in crate::topo_tools_full::vertices_of(s) {
            let mut bbv = occt_core::bnd::BndBox::new();
            bbv.add_point(&BRepTool::vertex_point(&v));
            bbv.enlarge(BRepTool::vertex_tolerance(&v));
            if solid_box.is_out_box(&bbv) {
                return true;
            }
        }
    }
    false
}
