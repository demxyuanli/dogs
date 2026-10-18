//! `BOPAlgo_Tools::ClassifyFaces` / `FillInternals` and `FillIn3DParts`.
//!
//! Source: `BOPAlgo_Tools.cxx` (`BOPAlgo_ShapeBox` at 1215,
//! `BOPAlgo_FillIn3DParts` at 1245, `Perform` at 1334,
//! `MapEdgesAndFaces` at 1523, `MakeConnexityBlock` at 1555,
//! `ClassifyFaces` at 1622, `FillInternals` at 1751).
//! Parallel `BOPTools_Parallel::Perform` is sequential. `IsInternalFace` is
//! `AlgoTools::compute_state` of a representative face point against the solid.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::int_tools_full::IntToolsContext;
use crate::iterator::ShapeIterator;
use crate::bop_occt_util::shape_key;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools_full::{edges_of, map_shapes_and_ancestors, vertices_of, ShapeMapKey};

/// `BOPAlgo_ShapeBox`.
#[derive(Clone)]
pub struct ShapeBox {
    pub shape: TopoShape,
    pub box_: BndBox,
}

impl ShapeBox {
    pub fn new(shape: TopoShape, box_: BndBox) -> Self {
        Self { shape, box_ }
    }
}

fn map_edges_and_faces(f: &TopoShape, ef_map: &mut HashMap<ShapeMapKey, (TopoShape, Vec<TopoShape>)>) {
    for w in ShapeIterator::of_shape(f) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        for e in ShapeIterator::of_shape(&w) {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            let k = shape_key(&e);
            let entry = ef_map
                .entry(k)
                .or_insert_with(|| (e.clone(), Vec::new()));
            // `BOPAlgo_Tools.cxx:1544-1549`: `pLF->Append(theF)` with no
            // deduplication. `IsInternalFace` reads `aLF.Extent() == 1` / `== 2`
            // (`BOPTools_AlgoTools.cxx:849-869`), so a `same_tshape` fence here
            // would collapse repeated ancestors and pick the wrong arm.
            entry.1.push(f.clone());
        }
    }
}

fn make_connexity_block(
    f_start: &TopoShape,
    me_avoid: &HashSet<ShapeMapKey>,
    ef_map: &HashMap<ShapeMapKey, (TopoShape, Vec<TopoShape>)>,
    mf_done: &mut HashSet<ShapeMapKey>,
    lcb: &mut Vec<TopoShape>,
    face_to_classify: &mut Option<TopoShape>,
) {
    lcb.push(f_start.clone());
    if ef_map.is_empty() {
        return;
    }
    let mut i = 0;
    while i < lcb.len() {
        let a_f = lcb[i].clone();
        i += 1;
        for w in ShapeIterator::of_shape(&a_f) {
            if w.shape_type() != ShapeType::Wire {
                continue;
            }
            for e in ShapeIterator::of_shape(&w) {
                if e.shape_type() != ShapeType::Edge {
                    continue;
                }
                let k = shape_key(&e);
                if me_avoid.contains(&k) || BRepTool::is_degenerated(&Edge(e.clone())) {
                    if face_to_classify.is_none() {
                        *face_to_classify = Some(a_f.clone());
                    }
                    continue;
                }
                let Some((_, p_lf)) = ef_map.get(&k) else {
                    continue;
                };
                for f_to_add in p_lf {
                    if mf_done.insert(shape_key(f_to_add)) {
                        lcb.push(f_to_add.clone());
                    }
                }
            }
        }
    }
}

/// `BOPTools_AlgoTools::IsInternalFace` via a representative point.
fn is_internal_face(
    face: &Face,
    solid: &TopoShape,
    mefds: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    tol: f64,
    ctx: &IntToolsContext,
) -> bool {
    crate::algo_tools_face::is_internal_face(
        face,
        solid,
        mefds,
        tol,
        ctx,
    )
}

/// `BOPTools_AlgoTools::ComputeStateByOnePoint` (`BOPAlgo_Tools.cxx` uses
/// this from `FillInternalShapes`). Public so the Builder-level
/// `FillInternalShapes` translation can call it.
pub fn compute_state_by_one_point(
    shape: &TopoShape,
    solid: &TopoShape,
    tol: f64,
) -> FaceState {
    match shape.shape_type() {
        ShapeType::Vertex => {
            let p = BRepTool::vertex_point(&Vertex(shape.clone()));
            AlgoTools::compute_state(solid, &p, tol).unwrap_or(FaceState::Unknown)
        }
        ShapeType::Edge => {
            let e = Edge(shape.clone());
            let (a, b) = BRepTool::edge_parameters(&e);
            let Some(curve) = BRepTool::edge_curve(&e) else {
                return FaceState::Unknown;
            };
            if !a.is_finite() || !b.is_finite() {
                return FaceState::Unknown;
            }
            let p = curve.d0(0.5 * (a + b));
            AlgoTools::compute_state(solid, &p, tol).unwrap_or(FaceState::Unknown)
        }
        ShapeType::Face => {
            // `BOPTools_AlgoTools::ComputeState` (face, solid, bounds, context)
            // (`BOPTools_AlgoTools.cxx:666-715`).
            let f = Face(shape.clone());
            crate::algo_tools_face::compute_state_face(&f, solid, tol)
        }
        _ => {
            let kids = shape.tshape.read().unwrap().children.clone();
            for sub in kids {
                let st = compute_state_by_one_point(&sub, solid, tol);
                if st != FaceState::Unknown {
                    return st;
                }
            }
            FaceState::Unknown
        }
    }
}

/// `BOPAlgo_FillIn3DParts::Perform` (`BOPAlgo_Tools.cxx:1334`).
pub fn fill_in_3d_parts(
    solid: &TopoShape,
    box_s: &BndBox,
    own_if: &[TopoShape],
    v_shape_box: &[ShapeBox],
    ctx: &IntToolsContext,
) -> Vec<TopoShape> {
    let mut in_faces: Vec<TopoShape> = Vec::new();
    let mut lifp: Vec<usize> = Vec::new();
    for (i, sb) in v_shape_box.iter().enumerate() {
        if !box_s.is_out_box(&sb.box_) {
            lifp.push(i);
        }
    }
    if lifp.is_empty() {
        return in_faces;
    }
    let mut mse: HashSet<ShapeMapKey> = HashSet::new();
    let mut msf: HashSet<ShapeMapKey> = HashSet::new();
    for e in edges_of(solid) {
        mse.insert(shape_key(&e.0));
    }
    for f in ShapeIterator::of_shape(solid) {
        if f.shape_type() == ShapeType::Face {
            msf.insert(shape_key(&f));
        }
    }
    let b_is_empty = msf.is_empty();
    for s in own_if {
        msf.insert(shape_key(s));
    }
    let mut ivec: Vec<usize> = Vec::new();
    for n_fp in lifp {
        let fp = &v_shape_box[n_fp].shape;
        if !msf.contains(&shape_key(fp)) {
            ivec.push(n_fp);
        }
    }
    ivec.sort_unstable();
    if b_is_empty {
        for k in ivec {
            in_faces.push(v_shape_box[k].shape.clone());
        }
        return in_faces;
    }
    let mut mefp: HashMap<ShapeMapKey, (TopoShape, Vec<TopoShape>)> = HashMap::new();
    if ivec.len() > 1 {
        for &k in &ivec {
            map_edges_and_faces(&v_shape_box[k].shape, &mut mefp);
        }
    }
    // `_cxx:1498-1501`: `aMEFDS` is filled by
    // `TopExp::MapShapesAndAncestors(mySolid, EDGE, FACE)`, keyed by
    // `TopTools_ShapeMapHasher` and built once.
    let mut a_mefds: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let mut a_mf_done: HashSet<ShapeMapKey> = HashSet::new();
    for &n_fp in &ivec {
        let fp = v_shape_box[n_fp].shape.clone();
        if !a_mf_done.insert(shape_key(&fp)) {
            continue;
        }
        let mut lcbf: Vec<TopoShape> = Vec::new();
        let mut a_face_to_classify: Option<TopoShape> = None;
        make_connexity_block(
            &fp,
            &mse,
            &mefp,
            &mut a_mf_done,
            &mut lcbf,
            &mut a_face_to_classify,
        );
        if !box_s.is_whole() {
            let mut b_out = false;
            for s in &lcbf {
                for v in vertices_of(s) {
                    let mut bbv = BndBox::new();
                    bbv.add_point(&BRepTool::vertex_point(&v));
                    bbv.enlarge(BRepTool::vertex_tolerance(&v));
                    if box_s.is_out_box(&bbv) {
                        b_out = true;
                        break;
                    }
                }
                if b_out {
                    break;
                }
            }
            if b_out {
                continue;
            }
        }
        let a_face = a_face_to_classify.unwrap_or(fp);
        if a_mefds.is_empty() {
            a_mefds = map_shapes_and_ancestors(solid, ShapeType::Edge, ShapeType::Face);
        }
        if is_internal_face(&Face(a_face), solid, &a_mefds, CONFUSION, ctx) {
            in_faces.extend(lcbf);
        }
    }
    in_faces
}

/// `BOPAlgo_Tools::ClassifyFaces` (`BOPAlgo_Tools.cxx:1622`).
pub fn classify_faces(
    the_faces: &[TopoShape],
    solids: &[TopoShape],
    ctx: &IntToolsContext,
    shape_box_map: &HashMap<ShapeMapKey, BndBox>,
    solids_if: &HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> HashMap<ShapeMapKey, Vec<TopoShape>> {
    let mut vsb: Vec<ShapeBox> = Vec::new();
    for a_f in the_faces {
        let bnd_box = shape_box_map
            .get(&shape_key(a_f))
            .cloned()
            .unwrap_or_else(|| shape_bbox(a_f));
        vsb.push(ShapeBox::new(a_f.clone(), bnd_box));
    }
    let mut in_parts: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    for a_solid in solids {
        let bnd_box = shape_box_map
            .get(&shape_key(a_solid))
            .cloned()
            .unwrap_or_else(|| {
                let mut b = shape_bbox(a_solid);
                if !b.is_whole() && AlgoTools::is_inverted_solid(a_solid) {
                    b.set_whole();
                }
                b
            });
        let own = solids_if
            .get(&shape_key(a_solid))
            .cloned()
            .unwrap_or_default();
        let in_faces = fill_in_3d_parts(a_solid, &bnd_box, &own, &vsb, ctx);
        in_parts.insert(shape_key(a_solid), in_faces);
    }
    in_parts
}

/// `BOPAlgo_Tools::FillInternals` (`BOPAlgo_Tools.cxx:1751`).
pub fn fill_internals(
    solids: &mut [TopoShape],
    parts: &[TopoShape],
    images: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    _the_ctx: &IntToolsContext,
) {
    if solids.is_empty() || parts.is_empty() {
        return;
    }
    let mut ms_solids: HashSet<ShapeMapKey> = HashSet::new();
    for a_solid in solids.iter() {
        if a_solid.shape_type() != ShapeType::Solid {
            continue;
        }
        for v in vertices_of(a_solid) {
            ms_solids.insert(shape_key(&v.0));
        }
        for e in edges_of(a_solid) {
            ms_solids.insert(shape_key(&e.0));
        }
        for f in ShapeIterator::of_shape(a_solid) {
            if f.shape_type() == ShapeType::Face {
                ms_solids.insert(shape_key(&f));
            }
        }
    }
    let mut l_parts: Vec<TopoShape> = Vec::new();
    let mut l_input: Vec<TopoShape> = parts.to_vec();
    let mut ii = 0;
    while ii < l_input.len() {
        let part = l_input[ii].clone();
        ii += 1;
        match part.shape_type() {
            ShapeType::Vertex | ShapeType::Edge | ShapeType::Face => {
                if let Some(p_im) = images.get(&shape_key(&part)) {
                    for part_im in p_im {
                        if !ms_solids.contains(&shape_key(part_im)) {
                            l_parts.push(part_im.clone());
                        }
                    }
                } else if !ms_solids.contains(&shape_key(&part)) {
                    l_parts.push(part);
                }
            }
            _ => {
                for c in ShapeIterator::of_shape(&part) {
                    l_input.push(c);
                }
            }
        }
    }
    let mut an_in_faces: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let b = TopoBuilder::new();
    for sd in solids.iter_mut() {
        if sd.shape_type() != ShapeType::Solid {
            continue;
        }
        let mut keep: Vec<TopoShape> = Vec::new();
        for part in l_parts.drain(..) {
            let state = compute_state_by_one_point(&part, sd, CONFUSION);
            if state == FaceState::In {
                if part.shape_type() == ShapeType::Face {
                    an_in_faces
                        .entry(shape_key(sd))
                        .or_default()
                        .push(part);
                } else {
                    let mut p = part;
                    p.set_orientation(Orientation::Internal);
                    b.add(sd, &p);
                }
            } else {
                keep.push(part);
            }
        }
        l_parts = keep;
    }
    for (k, faces) in an_in_faces {
        let Some(sd) = solids.iter_mut().find(|s| shape_key(s) == k) else {
            continue;
        };
        let blocks = connexity_face_blocks(&faces);
        for block in blocks {
        let faces: Vec<crate::shape::Face> = block
            .into_iter()
            .map(|mut f| {
                f.set_orientation(Orientation::Internal);
                crate::shape::Face(f)
            })
            .collect();
        let sh = b.make_shell(&faces);
        b.add(sd, &sh.0);
        }
    }
}

fn connexity_face_blocks(faces: &[TopoShape]) -> Vec<Vec<TopoShape>> {
    let mut e2f: HashMap<ShapeMapKey, Vec<usize>> = HashMap::new();
    for (i, f) in faces.iter().enumerate() {
        for e in edges_of(f) {
            e2f.entry(shape_key(&e.0)).or_default().push(i);
        }
    }
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for idxs in e2f.values() {
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
    for i in 0..faces.len() {
        if !seen.insert(i) {
            continue;
        }
        let mut stack = vec![i];
        let mut block = Vec::new();
        while let Some(n) = stack.pop() {
            block.push(faces[n].clone());
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
