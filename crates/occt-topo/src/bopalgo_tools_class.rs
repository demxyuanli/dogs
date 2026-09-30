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
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of, vertices_of};

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

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

fn map_edges_and_faces(the_f: &TopoShape, the_ef_map: &mut HashMap<usize, (TopoShape, Vec<TopoShape>)>) {
    for w in ShapeIterator::of_shape(the_f) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        for e in ShapeIterator::of_shape(&w) {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            let k = shape_key(&e);
            let entry = the_ef_map
                .entry(k)
                .or_insert_with(|| (e.clone(), Vec::new()));
            if !entry.1.iter().any(|f| f.same_tshape(the_f)) {
                entry.1.push(the_f.clone());
            }
        }
    }
}

fn make_connexity_block(
    the_f_start: &TopoShape,
    the_me_avoid: &HashSet<usize>,
    the_ef_map: &HashMap<usize, (TopoShape, Vec<TopoShape>)>,
    the_mf_done: &mut HashSet<usize>,
    the_lcb: &mut Vec<TopoShape>,
    the_face_to_classify: &mut Option<TopoShape>,
) {
    the_lcb.push(the_f_start.clone());
    if the_ef_map.is_empty() {
        return;
    }
    let mut i = 0;
    while i < the_lcb.len() {
        let a_f = the_lcb[i].clone();
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
                if the_me_avoid.contains(&k) || BRepTool::is_degenerated(&Edge(e.clone())) {
                    if the_face_to_classify.is_none() {
                        *the_face_to_classify = Some(a_f.clone());
                    }
                    continue;
                }
                let Some((_, p_lf)) = the_ef_map.get(&k) else {
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
}

fn face_sample_point(face: &Face) -> Option<occt_core::gp::GpPnt> {
    for e in edges_of(&face.0) {
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        let (a, b) = BRepTool::edge_parameters(&e);
        if a.is_finite() && b.is_finite() {
            if let Some(c) = BRepTool::edge_curve(&e) {
                return Some(c.d0(0.5 * (a + b)));
            }
        }
    }
    let surf = BRepTool::face_surface(face)?;
    let (u1, u2, v1, v2) = crate::int_tools_full::IntToolsContext::new().uv_bounds(face);
    if u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite() {
        return Some(surf.d0(0.5 * (u1 + u2), 0.5 * (v1 + v2)));
    }
    None
}

/// `BOPTools_AlgoTools::IsInternalFace` via a representative point.
fn is_internal_face(
    the_face: &Face,
    the_solid: &TopoShape,
    the_mefds: &HashMap<usize, (TopoShape, Vec<TopoShape>)>,
    the_tol: f64,
    the_ctx: &IntToolsContext,
) -> bool {
    crate::algo_tools_face::is_internal_face(
        the_face,
        the_solid,
        the_mefds,
        the_tol,
        the_ctx,
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
            let f = Face(shape.clone());
            let Some(p) = face_sample_point(&f) else {
                return FaceState::Unknown;
            };
            AlgoTools::compute_state(solid, &p, tol).unwrap_or(FaceState::Unknown)
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
    the_solid: &TopoShape,
    the_box_s: &BndBox,
    the_own_if: &[TopoShape],
    the_v_shape_box: &[ShapeBox],
    the_ctx: &IntToolsContext,
) -> Vec<TopoShape> {
    let mut my_in_faces: Vec<TopoShape> = Vec::new();
    let mut a_lifp: Vec<usize> = Vec::new();
    for (i, sb) in the_v_shape_box.iter().enumerate() {
        if !the_box_s.is_out_box(&sb.box_) {
            a_lifp.push(i);
        }
    }
    if a_lifp.is_empty() {
        return my_in_faces;
    }
    let mut a_mse: HashSet<usize> = HashSet::new();
    let mut a_msf: HashSet<usize> = HashSet::new();
    for e in edges_of(the_solid) {
        a_mse.insert(shape_key(&e.0));
    }
    for f in ShapeIterator::of_shape(the_solid) {
        if f.shape_type() == ShapeType::Face {
            a_msf.insert(shape_key(&f));
        }
    }
    let b_is_empty = a_msf.is_empty();
    for s in the_own_if {
        a_msf.insert(shape_key(s));
    }
    let mut a_ivec: Vec<usize> = Vec::new();
    for n_fp in a_lifp {
        let a_fp = &the_v_shape_box[n_fp].shape;
        if !a_msf.contains(&shape_key(a_fp)) {
            a_ivec.push(n_fp);
        }
    }
    a_ivec.sort_unstable();
    if b_is_empty {
        for k in a_ivec {
            my_in_faces.push(the_v_shape_box[k].shape.clone());
        }
        return my_in_faces;
    }
    let mut a_mefp: HashMap<usize, (TopoShape, Vec<TopoShape>)> = HashMap::new();
    if a_ivec.len() > 1 {
        for &k in &a_ivec {
            map_edges_and_faces(&the_v_shape_box[k].shape, &mut a_mefp);
        }
    }
    let mut a_mefds: HashMap<usize, (TopoShape, Vec<TopoShape>)> = HashMap::new();
    let mut a_mf_done: HashSet<usize> = HashSet::new();
    for &n_fp in &a_ivec {
        let a_fp = the_v_shape_box[n_fp].shape.clone();
        if !a_mf_done.insert(shape_key(&a_fp)) {
            continue;
        }
        let mut a_lcbf: Vec<TopoShape> = Vec::new();
        let mut a_face_to_classify: Option<TopoShape> = None;
        make_connexity_block(
            &a_fp,
            &a_mse,
            &a_mefp,
            &mut a_mf_done,
            &mut a_lcbf,
            &mut a_face_to_classify,
        );
        if !the_box_s.is_whole() {
            let mut b_out = false;
            for s in &a_lcbf {
                for v in vertices_of(s) {
                    let mut a_bbv = BndBox::new();
                    a_bbv.add_point(&BRepTool::vertex_point(&v));
                    a_bbv.enlarge(BRepTool::vertex_tolerance(&v));
                    if the_box_s.is_out_box(&a_bbv) {
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
        let a_face = a_face_to_classify.unwrap_or(a_fp);
        if a_mefds.is_empty() {
            for e in edges_of(the_solid) {
                let k = shape_key(&e.0);
                let faces: Vec<TopoShape> = ShapeIterator::of_shape(the_solid)
                    .filter(|f| {
                        f.shape_type() == ShapeType::Face
                            && edges_of(f).iter().any(|fe| fe.0.same_tshape(&e.0))
                    })
                    .collect();
                a_mefds.insert(k, (e.0, faces));
            }
        }
        if is_internal_face(&Face(a_face), the_solid, &a_mefds, CONFUSION, the_ctx) {
            my_in_faces.extend(a_lcbf);
        }
    }
    my_in_faces
}

/// `BOPAlgo_Tools::ClassifyFaces` (`BOPAlgo_Tools.cxx:1622`).
pub fn classify_faces(
    the_faces: &[TopoShape],
    the_solids: &[TopoShape],
    the_ctx: &IntToolsContext,
    the_shape_box_map: &HashMap<usize, BndBox>,
    the_solids_if: &HashMap<usize, Vec<TopoShape>>,
) -> HashMap<usize, Vec<TopoShape>> {
    let mut a_vsb: Vec<ShapeBox> = Vec::new();
    for a_f in the_faces {
        let a_box = the_shape_box_map
            .get(&shape_key(a_f))
            .cloned()
            .unwrap_or_else(|| shape_bbox(a_f));
        a_vsb.push(ShapeBox::new(a_f.clone(), a_box));
    }
    let mut the_in_parts: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for a_solid in the_solids {
        let a_box = the_shape_box_map
            .get(&shape_key(a_solid))
            .cloned()
            .unwrap_or_else(|| {
                let mut b = shape_bbox(a_solid);
                if !b.is_whole() && AlgoTools::is_inverted_solid(a_solid) {
                    b.set_whole();
                }
                b
            });
        let own = the_solids_if
            .get(&shape_key(a_solid))
            .cloned()
            .unwrap_or_default();
        let in_faces = fill_in_3d_parts(a_solid, &a_box, &own, &a_vsb, the_ctx);
        the_in_parts.insert(shape_key(a_solid), in_faces);
    }
    the_in_parts
}

/// `BOPAlgo_Tools::FillInternals` (`BOPAlgo_Tools.cxx:1751`).
pub fn fill_internals(
    the_solids: &mut [TopoShape],
    the_parts: &[TopoShape],
    the_images: &HashMap<usize, Vec<TopoShape>>,
    _the_ctx: &IntToolsContext,
) {
    if the_solids.is_empty() || the_parts.is_empty() {
        return;
    }
    let mut a_ms_solids: HashSet<usize> = HashSet::new();
    for a_solid in the_solids.iter() {
        if a_solid.shape_type() != ShapeType::Solid {
            continue;
        }
        for v in vertices_of(a_solid) {
            a_ms_solids.insert(shape_key(&v.0));
        }
        for e in edges_of(a_solid) {
            a_ms_solids.insert(shape_key(&e.0));
        }
        for f in ShapeIterator::of_shape(a_solid) {
            if f.shape_type() == ShapeType::Face {
                a_ms_solids.insert(shape_key(&f));
            }
        }
    }
    let mut a_l_parts: Vec<TopoShape> = Vec::new();
    let mut a_l_input: Vec<TopoShape> = the_parts.to_vec();
    let mut ii = 0;
    while ii < a_l_input.len() {
        let a_part = a_l_input[ii].clone();
        ii += 1;
        match a_part.shape_type() {
            ShapeType::Vertex | ShapeType::Edge | ShapeType::Face => {
                if let Some(p_im) = the_images.get(&shape_key(&a_part)) {
                    for a_part_im in p_im {
                        if !a_ms_solids.contains(&shape_key(a_part_im)) {
                            a_l_parts.push(a_part_im.clone());
                        }
                    }
                } else if !a_ms_solids.contains(&shape_key(&a_part)) {
                    a_l_parts.push(a_part);
                }
            }
            _ => {
                for c in ShapeIterator::of_shape(&a_part) {
                    a_l_input.push(c);
                }
            }
        }
    }
    let mut an_in_faces: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    let b = TopoBuilder::new();
    for a_sd in the_solids.iter_mut() {
        if a_sd.shape_type() != ShapeType::Solid {
            continue;
        }
        let mut keep: Vec<TopoShape> = Vec::new();
        for a_part in a_l_parts.drain(..) {
            let a_state = compute_state_by_one_point(&a_part, a_sd, CONFUSION);
            if a_state == FaceState::In {
                if a_part.shape_type() == ShapeType::Face {
                    an_in_faces
                        .entry(shape_key(a_sd))
                        .or_default()
                        .push(a_part);
                } else {
                    let mut p = a_part;
                    p.set_orientation(Orientation::Internal);
                    b.add(a_sd, &p);
                }
            } else {
                keep.push(a_part);
            }
        }
        a_l_parts = keep;
    }
    for (k, faces) in an_in_faces {
        let Some(a_sd) = the_solids.iter_mut().find(|s| shape_key(s) == k) else {
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
        b.add(a_sd, &sh.0);
        }
    }
}

fn connexity_face_blocks(faces: &[TopoShape]) -> Vec<Vec<TopoShape>> {
    let mut e2f: HashMap<usize, Vec<usize>> = HashMap::new();
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
