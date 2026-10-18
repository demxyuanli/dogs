//! `BOPAlgo_Tools::ClassifyFaces` with AABB culling in place of BVH.
//!
//! Source: `BOPAlgo_Tools.cxx:1622-1747` (`ClassifyFaces`) and
//! `BOPAlgo_FillIn3DParts::Perform` at 1334. The OCCT implementation builds
//! a `BOPTools_BoxTree` over face boxes and runs one `FillIn3DParts` task
//! per solid in parallel. This port:
//! * selects candidate face boxes with `Bnd_Box::IsOut` (AABB, not BVH);
//! * runs the per-solid classification sequentially
//!   (`BOPTools_Parallel::Perform` with `theRunParallel` ignored);
//! * reuses the connexity-block / `IsInternalFace` path already translated
//!   in [`crate::bopalgo_tools_class`].
//!
//! The Builder-level `FillIn3DParts` (`_3.cxx:97`) lives in
//! [`crate::bop_fill_in3d`] and calls this module for the classification step.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::HashMap;

use occt_core::bnd::BndBox;

use crate::algo_tools::AlgoTools;
use crate::bop_occt_util::shape_key;
use crate::bopalgo_tools_class::{fill_in_3d_parts, ShapeBox};
use crate::bbox_from_geometry::shape_bbox;
use crate::int_tools_full::IntToolsContext;
use crate::shape::TopoShape;

use crate::topo_tools::ShapeMapKey;

/// A face plus its bounding box, the `BOPAlgo_ShapeBox` vector entry.
pub type FaceBox = ShapeBox;

/// Build the `BOPAlgo_VectorOfShapeBox` from `the_faces` and an optional
/// precomputed box map (`theShapeBoxMap`). Missing entries are computed
/// with `BRepBndLib::Add` (`_cxx:1660-1666`).
pub fn make_shape_box_vector(
    faces: &[TopoShape],
    shape_box_map: &HashMap<ShapeMapKey, BndBox>,
) -> Vec<ShapeBox> {
    let mut vsb = Vec::with_capacity(faces.len());
    for f in faces {
        let bnd_box = shape_box_map
            .get(&shape_key(f))
            .copied()
            .filter(|b| !b.is_void())
            .unwrap_or_else(|| shape_bbox(f));
        vsb.push(ShapeBox::new(f.clone(), bnd_box));
    }
    vsb
}

/// AABB tree stand-in: indices of `v_sb` whose boxes are not out of
/// `box_s` (`BOPTools_BoxTreeSelector::Select`, `_cxx:1345-1354`).
pub fn aabb_select(box_s: &BndBox, v_sb: &[ShapeBox]) -> Vec<usize> {
    let mut idx = Vec::new();
    for (i, sb) in v_sb.iter().enumerate() {
        if !box_s.is_out_box(&sb.box_) {
            idx.push(i);
        }
    }
    idx
}

/// Solid box used by ClassifyFaces (`_cxx:1693-1712`): prefer the map as-is.
/// Invert-to-whole only when the box is computed here (`BRepBndLib::Add`).
pub fn solid_box_for_classify(
    solid: &TopoShape,
    shape_box_map: &HashMap<ShapeMapKey, BndBox>,
) -> BndBox {
    if let Some(b) = shape_box_map.get(&shape_key(solid)).copied() {
        if !b.is_void() {
            return b;
        }
    }
    let mut bnd_box = shape_bbox(solid);
    if !bnd_box.is_whole() && AlgoTools::is_inverted_solid(solid) {
        bnd_box.set_whole();
    }
    bnd_box
}

/// One `BOPAlgo_FillIn3DParts` task (`_cxx:1688-1722`).
#[derive(Clone)]
pub struct FillIn3dTask {
    pub solid: TopoShape,
    pub box_s: BndBox,
    pub own_if: Vec<TopoShape>,
}

impl FillIn3dTask {
    pub fn new(solid: TopoShape, box_s: BndBox, own_if: Vec<TopoShape>) -> Self {
        Self {
            solid,
            box_s,
            own_if,
        }
    }

    /// `BOPAlgo_FillIn3DParts::Perform`.
    pub fn perform(&self, v_sb: &[ShapeBox], ctx: &IntToolsContext) -> Vec<TopoShape> {
        fill_in_3d_parts(&self.solid, &self.box_s, &self.own_if, v_sb, ctx)
    }
}

/// `BOPAlgo_Tools::ClassifyFaces` (`_cxx:1622`).
///
/// Returns `theInParts`: draft-solid TShape key → faces classified IN.
/// The OCCT map is keyed by the draft solid shape itself; this port uses
/// [`shape_key`] of that solid so callers can look up by the same identity
/// FillIn3DParts stored in `aDraftSolid`.
pub fn classify_faces_occt(
    faces: &[TopoShape],
    solids: &[TopoShape],
    ctx: &IntToolsContext,
    shape_box_map: &HashMap<ShapeMapKey, BndBox>,
    solids_if: &HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> HashMap<ShapeMapKey, Vec<TopoShape>> {
    classify_faces_aabb(
        faces,
        solids,
        ctx,
        shape_box_map,
        solids_if,
    )
}

/// ClassifyFaces that also returns the per-solid tasks (diagnostics).
pub fn classify_faces_with_tasks(
    faces: &[TopoShape],
    solids: &[TopoShape],
    ctx: &IntToolsContext,
    shape_box_map: &HashMap<ShapeMapKey, BndBox>,
    solids_if: &HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> (HashMap<ShapeMapKey, Vec<TopoShape>>, Vec<FillIn3dTask>) {
    let v_sb = make_shape_box_vector(faces, shape_box_map);
    let mut tasks = Vec::new();
    let mut in_parts: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    for solid in solids {
        let box_s = solid_box_for_classify(solid, shape_box_map);
        let own = solids_if
            .get(&shape_key(solid))
            .cloned()
            .unwrap_or_default();
        let task = FillIn3dTask::new(solid.clone(), box_s, own);
        let in_faces = task.perform(&v_sb, ctx);
        in_parts.insert(shape_key(solid), in_faces);
        tasks.push(task);
    }
    let _ = v_sb;
    (in_parts, tasks)
}

/// Select-and-classify a single solid (the sequential body of one
/// `FillIn3DParts` parallel task).
pub fn classify_one_solid(
    solid: &TopoShape,
    faces: &[TopoShape],
    own_if: &[TopoShape],
    shape_box_map: &HashMap<ShapeMapKey, BndBox>,
    ctx: &IntToolsContext,
) -> Vec<TopoShape> {
    let v_sb = make_shape_box_vector(faces, shape_box_map);
    let box_s = solid_box_for_classify(solid, shape_box_map);
    fill_in_3d_parts(solid, &box_s, own_if, &v_sb, ctx)
}

/// Merge `the_in_faces` with `internal` as FillIn3DParts does when
/// writing `myInParts` (`_3.cxx:243-260`).
pub fn combine_in_and_internal(
    in_faces: &[TopoShape],
    internal: &[TopoShape],
) -> Vec<TopoShape> {
    let mut out = Vec::with_capacity(in_faces.len() + internal.len());
    out.extend(in_faces.iter().cloned());
    out.extend(internal.iter().cloned());
    out
}

/// True when a solid needs a draft binding after classification
/// (`_3.cxx:225-241`): either it has IN faces, or one of its shells has an
/// image (the solid was split even without IN faces).
pub fn solid_needs_draft(
    history: &crate::bop_hist::BopHistory,
    source_solid: &TopoShape,
    nb_in: usize,
) -> bool {
    if nb_in != 0 {
        return true;
    }
    crate::bop_occt_util::iter_children(source_solid)
        .into_iter()
        .any(|sh| history.has_image(&sh))
}

/// `BOPAlgo_FillIn3DParts::Perform` using the AABB tree (`_cxx:1334`).
///
/// Steps match OCCT:
/// 1. Select face boxes not out of `box_s` (`BoxTreeSelector`);
/// 2. Map solid edges/faces (`aMSE`, `aMSF`) and add own INTERNAL faces;
/// 3. Filter selected faces that are not already of the solid (`aIVec`);
/// 4. Empty solid → every remaining face is IN;
/// 5. Connexity blocks that do not cross solid edges;
/// 6. Vertex-box OUT cull;
/// 7. `IsInternalFace` on the representative, then the whole block is IN.
pub fn perform_fill_in_3d_parts_tree(
    solid: &TopoShape,
    box_s: &BndBox,
    own_if: &[TopoShape],
    v_sb: &[ShapeBox],
    ctx: &IntToolsContext,
) -> Vec<TopoShape> {
    use crate::algo_tools_face::is_internal_face;
    use crate::bop_aabb_faces::AabbTree;
    use crate::bop_connexity_faces::{
        block_vertices_out_of_solid, classification_face, make_connexity_block,
        map_edges_and_faces, solid_edge_avoid_set, solid_face_avoid_set, EdgeFaceMap,
    };
    use crate::bop_occt_util::origin_map_shapes_and_ancestors;
    use crate::shape::Face;
    use occt_core::precision::CONFUSION;
    use std::collections::HashSet;

    let mut in_faces: Vec<TopoShape> = Vec::new();
    let tree = AabbTree::from_shape_boxes(v_sb);
    let lifp = if tree.is_empty() {
        Vec::new()
    } else {
        tree.select(box_s)
    };
    if lifp.is_empty() {
        return in_faces;
    }

    let mse = solid_edge_avoid_set(solid);
    let msf = solid_face_avoid_set(solid, own_if);
    let b_is_empty = crate::topo_tools::faces_of(solid).is_empty();

    let mut ivec: Vec<usize> = Vec::new();
    for n_fp in lifp {
        if n_fp >= v_sb.len() {
            continue;
        }
        let fp = &v_sb[n_fp].shape;
        if !msf.contains(&shape_key(fp)) {
            ivec.push(n_fp);
        }
    }
    ivec.sort_unstable();
    if ivec.is_empty() {
        return in_faces;
    }

    if b_is_empty {
        for k in ivec {
            in_faces.push(v_sb[k].shape.clone());
        }
        return in_faces;
    }

    let mut mefp: EdgeFaceMap = EdgeFaceMap::new();
    if ivec.len() > 1 {
        for &k in &ivec {
            map_edges_and_faces(&v_sb[k].shape, &mut mefp);
        }
    }

    // `BOPTools_AlgoTools::IsInternalFace` reads `theMEF.Contains` /
    // `ChangeFromKey` only, so the node lists suffice; the map is
    // `TopExp::MapShapesAndAncestors(theSolid, EDGE, FACE)` built once
    // (`BOPTools_AlgoTools.cxx:811`, keyed by `TopTools_ShapeMapHasher`).
    let mut mefds: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let mut mf_done: HashSet<ShapeMapKey> = HashSet::new();
    for &n_fp in &ivec {
        let fp = v_sb[n_fp].shape.clone();
        if !mf_done.insert(shape_key(&fp)) {
            continue;
        }
        let block = make_connexity_block(&fp, &mse, &mefp, &mut mf_done);
        if block_vertices_out_of_solid(&block.faces, box_s) {
            continue;
        }
        let face = classification_face(&block, &fp);
        if mefds.is_empty() {
            mefds = origin_map_shapes_and_ancestors(solid, crate::abs::ShapeType::Edge, crate::abs::ShapeType::Face);
        }
        let is_in = is_internal_face(&Face(face.clone()), solid, &mefds, CONFUSION, ctx);
        if is_in {
            in_faces.extend(block.faces);
        }
    }
    in_faces
}

/// ClassifyFaces that uses the AABB-tree Perform instead of the linear
/// `fill_in_3d_parts` wrapper.
pub fn classify_faces_aabb(
    faces: &[TopoShape],
    solids: &[TopoShape],
    ctx: &IntToolsContext,
    shape_box_map: &HashMap<ShapeMapKey, BndBox>,
    solids_if: &HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> HashMap<ShapeMapKey, Vec<TopoShape>> {
    let v_sb = make_shape_box_vector(faces, shape_box_map);
    let mut in_parts: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    for solid in solids {
        let box_s = solid_box_for_classify(solid, shape_box_map);
        let own = solids_if
            .get(&shape_key(solid))
            .cloned()
            .unwrap_or_default();
        let in_faces = perform_fill_in_3d_parts_tree(solid, &box_s, &own, &v_sb, ctx);
        in_parts.insert(shape_key(solid), in_faces);
    }
    in_parts
}
