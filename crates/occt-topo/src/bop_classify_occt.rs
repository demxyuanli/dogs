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

use std::collections::HashMap;

use occt_core::bnd::BndBox;

use crate::algo_tools::AlgoTools;
use crate::bop_occt_util::shape_key;
use crate::bopalgo_tools_class::{fill_in_3d_parts, ShapeBox};
use crate::bbox_from_geometry::shape_bbox;
use crate::int_tools_full::IntToolsContext;
use crate::shape::TopoShape;

/// A face plus its bounding box, the `BOPAlgo_ShapeBox` vector entry.
pub type FaceBox = ShapeBox;

/// Build the `BOPAlgo_VectorOfShapeBox` from `the_faces` and an optional
/// precomputed box map (`theShapeBoxMap`). Missing entries are computed
/// with `BRepBndLib::Add` (`_cxx:1660-1666`).
pub fn make_shape_box_vector(
    the_faces: &[TopoShape],
    the_shape_box_map: &HashMap<usize, BndBox>,
) -> Vec<ShapeBox> {
    let mut a_vsb = Vec::with_capacity(the_faces.len());
    for a_f in the_faces {
        let a_box = the_shape_box_map
            .get(&shape_key(a_f))
            .copied()
            .filter(|b| !b.is_void())
            .unwrap_or_else(|| shape_bbox(a_f));
        a_vsb.push(ShapeBox::new(a_f.clone(), a_box));
    }
    a_vsb
}

/// AABB tree stand-in: indices of `the_v_sb` whose boxes are not out of
/// `the_box_s` (`BOPTools_BoxTreeSelector::Select`, `_cxx:1345-1354`).
pub fn aabb_select(the_box_s: &BndBox, the_v_sb: &[ShapeBox]) -> Vec<usize> {
    let mut idx = Vec::new();
    for (i, sb) in the_v_sb.iter().enumerate() {
        if !the_box_s.is_out_box(&sb.box_) {
            idx.push(i);
        }
    }
    idx
}

/// Solid box used by ClassifyFaces (`_cxx:1693-1712`): prefer the map as-is.
/// Invert-to-whole only when the box is computed here (`BRepBndLib::Add`).
pub fn solid_box_for_classify(
    the_solid: &TopoShape,
    the_shape_box_map: &HashMap<usize, BndBox>,
) -> BndBox {
    if let Some(b) = the_shape_box_map.get(&shape_key(the_solid)).copied() {
        if !b.is_void() {
            return b;
        }
    }
    let mut a_box = shape_bbox(the_solid);
    if !a_box.is_whole() && AlgoTools::is_inverted_solid(the_solid) {
        a_box.set_whole();
    }
    a_box
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
    the_faces: &[TopoShape],
    the_solids: &[TopoShape],
    the_ctx: &IntToolsContext,
    the_shape_box_map: &HashMap<usize, BndBox>,
    the_solids_if: &HashMap<usize, Vec<TopoShape>>,
) -> HashMap<usize, Vec<TopoShape>> {
    classify_faces_aabb(
        the_faces,
        the_solids,
        the_ctx,
        the_shape_box_map,
        the_solids_if,
    )
}

/// ClassifyFaces that also returns the per-solid tasks (diagnostics).
pub fn classify_faces_with_tasks(
    the_faces: &[TopoShape],
    the_solids: &[TopoShape],
    the_ctx: &IntToolsContext,
    the_shape_box_map: &HashMap<usize, BndBox>,
    the_solids_if: &HashMap<usize, Vec<TopoShape>>,
) -> (HashMap<usize, Vec<TopoShape>>, Vec<FillIn3dTask>) {
    let v_sb = make_shape_box_vector(the_faces, the_shape_box_map);
    let mut tasks = Vec::new();
    let mut in_parts: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for a_solid in the_solids {
        let box_s = solid_box_for_classify(a_solid, the_shape_box_map);
        let own = the_solids_if
            .get(&shape_key(a_solid))
            .cloned()
            .unwrap_or_default();
        let task = FillIn3dTask::new(a_solid.clone(), box_s, own);
        let in_faces = task.perform(&v_sb, the_ctx);
        in_parts.insert(shape_key(a_solid), in_faces);
        tasks.push(task);
    }
    let _ = v_sb;
    (in_parts, tasks)
}

/// Select-and-classify a single solid (the sequential body of one
/// `FillIn3DParts` parallel task).
pub fn classify_one_solid(
    the_solid: &TopoShape,
    the_faces: &[TopoShape],
    the_own_if: &[TopoShape],
    the_shape_box_map: &HashMap<usize, BndBox>,
    the_ctx: &IntToolsContext,
) -> Vec<TopoShape> {
    let v_sb = make_shape_box_vector(the_faces, the_shape_box_map);
    let box_s = solid_box_for_classify(the_solid, the_shape_box_map);
    fill_in_3d_parts(the_solid, &box_s, the_own_if, &v_sb, the_ctx)
}

/// Merge `the_in_faces` with `the_internal` as FillIn3DParts does when
/// writing `myInParts` (`_3.cxx:243-260`).
pub fn combine_in_and_internal(
    the_in_faces: &[TopoShape],
    the_internal: &[TopoShape],
) -> Vec<TopoShape> {
    let mut out = Vec::with_capacity(the_in_faces.len() + the_internal.len());
    out.extend(the_in_faces.iter().cloned());
    out.extend(the_internal.iter().cloned());
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
/// 1. Select face boxes not out of `the_box_s` (`BoxTreeSelector`);
/// 2. Map solid edges/faces (`aMSE`, `aMSF`) and add own INTERNAL faces;
/// 3. Filter selected faces that are not already of the solid (`aIVec`);
/// 4. Empty solid → every remaining face is IN;
/// 5. Connexity blocks that do not cross solid edges;
/// 6. Vertex-box OUT cull;
/// 7. `IsInternalFace` on the representative, then the whole block is IN.
pub fn perform_fill_in_3d_parts_tree(
    the_solid: &TopoShape,
    the_box_s: &BndBox,
    the_own_if: &[TopoShape],
    the_v_sb: &[ShapeBox],
    the_ctx: &IntToolsContext,
) -> Vec<TopoShape> {
    use crate::algo_tools_face::is_internal_face;
    use crate::bop_aabb_faces::AabbTree;
    use crate::bop_connexity_faces::{
        block_vertices_out_of_solid, classification_face, make_connexity_block,
        map_edges_and_faces, solid_edge_avoid_set, solid_face_avoid_set, EdgeFaceMap,
    };
    use crate::bop_occt_util::map_shapes_and_ancestors;
    use crate::shape::Face;
    use occt_core::precision::CONFUSION;
    use std::collections::HashSet;

    let mut my_in_faces: Vec<TopoShape> = Vec::new();
    let tree = AabbTree::from_shape_boxes(the_v_sb);
    let a_lifp = if tree.is_empty() {
        Vec::new()
    } else {
        tree.select(the_box_s)
    };
    if a_lifp.is_empty() {
        return my_in_faces;
    }

    let a_mse = solid_edge_avoid_set(the_solid);
    let a_msf = solid_face_avoid_set(the_solid, the_own_if);
    let b_is_empty = crate::topo_tools_full::faces_of(the_solid).is_empty();

    let mut a_ivec: Vec<usize> = Vec::new();
    for n_fp in a_lifp {
        if n_fp >= the_v_sb.len() {
            continue;
        }
        let a_fp = &the_v_sb[n_fp].shape;
        if !a_msf.contains(&shape_key(a_fp)) {
            a_ivec.push(n_fp);
        }
    }
    a_ivec.sort_unstable();
    if a_ivec.is_empty() {
        return my_in_faces;
    }

    if b_is_empty {
        for k in a_ivec {
            my_in_faces.push(the_v_sb[k].shape.clone());
        }
        return my_in_faces;
    }

    let mut a_mefp: EdgeFaceMap = EdgeFaceMap::new();
    if a_ivec.len() > 1 {
        for &k in &a_ivec {
            map_edges_and_faces(&the_v_sb[k].shape, &mut a_mefp);
        }
    }

    let mut a_mefds: std::collections::HashMap<usize, (TopoShape, Vec<TopoShape>)> =
        HashMap::new();
    let mut a_mf_done: HashSet<usize> = HashSet::new();
    for &n_fp in &a_ivec {
        let a_fp = the_v_sb[n_fp].shape.clone();
        if !a_mf_done.insert(shape_key(&a_fp)) {
            continue;
        }
        let block = make_connexity_block(&a_fp, &a_mse, &a_mefp, &mut a_mf_done);
        if block_vertices_out_of_solid(&block.faces, the_box_s) {
            continue;
        }
        let a_face = classification_face(&block, &a_fp);
        if a_mefds.is_empty() {
            a_mefds = map_shapes_and_ancestors(the_solid, crate::abs::ShapeType::Edge, crate::abs::ShapeType::Face)
                .into_iter()
                .map(|(k, faces)| {
                    let e = crate::bop_occt_util::explore(the_solid, crate::abs::ShapeType::Edge)
                        .into_iter()
                        .find(|e| shape_key(e) == k)
                        .unwrap_or_else(|| faces.first().cloned().unwrap_or_else(|| a_fp.clone()));
                    (k, (e, faces))
                })
                .collect();
        }
        let is_in = is_internal_face(&Face(a_face.clone()), the_solid, &a_mefds, CONFUSION, the_ctx);
        if is_in {
            my_in_faces.extend(block.faces);
        }
    }
    my_in_faces
}

/// ClassifyFaces that uses the AABB-tree Perform instead of the linear
/// `fill_in_3d_parts` wrapper.
pub fn classify_faces_aabb(
    the_faces: &[TopoShape],
    the_solids: &[TopoShape],
    the_ctx: &IntToolsContext,
    the_shape_box_map: &HashMap<usize, BndBox>,
    the_solids_if: &HashMap<usize, Vec<TopoShape>>,
) -> HashMap<usize, Vec<TopoShape>> {
    let v_sb = make_shape_box_vector(the_faces, the_shape_box_map);
    let mut the_in_parts: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for a_solid in the_solids {
        let box_s = solid_box_for_classify(a_solid, the_shape_box_map);
        let own = the_solids_if
            .get(&shape_key(a_solid))
            .cloned()
            .unwrap_or_default();
        let in_faces = perform_fill_in_3d_parts_tree(a_solid, &box_s, &own, &v_sb, the_ctx);
        the_in_parts.insert(shape_key(a_solid), in_faces);
    }
    the_in_parts
}
