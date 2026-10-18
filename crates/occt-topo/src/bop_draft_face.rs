//! `BOPAlgo_Builder::BuildDraftFace` and `HasMultiConnected`.
//!
//! Source: `BOPAlgo_Builder_2.cxx` (`HasMultiConnected` at 1014,
//! `BuildDraftFace` at 1052). BuildSplitFaces takes this fast path when a
//! face has no IN or section pave blocks (`_2.cxx:298-350`): the wires are
//! rebuilt from edge images, the original surface is kept, and the result is
//! stored as the face image without running `BOPAlgo_BuilderFace`.
//!
//! The draft is null (caller falls back to BuilderFace) when:
//! * a bounding edge is INTERNAL (it may split the face);
//! * a vertex is multi-connected (`HasMultiConnected`, more than two edges);
//! * two distinct non-closed edges unify onto the same image (the new wire
//!   would need a validity check the draft path does not perform);
//! * a closed original edge has a split that is not yet a seam and
//!   `DoSplitSEAMOnFace` is attempted (`_2.cxx:1163`).
//!
//! Same-domain / reverse: each split is oriented with
//! [`crate::bop_split_to_reverse::orient_split_with_warn`].
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::bop_hist::BopHistory;
use crate::bop_occt_util::{
    add_edge_to_wire, add_wire_to_face, as_edge, iter_children, make_empty_wire, shape_key,
    bop_wires_of_face,
};
use crate::bop_split_seam::{do_split_seam_on_face, is_closed_on_face};
use crate::bop_split_to_reverse::{emit_warnings, orient_split_with_warn, ReverseWarn};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{edge_vertices, ShapeMapKey};

/// Vertex → incident edges map used by `HasMultiConnected`.
///
/// OCCT stores `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// keyed by the vertex. The list is expected to stay 1–2 elements long; a
/// third distinct edge means the vertex can split a thin face, so the draft
/// path is abandoned.
pub type VertexEdgeMap = HashMap<ShapeMapKey, Vec<ShapeMapKey>>;

/// `HasMultiConnected` (`_2.cxx:1014`).
///
/// Walks the vertices of `edge`. For each vertex, records the edge in
/// `the_map` (by TShape key). Returns true as soon as any vertex lists more
/// than two distinct edges.
pub fn has_multi_connected(edge: &Edge, the_map: &mut VertexEdgeMap) -> bool {
    let e_key = shape_key(&edge.0);
    let (va, vb) = edge_vertices(edge);
    for v in [va, vb].into_iter().flatten() {
        let vk = shape_key(&v.0);
        let list = the_map.entry(vk).or_default();
        if !list.contains(&e_key) {
            list.push(e_key.clone());
        }
        if list.len() > 2 {
            return true;
        }
    }
    false
}

/// Result of a draft-face attempt.
#[derive(Debug, Clone)]
pub enum DraftFaceResult {
    /// A valid draft face ready to become the image of the original.
    Draft(Face),
    /// The draft path cannot be used; BuildSplitFaces must run BuilderFace.
    Fallback,
}

/// Warnings collected while building a draft face.
#[derive(Debug, Clone, Default)]
pub struct DraftFaceReport {
    pub warnings: Vec<ReverseWarn>,
}

/// `BuildDraftFace` (`_2.cxx:1052`).
///
/// `images` is `myImages`. `ctx` is the shared `IntTools_Context`.
/// `report` collects WithWarn alerts.
pub fn build_draft_face(
    face: &Face,
    images: &BopHistory,
    ctx: &mut IntToolsContext,
    report: &mut DraftFaceReport,
) -> DraftFaceResult {
    let Some(s) = BRepTool::face_surface(face) else {
        return DraftFaceResult::Fallback;
    };
    let bld = TopoBuilder::new();
    let mut draft_face = bld.make_face(s, &[]);
    draft_face.0.set_location(face.0.location());

    let mut vertices_counter: VertexEdgeMap = HashMap::new();
    let mut m_edges: HashSet<ShapeMapKey> = HashSet::new();

    let mut ff = face.clone();
    ff.0.set_orientation(Orientation::Forward);

    for a_w in bop_wires_of_face(&ff) {
        if a_w.0.shape_type() != ShapeType::Wire {
            continue;
        }
        let mut w_fwd = a_w.clone();
        w_fwd.0.set_orientation(Orientation::Forward);
        let it_e: Vec<TopoShape> = iter_children(&w_fwd.0);
        if it_e.is_empty() {
            continue;
        }

        let mut new_wire = make_empty_wire(Orientation::Forward);
        for e_shape in it_e {
            let Some(e) = as_edge(&e_shape) else {
                continue;
            };
            let an_ori_e = e.0.orientation();
            if an_ori_e == Orientation::Internal {
                // Internal edges could split the original face on halves.
                return DraftFaceResult::Fallback;
            }
            let b_is_degenerated = BRepTool::is_degenerated(&e);
            let b_is_closed = is_closed_on_face(&e, face);

            if !images.has_image(&e.0) {
                if !b_is_degenerated && has_multi_connected(&e, &mut vertices_counter) {
                    return DraftFaceResult::Fallback;
                }
                if !b_is_closed && !m_edges.insert(shape_key(&e.0)) {
                    return DraftFaceResult::Fallback;
                }
                add_edge_to_wire(&mut new_wire, &e);
                continue;
            }

            let splits = images.image(&e.0).unwrap_or(&[]).to_vec();
            for sp_shape in splits {
                let Some(mut sp) = as_edge(&sp_shape) else {
                    continue;
                };
                if !b_is_degenerated && has_multi_connected(&sp, &mut vertices_counter) {
                    return DraftFaceResult::Fallback;
                }
                if !b_is_closed && !m_edges.insert(shape_key(&sp.0)) {
                    return DraftFaceResult::Fallback;
                }
                sp.0.set_orientation(an_ori_e);
                if b_is_degenerated {
                    add_edge_to_wire(&mut new_wire, &sp);
                    continue;
                }
                if b_is_closed && !is_closed_on_face(&sp, face) {
                    let _ = do_split_seam_on_face(&sp, face);
                }
                if let Some(w) = orient_split_with_warn(&mut sp.0, &e.0, ctx) {
                    report.warnings.push(w);
                }
                add_edge_to_wire(&mut new_wire, &sp);
            }
        }
        new_wire.0.set_orientation(a_w.0.orientation());
        new_wire.0.set_closed(crate::bop_occt_util::bop_wire_is_closed(&new_wire));
        add_wire_to_face(&mut draft_face, &new_wire);
    }

    if face.0.orientation() == Orientation::Reversed {
        draft_face.0.reverse();
    }
    DraftFaceResult::Draft(draft_face)
}

/// Build a draft face and emit warnings through `add_warning`.
pub fn build_draft_face_reported<F: FnMut(String)>(
    face: &Face,
    images: &BopHistory,
    ctx: &mut IntToolsContext,
    mut add_warning: F,
) -> DraftFaceResult {
    let mut report = DraftFaceReport::default();
    let result = build_draft_face(face, images, ctx, &mut report);
    emit_warnings(&report.warnings, &mut add_warning);
    result
}

/// True when the face has any INTERNAL bounding edge (the check BuildSplitFaces
/// does before calling BuildDraftFace when there are no alone vertices).
pub fn face_has_internal_edges(face: &Face) -> bool {
    for w in bop_wires_of_face(face) {
        for e in iter_children(&w.0) {
            if e.shape_type() == ShapeType::Edge && e.orientation() == Orientation::Internal {
                return true;
            }
        }
    }
    false
}

/// True when any wire of the face has an image (`hasModified` in `_2.cxx:327`).
pub fn face_has_modified_wires(face: &Face, images: &BopHistory) -> bool {
    for w in bop_wires_of_face(face) {
        if images.has_image(&w.0) {
            return true;
        }
        for e in iter_children(&w.0) {
            if images.has_image(&e) {
                return true;
            }
        }
    }
    false
}

/// Decide whether BuildSplitFaces should try the draft-face fast path
/// (`_2.cxx:298-350`). Returns `true` when there are no IN/section paves
/// and the face is a candidate for `BuildDraftFace`.
pub fn should_try_draft_face(
    nb_pb_in: usize,
    nb_pb_sc: usize,
    nb_av: usize,
    has_internals: bool,
    has_modified: bool,
) -> bool {
    if nb_pb_in != 0 || nb_pb_sc != 0 {
        return false;
    }
    if nb_av == 0 && !has_internals && !has_modified {
        return false;
    }
    !has_internals
}

/// Combined predicate used by BuildSplitFaces: given the face, its pave
/// counts and the images table, either produce a draft, skip the face, or
/// fall through to BuilderFace.
#[derive(Debug, Clone)]
pub enum SplitFaceFastPath {
    Skip,
    Draft(Face),
    BuilderFace,
}

/// Evaluate the `_2.cxx:298-350` fast path for one source face.
pub fn split_face_fast_path<F: FnMut(String)>(
    face: &Face,
    nb_pb_in: usize,
    nb_pb_sc: usize,
    nb_av: usize,
    images: &BopHistory,
    ctx: &mut IntToolsContext,
    add_warning: F,
) -> SplitFaceFastPath {
    if nb_pb_in != 0 || nb_pb_sc != 0 {
        return SplitFaceFastPath::BuilderFace;
    }
    let has_internals = face_has_internal_edges(face);
    if nb_av == 0 {
        let has_modified = face_has_modified_wires(face, images);
        if !has_internals && !has_modified {
            return SplitFaceFastPath::Skip;
        }
    }
    if has_internals {
        return SplitFaceFastPath::BuilderFace;
    }
    match build_draft_face_reported(face, images, ctx, add_warning) {
        DraftFaceResult::Draft(fd) => SplitFaceFastPath::Draft(fd),
        DraftFaceResult::Fallback => SplitFaceFastPath::BuilderFace,
    }
}

/// Count the vertices in `the_map` that already have more than two edges.
pub fn multi_connected_vertex_count(the_map: &VertexEdgeMap) -> usize {
    the_map.values().filter(|l| l.len() > 2).count()
}

/// Reset helper used by tests and by callers that reuse a counter map.
pub fn clear_vertex_edge_map(the_map: &mut VertexEdgeMap) {
    the_map.clear();
}

/// Whether `edge` currently appears more than twice at any of its vertices
/// in `the_map` without mutating the map (a read-only probe).
pub fn vertex_is_multi_connected(the_map: &VertexEdgeMap, edge: &Edge) -> bool {
    let (va, vb) = edge_vertices(edge);
    for v in [va, vb].into_iter().flatten() {
        if the_map
            .get(&shape_key(&v.0))
            .map(|l| l.len() > 2)
            .unwrap_or(false)
        {
            return true;
        }
    }
    false
}

/// Edges of a draft face, in wire order, for diagnostics.
pub fn draft_face_edges(face: &Face) -> Vec<Edge> {
    let mut out = Vec::new();
    for w in bop_wires_of_face(face) {
        for e in iter_children(&w.0) {
            if let Some(ed) = as_edge(&e) {
                out.push(ed);
            }
        }
    }
    out
}

/// Number of wires on a draft face.
pub fn draft_face_wire_count(face: &Face) -> usize {
    bop_wires_of_face(face).len()
}
