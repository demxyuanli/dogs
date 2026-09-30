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
    wires_of_face,
};
use crate::bop_split_seam::{do_split_seam_on_face, is_closed_on_face};
use crate::bop_split_to_reverse::{emit_warnings, orient_split_with_warn, ReverseWarn};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::edge_vertices;

/// Vertex → incident edges map used by `HasMultiConnected`.
///
/// OCCT stores `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// keyed by the vertex. The list is expected to stay 1–2 elements long; a
/// third distinct edge means the vertex can split a thin face, so the draft
/// path is abandoned.
pub type VertexEdgeMap = HashMap<usize, Vec<usize>>;

/// `HasMultiConnected` (`_2.cxx:1014`).
///
/// Walks the vertices of `the_edge`. For each vertex, records the edge in
/// `the_map` (by TShape key). Returns true as soon as any vertex lists more
/// than two distinct edges.
pub fn has_multi_connected(the_edge: &Edge, the_map: &mut VertexEdgeMap) -> bool {
    let e_key = shape_key(&the_edge.0);
    let (va, vb) = edge_vertices(the_edge);
    for v in [va, vb].into_iter().flatten() {
        let vk = shape_key(&v.0);
        let list = the_map.entry(vk).or_default();
        if !list.contains(&e_key) {
            list.push(e_key);
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
/// `the_images` is `myImages`. `the_ctx` is the shared `IntTools_Context`.
/// `the_report` collects WithWarn alerts.
pub fn build_draft_face(
    the_face: &Face,
    the_images: &BopHistory,
    the_ctx: &mut IntToolsContext,
    the_report: &mut DraftFaceReport,
) -> DraftFaceResult {
    let Some(a_s) = BRepTool::face_surface(the_face) else {
        return DraftFaceResult::Fallback;
    };
    let bld = TopoBuilder::new();
    let mut a_draft_face = bld.make_face(a_s, &[]);
    a_draft_face.0.set_location(the_face.0.location());

    let mut a_vertices_counter: VertexEdgeMap = HashMap::new();
    let mut a_m_edges: HashSet<usize> = HashSet::new();

    let mut a_ff = the_face.clone();
    a_ff.0.set_orientation(Orientation::Forward);

    for a_w in wires_of_face(&a_ff) {
        if a_w.0.shape_type() != ShapeType::Wire {
            continue;
        }
        let mut a_w_fwd = a_w.clone();
        a_w_fwd.0.set_orientation(Orientation::Forward);
        let a_it_e: Vec<TopoShape> = iter_children(&a_w_fwd.0);
        if a_it_e.is_empty() {
            continue;
        }

        let mut a_new_wire = make_empty_wire(Orientation::Forward);
        for a_e_shape in a_it_e {
            let Some(a_e) = as_edge(&a_e_shape) else {
                continue;
            };
            let an_ori_e = a_e.0.orientation();
            if an_ori_e == Orientation::Internal {
                // Internal edges could split the original face on halves.
                return DraftFaceResult::Fallback;
            }
            let b_is_degenerated = BRepTool::is_degenerated(&a_e);
            let b_is_closed = is_closed_on_face(&a_e, the_face);

            if !the_images.has_image(&a_e.0) {
                if !b_is_degenerated && has_multi_connected(&a_e, &mut a_vertices_counter) {
                    return DraftFaceResult::Fallback;
                }
                if !b_is_closed && !a_m_edges.insert(shape_key(&a_e.0)) {
                    return DraftFaceResult::Fallback;
                }
                add_edge_to_wire(&mut a_new_wire, &a_e);
                continue;
            }

            let splits = the_images.image(&a_e.0).unwrap_or(&[]).to_vec();
            for sp_shape in splits {
                let Some(mut a_sp) = as_edge(&sp_shape) else {
                    continue;
                };
                if !b_is_degenerated && has_multi_connected(&a_sp, &mut a_vertices_counter) {
                    return DraftFaceResult::Fallback;
                }
                if !b_is_closed && !a_m_edges.insert(shape_key(&a_sp.0)) {
                    return DraftFaceResult::Fallback;
                }
                a_sp.0.set_orientation(an_ori_e);
                if b_is_degenerated {
                    add_edge_to_wire(&mut a_new_wire, &a_sp);
                    continue;
                }
                if b_is_closed && !is_closed_on_face(&a_sp, the_face) {
                    let _ = do_split_seam_on_face(&a_sp, the_face);
                }
                if let Some(w) = orient_split_with_warn(&mut a_sp.0, &a_e.0, the_ctx) {
                    the_report.warnings.push(w);
                }
                add_edge_to_wire(&mut a_new_wire, &a_sp);
            }
        }
        a_new_wire.0.set_orientation(a_w.0.orientation());
        a_new_wire.0.set_closed(crate::bop_occt_util::wire_is_closed(&a_new_wire));
        add_wire_to_face(&mut a_draft_face, &a_new_wire);
    }

    if the_face.0.orientation() == Orientation::Reversed {
        a_draft_face.0.reverse();
    }
    DraftFaceResult::Draft(a_draft_face)
}

/// Build a draft face and emit warnings through `add_warning`.
pub fn build_draft_face_reported<F: FnMut(String)>(
    the_face: &Face,
    the_images: &BopHistory,
    the_ctx: &mut IntToolsContext,
    mut add_warning: F,
) -> DraftFaceResult {
    let mut report = DraftFaceReport::default();
    let result = build_draft_face(the_face, the_images, the_ctx, &mut report);
    emit_warnings(&report.warnings, &mut add_warning);
    result
}

/// True when the face has any INTERNAL bounding edge (the check BuildSplitFaces
/// does before calling BuildDraftFace when there are no alone vertices).
pub fn face_has_internal_edges(face: &Face) -> bool {
    for w in wires_of_face(face) {
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
    for w in wires_of_face(face) {
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
    for w in wires_of_face(face) {
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
    wires_of_face(face).len()
}
