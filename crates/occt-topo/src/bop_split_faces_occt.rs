//! `BOPAlgo_Builder::BuildSplitFaces` — full collection + rebuild.
//!
//! Source: `BOPAlgo_Builder_2.cxx:233-555`. Area construction is
//! `BOPAlgo_BuilderFace::Perform` (`SetFace` / `SetShapes` / `Perform`),
//! which runs ShapesToAvoid → Loops → Areas → InternalShapes.
//!
//! Flow (`_2.cxx`):
//! 1. For every source FACE with a face-info record:
//!    * if no IN, ON, Sc paves and no alone vertices → skip;
//!    * if no IN/Sc → try [`crate::bop_draft_face::split_face_fast_path`];
//!    * otherwise collect bounding edges (images, INTERNAL doubled, closed
//!      seams doubled, open splits WithWarn-reversed) plus IN/Sc both ways;
//! 2. `BOPAlgo_BuilderFace::Perform` (ShapesToAvoid → Loops → Areas → Internal);
//! 3. Record images, re-applying the source face orientation.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::bop_build_faces::{append_on_face_edges, attach_pcurves, BopBuilderLike};
use crate::bop_draft_face::split_face_fast_path;
use crate::bop_geomlib_closed::split_face_surface_closed;
use crate::bop_occt_util::{
    alone_vertices, face_info_of, has_face_info, nb_source, unique_edge_paves,
};
use crate::bop_split_seam::{
    do_split_seam_on_face, do_split_seam_on_face_origin, is_closed_on_face, isoline_uv,
};
use crate::bop_split_to_reverse::orient_split_from_base_with_warn;
use crate::brep_tool::BRepTool;
use crate::builder_area::AreaBuilder;
use crate::builder_face::FaceBuilder;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::edges_of;

/// `GeomLib::IsClosed` stand-in: a surface is U/V closed when it is periodic
/// Whether a source edge is a seam on `face` that must be doubled on the
/// split (`_2.cxx:395-404`): the face is U or V closed, the edge is closed
/// on the face, and the p-curve is the matching isoline.
fn bounding_edge_is_closed_seam(edge: &Edge, face: &Face, is_u_closed: bool, is_v_closed: bool) -> bool {
    if !(is_u_closed || is_v_closed) {
        return false;
    }
    if !is_closed_on_face(edge, face) {
        return false;
    }
    let (is_u_iso, is_v_iso) = isoline_uv(edge, face);
    (is_u_closed && is_u_iso) || (is_v_closed && is_v_iso)
}

/// Make a split of a closed edge into a seam on `face` (`_2.cxx:429-446`).
/// Tries `DoSplitSEAMOnFace(split, face)` then the origin/split overload.
/// Returns true when the split is (now) closed on the face.
fn ensure_closed_split_on_face(origin: &Edge, split: &Edge, face: &Face) -> bool {
    if is_closed_on_face(split, face) {
        return true;
    }
    if do_split_seam_on_face(split, face) {
        return true;
    }
    if do_split_seam_on_face_origin(origin, split, face) {
        return true;
    }
    false
}

/// One BuilderFace task: the source face and the unordered edge set `LE`.
struct SplitFaceTask {
    face_index: usize,
    face: Face,
    edges: Vec<Edge>,
}

/// Collect the bounding + IN + Sc edges of `face` (`_2.cxx:353-494`).
fn collect_split_edge_set<B: BopBuilderLike>(
    f: &B,
    face: &Face,
    paves_in: &[(usize, f64, f64)],
    paves_sc: &[(usize, f64, f64)],
    ctx: &mut IntToolsContext,
    add_warning: &mut dyn FnMut(String),
) -> Vec<Edge> {
    let mut le: Vec<Edge> = Vec::new();
    let mut a_m_fence: HashSet<usize> = HashSet::new();
    let mut is_u_closed = false;
    let mut is_v_closed = false;
    let mut is_checked = false;

    let mut a_ff = face.clone();
    a_ff.0.set_orientation(Orientation::Forward);

    for be in edges_of(&a_ff.0) {
        let an_ori_e = be.0.orientation();
        if !f.history().has_image(&be.0) {
            if an_ori_e == Orientation::Internal {
                let mut fwd = be.clone();
                fwd.0.set_orientation(Orientation::Forward);
                le.push(fwd);
                let mut rev = be.clone();
                rev.0.set_orientation(Orientation::Reversed);
                le.push(rev);
            } else {
                le.push(be);
            }
            continue;
        }

        if !is_checked {
            let closed = split_face_surface_closed(face, &be);
            is_u_closed = closed.0;
            is_v_closed = closed.1;
            is_checked = true;
        }

        let b_is_closed = bounding_edge_is_closed_seam(&be, face, is_u_closed, is_v_closed);
        let b_is_degenerated = BRepTool::is_degenerated(&be);
        let splits = f.history().image(&be.0).map(|s| s.to_vec()).unwrap_or_default();
        for sp_shape in splits {
            let mut a_sp = sp_shape;
            if b_is_degenerated {
                a_sp.set_orientation(an_ori_e);
                le.push(Edge(a_sp));
                continue;
            }
            if an_ori_e == Orientation::Internal {
                a_sp.set_orientation(Orientation::Forward);
                le.push(Edge(a_sp.clone()));
                a_sp.set_orientation(Orientation::Reversed);
                le.push(Edge(a_sp));
                continue;
            }
            if b_is_closed {
                let sp_key = GeometryRegistry::shape_key(&a_sp);
                if a_m_fence.insert(sp_key) {
                    let sp_edge = Edge(a_sp.clone());
                    if !is_closed_on_face(&sp_edge, face)
                        && !ensure_closed_split_on_face(&be, &sp_edge, face)
                    {
                        add_warning(format!(
                            "BOPAlgo_AlertUnableToMakeClosedEdgeOnFace: split of closed edge on face {}",
                            GeometryRegistry::shape_key(&face.0)
                        ));
                    }
                    a_sp.set_orientation(Orientation::Forward);
                    le.push(Edge(a_sp.clone()));
                    a_sp.set_orientation(Orientation::Reversed);
                    le.push(Edge(a_sp));
                }
                continue;
            }
            let (oriented, warn) =
                orient_split_from_base_with_warn(&a_sp, &be.0, an_ori_e, ctx);
            if let Some(w) = warn {
                add_warning(w.message);
            }
            le.push(Edge(oriented));
        }
    }

    let mut on_le: Vec<Edge> = Vec::new();
    append_on_face_edges(f, face, ctx, paves_in, &mut le, &mut on_le);
    append_on_face_edges(f, face, ctx, paves_sc, &mut le, &mut on_le);
    let _ = on_le;
    le
}

fn rebuild_split_areas(face: &Face, edges: &[Edge]) -> Result<Vec<TopoShape>, String> {
    attach_pcurves(face, edges);
    let mut ff = face.clone();
    ff.0.set_orientation(Orientation::Forward);
    let shapes: Vec<TopoShape> = edges.iter().map(|e| e.0.clone()).collect();
    let mut fb = FaceBuilder::new();
    fb.set_face(&ff);
    fb.set_shapes(&shapes);
    crate::builder_face_occt::perform(&mut fb)?;
    for area in fb.areas() {
        attach_pcurves(&Face(area.clone()), edges);
    }
    Ok(fb.areas().to_vec())
}

/// Record split faces as images of the source, applying `anOriF` (`_2.cxx:534-552`).
fn bind_face_images<B: BopBuilderLike>(f: &mut B, face_idx: usize, areas: Vec<TopoShape>) {
    let Some(orig) = f.ds().shape(face_idx).cloned() else {
        return;
    };
    let an_ori = orig.orientation();
    for sf in areas {
        let mut sf = sf;
        if an_ori == Orientation::Reversed {
            sf.set_orientation(Orientation::Reversed);
        }
        f.history_mut().add_image(&orig, sf);
    }
}

/// `BOPAlgo_Builder::BuildSplitFaces`.
pub fn build_split_faces_occt<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    let n = nb_source(f.ds());
    let mut ctx = IntToolsContext::new();
    let mut tasks: Vec<SplitFaceTask> = Vec::new();
    let mut draft_images: HashMap<usize, Vec<TopoShape>> = HashMap::new();

    for i in 0..n {
        let is_face = f
            .ds()
            .shape_info(i)
            .map(|s| s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        if !is_face {
            continue;
        }
        if !has_face_info(f.ds(), i) {
            continue;
        }
        let Some(face_shape) = f.ds().shape(i).cloned() else {
            continue;
        };
        let face = Face(face_shape);
        let Some(fi) = face_info_of(f.ds(), i) else {
            continue;
        };

        let paves_in = unique_edge_paves(f.ds(), fi.paves_in());
        let paves_on = unique_edge_paves(f.ds(), fi.paves_on());
        let paves_sc = unique_edge_paves(f.ds(), fi.paves());
        let a_liav = alone_vertices(f.ds(), i);
        let a_nb_pb_in = paves_in.len();
        let a_nb_pb_on = paves_on.len();
        let a_nb_pb_sc = paves_sc.len();
        let a_nb_av = a_liav.len();
        if a_nb_pb_in == 0 && a_nb_pb_on == 0 && a_nb_pb_sc == 0 && a_nb_av == 0 {
            continue;
        }

        if a_nb_pb_in == 0 && a_nb_pb_sc == 0 {
            let mut warnings: Vec<String> = Vec::new();
            let path = split_face_fast_path(
                &face,
                a_nb_pb_in,
                a_nb_pb_sc,
                a_nb_av,
                f.history(),
                &mut ctx,
                |msg| warnings.push(msg),
            );
            for w in warnings {
                f.add_warning(w);
            }
            match path {
                crate::bop_draft_face::SplitFaceFastPath::Skip => continue,
                crate::bop_draft_face::SplitFaceFastPath::Draft(fd) => {
                    draft_images.entry(i).or_default().push(fd.0);
                    continue;
                }
                crate::bop_draft_face::SplitFaceFastPath::BuilderFace => {}
            }
        }

        let mut warnings: Vec<String> = Vec::new();
        let le = collect_split_edge_set(
            f,
            &face,
            &paves_in,
            &paves_sc,
            &mut ctx,
            &mut |msg| warnings.push(msg),
        );
        for w in warnings {
            f.add_warning(w);
        }
        tasks.push(SplitFaceTask {
            face_index: i,
            face,
            edges: le,
        });
    }

    let mut faces_im: HashMap<usize, Vec<TopoShape>> = draft_images;
    for task in tasks {
        match rebuild_split_areas(&task.face, &task.edges) {
            Ok(areas) => {
                faces_im.entry(task.face_index).or_default().extend(areas);
            }
            Err(e) => {
                f.add_error(format!(
                    "BuildSplitFaces: face {} rebuild failed: {e}",
                    task.face_index
                ));
            }
        }
    }

    for (face_idx, areas) in faces_im {
        bind_face_images(f, face_idx, areas);
    }
    Ok(())
}
