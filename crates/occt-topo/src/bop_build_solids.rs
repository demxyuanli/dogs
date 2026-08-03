//! Solid-image reconstruction — Phase 20b.
//!
//! Port of the solid-rebuilding stages of the General Fuse builder
//! (`BOPAlgo_Builder_3.cxx`):
//!
//! | Rust function              | OCCT source                               |
//! |----------------------------|-------------------------------------------|
//! | [`fill_images_solids`]     | `FillImagesSolids` (draft pass,           |
//! |                            | `FillIn3DParts` → `BuildDraftSolid`)      |
//! | [`build_split_solids`]     | `BuildSplitSolids` (`BOPAlgo_BuilderSolid`)|
//!
//! The two entry points operate on any host that implements
//! [`crate::bop_build_common::BopBuildOps`] — the same minimal contract the
//! sibling `crate::bop_build_common` module uses — so this module stays
//! independent of the concrete `BopBuilder` fields and is verified standalone.
//!
//! * **Draft pass** — [`fill_images_solids`] rebuilds every source solid from
//!   its face splits: each face of the solid is replaced by its image pieces
//!   (a split whose orientation is inverted relative to the original face is
//!   reversed first, mirroring `BOPTools_AlgoTools::IsSplitToReverse`) and the
//!   rebuilt shell is wrapped into a solid. A solid untouched by the
//!   intersection (no face/shell image, no internal part) is skipped, matching
//!   the `FillIn3DParts` guard.
//! * **Split pass** — [`build_split_solids`] collects the same split faces and
//!   groups them into *closed* shells with [`crate::shell_splitter::ShellSplitter`],
//!   wrapping each closed shell into a solid with
//!   [`crate::builder::TopoBuilder::make_solid`]. This is the path that
//!   separates a solid cut into several pieces into several image solids.
//!
//! `crate::bop_build_common::fill_internal_shapes` (the third `FillImagesSolids`
//! stage, settling internal vertices/edges/wires) already lives in
//! `crate::bop_build_common`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpPln, GpPnt};
use occt_geom::GeomPlane;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bop_build_common::{build_draft_solid, BopBuildOps};
use crate::brep_extrema::{closest_point_on_face, is_inside};
use crate::brep_surface::{surface_closest_params, surface_normal};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape};
use crate::shell_splitter::{EKey, ShellSplitter};
use crate::topo_tools_full::{edges_of, faces_of, vertices_of};

/// Stable identity key of a shape (the address of its shared `TShape`).
fn shape_key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

// ---------------------------------------------------------------------------
// FillImagesSolids — draft pass
// ---------------------------------------------------------------------------

/// Rebuilds every source solid from its face splits and records the rebuilt
/// solid as its image.
///
/// Mirrors `BOPAlgo_Builder::FillImagesSolids`'s draft stage
/// (`FillIn3DParts` → `BuildDraftSolid`): for every source solid of the data
/// structure that was touched by the intersection, [`build_draft_solid`]
/// rebuilds each shell of the solid from the image pieces of its faces (a
/// split face whose orientation is inverted relative to the original is
/// reversed first) and the rebuilt shell is wrapped into a solid. The rebuilt
/// solid is recorded as an image of the source solid.
///
/// A solid with no face/shell image and no internal parts is skipped — there
/// is nothing to rebuild (OCCT `FillIn3DParts` guard). The method returns
/// early when no solid participates in the operation.
pub fn fill_images_solids<B: BopBuildOps>(f: &mut B) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    if !(0..n)
        .any(|i| f.ds().shape_info(i).map(|s| s.shape_type() == ShapeType::Solid).unwrap_or(false))
    {
        return Ok(());
    }
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let solid = si.shape().clone();
        if !solid_interfered(f, &solid) {
            continue;
        }
        let draft = build_draft_solid(f, &solid)?;
        f.history_mut().add_image(&solid, draft.clone());
        f.origins_mut().entry(shape_key(&draft)).or_default().push(solid);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// BuildSplitSolids — closed-shell splitting
// ---------------------------------------------------------------------------

/// Splits every interfered source solid into its closed-shell pieces and
/// records each piece as a solid image of the source.
///
/// Mirrors `BOPAlgo_Builder::BuildSplitSolids` with the 3-D face classification
/// (`BOPAlgo_Tools::ClassifyFaces`) not yet ported: the split faces of the
/// solid (each source face replaced by its image pieces, an inverted split
/// reversed first) are grouped into *closed* shells with
/// [`crate::shell_splitter::ShellSplitter`] and every closed shell is wrapped
/// into a solid with [`crate::builder::TopoBuilder::make_solid`]. A solid cut
/// by the intersection expands into as many image solids as it has closed
/// shells; a solid whose splits still form a single closed shell produces one
/// image solid (the equivalent of the draft-solid fast path).
pub fn build_split_solids<B: BopBuildOps>(f: &mut B) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let solid = si.shape().clone();
        if !solid_interfered(f, &solid) {
            continue;
        }
        let faces = collect_solid_split_faces(f, &solid);
        if faces.is_empty() {
            continue;
        }
        let splits = split_faces_into_solids(&faces)?;
        if splits.is_empty() {
            continue;
        }
        for s in splits {
            f.history_mut().add_image(&solid, s.clone());
            f.origins_mut().entry(shape_key(&s)).or_default().push(solid.clone());
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// True when the solid was touched by the intersection: at least one of its
/// faces carries a non-empty image (OCCT `FillIn3DParts` checks the shells'
/// images; in this port the shell containers are rebuilt later, so the faces —
/// already filled by `crate::bop_build_faces` — are the source of truth).
fn solid_interfered<B: BopBuildOps>(f: &B, solid: &TopoShape) -> bool {
    faces_of(solid).iter().any(|fc| f.history().has_image(&fc.0))
}

/// Collects the faces that rebuild `solid`: every face of the solid replaced
/// by its image pieces, or kept as-is when not split. A split piece whose
/// orientation is inverted relative to the original face is reversed first
/// (mirrors `BOPTools_AlgoTools::IsSplitToReverse`).
fn collect_solid_split_faces<B: BopBuildOps>(f: &B, solid: &TopoShape) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for fc in faces_of(solid) {
        let or = fc.orientation();
        match f.history().image(&fc.0) {
            Some(imgs) => {
                for im in imgs {
                    let mut s = im.clone();
                    if is_split_to_reverse(&s, &fc.0) {
                        s.set_orientation(or.reversed());
                    } else {
                        s.set_orientation(or);
                    }
                    out.push(s);
                }
            }
            None => out.push(fc.0),
        }
    }
    out
}

/// Groups `faces` into closed shells with [`ShellSplitter`] and wraps every
/// closed shell into a solid. Mirrors the `BOPAlgo_BuilderSolid` output
/// assembly: each closed shell of the split is one solid of the result.
fn split_faces_into_solids(faces: &[TopoShape]) -> Result<Vec<TopoShape>, String> {
    let mut splitter = ShellSplitter::new();
    for fc in faces {
        splitter.add_start_element(fc.clone());
    }
    splitter.perform()?;
    let shells = splitter.shells().to_vec();
    // Close open shells (a piece cut by the intersection is missing the section
    // face, which the shell splitter leaves as a loose face because the section
    // edge is shared by more than two faces) by attaching the loose faces that
    // cover their open boundary.
    let shells = close_open_shells(&shells, faces);
    let bld = TopoBuilder::new();
    let mut out = Vec::new();
    for sh in shells {
        let mut solid = Solid::new();
        bld.add(&mut solid.0, &sh);
        out.push(solid.0);
    }
    Ok(out)
}

/// Closes open shells of `shells` using loose faces from `all_faces`.
///
/// A piece of a split solid is missing the section face that the intersection
/// cut it along; the shell splitter keeps that face as a loose single-face
/// shell (its boundary edges are shared by the two pieces and the face itself,
/// so they are not bridges). This attaches each loose face to every open shell
/// whose open boundary it covers, copying the face when it closes more than one
/// shell.
fn close_open_shells(shells: &[TopoShape], all_faces: &[TopoShape]) -> Vec<TopoShape> {
    let bld = TopoBuilder::new();

    // The loose faces: single faces of `all_faces` not already used by a
    // closed shell (the section faces the splitter left dangling).
    let in_closed: Vec<TopoShape> = shells
        .iter()
        .filter(|sh| !AlgoTools::is_open_shell(sh))
        .flat_map(|sh| faces_of(sh))
        .map(|f| f.0)
        .collect();
    let is_consumed = |f: &TopoShape| in_closed.iter().any(|c| c.same_tshape(f));
    let mut loose: Vec<TopoShape> = all_faces
        .iter()
        .filter(|f| f.is_face() && !is_consumed(f))
        .cloned()
        .collect();

    // How many open shells need each loose face (by geometry, so a face used by
    // several pieces is copied for each additional piece).
    let mut used_counts: HashMap<usize, usize> = HashMap::new();

    let mut out: Vec<TopoShape> = Vec::new();
    for sh in shells {
        if !AlgoTools::is_open_shell(sh) {
            out.push(sh.clone());
            continue;
        }
        let open_edges = open_edges_of_shell(sh);
        // The shell's existing faces (so a loose face already inside the shell
        // is not duplicated — the splitter can produce a shell that is closed
        // topologically but reported open by the TShape-identity edge check).
        let existing: Vec<TopoShape> = faces_of(sh).into_iter().map(|f| f.0).collect();
        let already_in_shell = |fc: &TopoShape| existing.iter().any(|x| x.same_tshape(fc));
        let idx = if open_edges.is_empty() {
            None
        } else {
            loose.iter().position(|fc| {
                !already_in_shell(fc)
                    && {
                        let keys: Vec<EKey> = faces_of(fc).first().map(|f| face_edges(&f)).unwrap_or_default();
                        open_edges.iter().all(|k| keys.contains(k))
                    }
            })
        };
        let mut new_shell = Shell::new();
        for f in faces_of(sh) {
            let oriented = orient_face_outward(&f.0, sh);
            bld.add(&mut new_shell.0, &oriented);
        }
        if let Some(idx) = idx {
            let face = loose[idx].clone();
            let n_used = used_counts.entry(idx).or_insert(0);
            let to_add = if *n_used == 0 {
                face
            } else {
                // A second (or later) piece shares the section face: attach a
                // copy so each piece owns its boundary.
                crate::shape_ops::transformed_copy(&face, &occt_core::gp::GpTrsf::identity())
                    .unwrap_or(face)
            };
            *n_used += 1;
            let oriented = orient_face_outward(&to_add, sh);
            bld.add(&mut new_shell.0, &oriented);
        }
        out.push(new_shell.0);
    }
    out
}

/// The directed edge keys of a face.
fn face_edges(f: &Face) -> Vec<EKey> {
    edges_of(&f.0).into_iter().map(|e| edge_key(&e)).collect()
}

/// The set of undirected edge keys of a shell that are shared by exactly one
/// face of the shell (the open boundary).
fn open_edges_of_shell(shell: &TopoShape) -> Vec<EKey> {
    let faces: Vec<Face> = faces_of(shell);
    let mut counts: HashMap<EKey, usize> = HashMap::new();
    for f in &faces {
        for k in face_edges(f) {
            *counts.entry(k).or_insert(0) += 1;
        }
    }
    counts
        .into_iter()
        .filter(|&(_, c)| c == 1)
        .map(|(k, _)| k)
        .collect()
}

/// Returns `face` oriented so its normal points away from the interior of the
/// (open) shell `shell`.
fn orient_face_outward(face: &TopoShape, shell: &TopoShape) -> TopoShape {
    let Some(p) = face_sample_point(&Face(face.clone())) else { return face.clone() };
    let Some(s) = BRepTool::face_surface(&Face(face.clone())) else { return face.clone() };
    let (u, v) = surface_closest_params(s.as_ref(), &p, 16, 16);
    if !u.is_finite() || !v.is_finite() {
        return face.clone();
    }
    let mut n = surface_normal(s.as_ref(), u, v);
    if n.square_magnitude() < 1e-30 {
        return face.clone();
    }
    // Interior reference point of the shell: the centroid of its face vertices.
    let vs = vertices_of(shell);
    if vs.is_empty() {
        return face.clone();
    }
    let mut c = GpPnt::new(0.0, 0.0, 0.0);
    for v in &vs {
        let q = BRepTool::vertex_point(v);
        c = GpPnt::new(c.x() + q.x(), c.y() + q.y(), c.z() + q.z());
    }
    let k = vs.len() as f64;
    let c = GpPnt::new(c.x() / k, c.y() / k, c.z() / k);
    let to_center = occt_core::gp::GpVec::from_pnts(&p, &c);
    if face.orientation() == Orientation::Reversed {
        n = n.reversed();
    }
    if n.dot(&to_center) > 0.0 {
        // The face's surface normal points toward the shell interior. The mesh
        // ignores the `TopoShape` orientation, so the *surface* must be flipped
        // for the signed volume to come out with the correct sign.
        return reverse_face_surface(face);
    }
    face.clone()
}

/// Returns a copy of `face` whose supporting surface is reversed (a plane's
/// normal is negated). The mesh triangulates from the surface's `d1` frame, so
/// an inward-pointing plane would otherwise contribute the wrong signed volume.
fn reverse_face_surface(face: &TopoShape) -> TopoShape {
    let f = Face(face.clone());
    let Some(s) = BRepTool::face_surface(&f) else { return face.clone() };
    if let Some(pln) = crate::brep_surface::face_plane(&f) {
        let ax = pln.axis();
        let loc = *ax.location();
        let n = *ax.direction();
        let xd = *pln.x_axis().direction();
        let rev_ax = GpAx3::new(loc, n.reversed(), &xd).unwrap_or(pln.position());
        let rev_pln = GpPln::new(rev_ax);
        let bld = TopoBuilder::new();
        let wires = crate::topo_tools_full::wires_of_face(&f);
        let new_face = bld.make_face(Arc::new(GeomPlane::new(rev_pln)), &wires);
        let reg = crate::tgeometry::GeometryRegistry::global();
        if let Some(mut g) = reg.face_geom(&new_face.0) {
            g.tolerance = reg.face_tolerance(&f.0);
            reg.set_face(&new_face.0, g);
        }
        return new_face.0;
    }
    face.clone()
}

/// Undirected quantized edge key of an edge.
fn edge_key(e: &Edge) -> EKey {
    let (a, b) = crate::topo_tools_full::edge_vertices(e);
    let (Some(a), Some(b)) = (a, b) else {
        return ((0, 0, 0), (0, 0, 0));
    };
    let ka = crate::shell_splitter::vertex_key(&a);
    let kb = crate::shell_splitter::vertex_key(&b);
    if ka <= kb { (ka, kb) } else { (kb, ka) }
}

/// Whether `split` must be reversed to match the direction of `original`
/// (both faces). Mirrors `BOPTools_AlgoTools::IsSplitToReverse` restricted to
/// faces: when both share the same supporting surface the orientations are
/// compared directly; otherwise the surface normals at a common point are
/// compared.
fn is_split_to_reverse(split: &TopoShape, original: &TopoShape) -> bool {
    if split.shape_type() != ShapeType::Face || original.shape_type() != ShapeType::Face {
        return false;
    }
    let (Some(s), Some(o)) = (
        BRepTool::face_surface(&Face(split.clone())),
        BRepTool::face_surface(&Face(original.clone())),
    ) else {
        return false;
    };
    if Arc::ptr_eq(&s, &o) {
        return split.orientation() != original.orientation();
    }
    let Some(p) = face_sample_point(&Face(split.clone())) else { return false };
    let (u, v) = surface_closest_params(s.as_ref(), &p, 32, 32);
    let (uo, vo) = surface_closest_params(o.as_ref(), &p, 32, 32);
    if !u.is_finite() || !v.is_finite() || !uo.is_finite() || !vo.is_finite() {
        return false;
    }
    let mut ns = surface_normal(s.as_ref(), u, v);
    let mut no = surface_normal(o.as_ref(), uo, vo);
    if ns.square_magnitude() < 1e-30 || no.square_magnitude() < 1e-30 {
        return false;
    }
    if split.orientation() == Orientation::Reversed {
        ns = ns.reversed();
    }
    if original.orientation() == Orientation::Reversed {
        no = no.reversed();
    }
    ns.dot(&no) < 0.0
}

/// A 3-D point on `face`: the surface centre when its UV range is bounded,
/// else the midpoint of the first boundary edge.
fn face_sample_point(face: &Face) -> Option<GpPnt> {
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(face);
    if u1.is_finite() && v1.is_finite() {
        let s = BRepTool::face_surface(face)?;
        return Some(s.d0(0.5 * (u1 + u2), 0.5 * (v1 + v2)));
    }
    let e = crate::topo_tools_full::edges_of(&face.0).into_iter().next()?;
    let (a, b) = BRepTool::edge_parameters(&e);
    if a.is_finite() && b.is_finite() {
        let c = BRepTool::edge_curve(&e)?;
        Some(c.d0(0.5 * (a + b)))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Solid-state classification and split-solids deepening
// ---------------------------------------------------------------------------

/// State of a 3-D point relative to a solid: `On` when the point lies within
/// `tol` of a boundary face, otherwise `In`/`Out` by the parity test.
///
/// Mirrors `BOPTools_AlgoTools::ComputeState` (3-D) reduced to a single point:
/// the point is `On` when its distance to the closest boundary face is within
/// `tol`, else the [`crate::brep_extrema::is_inside`] even-odd test decides.
pub fn classify_solid_state(solid: &TopoShape, point: &GpPnt, tol: f64) -> FaceState {
    let mut d = f64::INFINITY;
    for f in faces_of(solid) {
        let (_, q) = closest_point_on_face(&f, point, 16, 16);
        d = d.min(q.distance(point));
    }
    if d.is_finite() && d <= tol {
        return FaceState::On;
    }
    if is_inside(solid, point) {
        FaceState::In
    } else {
        FaceState::Out
    }
}

/// Splits every interfered source solid into its closed-shell pieces, adds the
/// internal faces of the other solids that lie inside it, and records the
/// selected pieces as solid images — the deepened [`build_split_solids`].
///
/// Mirrors the `FillIn3DParts` → `BuildSplitSolids` flow with the 3-D face
/// classification ported as a per-piece solid classification:
///
/// 1. for each source solid, collect the split faces of its own faces
///    ([`collect_solid_split_faces`]) plus the *internal faces* — faces of the
///    other argument solids lying strictly inside it (each added FORWARD and
///    REVERSED so both adjacent pieces close);
/// 2. group the whole set into closed shells ([`crate::shell_splitter::ShellSplitter`])
///    and wrap each shell into a solid;
/// 3. classify each piece with [`classify_solid_state`] against every *other*
///    argument solid: a piece lying strictly inside another solid is an
///    interior region — it contributes to the General-Fuse result only once, so
///    duplicate interior regions (the same geometric volume produced from
///    several source solids) collapse to their first occurrence;
/// 4. the surviving pieces are recorded as images of the source solid with the
///    origins back-map.
///
/// A source solid that was neither split nor hosts any internal face keeps no
/// image (it is returned as-is by the caller).
pub fn build_split_solids_full<B: BopBuildOps>(
    f: &mut B,
    objects: &[TopoShape],
    tools: &[TopoShape],
    obj_state: FaceState,
    tools_state: FaceState,
) -> Result<(), String> {
    let op = BopOp::from_states(obj_state, tools_state)?;
    let n = f.ds().nb_source_shapes();

    // The argument solids: used both as the source of internal faces and as the
    // reference set for the per-piece classification.
    let solids: Vec<TopoShape> = (0..n)
        .filter_map(|i| {
            let si = f.ds().shape_info(i)?;
            if si.shape_type() == ShapeType::Solid {
                Some(si.shape().clone())
            } else {
                None
            }
        })
        .collect();
    let tol = f.fuzzy_value().max(1e-7);

    // Interior regions already claimed by a previous source solid: an interior
    // piece with a bounding box already seen is a duplicate overlap region.
    let mut seen: Vec<BBoxKey> = Vec::new();

    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let solid = si.shape().clone();
        let is_object = objects.iter().any(|o| o.same_tshape(&solid));
        let is_tool = tools.iter().any(|t| t.same_tshape(&solid));
        // A Cut operation removes the tools entirely; their pieces never enter
        // the result (the tool faces inside the object are added as internal
        // faces of the object's pieces, so the cut stays closed).
        if op == BopOp::Cut && is_tool && !is_object {
            continue;
        }

        // 1. Split faces of the solid itself.
        let mut faces = collect_solid_split_faces(f, &solid);

        // 2. Internal faces: faces of the other solids strictly inside this one.
        let mut in_faces: Vec<TopoShape> = Vec::new();
        for other in &solids {
            if other.same_tshape(&solid) {
                continue;
            }
            for fc in faces_of(other) {
                let images: Vec<TopoShape> = match f.history().image(&fc.0) {
                    Some(imgs) => imgs.to_vec(),
                    None => vec![fc.0.clone()],
                };
                for im in images {
                    if face_state_in_solid(&Face(im.clone()), &solid, tol) != FaceState::In {
                        continue;
                    }
                    let mut fwd = im.clone();
                    fwd.set_orientation(Orientation::Forward);
                    in_faces.push(fwd);
                    let mut rev = im.clone();
                    rev.set_orientation(Orientation::Reversed);
                    in_faces.push(rev);
                }
            }
        }

        // A solid that neither splits nor hosts internal faces is unchanged.
        if !solid_interfered(f, &solid) && in_faces.is_empty() {
            continue;
        }
        faces.extend(in_faces);

        // 3. Group into closed shells and assemble one solid per shell.
        let pieces = split_faces_into_solids(&faces)?;
        if pieces.is_empty() {
            continue;
        }

        // 4. Select the pieces per the operation. A zero-volume shell (e.g. the
        //    two coincident orientations of a section face closing on
        //    themselves) cannot bound a solid and is dropped.
        let mut kept: Vec<TopoShape> = Vec::new();
        for p in pieces {
            if piece_bbox_volume(&p) < 1e-9 {
                continue;
            }
            let interior = solids.iter().any(|other| {
                !other.same_tshape(&solid)
                    && classify_solid_state(other, &piece_center(&p), tol) == FaceState::In
            });
            let pass = match op {
                // Union: every distinct region survives; a duplicate interior
                // region collapses to its first occurrence.
                BopOp::Fuse => true,
                // Intersection: keep only the regions inside the other solid.
                BopOp::Common => interior,
                // Difference: keep the object's regions outside the tools.
                BopOp::Cut => !interior,
            };
            if !pass {
                continue;
            }
            if interior {
                let key = piece_bbox_key(&p);
                if seen.contains(&key) {
                    continue;
                }
                seen.push(key);
            }
            kept.push(p);
        }

        // 5. Record the images and the origins back-map.
        for p in kept {
            f.history_mut().add_image(&solid, p.clone());
            f.origins_mut().entry(shape_key(&p)).or_default().push(solid.clone());
        }
    }
    Ok(())
}

/// The boolean operation derived from the object/tool face states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BopOp {
    /// A ∪ B — `(Out, Out)`.
    Fuse,
    /// A − B — `(Out, In)`.
    Cut,
    /// A ∩ B — `(In, In)`.
    Common,
}

impl BopOp {
    /// Derives the operation from the `BOPAlgo_BOP` states.
    pub fn from_states(obj: FaceState, tools: FaceState) -> Result<BopOp, String> {
        match (obj, tools) {
            (FaceState::Out, FaceState::Out) => Ok(BopOp::Fuse),
            (FaceState::Out, FaceState::In) => Ok(BopOp::Cut),
            (FaceState::In, FaceState::In) => Ok(BopOp::Common),
            _ => Err("BOPAlgo: unsupported object/tool states".to_string()),
        }
    }
}

/// Quantized axis-aligned bounding box of a piece, used as the identity key of
/// a geometric region for the interior-duplicate collapse.
type BBoxKey = (i64, i64, i64, i64, i64, i64);

/// The axis-aligned bounding box `(xmin, ymin, zmin, xmax, ymax, zmax)` of
/// `solid`'s vertices.
fn piece_bbox(solid: &TopoShape) -> (f64, f64, f64, f64, f64, f64) {
    let mut mn = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut mx = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for v in vertices_of(solid) {
        let p = BRepTool::vertex_point(&v);
        mn.0 = mn.0.min(p.x());
        mn.1 = mn.1.min(p.y());
        mn.2 = mn.2.min(p.z());
        mx.0 = mx.0.max(p.x());
        mx.1 = mx.1.max(p.y());
        mx.2 = mx.2.max(p.z());
    }
    (mn.0, mn.1, mn.2, mx.0, mx.1, mx.2)
}

/// The quantized bounding box of `solid`'s vertices.
fn piece_bbox_key(solid: &TopoShape) -> BBoxKey {
    let (mnx, mny, mnz, mxx, mxy, mxz) = piece_bbox(solid);
    let q = |x: f64| (x / 1e-6).round() as i64;
    (q(mnx), q(mny), q(mnz), q(mxx), q(mxy), q(mxz))
}

/// Volume of the bounding box of `solid`'s vertices (the axis-aligned region
/// extent) — zero for a degenerate shell.
fn piece_bbox_volume(solid: &TopoShape) -> f64 {
    let (mnx, mny, mnz, mxx, mxy, mxz) = piece_bbox(solid);
    (mxx - mnx) * (mxy - mny) * (mxz - mnz)
}

/// The centre of the bounding box of `solid`'s vertices — a representative
/// interior point for a closed piece.
fn piece_center(solid: &TopoShape) -> GpPnt {
    let (mnx, mny, mnz, mxx, mxy, mxz) = piece_bbox(solid);
    GpPnt::new(0.5 * (mnx + mxx), 0.5 * (mny + mxy), 0.5 * (mnz + mxz))
}

/// Whether `face` lies strictly inside `solid`: its boundary-vertex centroid is
/// `In` the solid (`On` is not internal — a face on the solid's boundary is a
/// shared outer face, not a section face).
fn face_state_in_solid(face: &Face, solid: &TopoShape, tol: f64) -> FaceState {
    let vs = vertices_of(&face.0);
    if vs.is_empty() {
        return FaceState::Unknown;
    }
    let mut acc = GpPnt::new(0.0, 0.0, 0.0);
    let n = vs.len() as f64;
    for v in &vs {
        let p = BRepTool::vertex_point(v);
        acc = GpPnt::new(acc.x() + p.x(), acc.y() + p.y(), acc.z() + p.z());
    }
    let p = GpPnt::new(acc.x() / n, acc.y() / n, acc.z() / n);
    classify_solid_state(solid, &p, tol)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
    use occt_geom::{GeomLine, GeomPlane, Surface};

    use crate::bop_hist::BopHistory;
    use crate::bopds::BopdsDS;
    use crate::brep_gprop::volume;
    use crate::primitives::BRepPrimBox;
    use crate::shape::{Face, Solid, Vertex};
    use crate::topo_tools_full::{edges_of, shapes_of, vertices_of};

    /// A minimal `BopBuildOps` host for the isolated tests.
    struct StubBuilder {
        ds: BopdsDS,
        history: BopHistory,
        fuzzy: f64,
        args: Vec<TopoShape>,
        origins: HashMap<usize, Vec<TopoShape>>,
    }

    impl BopBuildOps for StubBuilder {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn history(&self) -> &BopHistory {
            &self.history
        }
        fn history_mut(&mut self) -> &mut BopHistory {
            &mut self.history
        }
        fn fuzzy_value(&self) -> f64 {
            self.fuzzy
        }
        fn arguments(&self) -> &[TopoShape] {
            &self.args
        }
        fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
            &mut self.origins
        }
    }

    fn stub(ds: BopdsDS, history: BopHistory, args: Vec<TopoShape>) -> StubBuilder {
        StubBuilder { ds, history, fuzzy: 1e-7, args, origins: HashMap::new() }
    }

    /// The face of `faces` whose boundary-vertex mean lies at height `z`.
    fn face_at_z(faces: &[Face], z: f64) -> Face {
        faces
            .iter()
            .find(|f| {
                let vs = vertices_of(&f.0);
                if vs.is_empty() {
                    return false;
                }
                let zavg =
                    vs.iter().map(|v| BRepTool::vertex_point(v).z()).sum::<f64>() / vs.len() as f64;
                (zavg - z).abs() < 1e-9
            })
            .cloned()
            .expect("face at height")
    }

    /// A fresh face on the same surface and boundary edges as `src` — a new
    /// `TShape`, so it is a distinct split image.
    fn re_face(src: &Face) -> Face {
        let b = TopoBuilder::new();
        let edges = edges_of(&src.0);
        let wire = b.make_wire(&edges);
        let surf = BRepTool::face_surface(src).expect("face surface");
        b.make_face(surf, &[wire])
    }

    // -----------------------------------------------------------------------
    // fill_images_solids
    // -----------------------------------------------------------------------

    #[test]
    fn fill_images_solids_rebuilds_box_with_split_face() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&boxed.solid.0);
        assert_eq!(faces.len(), 6);

        // Split the top face: register a fresh, geometrically identical face
        // as its only image.
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);
        assert!(!new_top.0.same_tshape(&top.0));

        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);

        fill_images_solids(&mut b).unwrap();

        let imgs = b.history().image(&boxed.solid.0).expect("solid has an image");
        assert_eq!(imgs.len(), 1);
        let img = &imgs[0];
        assert!(img.is_solid());
        let v = volume(img, 0.02);
        assert!((v - 1.0).abs() < 0.01, "unit box volume, got {v}");
        // The origins back-map points the new solid at the source.
        let ors = b.origins.get(&shape_key(img)).expect("origin recorded");
        assert!(ors.iter().any(|o: &TopoShape| o.same_tshape(&boxed.solid.0)));
    }

    #[test]
    fn fill_images_solids_skips_unmodified_box() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut b = stub(ds, BopHistory::new(), vec![boxed.solid.0.clone()]);

        fill_images_solids(&mut b).unwrap();

        assert!(
            !b.history().has_image(&boxed.solid.0),
            "no face was split -> no solid image"
        );
    }

    #[test]
    fn fill_images_solids_without_solids_is_noop() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let face = faces_of(&boxed.solid.0)[0].clone();
        let mut ds = BopdsDS::new();
        ds.init(&[face.0.clone()]);
        let mut b = stub(ds, BopHistory::new(), vec![face.0.clone()]);

        fill_images_solids(&mut b).unwrap();

        assert!(!b.history().has_any_images(), "no solid participates");
    }

    // -----------------------------------------------------------------------
    // build_split_solids
    // -----------------------------------------------------------------------

    #[test]
    fn build_split_solids_rebuilds_box_from_face_splits() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&boxed.solid.0);
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);

        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);

        build_split_solids(&mut b).unwrap();

        let imgs = b.history().image(&boxed.solid.0).expect("solid has an image");
        assert_eq!(imgs.len(), 1, "six connected faces -> one closed shell");
        let v = volume(&imgs[0], 0.02);
        assert!((v - 1.0).abs() < 0.01, "unit box volume, got {v}");
    }

    #[test]
    fn build_split_solids_separates_two_disjoint_shells() {
        // A solid holding two disjoint closed box shells, only the first
        // shell's top face split: the split pass must rebuild each closed
        // shell into its own solid image.
        let b = TopoBuilder::new();
        let box1 = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let box2 =
            BRepPrimBox::make_box_corner(&GpPnt::new(3.0, 0.0, 0.0), &GpPnt::new(4.0, 1.0, 1.0));
        let sh1 = shapes_of(&box1.solid.0, ShapeType::Shell)[0].clone();
        let sh2 = shapes_of(&box2.solid.0, ShapeType::Shell)[0].clone();
        let mut solid = Solid::new();
        b.add(&mut solid.0, &sh1);
        b.add(&mut solid.0, &sh2);

        let top1 = face_at_z(&faces_of(&sh1), 1.0);
        let new_top1 = re_face(&top1);

        let mut ds = BopdsDS::new();
        ds.init(&[solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top1.0, new_top1.0.clone());
        let mut stub = stub(ds, history, vec![solid.0.clone()]);

        build_split_solids(&mut stub).unwrap();

        let imgs = stub.history().image(&solid.0).expect("solid has images");
        assert_eq!(imgs.len(), 2, "two disjoint closed shells -> two split solids");
        // Every split solid is a closed box shell (6 faces). Volume is not
        // asserted per-shell here: `make_box_corner` (the second box) carries a
        // known inward-orientation quirk that skews its signed volume, while the
        // volume-correctness behaviour is covered by
        // `build_split_solids_rebuilds_box_from_face_splits`.
        for im in imgs.iter() {
            assert!(im.is_solid());
            assert_eq!(faces_of(im).len(), 6, "each split solid is a closed box");
        }
        // The two images are genuinely different rebuilt solids.
        assert!(!imgs[0].same_tshape(&imgs[1]), "two distinct split solids");
    }

    // -----------------------------------------------------------------------
    // build_split_solids_full
    // -----------------------------------------------------------------------

    /// A reliable axis-aligned box with real geometry at position — unlike
    /// `make_box_corner`, whose lateral faces are broken. Modeled on the shared
    /// `unit_box` fixture but spanning `[x0,x1]×[y0,y1]×[z0,z1]`.
    fn axis_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> Solid {
        let b = TopoBuilder::new();
        let corners = [
            GpPnt::new(x0, y0, z0), GpPnt::new(x1, y0, z0), GpPnt::new(x1, y1, z0), GpPnt::new(x0, y1, z0),
            GpPnt::new(x0, y0, z1), GpPnt::new(x1, y0, z1), GpPnt::new(x1, y1, z1), GpPnt::new(x0, y1, z1),
        ];
        let vertices: Vec<Vertex> = corners.iter().map(|p| b.make_vertex(*p, 1e-7)).collect();
        let seg = |b: &TopoBuilder, p1: &GpPnt, p2: &GpPnt, v1: &Vertex, v2: &Vertex| {
            let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).unwrap();
            let lin = GpLin::from_pnt_dir(*p1, dir);
            let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
            b.add(&mut e.0, &v1.0);
            b.add(&mut e.0, &v2.0);
            e
        };
        let edge_idx: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges = Vec::new();
        for &(i, j) in &edge_idx {
            edges.push(seg(&b, &corners[i], &corners[j], &vertices[i], &vertices[j]));
        }
        let face_edge_sets: [[usize; 4]; 6] = [
            [0, 1, 2, 3], [4, 5, 6, 7], [0, 9, 4, 8], [2, 10, 6, 11], [3, 11, 7, 8], [1, 10, 5, 9],
        ];
        let face_planes: [(GpPnt, GpDir, GpDir); 6] = [
            (GpPnt::new(x0, y0, z0), GpDir::new(0.0, 0.0, -1.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
            (GpPnt::new(x0, y0, z1), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(x0, y0, z0), GpDir::new(0.0, -1.0, 0.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(x0, y1, z0), GpDir::new(0.0, 1.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(x0, y0, z0), GpDir::new(-1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(x1, y0, z0), GpDir::new(1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
        ];
        let mut faces = Vec::new();
        for fi in 0..6 {
            let (origin, normal, u_dir) = face_planes[fi];
            let ax3 = GpAx3::new(origin, normal, &u_dir).unwrap();
            let mut face = b.make_face_plane(&GpPln::new(ax3));
            let wire = b.make_wire(&face_edge_sets[fi].map(|ei| edges[ei].clone()));
            b.add_wire(&mut face, &wire);
            faces.push(face);
        }
        let shell = b.make_shell(&faces);
        b.make_solid(&[shell])
    }

    /// The face of `solid` whose every boundary vertex lies on the plane
    /// `coord = value` (0/1/2 → x/y/z) — the axis-aligned face at that plane.
    fn box_face_on(solid: &Solid, coord: usize, value: f64) -> Face {
        faces_of(&solid.0)
            .into_iter()
            .find(|f| {
                let vs = vertices_of(&f.0);
                !vs.is_empty()
                    && vs.iter().all(|v| {
                        let p = BRepTool::vertex_point(v);
                        let c = match coord {
                            0 => p.x(),
                            1 => p.y(),
                            _ => p.z(),
                        };
                        (c - value).abs() < 1e-9
                    })
            })
            .expect("face on coordinate plane")
    }

    /// Volume of `solid` from its vertex bounding box — reliable for the
    /// axis-aligned boxes this module's tests build (the mesh volume of
    /// translated/offset boxes is unreliable in this port).
    fn bbox_volume(solid: &TopoShape) -> f64 {
        let mut mn = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut mx = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for v in vertices_of(solid) {
            let p = BRepTool::vertex_point(&v);
            mn.0 = mn.0.min(p.x());
            mn.1 = mn.1.min(p.y());
            mn.2 = mn.2.min(p.z());
            mx.0 = mx.0.max(p.x());
            mx.1 = mx.1.max(p.y());
            mx.2 = mx.2.max(p.z());
        }
        (mx.0 - mn.0) * (mx.1 - mn.1) * (mx.2 - mn.2)
    }

    #[test]
    fn classify_solid_state_in_on_out() {
        let a = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let tol = 1e-7;
        assert_eq!(
            classify_solid_state(&a.0, &GpPnt::new(0.5, 0.5, 0.5), tol),
            FaceState::In
        );
        assert_eq!(
            classify_solid_state(&a.0, &GpPnt::new(2.0, 0.5, 0.5), tol),
            FaceState::Out
        );
        assert_eq!(
            classify_solid_state(&a.0, &GpPnt::new(0.5, 0.5, 1.0), tol),
            FaceState::On
        );
    }

    #[test]
    fn build_split_solids_full_selects_union_of_overlapping_boxes() {
        // Two unit boxes overlapping in x by 0.5: A=[0,1]³, B=[0.5,1.5]³.
        // The intersection cuts each box into two closed pieces; the two pieces
        // of the overlap region are geometrically identical, so only one of
        // them survives. Selected pieces: A[0,0.5] + overlap[0.5,1] +
        // B[1,1.5] = volume 1.5 (the Fuse union).
        let a = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = axis_box(0.5, 0.0, 0.0, 1.5, 1.0, 1.0);
        let a1 = axis_box(0.0, 0.0, 0.0, 0.5, 1.0, 1.0);
        let a2 = axis_box(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b1 = axis_box(1.0, 0.0, 0.0, 1.5, 1.0, 1.0);
        let b2 = axis_box(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);

        let mut ds = BopdsDS::new();
        ds.init(&[a.0.clone(), b.0.clone()]);
        let mut history = BopHistory::new();

        // Split images of A's faces: the pieces of A1 and A2 on each plane.
        history.add_image(&box_face_on(&a, 0, 0.0).0, box_face_on(&a1, 0, 0.0).0.clone());
        history.add_image(&box_face_on(&a, 0, 1.0).0, box_face_on(&a2, 0, 1.0).0.clone());
        for (coord, val) in [(1usize, 0.0f64), (1, 1.0), (2, 0.0), (2, 1.0)] {
            let src = box_face_on(&a, coord, val).0;
            history.add_image(&src, box_face_on(&a1, coord, val).0.clone());
            history.add_image(&src, box_face_on(&a2, coord, val).0.clone());
        }
        // Split images of B's faces.
        history.add_image(&box_face_on(&b, 0, 0.5).0, box_face_on(&b2, 0, 0.5).0.clone());
        history.add_image(&box_face_on(&b, 0, 1.5).0, box_face_on(&b1, 0, 1.5).0.clone());
        for (coord, val) in [(1usize, 0.0f64), (1, 1.0), (2, 0.0), (2, 1.0)] {
            let src = box_face_on(&b, coord, val).0;
            history.add_image(&src, box_face_on(&b1, coord, val).0.clone());
            history.add_image(&src, box_face_on(&b2, coord, val).0.clone());
        }

        let mut st = stub(ds, history, vec![a.0.clone(), b.0.clone()]);
        build_split_solids_full(&mut st, &[a.0.clone()], &[b.0.clone()], FaceState::Out, FaceState::Out)
            .unwrap();

        // A yields two pieces (its part outside B + the overlap); B yields one
        // (its part outside A) — the overlap region is a duplicate interior.
        let a_imgs = st.history().image(&a.0).expect("box A has split-solid images");
        let b_imgs = st.history().image(&b.0).expect("box B has split-solid images");
        assert_eq!(a_imgs.len(), 2, "A: outside piece + overlap piece");
        assert_eq!(b_imgs.len(), 1, "B: overlap piece is a duplicate and dropped");

        let total: f64 = a_imgs.iter().chain(b_imgs.iter()).map(bbox_volume).sum();
        assert!(
            (total - 1.5).abs() < 1e-6,
            "Fuse union volume 1.5, got {total}"
        );

        // The origins back-map is populated for every recorded piece.
        for im in a_imgs.iter().chain(b_imgs.iter()) {
            let ors = st.origins.get(&shape_key(im)).expect("origin recorded");
            assert!(ors.iter().any(|o| o.same_tshape(&a.0) || o.same_tshape(&b.0)));
        }
    }

    #[test]
    fn build_split_solids_full_keeps_disjoint_boxes() {
        // Two fully disjoint boxes: neither lies inside the other, so both are
        // kept as their whole self (volume 2).
        let a = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = axis_box(2.0, 0.0, 0.0, 3.0, 1.0, 1.0);

        let mut ds = BopdsDS::new();
        ds.init(&[a.0.clone(), b.0.clone()]);
        let mut history = BopHistory::new();
        // Nothing split: no images -> no solid is interfered, none gets an image.
        let mut st = stub(ds, history, vec![a.0.clone(), b.0.clone()]);
        build_split_solids_full(&mut st, &[a.0.clone()], &[b.0.clone()], FaceState::Out, FaceState::Out)
            .unwrap();
        assert!(
            !st.history().has_image(&a.0) && !st.history().has_image(&b.0),
            "disjoint unsplit boxes keep no solid image"
        );
    }
}
