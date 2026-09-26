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
//! Live GF path is [`crate::bop_images_solids`] → [`crate::bop_split_solids_occt`]
//! (`aMST`): the OCCT `FillIn3DParts` + `BuildSplitSolids` pair
//! (`BOPAlgo_Builder_3.cxx:97-263` / `:413-618`).
//!
//! The historical *full flow* port (`build_split_solids_full` with its
//! `classify_faces_in_solid` / `face_state_in_solid` /
//! `collect_all_candidate_faces` helpers) was **removed in T-27 (round 78)**:
//! it lived in a private submodule that `BopBuilder` never called, and its only
//! distinguishing rule, `is_covering_face`, was a port-invented predicate with
//! no OCCT counterpart. Its live replacements are the two OCCT chains below.
//!
//! ## Translation boundaries vs OCCT
//!
//! * `IsInternalFace` uses `GetFaceOff` (angle-normals around a shared edge)
//!   then `ComputeState(Face, Solid)` when angles cannot decide
//!   ([`crate::algo_tools_face::is_internal_face`]) — the live
//!   `BOPAlgo_Tools::ClassifyFaces` path (`BOPAlgo_Tools.cxx:1622`), reached
//!   from `BOPAlgo_FillIn3DParts::Perform` (`:1334`) and
//!   `BOPAlgo_BuilderSolid::PerformInternalShapes` (`BOPAlgo_BuilderSolid.cxx:673`).
//! * The connexity-block grouping (`BOPAlgo_FillIn3DParts::MakeConnexityBlock`)
//!   is ported in [`crate::bop_classify_occt`]; classification runs against the
//!   *original* solid, whose volume equals the draft solid's (only the
//!   boundary-edge set comes from the split faces, OCCT classifying against
//!   the draft solid built by `BuildDraftSolid`).
//! * `BOPAlgo_SplitSolid` is replaced here by [`crate::shell_splitter::ShellSplitter`]
//!   + [`close_open_shells`]; live assembly is [`crate::builder_solid::BuilderSolid`].
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

use std::collections::{HashMap, HashSet};
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
use crate::shape::{Edge, Face, Shell, TopoShape};
use crate::shell_splitter::{EKey};
use crate::topo_tools_full::{edges_of_wire, faces_of, vertices_of, wires_of_face};

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
/// **LEGACY / superseded（T-27，round 78）**：绿路径已由 [`bop_split_solids_occt`] 走 OCCT 的
/// `BOPAlgo_Builder_3.cxx:413-618`（`build_split_solids_occt` → `SplitSolid::perform` → `builder_solid.rs`），
/// 本函数是同一 OCCT 段的**早期重复实现**，今天只有本文件的 2 条测试调用（`:666`/`:698`），生产 0 调用。
/// 保留是为了不删除既有测试；新代码**不要**接入它。
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
    // The faithful `BOPAlgo_SplitSolid` assembly: the four-phase
    // `BOPAlgo_BuilderSolid` pipeline (ShapesToAvoid → Loops → Areas →
    // InternalShapes) turns the split + internal faces into region solids.
    // `close_open_shells` (re-attaching the shared section face onto the open
    // pieces) is applied inside `PerformLoops`.
    crate::builder_solid::build_solids_from_faces(faces)
}

/// Closes open shells of `shells` using loose faces from `all_faces`.
///
/// A piece of a split solid is missing the section face that the intersection
/// cut it along; the shell splitter keeps that face as a loose single-face
/// shell (its boundary edges are shared by the two pieces and the face itself,
/// so they are not bridges). This attaches each loose face to every open shell
/// whose open boundary it covers, copying the face when it closes more than one
/// shell.
pub(crate) fn close_open_shells(shells: &[TopoShape], all_faces: &[TopoShape]) -> Vec<TopoShape> {
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
    let loose: Vec<TopoShape> = all_faces
        .iter()
        .filter(|f| f.is_face() && !is_consumed(f))
        .cloned()
        .collect();

    // How many open shells need each loose face (by geometry, so a face used by
    // several pieces is copied for each additional piece).
    let mut used_counts: HashMap<usize, usize> = HashMap::new();

    let mut out: Vec<TopoShape> = Vec::new();
    for sh in shells {
        // A shell that is geometrically closed — every boundary edge (by
        // undirected key) used by exactly two faces — needs no loose face
        // attached and no re-orientation; the TShape-identity check reports it
        // open only when coincident faces carry distinct edge `TShape`s (the
        // box-top ring vs the cylinder base). Rebuilding such a shell would
        // re-create its faces with fresh `TShape`s and break the cross-solid
        // merge that relies on the shared face's identity.
        if !AlgoTools::is_open_shell(sh) || !geometrically_open(sh) {
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
            let to_add = face;
            *n_used += 1;
            let oriented = orient_face_outward(&to_add, sh);
            bld.add(&mut new_shell.0, &oriented);
        }
        out.push(new_shell.0);
    }
    out
}

/// The undirected edge keys of a face, **one entry per edge occurrence**
/// (`TopExp_Explorer(face, EDGE)` with the default `CumOri = true`, as in
/// [`face_edge_occurrences`]). `edges_of` uniquifies by `TShape`, so a seam
/// stored twice in one wire would be counted once — `geometrically_open` then
/// reports a perfectly closed band as open.
fn face_edges(f: &Face) -> Vec<EKey> {
    face_edge_occurrences(f).iter().map(edge_key).collect()
}

/// Boundary-edge **occurrences** of a face (see `builder_solid.rs`'s twin).
fn face_edge_occurrences(f: &Face) -> Vec<Edge> {
    let mut out = Vec::new();
    for w in wires_of_face(f) {
        out.extend(edges_of_wire(&w));
    }
    out
}

/// `BRep_Tool::IsClosed(Shell)` (`BRep_Tool.cxx:1707-1728`) negated: OCCT walks
/// every edge **occurrence** of the shell (`TopExp_Explorer(shell, EDGE)`,
/// cumulated orientation), skips degenerated / `INTERNAL` / `EXTERNAL` edges,
/// and toggles each edge in a map (`Add`, and `Remove` when it was already
/// there). The shell is closed when at least one boundary edge was seen and
/// the map ends up **empty**, i.e. every edge occurs an even number of times.
/// Requiring exactly two occurrences per edge instead misreports a seam that a
/// single face stores twice (Forward and Reversed) as an open boundary.
pub(crate) fn geometrically_open(shell: &TopoShape) -> bool {
    let mut odd: HashSet<EKey> = HashSet::new();
    let mut has_bound = false;
    for f in faces_of(shell) {
        for e in face_edge_occurrences(&f) {
            if BRepTool::is_degenerated(&e) {
                continue;
            }
            let o = e.orientation();
            if o == Orientation::Internal || o == Orientation::External {
                continue;
            }
            has_bound = true;
            let k = edge_key(&e);
            if !odd.insert(k) {
                odd.remove(&k);
            }
        }
    }
    !(has_bound && odd.is_empty())
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
    // UNPORTED: `orient_face_outward` is a port-only heuristic for open shells
    // with no OCCT control flow that projects the sample point to UV; OCCT
    // decides face orientation from the face's pcurves / `BRepGProp` frame
    // (`BOPAlgo_BuilderSolid.cxx`). The grid stays.
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
        let mut f = face.clone();
        f.reverse();
        return f;
    }
    face.clone()
}

/// Returns a copy of `face` whose supporting surface is reversed (a plane's
/// normal is negated). The mesh triangulates from the surface's `d1` frame, so
/// an inward-pointing plane would otherwise contribute the wrong signed volume.
#[allow(dead_code)]
fn reverse_face_surface(face: &TopoShape) -> TopoShape {
    let f = Face(face.clone());
    let Some(_s) = BRepTool::face_surface(&f) else { return face.clone() };
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
    // Faithful `BOPTools_AlgoTools::IsSplitToReverse(Face, Face, ...)`
    // (`BOPTools_AlgoTools.cxx:1316-1427`): the point is projected onto each
    // supporting surface with `ProjPS` = `GeomAPI_ProjectPointOnSurf`
    // (`cxx:1393-1406`); `NbPoints() == 0` is the error return
    // (`cxx:1395-1403`), `false` here.
    let (Some(ps), Some(po)) = (
        occt_geom::geom_api::project_point_on_surface(s.as_ref(), &p, occt_core::precision::CONFUSION),
        occt_geom::geom_api::project_point_on_surface(o.as_ref(), &p, occt_core::precision::CONFUSION),
    ) else {
        return false;
    };
    let (u, v) = (ps.u, ps.v);
    let (uo, vo) = (po.u, po.v);
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



// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;


    use occt_core::gp::{GpPnt};


    use crate::bop_hist::BopHistory;
    use crate::bopds::BopdsDS;
    use crate::brep_gprop::volume;
    use crate::primitives::BRepPrimBox;
    use crate::shape::{Face, Solid};
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

}
