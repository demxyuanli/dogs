//! Face-image reconstruction — Phase 20 wave C2b-2a.
//!
//! Port of `BOPAlgo_Builder::FillImagesFaces` (with its two stages
//! `BuildSplitFaces` and `FillSameDomainFaces`) from
//! `BOPAlgo_Builder_2.cxx`:
//!
//! | OCCT method               | Rust function                         |
//! |---------------------------|---------------------------------------|
//! | `FillImagesFaces`         | [`fill_images_faces`]                 |
//! | `BuildSplitFaces`         | [`build_split_faces`]                 |
//! | `FillSameDomainFaces`     | [`fill_same_domain_faces`]            |
//!
//! The three entry points are generic over a [`BopBuilderLike`] host — the
//! minimal contract `crate::bop_builder2::BopBuilder` (General Fuse builder,
//! Phase 20) satisfies. The host exposes the BOPDS, the images history, the
//! same-domain map and the option flags the OCCT builder stores as members
//! (`myDS`, `myImages`, `myShapesSD`, `myPaveFiller`, …). Keeping the logic
//! decoupled from the concrete type follows the [`crate::pave_blocks::PaveFillerLike`]
//! / [`crate::bop_build_common::BopBuildOps`] pattern: the module is verifiable
//! in isolation with a local test stub, and the builder main class adapts by
//! implementing the trait (the only addition it needs is a mutable accessor to
//! its same-domain map).
//!
//! Algorithm overview (matching the C++ flow):
//!
//! 1. [`build_split_faces`] — for every source face carrying on-face
//!    (IN/section) edges in the BOPDS face-info pool, collect the face's
//!    bounding edges (substituting split images) plus the on-face edges, close
//!    the set into wires with [`crate::wire_splitter::WireSplitter`], then
//!    rebuild one closed face per wire on the original surface with
//!    [`crate::builder_face::FaceBuilder`]. The resulting faces are recorded as
//!    images of the source face.
//! 2. [`fill_same_domain_faces`] — group the result faces by (boundary edge
//!    signature, geometric surface); each group of coincident faces collapses
//!    to one representative recorded in the same-domain map
//!    ([`BopBuilderLike::bind_shapes_sd`]).
//!
//! Dependencies: [`crate::builder_face`], [`crate::wire_splitter`],
//! [`crate::bopds`], [`crate::pave_filler`], [`crate::algo_tools`],
//! [`crate::bop_builder2`].

use std::collections::HashMap;

use occt_core::gp::{GpPnt2d, GpVec};
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bop_hist::BopHistory;
use crate::bopds::BopdsDS;
use crate::boptools_2d;
use crate::brep_surface::{classify_surface, surface_closest_params, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::builder_area::AreaBuilder;
use crate::builder_face::FaceBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, vertex_position};
use crate::wire_splitter::{WireEdgeSet, WireSplitter};

/// The context a face-reconstruction stage exposes to the BOPDS + history.
///
/// Implemented by `crate::bop_builder2::BopBuilder` once it gains a mutable
/// same-domain accessor; the local test stub provides the same surface so this
/// module is verifiable in isolation. Mirrors the members of `BOPAlgo_Builder`
/// that `FillImagesFaces`/`BuildSplitFaces`/`FillSameDomainFaces` touch.
pub trait BopBuilderLike {
    /// Read access to the BOPDS.
    fn ds(&self) -> &BopdsDS;
    /// Mutable access to the BOPDS.
    fn ds_mut(&mut self) -> &mut BopdsDS;
    /// The images naming table (`BOPAlgo_Builder::myImages`).
    fn history(&self) -> &BopHistory;
    /// Mutable access to the images naming table.
    fn history_mut(&mut self) -> &mut BopHistory;
    /// True when the builder has accumulated a fatal error.
    fn has_errors(&self) -> bool;
    /// Record a fatal error.
    fn add_error(&mut self, msg: String);
    /// Record a non-fatal warning.
    fn add_warning(&mut self, msg: String);
    /// Whether the pave-filler runs in non-destructive mode
    /// (`myPaveFiller->NonDestructive()`).
    fn non_destructive(&self) -> bool;
    /// The fuzzy tolerance of the operation (`myFuzzyValue`).
    fn fuzzy_value(&self) -> f64;
    /// Bind `shape` as same-domain with `sd` (`myShapesSD.Bind`).
    fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape);
    /// The same-domain representative of `shape`, when bound (`myShapesSD.Seek`).
    fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape>;
}

// ---------------------------------------------------------------------------
// FillImagesFaces
// ---------------------------------------------------------------------------

/// Run the whole face-image construction: [`build_split_faces`] then
/// [`fill_same_domain_faces`], aborting when either stage reported errors.
///
/// Mirrors `BOPAlgo_Builder::FillImagesFaces` (the `FillInternalVertices`
/// epilogue lives in `crate::bop_build_common::fill_internal_vertices`).
pub fn fill_images_faces<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    build_split_faces(f)?;
    if f.has_errors() {
        return Err("FillImagesFaces: BuildSplitFaces reported errors".to_string());
    }
    fill_same_domain_faces(f)?;
    if f.has_errors() {
        return Err("FillImagesFaces: FillSameDomainFaces reported errors".to_string());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// BuildSplitFaces
// ---------------------------------------------------------------------------

/// Rebuild the split pieces of every source face that was cut by the
/// intersection and record them as images of that face.
///
/// For each source face with a non-empty face-info record (i.e. with on-face
/// IN/section edges), the boundary edge set is collected — each bounding edge
/// replaced by its split images when present — and the on-face edges are
/// appended both FORWARD and REVERSED so the wire splitter can close loops on
/// both sides. The set is closed into wires by
/// [`crate::wire_splitter::WireSplitter`]; each closed wire is rebuilt into a
/// face on the *original* surface by [`FaceBuilder`]. Every resulting face is
/// appended to the images of the source face (re-applying the source face
/// orientation). Mirrors `BOPAlgo_Builder::BuildSplitFaces`.
///
/// Per-face failures are recorded with [`BopBuilderLike::add_error`] and the
/// remaining faces are still processed; the caller aborts on
/// [`BopBuilderLike::has_errors`].
///
/// *Simplifications vs OCCT:* the draft-face fast path (`BuildDraftFace`), the
/// seam-closing step for periodic surfaces (`DoSplitSEAMOnFace`), the
/// `IsSplitToReverseWithWarn` re-orientation and the alone-vertex handling are
/// not ported; the planar-face p-curve speed-up is replaced by a best-effort
/// p-curve attachment on the split faces.
pub fn build_split_faces<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();

    /// One face to split: the source face and the edge set that bounds its
    /// split pieces. `on_edges` are the on-face (section) edges separately,
    /// for the planar arrangement fast path.
    struct SplitTask {
        face_index: usize,
        face: Face,
        edges: Vec<Edge>,
        on_edges: Vec<Edge>,
    }
    let mut tasks: Vec<SplitTask> = Vec::new();

    // 1. Prepare the per-face tasks.
    for i in 0..n {
        let is_face = f
            .ds()
            .shape_info(i)
            .map(|s| s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        if !is_face {
            continue;
        }
        let Some(face_shape) = f.ds().shape(i).cloned() else { continue };
        let face = Face(face_shape);

        // The face participates only when the intersection recorded on-face
        // (IN/ON/section) edges for it. The face-info `paves` list mixes edge
        // indices (E/F, F/F section edges) with vertex indices (V/F vertices on
        // the face); only the edges are on-face section edges. Each pave carries
        // the `(first, last)` parameter range of the edge lying on the face, so
        // the split face is bounded by the *trimmed* sub-edge.
        let on_edges: Vec<(usize, f64, f64)> = {
            let pool = f.ds().face_info_pool();
            let mut v: Vec<(usize, f64, f64)> = Vec::new();
            if let Some(fi) = pool.iter().find(|fi| fi.face_index == i) {
                for &(e, fl, ll) in &fi.paves {
                    let is_edge = f
                        .ds()
                        .shape(e)
                        .map(|s| s.shape_type() == ShapeType::Edge)
                        .unwrap_or(false);
                    if is_edge && !v.iter().any(|&(x, _, _)| x == e) {
                        v.push((e, fl, ll));
                    }
                }
            }
            v
        };
        if on_edges.is_empty() {
            continue;
        }

        // Normalize the face to FORWARD for edge collection; the original
        // orientation is re-applied to the split faces at the end.
        let mut ff = face.clone();
        ff.0.set_orientation(Orientation::Forward);

        // 1.1 Bounding edges of the face, substituting split images.
        let mut le: Vec<Edge> = Vec::new();
        for be in edges_of(&ff.0) {
            let an_ori = be.0.orientation();
            if !f.history().has_image(&be.0) {
                if an_ori == Orientation::Internal {
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

            // The edge was split during intersection: use its image pieces.
            let splits = f.history().image(&be.0).map(|s| s.to_vec()).unwrap_or_default();
            for mut sp in splits {
                if BRepTool::is_degenerated(&be) {
                    sp.set_orientation(an_ori);
                    le.push(Edge(sp));
                    continue;
                }
                if an_ori == Orientation::Internal {
                    sp.set_orientation(Orientation::Forward);
                    le.push(Edge(sp.clone()));
                    sp.set_orientation(Orientation::Reversed);
                    le.push(Edge(sp));
                } else {
                    sp.set_orientation(an_ori);
                    le.push(Edge(sp));
                }
            }
        }

        // 1.2 + 1.3 On-face (IN/section) edges — both orientations. A pave
        // block of the on-face edge covers exactly the part of that edge lying
        // on the face; its split edge (when the block was split) is the
        // properly-trimmed sub-edge. Fall back to the whole edge otherwise.
        let mut on_le: Vec<Edge> = Vec::new();
        for &(e_idx, fl, ll) in &on_edges {
            let on_edge = on_face_split_edge(f, e_idx, fl, ll);
            let Some(mut sp) = on_edge else { continue };
            sp.set_orientation(Orientation::Forward);
            le.push(Edge(sp.clone()));
            on_le.push(Edge(sp.clone()));
            sp.set_orientation(Orientation::Reversed);
            le.push(Edge(sp));
        }

        tasks.push(SplitTask { face_index: i, face, edges: le, on_edges: on_le });
    }

    // 2. Execute the tasks: close each edge set into wires, then build one
    //    closed face per wire on the original surface.
    let mut faces_im: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for task in tasks {
        // The on-face edges are added FORWARD and REVERSED, and the boundary
        // edges of coincident faces coincide geometrically (A's edge and B's
        // edge over the shared segment). The duplicates are needed: an interior
        // section edge is shared by two split pieces, so the wire splitter must
        // see it in each loop.
        let mut wes = WireEdgeSet::new();
        wes.set_face(task.face.clone());
        wes.add_edges(&task.edges);
        let mut ws = WireSplitter::new();
        ws.set_wes(wes);
        if let Err(e) = ws.perform() {
            f.add_error(format!(
                "BuildSplitFaces: face {} wire splitting failed: {e}",
                task.face_index
            ));
            continue;
        }
        let raw_wires = ws.wires().to_vec();
        let wires = filter_wires(&raw_wires, &task.face);
        if wires.is_empty() {
            continue;
        }

        let mut ff = task.face.clone();
        ff.0.set_orientation(Orientation::Forward);
        // Group the closed wires into (growth outer loop, hole loops) per
        // `BOPAlgo_BuilderFace::PerformAreas`: a wire strictly containing
        // another wire's region is a growth, the inner wire is its hole. Each
        // growth becomes a face bounded by its outer loop plus the holes that
        // fall inside it. Every hole region is also a split piece on its own
        // (the disk cut out of the outer loop), rebuilt as an independent face.
        let wire_edges_list: Vec<Vec<Edge>> = wires.iter().map(|w| edges_of_wire(&Wire(w.clone()))).collect();
        let mut used = vec![false; wire_edges_list.len()];
        for i in 0..wire_edges_list.len() {
            if used[i] {
                continue;
            }
            used[i] = true;
            let mut holes: Vec<Vec<Edge>> = Vec::new();
            for j in (i + 1)..wire_edges_list.len() {
                if used[j] {
                    continue;
                }
                if loop_contains(&wire_edges_list[i], &wire_edges_list[j]) {
                    used[j] = true;
                    holes.push(wire_edges_list[j].clone());
                }
            }
            let mut all_edges = wire_edges_list[i].clone();
            for h in &holes {
                all_edges.extend(h.iter().cloned());
            }
            let shapes: Vec<TopoShape> = all_edges.iter().map(|e| e.0.clone()).collect();
            match crate::builder_face::build_face_with_holes(&wire_edges_list[i], &holes) {
                Ok(face) => {
                    attach_pcurves(&face, &all_edges);
                    faces_im.entry(task.face_index).or_default().push(face.0);
                }
                Err(_) => {
                    // Fall back to a plain single-wire face (no holes attach).
                    let mut fb = FaceBuilder::new();
                    fb.set_face(&ff);
                    fb.set_shapes(&shapes);
                    if let Err(e2) = fb.perform() {
                        f.add_error(format!(
                            "BuildSplitFaces: face {} area construction failed: {e2}",
                            task.face_index
                        ));
                        continue;
                    }
                    for area in fb.areas() {
                        faces_im.entry(task.face_index).or_default().push(area.clone());
                    }
                }
            }
        }
    }

    // 3. Record the split faces as images, preserving the source orientation.
    for (face_idx, split_faces) in faces_im {
        let Some(orig) = f.ds().shape(face_idx).cloned() else { continue };
        let an_ori = orig.orientation();
        for sf in split_faces {
            let mut sf = sf;
            if an_ori == Orientation::Reversed {
                sf.set_orientation(Orientation::Reversed);
            }
            f.history_mut().add_image(&orig, sf);
        }
    }
    Ok(())
}

/// The sub-edge of the on-face edge `e_idx` covering the pave range
/// `[fl, ll]`: the pave block of `e_idx` whose range matches `[fl, ll]` has a
/// split edge (`pb.edge()`) that is exactly the trimmed on-face part. When no
/// block matches (or it was not split), the whole edge is used.
fn on_face_split_edge<B: BopBuilderLike>(
    f: &B,
    e_idx: usize,
    fl: f64,
    ll: f64,
) -> Option<TopoShape> {
    let tol = 1e-7;
    let blocks = f.ds().pave_blocks(e_idx);
    if let Some(pb) = blocks.iter().find(|pb| {
        let (a, b) = pb.range();
        ((a - fl).abs() <= tol && (b - ll).abs() <= tol)
            || ((a - ll).abs() <= tol && (b - fl).abs() <= tol)
    }) {
        if let Some(sp) = f.ds().shape(pb.edge()).cloned() {
            return Some(sp);
        }
    }
    f.ds().shape(e_idx).cloned()
}

/// Whether the region bounded by `outer` strictly contains the region bounded
/// by `inner`: `inner` has the smaller |area| and one of its vertices lies
/// strictly inside `outer`. Mirrors the hole→growth attachment of
/// `BOPAlgo_BuilderFace::PerformAreas` (`IsInside` on the 2-D classification).
fn loop_contains(outer: &[Edge], inner: &[Edge]) -> bool {
    let op = match crate::builder_area::plane_from_loop(outer) {
        Ok(p) => p,
        Err(_) => return false,
    };
    let a_outer = crate::builder_face::loop_signed_area(outer, &op).abs();
    let a_inner = crate::builder_face::loop_signed_area(inner, &op).abs();
    if a_inner >= a_outer {
        return false;
    }
    // Project `outer` into its own plane frame and test one vertex of `inner`
    // with the ray-crossing test (strictly inside: on-boundary → false).
    let outer_poly = project_loop_2d(outer, &op);
    let inner_poly = project_loop_2d(inner, &op);
    let Some(&p) = inner_poly.first() else { return false };
    point_in_polygon(&outer_poly, &p, false)
}

/// Project an edge loop onto the plane frame `(u, v)` of `pln`.
fn project_loop_2d(edges: &[Edge], pln: &occt_core::gp::GpPln) -> Vec<occt_core::gp::GpPnt2d> {
    let xd = *pln.position().x_direction().xyz();
    let yd = *pln.position().y_direction().xyz();
    let loc = pln.position().location();
    let mut out = Vec::new();
    for e in edges {
        let (a, _) = edge_vertices(e);
        let Some(a) = a else { continue };
        let p = vertex_position(&a);
        let v = p.coord.subtracted(&loc.coord);
        out.push(occt_core::gp::GpPnt2d::new(v.dot(&xd), v.dot(&yd)));
    }
    out
}

/// Ray-crossing point-in-polygon test on the plane; `on_edge_is_inside=false`
/// makes a boundary point count as outside.
fn point_in_polygon(poly: &[occt_core::gp::GpPnt2d], p: &occt_core::gp::GpPnt2d, on_edge_is_inside: bool) -> bool {
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        let ((ax, ay), (bx, by)) = ((a.x(), a.y()), (b.x(), b.y()));
        // Boundary test (point on segment).
        let cross = (p.x() - ax) * (by - ay) - (p.y() - ay) * (bx - ax);
        let seg_len2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
        if seg_len2 > 1e-24 {
            let t = (((p.x() - ax) * (bx - ax) + (p.y() - ay) * (by - ay)) / seg_len2).clamp(0.0, 1.0);
            let (qx, qy) = (ax + (bx - ax) * t, ay + (by - ay) * t);
            if (p.x() - qx).hypot(p.y() - qy) < 1e-9 {
                return on_edge_is_inside;
            }
        }
        // Ray-crossing.
        if ((ay > p.y()) != (by > p.y()))
            && (p.x() < (bx - ax) * (p.y() - ay) / (by - ay + 1e-30) + ax)
        {
            inside = !inside;
        }
    }
    inside
}

/// Keeps the closed wires that bound genuine split pieces of `face`: a
/// degenerate wire (zero area) cannot bound a face, and duplicate wires that
/// bound the same region collapse to one representative.
fn filter_wires(wires: &[TopoShape], face: &Face) -> Vec<TopoShape> {
    let pln = crate::brep_surface::face_plane(face);
    let mut seen: Vec<Vec<((i64, i64, i64), (i64, i64, i64))>> = Vec::new();
    let mut out: Vec<TopoShape> = Vec::new();
    for w in wires {
        let we = edges_of_wire(&Wire(w.clone()));
        let area = match &pln {
            Some(p) => wire_signed_area(&we, p).abs(),
            None => 0.0,
        };
        if area <= 1e-7 {
            continue;
        }
        let sig = wire_signature(&we);
        if sig.is_empty() || seen.contains(&sig) {
            continue;
        }
        seen.push(sig);
        out.push(w.clone());
    }
    out
}

/// Signed area of a closed wire on the plane, from the unique vertex positions
/// of its edges (shoelace over the polygon vertices). Unlike
/// `builder_face::loop_signed_area`, this does not require the wire edges to be
/// stored in chain order — it collects the distinct endpoint points and walks
/// them as the polygon. Zero for a degenerate (self-coincident) wire.
fn wire_signed_area(edges: &[Edge], pln: &occt_core::gp::GpPln) -> f64 {
    use std::collections::BTreeSet;
    let mut pts: Vec<occt_core::gp::GpPnt> = Vec::new();
    for e in edges {
        let (a, b) = edge_vertices(e);
        for v in [a, b].into_iter().flatten() {
            pts.push(vertex_position(&v));
        }
    }
    if pts.len() < 3 {
        return 0.0;
    }
    // Project onto the plane frame; collect distinct points (a closed ring has
    // n edges → n distinct vertices; the closing vertex equals the first).
    let xd = *pln.position().x_direction().xyz();
    let yd = *pln.position().y_direction().xyz();
    let loc = pln.position().location();
    let mut ring: Vec<GpPnt2d> = Vec::new();
    let mut seen: BTreeSet<(i64, i64)> = BTreeSet::new();
    for p in pts {
        let v = p.coord.subtracted(&loc.coord);
        let q = GpPnt2d::new(v.dot(&xd), v.dot(&yd));
        let k = ((q.x() / 1e-6).round() as i64, (q.y() / 1e-6).round() as i64);
        if seen.insert(k) {
            ring.push(q);
        }
    }
    if ring.len() < 3 {
        return 0.0;
    }
    let mut acc = 0.0;
    let n = ring.len();
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        acc += a.x() * b.y() - b.x() * a.y();
    }
    0.5 * acc
}

/// Sorted multiset of an edge loop's quantized endpoint pairs — the identity
/// of the bounded region.
fn wire_signature(edges: &[Edge]) -> Vec<((i64, i64, i64), (i64, i64, i64))> {
    const TOL: f64 = 1e-6;
    let key = |p: &occt_core::gp::GpPnt| -> (i64, i64, i64) {
        (
            (p.x() / TOL).round() as i64,
            (p.y() / TOL).round() as i64,
            (p.z() / TOL).round() as i64,
        )
    };
    let mut sig: Vec<((i64, i64, i64), (i64, i64, i64))> = Vec::new();
    for e in edges {
        let (a, b) = edge_vertices(e);
        let (Some(a), Some(b)) = (a, b) else { return Vec::new() };
        let (ka, kb) = (key(&vertex_position(&a)), key(&vertex_position(&b)));
        sig.push(if ka <= kb { (ka, kb) } else { (kb, ka) });
    }
    sig.sort_unstable();
    sig
}

/// Attach the p-curve of every edge of a freshly built split face onto that
/// face, when a p-curve is not already registered. Mirrors the OCCT
/// `BRepLib::BuildPCurveForEdgesOnPlane` step; a failure only warns.
fn attach_pcurves(face: &Face, edges: &[Edge]) {
    for e in edges {
        if boptools_2d::curve_on_surface(e, face).is_some() {
            continue;
        }
        if let Ok(pc) = AlgoTools::make_pcurve(e, face) {
            let face_key = GeometryRegistry::shape_key(&face.0);
            GeometryRegistry::global().set_edge_pcurve(&e.0, face_key, pc);
        }
    }
}

// ---------------------------------------------------------------------------
// FillSameDomainFaces
// ---------------------------------------------------------------------------

/// Collapse coincident faces of the result into a single representative.
///
/// Collects every result face (each source face or, when split, its images),
/// groups them by (boundary edge signature, geometric surface) and, for each
/// group of two or more coincident faces, binds every non-representative face
/// to the representative in the same-domain map. The representative is an
/// original source face with the smallest BOPDS index when one exists, else
/// the first face of the group. Original faces that acquire a same-domain
/// twin are appended to their own image list, mirroring the OCCT step that
/// marks an unchanged face as "split" once it has a same-domain companion.
///
/// Mirrors `BOPAlgo_Builder::FillSameDomainFaces` with the geometric
/// same-domain test reduced to [`brep_surface::classify_surface`] + surface
/// sampling (planar faces are compared exactly by plane coincidence).
///
/// *Note:* the caller resolves images through [`BopBuilderLike::seek_shapes_sd`]
/// when assembling the final result, so each SD group contributes only its
/// representative.
pub fn fill_same_domain_faces<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();

    // Collect the faces that participate in the result.
    let mut all_faces: Vec<TopoShape> = Vec::new();
    for i in 0..n {
        let is_face = f
            .ds()
            .shape_info(i)
            .map(|s| s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        if !is_face {
            continue;
        }
        let Some(orig) = f.ds().shape(i).cloned() else { continue };
        match f.history().image(&orig) {
            Some(imgs) => {
                all_faces.extend(imgs.iter().cloned())
            }
            None => {
                all_faces.push(orig)
            }
        }
    }
    if all_faces.len() < 2 {
        return Ok(());
    }

    // Group by (edge signature, geometric surface).
    let tol = f.fuzzy_value().max(1e-7);
    let mut groups: Vec<Vec<TopoShape>> = Vec::new();
    let mut used = vec![false; all_faces.len()];
    for i in 0..all_faces.len() {
        if used[i] {
            continue;
        }
        let mut group = vec![all_faces[i].clone()];
        used[i] = true;
        for j in (i + 1)..all_faces.len() {
            if used[j] {
                continue;
            }
            let sd_ij = faces_same_domain(&Face(all_faces[i].clone()), &Face(all_faces[j].clone()), tol);
            if i == 0 && j == 6 {
            }
            if sd_ij {
                group.push(all_faces[j].clone());
                used[j] = true;
            }
        }
        groups.push(group);
    }

    // Fill the same-domain map.
    for group in groups {
        if group.len() < 2 {
            continue;
        }
        // Representative: an original face (present in the DS) with the
        // smallest index, else the first face of the group.
        let mut rep: Option<TopoShape> = None;
        let mut rep_idx = usize::MAX;
        for gf in &group {
            if let Some(idx) = f.ds().index(gf) {
                if idx < rep_idx {
                    rep_idx = idx;
                    rep = Some(gf.clone());
                }
            }
        }
        let rep = rep.unwrap_or_else(|| group[0].clone());

        // Bind every non-representative face to the representative.
        for gf in &group {
            if !gf.same_tshape(&rep) {
                let gi = f.ds().index(gf).map(|x| x.to_string()).unwrap_or("?".into());
                let ri = f.ds().index(&rep).map(|x| x.to_string()).unwrap_or("?".into());
                let is_src = f.ds().nb_source_shapes();
                f.bind_shapes_sd(gf.clone(), rep.clone());
            }
        }
        // Original faces in the group get themselves as images (OCCT marks an
        // unchanged face that has a same-domain twin as split into itself).
        // A face that already carries split images is not additionally marked
        // as its own image: its splits replace it in the result, and adding a
        // self-image on top would make the draft-solid assembly count the face
        // both as itself and as its pieces.
        for gf in &group {
            if f.ds().index(gf).is_some() && !f.history().has_image(gf) {
                f.history_mut().add_image(gf, gf.clone());
            }
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Same-domain geometric helpers
// ---------------------------------------------------------------------------

/// Grid size used to quantize vertex coordinates for the edge signature.
const VERTEX_TOL: f64 = 1e-6;

/// Quantized identity of a vertex position (matches `wire_splitter`'s grid).
fn vkey(v: &Vertex) -> (i64, i64, i64) {
    let p = vertex_position(v);
    (
        (p.x() / VERTEX_TOL).round() as i64,
        (p.y() / VERTEX_TOL).round() as i64,
        (p.z() / VERTEX_TOL).round() as i64,
    )
}

/// Sorted multiset of the face's boundary edges, each identified by its
/// endpoint vertex positions. Two coincident faces of the same shape (from
/// different operands, hence different `TShape`s) produce equal signatures.
fn edge_signature(face: &Face) -> Vec<((i64, i64, i64), (i64, i64, i64))> {
    let mut sig: Vec<((i64, i64, i64), (i64, i64, i64))> = edges_of(&face.0)
        .iter()
        .filter_map(|e| {
            let (a, b) = edge_vertices(e);
            match (a, b) {
                (Some(a), Some(b)) => Some((vkey(&a), vkey(&b))),
                _ => None,
            }
        })
        .collect();
    sig.sort_unstable();
    sig
}

/// Whether two faces are Same-Domain: identical boundary and coincident
/// supporting surfaces.
fn faces_same_domain(f1: &Face, f2: &Face, tol: f64) -> bool {
    if edge_signature(f1) != edge_signature(f2) {
        return false;
    }
    let (Some(s1), Some(s2)) = (
        BRepTool::face_surface_world(f1),
        BRepTool::face_surface_world(f2),
    ) else {
        return false;
    };
    surfaces_match(s1.as_ref(), s2.as_ref(), tol)
}

/// Whether two surfaces describe the same geometric surface within `tol`.
///
/// Planes are compared exactly (parallel normals and coincident location);
/// other surface types are compared by sampling the first surface on a grid
/// and projecting each sample onto the second.
fn surfaces_match(s1: &dyn Surface, s2: &dyn Surface, tol: f64) -> bool {
    let k1 = classify_surface(s1);
    if k1 != classify_surface(s2) {
        return false;
    }
    match k1 {
        SurfaceKind::Plane => {
            let n1 = crate::brep_surface::surface_normal(s1, 0.0, 0.0);
            let n2 = crate::brep_surface::surface_normal(s2, 0.0, 0.0);
            // Parallel planes: the normals are collinear (either direction).
            if n1.xyz().crossed(n2.xyz()).modulus() > tol {
                return false;
            }
            // Coincident planes: the location of one lies on the other.
            let p1 = s1.d0(0.0, 0.0);
            let p2 = s2.d0(0.0, 0.0);
            let offset = GpVec::from_pnts(&p1, &p2).dot(&n1).abs();
            offset <= tol
        }
        _ => {
            // Sampled comparison: every sample point of s1 lies on s2. The
            // tolerance is floored at 1e-3 to absorb the closest-parameter
            // grid-search error.
            let stol = tol.max(1e-3);
            let (u0, u1, v0, v1) = sample_bounds(s1);
            let (nu, nv) = (8, 8);
            for i in 0..nu {
                for j in 0..nv {
                    let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
                    let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
                    let p = s1.d0(u, v);
                    let (pu, pv) = surface_closest_params(s2, &p, 16, 16);
                    if p.distance(&s2.d0(pu, pv)) > stol {
                        return false;
                    }
                }
            }
            true
        }
    }
}

/// Finite sampling bounds of a surface (unbounded ranges clamp to `[-1, 1]`).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bopds::BopdsFaceInfo;
    use crate::brep_surface::face_plane;
    use crate::builder::TopoBuilder;
    use crate::builder_face::loop_signed_area;
    use crate::brep_extrema::test_box::unit_box;
    use crate::topo_tools_full::wires_of_face;
    use occt_core::gp::GpPnt;

    /// A minimal `BopBuilderLike` host for the isolated tests.
    #[derive(Debug)]
    struct StubBopBuilder {
        ds: BopdsDS,
        history: BopHistory,
        errors: Vec<String>,
        warnings: Vec<String>,
        shapes_sd: Vec<(TopoShape, TopoShape)>,
        fuzzy: f64,
    }

    impl StubBopBuilder {
        fn new() -> Self {
            Self {
                ds: BopdsDS::new(),
                history: BopHistory::new(),
                errors: Vec::new(),
                warnings: Vec::new(),
                shapes_sd: Vec::new(),
                fuzzy: 1e-7,
            }
        }
    }

    impl BopBuilderLike for StubBopBuilder {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn ds_mut(&mut self) -> &mut BopdsDS {
            &mut self.ds
        }
        fn history(&self) -> &BopHistory {
            &self.history
        }
        fn history_mut(&mut self) -> &mut BopHistory {
            &mut self.history
        }
        fn has_errors(&self) -> bool {
            !self.errors.is_empty()
        }
        fn add_error(&mut self, msg: String) {
            self.errors.push(msg);
        }
        fn add_warning(&mut self, msg: String) {
            self.warnings.push(msg);
        }
        fn non_destructive(&self) -> bool {
            false
        }
        fn fuzzy_value(&self) -> f64 {
            self.fuzzy
        }
        fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
            if let Some(e) = self.shapes_sd.iter_mut().find(|(s, _)| s.same_tshape(&shape)) {
                e.1 = sd;
            } else {
                self.shapes_sd.push((shape, sd));
            }
        }
        fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
            self.shapes_sd
                .iter()
                .find(|(s, _)| s.same_tshape(shape))
                .map(|(_, sd)| sd.clone())
        }
    }

    /// True polygon area of a planar face: the sum of |signed loop areas| of
    /// its wires projected on the supporting plane. Unlike `face_area_exact`
    /// (which triangulates the whole UV rectangle), this gives the trimmed
    /// area of an unbounded-plane face.
    fn polygon_area_on_plane(face: &Face) -> f64 {
        let Some(pln) = face_plane(face) else { return f64::NAN };
        wires_of_face(face)
            .iter()
            .map(|w| loop_signed_area(&edges_of_wire(w), &pln).abs())
            .sum()
    }

    /// Record the diagonal `(0,0,0)-(1,1,0)` as an on-face edge of the bottom
    /// face (`z = 0`) of `ub`, in the stub's BOPDS.
    fn stub_with_split_bottom_face(ub: &crate::brep_extrema::test_box::UnitBox) -> (StubBopBuilder, usize) {
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        let bottom = &ub.faces[0];
        let face_idx = stub.ds().index(&bottom.0).expect("bottom face indexed");

        let b = TopoBuilder::new();
        let diag = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let diag_idx = stub.ds_mut().append(diag.0.clone()).expect("append diagonal edge");

        let pool = stub.ds_mut().change_face_info_pool();
        let fi = match pool.iter_mut().find(|fi| fi.face_index == face_idx) {
            Some(fi) => fi,
            None => {
                pool.push(BopdsFaceInfo::new(face_idx));
                pool.last_mut().expect("just pushed")
            }
        };
        fi.add_pave(diag_idx, 0.0, 1.0);
        (stub, face_idx)
    }

    #[test]
    fn box_face_split_by_diagonal_builds_two_closed_faces() {
        let ub = unit_box();
        let (mut stub, face_idx) = stub_with_split_bottom_face(&ub);

        build_split_faces(&mut stub).expect("build split faces");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);
        assert!(stub.warnings.is_empty(), "warnings: {:?}", stub.warnings);

        // The bottom face must have exactly two images (the two triangles).
        let orig = stub.ds().shape(face_idx).cloned().expect("face shape");
        let imgs = stub.history().image(&orig).expect("images recorded");
        assert_eq!(imgs.len(), 2, "two split faces, got {}", imgs.len());

        let a1 = polygon_area_on_plane(&Face(imgs[0].clone()));
        let a2 = polygon_area_on_plane(&Face(imgs[1].clone()));
        assert!(
            (a1 + a2 - 1.0).abs() < 1e-6,
            "split areas must sum to the original face area: {a1} + {a2}"
        );
        assert!(
            (a1 - 0.5).abs() < 1e-6 && (a2 - 0.5).abs() < 1e-6,
            "each diagonal half is a triangle of area 0.5, got {a1}, {a2}"
        );

        // Both split faces are closed (single wire, closed loop).
        for im in imgs {
            let f = Face(im.clone());
            assert_eq!(wires_of_face(&f).len(), 1, "one wire per split face");
        }
    }

    #[test]
    fn fill_same_domain_merges_coincident_keeps_parallel() {
        let ub1 = unit_box();
        let ub2 = unit_box();
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub1.solid.0.clone(), ub2.solid.0.clone()]);

        fill_same_domain_faces(&mut stub).expect("fill same domain");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);

        // Two coincident bottom faces merge: the second box's bottom is SD of
        // the first box's bottom (the smaller DS index is the representative).
        let f1 = &ub1.faces[0]; // bottom z = 0
        let f2 = &ub2.faces[0]; // bottom z = 0, same position
        let sd = stub.seek_shapes_sd(&f2.0);
        assert!(sd.is_some(), "coincident bottom faces must merge");
        assert!(sd.unwrap().same_tshape(&f1.0), "representative is the first box's bottom");

        // Parallel but distinct faces of one box do not merge: the bottom and
        // the top of the same box are never bound together.
        let f_top = &ub1.faces[1]; // top z = 1
        let sd_bottom = stub.seek_shapes_sd(&f1.0);
        let sd_top = stub.seek_shapes_sd(&f_top.0);
        assert!(!sd_bottom.map_or(false, |s| s.same_tshape(&f_top.0)));
        assert!(!sd_top.map_or(false, |s| s.same_tshape(&f1.0)));

        // Original faces that acquired a same-domain twin got themselves as
        // images (they are unchanged but now have a representative).
        assert!(stub.history().has_image(&f1.0));
        assert!(stub.history().has_image(&f2.0));
    }

    #[test]
    fn fill_same_domain_single_box_produces_no_merge() {
        let ub = unit_box();
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        fill_same_domain_faces(&mut stub).expect("fill same domain");
        assert!(!stub.has_errors());

        // None of the six distinct box faces is coincident with another.
        assert!(stub.shapes_sd.is_empty(), "no same-domain bindings expected");
        assert!(!stub.history().has_any_images(), "no face should become its own image");
    }

    #[test]
    fn fill_images_faces_runs_clean_without_face_info() {
        let ub = unit_box();
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        // No face carries face-info: nothing to split, nothing to merge.
        fill_images_faces(&mut stub).expect("fill images faces");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);
        assert!(!stub.history().has_any_images(), "no split faces recorded");
        assert!(stub.shapes_sd.is_empty(), "no same-domain bindings");
    }

    #[test]
    fn face_without_face_info_is_skipped() {
        let ub = unit_box();
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        // Only the top face carries a face-info entry — and it is empty.
        let top = &ub.faces[1];
        let top_idx = stub.ds().index(&top.0).expect("top face indexed");
        stub.ds_mut().change_face_info_pool().push(BopdsFaceInfo::new(top_idx));

        build_split_faces(&mut stub).expect("build split faces");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);
        // An empty face-info record means the face was not split.
        assert!(!stub.history().has_any_images(), "no split faces recorded");
    }
}
