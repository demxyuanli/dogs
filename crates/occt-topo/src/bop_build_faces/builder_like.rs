use super::prelude::*;
use super::*;

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
    /// Origins back-map (`myOrigins`) keyed by TShape identity.
    fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>>;
}

// ---------------------------------------------------------------------------
// FillImagesFaces
// ---------------------------------------------------------------------------

/// Run the whole face-image construction: [`build_split_faces`],
/// [`fill_same_domain_faces`], then
/// [`crate::bop_fill_internal_verts::fill_internal_vertices_occt`].
///
/// Mirrors `BOPAlgo_Builder::FillImagesFaces` (`_2.cxx:215-229`).
pub fn fill_images_faces<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    build_split_faces(f)?;
    if f.has_errors() {
        return Err("FillImagesFaces: BuildSplitFaces reported errors".to_string());
    }
    fill_same_domain_faces(f)?;
    if f.has_errors() {
        return Err("FillImagesFaces: FillSameDomainFaces reported errors".to_string());
    }
    crate::bop_fill_internal_verts::fill_internal_vertices_occt(f)?;
    if f.has_errors() {
        return Err("FillImagesFaces: FillInternalVertices reported errors".to_string());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// BuildSplitFaces
// ---------------------------------------------------------------------------

/// Rebuild the split pieces of every source face that was cut by the
/// intersection and record them as images of that face.
///
/// Delegates to [`crate::bop_split_faces_occt::build_split_faces_occt`], the
/// full `BOPAlgo_Builder::BuildSplitFaces` translation (draft-face fast path,
/// `DoSplitSEAMOnFace`, `IsSplitToReverseWithWarn`, IN/Sc both ways).
pub fn build_split_faces<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    crate::bop_split_faces_occt::build_split_faces_occt(f)
}

pub(crate) fn collect_unique_edge_paves<B: BopBuilderLike>(
    f: &B,
    src: &[(usize, f64, f64)],
    dst: &mut Vec<(usize, f64, f64)>,
) {
    for &(e, fl, ll) in src {
        let is_edge = f
            .ds()
            .shape(e)
            .map(|s| s.shape_type() == ShapeType::Edge)
            .unwrap_or(false);
        if is_edge && !dst.iter().any(|&(x, _, _)| x == e) {
            dst.push((e, fl, ll));
        }
    }
}

pub(crate) fn append_on_face_edges<B: BopBuilderLike>(
    f: &B,
    _face: &Face,
    _ctx: &mut IntToolsContext,
    paves: &[(usize, f64, f64)],
    le: &mut Vec<Edge>,
    on_le: &mut Vec<Edge>,
) {
    for &(e_idx, fl, ll) in paves {
        for mut sp in on_face_split_edges(f, e_idx, fl, ll) {
            sp.set_orientation(Orientation::Forward);
            le.push(Edge(sp.clone()));
            on_le.push(Edge(sp.clone()));
            sp.set_orientation(Orientation::Reversed);
            le.push(Edge(sp));
        }
    }
}

/// `BOPAlgo_Builder_2.cxx:468-479`: `aPB->Edge()` is the split (or the
/// original when the block was not split). FaceInfo stores `(edge, first,
/// last)`; overlapping pave blocks of `e_idx` supply those splits. A range
/// miss used to fall back to the whole original, which can be a 3D edge that
/// only *touches* the face at a vertex and then poisons `ShapesToAvoid`.
pub(super) fn on_face_split_edges<B: BopBuilderLike>(
    f: &B,
    e_idx: usize,
    fl: f64,
    ll: f64,
) -> Vec<TopoShape> {
    let tol = 1e-7;
    let lo = fl.min(ll);
    let hi = fl.max(ll);
    let blocks = f.ds().pave_blocks(e_idx);
    let mut out = Vec::new();
    for pb in blocks {
        let (a, b) = pb.range();
        let (a, b) = (a.min(b), a.max(b));
        // Positive-length overlap with the FaceInfo range, not a shared end.
        // `_2.cxx` takes that one `aPB->Edge()`; a neighbour whose Edge() is
        // still the unsplit original poisons ShapesToAvoid.
        if a >= hi - tol || b <= lo + tol {
            continue;
        }
        let n_sp = pb.edge();
        if n_sp == 0 {
            continue;
        }
        if let Some(sp) = f.ds().shape(n_sp).cloned() {
            if sp.shape_type() == ShapeType::Edge
                && !out.iter().any(|s: &TopoShape| s.same_tshape(&sp))
            {
                out.push(sp);
            }
        }
    }
    if out.is_empty() && blocks.is_empty() {
        if let Some(s) = f.ds().shape(e_idx).cloned() {
            if s.shape_type() == ShapeType::Edge {
                out.push(s);
            }
        }
    }
    out
}

/// Whether the region bounded by `outer` strictly contains the region bounded
/// by `inner`: `inner` has the smaller |area| and one of its vertices lies
/// strictly inside `outer`. Mirrors the hole→growth attachment of
/// `BOPAlgo_BuilderFace::PerformAreas` (`IsInside` on the 2-D classification).
pub(super) fn loop_contains(outer: &[Edge], inner: &[Edge]) -> bool {
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

/// Group the closed split wires into (growth outer loop, hole loops) per
/// `BOPAlgo_BuilderFace::PerformAreas` (BOPAlgo_BuilderFace.cxx:387).
///
/// Mirrors the OCCT flow:
/// 1. classify each wire as a growth or a hole — the `IsGrowthWire` fast path
///    (a wire sharing an edge with a known hole wire is a growth) plus a
///    winding-independent geometric fallback for the remaining wires
///    (BOPAlgo_BuilderFace.cxx:441-458);
/// 2. attach every hole to the *nearest* growth that contains it — the
///    `IntTools_FClass2d` point-in-region test over the growth's draft face
///    (`IsInside`, BOPAlgo_BuilderFace.cxx:842-894), keeping the most-internal
///    containing growth (the `IsInside(aFace, *pFaceWas)` owner resolution,
///    BOPAlgo_BuilderFace.cxx:524-536).
///
/// `FClass2d::IsHole` is deliberately *not* used for the growth/hole decision:
/// the port's split wires are emitted in the surface frame regardless of the
/// material side, so a winding test misclassifies (the outer loop of a
/// box-bottom face reports `is_hole == true`). The classification stays
/// geometric; `FClass2d` is used only for the containment (`IsInside`) test,
/// which is winding-independent. The `Box2dTree` bbox cull of
/// BOPAlgo_BuilderFace.cxx:471-481 is skipped: the wire count per face is tiny,
/// so a plain scan over every growth is equivalent.
///
/// Returns, in wire order, the index of each growth wire and the indices of
/// the hole wires attached to it.
pub(crate) fn group_wires_as_areas(
    wire_edges: &[Vec<Edge>],
    surface: &Option<Arc<dyn Surface>>,
) -> Vec<(usize, Vec<usize>)> {
    let n = wire_edges.len();
    let edge_key = |e: &Edge| GeometryRegistry::shape_key(&e.0);

    // Phase 1 — growth/hole classification (PerformAreas lines 425-459).
    let mut mhe: HashSet<usize> = HashSet::new();
    let mut is_growth = vec![false; n];
    for i in 0..n {
        // IsGrowthWire fast path: the wire contains an edge of a known hole.
        if wire_edges[i].iter().any(|e| mhe.contains(&edge_key(e))) {
            is_growth[i] = true;
            continue;
        }
        // Geometric fallback: strictly contained in another wire → a hole.
        let mut hole = false;
        for (j, oes) in wire_edges.iter().enumerate() {
            if j != i && loop_contains(oes, &wire_edges[i]) {
                hole = true;
                break;
            }
        }
        if hole {
            for e in &wire_edges[i] {
                mhe.insert(edge_key(e));
            }
        } else {
            is_growth[i] = true;
        }
    }

    // Draft-face classifier of every growth (a single-wire face on the
    // original surface), used for the hole containment test.
    let mut cls: Vec<Option<FClass2d>> = vec![None; n];
    for i in 0..n {
        if is_growth[i] {
            if let Some(surf) = surface {
                if let Ok(face) = make_face_from_wire(&wire_edges[i], Some(surf.clone())) {
                    if let Ok(c) = FClass2d::new(&face, 1e-7) {
                        cls[i] = Some(c);
                    }
                }
            }
        }
    }

    // Phase 2 — attach each hole to the nearest (most internal) growth
    // (PerformAreas lines 486-537).
    let mut owners: Vec<Option<usize>> = vec![None; n];
    for h in 0..n {
        if is_growth[h] {
            continue;
        }
        for g in 0..n {
            if !is_growth[g] {
                continue;
            }
            if !hole_inside(&wire_edges[h], g, &wire_edges, &cls, surface) {
                continue;
            }
            match owners[h] {
                None => owners[h] = Some(g),
                Some(g0) => {
                    // Most-internal growth wins (OCCT `IsInside(aFace, *pFaceWas)`).
                    let new_inside_old = match &cls[g0] {
                        Some(co) => {
                            wire_point_state(&wire_edges[g], co, surface) == Some(FaceState::In)
                        }
                        None => loop_contains(&wire_edges[g0], &wire_edges[g]),
                    };
                    if new_inside_old {
                        owners[h] = Some(g);
                    }
                }
            }
        }
    }

    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for g in 0..n {
        if is_growth[g] {
            groups.push((g, (0..n).filter(|&h| owners[h] == Some(g)).collect()));
        }
    }
    groups
}

/// Whether the hole wire `hole` lies inside the growth wire `g` — the
/// `IntTools_FClass2d` `IsInside` test of PerformAreas (a hole edge midpoint
/// classified against the growth's draft-face region), falling back to the
/// geometric containment when no classifier is available. A hole sharing an
/// edge with the growth is never inside it (OCCT returns early at
/// BOPAlgo_BuilderFace.cxx:870).
pub(super) fn hole_inside(
    hole: &[Edge],
    g: usize,
    wire_edges: &[Vec<Edge>],
    cls: &[Option<FClass2d>],
    surface: &Option<Arc<dyn Surface>>,
) -> bool {
    let geometric = loop_contains(&wire_edges[g], hole);
    match (&cls[g], surface) {
        (Some(cl), Some(_)) => {
            let gk: HashSet<usize> = wire_edges[g]
                .iter()
                .map(|e| GeometryRegistry::shape_key(&e.0))
                .collect();
            if hole.iter().any(|e| gk.contains(&GeometryRegistry::shape_key(&e.0))) {
                return false;
            }
            geometric || wire_point_state(hole, cl, surface) == Some(FaceState::In)
        }
        _ => geometric,
    }
}

/// The `IntTools_FClass2d` state of a point on `wire` (one edge midpoint
/// projected to the surface UV), classified against the classifier `cl`.
/// `None` when no edge yields a usable point. Mirrors the OCCT `IsInside`
/// loop over the wire's edges (BOPAlgo_BuilderFace.cxx:859-892).
pub(super) fn wire_point_state(
    wire: &[Edge],
    cl: &FClass2d,
    surface: &Option<Arc<dyn Surface>>,
) -> Option<FaceState> {
    let surf = surface.as_ref()?;
    for e in wire {
        let (a, b) = edge_vertices(e);
        let (Some(va), Some(vb)) = (a, b) else { continue };
        let pa = vertex_position(&va);
        let pb = vertex_position(&vb);
        let mid = GpPnt::new(
            (pa.x() + pb.x()) / 2.0,
            (pa.y() + pb.y()) / 2.0,
            (pa.z() + pb.z()) / 2.0,
        );
        let (u, v) = surface_closest_params(surf.as_ref(), &mid, 16, 16);
        let st = cl.perform(GpPnt2d::new(u, v));
        if st != FaceState::Unknown {
            return Some(st);
        }
    }
    None
}

/// Project an edge loop onto the plane frame `(u, v)` of `pln`.
pub(super) fn project_loop_2d(edges: &[Edge], pln: &occt_core::gp::GpPln) -> Vec<occt_core::gp::GpPnt2d> {
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
pub(super) fn point_in_polygon(poly: &[occt_core::gp::GpPnt2d], p: &occt_core::gp::GpPnt2d, on_edge_is_inside: bool) -> bool {
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        let ((ax, ay), (bx, by)) = ((a.x(), a.y()), (b.x(), b.y()));
        // Boundary test (point on segment).
        let _cross = (p.x() - ax) * (by - ay) - (p.y() - ay) * (bx - ax);
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
/// Reverse the geometric direction of every edge of a loop, so the loop winds
/// opposite its original direction (`TopoDS::Reverse` on an edge list).
///
/// Each edge's two child vertices are swapped and its parameter range is
/// mirrored, so the child-order chaining used by
/// `builder_face::{build_loops, loop_signed_area}` walks the loop the other
/// way — the geometric winding flips. A B-Rep hole must wind opposite its
/// outer loop; the wire splitter emits every bounded region CCW, so the hole
/// attached to a growth face needs this reversal.
pub(crate) fn reverse_loop(edges: &[Edge]) -> Vec<Edge> {
    let reg = GeometryRegistry::global();
    edges
        .iter()
        .map(|e| {
            let mut nt = crate::tshape::TShape::new(ShapeType::Edge);
            for k in e.0.tshape.read().unwrap().children.iter().rev() {
                nt.add_child(k.clone());
            }
            let mut ne = TopoShape::from_handle(std::sync::Arc::new(std::sync::RwLock::new(nt)));
            ne.set_orientation(e.0.orientation().reversed());
            let ne = Edge(ne);
            if let Some(mut g) = reg.edge_geom(&e.0) {
                std::mem::swap(&mut g.first, &mut g.last);
                reg.set_edge(&ne.0, g);
            }
            ne
        })
        .collect()
}

pub(crate) fn filter_wires(wires: &[TopoShape], face: &Face) -> Vec<TopoShape> {
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
pub(super) fn wire_signed_area(edges: &[Edge], pln: &occt_core::gp::GpPln) -> f64 {
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
pub(super) fn wire_signature(edges: &[Edge]) -> Vec<((i64, i64, i64), (i64, i64, i64))> {
    pub(super) const TOL: f64 = 1e-6;
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
pub(crate) fn attach_pcurves(face: &Face, edges: &[Edge]) {
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
/// Delegates to [`crate::bop_same_domain_faces::fill_same_domain_faces_occt`],
/// the full `BOPAlgo_Builder::FillSameDomainFaces` translation
/// (`_2.cxx:580-925`): F/F + FaceInfo filter, parent-solid exclusion,
/// `BOPTools_Set` of edges, planar-bounded shortcut, `AreFacesSameDomain`,
/// `FillMap` / `MakeBlocks`, `myShapesSD` + image rewrite + `myOrigins`.
pub fn fill_same_domain_faces<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    crate::bop_same_domain_faces::fill_same_domain_faces_occt(f)
}
