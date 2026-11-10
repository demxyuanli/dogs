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
    /// Origins back-map (`myOrigins`) keyed by `TopTools_ShapeMapHasher`
    /// (`BOPAlgo_Builder.hxx:501`): `TShape` + `Location`.
    fn origins_mut(&mut self) -> &mut HashMap<ShapeMapKey, Vec<TopoShape>>;
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
/// Delegates to [`crate::bop_builder_split_faces::build_split_faces_occt`], the
/// full `BOPAlgo_Builder::BuildSplitFaces` translation (draft-face fast path,
/// `DoSplitSEAMOnFace`, `IsSplitToReverseWithWarn`, IN/Sc both ways).
pub fn build_split_faces<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    crate::bop_builder_split_faces::build_split_faces_occt(f)
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
        // `_2.cxx` takes that one `aPB->Edge()`; a neighbor whose Edge() is
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
