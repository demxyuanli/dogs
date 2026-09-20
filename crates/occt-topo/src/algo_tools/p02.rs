use super::prelude::*;
use super::*;

impl AlgoTools {

    /// `CorrectTolerances` with an avoid-map of source V/E/F keys.
    pub fn correct_tolerances_avoid(shape: &TopoShape, avoid: &HashSet<usize>, tol: f64) {
        for edge in edges_of(shape) {
            if avoid.contains(&GeometryRegistry::shape_key(&edge.0)) {
                continue;
            }
            let e_tol = BRepTool::edge_tolerance(&edge);
            let curve = BRepTool::edge_curve(&edge);
            let (a, b) = BRepTool::edge_parameters(&edge);
            let kids: Vec<TopoShape> = edge
                .0
                .tshape
                .read()
                .unwrap()
                .children
                .iter()
                .filter(|h| h.shape_type() == ShapeType::Vertex)
                .cloned()
                .collect();
            if let (Some(curve), true, true) = (curve, a.is_finite(), b.is_finite()) {
                let pa = curve.d0(a);
                let pb = curve.d0(b);
                for (i, p_end) in [(0usize, &pa), (1, &pb)] {
                    if let Some(k) = kids.get(i) {
                        if avoid.contains(&GeometryRegistry::shape_key(k)) {
                            continue;
                        }
                        let v = Vertex(k.clone());
                        let want = BRepTool::vertex_point(&v).distance(p_end) + D_TOLERANCE;
                        if want > BRepTool::vertex_tolerance(&v) && want <= tol {
                            v.set_tolerance(want);
                        }
                    }
                }
            }
            for k in &kids {
                if avoid.contains(&GeometryRegistry::shape_key(k)) {
                    continue;
                }
                let v = Vertex(k.clone());
                let cur = BRepTool::vertex_tolerance(&v);
                if e_tol > cur && e_tol <= tol {
                    v.set_tolerance(e_tol);
                }
            }
        }
    }

    /// Find a parameter of `edge` at which the edge's 3-D curve is *off* `face`
    /// — a sampled point whose UV classification on `face` is `Out`
    /// (`BOPTools_AlgoTools::GetEdgeOff`, geometric form).
    ///
    /// Returns `None` when the whole edge lies on or inside the face.
    pub fn get_edge_off(edge: &Edge, face: &Face) -> Option<f64> {
        let surf = BRepTool::face_surface(face)?;
        let cl = FClass2d::new(face, 1e-7).ok()?;
        let (a, b) = BRepTool::edge_parameters(edge);
        if !a.is_finite() || !b.is_finite() {
            return None;
        }
        for i in 0..=32 {
            let t = a + (b - a) * i as f64 / 32.0;
            let p = AlgoTools::point_on_edge(edge, t).ok()?;
            // Faithful `GeomAPI_ProjectPointOnSurf` (`Extrema_ExtPS`, T-52); the
            // previous 32x32 grid was A1's substitute. A sample whose projection
            // is not done contributes nothing (OCCT's classifier has no (u, v)
            // to work with there).
            let Some(ps) = occt_geom::geom_api::project_point_on_surface(surf.as_ref(), &p, 1e-7)
            else {
                continue;
            };
            if cl.perform(GpPnt2d::new(ps.u, ps.v)) == FaceState::Out {
                return Some(t);
            }
        }
        None
    }

    // ---- private helpers ---------------------------------------------------

    /// Distinct `TShape` identity keys of the edges of `shape`. A face yields
    /// the edges of its boundary wires; an edge yields itself; any other shape
    /// yields its distinct edge sub-shapes.
    pub(super) fn shape_edge_keys(shape: &TopoShape) -> Vec<usize> {
        match shape.shape_type() {
            ShapeType::Edge => vec![GeometryRegistry::shape_key(shape)],
            ShapeType::Face => {
                let f = Face(shape.clone());
                let mut out = Vec::new();
                for w in wires_of_face(&f) {
                    for e in edges_of_wire(&w) {
                        out.push(GeometryRegistry::shape_key(&e.0));
                    }
                }
                out
            }
            _ => edges_of(shape).into_iter().map(|e| GeometryRegistry::shape_key(&e.0)).collect(),
        }
    }

    /// Wrap a connected group into a [`ConnexityBlock`] and compute its
    /// regularity (every shared edge used by exactly two shapes).
    pub(super) fn connexity_block(group: Vec<TopoShape>) -> ConnexityBlock {
        let mut block = ConnexityBlock::new();
        block.change_shapes_mut().extend(group);
        let edge_sets: Vec<Vec<usize>> =
            block.shapes().iter().map(AlgoTools::shape_edge_keys).collect();
        let mut counts: HashMap<usize, usize> = HashMap::new();
        for es in &edge_sets {
            for e in es {
                *counts.entry(*e).or_insert(0) += 1;
            }
        }
        block.set_regular(counts.values().all(|&c| c == 2));
        block
    }

    /// All boundary edges of a face (its wires' edges).
    pub(super) fn face_edges(face: &Face) -> Vec<Edge> {
        let mut out = Vec::new();
        for w in wires_of_face(face) {
            for e in edges_of_wire(&w) {
                out.push(e);
            }
        }
        out
    }

    /// Unit surface normal of `face` at the parameter point closest to `p`.
    pub(super) fn face_normal_at_point(face: &Face, p: &GpPnt) -> Option<GpVec> {
        let surf = BRepTool::face_surface(face)?;
        // Faithful `GeomAPI_ProjectPointOnSurf` (`Extrema_ExtPS`, T-52).
        let (u, v) = occt_geom::geom_api::project_point_on_surface(surf.as_ref(), p, 1e-7)
            .map(|ps| (ps.u, ps.v))?;
        let n = surface_normal(surf.as_ref(), u, v);
        if n.square_magnitude() < 1e-30 { None } else { Some(n) }
    }
}
/// 2D (p-curve) helpers mirroring `BOPTools_AlgoTools2D`.
pub struct AlgoTools2D;

impl AlgoTools2D {
    /// The p-curve of `edge` on `face` — a convenience alias of
    /// [`AlgoTools::make_pcurve`] matching the OCCT `EdgeToFace` role of
    /// producing the edge's UV representation on a face.
    pub fn edge_to_face(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
        AlgoTools::make_pcurve(edge, face)
    }

    /// Bring `curve` inside the UV bounds of `face`
    /// (`BOPTools_AlgoTools2D::AdjustPCurveOnSurf`).
    ///
    /// Periodic dimensions are shifted by whole periods so the curve midpoint
    /// lands in range; a non-periodic dimension that still exits the bounds is
    /// trimmed to the in-bounds parameter subrange. Delegates to
    /// [`crate::pcurve_full::trim_pcurve_to_face`].
    pub fn adjust_pcurve_on_surf(
        curve: &Arc<dyn Curve2d>,
        face: &Face,
        tol: f64,
    ) -> Result<Arc<dyn Curve2d>, String> {
        pcurve_full::trim_pcurve_to_face(curve, face, tol)
    }
}

/// Shape-set utilities mirroring `BOPTools_Set`.
pub struct BOPToolsSet;

impl BOPToolsSet {
    /// Deduplicate `shapes`, keeping the first occurrence of each distinct
    /// shape. Two shapes are the same when they share the same `TShape` data
    /// (`TopoShape::same_tshape` — pointer identity), matching OCCT's shape-map
    /// hashing by `TopoDS_Shape`.
    pub fn shape_list(shapes: &[TopoShape]) -> Vec<TopoShape> {
        let mut out: Vec<TopoShape> = Vec::new();
        for s in shapes {
            if !out.iter().any(|o| o.same_tshape(s)) {
                out.push(s.clone());
            }
        }
        out
    }

    /// Count how many shapes in `shapes` have type `kind` (`TopAbs_ShapeEnum`).
    pub fn type_count(shapes: &[TopoShape], kind: ShapeType) -> usize {
        shapes.iter().filter(|s| s.shape_type() == kind).count()
    }
}
