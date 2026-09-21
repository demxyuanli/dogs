use super::prelude::*;
use super::*;

/// OCCT `BOPTools_AlgoTools::DTolerance()` — the extra tolerance added to new
/// or updated entities so their tolerance is *slightly* larger than the actual
/// distance, protecting against numerical instability.

pub const D_TOLERANCE: f64 = 1e-12;

/// Static build helpers mirroring `BOPTools_AlgoTools`.
pub struct AlgoTools;
impl AlgoTools {

    /// Make a vertex at the 3D point `p` with tolerance `tol`
    /// (`BOPTools_AlgoTools::MakeNewVertex(gp_Pnt, tol)`).
    ///
    /// The returned shape is a `Vertex` whose geometry is registered in the
    /// `GeometryRegistry` side-table (`BRep_TVertex`).
    pub fn make_new_vertex(p: &GpPnt, tol: f64) -> Result<TopoShape, String> {
        let v = TopoBuilder::new().make_vertex(*p, tol);
        Ok(v.into())
    }

    /// Make an edge from the intersection curve `curve` bounded by the two
    /// optional vertices `v1` / `v2` at parameters `t1` / `t2`
    /// (`BOPTools_AlgoTools::MakeEdge` + `MakeSectEdge`).
    ///
    /// * The endpoint vertex tolerances are raised to `tol + DTolerance()` so
    ///   each vertex covers the curve point at its parameter.
    /// * The edge's range is stored exactly as given (`(t1, t2)`), preserving
    ///   the intersection range. A reversed range (`t1 > t2`) adds `v1`
    ///   reversed / `v2` forward.
    /// * The edge's own tolerance is set to `tol`.
    pub fn make_edge(
        curve: Arc<dyn Curve>,
        v1: Option<&TopoShape>,
        t1: f64,
        v2: Option<&TopoShape>,
        t2: f64,
        tol: f64,
    ) -> Result<TopoShape, String> {
        let b = TopoBuilder::new();
        let need_tol = tol + D_TOLERANCE;
        if let Some(v) = v1 {
            Vertex(v.clone()).set_tolerance(need_tol);
        }
        if let Some(v) = v2 {
            Vertex(v.clone()).set_tolerance(need_tol);
        }

        let mut e = b.make_edge(curve, t1, t2);
        if let Some(v) = v1 {
            let o = if t1 < t2 { Orientation::Forward } else { Orientation::Reversed };
            b.add(&mut e.0, &v.oriented(o));
        }
        if let Some(v) = v2 {
            let o = if t1 < t2 { Orientation::Reversed } else { Orientation::Forward };
            b.add(&mut e.0, &v.oriented(o));
        }

        // `BRep_Builder::UpdateEdge(theE, theTolR3D)` — record the edge
        // tolerance in the side-table.
        if let Some(mut g) = GeometryRegistry::global().edge_geom(&e.0) {
            g.tolerance = tol;
            GeometryRegistry::global().set_edge(&e.0, g);
        }
        Ok(e.0)
    }

    /// Build the p-curve of `edge` on `face` (`BOPTools_AlgoTools2D::Make2D`).
    ///
    /// Delegates to [`crate::pcurve_full::make_pcurve_full`], which gives the
    /// analytic isoparametric projection on plane/cylinder/cone/sphere/torus
    /// faces and a sampled B-spline on general faces.
    pub fn make_pcurve(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
        pcurve_full::make_pcurve_full(edge, face)
    }

    /// Decide whether the vertex `v` coincides with the point `p` within the
    /// summed tolerance `tol_v + tol_p + Precision::Confusion()`
    /// (`BOPTools_AlgoTools::ComputeVV(vertex, point, tol)`).
    ///
    /// Returns `1` when the vertex interferes with the point (distance is not
    /// greater than the tolerance sum), `0` when the point is separated. This
    /// is the predicate convention used by the pave-filler builder in this
    /// port; the OCCT `ComputeVV` reports the inverse error code for the same
    /// distance comparison.
    pub fn compute_vv(v: &TopoShape, p: &GpPnt, tol: f64) -> i32 {
        let vtx = Vertex(v.clone());
        let tol_v = BRepTool::vertex_tolerance(&vtx);
        let tol_sum = tol_v + tol + CONFUSION;
        let tol_sum2 = tol_sum * tol_sum;
        let p1 = BRepTool::vertex_point(&vtx);
        let d2 = p1.square_distance(p);
        if d2 > tol_sum2 {
            0
        } else {
            1
        }
    }

    /// The 3D point of the edge's curve at parameter `t`
    /// (`BOPTools_AlgoTools::PointOnEdge`).
    ///
    /// Errors when the edge has no registered 3D curve.
    pub fn point_on_edge(edge: &Edge, t: f64) -> Result<GpPnt, String> {
        let curve =
            BRepTool::edge_curve(edge).ok_or("point_on_edge: edge has no 3D curve")?;
        Ok(curve.d0(t))
    }

    /// Make the vertex `v` cover the point of `edge` at parameter `t`
    /// (`BOPTools_AlgoTools::UpdateVertex(edge, t, vertex)`).
    ///
    /// Computes the distance between the vertex point and the edge point at
    /// `t`; when it exceeds the vertex's current tolerance, the tolerance is
    /// raised to `dist + DTolerance()`. The vertex position itself is left
    /// unchanged (the tolerance sphere is what grows to absorb the edge point).
    /// A no-op when the edge has no 3D curve.
    pub fn update_vertex(edge: &Edge, t: f64, v: &mut TopoShape) {
        let Some(curve) = BRepTool::edge_curve(edge) else { return };
        let vtx = Vertex(v.clone());
        let p_v = BRepTool::vertex_point(&vtx);
        let tol_v = BRepTool::vertex_tolerance(&vtx);
        let dist = p_v.distance(&curve.d0(t));
        if dist > tol_v {
            Vertex(v.clone()).set_tolerance(dist + D_TOLERANCE);
        }
    }

    /// Unit normal of the surface `s` at `(u, v)`
    /// (`BOPTools_AlgoTools3D::GetNormalToSurface`).
    ///
    /// Errors when the surface has a degenerate parameterization at that point
    /// (both partial derivatives collapse to zero).
    pub fn get_normal_to_surface(s: &dyn Surface, u: f64, v: f64) -> Result<GpVec, String> {
        let n = surface_normal(s, u, v);
        if n.square_magnitude() < 1e-30 {
            return Err(format!(
                "get_normal_to_surface: degenerate surface at (u={u}, v={v})"
            ));
        }
        Ok(n)
    }

    // ----------------------------------------------------------------------
    // Phase 18b — BOPTools_AlgoTools remainder
    // ----------------------------------------------------------------------

    /// Classify the 3-D point `p` relative to `shape`
    /// (`BOPTools_AlgoTools::ComputeState`).
    ///
    /// * `Vertex` — `On` when `p` is within `tol` of the vertex point, else `Out`;
    /// * `Edge`   — `On` when the distance from `p` to the edge is within `tol`,
    ///   else `Out` (a 1-D edge has no inside/outside region);
    /// * `Face`   — `p` is projected onto the face surface and classified in the
    ///   face's UV with [`crate::fclass2d::FClass2d`]; a projection whose 3-D
    ///   distance exceeds `tol` reports `Out`;
    /// * `Solid`  — [`crate::brep_class3d::SolidClassifier`] (`BRepClass3d`);
    /// * any container (`Compound`, `Wire`, `Shell`, …) — recurses into the
    ///   first child with a known state.
    ///
    /// The returned state reuses [`crate::fclass2d::FaceState`]
    /// (`TopAbs_State` semantics).
    pub fn compute_state(shape: &TopoShape, p: &GpPnt, tol: f64) -> Result<FaceState, String> {
        match shape.shape_type() {
            ShapeType::Vertex => {
                let v = Vertex(shape.clone());
                let d = BRepTool::vertex_point(&v).distance(p);
                Ok(if d <= tol { FaceState::On } else { FaceState::Out })
            }
            ShapeType::Edge => {
                let e = Edge(shape.clone());
                if BRepTool::edge_curve(&e).is_none() {
                    return Ok(FaceState::Out);
                }
                let (_, q) = closest_point_on_edge(&e, p, 32);
                Ok(if q.distance(p) <= tol { FaceState::On } else { FaceState::Out })
            }
            ShapeType::Face => {
                let f = Face(shape.clone());
                let surf = BRepTool::face_surface(&f)
                    .ok_or("compute_state: face has no registered surface")?;
                // Port of `BRepClass_FaceClassifier::Perform(F, P, Tol)`
                // (`BRepClass_FaceClassifier.cxx:76-125`): `Extrema_ExtPS` over
                // the face's `BRepTools::UVBounds` (`cxx:91-92`), the
                // `IsDone()`/`NbExt()` guards (`cxx:98-107`, leaving the state
                // `UNKNOWN` with `rejected = true`), the solution of smallest
                // square distance (`cxx:109-117`), then the 2D classifier on
                // that `(u, v)` (`cxx:121-123`).
                //
                // OCCT does **not** test the 3-D projection distance here: the
                // previous body's `surf.d0(u,v).distance(p) > tol → Out` gate
                // was port-invented (audit A4/T-40), and the `32×32` grid
                // search was A1's substitute. `BRepAdaptor_Surface(theF, false)`
                // (`cxx:90`) is the *unrestricted* face surface, so the window
                // comes from the explicit `UVBounds` argument, as here.
                let uv = crate::brep_uv_bounds::uv_box_of_face(&f);
                let ex = occt_geom::extrema_surf::ExtPs::with_window(
                    p,
                    surf.as_ref(),
                    uv.xmin(),
                    uv.xmax(),
                    uv.ymin(),
                    uv.ymax(),
                    tol,
                    tol,
                );
                if !ex.is_done() || ex.nb_ext() == 0 {
                    return Ok(FaceState::Unknown);
                }
                let mut best = 1usize;
                for i in 2..=ex.nb_ext() {
                    if ex.square_distance(i) < ex.square_distance(best) {
                        best = i;
                    }
                }
                let (u, v, _) = ex.point(best);
                let cl = FClass2d::new(&f, tol)?;
                Ok(cl.perform(GpPnt2d::new(u, v)))
            }
            ShapeType::Solid | ShapeType::CompSolid => {
                Ok(SolidClassifier::classify(shape, p, tol))
            }
            _ => {
                let kids = shape.tshape.read().unwrap().children.clone();
                for sub in kids {
                    let st = AlgoTools::compute_state(&sub, p, tol)?;
                    if st != FaceState::Unknown {
                        return Ok(st);
                    }
                }
                Ok(FaceState::Unknown)
            }
        }
    }

    /// Build a single connexity block starting from the first shape of
    /// `shapes`, following shared-edge adjacency
    /// (`BOPTools_AlgoTools::MakeConnexityBlock`).
    ///
    /// Two shapes are connected when they share an edge (`TShape` identity).
    /// The block's `shapes` are the connected component reachable from the
    /// first input shape; [`ConnexityBlock::is_regular`] is `true` when every
    /// shared edge is used by exactly two shapes of the block.
    pub fn make_connexity_block(shapes: &[TopoShape]) -> ConnexityBlock {
        if shapes.is_empty() {
            return ConnexityBlock::new();
        }
        let edge_sets: Vec<Vec<usize>> = shapes.iter().map(AlgoTools::shape_edge_keys).collect();
        let mut in_block = vec![false; shapes.len()];
        let mut stack = vec![0usize];
        in_block[0] = true;
        while let Some(i) = stack.pop() {
            for (j, ej) in edge_sets.iter().enumerate() {
                if in_block[j] {
                    continue;
                }
                if edge_sets[i].iter().any(|e| ej.contains(e)) {
                    in_block[j] = true;
                    stack.push(j);
                }
            }
        }
        let mut component: Vec<TopoShape> = Vec::new();
        for (i, s) in shapes.iter().enumerate() {
            if in_block[i] {
                component.push(s.clone());
            }
        }
        AlgoTools::connexity_block(component)
    }

    /// Group `shapes` into connexity blocks by shared-edge connectivity
    /// (`BOPTools_AlgoTools::MakeConnexityBlocks`).
    ///
    /// Two shapes are in the same block when they are connected through a chain
    /// of shapes pairwise sharing an edge. The returned blocks are ordered by
    /// the index of their first shape.
    pub fn make_connexity_blocks(shapes: &[TopoShape]) -> Vec<ConnexityBlock> {
        let n = shapes.len();
        if n == 0 {
            return Vec::new();
        }
        let edge_sets: Vec<Vec<usize>> = shapes.iter().map(AlgoTools::shape_edge_keys).collect();
        let mut parent: Vec<usize> = (0..n).collect();
        pub(super) fn find(parent: &mut [usize], x: usize) -> usize {
            if parent[x] != x {
                let r = find(parent, parent[x]);
                parent[x] = r;
            }
            parent[x]
        }
        pub(super) fn union(parent: &mut [usize], a: usize, b: usize) {
            let (ra, rb) = (find(parent, a), find(parent, b));
            if ra != rb {
                parent[ra] = rb;
            }
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if edge_sets[i].iter().any(|e| edge_sets[j].contains(e)) {
                    union(&mut parent, i, j);
                }
            }
        }
        let mut groups: HashMap<usize, Vec<TopoShape>> = HashMap::new();
        for i in 0..n {
            let r = find(&mut parent, i);
            groups.entry(r).or_default().push(shapes[i].clone());
        }
        let mut out: Vec<ConnexityBlock> = Vec::new();
        for i in 0..n {
            let r = find(&mut parent, i);
            if let Some(g) = groups.remove(&r) {
                out.push(AlgoTools::connexity_block(g));
            }
        }
        out
    }

    /// Splits the keys of the `adjacency` map into connected blocks.
    ///
    /// Port of `BOPAlgo_Tools::MakeBlocks` (BOPAlgo_Tools.hxx): each key of
    /// `adjacency` starts a chain; the chain is grown by repeatedly appending
    /// the neighbours of its members (a breadth-first walk over the symmetric
    /// adjacency map — the map `FillMap` produces), so every element belongs to
    /// exactly one block. Blocks are returned in the order of their first key
    /// (keys visited in ascending order), matching the OCCT insertion order.
    pub fn make_blocks<K>(adjacency: &HashMap<K, Vec<K>>) -> Vec<Vec<K>>
    where
        K: Clone + Eq + std::hash::Hash + Ord,
    {
        let mut keys: Vec<&K> = adjacency.keys().collect();
        keys.sort();
        let mut used: HashSet<K> = HashSet::new();
        let mut blocks: Vec<Vec<K>> = Vec::new();
        for &key in &keys {
            if !used.insert(key.clone()) {
                continue;
            }
            // Start the chain.
            let mut chain: Vec<K> = vec![key.clone()];
            // Grow it: the neighbours of every member join.
            let mut i = 0;
            while i < chain.len() {
                if let Some(neighbours) = adjacency.get(&chain[i]) {
                    for n in neighbours {
                        if used.insert(n.clone()) {
                            chain.push(n.clone());
                        }
                    }
                }
                i += 1;
            }
            blocks.push(chain);
        }
        blocks
    }

    /// Reorder the edges of `wire` so they form a closed chain: each edge's end
    /// vertex coincides with the next edge's start vertex
    /// (`BOPTools_AlgoTools::OrientEdgesOnWire`).
    ///
    /// The wire's child list is rewritten in place; edges that cannot be chained
    /// are appended at the tail in their original order. Edges are matched by
    /// their endpoint 3-D positions (the flat model drops per-edge orientation,
    /// so the chain is defined geometrically).
    pub fn orient_edges_on_wire(wire: &mut TopoShape) {
        let edges: Vec<TopoShape> = wire
            .tshape
            .read()
            .unwrap()
            .children
            .iter()
            .filter(|h| h.shape_type() == ShapeType::Edge)
            .cloned()
            .collect();
        let mut unique: Vec<TopoShape> = Vec::new();
        for e in edges {
            // Collapse only exact duplicate references (same TShape *and* same
            // orientation). A seam edge of a periodic surface appears twice with
            // opposite orientations (its two sides); both must be kept or the
            // wire's UV polygon degenerates to one seam side.
            if !unique.iter().any(|u| u.same_tshape(&e) && u.orientation() == e.orientation()) {
                unique.push(e);
            }
        }
        if unique.is_empty() {
            return;
        }
        let endpoints: Vec<(GpPnt, GpPnt)> = unique
            .iter()
            .map(|e| match BRepTool::edge_vertices(&Edge(e.clone())) {
                Some((a, b)) => (a, b),
                None => (GpPnt::zero(), GpPnt::zero()),
            })
            .collect();
        let n = unique.len();
        let mut order: Vec<usize> = Vec::with_capacity(n);
        let mut used = vec![false; n];
        order.push(0);
        used[0] = true;
        let mut prev_end = endpoints[0].1;
        for _ in 1..n {
            let mut next: Option<(usize, bool)> = None; // (index, reversed?)
            for j in 0..n {
                if used[j] {
                    continue;
                }
                if endpoints[j].0.distance(&prev_end) <= 1e-9 {
                    next = Some((j, false));
                    break;
                }
                if endpoints[j].1.distance(&prev_end) <= 1e-9 {
                    next = Some((j, true));
                    break;
                }
            }
            match next {
                Some((j, rev)) => {
                    used[j] = true;
                    order.push(j);
                    prev_end = if rev { endpoints[j].0 } else { endpoints[j].1 };
                }
                None => break,
            }
        }
        for j in 0..n {
            if !used[j] {
                order.push(j);
            }
        }
        let ordered: Vec<TopoShape> =
            order.into_iter().map(|i| unique[i].clone()).collect();
        if let Ok(mut t) = wire.tshape.write() {
            t.children = ordered;
        }
    }

    /// Rebuild `shell` from its distinct faces, ordered by adjacency along
    /// shared edges (`BOPTools_AlgoTools::OrientFacesOnShell`).
    ///
    /// A face sharing an edge with an already-placed face is added next; the
    /// remaining faces follow in their original order. Duplicate references to
    /// the same face (e.g. a seam edge counted twice) collapse to one.
    pub fn orient_faces_on_shell(shell: &mut TopoShape) {
        let faces: Vec<TopoShape> = shell
            .tshape
            .read()
            .unwrap()
            .children
            .iter()
            .filter(|h| h.shape_type() == ShapeType::Face)
            .cloned()
            .collect();
        let mut unique: Vec<TopoShape> = Vec::new();
        for f in faces {
            if !unique.iter().any(|u| u.same_tshape(&f)) {
                unique.push(f);
            }
        }
        let n = unique.len();
        let edge_sets: Vec<Vec<usize>> = unique.iter().map(AlgoTools::shape_edge_keys).collect();
        let mut order: Vec<usize> = Vec::with_capacity(n);
        let mut used = vec![false; n];
        if n > 0 {
            order.push(0);
            used[0] = true;
        }
        while order.len() < n {
            let mut progressed = false;
            for i in 0..n {
                if used[i] {
                    continue;
                }
                if order
                    .iter()
                    .any(|&k| edge_sets[k].iter().any(|e| edge_sets[i].contains(e)))
                {
                    used[i] = true;
                    order.push(i);
                    progressed = true;
                }
            }
            if !progressed {
                break;
            }
        }
        for i in 0..n {
            if !used[i] {
                order.push(i);
            }
        }
        let ordered: Vec<TopoShape> =
            order.into_iter().map(|i| unique[i].clone()).collect();
        if let Ok(mut t) = shell.tshape.write() {
            t.children = ordered;
        }
        crate::brep_class3d::apply_orient_faces_on_shell(shell);
    }

    /// Copy `edge` into a new edge sharing the same 3-D curve, parameter range,
    /// tolerance, p-curves and vertex children
    /// (`BOPTools_AlgoTools::CopyEdge`).
    pub fn copy_edge(edge: &Edge) -> Result<Edge, String> {
        let reg = GeometryRegistry::global();
        let geom = reg
            .edge_geom(&edge.0)
            .ok_or("copy_edge: edge has no registered geometry")?;
        let b = TopoBuilder::new();
        let mut e = b.make_edge(geom.curve.clone(), geom.first, geom.last);
        if let Some(mut g) = reg.edge_geom(&e.0) {
            g.tolerance = geom.tolerance;
            g.same_parameter = geom.same_parameter;
            g.same_range = geom.same_range;
            g.degenerated = geom.degenerated;
            g.pcurves = geom.pcurves.clone();
            reg.set_edge(&e.0, g);
        }
        let kids = edge.0.tshape.read().unwrap().children.clone();
        for k in kids {
            b.add(&mut e.0, &k);
        }
        Ok(e)
    }

    /// Make a new edge from the base `edge`'s curve bounded by the optional
    /// vertices `v1` / `v2` at parameters `t1` / `t2`
    /// (`BOPTools_AlgoTools::MakeSplitEdge`).
    ///
    /// The stored range is `(min(t1,t2), max(t1,t2))`; a reversed parameter
    /// order adds `v1` reversed / `v2` forward, mirroring the OCCT convention.
    pub fn make_split_edge(
        edge: &Edge,
        v1: Option<&TopoShape>,
        t1: f64,
        v2: Option<&TopoShape>,
        t2: f64,
    ) -> Result<Edge, String> {
        let reg = GeometryRegistry::global();
        let geom = reg
            .edge_geom(&edge.0)
            .ok_or("make_split_edge: edge has no registered geometry")?;
        let (lo, hi) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
        let b = TopoBuilder::new();
        let mut e = b.make_edge(geom.curve.clone(), lo, hi);
        if let Some(mut g) = reg.edge_geom(&e.0) {
            g.tolerance = geom.tolerance;
            g.same_parameter = geom.same_parameter;
            g.same_range = geom.same_range;
            g.degenerated = geom.degenerated;
            reg.set_edge(&e.0, g);
        }
        if let Some(v) = v1 {
            let o = if t1 < t2 { Orientation::Forward } else { Orientation::Reversed };
            b.add(&mut e.0, &v.oriented(o));
        }
        if let Some(v) = v2 {
            let o = if t1 < t2 { Orientation::Reversed } else { Orientation::Forward };
            b.add(&mut e.0, &v.oriented(o));
        }
        Ok(e)
    }

    /// Whether `edge` is a micro-edge: a degenerated edge or an edge whose
    /// 3-D curve length is below `tol`
    /// (`BOPTools_AlgoTools::IsMicroEdge`).
    pub fn is_micro_edge(edge: &Edge, tol: f64) -> bool {
        if BRepTool::is_degenerated(edge) {
            return true;
        }
        let (a, b) = BRepTool::edge_parameters(edge);
        if !a.is_finite() || !b.is_finite() || (a - b).abs() < 1e-15 {
            return true;
        }
        edge_length(edge, 32) < tol
    }

    /// Whether `shell` is open: it has a non-degenerated boundary edge shared
    /// by exactly one face (`BOPTools_AlgoTools::IsOpenShell`). A `Solid` is
    /// open when any of its shells is open.
    pub fn is_open_shell(shell: &TopoShape) -> bool {
        let shells: Vec<TopoShape> = match shell.shape_type() {
            ShapeType::Shell => vec![shell.clone()],
            ShapeType::Solid => shell
                .tshape
                .read()
                .unwrap()
                .children
                .iter()
                .filter(|h| h.shape_type() == ShapeType::Shell)
                .cloned()
                .collect(),
            _ => return false,
        };
        shells.iter().any(|sh| {
            let faces: Vec<TopoShape> = sh
                .tshape
                .read()
                .unwrap()
                .children
                .iter()
                .filter(|h| h.shape_type() == ShapeType::Face)
                .cloned()
                .collect();
            let mut counts: HashMap<usize, usize> = HashMap::new();
            for f in &faces {
                for e in AlgoTools::shape_edge_keys(f) {
                    *counts.entry(e).or_insert(0) += 1;
                }
            }
            counts.values().any(|&c| c == 1)
        })
    }

    /// `BOPTools_AlgoTools::IsInvertedSolid`: the infinite point classifies as
    /// `In` (`BRepClass3d_SolidClassifier::PerformInfinitePoint` at `1.e-7`).
    pub fn is_inverted_solid(solid: &TopoShape) -> bool {
        if solid.shape_type() != ShapeType::Solid && solid.shape_type() != ShapeType::CompSolid {
            return false;
        }
        let mut sc = SolidClassifier::new();
        sc.load(solid.clone());
        sc.perform_infinite_point(1e-7);
        sc.state() == FaceState::In
    }

    /// Sense of two faces sharing an edge: `1` when their normals near the
    /// shared edge point the same way, `-1` when opposite, `0` when no shared
    /// edge exists or the normals are not collinear
    /// (`BOPTools_AlgoTools::Sense`).
    pub fn sense(f1: &Face, f2: &Face) -> i32 {
        let e1s = AlgoTools::face_edges(f1);
        let e2s = AlgoTools::face_edges(f2);
        let shared = e1s.iter().find(|e| {
            !BRepTool::is_degenerated(e) && e2s.iter().any(|e2| e2.same_tshape(&e.0))
        });
        let Some(e) = shared else { return 0 };
        let (a, b) = BRepTool::edge_parameters(e);
        if !a.is_finite() || !b.is_finite() {
            return 0;
        }
        let t = 0.5 * (a + b);
        let p = BRepTool::edge_curve(e).map(|c| c.d0(t)).unwrap_or(GpPnt::zero());
        let (Some(n1), Some(n2)) = (
            AlgoTools::face_normal_at_point(f1, &p),
            AlgoTools::face_normal_at_point(f2, &p),
        ) else {
            return 0;
        };
        let dot = n1.xyz().dot(n2.xyz());
        if dot > 1e-9 {
            1
        } else if dot < -1e-9 {
            -1
        } else {
            0
        }
    }

    /// Whether the wire `w` is a hole of `face`: its UV pcurve winds with
    /// positive signed area on the face surface
    /// (`BOPTools_AlgoTools::IsHole`).
    pub fn is_hole(w: &TopoShape, face: &Face) -> bool {
        let mut s = 0.0;
        for e in edges_of_wire(&Wire(w.clone())) {
            let or = e.orientation();
            if or != Orientation::Forward && or != Orientation::Reversed {
                continue;
            }
            let (a, b) = BRepTool::edge_parameters(&e);
            if !a.is_finite() || !b.is_finite() {
                continue;
            }
            let Ok(pc) = AlgoTools::make_pcurve(&e, face) else { continue };
            let dir = if or == Orientation::Reversed { -1.0 } else { 1.0 };
            let mut prev = if dir > 0.0 { pc.d0(a) } else { pc.d0(b) };
            for i in 1..=16 {
                let t = a + (b - a) * i as f64 / 16.0;
                let cur = pc.d0(t);
                s += (prev.y() + cur.y()) * (cur.x() - prev.x()) * dir;
                prev = cur;
            }
        }
        s > 0.0
    }

    /// Dimension of `shape`: `0` for a vertex, `1` for an edge or wire, `2` for
    /// a face or shell, `3` for a solid or comp-solid
    /// (`BOPTools_AlgoTools::Dimension`). A compound reports the maximum
    /// dimension of its children.
    pub fn dimension(shape: &TopoShape) -> usize {
        match shape.shape_type() {
            ShapeType::Vertex => 0,
            ShapeType::Edge | ShapeType::Wire => 1,
            ShapeType::Face | ShapeType::Shell => 2,
            ShapeType::Solid | ShapeType::CompSolid => 3,
            ShapeType::Compound => {
                let kids = shape.tshape.read().unwrap().children.clone();
                kids.iter()
                    .map(|h| AlgoTools::dimension(h))
                    .max()
                    .unwrap_or(0)
            }
            ShapeType::Shape => 0,
        }
    }

    /// Grow insufficient vertex tolerances of `shape` so every vertex on an
    /// edge covers the edge's curve endpoints and at least the edge tolerance
    /// (`BOPTools_AlgoTools::CorrectTolerances` / `CorrectShapeTolerances`).
    ///
    /// Shapes whose TShape key is in `avoid` are skipped (`MapToAvoid` / `aMA`
    /// in `BOPAlgo_Builder::PostTreat` when non-destructive).
    ///
    /// Tolerances are only *raised*, never lowered, and never exceed `tol`.
    pub fn correct_tolerances(shape: &TopoShape, tol: f64) {
        Self::correct_tolerances_avoid(shape, &HashSet::new(), tol);
    }
}
