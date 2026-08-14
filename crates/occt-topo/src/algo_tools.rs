//! Build-time geometry helpers for the boolean PaveFiller — a port of the
//! `BOPTools_AlgoTools` / `BOPTools_AlgoTools2D` / `BOPTools_AlgoTools3D` /
//! `BOPTools_Set` utility families (TKBO).
//!
//! These are the low-level construction primitives the pave-filler uses while
//! building a boolean result:
//!
//! * **Vertex construction** — [`AlgoTools::make_new_vertex`] builds a vertex
//!   at a 3D point with a tolerance;
//! * **Edge construction** — [`AlgoTools::make_edge`] builds an edge from an
//!   intersection curve plus two endpoint vertices and parameters, bumping the
//!   vertex tolerances so they cover the curve ends (`BOPTools_AlgoTools::MakeEdge`);
//! * **Vertex coincidence** — [`AlgoTools::compute_vv`] decides whether a vertex
//!   and a point interfere within the summed tolerances (`ComputeVV`);
//! * **Point-on-edge** — [`AlgoTools::point_on_edge`] evaluates the edge's 3D
//!   curve at a parameter; [`AlgoTools::update_vertex`] grows a vertex's
//!   tolerance to cover such a point (`UpdateVertex`);
//! * **P-curves** — [`AlgoTools::make_pcurve`] / [`AlgoTools2D::edge_to_face`]
//!   build the edge→face UV curve by delegating to
//!   [`crate::pcurve_full::make_pcurve_full`];
//!   [`AlgoTools2D::adjust_pcurve_on_surf`] brings a pcurve inside the face UV
//!   bounds via [`crate::pcurve_full::trim_pcurve_to_face`]
//!   (`AdjustPCurveOnSurf`);
//! * **Surface normal** — [`AlgoTools::get_normal_to_surface`] computes the
//!   unit normal of a surface at `(u, v)` (`GetNormalToSurface`);
//! * **Shape sets** — [`BOPToolsSet::shape_list`] dedupes a shape collection
//!   and [`BOPToolsSet::type_count`] counts shapes of a given type.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_extrema::{closest_point_on_edge, is_inside, point_shape_distance};
use crate::brep_measure::edge_length;
use crate::brep_surface::{surface_closest_params, surface_normal};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::connexity_block::ConnexityBlock;
use crate::fclass2d::{FClass2d, FaceState};
use crate::pcurve_full;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, wires_of_face};
use crate::tshape::HandleTShape;

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
    /// * `Solid`  — `On` when `p` is within `tol` of the boundary, otherwise the
    ///   parity test [`crate::brep_extrema::is_inside`];
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
                let (u, v) = surface_closest_params(surf.as_ref(), p, 32, 32);
                if surf.d0(u, v).distance(p) > tol {
                    return Ok(FaceState::Out);
                }
                let cl = FClass2d::new(&f, tol)?;
                Ok(cl.perform(GpPnt2d::new(u, v)))
            }
            ShapeType::Solid | ShapeType::CompSolid => {
                if point_shape_distance(p, shape) <= tol {
                    return Ok(FaceState::On);
                }
                Ok(if is_inside(shape, p) { FaceState::In } else { FaceState::Out })
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
        fn find(parent: &mut [usize], x: usize) -> usize {
            if parent[x] != x {
                let r = find(parent, parent[x]);
                parent[x] = r;
            }
            parent[x]
        }
        fn union(parent: &mut [usize], a: usize, b: usize) {
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
            if !unique.iter().any(|u| u.same_tshape(&e)) {
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

    /// Whether `solid` is inverted: the natural surface normals of its boundary
    /// faces point inward, i.e. the signed boundary volume is negative. This is
    /// the 3-D analogue of `BOPTools_AlgoTools::IsInvertedSolid`, which
    /// classifies the infinite point as `In`.
    pub fn is_inverted_solid(solid: &TopoShape) -> bool {
        if solid.shape_type() != ShapeType::Solid && solid.shape_type() != ShapeType::CompSolid {
            return false;
        }
        let mut signed_vol = 0.0;
        for face in faces_of(solid) {
            let Some(surf) = BRepTool::face_surface(&face) else { continue };
            let (u1, u2, v1, v2) = BRepTool::uv_bounds(&face);
            let (u1, u2, v1, v2) = (
                if u1.is_finite() { u1 } else { -1.0 },
                if u2.is_finite() { u2 } else { 1.0 },
                if v1.is_finite() { v1 } else { -1.0 },
                if v2.is_finite() { v2 } else { 1.0 },
            );
            const N: usize = 8;
            for i in 0..N {
                for j in 0..N {
                    let ua = u1 + (u2 - u1) * i as f64 / N as f64;
                    let ub = u1 + (u2 - u1) * (i + 1) as f64 / N as f64;
                    let va = v1 + (v2 - v1) * j as f64 / N as f64;
                    let vb = v1 + (v2 - v1) * (j + 1) as f64 / N as f64;
                    let p00 = surf.d0(ua, va);
                    let p10 = surf.d0(ub, va);
                    let p01 = surf.d0(ua, vb);
                    let p11 = surf.d0(ub, vb);
                    let n_surf = surface_normal(surf.as_ref(), 0.5 * (ua + ub), 0.5 * (va + vb));
                    signed_vol += AlgoTools::tet_signed_volume(&p00, &p10, &p01, &n_surf);
                    signed_vol += AlgoTools::tet_signed_volume(&p00, &p01, &p11, &n_surf);
                }
            }
        }
        signed_vol < 0.0
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
    /// Tolerances are only *raised*, never lowered, and never exceed `tol`.
    pub fn correct_tolerances(shape: &TopoShape, tol: f64) {
        for edge in edges_of(shape) {
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
                        let v = Vertex(k.clone());
                        let want = BRepTool::vertex_point(&v).distance(p_end) + D_TOLERANCE;
                        if want > BRepTool::vertex_tolerance(&v) && want <= tol {
                            v.set_tolerance(want);
                        }
                    }
                }
            }
            for k in &kids {
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
            let (u, v) = surface_closest_params(surf.as_ref(), &p, 32, 32);
            if cl.perform(GpPnt2d::new(u, v)) == FaceState::Out {
                return Some(t);
            }
        }
        None
    }

    // ---- private helpers ---------------------------------------------------

    /// Distinct `TShape` identity keys of the edges of `shape`. A face yields
    /// the edges of its boundary wires; an edge yields itself; any other shape
    /// yields its distinct edge sub-shapes.
    fn shape_edge_keys(shape: &TopoShape) -> Vec<usize> {
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
    fn connexity_block(group: Vec<TopoShape>) -> ConnexityBlock {
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
    fn face_edges(face: &Face) -> Vec<Edge> {
        let mut out = Vec::new();
        for w in wires_of_face(face) {
            for e in edges_of_wire(&w) {
                out.push(e);
            }
        }
        out
    }

    /// Unit surface normal of `face` at the parameter point closest to `p`.
    fn face_normal_at_point(face: &Face, p: &GpPnt) -> Option<GpVec> {
        let surf = BRepTool::face_surface(face)?;
        let (u, v) = surface_closest_params(surf.as_ref(), p, 32, 32);
        let n = surface_normal(surf.as_ref(), u, v);
        if n.square_magnitude() < 1e-30 { None } else { Some(n) }
    }

    /// Signed tetrahedron volume `(1/6)·p0·(p1×p2)` of a triangle `(p0,p1,p2)`,
    /// sign-corrected so the triangle winds with the surface normal `n_surf`.
    fn tet_signed_volume(p0: &GpPnt, p1: &GpPnt, p2: &GpPnt, n_surf: &GpVec) -> f64 {
        let e1 = GpVec::from_pnts(p0, p1);
        let e2 = GpVec::from_pnts(p0, p2);
        let tri_n = e1.xyz().crossed(e2.xyz());
        let sign = if tri_n.dot(n_surf.xyz()) >= 0.0 { 1.0 } else { -1.0 };
        sign * p0.coord.dot(&tri_n) / 6.0
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::builder_face::build_face_with_holes;
    use crate::topo_tools_full::wire_is_closed;
    use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt};
    use occt_geom::GeomLine;

    #[test]
    fn make_new_vertex_at_box_corner() {
        let b = unit_box();
        let v = AlgoTools::make_new_vertex(&b.corners[0], 1e-7).expect("vertex");
        assert!(v.is_vertex());
        let vtx = Vertex(v);
        assert!(BRepTool::vertex_point(&vtx).is_equal(&b.corners[0]), "point must match corner");
        assert_eq!(BRepTool::vertex_tolerance(&vtx), 1e-7);
    }

    #[test]
    fn make_edge_from_box_edge_curve_connects_endpoints() {
        let b = unit_box();
        let curve = BRepTool::edge_curve(&b.edges[0]).expect("box edge curve");
        let e = AlgoTools::make_edge(
            curve,
            Some(&b.vertices[0].0),
            0.0,
            Some(&b.vertices[1].0),
            1.0,
            1e-7,
        )
        .expect("edge");
        assert!(e.is_edge());
        assert_eq!(e.tshape.read().unwrap().children.len(), 2);

        let p0 = AlgoTools::point_on_edge(&Edge(e.clone()), 0.0).expect("start point");
        let p1 = AlgoTools::point_on_edge(&Edge(e.clone()), 1.0).expect("end point");
        assert!(p0.is_equal(&b.corners[0]), "start {:?}", p0);
        assert!(p1.is_equal(&b.corners[1]), "end {:?}", p1);
        assert_eq!(BRepTool::edge_tolerance(&Edge(e)), 1e-7);

        // Vertex tolerances were bumped to cover the curve ends.
        assert!(
            BRepTool::vertex_tolerance(&b.vertices[0]) >= 1e-7,
            "v0 tolerance grew"
        );
    }

    #[test]
    fn make_pcurve_agrees_with_pcurve_full() {
        let b = unit_box();
        let edge = &b.edges[0];
        let face = &b.faces[0];
        let a = AlgoTools::make_pcurve(edge, face).expect("make_pcurve");
        let full = pcurve_full::make_pcurve_full(edge, face).expect("pcurve_full");
        for i in 0..=8 {
            let t = i as f64 / 8.0;
            let pa = a.d0(t);
            let pf = full.d0(t);
            assert!(
                (pa.x() - pf.x()).abs() < 1e-9 && (pa.y() - pf.y()).abs() < 1e-9,
                "sample t={t}"
            );
        }
    }

    #[test]
    fn compute_vv_coincident_returns_one_separated_returns_zero() {
        let b = unit_box();
        // The box corner vertex coincides with its own point → 1.
        assert_eq!(AlgoTools::compute_vv(&b.vertices[0].0, &b.corners[0], 1e-7), 1);
        // A far-away point is separated → 0.
        assert_eq!(
            AlgoTools::compute_vv(&b.vertices[0].0, &GpPnt::new(10.0, 10.0, 10.0), 1e-7),
            0
        );
    }

    #[test]
    fn point_on_edge_at_endpoints() {
        let b = unit_box();
        let p0 = AlgoTools::point_on_edge(&b.edges[0], 0.0).expect("t=0");
        let p1 = AlgoTools::point_on_edge(&b.edges[0], 1.0).expect("t=1");
        assert!(p0.is_equal(&b.corners[0]), "t=0 {:?}", p0);
        assert!(p1.is_equal(&b.corners[1]), "t=1 {:?}", p1);
    }

    #[test]
    fn update_vertex_grows_tolerance_to_cover_edge_point() {
        let b = unit_box();
        // Vertex near (0,0,0.5): edge 0 runs (0,0,0)→(1,0,0), so the point at
        // t=0.5 is (0.5,0,0) — distance ≈ 0.707 → tolerance must grow above it.
        let mut v = TopoBuilder::new().make_vertex(GpPnt::new(0.0, 0.0, 0.5), 1e-7);
        AlgoTools::update_vertex(&b.edges[0], 0.5, &mut v.0);
        let tol = BRepTool::vertex_tolerance(&v);
        assert!(tol > 0.7, "tolerance {tol} must cover the edge point");

        // A vertex already covering the edge point is left untouched.
        let mut v2 = TopoBuilder::new().make_vertex(GpPnt::new(0.5, 0.0, 0.0), 0.1);
        AlgoTools::update_vertex(&b.edges[0], 0.5, &mut v2.0);
        assert_eq!(BRepTool::vertex_tolerance(&v2), 0.1);
    }

    #[test]
    fn get_normal_to_surface_plane_is_z() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let s = BRepTool::face_surface(&face).expect("plane surface");
        let n = AlgoTools::get_normal_to_surface(s.as_ref(), 0.5, 0.5).expect("normal");
        // Standard plane: unit normal along ±Z.
        assert!(n.xyz().z.abs() > 0.99, "normal {:?}", n);
        assert!(n.square_magnitude() - 1.0 < 1e-6, "unit normal");
    }

    #[test]
    fn edge_to_face_and_adjust_pcurve_on_surf() {
        let b = unit_box();
        let edge = &b.edges[0];
        let face = &b.faces[0];
        let pc = AlgoTools2D::edge_to_face(edge, face).expect("edge to face");
        // Edge 0 (0,0,0)→(1,0,0) on the bottom face maps to a (0, v) isoline.
        let q0 = pc.d0(0.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {:?}", q0);

        // Adjusting keeps the pcurve inside the face UV bounds.
        let adj = AlgoTools2D::adjust_pcurve_on_surf(&pc, face, 1e-7).expect("adjusted");
        let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(face);
        for i in 0..=8 {
            let q = adj.d0(i as f64 / 8.0);
            assert!(q.x() >= umin - 1e-6 && q.x() <= umax + 1e-6, "u {} at i={}", q.x(), i);
            assert!(q.y() >= vmin - 1e-6 && q.y() <= vmax + 1e-6, "v {} at i={}", q.y(), i);
        }
    }

    #[test]
    fn shape_list_dedupes_and_type_count() {
        let b = unit_box();
        let faces: Vec<TopoShape> = b.faces.iter().map(|f| f.0.clone()).collect();
        assert_eq!(BOPToolsSet::type_count(&faces, ShapeType::Face), 6);
        assert_eq!(BOPToolsSet::type_count(&faces, ShapeType::Edge), 0);

        // Duplicated faces collapse to the 6 distinct ones.
        let mut dup = faces.clone();
        dup.extend(faces.iter().cloned());
        assert_eq!(BOPToolsSet::shape_list(&dup).len(), 6);

        // Repeated identical edges collapse to one.
        let e0 = b.edges[0].0.clone();
        let lst = vec![e0.clone(), e0.clone(), b.edges[1].0.clone()];
        assert_eq!(BOPToolsSet::shape_list(&lst).len(), 2);
    }

    // ---- Phase 18b tests ----------------------------------------------------

    #[test]
    fn compute_state_face_in_out_on() {
        let b = unit_box();
        let face = &b.faces[0]; // bottom (z = 0)
        // Interior point of the face → In.
        assert_eq!(
            AlgoTools::compute_state(&face.0, &GpPnt::new(0.5, 0.5, 0.0), 1e-6).unwrap(),
            FaceState::In
        );
        // Point above the face plane → Out.
        assert_eq!(
            AlgoTools::compute_state(&face.0, &GpPnt::new(0.5, 0.5, 0.1), 1e-6).unwrap(),
            FaceState::Out
        );
        // Point on the face boundary edge → On.
        assert_eq!(
            AlgoTools::compute_state(&face.0, &GpPnt::new(0.5, 0.0, 0.0), 1e-6).unwrap(),
            FaceState::On
        );
    }

    #[test]
    fn compute_state_solid_in_out_on() {
        let b = unit_box();
        assert_eq!(
            AlgoTools::compute_state(&b.solid.0, &GpPnt::new(0.5, 0.5, 0.5), 1e-6).unwrap(),
            FaceState::In
        );
        assert_eq!(
            AlgoTools::compute_state(&b.solid.0, &GpPnt::new(2.0, 0.0, 0.0), 1e-6).unwrap(),
            FaceState::Out
        );
        // On the front face (y = 0).
        assert_eq!(
            AlgoTools::compute_state(&b.solid.0, &GpPnt::new(0.5, 0.0, 0.5), 1e-6).unwrap(),
            FaceState::On
        );
    }

    #[test]
    fn compute_state_edge_and_vertex() {
        let b = unit_box();
        // Edge 0 runs (0,0,0)→(1,0,0).
        assert_eq!(
            AlgoTools::compute_state(&b.edges[0].0, &GpPnt::new(0.5, 0.0, 0.0), 1e-6).unwrap(),
            FaceState::On
        );
        assert_eq!(
            AlgoTools::compute_state(&b.edges[0].0, &GpPnt::new(0.5, 1.0, 0.0), 1e-6).unwrap(),
            FaceState::Out
        );
        // Vertex at a corner: coincident → On, far away → Out.
        assert_eq!(
            AlgoTools::compute_state(&b.vertices[0].0, &b.corners[0], 1e-6).unwrap(),
            FaceState::On
        );
        assert_eq!(
            AlgoTools::compute_state(&b.vertices[0].0, &GpPnt::new(5.0, 5.0, 5.0), 1e-6).unwrap(),
            FaceState::Out
        );
    }

    #[test]
    fn make_connexity_block_box_faces_single() {
        let b = unit_box();
        let faces: Vec<TopoShape> = b.faces.iter().map(|f| f.0.clone()).collect();
        let block = AlgoTools::make_connexity_block(&faces);
        assert_eq!(block.shapes().len(), 6);
        assert!(block.is_regular());
    }

    #[test]
    fn make_connexity_blocks_box_one_block() {
        let b = unit_box();
        let faces: Vec<TopoShape> = b.faces.iter().map(|f| f.0.clone()).collect();
        let blocks = AlgoTools::make_connexity_blocks(&faces);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].shapes().len(), 6);
    }

    #[test]
    fn make_connexity_blocks_separates_disjoint_faces() {
        let tb = TopoBuilder::new();
        let sq = |tb: &TopoBuilder, off: f64| -> TopoShape {
            let pts = [
                GpPnt::new(off, 0.0, 0.0),
                GpPnt::new(off + 1.0, 0.0, 0.0),
                GpPnt::new(off + 1.0, 1.0, 0.0),
                GpPnt::new(off, 1.0, 0.0),
            ];
            let edges: Vec<Edge> = (0..4)
                .map(|i| tb.make_edge_segment(&pts[i], &pts[(i + 1) % 4]))
                .collect();
            build_face_with_holes(&edges, &[]).expect("square face").0
        };
        let f1 = sq(&tb, 0.0);
        let f2 = sq(&tb, 10.0);
        let blocks = AlgoTools::make_connexity_blocks(&[f1, f2]);
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].shapes().len(), 1);
        assert_eq!(blocks[1].shapes().len(), 1);
    }

    #[test]
    fn make_blocks_connects_transitive_chain() {
        // 0~1, 1~2 (no direct 0~2): MakeBlocks joins the chain into one block.
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        adj.entry(0).or_default().extend([1]);
        adj.entry(1).or_default().extend([0, 2]);
        adj.entry(2).or_default().extend([1]);
        let blocks = AlgoTools::make_blocks(&adj);
        assert_eq!(blocks.len(), 1, "transitive chain must collapse into one block");
        assert_eq!(blocks[0].len(), 3);
        assert!(blocks[0].contains(&0) && blocks[0].contains(&1) && blocks[0].contains(&2));
    }

    #[test]
    fn make_blocks_separates_disconnected() {
        // Two disconnected pairs and one isolated key.
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        adj.entry(0).or_default().extend([1]);
        adj.entry(1).or_default().extend([0]);
        adj.entry(2).or_default().extend([3]);
        adj.entry(3).or_default().extend([2]);
        adj.entry(4).or_default(); // isolated: present as a key, no neighbours
        let blocks = AlgoTools::make_blocks(&adj);
        assert_eq!(blocks.len(), 3);
        let mut sizes: Vec<usize> = blocks.iter().map(|b| b.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![1, 2, 2]);
    }

    #[test]
    fn orient_edges_on_wire_chains_shuffled_square() {
        let tb = TopoBuilder::new();
        let p = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let e1 = tb.make_edge_segment(&p[0], &p[1]);
        let e2 = tb.make_edge_segment(&p[1], &p[2]);
        let e3 = tb.make_edge_segment(&p[2], &p[3]);
        let e4 = tb.make_edge_segment(&p[3], &p[0]);
        let mut wire = tb.make_wire(&[e2, e4, e1, e3]);
        assert!(!wire_is_closed(&wire));
        AlgoTools::orient_edges_on_wire(&mut wire.0);
        assert!(wire_is_closed(&wire));
    }

    #[test]
    fn orient_faces_on_shell_dedupes_and_orders() {
        let b = unit_box();
        let shell = TopoBuilder::new().make_shell(&b.faces);
        let mut sh = shell.0;
        AlgoTools::orient_faces_on_shell(&mut sh);
        assert_eq!(sh.tshape.read().unwrap().children.len(), 6);
        for h in sh.tshape.read().unwrap().children.iter() {
            assert_eq!(h.shape_type(), ShapeType::Face);
        }
    }

    #[test]
    fn copy_edge_preserves_geometry() {
        let b = unit_box();
        let copy = AlgoTools::copy_edge(&b.edges[0]).expect("copy edge");
        assert!(!copy.0.same_tshape(&b.edges[0].0));
        assert_eq!(BRepTool::edge_parameters(&copy), BRepTool::edge_parameters(&b.edges[0]));
        let (a1, b1) = BRepTool::edge_vertices(&b.edges[0]).expect("orig endpoints");
        let (a2, b2) = BRepTool::edge_vertices(&copy).expect("copy endpoints");
        assert!(a1.is_equal(&a2) && b1.is_equal(&b2));
        assert_eq!(copy.0.tshape.read().unwrap().children.len(), 2);
    }

    #[test]
    fn make_split_edge_bounds_at_params() {
        let b = unit_box();
        let p1 = AlgoTools::point_on_edge(&b.edges[0], 0.25).expect("p1");
        let p2 = AlgoTools::point_on_edge(&b.edges[0], 0.75).expect("p2");
        let v1 = AlgoTools::make_new_vertex(&p1, 1e-7).expect("v1");
        let v2 = AlgoTools::make_new_vertex(&p2, 1e-7).expect("v2");
        let e = AlgoTools::make_split_edge(&b.edges[0], Some(&v1), 0.25, Some(&v2), 0.75)
            .expect("split edge");
        assert_eq!(BRepTool::edge_parameters(&e), (0.25, 0.75));
        let q1 = AlgoTools::point_on_edge(&e, 0.25).expect("start");
        let q2 = AlgoTools::point_on_edge(&e, 0.75).expect("end");
        assert!(q1.distance(&GpPnt::new(0.25, 0.0, 0.0)) < 1e-9, "start {q1:?}");
        assert!(q2.distance(&GpPnt::new(0.75, 0.0, 0.0)) < 1e-9, "end {q2:?}");
    }

    #[test]
    fn is_micro_edge_short_true_long_false() {
        let b = unit_box();
        assert!(!AlgoTools::is_micro_edge(&b.edges[0], 0.5));
        let tb = TopoBuilder::new();
        let micro = tb.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1e-6, 0.0, 0.0));
        assert!(AlgoTools::is_micro_edge(&micro, 0.5));
    }

    #[test]
    fn is_open_shell_closed_false_open_true() {
        let b = unit_box();
        let shell = b.solid.0.tshape.read().unwrap().children[0].clone();
        assert!(!AlgoTools::is_open_shell(&shell));

        // Five faces of the box form an open shell (the right face is missing).
        let open = TopoBuilder::new().make_shell(&b.faces[0..5]);
        assert!(AlgoTools::is_open_shell(&open.0));
    }

    #[test]
    fn is_inverted_solid_normal_box_false() {
        let b = unit_box();
        assert!(!AlgoTools::is_inverted_solid(&b.solid.0));
    }

    #[test]
    fn sense_perpendicular_box_faces_is_zero() {
        let b = unit_box();
        // Bottom (normal −Z) and front (normal −Y) share edge 0; the normals
        // are perpendicular so the sense is 0.
        assert_eq!(AlgoTools::sense(&b.faces[0], &b.faces[2]), 0);
        // Two faces that do not share an edge also report 0.
        assert_eq!(AlgoTools::sense(&b.faces[0], &b.faces[1]), 0);
    }

    #[test]
    fn is_hole_detects_inner_wire() {
        let tb = TopoBuilder::new();
        let outer = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let hole = [
            GpPnt::new(0.75, 0.75, 0.0),
            GpPnt::new(0.75, 0.25, 0.0),
            GpPnt::new(0.25, 0.25, 0.0),
            GpPnt::new(0.25, 0.75, 0.0),
        ];
        let square = |pts: &[GpPnt; 4]| -> Vec<Edge> {
            (0..4).map(|i| tb.make_edge_segment(&pts[i], &pts[(i + 1) % 4])).collect()
        };
        let face = build_face_with_holes(&square(&outer), &[square(&hole)]).expect("face with hole");
        let wires = wires_of_face(&face);
        assert_eq!(wires.len(), 2);
        assert!(!AlgoTools::is_hole(&wires[0].0, &face), "outer wire is not a hole");
        assert!(AlgoTools::is_hole(&wires[1].0, &face), "inner wire is a hole");
    }

    #[test]
    fn dimension_by_shape_type() {
        let b = unit_box();
        assert_eq!(AlgoTools::dimension(&b.vertices[0].0), 0);
        assert_eq!(AlgoTools::dimension(&b.edges[0].0), 1);
        assert_eq!(AlgoTools::dimension(&b.faces[0].0), 2);
        assert_eq!(AlgoTools::dimension(&b.solid.0), 3);
        let wire = TopoBuilder::new().make_wire(&b.edges[0..4]);
        assert_eq!(AlgoTools::dimension(&wire.0), 1);
    }

    #[test]
    fn correct_tolerances_grows_vertex_tolerance() {
        let tb = TopoBuilder::new();
        let v = tb.make_vertex(GpPnt::new(0.0, 0.0, 0.5), 1e-7);
        let lin = GpLin::from_pnt_dir(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let mut e = tb.make_edge(Arc::new(GeomLine::new(lin)), 0.0, 1.0);
        tb.add(&mut e.0, &v.0);
        // The vertex sits 0.5 above the curve start; correct_tolerances must
        // grow its tolerance to cover the endpoint (capped by the 1.0 budget).
        AlgoTools::correct_tolerances(&e.0, 1.0);
        assert!(
            BRepTool::vertex_tolerance(&v) >= 0.5,
            "tol {}",
            BRepTool::vertex_tolerance(&v)
        );
    }

    #[test]
    fn get_edge_off_none_for_edge_on_face() {
        let b = unit_box();
        // Edge 0 lies on the bottom face → no off-face parameter.
        assert!(AlgoTools::get_edge_off(&b.edges[0], &b.faces[0]).is_none());
    }
}
