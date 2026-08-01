//! BVH-accelerated mesh queries and `Poly` mesh topology tools.
//!
//! This module layers the standard OCCT query surface on top of the triangle
//! BVH built by [`crate::bvh::builder_tri::build_tri_bvh`], and adds a set of
//! mesh utility functions that OCCT exposes through the `Poly` package and the
//! `BRepMesh` meshing algorithms:
//!
//! * Ray, segment and box queries against a triangle mesh (pruned by the BVH).
//! * Point containment and closest-point queries.
//! * Integral quantities: signed volume (divergence theorem) and surface area.
//! * Edge/vertex topology: boundary edges, manifold-closed test, connected
//!   components, boundary loops and per-vertex incident triangle lists.
//! * 2D polygon ear-clipping triangulation.
//!
//! Source: `BVH_Tools.hxx`, `BVH_BoxSet.hxx`, `Poly_Connect.hxx`,
//! `Poly_Triangulation.hxx`, `BRepMesh`-style acceleration.
//!
//! All mesh-level queries take a [`TriBvh`] built over the *same* triangle
//! slice they are called with — the leaf ranges stored in the BVH nodes index
//! directly into that slice. The module therefore speaks two triangle
//! representations and bridges them:
//!
//! * `&[(GpPnt, GpPnt, GpPnt)]` — a flat triangle soup; this is what the BVH
//!   builder consumes and what every query here accepts.
//! * `&[(usize, usize, usize)]` — an index mesh referencing a vertex array;
//!   this is the form used by the topology helpers ([`mesh_edge_topology`],
//!   [`mesh_connected_components`], [`mesh_boundary_loops`], ...).
//!
//! [`triangulation_triangles`] converts a [`crate::poly::Triangulation`] into
//! the flat soup form so meshes built by the `Poly` package can be queried
//! directly.
//!
//! ## Relationship to OCCT classes
//!
//! | This module | OCCT |
//! |---|---|
//! | [`ray_cast_mesh`], [`segment_query_mesh`] | `BVH_Tools::RaySegmentIntersection` |
//! | [`box_query_union`] | `BVH_BoxSet::Select` / `BVH_Tree` box query |
//! | [`point_inside_box`] | `BRepMesh` even-odd containment |
//! | [`closest_point_mesh`] | `Poly_Connect` + point-triangle distance |
//! | [`mesh_edge_topology`] | `Poly_Connect` edge counts |
//! | [`mesh_to_polygon_indices`] | `Poly_Connect::TriangleToNodes` |
//!
//! ## Performance notes
//!
//! The ray, segment and box queries prune with the BVH first and only test the
//! triangles that land in candidate leaves, so their cost is logarithmic in
//! the leaf count for meshes with a reasonable spatial distribution. The
//! topology helpers ([`mesh_edge_topology`], [`mesh_connected_components`],
//! [`mesh_boundary_loops`], [`mesh_is_manifold_closed`]) operate directly on
//! index triples and do **not** require a BVH at all — they are linear in the
//! number of edges. [`closest_point_mesh`] combines a BVH candidate lookup
//! with a full scan fallback (see its documentation), so its worst case is
//! linear in the triangle count; the BVH only short-circuits on-surface
//! queries.

use std::collections::HashMap;

use crate::bnd::BndBox;
use crate::bvh::builder_tri::{closest_triangle_to_point, query_triangles, TriBvh};
use crate::bvh::bvh_ops::{bvh_point_in_mesh, ray_hits_bbox, ray_triangle_t};
use crate::bvh::BvhNode;
use crate::gp::{GpPnt, GpPnt2d, GpVec};

// ---------------------------------------------------------------------------
// Poly bridge
// ---------------------------------------------------------------------------

/// Convert a [`crate::poly::Triangulation`] into a flat triangle soup.
///
/// `Poly_Triangulation` stores a node array and per-triangle index triples;
/// every BVH query in this module works on `(GpPnt, GpPnt, GpPnt)` slices, so
/// this is the one-line adapter that lets a `Poly` mesh be handed straight to
/// [`build_tri_bvh`](crate::bvh::builder_tri::build_tri_bvh) and then to any
/// of the ray / box / closest-point queries.
///
/// Triangles whose indices fall outside the node array are skipped, matching
/// the bounds-checking behaviour used throughout this crate's mesh utilities.
///
/// * `tri` — a `Poly` triangulation (nodes + index triples).
///
/// Returns the point-triangle soup in the same triangle order.
pub fn triangulation_triangles(
    tri: &crate::poly::Triangulation,
) -> Vec<(GpPnt, GpPnt, GpPnt)> {
    let nodes = &tri.nodes;
    tri.triangles
        .iter()
        .filter_map(|t| {
            let (a, b, c) = (t.n0, t.n1, t.n2);
            if a < nodes.len() && b < nodes.len() && c < nodes.len() {
                Some((nodes[a], nodes[b], nodes[c]))
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Ray / segment queries
// ---------------------------------------------------------------------------

/// The result of a ray or segment intersection against a triangle mesh.
///
/// Mirrors the payload OCCT's `BVH_Tools::RaySegmentIntersection` helpers
/// produce: a boolean hit flag, the intersection distance `t` along the ray
/// (or segment) direction, the interpolated hit point, and the index of the
/// triangle that was struck.
///
/// For a ray query the parameterisation is `point = origin + t * dir` with
/// `dir` used exactly as passed to [`ray_cast_mesh`] (not necessarily unit
/// length). For a segment query `t` lies in `[0, length]` after the direction
/// is normalised internally by [`segment_query_mesh`].
#[derive(Debug, Clone, Copy)]
pub struct RayQueryResult {
    /// Whether an intersection was found. Always `true` when the value is
    /// produced by [`ray_cast_mesh`] / [`segment_query_mesh`]; kept as a field
    /// for structural compatibility with the OCCT query result.
    pub hit: bool,
    /// The ray parameter at the closest intersection.
    pub t: f64,
    /// The world-space hit point `origin + t * dir`.
    pub point: GpPnt,
    /// Index into the triangle slice of the hit triangle.
    pub triangle_index: usize,
}

/// Collect the leaf candidate index ranges a ray passes through.
///
/// Walks the BVH from the root, pruning any subtree whose node bounding box
/// does not intersect the ray (slab test via [`ray_hits_bbox`]). Each leaf is
/// recorded as a `(start_idx, end_idx)` range into the original triangle
/// slice — the same convention used by [`crate::bvh::builder_tri::TriBvh`]
/// leaves.
fn ray_leaves(node: &BvhNode, origin: &GpPnt, dir: &GpVec, out: &mut Vec<(usize, usize)>) {
    if !ray_hits_bbox(origin, dir, &node.bbox) {
        return;
    }
    if node.is_leaf() {
        out.push((node.start_idx, node.end_idx));
        return;
    }
    if let Some(l) = &node.left {
        ray_leaves(l, origin, dir, out);
    }
    if let Some(r) = &node.right {
        ray_leaves(r, origin, dir, out);
    }
}

/// Cast a ray against a triangle mesh and return the *nearest* hit.
///
/// The BVH is marched by collecting the leaves whose bounding box the ray
/// crosses ([`ray_leaves`]); every candidate triangle in those leaves is then
/// tested with the Möller–Trumbore primitive [`ray_triangle_t`]. Among the
/// hits with `t <= max_t` the closest one is kept, so this returns the first
/// surface point the ray reaches *before* the caller's distance budget, not
/// merely the first triangle the traversal happens to visit.
///
/// The direction `dir` is used unnormalised: `t` is measured in the same units
/// as `dir`, and the returned [`RayQueryResult::point`] is
/// `origin + t * dir`.
///
/// * `bvh` — BVH built over `triangles`.
/// * `triangles` — the triangle soup being queried.
/// * `origin` — ray start point.
/// * `dir` — ray direction (need not be unit length; a zero vector yields
///   `None`).
/// * `max_t` — distance budget; only hits with `0 < t <= max_t` are reported.
///
/// Returns `None` when the ray misses the mesh or when every hit lies beyond
/// `max_t`.
///
/// # Example
/// ```ignore
/// let bvh = build_tri_bvh(&tris, 4);
/// if let Some(r) = ray_cast_mesh(&bvh, &tris, origin, dir, 1e3) {
///     println!("hit triangle {} at t = {}", r.triangle_index, r.t);
/// }
/// ```
pub fn ray_cast_mesh(
    bvh: &TriBvh,
    triangles: &[(GpPnt, GpPnt, GpPnt)],
    origin: GpPnt,
    dir: GpVec,
    max_t: f64,
) -> Option<RayQueryResult> {
    if dir.xyz().square_modulus() <= 1e-30 {
        return None;
    }
    let root = bvh.root.as_ref()?;
    let mut ranges = Vec::new();
    ray_leaves(root, &origin, &dir, &mut ranges);
    let mut best: Option<(usize, f64)> = None;
    for (s, e) in ranges {
        for i in s..e {
            if i >= triangles.len() {
                continue;
            }
            let (a, b, c) = triangles[i];
            if let Some(t) = ray_triangle_t(&origin, &dir, &a, &b, &c) {
                if t <= max_t && best.map_or(true, |(_, bt)| t < bt) {
                    best = Some((i, t));
                }
            }
        }
    }
    best.map(|(ti, t)| RayQueryResult {
        hit: true,
        t,
        point: origin.translated_vec(&dir.multiplied_scalar(t)),
        triangle_index: ti,
    })
}

/// Query a finite segment `[a, b]` against a triangle mesh.
///
/// A segment is a clamped ray: the direction is normalised internally and the
/// distance budget is exactly the segment length. Any hit with
/// `0 < t <= |b - a|` is reported, and the returned [`RayQueryResult::t`] is
/// the distance from `a` along the segment. Degenerate segments (`a == b`)
/// yield `None`.
///
/// * `bvh` — BVH built over `triangles`.
/// * `triangles` — the triangle soup being queried.
/// * `a` — segment start.
/// * `b` — segment end.
///
/// Returns the nearest intersection, or `None` if the segment does not touch
/// the mesh.
pub fn segment_query_mesh(
    bvh: &TriBvh,
    triangles: &[(GpPnt, GpPnt, GpPnt)],
    a: GpPnt,
    b: GpPnt,
) -> Option<RayQueryResult> {
    let dir = GpVec::from_pnts(&a, &b);
    let len = dir.xyz().modulus();
    if len <= 1e-30 {
        return None;
    }
    let unit = dir.divided(len);
    ray_cast_mesh(bvh, triangles, a, unit, len)
}

// ---------------------------------------------------------------------------
// Bounding boxes
// ---------------------------------------------------------------------------

/// Compute the axis-aligned bounding box of a single triangle.
///
/// The box is grown from the three corner points; it is the same per-triangle
/// AABB that the BVH builder stores for every triangle, and it is the shape
/// that box queries and ray-pruning tests actually operate on.
///
/// * `a`, `b`, `c` — the triangle corners.
///
/// Returns a finite [`BndBox`] enclosing the triangle.
pub fn triangle_bounding_box(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> BndBox {
    let mut bb = BndBox::new();
    bb.add_point(a);
    bb.add_point(b);
    bb.add_point(c);
    bb
}

/// Compute the axis-aligned bounding box of a whole triangle mesh.
///
/// The box is the union of every triangle's bounding box (equivalently, the
/// union of all triangle corners). For an empty mesh a void box is returned.
///
/// * `triangles` — the triangle soup.
///
/// Returns the mesh AABB.
pub fn mesh_bounding_box(triangles: &[(GpPnt, GpPnt, GpPnt)]) -> BndBox {
    let mut bb = BndBox::new();
    for (a, b, c) in triangles {
        bb.add_point(a);
        bb.add_point(b);
        bb.add_point(c);
    }
    bb
}

// ---------------------------------------------------------------------------
// Box queries
// ---------------------------------------------------------------------------

/// Return the indices of every triangle whose bounding box overlaps a query
/// box.
///
/// This is a straight BVH box query via the builder's [`query_triangles`] (the
/// same primitive exposed as `crate::bvh::bvh_ops::bvh_box_query`) — the box
/// overlap test runs on the per-triangle bounding boxes stored in the BVH, so
/// the result is a conservative superset of the triangles whose *interior*
/// intersects the box. If an exact, triangle-by-triangle test is needed the
/// returned candidates must be filtered further.
///
/// * `bvh` — BVH built over the mesh.
/// * `box3d` — the axis-aligned query box.
///
/// Returns the triangle indices (in arbitrary traversal order, not sorted).
pub fn box_query_union(bvh: &TriBvh, box3d: &BndBox) -> Vec<usize> {
    let mut out = Vec::new();
    query_triangles(bvh, box3d, &mut out);
    out
}

// ---------------------------------------------------------------------------
// Point containment
// ---------------------------------------------------------------------------

/// Test whether a point lies inside a closed triangle mesh.
///
/// Delegates to [`bvh_point_in_mesh`], which casts a jittered ray in +X and
/// counts parity (even-odd rule) with coincident-hit de-duplication. The
/// jitter moves the ray origin off any coplanar geometry so that edge/grazing
/// hits do not corrupt the count.
///
/// * `bvh` — BVH built over `triangles`.
/// * `triangles` — the mesh (assumed closed and manifold).
/// * `p` — the query point.
///
/// Returns `true` for interior points, `false` for exterior points or empty
/// meshes.
pub fn point_inside_box(bvh: &TriBvh, triangles: &[(GpPnt, GpPnt, GpPnt)], p: GpPnt) -> bool {
    bvh_point_in_mesh(bvh, triangles, &p)
}

// ---------------------------------------------------------------------------
// Closest-point queries
// ---------------------------------------------------------------------------

/// Compute the closest point on a triangle to a query point.
///
/// Implements the region-based closest-point algorithm from Ericson,
/// *Real-Time Collision Detection* (Section 5.1.5): the barycentric
/// coordinates of the query point are classified into the seven Voronoi
/// regions of the triangle (three vertex regions, three edge regions and the
/// interior face region) and the projection onto the appropriate feature is
/// returned.
///
/// * `p` — the query point (anywhere in space).
/// * `a`, `b`, `c` — the triangle corners.
///
/// Returns `(closest_point, distance)` where `distance` is the Euclidean
/// distance from `p` to the triangle.
pub fn closest_point_on_triangle(
    p: &GpPnt,
    a: &GpPnt,
    b: &GpPnt,
    c: &GpPnt,
) -> (GpPnt, f64) {
    let ab = GpVec::from_pnts(a, b);
    let ac = GpVec::from_pnts(a, c);
    let ap = GpVec::from_pnts(a, p);

    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);

    // Vertex region A: p projects outside the corner at a.
    if d1 <= 0.0 && d2 <= 0.0 {
        return (*a, p.distance(a));
    }

    let bp = GpVec::from_pnts(b, p);
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);

    // Vertex region B: p projects outside the corner at b.
    if d3 >= 0.0 && d4 <= d3 {
        return (*b, p.distance(b));
    }

    // Edge region AB: project p onto the segment a-b.
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let pt = GpPnt::from_xyz(&a.coord.added(&ab.xyz().multiplied(v)));
        return (pt, p.distance(&pt));
    }

    let cp = GpVec::from_pnts(c, p);
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);

    // Vertex region C: p projects outside the corner at c.
    if d6 >= 0.0 && d5 <= d6 {
        return (*c, p.distance(c));
    }

    // Edge region AC: project p onto the segment a-c.
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let pt = GpPnt::from_xyz(&a.coord.added(&ac.xyz().multiplied(w)));
        return (pt, p.distance(&pt));
    }

    // Edge region BC: project p onto the segment b-c.
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let pt = GpPnt::from_xyz(&b.coord.added(&c.coord.subtracted(&b.coord).multiplied(w)));
        return (pt, p.distance(&pt));
    }

    // Face region: p projects inside the triangle; barycentric interpolation
    // of the projection is the closest point.
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let u = 1.0 - v - w;
    let pt = GpPnt::from_xyz(
        &a.coord.multiplied(u).added(&b.coord.multiplied(v)).added(&c.coord.multiplied(w)),
    );
    (pt, p.distance(&pt))
}

/// Find the closest point on a triangle mesh to a query point.
///
/// The BVH-accelerated candidate lookup from
/// [`closest_triangle_to_point`] is used as a fast path: it returns the
/// triangle whose *bounding box* is nearest to the point. The exact closest
/// point is then computed on that triangle with
/// [`closest_point_on_triangle`].
///
/// Because the BVH candidate probe uses a zero-size query box around `p`, it
/// only matches triangles whose bounding box actually contains the point; an
/// off-surface point would otherwise miss every candidate. A full scan over
/// the triangle slice is therefore performed as a correctness fallback, so
/// this function always returns the true closest point of the mesh.
///
/// * `bvh` — BVH built over `triangles`.
/// * `triangles` — the triangle soup.
/// * `p` — the query point.
///
/// Returns `(closest_point, distance)`; for an empty mesh the query point
/// itself is returned with an infinite distance.
pub fn closest_point_mesh(
    bvh: &TriBvh,
    triangles: &[(GpPnt, GpPnt, GpPnt)],
    p: GpPnt,
) -> (GpPnt, f64) {
    let mut best_pt = p;
    let mut best_d = f64::INFINITY;

    // Fast path: reuse the BVH triangle lookup.
    if let Some(i) = closest_triangle_to_point(bvh, &p, triangles) {
        let (a, b, c) = triangles[i];
        let (pt, d) = closest_point_on_triangle(&p, &a, &b, &c);
        if d < best_d {
            best_d = d;
            best_pt = pt;
        }
    }

    // ponytail: the BVH probe misses off-surface points, so a full scan is the
    // only guarantee of correctness. Replace with a growing-radius BVH query
    // if meshes get large.
    for &(a, b, c) in triangles {
        let (pt, d) = closest_point_on_triangle(&p, &a, &b, &c);
        if d < best_d {
            best_d = d;
            best_pt = pt;
        }
    }
    (best_pt, best_d)
}

// ---------------------------------------------------------------------------
// Integral quantities
// ---------------------------------------------------------------------------

/// Compute the enclosed volume of a closed triangle mesh.
///
/// Uses the divergence theorem: the signed volume is the sum over triangles
/// of the origin-based tetrahedron volume `v0 · (v1 × v2) / 6`. For a closed,
/// consistently oriented surface the origin terms cancel and the sum is the
/// enclosed volume, independent of the choice of origin. The result is
/// returned as an absolute value so that either winding direction yields a
/// positive volume.
///
/// The `bvh` argument is accepted for API symmetry with the other mesh
/// queries but is not needed by the computation.
///
/// * `bvh` — ignored (BVH is not required for the closed-surface integral).
/// * `triangles` — the closed triangle mesh.
///
/// Returns the enclosed volume. For non-closed meshes the value is the signed
/// "volume" of the open surface and has no strict geometric meaning.
pub fn mesh_volume(bvh: &TriBvh, triangles: &[(GpPnt, GpPnt, GpPnt)]) -> f64 {
    let _ = bvh;
    let mut vol = 0.0;
    for (a, b, c) in triangles {
        vol += a.coord.dot(&b.coord.crossed(&c.coord)) / 6.0;
    }
    vol.abs()
}

/// Compute the total surface area of a triangle mesh.
///
/// Each triangle contributes half the magnitude of the cross product of its
/// two edge vectors; the contributions are summed across the mesh. The
/// winding (orientation) of the triangles is irrelevant — only the areas are
/// accumulated.
///
/// * `triangles` — the triangle soup.
///
/// Returns the summed triangle areas.
pub fn mesh_surface_area(triangles: &[(GpPnt, GpPnt, GpPnt)]) -> f64 {
    let mut area = 0.0;
    for (a, b, c) in triangles {
        let ab = GpVec::from_pnts(a, b);
        let ac = GpVec::from_pnts(a, c);
        area += 0.5 * ab.xyz().crossed(ac.xyz()).modulus();
    }
    area
}

// ---------------------------------------------------------------------------
// Edge topology
// ---------------------------------------------------------------------------

/// Count how many times each undirected edge appears in a triangle mesh.
///
/// Every triangle contributes its three edges `(i, j), (j, k), (k, i)` with
/// the endpoints sorted so that the undirected edge `{a, b}` has a single
/// canonical key `(min, max)`. The returned map is the raw material for the
/// topology queries below.
fn edge_usage_counts(triangles: &[(usize, usize, usize)]) -> HashMap<(usize, usize), usize> {
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for &(i, j, k) in triangles {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            let key = if a < b { (a, b) } else { (b, a) };
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    counts
}

/// Map each undirected mesh edge to its incident triangles.
///
/// This is the dual adjacency table that OCCT's `Poly_Connect::EdgeTriangles`
/// builds: for every canonical edge key `(min, max)` it lists the triangle
/// indices that use that edge. A well-formed manifold edge has exactly two
/// entries; boundary edges have one; non-manifold edges have three or more.
/// It is the natural bridge between [`mesh_edge_topology`] and per-edge mesh
/// editing algorithms.
///
/// * `triangles` — the mesh as index triples.
///
/// Returns the edge → incident triangle map (only edges that appear in the
/// mesh are present).
pub fn edge_to_triangles(triangles: &[(usize, usize, usize)]) -> HashMap<(usize, usize), Vec<usize>> {
    let mut map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (ti, &(i, j, k)) in triangles.iter().enumerate() {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            let key = if a < b { (a, b) } else { (b, a) };
            map.entry(key).or_default().push(ti);
        }
    }
    map
}

/// Analyse the edge topology of a triangle mesh.
///
/// This is a thin re-export of the shared implementation in
/// [`crate::geom::mesh_analysis::mesh_edge_topology`] (kept here so the `Poly`
/// query surface is available from the BVH module as well):
///
/// * a *boundary* edge is used by exactly one triangle;
/// * a *non-manifold* edge is used by more than two triangles.
///
/// * `triangles` — the mesh as index triples.
///
/// Returns `(boundary_edge_count, non_manifold_edge_count)`.
pub fn mesh_edge_topology(triangles: &[(usize, usize, usize)]) -> (usize, usize) {
    crate::geom::mesh_analysis::mesh_edge_topology(triangles)
}

/// Whether a triangle mesh is a closed 2-manifold.
///
/// A mesh is a closed manifold exactly when every undirected edge is shared by
/// precisely two triangles (each edge has two incident faces and no dangling
/// boundary edges). This is the discrete analogue of "every point has a
/// neighbourhood homeomorphic to a disk" with the boundary empty.
///
/// * `triangles` — the mesh as index triples.
///
/// Returns `true` when every edge is used exactly twice. An empty mesh
/// vacuously satisfies the property.
pub fn mesh_is_manifold_closed(triangles: &[(usize, usize, usize)]) -> bool {
    edge_usage_counts(triangles).values().all(|&c| c == 2)
}

// ---------------------------------------------------------------------------
// Connected components
// ---------------------------------------------------------------------------

/// Find the root of the union-find set containing `x`, with path compression.
///
/// Iterative two-pass variant: the first pass walks to the root, the second
/// pass re-parents every node on the path directly under the root so that
/// subsequent lookups are amortised near-constant time (inverse Ackermann).
fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
    let mut root = x;
    while parent[root] != root {
        root = parent[root];
    }
    while parent[x] != x {
        let next = parent[x];
        parent[x] = root;
        x = next;
    }
    root
}

/// Merge the union-find sets containing `a` and `b`.
///
/// Union by attaching the root of `a` under the root of `b`. Sizes are not
/// tracked (union-by-rank is unnecessary for the small meshes this port
/// targets); the `find` path compression keeps the structure shallow.
fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let ra = uf_find(parent, a);
    let rb = uf_find(parent, b);
    if ra != rb {
        parent[ra] = rb;
    }
}

/// Partition the vertices of a triangle mesh into connected components.
///
/// Vertices are unioned together whenever they co-occur in a triangle, using
/// a union-find structure sized to the maximum vertex index. Only vertices
/// that actually appear in at least one triangle are reported (isolated
/// vertices are not emitted as singleton components).
///
/// * `triangles` — the mesh as index triples.
///
/// Returns a list of vertex-index sets, one per connected component. Each
/// component's vertex list is sorted ascending; components are ordered by
/// their smallest vertex index for deterministic output.
pub fn mesh_connected_components(triangles: &[(usize, usize, usize)]) -> Vec<Vec<usize>> {
    if triangles.is_empty() {
        return Vec::new();
    }
    let max_v = triangles.iter().flat_map(|&(a, b, c)| [a, b, c]).max().unwrap();
    let n = max_v + 1;
    let mut parent: Vec<usize> = (0..n).collect();
    let mut used = vec![false; n];
    for &(a, b, c) in triangles {
        used[a] = true;
        used[b] = true;
        used[c] = true;
        uf_union(&mut parent, a, b);
        uf_union(&mut parent, a, c);
    }
    let mut comps: HashMap<usize, Vec<usize>> = HashMap::new();
    for v in 0..n {
        if used[v] {
            let r = uf_find(&mut parent, v);
            comps.entry(r).or_default().push(v);
        }
    }
    let mut out: Vec<Vec<usize>> = comps.into_values().collect();
    for c in out.iter_mut() {
        c.sort_unstable();
    }
    out.sort_by(|a, b| a[0].cmp(&b[0]));
    out
}

// ---------------------------------------------------------------------------
// Boundary loops
// ---------------------------------------------------------------------------

/// Chain the boundary edges of a triangle mesh into closed loops.
///
/// The boundary of a manifold-with-boundary mesh consists of one or more
/// disjoint closed polylines. Each boundary edge (an edge used by exactly one
/// triangle) is followed to its far endpoint, where the walk continues along
/// the other incident boundary edge until the start vertex is revisited.
///
/// * `triangles` — the mesh as index triples.
///
/// Returns one `Vec<usize>` per boundary loop. Each loop is a closed cycle of
/// vertex indices, i.e. the first vertex is repeated as the last element (a
/// square's boundary is `[0, 1, 2, 3, 0]`, four edges). If a vertex is
/// incident to more than two boundary edges (a non-manifold boundary), the
/// chain follows the first unvisited edge at each step; the result may then
/// decompose a boundary into multiple open chains rather than closed loops.
pub fn mesh_boundary_loops(triangles: &[(usize, usize, usize)]) -> Vec<Vec<usize>> {
    let counts = edge_usage_counts(triangles);
    let mut edges: Vec<(usize, usize)> = Vec::new();
    let mut edge_id: HashMap<(usize, usize), usize> = HashMap::new();
    for (&e, &c) in &counts {
        if c == 1 {
            edge_id.insert(e, edges.len());
            edges.push(e);
        }
    }
    if edges.is_empty() {
        return Vec::new();
    }
    // Vertex -> incident boundary edge ids.
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for (id, &(u, v)) in edges.iter().enumerate() {
        adj.entry(u).or_default().push(id);
        adj.entry(v).or_default().push(id);
    }
    let mut visited = vec![false; edges.len()];
    let mut loops: Vec<Vec<usize>> = Vec::new();
    for start in 0..edges.len() {
        if visited[start] {
            continue;
        }
        let mut loop_verts: Vec<usize> = Vec::new();
        let mut cur_edge = start;
        let mut cur_vert = edges[start].0;
        visited[cur_edge] = true;
        loop_verts.push(cur_vert);
        loop {
            let (u, v) = edges[cur_edge];
            let next_vert = if u == cur_vert { v } else { u };
            loop_verts.push(next_vert);
            let next_edge = adj[&next_vert]
                .iter()
                .copied()
                .find(|&e| e != cur_edge && !visited[e]);
            match next_edge {
                Some(e) => {
                    visited[e] = true;
                    cur_edge = e;
                    cur_vert = next_vert;
                }
                None => break,
            }
            if cur_edge == start {
                break;
            }
        }
        loops.push(loop_verts);
    }
    loops
}

// ---------------------------------------------------------------------------
// 2D polygon triangulation
// ---------------------------------------------------------------------------

/// Triangulate a simple 2D polygon by ear clipping.
///
/// Wraps the shared ear-clipping implementation in
/// [`crate::geom::polygon_ops::triangulate_polygon2d`], which accepts a simple
/// polygon in either winding order and repeatedly clips convex "ears"
/// (corners that contain no other polygon vertex) until a single triangle
/// remains. The resulting triangle index triples reference positions in `pts`.
///
/// * `pts` — polygon vertices in boundary order (CCW or CW).
///
/// Returns `Some(triangles)` on success, or `None` for fewer than three
/// vertices, degenerate (zero-area) or self-intersecting input.
pub fn polygon_triangulate_ear(pts: &[GpPnt2d]) -> Option<Vec<(usize, usize, usize)>> {
    crate::geom::polygon_ops::triangulate_polygon2d(pts)
}

// ---------------------------------------------------------------------------
// Vertex-connectivity
// ---------------------------------------------------------------------------

/// Compute the per-vertex incident triangle lists of a mesh.
///
/// For each vertex index this returns the list of triangles that reference it,
/// i.e. the 1-ring of faces around the vertex. This is the data structure
/// `Poly_Connect::TriangleToNodes` / `Poly_Triangulation` adjacency tables
/// expose, and it is the starting point for vertex-normal averaging, boundary
/// extraction and local mesh editing.
///
/// * `verts` — the mesh vertex array; its length bounds the result.
/// * `triangles` — the mesh as index triples.
///
/// Returns one `Vec<usize>` per vertex; vertices that are not referenced by
/// any triangle get an empty list. Triangle indices are pushed in ascending
/// triangle order, so each list is naturally sorted.
pub fn mesh_to_polygon_indices(
    verts: &[GpPnt],
    triangles: &[(usize, usize, usize)],
) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); verts.len()];
    for (ti, &(a, b, c)) in triangles.iter().enumerate() {
        if a < verts.len() {
            out[a].push(ti);
        }
        if b < verts.len() {
            out[b].push(ti);
        }
        if c < verts.len() {
            out[c].push(ti);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit cube centred at the origin: 8 vertices, 12 triangles.
    ///
    /// The winding matches `geom::mesh_analysis::box_mesh` (outward CCW) so
    /// the divergence-theorem volume comes out positive.
    fn cube_triangles() -> Vec<(GpPnt, GpPnt, GpPnt)> {
        let v: Vec<GpPnt> = vec![
            GpPnt::new(-0.5, -0.5, -0.5), GpPnt::new(0.5, -0.5, -0.5),
            GpPnt::new(0.5, 0.5, -0.5), GpPnt::new(-0.5, 0.5, -0.5),
            GpPnt::new(-0.5, -0.5, 0.5), GpPnt::new(0.5, -0.5, 0.5),
            GpPnt::new(0.5, 0.5, 0.5), GpPnt::new(-0.5, 0.5, 0.5),
        ];
        let t: Vec<(usize, usize, usize)> = vec![
            (0, 3, 2), (0, 2, 1), // -Z
            (4, 5, 6), (4, 6, 7), // +Z
            (0, 1, 5), (0, 5, 4), // -Y
            (3, 7, 6), (3, 6, 2), // +Y
            (0, 4, 7), (0, 7, 3), // -X
            (1, 2, 6), (1, 6, 5), // +X
        ];
        t.into_iter()
            .map(|(i, j, k)| (v[i], v[j], v[k]))
            .collect()
    }

    /// Cube corner-index mesh (parallel to `cube_triangles`).
    fn cube_index_mesh() -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
        let v: Vec<GpPnt> = vec![
            GpPnt::new(-0.5, -0.5, -0.5), GpPnt::new(0.5, -0.5, -0.5),
            GpPnt::new(0.5, 0.5, -0.5), GpPnt::new(-0.5, 0.5, -0.5),
            GpPnt::new(-0.5, -0.5, 0.5), GpPnt::new(0.5, -0.5, 0.5),
            GpPnt::new(0.5, 0.5, 0.5), GpPnt::new(-0.5, 0.5, 0.5),
        ];
        let t: Vec<(usize, usize, usize)> = vec![
            (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
            (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
        ];
        (v, t)
    }

    #[test]
    fn ray_hits_triangle_nearest() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Ray from below the cube along +Z: front face at z = -0.5, so t = 4.5.
        // Offset in xy keeps the hit strictly inside a face triangle (off the
        // shared diagonal).
        let hit = ray_cast_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, -5.0), GpVec::new(0.0, 0.0, 1.0), 100.0)
            .expect("ray hits the cube");
        assert!((hit.t - 4.5).abs() < 1e-9, "t = {}", hit.t);
        assert!((hit.point.z() + 0.5).abs() < 1e-9, "z = {}", hit.point.z());
        assert!(hit.triangle_index < tris.len());
    }

    #[test]
    fn ray_misses_mesh() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Pointing away from the cube (downward, cube is above origin).
        assert!(ray_cast_mesh(&bvh, &tris, GpPnt::new(0.0, 0.0, -5.0), GpVec::new(0.0, 0.0, -1.0), 100.0).is_none());
        // Direction is zero.
        assert!(ray_cast_mesh(&bvh, &tris, GpPnt::new(0.0, 0.0, -5.0), GpVec::zero(), 100.0).is_none());
    }

    #[test]
    fn segment_clamped_hit() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Short segment stopping before the front face (z = -1) -> miss.
        let short = segment_query_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, -5.0), GpPnt::new(0.2, 0.1, -1.0));
        assert!(short.is_none(), "short segment must not reach the cube");
        // Long segment through the cube -> hit at distance 4.5 from a.
        let long = segment_query_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, -5.0), GpPnt::new(0.2, 0.1, 5.0))
            .expect("segment through the cube hits");
        assert!((long.t - 4.5).abs() < 1e-9, "t = {}", long.t);
    }

    #[test]
    fn box_query_overlap() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // A box fully inside the cube overlaps every face's triangle bbox.
        let mut q = BndBox::from_corners(&GpPnt::new(-0.2, -0.2, -0.2), &GpPnt::new(0.2, 0.2, 0.2));
        let hits = box_query_union(&bvh, &q);
        assert!(!hits.is_empty(), "interior box must overlap triangles");
        q = BndBox::from_corners(&GpPnt::new(10.0, 10.0, 10.0), &GpPnt::new(11.0, 11.0, 11.0));
        assert!(box_query_union(&bvh, &q).is_empty(), "distant box must miss");
    }

    #[test]
    fn point_inside_closed() {
        // Tetrahedron (0,0,0),(1,0,0),(0,1,0),(0,0,1) — 4 triangles, closed.
        let o = GpPnt::new(0.0, 0.0, 0.0);
        let x = GpPnt::new(1.0, 0.0, 0.0);
        let y = GpPnt::new(0.0, 1.0, 0.0);
        let z = GpPnt::new(0.0, 0.0, 1.0);
        let tris = vec![(o, x, y), (o, z, x), (o, y, z), (x, z, y)];
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        assert!(point_inside_box(&bvh, &tris, GpPnt::new(0.1, 0.1, 0.1)), "interior point");
        assert!(!point_inside_box(&bvh, &tris, GpPnt::new(10.0, 10.0, 10.0)), "exterior point");
    }

    #[test]
    fn closest_point_on_mesh() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Point above the +Z face (z = 0.5): distance should be 1.5.
        let (pt, d) = closest_point_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, 2.0));
        assert!((d - 1.5).abs() < 1e-9, "distance = {d}");
        assert!((pt.z() - 0.5).abs() < 1e-9, "closest z = {}", pt.z());
        // Point exactly on a face: distance ~ 0.
        let (_, d0) = closest_point_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, 0.5));
        assert!(d0 < 1e-9, "on-surface distance = {d0}");
    }

    #[test]
    fn mesh_volume_box() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        let vol = mesh_volume(&bvh, &tris);
        assert!((vol - 1.0).abs() < 0.05, "unit cube volume = {vol}");
    }

    #[test]
    fn surface_area_unit_square() {
        // Two triangles on a unit square in the z=0 plane.
        let tris = vec![
            (GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)),
            (GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)),
        ];
        assert!((mesh_surface_area(&tris) - 1.0).abs() < 1e-12, "area = {}", mesh_surface_area(&tris));
    }

    #[test]
    fn edge_topology_boundary_and_manifold() {
        // Single triangle: 3 boundary edges, 0 non-manifold.
        assert_eq!(mesh_edge_topology(&[(0, 1, 2)]), (3, 0));
        // Two triangles sharing edge (0,2): 4 boundary edges, 1 shared edge.
        assert_eq!(mesh_edge_topology(&[(0, 1, 2), (0, 2, 3)]), (4, 0));
        // Closed cube: 0 boundary.
        let (_, t) = cube_index_mesh();
        assert_eq!(mesh_edge_topology(&t), (0, 0));
    }

    #[test]
    fn manifold_closed_cube() {
        let (_, t) = cube_index_mesh();
        assert!(mesh_is_manifold_closed(&t), "cube is closed and manifold");
        // A single triangle is not closed (boundary edges).
        assert!(!mesh_is_manifold_closed(&[(0, 1, 2)]));
    }

    #[test]
    fn connected_components_two_boxes() {
        // Two disjoint cubes: vertices 0..7 and 8..15.
        let v: Vec<GpPnt> = (0..16)
            .map(|i| {
                if i < 8 {
                    GpPnt::new(i as f64, 0.0, 0.0)
                } else {
                    GpPnt::new(100.0 + i as f64, 0.0, 0.0)
                }
            })
            .collect();
        let cube: Vec<(usize, usize, usize)> = vec![
            (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
            (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
        ];
        let shifted: Vec<(usize, usize, usize)> =
            cube.iter().map(|&(a, b, c)| (a + 8, b + 8, c + 8)).collect();
        let mut tris = cube;
        tris.extend(shifted);
        let comps = mesh_connected_components(&tris);
        assert_eq!(comps.len(), 2, "two disjoint cubes -> two components");
        // Sorted by smallest vertex: first component is 0..8.
        assert_eq!(comps[0].len(), 8);
        assert_eq!(comps[1].len(), 8);
        assert!(comps.iter().all(|c| !c.is_empty()));
        let _ = v;
    }

    #[test]
    fn boundary_loops_square() {
        // Square: (0,0)-(1,0)-(1,1)-(0,1), two triangles.
        let tris = vec![(0, 1, 2), (0, 2, 3)];
        let loops = mesh_boundary_loops(&tris);
        assert_eq!(loops.len(), 1, "a square has one boundary loop");
        // Closed loop: first == last, four unique vertices -> 5 entries.
        let l = &loops[0];
        assert_eq!(l.len(), 5, "loop = {l:?}");
        assert_eq!(l[0], l[4], "loop is closed");
        // Closed cube: no boundary loops at all.
        let (_, t) = cube_index_mesh();
        assert!(mesh_boundary_loops(&t).is_empty());
    }

    #[test]
    fn ear_clipping_square() {
        let sq = [
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let tris = polygon_triangulate_ear(&sq).expect("square triangulates");
        assert_eq!(tris.len(), 2);
        let mut area = 0.0;
        for (a, b, c) in &tris {
            let (pa, pb, pc) = (sq[*a], sq[*b], sq[*c]);
            area += 0.5
                * ((pb.x() - pa.x()) * (pc.y() - pa.y()) - (pb.y() - pa.y()) * (pc.x() - pa.x()))
                    .abs();
        }
        assert!((area - 1.0).abs() < 1e-12, "triangulated area = {area}");
    }

    #[test]
    fn incident_triangle_lists() {
        let (v, t) = cube_index_mesh();
        let incident = mesh_to_polygon_indices(&v, &t);
        assert_eq!(incident.len(), 8);
        for (i, list) in incident.iter().enumerate() {
            assert!(!list.is_empty(), "vertex {i} has no incident triangles");
            // Every triangle in the list references vertex i.
            for &ti in list {
                let (a, b, c) = t[ti];
                assert!(a == i || b == i || c == i);
            }
        }
    }
}
