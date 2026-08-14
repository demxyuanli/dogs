//! Rolling-ball fillet along an edge chain — a simplified port of
//! `BRepFilletAPI_MakeFillet` / `ChFi3d` (constant-radius fillet).
//!
//! A sharp convex edge of a solid is replaced by a cylindrical blend surface
//! tangent to both adjacent planar faces at radius `R` (a quarter cylinder for
//! a 90° corner). The two adjacent faces are trimmed to their tangency lines,
//! and the two end faces (perpendicular to the edge) are rebuilt with a
//! quarter-circular arc replacing the corner vertex. Edge instances are shared
//! between faces so the rebuilt shell satisfies `shell_is_closed`.
//!
//! Non-planar adjacent faces and reflex (concave) corners are rejected with an
//! error. The corner patch (`fillet_corner_solid`) is implemented as a
//! spherical octant — see the note there.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx3, GpCylinder, GpDir, GpPnt, GpSphere, GpVec};
use occt_geom::{GeomCylinder, GeomSphere, Surface};

use crate::brep_surface::{classify_surface, face_plane, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, is_same, vertex_position, wires_of_face,
};

/// Fillet parameters for a rolling-ball blend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FilletSpec {
    pub radius: f64,
}

impl FilletSpec {
    /// Validate the radius is strictly positive.
    pub fn check(&self) -> Result<(), String> {
        if self.radius <= 0.0 || !self.radius.is_finite() {
            return Err("fillet: radius must be a positive finite value".to_string());
        }
        Ok(())
    }
}

/// Builds and caches topology edges so faces that must share a boundary edge
/// all receive the *same* `Edge` instance (required for `shell_is_closed`,
/// which counts edge usage by `TShape` identity).
///
/// Points are deduplicated by position; line segments are keyed by their
/// (unordered) endpoint pair; circular arcs by the same. Original edges that
/// are fully retained by a face are registered here too, so a rebuilt face
/// reuses the identical edge instance that a kept neighbour face still
/// references.
struct EdgeCache {
    b: TopoBuilder,
    pts: Vec<GpPnt>,
    segs: HashMap<(usize, usize), Edge>,
    arcs: HashMap<(usize, usize), Edge>,
}

impl EdgeCache {
    fn new() -> Self {
        EdgeCache {
            b: TopoBuilder::new(),
            pts: Vec::new(),
            segs: HashMap::new(),
            arcs: HashMap::new(),
        }
    }

    /// Index of an existing point within `tol` of `p`, or a fresh index.
    ///
    /// A 1e-6 tolerance: tangency points are recomputed independently from
    /// clipped polygon intersections, which can differ by ~1e-9 (bisection
    /// tolerance); the merged points must resolve to the same cache key so the
    /// shared edge instance is reused. Distinct feature points in the shapes
    /// handled here are orders of magnitude farther apart.
    fn point_index(&mut self, p: GpPnt) -> usize {
        for (i, q) in self.pts.iter().enumerate() {
            if p.distance(q) < 1e-6 {
                return i;
            }
        }
        self.pts.push(p);
        self.pts.len() - 1
    }

    /// Register an original edge so future `seg` requests between the same two
    /// points return this exact instance.
    fn register_original(&mut self, a: &GpPnt, b: &GpPnt, e: &Edge) {
        let ia = self.point_index(*a);
        let ib = self.point_index(*b);
        self.segs.insert((ia.min(ib), ia.max(ib)), e.clone());
    }

    /// A straight segment edge from `a` to `b` (shared by all requesters).
    fn seg(&mut self, a: &GpPnt, b: &GpPnt) -> Edge {
        if a.distance(b) < 1e-9 {
            panic!("fillet: degenerate segment at {:?}", a);
        }
        let ia = self.point_index(*a);
        let ib = self.point_index(*b);
        let key = (ia.min(ib), ia.max(ib));
        if let Some(e) = self.segs.get(&key) {
            return e.clone();
        }
        let e = self.b.make_edge_segment(&self.pts[ia], &self.pts[ib]);
        self.segs.insert(key, e.clone());
        e
    }

    /// A circular-arc edge from `a` to `b` on the circle of `radius` centred at
    /// `center` with plane normal `normal`, shared by all requesters. The arc
    /// is the minor arc that sweeps from `a` to `b` counter-clockwise in the
    /// frame whose X axis points at `a`.
    fn arc(
        &mut self,
        center: &GpPnt,
        normal: &GpDir,
        radius: f64,
        a: &GpPnt,
        b: &GpPnt,
    ) -> Result<Edge, String> {
        let ia = self.point_index(*a);
        let ib = self.point_index(*b);
        let key = (ia.min(ib), ia.max(ib));
        if let Some(e) = self.arcs.get(&key) {
            return Ok(e.clone());
        }
        let xd = GpDir::from_vec(&GpVec::from_pnts(center, a))
            .map_err(|_| "fillet: arc start coincides with the arc centre".to_string())?;
        let ax2 = GpAx2::new(*center, *normal, xd)
            .map_err(|e| format!("fillet: arc frame is degenerate: {e}"))?;
        let xdir = *ax2.x_direction();
        let ydir = *ax2.y_direction();
        let va = GpVec::from_pnts(center, a);
        let vb = GpVec::from_pnts(center, b);
        let a1 = va.coord.dot(&ydir.xyz()).atan2(va.coord.dot(&xdir.xyz()));
        let a2 = vb.coord.dot(&ydir.xyz()).atan2(vb.coord.dot(&xdir.xyz()));
        let mut e = self.b.make_edge_circle(&ax2, radius, a1, a2);
        let v1 = self.b.make_vertex(*a, 0.0);
        let v2 = self.b.make_vertex(*b, 0.0);
        self.b.add(&mut e.0, &v1.0);
        self.b.add(&mut e.0, &v2.0);
        self.arcs.insert(key, e.clone());
        Ok(e)
    }
}

/// Compute the two in-wedge unit directions and the wedge angle for two planar
/// faces meeting at a straight edge.
///
/// `a` is the unit edge direction, `n1`/`n2` the two outward face normals, and
/// `p` any point on the edge. The returned `u1`/`u2` lie in the plane
/// perpendicular to the edge, point *into* the material corner, and are unit
/// vectors along each face. The wedge angle `theta ∈ (0, π)` is the interior
/// (material) angle between the faces.
fn in_wedge_directions(a: &GpVec, n1: &GpVec, n2: &GpVec) -> Result<(GpVec, GpVec, f64), String> {
    // b_i = a × n_i is a unit direction along face i in the cross-section.
    let b1 = a.crossed(n1).normalized();
    let b2 = a.crossed(n2).normalized();
    let mag = b1.xyz().crossed(&b2.xyz()).modulus();
    if mag < 1e-12 {
        return Err("fillet: adjacent faces are coplanar (no wedge)".to_string());
    }
    // Choose the sign of each so it points into the other face's material
    // half-space (dot with the other outward normal is negative).
    let s1 = if b1.dot(n2) < 0.0 { 1.0 } else { -1.0 };
    let s2 = if b2.dot(n1) < 0.0 { 1.0 } else { -1.0 };
    let u1 = b1.multiplied_scalar(s1).normalized();
    let u2 = b2.multiplied_scalar(s2).normalized();
    let cos_t = u1.dot(&u2).clamp(-1.0, 1.0);
    let sin_t = u1.xyz().crossed(&u2.xyz()).modulus();
    if sin_t < 1e-9 {
        return Err("fillet: degenerate wedge angle".to_string());
    }
    // For a convex corner the wedge angle is in (0, π).
    if cos_t <= 0.0 && sin_t < 1e-6 {
        return Err("fillet: reflex (concave) corner unsupported".to_string());
    }
    let theta = cos_t.acos();
    Ok((u1, u2, theta))
}

/// Geometric description of a constant-radius fillet along one straight edge.
#[derive(Debug, Clone)]
struct BlendGeom {
    /// Tangency point on face 1 at the two edge endpoints.
    t1_0: GpPnt,
    t1_1: GpPnt,
    /// Tangency point on face 2 at the two edge endpoints.
    t2_0: GpPnt,
    t2_1: GpPnt,
    /// Cylinder axis (a point on the axis at the p0 cross-section, direction).
    axis_origin: GpPnt,
    axis_dir: GpVec,
    /// Cylinder frame: X points at t1_0 (so u = 0 → t1_0), Y chosen so the
    /// angular sweep to t2_0 is positive.
    frame: GpAx3,
    /// Face 1 in-wedge direction (into the material corner).
    u1_dir: GpVec,
    /// Face 2 in-wedge direction.
    u2_dir: GpVec,
    radius: f64,
}

/// Build the blend geometry for a fillet of `radius` along the straight edge
/// `p0 → p1`, with the two adjacent planar faces' outward normals `n1`/`n2`.
fn blend_geometry(
    p0: &GpPnt,
    p1: &GpPnt,
    n1: &GpVec,
    n2: &GpVec,
    radius: f64,
) -> Result<BlendGeom, String> {
    let a = GpVec::from_pnts(p0, p1).normalized();
    let a_mag = GpVec::from_pnts(p0, p1).magnitude();
    if a_mag < 1e-12 {
        return Err("fillet: edge has zero length".to_string());
    }
    let (u1, u2, theta) = in_wedge_directions(&a, n1, n2)?;
    let sin_t = theta.sin();
    if sin_t < 1e-9 {
        return Err("fillet: degenerate wedge".to_string());
    }
    let cot_half = (1.0 + theta.cos()) / sin_t;
    let off = radius * cot_half; // tangency distance from the edge along each face
    let t1_0 = p0.translated_vec(&u1.multiplied_scalar(off));
    let t1_1 = p1.translated_vec(&u1.multiplied_scalar(off));
    let t2_0 = p0.translated_vec(&u2.multiplied_scalar(off));
    let t2_1 = p1.translated_vec(&u2.multiplied_scalar(off));
    let center_off = u1.added(&u2).multiplied_scalar(radius / sin_t);
    let c0 = p0.translated_vec(&center_off);
    let c1 = p1.translated_vec(&center_off);

    // Cylinder frame: X points at t1_0 so u = 0 lands on the face-1 tangency.
    let xd_vec = GpVec::from_pnts(&c0, &t1_0).normalized();
    let xd = GpDir::from_xyz(&xd_vec.xyz()).map_err(|_| "fillet: degenerate blend frame")?;
    let zd = GpDir::from_xyz(&a.xyz()).map_err(|_| "fillet: degenerate edge axis")?;
    let ax3 = GpAx3::new(c0, zd, &xd).map_err(|e| format!("fillet: frame error: {e}"))?;
    let xdir = *ax3.x_direction();
    let ydir = *ax3.y_direction();
    // Verify the angular sweep from t1_0 to t2_0 is sane (quarter turn for a
    // 90° corner); the arc edge builders re-derive the exact angles from the
    // frame, so this is only a sanity check on the geometry.
    let r2 = GpVec::from_pnts(&c0, &t2_0);
    let u2_angle = r2.coord.dot(&ydir.xyz()).atan2(r2.coord.dot(&xdir.xyz()));
    if u2_angle.abs() < 1e-6 || u2_angle.abs() > PI - 1e-6 {
        return Err("fillet: blend arc angle is degenerate".to_string());
    }

    let _ = (c1, a_mag);
    Ok(BlendGeom {
        t1_0,
        t1_1,
        t2_0,
        t2_1,
        axis_origin: c0,
        axis_dir: a,
        frame: ax3,
        u1_dir: u1,
        u2_dir: u2,
        radius,
    })
}

/// Faces of `shape` whose boundary wire contains the edge `edge` (same `TShape`).
pub fn faces_touching_edge(shape: &TopoShape, edge: &Edge) -> Vec<Face> {
    faces_of(shape)
        .into_iter()
        .filter(|f| {
            wires_of_face(f).iter().any(|w| {
                edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0))
            })
        })
        .collect()
}

/// Whether `p` coincides with any vertex of `face`'s outer boundary wire.
fn face_contains_point(face: &Face, p: &GpPnt) -> bool {
    wires_of_face(face).iter().any(|w| {
        edges_of_wire(w).iter().any(|e| {
            let (a, b) = edge_vertices(e);
            a.map_or(false, |v| vertex_position(&v).distance(p) < 1e-9)
                || b.map_or(false, |v| vertex_position(&v).distance(p) < 1e-9)
        })
    })
}

/// The ordered polygon of a face's outer wire: 3D vertex points around the
/// boundary, in a connected cycle. The wire's edges are chained by their
/// endpoint coordinates (shared edges may be stored in any order and any
/// orientation), so the result is the true boundary cycle regardless of how
/// the builder stored the wire. Returns an error for faces with no wire.
fn face_polygon(face: &Face) -> Result<Vec<GpPnt>, String> {
    let w = wires_of_face(face)
        .into_iter()
        .next()
        .ok_or("fillet: face has no boundary wire")?;
    let es = edges_of_wire(&w);
    if es.is_empty() {
        return Err("fillet: face wire has no edges".to_string());
    }
    let ends: Vec<(GpPnt, GpPnt)> = es
        .iter()
        .map(|e| {
            let (a, b) = edge_vertices(e);
            let a = a.ok_or("fillet: edge has no start vertex")?;
            let b = b.ok_or("fillet: edge has no end vertex")?;
            Ok((vertex_position(&a), vertex_position(&b)))
        })
        .collect::<Result<_, String>>()?;
    let mut used = vec![false; ends.len()];
    let mut poly: Vec<GpPnt> = Vec::new();
    let (a0, b0) = ends[0];
    poly.push(a0);
    poly.push(b0);
    used[0] = true;
    let mut cur = b0;
    let mut guard = 0;
    while guard < ends.len() {
        let mut found = false;
        for (i, (a, b)) in ends.iter().enumerate() {
            if used[i] {
                continue;
            }
            if a.distance(&cur) < 1e-9 {
                poly.push(*b);
                cur = *b;
                used[i] = true;
                found = true;
                break;
            } else if b.distance(&cur) < 1e-9 {
                poly.push(*a);
                cur = *a;
                used[i] = true;
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
        guard += 1;
    }
    if poly.len() > 1 && poly[0].distance(poly.last().unwrap()) < 1e-9 {
        poly.pop();
    }
    Ok(poly)
}

/// Sutherland–Hodgman clip of `poly` against the half-plane `inside(p)`.
/// Returns the clipped polygon (counter-clockwise or clockwise, matching input
/// winding). The polygon is treated as closed (last vertex wraps to first).
fn clip_polygon(poly: &[GpPnt], inside: &impl Fn(&GpPnt) -> bool) -> Vec<GpPnt> {
    if poly.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<GpPnt> = Vec::new();
    let n = poly.len();
    let mut prev = poly[n - 1];
    let mut prev_in = inside(&prev);
    for i in 0..n {
        let cur = poly[i];
        let cur_in = inside(&cur);
        if cur_in != prev_in {
            // Intersection of segment prev→cur with the clipping boundary.
            let d = GpVec::from_pnts(&prev, &cur);
            let t = intersect_t(&prev, &cur, inside);
            if let Some(tv) = t {
                let p = prev.translated_vec(&d.multiplied_scalar(tv));
                out.push(p);
            }
        }
        if cur_in {
            out.push(cur);
        }
        prev = cur;
        prev_in = cur_in;
    }
    dedupe_polygon(&out)
}

/// Remove consecutive (and wrap-around) duplicate vertices, keeping the first
/// of each run. The polygon must have at least 3 distinct vertices afterwards.
fn dedupe_polygon(poly: &[GpPnt]) -> Vec<GpPnt> {
    if poly.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<GpPnt> = Vec::new();
    for (i, p) in poly.iter().enumerate() {
        let next = &poly[(i + 1) % poly.len()];
        if p.distance(next) > 1e-9 {
            out.push(*p);
        }
    }
    // Also drop a trailing point equal to the first.
    while out.len() >= 2 && out[0].distance(out.last().unwrap()) < 1e-9 {
        out.pop();
    }
    out
}

/// Bisection search for the fraction `t ∈ [0,1]` where the segment `a→b`
/// crosses the `inside` boundary (used by the polygon clipper).
fn intersect_t(a: &GpPnt, b: &GpPnt, inside: &impl Fn(&GpPnt) -> bool) -> Option<f64> {
    let mut lo = 0.0f64;
    let mut hi = 1.0f64;
    let mut lo_in = inside(a);
    if lo_in == inside(b) {
        return None;
    }
    for _ in 0..48 {
        let mid = 0.5 * (lo + hi);
        let pm = a.translated_vec(&GpVec::from_pnts(a, b).multiplied_scalar(mid));
        let m_in = inside(&pm);
        if lo_in == m_in {
            lo = mid;
            lo_in = m_in;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// Build the wire for a trimmed planar face whose retained polygon is `poly`.
/// Every polygon segment is created through the shared `cache`, so edges that
/// coincide with a registered original edge (or another face's request) are
/// reused.
fn build_polygon_wire(cache: &mut EdgeCache, poly: &[GpPnt]) -> Result<Wire, String> {
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    let n = poly.len();
    if n < 3 {
        return Err("fillet: clipped face polygon has fewer than 3 vertices".to_string());
    }
    for i in 0..n {
        let a = &poly[i];
        let c = &poly[(i + 1) % n];
        edges.push(cache.seg(a, c));
    }
    Ok(b.make_wire(&edges))
}

/// Register the original edges of `face` with the cache so future segment
/// requests reuse them. Returns the face's polygon (in wire order) alongside.
fn register_face_edges(face: &Face, cache: &mut EdgeCache) -> Result<Vec<GpPnt>, String> {
    let poly = face_polygon(face)?;
    let w = wires_of_face(face)
        .into_iter()
        .next()
        .ok_or("fillet: face has no wire")?;
    let es = edges_of_wire(&w);
    // Register each original edge by its *own* endpoint pair (the polygon
    // segments chain through these edges, but the wire order need not match
    // the polygon order, and a shared edge may be stored with either
    // orientation).
    for e in &es {
        let (a, b) = edge_vertices(e);
        let pa = vertex_position(&a.ok_or("fillet: edge has no start vertex")?);
        let pb = vertex_position(&b.ok_or("fillet: edge has no end vertex")?);
        cache.register_original(&pa, &pb, e);
    }
    Ok(poly)
}

/// The plane normal (outward) of a planar face.
fn face_outward_normal(face: &Face) -> Result<GpVec, String> {
    let pln = face_plane(face).ok_or("fillet: face is not planar")?;
    let n = *pln.axis().direction().xyz();
    Ok(GpVec::from_xyz(&n))
}

/// Rebuild one planar face into its trimmed version: the retained region is
/// `inside(p)` (a half-plane for adjacent faces, the outside-of-sphere test for
/// corner faces). Original edges that survive untouched are reused so shared
/// boundaries with kept neighbour faces stay intact.
fn rebuild_face(face: &Face, inside: &impl Fn(&GpPnt) -> bool, cache: &mut EdgeCache) -> Result<Face, String> {
    let poly = register_face_edges(face, cache)?;
    let clipped = clip_polygon(&poly, inside);
    if clipped.len() < 3 {
        return Err("fillet: clipping a face left no retained region".to_string());
    }
    let wire = build_polygon_wire(cache, &clipped)?;
    let surf = BRepTool::face_surface(face).ok_or("fillet: face has no surface")?;
    Ok(cache.b.make_face(surf, &[wire]))
}

/// Rebuild the end face at an edge endpoint: replace the corner vertex `p`
/// (on the filleted edge) with the tangency points `t1`, `t2` and the blend
/// arc between them, all shared with the neighbouring rebuilt faces.
///
/// `center` is the blend-cylinder axis point at this endpoint (the arc's
/// centre), `n1`/`n2` are the adjacent faces' outward normals used to tell
/// which incident edge of the end face lies on which adjacent face.
fn rebuild_end_face(
    face: &Face,
    p: &GpPnt,
    t1: &GpPnt,
    t2: &GpPnt,
    center: &GpPnt,
    axis_dir: &GpVec,
    radius: f64,
    n1: &GpVec,
    cache: &mut EdgeCache,
) -> Result<Face, String> {
    let poly = register_face_edges(face, cache)?;
    let n = poly.len();
    let idx = poly
        .iter()
        .position(|q| q.distance(p) < 1e-9)
        .ok_or("fillet: end face does not contain the edge endpoint")?;
    let prev = &poly[(idx + n - 1) % n];
    let next = &poly[(idx + 1) % n];

    // The shared edge prev→p is on face 1 iff its direction is perpendicular
    // to face 1's outward normal (the edge lies in face 1's plane).
    let dir_prev = GpVec::from_pnts(prev, p).normalized();
    let dir_next = GpVec::from_pnts(p, next).normalized();
    let prev_is_face1 = dir_prev.dot(n1).abs() < 1e-6;
    let next_is_face1 = dir_next.dot(n1).abs() < 1e-6;

    let (first_seg_a, first_seg_b, last_seg_a, last_seg_b);
    if prev_is_face1 && !next_is_face1 {
        // prev→p is the face-1 edge (ends at t1); p→next is the face-2 edge.
        first_seg_a = prev;
        first_seg_b = t1;
        last_seg_a = t2;
        last_seg_b = next;
    } else if next_is_face1 && !prev_is_face1 {
        // p→next is the face-1 edge (ends at t1); prev→p is the face-2 edge
        // (ends at t2). The wire runs prev→t2, arc, t1→next.
        first_seg_a = prev;
        first_seg_b = t2;
        last_seg_a = t1;
        last_seg_b = next;
    } else {
        return Err("fillet: cannot classify end-face incident edges".to_string());
    }

    // The blend arc at this endpoint, shared with the blend cylinder face.
    let zd = GpDir::from_xyz(&axis_dir.xyz()).map_err(|e| e.to_string())?;
    let arc_edge = cache.arc(center, &zd, radius, t1, t2)?;

    // Build the new wire: [face-1 edge → t1, arc t1→t2, t2 → face-2 edge] then
    // the retained original edges continuing around the polygon.
    let b = TopoBuilder::new();
    let mut wire_edges: Vec<Edge> = Vec::new();
    wire_edges.push(cache.seg(first_seg_a, first_seg_b));
    wire_edges.push(arc_edge.clone());
    wire_edges.push(cache.seg(last_seg_a, last_seg_b));
    // Continue from `next` around the polygon back to `prev`, reusing original
    // edges.
    let mut j = (idx + 1) % n;
    let mut guard = 0;
    while j != idx && guard < n {
        if j == (idx + n - 1) % n {
            break; // the last segment returns to prev; already handled above.
        }
        let a = &poly[j];
        let c = &poly[(j + 1) % n];
        wire_edges.push(cache.seg(a, c));
        j = (j + 1) % n;
        guard += 1;
    }
    let wire = b.make_wire(&wire_edges);
    let surf = BRepTool::face_surface(face).ok_or("fillet: face has no surface")?;
    Ok(b.make_face(surf, &[wire]))
}

/// Build the blend cylinder face: a `GeomCylinder` surface with a wire of four
/// shared edges (two tangency lines, two quarter arcs).
fn build_blend_face(blend: &BlendGeom, cache: &mut EdgeCache) -> Result<Face, String> {
    let ax3 = blend.frame.clone();
    let gcyl = GpCylinder::new(ax3, blend.radius).map_err(|e| e.to_string())?;
    let surface: Arc<dyn Surface> = Arc::new(GeomCylinder::new(gcyl));
    let b = TopoBuilder::new();

    let t1_0 = &blend.t1_0;
    let t1_1 = &blend.t1_1;
    let t2_0 = &blend.t2_0;
    let t2_1 = &blend.t2_1;
    let zd = GpDir::from_xyz(&blend.axis_dir.xyz()).map_err(|e| e.to_string())?;

    let la = cache.seg(t1_0, t1_1); // tangency line on face 1
    let lb = cache.seg(t2_0, t2_1); // tangency line on face 2
    let ap0 = cache.arc(&blend.axis_origin, &zd, blend.radius, t1_0, t2_0)?; // arc at p0
    let ap1 = cache.arc(
        &blend.axis_origin.translated_vec(&blend.axis_dir.multiplied_scalar(
            blend.t1_1.distance(&blend.t1_0),
        )),
        &zd,
        blend.radius,
        t1_1,
        t2_1,
    )?; // arc at p1

    let wire = b.make_wire(&[la, ap1, lb, ap0]);
    Ok(b.make_face(surface, &[wire]))
}

/// Replace a sharp edge of `solid` with a constant-radius rolling-ball fillet.
///
/// The two adjacent faces must be planar; the end faces (perpendicular to the
/// edge, containing its endpoints) are rebuilt with the blend arc. All other
/// faces are kept untouched. Returns the rebuilt solid, or an error when the
/// edge is non-manifold, a face is non-planar, or the radius does not fit.
pub fn fillet_edge(solid: &TopoShape, edge: &Edge, radius: f64) -> Result<TopoShape, String> {
    FilletSpec { radius }.check()?;
    let adjacent = faces_touching_edge(solid, edge);
    if adjacent.len() != 2 {
        return Err(format!(
            "fillet: edge is adjacent to {} faces (expected 2)",
            adjacent.len()
        ));
    }
    let f1 = &adjacent[0];
    let f2 = &adjacent[1];
    let n1 = face_outward_normal(f1)?;
    let n2 = face_outward_normal(f2)?;

    let (p0, p1) = BRepTool::edge_vertices(edge)
        .ok_or("fillet: edge has no evaluable curve")?;
    let blend = blend_geometry(&p0, &p1, &n1, &n2, radius)?;

    // Retained half-planes for the two adjacent faces.
    let t1_0 = blend.t1_0;
    let t2_0 = blend.t2_0;
    let inside1 = |p: &GpPnt| GpVec::from_pnts(&t1_0, p).dot(&blend.u1_dir) >= -1e-12;
    let inside2 = |p: &GpPnt| GpVec::from_pnts(&t2_0, p).dot(&blend.u2_dir) >= -1e-12;

    // End faces: faces (other than the two adjacent) containing p0 / p1.
    let all_faces = faces_of(solid);
    let end0: Vec<Face> = all_faces
        .iter()
        .filter(|f| {
            face_contains_point(f, &p0)
                && !is_same(&f.0, &f1.0)
                && !is_same(&f.0, &f2.0)
        })
        .cloned()
        .collect();
    let end1: Vec<Face> = all_faces
        .iter()
        .filter(|f| {
            face_contains_point(f, &p1)
                && !is_same(&f.0, &f1.0)
                && !is_same(&f.0, &f2.0)
        })
        .cloned()
        .collect();
    if end0.len() != 1 || end1.len() != 1 {
        return Err(format!(
            "fillet: expected one end face per edge endpoint, got {} and {}",
            end0.len(),
            end1.len()
        ));
    }

    let mut cache = EdgeCache::new();
    let axis_dir = blend.axis_dir;
    let radius = blend.radius;

    let trimmed1 = rebuild_face(f1, &inside1, &mut cache)?;
    let trimmed2 = rebuild_face(f2, &inside2, &mut cache)?;
    let center0 = blend.axis_origin;
    let center1 = center0.translated_vec(&axis_dir.multiplied_scalar(blend.t1_1.distance(&blend.t1_0)));
    let end0_face = rebuild_end_face(
        &end0[0],
        &p0,
        &blend.t1_0,
        &blend.t2_0,
        &center0,
        &axis_dir,
        radius,
        &n1,
        &mut cache,
    )?;
    let end1_face = rebuild_end_face(
        &end1[0],
        &p1,
        &blend.t1_1,
        &blend.t2_1,
        &center1,
        &axis_dir,
        radius,
        &n1,
        &mut cache,
    )?;
    let blend_face = build_blend_face(&blend, &mut cache)?;

    // Assemble the new shell: keep every face except the replaced ones.
    let mut faces: Vec<Face> = Vec::new();
    for f in all_faces {
        if is_same(&f.0, &f1.0)
            || is_same(&f.0, &f2.0)
            || is_same(&f.0, &end0[0].0)
            || is_same(&f.0, &end1[0].0)
        {
            continue;
        }
        faces.push(f);
    }
    faces.push(trimmed1);
    faces.push(trimmed2);
    faces.push(end0_face);
    faces.push(end1_face);
    faces.push(blend_face);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Fillet several edges of `solid` sequentially (each step feeds the next).
/// `edge_indices` indexes into `edges_of(solid)`. Edges are filleted in the
/// given order; the first failure aborts with an error.
pub fn fillet_edge_chain(
    solid: &TopoShape,
    edge_indices: &[usize],
    radius: f64,
) -> Result<TopoShape, String> {
    let mut current = solid.clone();
    for &i in edge_indices {
        let es = edges_of(&current);
        let e = es
            .get(i)
            .ok_or_else(|| format!("fillet_edge_chain: edge index {i} out of range"))?;
        current = fillet_edge(&current, e, radius)?;
    }
    Ok(current)
}

/// Build a spherical-octant blend face for a box corner. `p` is the corner
/// point, `dirs` the three unit edge directions from `p` into the material,
/// `radius` the blend radius. The sphere is centred at `p` (an octant of a
/// sphere centred on the corner), which is the smallest single patch that
/// closes the trimmed faces.
///
/// ponytail: a rolling-ball corner patch that is tangent to all three faces
/// (a sphere centred at p + R·Σdirᵢ) cannot be closed by a single spherical
/// face — it meets each plane at a single point and needs the three edge
/// fillet cylinders. The corner here is therefore a spherical octant cut,
/// which forms a closed shell with exactly one added face.
///
/// The octant boundary is three quarter-circle arcs, each shared with one of
/// the three trimmed planar faces: `t_i = p + R·dirs[i]` and the arc between
/// `t_j`/`t_k` lies in the plane of the face perpendicular to `dirs[i]`.
fn build_corner_blend(p: &GpPnt, dirs: &[GpVec; 3], radius: f64, cache: &mut EdgeCache) -> Result<Face, String> {
    let d = [
        dirs[0].normalized(),
        dirs[1].normalized(),
        dirs[2].normalized(),
    ];
    let t = [
        p.translated_vec(&d[0].multiplied_scalar(radius)),
        p.translated_vec(&d[1].multiplied_scalar(radius)),
        p.translated_vec(&d[2].multiplied_scalar(radius)),
    ];
    let nd0 = GpDir::from_xyz(&d[0].xyz()).map_err(|e| e.to_string())?;
    let nd1 = GpDir::from_xyz(&d[1].xyz()).map_err(|e| e.to_string())?;
    let nd2 = GpDir::from_xyz(&d[2].xyz()).map_err(|e| e.to_string())?;
    let a01 = cache.arc(p, &nd2, radius, &t[0], &t[1])?; // in plane of face ⊥ d2
    let a12 = cache.arc(p, &nd0, radius, &t[1], &t[2])?; // in plane of face ⊥ d0
    let a20 = cache.arc(p, &nd1, radius, &t[2], &t[0])?; // in plane of face ⊥ d1
    let wire = TopoBuilder::new().make_wire(&[a01, a12, a20]);
    let mut sph = GpSphere::new(GpAx3::standard(), radius).map_err(|e| e.to_string())?;
    sph.set_location(*p);
    Ok(TopoBuilder::new().make_face(Arc::new(GeomSphere::new(sph)), &[wire]))
}

/// Fillet the three edges meeting at a box corner: replaces the corner with a
/// spherical-quadrant blend and trims the three adjacent planar faces.
/// `corner_vertex` must be a vertex of `solid` with exactly three incident
/// edges (a box corner). Returns a closed solid with one additional face.
pub fn fillet_corner_solid(
    solid: &TopoShape,
    corner_vertex: &Vertex,
    radius: f64,
) -> Result<TopoShape, String> {
    FilletSpec { radius }.check()?;
    let p = BRepTool::vertex_point(corner_vertex);

    // The three incident edges and their unit directions into the material.
    let es = edges_of(solid);
    let mut incident: Vec<Edge> = Vec::new();
    for e in &es {
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            let pa = vertex_position(&va);
            let pb = vertex_position(&vb);
            if pa.distance(&p) < 1e-9 || pb.distance(&p) < 1e-9 {
                incident.push(e.clone());
            }
        }
    }
    if incident.len() != 3 {
        return Err(format!(
            "fillet_corner_solid: vertex has {} incident edges (expected 3)",
            incident.len()
        ));
    }
    let mut dirs: Vec<GpVec> = Vec::with_capacity(3);
    for e in &incident {
        let (a, b) = edge_vertices(e);
        let (pa, pb) = (
            vertex_position(&a.ok_or("fillet: edge has no start vertex")?),
            vertex_position(&b.ok_or("fillet: edge has no end vertex")?),
        );
        let other = if pa.distance(&p) < 1e-9 { pb } else { pa };
        dirs.push(GpVec::from_pnts(&p, &other).normalized());
    }
    let dirs = [dirs[0], dirs[1], dirs[2]];

    // The three faces containing the corner.
    let all_faces = faces_of(solid);
    let corner_faces: Vec<Face> = all_faces
        .iter()
        .filter(|f| face_contains_point(f, &p))
        .cloned()
        .collect();
    if corner_faces.len() != 3 {
        return Err(format!(
            "fillet_corner_solid: corner touches {} faces (expected 3)",
            corner_faces.len()
        ));
    }

    // For each corner face, find its outward normal and its polygon. The two
    // in-plane directions from the corner are the two edges on the face.
    let mut cache = EdgeCache::new();
    let mut rebuilt: Vec<Face> = Vec::new();
    let mut idx_by_normal: Vec<(usize, GpVec)> = Vec::new();
    for f in &corner_faces {
        let n = face_outward_normal(f)?;
        // Match this face to the axis direction it is perpendicular to.
        let mut match_idx = usize::MAX;
        for (k, d) in dirs.iter().enumerate() {
            if n.xyz().crossed(&d.xyz()).modulus() < 1e-6 {
                match_idx = k;
                break;
            }
        }
        if match_idx == usize::MAX {
            return Err("fillet_corner_solid: cannot match corner face to axis".to_string());
        }
        idx_by_normal.push((match_idx, n));
    }
    // Build the retained polygon for each face: outside the sphere centred at p
    // with radius R (in the two in-plane coordinates).
    for (f, (k, n)) in corner_faces.iter().zip(idx_by_normal.iter()) {
        let poly = register_face_edges(f, &mut cache)?;
        let d_j = dirs[(k + 1) % 3];
        let d_k = dirs[(k + 2) % 3];
        // Retained: the point is outside the radius-R disk in the (d_j, d_k)
        // plane (i.e. its squared distance from p along the two in-plane axes
        // is at least R²).
        let inside = |q: &GpPnt| {
            let v = GpVec::from_pnts(&p, q);
            let a = v.dot(&d_j);
            let b = v.dot(&d_k);
            a * a + b * b >= radius * radius - 1e-9
        };
        let clipped = clip_polygon(&poly, &inside);
        if clipped.len() < 3 {
            return Err("fillet_corner_solid: clipping left no face".to_string());
        }
        // The clipped polygon has a chamfer corner at the two tangency points;
        // replace it with the spherical arc shared with the blend face.
        let wire = build_corner_face_wire(&p, &dirs, radius, &clipped, &mut cache, &n)?;
        let surf = BRepTool::face_surface(f).ok_or("fillet: face has no surface")?;
        rebuilt.push(cache.b.make_face(surf, &[wire]));
    }

    let blend_face = build_corner_blend(&p, &dirs, radius, &mut cache)?;

    let mut faces: Vec<Face> = Vec::new();
    for f in all_faces {
        if corner_faces.iter().any(|cf| is_same(&f.0, &cf.0)) {
            continue;
        }
        faces.push(f);
    }
    faces.extend(rebuilt);
    faces.push(blend_face);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Build the wire of a corner-trimmed planar face: the clipped polygon with its
/// straight chamfer edge (between the two tangency points) replaced by the
/// spherical arc shared with the corner blend face.
///
/// `n` is the face's outward normal; the face is the one perpendicular to the
/// axis direction `dirs[k]` that is parallel to `n`. `dirs` are the three
/// in-wedge directions from the corner.
fn build_corner_face_wire(
    p: &GpPnt,
    dirs: &[GpVec; 3],
    radius: f64,
    clipped: &[GpPnt],
    cache: &mut EdgeCache,
    n: &GpVec,
) -> Result<Wire, String> {
    let k = (0..3)
        .find(|&i| dirs[i].xyz().crossed(&n.xyz()).modulus() < 1e-6)
        .ok_or("fillet: cannot classify corner face")?;
    let d = [
        dirs[0].normalized(),
        dirs[1].normalized(),
        dirs[2].normalized(),
    ];
    let t_a = p.translated_vec(&d[(k + 1) % 3].multiplied_scalar(radius));
    let t_b = p.translated_vec(&d[(k + 2) % 3].multiplied_scalar(radius));
    // The spherical arc shared with the blend face: in this face's plane
    // (normal d[k]), from t_a to t_b.
    let nd = GpDir::from_xyz(&d[k].xyz()).map_err(|e| e.to_string())?;
    let arc = cache.arc(p, &nd, radius, &t_a, &t_b)?;

    // Find the chamfer edge of the clipped polygon: the straight segment
    // connecting t_a and t_b (either orientation).
    let m = clipped.len();
    let mut chamfer = None;
    for i in 0..m {
        let a = &clipped[i];
        let c = &clipped[(i + 1) % m];
        let is_ab = (a.distance(&t_a) < 1e-6 && c.distance(&t_b) < 1e-6)
            || (a.distance(&t_b) < 1e-6 && c.distance(&t_a) < 1e-6);
        if is_ab {
            chamfer = Some(i);
            break;
        }
    }
    let ci = chamfer.ok_or("fillet: corner polygon has no chamfer edge")?;

    // Replace the chamfer edge with the arc (which shares its endpoint pair,
    // so the cache returns the same instance the blend face uses).
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    for i in 0..m {
        if i == ci {
            edges.push(arc.clone());
        } else {
            edges.push(cache.seg(&clipped[i], &clipped[(i + 1) % m]));
        }
    }
    Ok(b.make_wire(&edges))
}

/// Classify a surface, adding cylinder detection on top of
/// `brep_surface::classify_surface` (which only distinguishes Plane / Sphere /
/// Other). A cylinder is recognised by constant distance from an axis line.
pub fn classify_surface_full(s: &dyn Surface) -> SurfaceKind {
    let k = classify_surface(s);
    if k != SurfaceKind::Other {
        return k;
    }
    if cylinder_radius(s).is_some() {
        SurfaceKind::Cylinder
    } else {
        SurfaceKind::Other
    }
}

/// If `s` is a cylinder, return its radius. All sampled surface normals must be
/// perpendicular to one axis direction, and all sampled points equidistant
/// from that axis line.
pub fn cylinder_radius(s: &dyn Surface) -> Option<f64> {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let nu = 6;
    let nv = 6;
    let hu = ((u1 - u0) * 1e-4).max(1e-7);
    let hv = ((v1 - v0) * 1e-4).max(1e-7);
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut nrm: Vec<GpVec> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            pts.push(s.d0(u, v));
            // Local finite-difference normal: `brep_surface::surface_normal`
            // falls back to an infinite step on unbounded parameter ranges
            // (e.g. the cylinder's infinite v range), so compute it here with
            // the clamped sampling bounds.
            let p0 = s.d0(u, v);
            let du = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
            let dv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
            let nn = du.xyz().crossed(&dv.xyz());
            let m = nn.modulus();
            nrm.push(if m > 1e-30 { GpVec::from_xyz(&nn.divided(m)) } else { GpVec::zero() });
        }
    }
    // Axis direction: perpendicular to all normals. Pick two non-parallel
    // normals and cross them.
    let mut axis = None;
    for a in 0..nrm.len() {
        for b in (a + 1)..nrm.len() {
            let c = nrm[a].xyz().crossed(&nrm[b].xyz());
            if c.modulus() > 1e-6 {
                axis = Some(GpVec::from_xyz(&c).normalized());
                break;
            }
        }
        if axis.is_some() {
            break;
        }
    }
    let axis = axis?;
    // All normals must be perpendicular to the axis (dot product ≈ 0).
    for nn in &nrm {
        if nn.dot(&axis).abs() > 1e-3 {
            return None;
        }
    }
    // Project points onto the plane perpendicular to the axis; the axis line
    // passes through the 2D circumcenter of three non-collinear projections.
    let (x2, y2) = project_basis(&axis);
    let proj: Vec<(f64, f64)> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&pts[0], q);
            (v.dot(&x2), v.dot(&y2))
        })
        .collect();
    let (mut cx, mut cy) = (f64::NAN, f64::NAN);
    'outer: for a in 0..proj.len() {
        for b in (a + 1)..proj.len() {
            for c in (b + 1)..proj.len() {
                let (pa, pb, pc) = (proj[a], proj[b], proj[c]);
                let d = 2.0 * (pa.0 * (pb.1 - pc.1) + pb.0 * (pc.1 - pa.1) + pc.0 * (pa.1 - pb.1));
                if d.abs() < 1e-12 {
                    continue;
                }
                let ux = ((pa.0 * pa.0 + pa.1 * pa.1) * (pb.1 - pc.1)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pc.1 - pa.1)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pa.1 - pb.1))
                    / d;
                let uy = ((pa.0 * pa.0 + pa.1 * pa.1) * (pc.0 - pb.0)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pa.0 - pc.0)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pb.0 - pa.0))
                    / d;
                cx = ux;
                cy = uy;
                break 'outer;
            }
        }
    }
    if !cx.is_finite() {
        return None;
    }
    let axis_pt = pts[0].translated_vec(&x2.multiplied_scalar(cx).added(&y2.multiplied_scalar(cy)));
    let mut radii: Vec<f64> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&axis_pt, q);
            let along = v.dot(&axis);
            let perp = v.subtracted(&axis.multiplied_scalar(along));
            perp.magnitude()
        })
        .collect();
    if radii.is_empty() {
        return None;
    }
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = radii[radii.len() / 2];
    if med < 1e-9 {
        return None;
    }
    for r in &radii {
        if (r - med).abs() > 1e-3 * med.abs().max(1.0) {
            return None;
        }
    }
    Some(med)
}

/// Two orthonormal vectors spanning the plane perpendicular to `axis`.
fn project_basis(axis: &GpVec) -> (GpVec, GpVec) {
    let ref_v = if axis.x().abs() < 0.9 {
        GpVec::new(1.0, 0.0, 0.0)
    } else {
        GpVec::new(0.0, 1.0, 0.0)
    };
    let x2 = axis.crossed(&ref_v).normalized();
    let y2 = axis.crossed(&x2).normalized();
    (x2, y2)
}

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_surface::classify_surface;
    use crate::primitives::BRepPrimBox;
    use crate::shape::Shell;
    use crate::shell_check::shell_is_closed;
    use crate::topo_tools_full::vertices_of;

    fn box_edges(b: &BRepPrimBox) -> Vec<Edge> {
        edges_of(&b.solid.0)
    }

    fn find_face_by_surface<'a>(
        faces: &'a [Face],
        f: impl Fn(&dyn Surface) -> bool,
    ) -> Option<&'a Face> {
        faces
            .iter()
            .find(|fa| BRepTool::face_surface(fa).map(|s| f(s.as_ref())).unwrap_or(false))
    }

    #[test]
    fn corner_wire_fillet_alias() {
        // The 2D corner fillet module and this 3D edge fillet coexist without
        // a name clash: `fillet::fillet_corner` (wire) vs `fillet_edge` here.
        assert!(FilletSpec { radius: 0.5 }.check().is_ok());
        assert!(FilletSpec { radius: 0.0 }.check().is_err());
    }

    #[test]
    fn fillet_box_edge_90deg() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        // Edge from (0,0,0) → (2,0,0): the X edge at y=0, z=0.
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .expect("find X edge");
        let out = fillet_edge(&bx.solid.0, e, 0.4).expect("fillet");
        let faces = faces_of(&out);
        // 6 → 7 faces (blend cylinder added).
        assert_eq!(faces.len(), 7, "face count {}", faces.len());
        // Blend face classifies as a cylinder.
        let blend = find_face_by_surface(&faces, |s| classify_surface_full(s) == SurfaceKind::Cylinder)
            .expect("blend face");
        let r = cylinder_radius(BRepTool::face_surface(blend).unwrap().as_ref()).unwrap();
        assert!((r - 0.4).abs() < 1e-6, "blend radius {r}");
        // Closed shell.
        let shell = Shell(out.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "filleted box is closed");
        // Every edge endpoint of the result lies on the original box surface
        // (within tolerance). The box spans [0,2]³; the surface is the union
        // of the six coordinate planes at 0 and 2.
        let on_box_surface = |p: &GpPnt| {
            let d = [
                p.x(), p.x() - 2.0, p.y(), p.y() - 2.0, p.z(), p.z() - 2.0,
            ];
            d.iter().map(|v| v.abs()).fold(f64::INFINITY, f64::min) < 1e-6
        };
        for edge in edges_of(&out) {
            let (a, b) = BRepTool::edge_vertices(&edge).unwrap();
            assert!(on_box_surface(&a), "endpoint {a:?} off box surface");
            assert!(on_box_surface(&b), "endpoint {b:?} off box surface");
        }
    }

    #[test]
    fn fillet_increases_face_count() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let out = fillet_edge(&bx.solid.0, e, 0.3).expect("fillet");
        assert_eq!(faces_of(&out).len(), 7);
    }

    #[test]
    fn fillet_blend_surface_is_cylinder() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let out = fillet_edge(&bx.solid.0, e, 0.4).unwrap();
        let faces = faces_of(&out);
        let blend = find_face_by_surface(&faces, |s| classify_surface_full(s) == SurfaceKind::Cylinder)
            .expect("blend face");
        // The vanilla classifier cannot see cylinders; the extended one can.
        let vanilla = classify_surface(BRepTool::face_surface(blend).unwrap().as_ref());
        assert_eq!(vanilla, SurfaceKind::Other);
        assert_eq!(classify_surface_full(BRepTool::face_surface(blend).unwrap().as_ref()), SurfaceKind::Cylinder);
    }

    #[test]
    fn fillet_chain_two_edges() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        // Edge 0: (0,0,0)→(2,0,0). Edge at index of the opposite top edge
        // (2,2,2)→(0,2,2): find by endpoints.
        let e0 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let e1 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                (a.is_equal(&GpPnt::new(2.0, 2.0, 2.0)) && b.is_equal(&GpPnt::new(0.0, 2.0, 2.0)))
                    || (b.is_equal(&GpPnt::new(2.0, 2.0, 2.0)) && a.is_equal(&GpPnt::new(0.0, 2.0, 2.0)))
            })
            .unwrap();
        let out = fillet_edge_chain(&bx.solid.0, &[e0, e1], 0.3).expect("chain");
        assert_eq!(faces_of(&out).len(), 8, "two fillets add two faces");
        let shell = Shell(out.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "chained fillet is closed");
    }

    #[test]
    fn fillet_corner_solid_sphere_blend() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let verts = vertices_of(&bx.solid.0);
        let corner = verts
            .iter()
            .find(|v| BRepTool::vertex_point(v).is_equal(&GpPnt::new(0.0, 0.0, 0.0)))
            .expect("corner vertex");
        let out = fillet_corner_solid(&bx.solid.0, corner, 0.4).expect("corner fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 7, "corner fillet adds one face");
        let sphere = find_face_by_surface(&faces, |s| classify_surface(s) == SurfaceKind::Sphere)
            .expect("sphere blend face");
        let _ = sphere;
        let shell = Shell(out.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "corner fillet is closed");
    }

    #[test]
    fn nonplanar_edge_errors() {
        let cyl = crate::primitives::BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let edges = edges_of(&cyl.solid.0);
        // The cap circles are adjacent to a planar cap and the lateral face;
        // the seam is adjacent to the lateral (non-planar) face on both sides.
        // Filleting the seam must fail because both adjacent faces are the
        // same non-planar lateral surface.
        let seam = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.x().abs() > 0.99 && b.x().abs() > 0.99 && (a.z() - b.z()).abs() > 0.1
            })
            .unwrap();
        assert!(fillet_edge(&cyl.solid.0, seam, 0.2).is_err());
    }

    #[test]
    fn radius_positive_required() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = &edges[0];
        assert!(fillet_edge(&bx.solid.0, e, 0.0).is_err());
        assert!(fillet_edge(&bx.solid.0, e, -1.0).is_err());
    }

    #[test]
    fn vertex_endpoint_geometry_preserved() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let out = fillet_edge(&bx.solid.0, e, 0.4).unwrap();
        // The two far endpoints of the filleted edge are the +Y/+Z corners on
        // the original surface; the tangency vertices lie on the adjacent
        // planes. Check the far corner (2,2,0) and (0,2,0) remain present.
        let verts = vertices_of(&out);
        let ps: Vec<GpPnt> = verts.iter().map(BRepTool::vertex_point).collect();
        for want in [
            GpPnt::new(2.0, 2.0, 0.0),
            GpPnt::new(0.0, 2.0, 0.0),
            GpPnt::new(2.0, 0.0, 2.0),
            GpPnt::new(0.0, 0.0, 2.0),
        ] {
            assert!(
                ps.iter().any(|p| p.distance(&want) < 1e-6),
                "missing preserved vertex {want:?}"
            );
        }
        // The original sharp-edge endpoints are replaced by tangency points on
        // the adjacent planes (z=0 and y=0 for this edge).
        assert!(
            ps.iter().any(|p| (p.distance(&GpPnt::new(0.4, 0.0, 0.0)) < 1e-6)
                || (p.distance(&GpPnt::new(0.0, 0.4, 0.0)) < 1e-6)),
            "tangency vertices present"
        );
    }
}
