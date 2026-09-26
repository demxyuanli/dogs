use super::prelude::*;


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
pub(super) struct EdgeCache {
    pub(super) b: TopoBuilder,
    pub(super) pts: Vec<GpPnt>,
    pub(super) segs: HashMap<(usize, usize), Edge>,
    pub(super) arcs: HashMap<(usize, usize), Edge>,
}

impl EdgeCache {
    pub(super) fn new() -> Self {
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
    pub(super) fn point_index(&mut self, p: GpPnt) -> usize {
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
    pub(super) fn register_original(&mut self, a: &GpPnt, b: &GpPnt, e: &Edge) {
        let ia = self.point_index(*a);
        let ib = self.point_index(*b);
        self.segs.insert((ia.min(ib), ia.max(ib)), e.clone());
    }

    /// A straight segment edge from `a` to `b` (shared by all requesters).
    pub(super) fn seg(&mut self, a: &GpPnt, b: &GpPnt) -> Edge {
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
    pub(super) fn arc(
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
        self.b.add_edge_vertices(&mut e, &v1, &v2);
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
pub(super) fn in_wedge_directions(a: &GpVec, n1: &GpVec, n2: &GpVec) -> Result<(GpVec, GpVec, f64), String> {
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
pub(super) struct BlendGeom {
    /// Tangency point on face 1 at the two edge endpoints.
    pub(super) t1_0: GpPnt,
    pub(super) t1_1: GpPnt,
    /// Tangency point on face 2 at the two edge endpoints.
    pub(super) t2_0: GpPnt,
    pub(super) t2_1: GpPnt,
    /// Cylinder axis (a point on the axis at the p0 cross-section, direction).
    pub(super) axis_origin: GpPnt,
    pub(super) axis_dir: GpVec,
    /// Cylinder frame: X points at t1_0 (so u = 0 → t1_0), Y chosen so the
    /// angular sweep to t2_0 is positive.
    pub(super) frame: GpAx3,
    /// Face 1 in-wedge direction (into the material corner).
    pub(super) u1_dir: GpVec,
    /// Face 2 in-wedge direction.
    pub(super) u2_dir: GpVec,
    pub(super) radius: f64,
}

/// Build the blend geometry for a fillet of `radius` along the straight edge
/// `p0 → p1`, with the two adjacent planar faces' outward normals `n1`/`n2`.
pub(super) fn blend_geometry(
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
pub(super) fn face_contains_point(face: &Face, p: &GpPnt) -> bool {
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
pub(super) fn face_polygon(face: &Face) -> Result<Vec<GpPnt>, String> {
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
pub(super) fn clip_polygon(poly: &[GpPnt], inside: &impl Fn(&GpPnt) -> bool) -> Vec<GpPnt> {
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
pub(super) fn dedupe_polygon(poly: &[GpPnt]) -> Vec<GpPnt> {
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
pub(super) fn intersect_t(a: &GpPnt, b: &GpPnt, inside: &impl Fn(&GpPnt) -> bool) -> Option<f64> {
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
pub(super) fn build_polygon_wire(cache: &mut EdgeCache, poly: &[GpPnt]) -> Result<Wire, String> {
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
pub(super) fn register_face_edges(face: &Face, cache: &mut EdgeCache) -> Result<Vec<GpPnt>, String> {
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
pub(super) fn face_outward_normal(face: &Face) -> Result<GpVec, String> {
    let pln = face_plane(face).ok_or("fillet: face is not planar")?;
    let n = *pln.axis().direction().xyz();
    Ok(GpVec::from_xyz(&n))
}

/// Rebuild one planar face into its trimmed version: the retained region is
/// `inside(p)` (a half-plane for adjacent faces, the outside-of-sphere test for
/// corner faces). Original edges that survive untouched are reused so shared
/// boundaries with kept neighbour faces stay intact.
pub(super) fn rebuild_face(face: &Face, inside: &impl Fn(&GpPnt) -> bool, cache: &mut EdgeCache) -> Result<Face, String> {
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
pub(super) fn rebuild_end_face(
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
pub(super) fn build_blend_face(blend: &BlendGeom, cache: &mut EdgeCache) -> Result<Face, String> {
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
