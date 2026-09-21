use super::prelude::*;
use super::*;

/// Default number of samples along the filleted edge for the variable-radius
/// blend (mirrors the `NbSamples` default used by the OCCT fillet builder).

pub const DEFAULT_VAR_SAMPLES: usize = 17;

/// Cross-section arc samples per edge sample (the blend b-spline's v-direction).
pub(super) const ARC_SAMPLES: usize = 5;

/// Radius law along the filleted edge (`Law_Function`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RadiusLaw {
    /// r(t) = r0 + (r1 − r0)·t
    Linear,
    /// r(t) = r0 + (r1 − r0)·(2t − t²), monotone with zero slope at t = 1.
    Quadratic,
    /// r(t) = r0 + (r1 − r0)·(3t² − 2t³) — the smoothstep, C¹ with zero
    /// slope at both ends.
    Cubic,
}

/// Variable-radius fillet parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VarFilletSpec {
    /// Ball radius at the edge's first parameter (t = 0).
    pub r_start: f64,
    /// Ball radius at the edge's last parameter (t = 1).
    pub r_end: f64,
    /// Interpolation law between the two.
    pub law: RadiusLaw,
}

impl VarFilletSpec {
    /// Linear law between `r_start` and `r_end`.
    pub fn new(r_start: f64, r_end: f64) -> Self {
        Self { r_start, r_end, law: RadiusLaw::Linear }
    }

    /// Validate both radii are strictly positive and finite.
    pub fn check(&self) -> Result<(), String> {
        if self.r_start <= 0.0 || !self.r_start.is_finite() {
            return Err("fillet_var: r_start must be a positive finite value".to_string());
        }
        if self.r_end <= 0.0 || !self.r_end.is_finite() {
            return Err("fillet_var: r_end must be a positive finite value".to_string());
        }
        Ok(())
    }
}

/// Ball radius at edge fraction `t ∈ [0, 1]` under the spec's law.
pub fn radius_at(spec: &VarFilletSpec, t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    let d = spec.r_end - spec.r_start;
    match spec.law {
        RadiusLaw::Linear => spec.r_start + d * t,
        RadiusLaw::Quadratic => spec.r_start + d * (2.0 * t - t * t),
        RadiusLaw::Cubic => spec.r_start + d * (3.0 * t * t - 2.0 * t * t * t),
    }
}

/// Per-sample blend geometry along the filleted edge.
#[derive(Debug, Clone)]
pub(super) struct BlendSample {
    /// Ball radius at this sample.
    pub(super) radius: f64,
    /// Blend-circle centre (the ball centre at this cross-section).
    pub(super) center: GpPnt,
    /// Tangency point on adjacent face 1.
    pub(super) t1: GpPnt,
    /// Tangency point on adjacent face 2.
    pub(super) t2: GpPnt,
    /// Sampled points on the cross-section arc (t1 → t2), for the b-spline grid.
    pub(super) arc_points: Vec<GpPnt>,
}

/// Compute the tangent-plane in-wedge directions for two planar faces meeting
/// at a straight edge (copied from `fillet_edge`, which owns the canonical
/// implementation).
pub(super) fn in_wedge_directions(a: &GpVec, n1: &GpVec, n2: &GpVec) -> Result<(GpVec, GpVec, f64), String> {
    let b1 = a.crossed(n1).normalized();
    let b2 = a.crossed(n2).normalized();
    let mag = b1.xyz().crossed(&b2.xyz()).modulus();
    if mag < 1e-12 {
        return Err("fillet: adjacent faces are coplanar (no wedge)".to_string());
    }
    let s1 = if b1.dot(n2) < 0.0 { 1.0 } else { -1.0 };
    let s2 = if b2.dot(n1) < 0.0 { 1.0 } else { -1.0 };
    let u1 = b1.multiplied_scalar(s1).normalized();
    let u2 = b2.multiplied_scalar(s2).normalized();
    let cos_t = u1.dot(&u2).clamp(-1.0, 1.0);
    let sin_t = u1.xyz().crossed(&u2.xyz()).modulus();
    if sin_t < 1e-9 {
        return Err("fillet: degenerate wedge angle".to_string());
    }
    if cos_t <= 0.0 && sin_t < 1e-6 {
        return Err("fillet: reflex (concave) corner unsupported".to_string());
    }
    let theta = cos_t.acos();
    Ok((u1, u2, theta))
}

/// Builds and caches topology edges so faces that must share a boundary edge
/// all receive the *same* `Edge` instance (required for `shell_is_closed`).
/// Copied from `fillet_edge` (see the notes there).
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

    pub(super) fn point_index(&mut self, p: GpPnt) -> usize {
        for (i, q) in self.pts.iter().enumerate() {
            if p.distance(q) < 1e-6 {
                return i;
            }
        }
        self.pts.push(p);
        self.pts.len() - 1
    }

    pub(super) fn register_original(&mut self, a: &GpPnt, b: &GpPnt, e: &Edge) {
        let ia = self.point_index(*a);
        let ib = self.point_index(*b);
        self.segs.insert((ia.min(ib), ia.max(ib)), e.clone());
    }

    pub(super) fn seg(&mut self, a: &GpPnt, b: &GpPnt) -> Edge {
        if a.distance(b) < 1e-9 {
            panic!("fillet_var: degenerate segment at {:?}", a);
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
            .map_err(|_| "fillet_var: arc start coincides with the arc centre".to_string())?;
        let ax2 = GpAx2::new(*center, *normal, xd)
            .map_err(|e| format!("fillet_var: arc frame is degenerate: {e}"))?;
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

/// The ordered boundary polygon of a planar face's outer wire (copied from
/// `fillet_edge`).
pub(super) fn face_polygon(face: &Face) -> Result<Vec<GpPnt>, String> {
    let w = wires_of_face(face)
        .into_iter()
        .next()
        .ok_or("fillet_var: face has no boundary wire")?;
    let es = edges_of_wire(&w);
    if es.is_empty() {
        return Err("fillet_var: face wire has no edges".to_string());
    }
    let ends: Vec<(GpPnt, GpPnt)> = es
        .iter()
        .map(|e| {
            let (a, b) = edge_vertices(e);
            let a = a.ok_or("fillet_var: edge has no start vertex")?;
            let b = b.ok_or("fillet_var: edge has no end vertex")?;
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

/// Sutherland–Hodgman clip of `poly` against the half-plane `inside(p)`
/// (copied from `fillet_edge`; works for any predicate, including the curved
/// tangency boundary used here).
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
    while out.len() >= 2 && out[0].distance(out.last().unwrap()) < 1e-9 {
        out.pop();
    }
    out
}

/// Bisection search for the fraction `t ∈ [0, 1]` where the segment `a→b`
/// crosses the `inside` boundary.
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

/// Register the original edges of `face` with the cache so future segment
/// requests reuse them. Returns the face's polygon alongside.
pub(super) fn register_face_edges(face: &Face, cache: &mut EdgeCache) -> Result<Vec<GpPnt>, String> {
    let poly = face_polygon(face)?;
    let w = wires_of_face(face)
        .into_iter()
        .next()
        .ok_or("fillet_var: face has no wire")?;
    let es = edges_of_wire(&w);
    for e in &es {
        let (a, b) = edge_vertices(e);
        let pa = vertex_position(&a.ok_or("fillet_var: edge has no start vertex")?);
        let pb = vertex_position(&b.ok_or("fillet_var: edge has no end vertex")?);
        cache.register_original(&pa, &pb, e);
    }
    Ok(poly)
}

/// The outward unit normal of a planar face.
pub(super) fn face_outward_normal(face: &Face) -> Result<GpVec, String> {
    let pln = face_plane(face).ok_or("fillet_var: face is not planar")?;
    let n = *pln.axis().direction().xyz();
    Ok(GpVec::from_xyz(&n))
}

/// Whether `p` coincides with any vertex of `face`'s boundary wires.
pub(super) fn face_contains_point(face: &Face, p: &GpPnt) -> bool {
    wires_of_face(face).iter().any(|w| {
        edges_of_wire(w).iter().any(|e| {
            let (a, b) = edge_vertices(e);
            a.map_or(false, |v| vertex_position(&v).distance(p) < 1e-9)
                || b.map_or(false, |v| vertex_position(&v).distance(p) < 1e-9)
        })
    })
}

/// Sample the cross-section arc at sample `i`: points on the circle of
/// `radius` centred at `center`, in the plane perpendicular to `a`, from the
/// tangency point `t1` to `t2`. `m` points, uniform in arc angle.
pub(super) fn arc_sample_points(
    center: &GpPnt,
    a: &GpVec,
    radius: f64,
    t1: &GpPnt,
    t2: &GpPnt,
    m: usize,
) -> Result<Vec<GpPnt>, String> {
    if m < 2 {
        return Err("fillet_var: need at least 2 arc samples".to_string());
    }
    let xd_vec = GpVec::from_pnts(center, t1).normalized();
    let xd = GpDir::from_xyz(&xd_vec.xyz()).map_err(|_| "fillet_var: degenerate arc frame")?;
    let zd = GpDir::from_xyz(&a.xyz()).map_err(|_| "fillet_var: degenerate arc axis")?;
    let ax2 = GpAx2::new(*center, zd, xd).map_err(|e| format!("fillet_var: arc frame: {e}"))?;
    let xdir = *ax2.x_direction();
    let ydir = *ax2.y_direction();
    let va = GpVec::from_pnts(center, t1);
    let vb = GpVec::from_pnts(center, t2);
    let a1 = va.coord.dot(&ydir.xyz()).atan2(va.coord.dot(&xdir.xyz()));
    let a2 = vb.coord.dot(&ydir.xyz()).atan2(vb.coord.dot(&xdir.xyz()));
    let mut pts = Vec::with_capacity(m);
    for j in 0..m {
        let s = j as f64 / (m - 1) as f64;
        let phi = a1 + s * (a2 - a1);
        let xv = GpVec::from_xyz(&xdir.xyz()).multiplied_scalar(radius * phi.cos());
        let yv = GpVec::from_xyz(&ydir.xyz()).multiplied_scalar(radius * phi.sin());
        pts.push(center.translated_vec(&xv.added(&yv)));
    }
    Ok(pts)
}

/// Sample the variable-radius blend geometry along a straight edge `p0 → p1`
/// with adjacent planar-face normals `n1`/`n2`. Returns the samples plus the
/// in-wedge directions, the tangency cotangent factor and the edge axis.
#[allow(clippy::type_complexity)]
pub(super) fn sample_blend(
    p0: &GpPnt,
    p1: &GpPnt,
    n1: &GpVec,
    n2: &GpVec,
    spec: &VarFilletSpec,
    n: usize,
) -> Result<(Vec<BlendSample>, GpVec, GpVec, f64), String> {
    let n = n.max(2);
    let edge_len = p0.distance(p1);
    if edge_len < 1e-12 {
        return Err("fillet_var: edge has zero length".to_string());
    }
    let a = GpVec::from_pnts(p0, p1).normalized();
    let (u1, u2, theta) = in_wedge_directions(&a, n1, n2)?;
    let sin_t = theta.sin();
    if sin_t < 1e-9 {
        return Err("fillet_var: degenerate wedge".to_string());
    }
    let cot_half = (1.0 + theta.cos()) / sin_t;
    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let r = radius_at(spec, t);
        let p = p0.translated_vec(&a.multiplied_scalar(t * edge_len));
        let off = r * cot_half;
        let t1 = p.translated_vec(&u1.multiplied_scalar(off));
        let t2 = p.translated_vec(&u2.multiplied_scalar(off));
        let center = p.translated_vec(&u1.added(&u2).multiplied_scalar(r / sin_t));
        let arc_points = arc_sample_points(&center, &a, r, &t1, &t2, ARC_SAMPLES)?;
        samples.push(BlendSample { radius: r, center, t1, t2, arc_points });
    }
    Ok((samples, u1, u2, cot_half))
}

/// Rebuild one planar adjacent face into its variable-radius trimmed version:
/// the retained region is `inside(p)` (the far side of the curved tangency
/// boundary). The straight chamfer edge between the first and last tangency
/// points is replaced by the sampled tangency polyline, whose segments are
/// shared with the blend face.
pub(super) fn rebuild_var_face(
    face: &Face,
    inside: &impl Fn(&GpPnt) -> bool,
    polyline_first: &GpPnt,
    polyline_last: &GpPnt,
    polyline: &[GpPnt],
    tol: f64,
    cache: &mut EdgeCache,
) -> Result<Face, String> {
    let poly = register_face_edges(face, cache)?;
    let clipped = clip_polygon(&poly, inside);
    if clipped.len() < 3 {
        return Err("fillet_var: clipping a face left no retained region".to_string());
    }
    // The clip boundary is evaluated with an `eps` slack, so the crossing
    // points sit up to `eps` inside the true tangency curve. The chamfer
    // search tolerance must exceed that slack; it stays far below the spacing
    // between distinct tangency points.
    let chamfer_tol = (10.0 * tol).max(1e-4);
    let m = clipped.len();
    let mut ci = None;
    for i in 0..m {
        let a = &clipped[i];
        let c = &clipped[(i + 1) % m];
        let is_ab = (a.distance(polyline_first) < chamfer_tol && c.distance(polyline_last) < chamfer_tol)
            || (a.distance(polyline_last) < chamfer_tol && c.distance(polyline_first) < chamfer_tol);
        if is_ab {
            ci = Some(i);
            break;
        }
    }
    let ci = ci.ok_or("fillet_var: trimmed face has no tangency chamfer edge")?;
    let first_is_first = clipped[ci].distance(polyline_first) < chamfer_tol;
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    for i in 0..m {
        if i == ci {
            if first_is_first {
                for w in polyline.windows(2) {
                    edges.push(cache.seg(&w[0], &w[1]));
                }
            } else {
                for w in polyline.windows(2).rev() {
                    edges.push(cache.seg(&w[1], &w[0]));
                }
            }
        } else {
            edges.push(cache.seg(&clipped[i], &clipped[(i + 1) % m]));
        }
    }
    let wire = b.make_wire(&edges);
    let surf = BRepTool::face_surface(face).ok_or("fillet_var: face has no surface")?;
    Ok(b.make_face(surf, &[wire]))
}

/// Rebuild the end face at an edge endpoint: replace the corner vertex `p`
/// with the blend arc between the tangency points `t1`/`t2`, shared with the
/// blend face (copied from `fillet_edge`).
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
        .ok_or("fillet_var: end face does not contain the edge endpoint")?;
    let prev = &poly[(idx + n - 1) % n];
    let next = &poly[(idx + 1) % n];

    let dir_prev = GpVec::from_pnts(prev, p).normalized();
    let dir_next = GpVec::from_pnts(p, next).normalized();
    let prev_is_face1 = dir_prev.dot(n1).abs() < 1e-6;
    let next_is_face1 = dir_next.dot(n1).abs() < 1e-6;

    let (first_seg_a, first_seg_b, last_seg_a, last_seg_b);
    if prev_is_face1 && !next_is_face1 {
        first_seg_a = prev;
        first_seg_b = t1;
        last_seg_a = t2;
        last_seg_b = next;
    } else if next_is_face1 && !prev_is_face1 {
        first_seg_a = prev;
        first_seg_b = t2;
        last_seg_a = t1;
        last_seg_b = next;
    } else {
        return Err("fillet_var: cannot classify end-face incident edges".to_string());
    }

    let zd = GpDir::from_xyz(&axis_dir.xyz()).map_err(|e| e.to_string())?;
    let arc_edge = cache.arc(center, &zd, radius, t1, t2)?;

    let b = TopoBuilder::new();
    let mut wire_edges: Vec<Edge> = Vec::new();
    wire_edges.push(cache.seg(first_seg_a, first_seg_b));
    wire_edges.push(arc_edge.clone());
    wire_edges.push(cache.seg(last_seg_a, last_seg_b));
    let mut j = (idx + 1) % n;
    let mut guard = 0;
    while j != idx && guard < n {
        if j == (idx + n - 1) % n {
            break;
        }
        let a = &poly[j];
        let c = &poly[(j + 1) % n];
        wire_edges.push(cache.seg(a, c));
        j = (j + 1) % n;
        guard += 1;
    }
    let wire = b.make_wire(&wire_edges);
    let surf = BRepTool::face_surface(face).ok_or("fillet_var: face has no surface")?;
    Ok(b.make_face(surf, &[wire]))
}

/// Build the variable-radius blend face: a degree-1 tensor-product b-spline
/// through the sampled cross-section arcs, bounded by the two tangency
/// polylines and the two end arcs (all shared with the neighbouring faces).
pub(super) fn build_var_blend_face(samples: &[BlendSample], axis_dir: &GpVec, cache: &mut EdgeCache) -> Result<Face, String> {
    let n = samples.len();
    if n < 2 {
        return Err("fillet_var: need at least 2 edge samples".to_string());
    }
    let b = TopoBuilder::new();
    let zd = GpDir::from_xyz(&axis_dir.xyz()).map_err(|e| e.to_string())?;

    let mut edges: Vec<Edge> = Vec::new();
    for i in 0..n - 1 {
        edges.push(cache.seg(&samples[i].t1, &samples[i + 1].t1));
    }
    let last = n - 1;
    edges.push(cache.arc(
        &samples[last].center,
        &zd,
        samples[last].radius,
        &samples[last].t1,
        &samples[last].t2,
    )?);
    for i in (1..n).rev() {
        edges.push(cache.seg(&samples[i].t2, &samples[i - 1].t2));
    }
    edges.push(cache.arc(
        &samples[0].center,
        &zd,
        samples[0].radius,
        &samples[0].t1,
        &samples[0].t2,
    )?);
    let wire = b.make_wire(&edges);

    let poles: Vec<Vec<GpPnt>> = samples.iter().map(|s| s.arc_points.clone()).collect();
    let surf = fit_surface_grid(&poles, 1, 1)
        .map_err(|e| format!("fillet_var: blend surface fit failed: {e}"))?;
    let surface: Arc<dyn Surface> = Arc::new(surf);
    Ok(b.make_face(surface, &[wire]))
}
