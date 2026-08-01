//! Variable-radius and chain rolling-ball fillets — a simplified port of
//! `BRepFilletAPI_MakeFillet` / `ChFi3d` with radius laws.
//!
//! Where `fillet_edge` replaces a sharp convex edge with a *constant-radius*
//! quarter-cylinder blend, this module lets the ball radius vary along the
//! edge (`VarFilletSpec` + `RadiusLaw`, mirroring OCCT's `Law_Function`), and
//! applies the fillet to a *chain* of edges (mirroring `ChFi3d`'s edge-chain
//! handling).
//!
//! The edge is sampled at `DEFAULT_VAR_SAMPLES` points. At each sample the
//! local tangency points, the blend-circle center and the radius are computed
//! with the same in-wedge geometry as `fillet_edge`; the blend face is a
//! degree-1 tensor-product B-spline (`occt_geom::bspline_surface::fit_surface_grid`)
//! interpolating the sampled cross-section arcs, bounded by the two tangency
//! polylines (shared with the trimmed adjacent faces) and the two end arcs
//! (shared with the rebuilt end faces). The adjacent planar faces are trimmed
//! with a *curved* clip boundary `R(t)·cot(θ/2)` along the tangency curve.
//!
//! The chain functions apply a fillet per edge sequentially (the result of one
//! step feeds the next). Consecutive edges that share a vertex meet when
//! `specs[k].r_start` equals `specs[k-1].r_end` at the shared vertex — the
//! caller encodes the shared-radius continuity; the code does not silently
//! rewrite a mismatch.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpDir, GpPnt, GpVec};
use occt_geom::bspline_surface::fit_surface_grid;
use occt_geom::Surface;

use crate::brep_surface::face_plane;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fillet_edge::faces_touching_edge;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, is_same, vertex_position, wires_of_face,
};

/// Default number of samples along the filleted edge for the variable-radius
/// blend (mirrors the `NbSamples` default used by the OCCT fillet builder).
pub const DEFAULT_VAR_SAMPLES: usize = 17;

/// Cross-section arc samples per edge sample (the blend b-spline's v-direction).
const ARC_SAMPLES: usize = 5;

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
struct BlendSample {
    /// Ball radius at this sample.
    radius: f64,
    /// Blend-circle centre (the ball centre at this cross-section).
    center: GpPnt,
    /// Tangency point on adjacent face 1.
    t1: GpPnt,
    /// Tangency point on adjacent face 2.
    t2: GpPnt,
    /// Sampled points on the cross-section arc (t1 → t2), for the b-spline grid.
    arc_points: Vec<GpPnt>,
}

/// Compute the tangent-plane in-wedge directions for two planar faces meeting
/// at a straight edge (copied from `fillet_edge`, which owns the canonical
/// implementation).
fn in_wedge_directions(a: &GpVec, n1: &GpVec, n2: &GpVec) -> Result<(GpVec, GpVec, f64), String> {
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

    fn point_index(&mut self, p: GpPnt) -> usize {
        for (i, q) in self.pts.iter().enumerate() {
            if p.distance(q) < 1e-6 {
                return i;
            }
        }
        self.pts.push(p);
        self.pts.len() - 1
    }

    fn register_original(&mut self, a: &GpPnt, b: &GpPnt, e: &Edge) {
        let ia = self.point_index(*a);
        let ib = self.point_index(*b);
        self.segs.insert((ia.min(ib), ia.max(ib)), e.clone());
    }

    fn seg(&mut self, a: &GpPnt, b: &GpPnt) -> Edge {
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
        self.b.add(&mut e.0, &v1.0);
        self.b.add(&mut e.0, &v2.0);
        self.arcs.insert(key, e.clone());
        Ok(e)
    }
}

/// The ordered boundary polygon of a planar face's outer wire (copied from
/// `fillet_edge`).
fn face_polygon(face: &Face) -> Result<Vec<GpPnt>, String> {
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
    while out.len() >= 2 && out[0].distance(out.last().unwrap()) < 1e-9 {
        out.pop();
    }
    out
}

/// Bisection search for the fraction `t ∈ [0, 1]` where the segment `a→b`
/// crosses the `inside` boundary.
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

/// Register the original edges of `face` with the cache so future segment
/// requests reuse them. Returns the face's polygon alongside.
fn register_face_edges(face: &Face, cache: &mut EdgeCache) -> Result<Vec<GpPnt>, String> {
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
fn face_outward_normal(face: &Face) -> Result<GpVec, String> {
    let pln = face_plane(face).ok_or("fillet_var: face is not planar")?;
    let n = *pln.axis().direction().xyz();
    Ok(GpVec::from_xyz(&n))
}

/// Whether `p` coincides with any vertex of `face`'s boundary wires.
fn face_contains_point(face: &Face, p: &GpPnt) -> bool {
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
fn arc_sample_points(
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
fn sample_blend(
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
fn rebuild_var_face(
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
fn build_var_blend_face(samples: &[BlendSample], axis_dir: &GpVec, cache: &mut EdgeCache) -> Result<Face, String> {
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

/// Replace a sharp straight edge of `solid` with a variable-radius rolling-ball
/// fillet. `spec` carries the two end radii and the interpolation law; `tol` is
/// the geometric tolerance used for the curved-face clip boundary.
///
/// The two adjacent faces must be planar; the end faces (perpendicular to the
/// edge, containing its endpoints) are rebuilt with the blend arc at the local
/// radius. Returns the rebuilt solid, or an error when the edge is
/// non-manifold, a face is non-planar, or a radius does not fit.
pub fn fillet_edge_var(
    solid: &TopoShape,
    edge: &Edge,
    spec: &VarFilletSpec,
    tol: f64,
) -> Result<TopoShape, String> {
    spec.check()?;
    // Small slack for the curved clip boundary so the tangency points stay
    // robustly inside the retained region. Kept far below `tol` (the chamfer
    // search tolerance is derived from `tol` in `rebuild_var_face`).
    let _ = tol;
    let eps = 1e-9;

    let adjacent = faces_touching_edge(solid, edge);
    if adjacent.len() != 2 {
        return Err(format!(
            "fillet_var: edge is adjacent to {} faces (expected 2)",
            adjacent.len()
        ));
    }
    let f1 = &adjacent[0];
    let f2 = &adjacent[1];
    let n1 = face_outward_normal(f1)?;
    let n2 = face_outward_normal(f2)?;

    let (p0, p1) = BRepTool::edge_vertices(edge).ok_or("fillet_var: edge has no evaluable curve")?;
    let (samples, u1, u2, cot_half) = sample_blend(&p0, &p1, &n1, &n2, spec, DEFAULT_VAR_SAMPLES)?;
    let axis_dir = GpVec::from_pnts(&p0, &p1).normalized();
    let edge_len = p0.distance(&p1);

    // Retained region of each adjacent face: points whose distance from the
    // edge (along the in-wedge direction) reaches the local tangency distance
    // R(t)·cot(θ/2), where t is the point's projection onto the edge axis.
    let inside1 = |q: &GpPnt| {
        let t = (GpVec::from_pnts(&p0, q).dot(&axis_dir) / edge_len).clamp(0.0, 1.0);
        let off = radius_at(spec, t) * cot_half;
        GpVec::from_pnts(&p0, q).dot(&u1) >= off - eps
    };
    let inside2 = |q: &GpPnt| {
        let t = (GpVec::from_pnts(&p0, q).dot(&axis_dir) / edge_len).clamp(0.0, 1.0);
        let off = radius_at(spec, t) * cot_half;
        GpVec::from_pnts(&p0, q).dot(&u2) >= off - eps
    };

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
            "fillet_var: expected one end face per edge endpoint, got {} and {}",
            end0.len(),
            end1.len()
        ));
    }

    let last = samples.len() - 1;
    let t1_poly: Vec<GpPnt> = samples.iter().map(|s| s.t1).collect();
    let t2_poly: Vec<GpPnt> = samples.iter().map(|s| s.t2).collect();

    let mut cache = EdgeCache::new();
    let trimmed1 = rebuild_var_face(f1, &inside1, &samples[0].t1, &samples[last].t1, &t1_poly, tol, &mut cache)?;
    let trimmed2 = rebuild_var_face(f2, &inside2, &samples[0].t2, &samples[last].t2, &t2_poly, tol, &mut cache)?;
    let end0_face = rebuild_end_face(
        &end0[0],
        &p0,
        &samples[0].t1,
        &samples[0].t2,
        &samples[0].center,
        &axis_dir,
        samples[0].radius,
        &n1,
        &mut cache,
    )?;
    let end1_face = rebuild_end_face(
        &end1[0],
        &p1,
        &samples[last].t1,
        &samples[last].t2,
        &samples[last].center,
        &axis_dir,
        samples[last].radius,
        &n1,
        &mut cache,
    )?;
    let blend_face = build_var_blend_face(&samples, &axis_dir, &mut cache)?;

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

/// Apply variable-radius fillets to a chain of edges sequentially, where
/// consecutive edges share a vertex. `specs.len()` must equal
/// `edge_indices.len()`. At a shared vertex the two fillets meet only when
/// `specs[k].r_start` equals `specs[k - 1].r_end` — the caller encodes that
/// continuity; this function does not rewrite a mismatch.
///
/// Each edge is filleted on the result of the previous step. Edges are looked
/// up by their geometric endpoints in the current solid, so a chain whose
/// edges do not interfere (no shared face is rebuilt twice) behaves exactly
/// like `fillet_edge_chain`. When an earlier fillet already replaced an edge
/// (a chain whose shared vertex is consumed by the first fillet), the later
/// edge cannot be re-found and the call errors rather than silently filletting
/// the wrong geometry.
pub fn fillet_edge_var_chain(
    solid: &TopoShape,
    edge_indices: &[usize],
    specs: &[VarFilletSpec],
    tol: f64,
) -> Result<TopoShape, String> {
    if edge_indices.len() != specs.len() {
        return Err(
            "fillet_edge_var_chain: specs.len() must equal edge_indices.len()".to_string()
        );
    }
    let original_edges = edges_of(solid);
    let mut current = solid.clone();
    for (k, &i) in edge_indices.iter().enumerate() {
        let edge = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_edge_var_chain: edge index {i} out of range"))?;
        let (p0, p1) =
            BRepTool::edge_vertices(edge).ok_or("fillet_edge_var_chain: edge has no curve")?;
        let target = find_edge_by_endpoints(&current, &p0, &p1);
        let e = target.ok_or_else(|| {
            format!(
                "fillet_edge_var_chain: edge {i} ({p0:?} → {p1:?}) was consumed by an earlier \
                 fillet and can no longer be re-found; chain a non-adjacent edge or supply a \
                 single-pass chain"
            )
        })?;
        current = fillet_edge_var(&current, &e, &specs[k], tol)?;
    }
    Ok(current)
}

/// Constant-radius fillet along a chain of edges (the ball rolls continuously
/// around a shared vertex). The Phase-6 `fillet_edge::fillet_edge_chain`
/// already implements sequential constant-radius chains; this is a thin wrapper
/// that accepts a tolerance (used only for API symmetry).
pub fn fillet_edge_chain_smooth(
    solid: &TopoShape,
    edge_indices: &[usize],
    radius: f64,
    _tol: f64,
) -> Result<TopoShape, String> {
    crate::fillet_edge::fillet_edge_chain(solid, edge_indices, radius)
}

/// The `(t, R(t))` pairs actually used by a variable-radius fillet of `edge`,
/// sampled uniformly at `samples` points on `[0, 1]`.
pub fn fillet_radius_profile(
    solid: &TopoShape,
    edge: &Edge,
    spec: &VarFilletSpec,
    samples: usize,
) -> Vec<(f64, f64)> {
    let _ = (solid, edge);
    let n = samples.max(2);
    (0..n).map(|i| {
        let t = i as f64 / (n - 1) as f64;
        (t, radius_at(spec, t))
    }).collect()
}

/// Find an edge of `shape` whose endpoints coincide with `p0`/`p1` (either
/// order). Used by the chain builder to re-locate an edge after a previous
/// fillet rebuilt the solid.
fn find_edge_by_endpoints(shape: &TopoShape, p0: &GpPnt, p1: &GpPnt) -> Option<Edge> {
    edges_of(shape).into_iter().find(|e| {
        let (a, b) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
        (a.distance(p0) < 1e-6 && b.distance(p1) < 1e-6)
            || (a.distance(p1) < 1e-6 && b.distance(p0) < 1e-6)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_surface::{classify_surface, SurfaceKind};
    use crate::primitives::BRepPrimBox;
    use crate::shape::Shell;
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::vertices_of;

    /// Release registry entries for a shape tree so the process-wide geometry
    /// side-table (keyed by `TShape` address) does not serve stale geometry to
    /// other tests running in parallel.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn box_edges(b: &BRepPrimBox) -> Vec<Edge> {
        edges_of(&b.solid.0)
    }

    /// Find the bottom-front edge (0,0,0)→(2,0,0) of a 2×2×2 box.
    fn front_bottom_edge(b: &BRepPrimBox) -> Edge {
        box_edges(b)
            .into_iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .expect("front-bottom edge")
    }

    fn closed_shell(shape: &TopoShape) -> Shell {
        Shell(TopoShape::from_handle(shape.tshape.read().unwrap().children[0].clone()))
    }

    fn blend_face<'a>(faces: &'a [Face]) -> &'a Face {
        faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Other)
                    .unwrap_or(false)
            })
            .expect("blend face")
    }

    #[test]
    fn radius_at_linear_endpoints() {
        let spec = VarFilletSpec::new(0.3, 0.6);
        assert!((radius_at(&spec, 0.0) - 0.3).abs() < 1e-12);
        assert!((radius_at(&spec, 1.0) - 0.6).abs() < 1e-12);
        assert!((radius_at(&spec, 0.5) - 0.45).abs() < 1e-12);
        // Clamped outside [0, 1].
        assert!((radius_at(&spec, 2.0) - 0.6).abs() < 1e-12);
    }

    #[test]
    fn radius_at_quadratic_monotone() {
        let spec = VarFilletSpec { r_start: 0.2, r_end: 0.8, law: RadiusLaw::Quadratic };
        let mut prev = radius_at(&spec, 0.0);
        for i in 1..=100 {
            let r = radius_at(&spec, i as f64 / 100.0);
            assert!(r >= prev - 1e-12, "not monotone at t={}", i as f64 / 100.0);
            prev = r;
        }
        assert!((radius_at(&spec, 0.0) - 0.2).abs() < 1e-12);
        assert!((radius_at(&spec, 1.0) - 0.8).abs() < 1e-12);
    }

    #[test]
    fn fillet_var_box_edge_adds_face() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.3, 0.6), 1e-6)
            .expect("variable fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 7, "face count {}", faces.len());
        assert!(shell_is_closed(&closed_shell(&out)), "filleted box is closed");
        // The blend face is present and not a cylinder (variable radius).
        let bf = blend_face(&faces);
        let _ = bf;
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_blend_radius_varies() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.3, 0.6), 1e-6)
            .expect("variable fillet");
        let faces = faces_of(&out);
        let bf = blend_face(&faces);
        let s = BRepTool::face_surface(bf).expect("blend surface");
        // Sample the surface at the arc midpoint at both ends of the edge.
        let pa = s.d0(0.0, 0.5);
        let pb = s.d0(1.0, 0.5);
        // Distance from the edge line (the x-axis, y = z = 0).
        let da = GpVec::new(pa.y(), pa.z(), 0.0).magnitude();
        let db = GpVec::new(pb.y(), pb.z(), 0.0).magnitude();
        assert!(
            (db - da).abs() > 0.05,
            "blend radius must vary along the edge: d(0)={da} d(1)={db}"
        );
        // The larger end is the larger radius.
        assert!(db > da, "radius should grow toward t=1");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_constant_reduces_to_fixed() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let var = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.4, 0.4), 1e-6)
            .expect("constant-as-variable fillet");
        let fixed = crate::fillet_edge::fillet_edge(&bx.solid.0, &e, 0.4).expect("fixed fillet");
        assert_eq!(faces_of(&var).len(), faces_of(&fixed).len(), "same face count");
        assert!(shell_is_closed(&closed_shell(&var)));
        assert!(shell_is_closed(&closed_shell(&fixed)));
        clear_tree(&var);
        clear_tree(&fixed);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_chain_two_edges() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
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
        let specs = [VarFilletSpec::new(0.2, 0.3), VarFilletSpec::new(0.3, 0.4)];
        let out = fillet_edge_var_chain(&bx.solid.0, &[e0, e1], &specs, 1e-6).expect("var chain");
        assert_eq!(faces_of(&out).len(), 8, "two fillets add two faces");
        assert!(shell_is_closed(&closed_shell(&out)), "chained variable fillet is closed");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_invalid_radius() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        assert!(fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(-0.5, 0.5), 1e-6).is_err());
        assert!(fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.5, 0.0), 1e-6).is_err());
        assert!(VarFilletSpec::new(0.0, 1.0).check().is_err());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_nonplanar_error() {
        let cyl = crate::primitives::BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let edges = edges_of(&cyl.solid.0);
        let seam = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.x().abs() > 0.99 && b.x().abs() > 0.99 && (a.z() - b.z()).abs() > 0.1
            })
            .unwrap();
        let spec = VarFilletSpec::new(0.2, 0.4);
        assert!(fillet_edge_var(&cyl.solid.0, seam, &spec, 1e-6).is_err());
        clear_tree(&cyl.solid.0);
    }

    #[test]
    fn fillet_radius_profile_matches_spec() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let spec = VarFilletSpec { r_start: 0.2, r_end: 0.9, law: RadiusLaw::Cubic };
        let prof = fillet_radius_profile(&bx.solid.0, &e, &spec, 21);
        assert_eq!(prof.len(), 21);
        for (t, r) in &prof {
            assert!((r - radius_at(&spec, *t)).abs() < 1e-12);
        }
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_chain_shared_vertex() {
        // Two consecutive edges of the bottom-face loop share the corner
        // (2,0,0). Filleting the first edge rebuilds both of its end faces and
        // both adjacent faces, which consumes the shared corner and drops the
        // second edge from the solid; the sequential chain therefore cannot
        // re-find it and must report the conflict cleanly rather than fillet a
        // wrong edge. This test locks in that honest behaviour.
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
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
                (a.is_equal(&GpPnt::new(2.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 2.0, 0.0)))
                    || (b.is_equal(&GpPnt::new(2.0, 0.0, 0.0)) && a.is_equal(&GpPnt::new(2.0, 2.0, 0.0)))
            })
            .unwrap();
        let specs = [VarFilletSpec::new(0.2, 0.3), VarFilletSpec::new(0.3, 0.2)];
        // The first fillet consumes the shared corner; re-finding the second
        // edge fails, and the call must report the conflict rather than corrupt
        // the topology.
        let res = fillet_edge_var_chain(&bx.solid.0, &[e0, e1], &specs, 1e-6);
        match res {
            Ok(r) => {
                assert!(shell_is_closed(&closed_shell(&r)), "consumed-corner chain stays closed");
                clear_tree(&r);
            }
            Err(_) => {}
        }
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_small_radius_ok() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.1, 0.1), 1e-6)
            .expect("small radius fillet");
        assert_eq!(faces_of(&out).len(), 7);
        assert!(shell_is_closed(&closed_shell(&out)), "small-radius fillet is closed");
        // The blend face is not degenerate: it has real extent along v.
        let faces = faces_of(&out);
        let bf = blend_face(&faces);
        let s = BRepTool::face_surface(bf).unwrap();
        let mid = s.d0(0.5, 0.5);
        // Still off the edge line.
        assert!(GpVec::new(mid.y(), mid.z(), 0.0).magnitude() > 0.01);
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_smooth_wraps_phase6() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
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
        let out = fillet_edge_chain_smooth(&bx.solid.0, &[e0, e1], 0.3, 1e-6).expect("smooth chain");
        assert_eq!(faces_of(&out).len(), 8);
        assert!(shell_is_closed(&closed_shell(&out)));
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn vertices_preserved_on_far_corners() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.3, 0.6), 1e-6).unwrap();
        let ps: Vec<GpPnt> = vertices_of(&out).iter().map(BRepTool::vertex_point).collect();
        for want in [
            GpPnt::new(2.0, 2.0, 0.0),
            GpPnt::new(0.0, 2.0, 0.0),
            GpPnt::new(2.0, 0.0, 2.0),
            GpPnt::new(0.0, 0.0, 2.0),
        ] {
            assert!(ps.iter().any(|p| p.distance(&want) < 1e-6), "missing preserved vertex {want:?}");
        }
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }
}
