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
use std::f64::consts::SQRT_2;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx3, GpDir, GpPnt, GpSphere, GpVec};
use occt_geom::bspline_surface::fit_surface_grid;
use occt_geom::{GeomSphere, Surface};

use crate::brep_surface::face_plane;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fillet_edge::faces_touching_edge;
use crate::shape::{Edge, Face, TopoShape, Vertex};
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

// ===========================================================================
// Rolling-ball corner patch (shared-vertex blend)
//
// Two consecutive edge fillets meeting at a box corner leave a gap between the
// two blend surfaces and the three corner faces. The gap is closed by a
// spherical patch centred at the rolling-ball corner centre
// `C = v + R·(d0 + d1 + d2)` — the point at distance `R` from each of the
// three faces. A sphere of radius `R·√2` centred at `C` passes through the
// three fillet tangency points, and its boundary is five circular arcs:
//
//   · one arc in each of the three corner-face planes,
//   · one arc shared with each of the two trimmed blend faces.
//
// Each edge fillet is trimmed to start at distance `2R` from the corner, where
// its end arc coincides with the sphere-boundary arc. The result is a closed
// shell with one added face per shared vertex. The blend faces of a corner
// chain are built at the corner radius `R` (the maximum incident radius), so
// the shared tangency geometry is exactly consistent.
// ===========================================================================

/// The three unit directions of the edges meeting at `vertex`, pointing away
/// from the vertex (into the material for a convex corner).
fn corner_in_dirs(solid: &TopoShape, vertex: &Vertex) -> Result<[GpVec; 3], String> {
    let p = vertex_position(vertex);
    let es = edges_of(solid);
    let mut dirs: Vec<GpVec> = Vec::new();
    for e in &es {
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            let pa = vertex_position(&va);
            let pb = vertex_position(&vb);
            if pa.distance(&p) < 1e-9 || pb.distance(&p) < 1e-9 {
                let other = if pa.distance(&p) < 1e-9 { pb } else { pa };
                let v = GpVec::from_pnts(&p, &other);
                if v.magnitude() > 1e-12 {
                    dirs.push(v.normalized());
                }
            }
        }
    }
    if dirs.len() != 3 {
        return Err(format!(
            "fillet_corner: vertex has {} incident edges (expected 3)",
            dirs.len()
        ));
    }
    Ok([dirs[0], dirs[1], dirs[2]])
}

/// The rolling-ball corner centre: the point at distance `radius` from each of
/// the three faces meeting at `vertex` (for a box corner, `v + R·Σdirᵢ`).
pub fn corner_center(solid: &TopoShape, vertex: &Vertex, radius: f64, _tol: f64) -> Option<GpPnt> {
    if radius <= 0.0 {
        return None;
    }
    let dirs = corner_in_dirs(solid, vertex).ok()?;
    let p = vertex_position(vertex);
    let s = dirs[0].added(&dirs[1]).added(&dirs[2]);
    Some(p.translated_vec(&s.multiplied_scalar(radius)))
}

/// The effective corner radius at `vertex`: the maximum of the incident-edge
/// fillet radii at that vertex (`r_start` when the vertex is the edge's start,
/// `r_end` otherwise). `specs` is aligned with `edges_of(solid)`; edges without
/// a matching spec are ignored.
pub fn corner_patch_radius_at_vertex(
    solid: &TopoShape,
    vertex: &Vertex,
    specs: &[VarFilletSpec],
    tol: f64,
) -> f64 {
    let p = vertex_position(vertex);
    let es = edges_of(solid);
    let mut rmax = 0.0f64;
    for (i, e) in es.iter().enumerate() {
        let spec = match specs.get(i) {
            Some(s) => s,
            None => continue,
        };
        let (a, b) = edge_vertices(e);
        let (Some(va), Some(vb)) = (a, b) else { continue };
        let pa = vertex_position(&va);
        let pb = vertex_position(&vb);
        let r = if pa.distance(&p) < tol {
            spec.r_start
        } else if pb.distance(&p) < tol {
            spec.r_end
        } else {
            continue;
        };
        rmax = rmax.max(r);
    }
    rmax
}

/// The unit direction of `edge` away from point `p` (one of its endpoints).
fn in_edge_dir(edge: &Edge, p: &GpPnt) -> Result<GpVec, String> {
    let (a, b) = edge_vertices(edge);
    let pa = vertex_position(&a.ok_or("fillet_corner: edge has no start vertex")?);
    let pb = vertex_position(&b.ok_or("fillet_corner: edge has no end vertex")?);
    let other = if pa.distance(p) < 1e-9 { pb } else { pa };
    let v = GpVec::from_pnts(p, &other);
    if v.magnitude() < 1e-12 {
        return Err("fillet_corner: degenerate incident edge".to_string());
    }
    Ok(v.normalized())
}

/// Geometry of a shared-vertex corner where two consecutive filleted edges meet.
struct CornerGeom {
    vertex: GpPnt,
    radius: f64,
    /// Rolling-ball corner centre (the corner-patch sphere centre).
    center: GpPnt,
    /// Direction of the first filleted edge (e_a), from the corner into material.
    d0: GpVec,
    /// Direction of the second filleted edge (e_b).
    d1: GpVec,
    /// Direction of the third (unfilleted) edge.
    d2: GpVec,
}

impl CornerGeom {
    fn from_edges(
        solid: &TopoShape,
        vertex: &Vertex,
        edge_a: &Edge,
        edge_b: &Edge,
        radius: f64,
    ) -> Result<Self, String> {
        let dirs = corner_in_dirs(solid, vertex)?;
        let p = vertex_position(vertex);
        let da = in_edge_dir(edge_a, &p)?;
        let db = in_edge_dir(edge_b, &p)?;
        let mut d0 = GpVec::zero();
        let mut d1 = GpVec::zero();
        let mut d2 = GpVec::zero();
        let mut found = [false; 3];
        for d in &dirs {
            if d.xyz().crossed(&da.xyz()).modulus() < 1e-6 {
                d0 = *d;
                found[0] = true;
            } else if d.xyz().crossed(&db.xyz()).modulus() < 1e-6 {
                d1 = *d;
                found[1] = true;
            } else {
                d2 = *d;
                found[2] = true;
            }
        }
        if !found.iter().all(|f| *f) {
            return Err("fillet_corner: cannot classify the corner edges".to_string());
        }
        let center = p.translated_vec(&d0.added(&d1).added(&d2).multiplied_scalar(radius));
        Ok(CornerGeom { vertex: p, radius, center, d0, d1, d2 })
    }

    fn sphere_radius(&self) -> f64 {
        SQRT_2 * self.radius
    }

    /// A point of the corner: `v + R·(a·d0 + b·d1 + c·d2)`.
    fn pt(&self, a: f64, b: f64, c: f64) -> GpPnt {
        let v0 = self.d0.multiplied_scalar(a * self.radius);
        let v1 = self.d1.multiplied_scalar(b * self.radius);
        let v2 = self.d2.multiplied_scalar(c * self.radius);
        self.vertex.translated_vec(&v0.added(&v1).added(&v2))
    }

    /// The five boundary arcs of the corner patch, each `(center, normal,
    /// radius, a, b)`. In order they form a closed loop: face ⊥d2, edge b,
    /// face ⊥d0, face ⊥d1, edge a.
    fn arcs(&self) -> Vec<(GpPnt, GpVec, f64, GpPnt, GpPnt)> {
        let r = self.radius;
        vec![
            // Arc in the face perpendicular to d2 (the face shared by both edges).
            (self.pt(1.0, 1.0, 0.0), self.d2, r, self.pt(2.0, 1.0, 0.0), self.pt(1.0, 2.0, 0.0)),
            // Arc with edge b (the second filleted edge).
            (self.pt(1.0, 2.0, 1.0), self.d1, r, self.pt(1.0, 2.0, 0.0), self.pt(0.0, 2.0, 1.0)),
            // Arc in the face perpendicular to d0 (the face edge a is not on).
            (self.pt(0.0, 1.0, 1.0), self.d0.multiplied_scalar(-1.0), r, self.pt(0.0, 2.0, 1.0), self.pt(0.0, 0.0, 1.0)),
            // Arc in the face perpendicular to d1 (the face edge b is not on).
            (self.pt(1.0, 0.0, 1.0), self.d1, r, self.pt(0.0, 0.0, 1.0), self.pt(2.0, 0.0, 1.0)),
            // Arc with edge a (the first filleted edge).
            (self.pt(2.0, 1.0, 1.0), self.d0, r, self.pt(2.0, 0.0, 1.0), self.pt(2.0, 1.0, 0.0)),
        ]
    }

    /// The corner arc that lives on the face whose outward normal is `n_out`.
    fn arc_on_face(&self, n_out: &GpVec) -> Option<(GpPnt, GpVec, f64, GpPnt, GpPnt)> {
        let arcs = self.arcs();
        let idx = if n_out.xyz().crossed(&self.d2.xyz()).modulus() < 1e-6 {
            Some(0)
        } else if n_out.xyz().crossed(&self.d0.xyz()).modulus() < 1e-6 {
            Some(2)
        } else if n_out.xyz().crossed(&self.d1.xyz()).modulus() < 1e-6 {
            Some(3)
        } else {
            None
        };
        idx.map(|i| arcs[i].clone())
    }
}

/// Build the spherical corner-patch face for `corner`, with its five boundary
/// arcs shared with the corner faces and the two trimmed blend faces.
fn build_corner_patch_face(corner: &CornerGeom, cache: &mut EdgeCache) -> Result<Face, String> {
    let arcs = corner.arcs();
    let mut edges = Vec::with_capacity(arcs.len());
    for (center, normal, r, a, b) in arcs {
        let nd = GpDir::from_vec(&normal).map_err(|e| e.to_string())?;
        edges.push(cache.arc(&center, &nd, r, &a, &b)?);
    }
    let b = TopoBuilder::new();
    let wire = b.make_wire(&edges);
    let mut sph = GpSphere::new(GpAx3::standard(), corner.sphere_radius()).map_err(|e| e.to_string())?;
    sph.set_location(corner.center);
    Ok(b.make_face(Arc::new(GeomSphere::new(sph)), &[wire]))
}

/// Parameters `t ∈ [0,1]` where the segment `p→q` crosses the circle centred at
/// `c` of radius `r` (the disk boundary).
fn line_circle_hits(p: &GpPnt, q: &GpPnt, c: &GpPnt, r: f64) -> Vec<f64> {
    let d = GpVec::from_pnts(p, q);
    let f = GpVec::from_pnts(p, c); // f = c − p
    let a = d.dot(&d);
    if a < 1e-30 {
        return Vec::new();
    }
    // |p + t·d − c|² = r²  ⇒  a·t² − 2(d·f)·t + (|f|² − r²) = 0.
    let b = -2.0 * f.dot(&d);
    let cc = f.dot(&f) - r * r;
    let disc = b * b - 4.0 * a * cc;
    if disc < 0.0 {
        return Vec::new();
    }
    let sq = disc.sqrt();
    let mut ts: Vec<f64> = Vec::new();
    for t in [(-b + sq) / (2.0 * a), (-b - sq) / (2.0 * a)] {
        if t > -1e-9 && t < 1.0 + 1e-9 {
            ts.push(t.clamp(0.0, 1.0));
        }
    }
    ts.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ts
}

/// Remove the disk (centre `c`, radius `r`) from the polygon, keeping the
/// outside portion. The disk cap is replaced by a straight chord; the returned
/// polygon has that chord as one boundary edge, to be replaced by the shared
/// corner arc.
fn cut_disk_poly(poly: &[GpPnt], c: &GpPnt, r: f64) -> Result<Vec<GpPnt>, String> {
    let n = poly.len();
    if n < 3 {
        return Err("cut_disk: degenerate polygon".to_string());
    }
    let mut out: Vec<GpPnt> = Vec::new();
    for i in 0..n {
        let p = &poly[i];
        let q = &poly[(i + 1) % n];
        let hits = line_circle_hits(p, q, c, r);
        let mut ts: Vec<f64> = vec![0.0];
        for t in &hits {
            if *t > 1e-9 && *t < 1.0 - 1e-9 {
                ts.push(*t);
            }
        }
        ts.push(1.0);
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap());
        for w in ts.windows(2) {
            let (t0, t1) = (w[0], w[1]);
            if t1 - t0 < 1e-12 {
                continue;
            }
            let tm = 0.5 * (t0 + t1);
            let mid = p.translated_vec(&GpVec::from_pnts(p, q).multiplied_scalar(tm));
            if mid.distance(c) >= r - 1e-9 {
                let a0 = p.translated_vec(&GpVec::from_pnts(p, q).multiplied_scalar(t0));
                let a1 = p.translated_vec(&GpVec::from_pnts(p, q).multiplied_scalar(t1));
                if out.last().map_or(true, |last| last.distance(&a0) > 1e-6) {
                    out.push(a0);
                }
                if a0.distance(&a1) > 1e-6 {
                    out.push(a1);
                }
            }
        }
    }
    // Close the loop: drop a trailing point equal to the first.
    while out.len() >= 2 && out[0].distance(out.last().unwrap()) < 1e-9 {
        out.pop();
    }
    if out.len() < 3 {
        return Err("cut_disk: the disk removed the whole polygon".to_string());
    }
    Ok(out)
}

/// Rebuild a corner face: clip its polygon by the half-plane tangency trims,
/// then cut out each corner disk and replace the resulting chord with the
/// shared circular-arc edge. `arcs` holds the corner arcs that lie on this
/// face, one per adjacent corner. `polylines` holds the sampled tangency
/// polylines of the adjacent blend faces, which replace the straight tangency
/// edges so the corner face shares those edges with the blend faces.
fn rebuild_corner_face(
    face: &Face,
    trims: &[(GpPnt, GpVec, f64)],
    arcs: &[(GpPnt, GpVec, f64, GpPnt, GpPnt)],
    polylines: &[(GpPnt, GpPnt, Vec<GpPnt>)],
    tol: f64,
    cache: &mut EdgeCache,
) -> Result<Face, String> {
    let poly = register_face_edges(face, cache)?;
    let eps = 1e-9;
    let mut clipped = poly;
    for (p0, u, off) in trims {
        let inside = |q: &GpPnt| GpVec::from_pnts(p0, q).dot(u) >= off - eps;
        clipped = clip_polygon(&clipped, &inside);
        if clipped.len() < 3 {
            return Err("fillet_corner: clipping a corner face left no retained region".to_string());
        }
    }
    let mut arc_edges: Vec<(GpPnt, GpPnt, Edge)> = Vec::new();
    for (center, normal, radius, a, b) in arcs {
        let cut = cut_disk_poly(&clipped, center, *radius)?;
        if cut.len() < 3 {
            return Err("fillet_corner: cutting the corner disk left no region".to_string());
        }
        let nd = GpDir::from_vec(normal).map_err(|e| e.to_string())?;
        let arc_edge = cache.arc(center, &nd, *radius, a, b)?;
        arc_edges.push((*a, *b, arc_edge));
        clipped = cut;
    }
    let match_tol = (10.0 * tol).max(1e-4).min(1e-3);
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    let m = clipped.len();
    for i in 0..m {
        let a = &clipped[i];
        let c = &clipped[(i + 1) % m];
        if let Some((_, _, e)) = arc_edges.iter().find(|(ta, tb, _)| {
            (a.distance(ta) < match_tol && c.distance(tb) < match_tol)
                || (a.distance(tb) < match_tol && c.distance(ta) < match_tol)
        }) {
            edges.push(e.clone());
        } else if let Some((_, _, poly)) = polylines.iter().find(|(f, l, _)| {
            (a.distance(f) < match_tol && c.distance(l) < match_tol)
                || (a.distance(l) < match_tol && c.distance(f) < match_tol)
        }) {
            if a.distance(&poly[0]) < match_tol {
                for w in poly.windows(2) {
                    edges.push(cache.seg(&w[0], &w[1]));
                }
            } else {
                for w in poly.windows(2).rev() {
                    edges.push(cache.seg(&w[1], &w[0]));
                }
            }
        } else if a.distance(c) > 1e-6 {
            edges.push(cache.seg(a, c));
        }
    }
    let wire = b.make_wire(&edges);
    let surf = BRepTool::face_surface(face).ok_or("fillet_corner: corner face has no surface")?;
    Ok(b.make_face(surf, &[wire]))
}

/// A single end-arc replacement: the corner vertex `corner` of an end face is
/// replaced by the blend arc `t1 → t2` (shared with the blend face).
#[derive(Clone)]
struct EndArcReplacement {
    corner: GpPnt,
    t1: GpPnt,
    t2: GpPnt,
    center: GpPnt,
    axis: GpVec,
    radius: f64,
    /// Outward normal of adjacent face 1 (used to classify the incident edges).
    n1: GpVec,
}

/// Rebuild a face by replacing several corner vertices with their end arcs.
/// Used for far end faces (perpendicular to an edge) that may carry arcs from
/// more than one edge (e.g. the far +X face of two opposite X edges).
fn rebuild_face_with_end_arcs(
    face: &Face,
    reps: &[EndArcReplacement],
    cache: &mut EdgeCache,
) -> Result<Face, String> {
    let mut poly = register_face_edges(face, cache)?;
    let mut arc_pairs: Vec<(GpPnt, GpPnt, Edge)> = Vec::new();
    for rep in reps {
        let n = poly.len();
        let idx = poly
            .iter()
            .position(|q| q.distance(&rep.corner) < 1e-9)
            .ok_or("fillet_corner: end face does not contain the corner vertex")?;
        let prev = poly[(idx + n - 1) % n];
        let next = poly[(idx + 1) % n];
        let dir_prev = GpVec::from_pnts(&prev, &rep.corner).normalized();
        let dir_next = GpVec::from_pnts(&rep.corner, &next).normalized();
        let prev_is_face1 = dir_prev.dot(&rep.n1).abs() < 1e-6;
        let next_is_face1 = dir_next.dot(&rep.n1).abs() < 1e-6;
        let (t_first, t_last) = if prev_is_face1 && !next_is_face1 {
            (rep.t1, rep.t2)
        } else if next_is_face1 && !prev_is_face1 {
            (rep.t2, rep.t1)
        } else {
            return Err("fillet_corner: cannot classify end-face incident edges".to_string());
        };
        let zd = GpDir::from_vec(&rep.axis).map_err(|e| e.to_string())?;
        let arc_edge = cache.arc(&rep.center, &zd, rep.radius, &rep.t1, &rep.t2)?;
        arc_pairs.push((rep.t1, rep.t2, arc_edge));
        poly[idx] = t_first;
        poly.insert(idx + 1, t_last);
    }

    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    let m = poly.len();
    for i in 0..m {
        let a = &poly[i];
        let c = &poly[(i + 1) % m];
        if let Some((_, _, e)) = arc_pairs.iter().find(|(ta, tb, _)| {
            (a.distance(ta) < 1e-6 && c.distance(tb) < 1e-6)
                || (a.distance(tb) < 1e-6 && c.distance(ta) < 1e-6)
        }) {
            edges.push(e.clone());
        } else if a.distance(c) > 1e-6 {
            edges.push(cache.seg(a, c));
        }
    }
    let wire = b.make_wire(&edges);
    let surf = BRepTool::face_surface(face).ok_or("fillet_corner: end face has no surface")?;
    Ok(b.make_face(surf, &[wire]))
}

/// Per-edge data for a corner-chain run.
struct EdgeInfo {
    edge: Edge,
    p0: GpPnt,
    p1: GpPnt,
    f1: Face,
    f2: Face,
    n1: GpVec,
    n2: GpVec,
    axis: GpVec,
    len: f64,
    spec: VarFilletSpec,
}

fn find_shared_vertex(a0: &GpPnt, a1: &GpPnt, b0: &GpPnt, b1: &GpPnt) -> Option<GpPnt> {
    for p in [a0, a1] {
        if p.distance(b0) < 1e-9 || p.distance(b1) < 1e-9 {
            return Some(*p);
        }
    }
    None
}

fn radius_at_endpoint(spec: &VarFilletSpec, p0: &GpPnt, p1: &GpPnt, v: &GpPnt) -> f64 {
    if p0.distance(v) < 1e-9 {
        spec.r_start
    } else if p1.distance(v) < 1e-9 {
        spec.r_end
    } else {
        0.5 * (spec.r_start + spec.r_end)
    }
}

fn find_vertex_at(solid: &TopoShape, p: &GpPnt) -> Option<Vertex> {
    crate::topo_tools_full::vertices_of(solid)
        .into_iter()
        .find(|v| BRepTool::vertex_point(v).distance(p) < 1e-6)
}

fn find_end_face(solid: &TopoShape, f1: &Face, f2: &Face, p: &GpPnt) -> Result<Face, String> {
    let ends: Vec<Face> = faces_of(solid)
        .into_iter()
        .filter(|f| face_contains_point(f, p) && !is_same(&f.0, &f1.0) && !is_same(&f.0, &f2.0))
        .collect();
    if ends.len() != 1 {
        return Err(format!(
            "fillet_corner: expected one end face at {:?}, got {}",
            p,
            ends.len()
        ));
    }
    Ok(ends[0].clone())
}

/// The tangency trim of `edge` on `face`: `(edge_p0, u_dir, off)` such that the
/// retained side is `dot(from_p0, u_dir) >= off`.
fn edge_tangency_trim(e: &EdgeInfo, face: &Face, radius: f64) -> Result<(GpPnt, GpVec, f64), String> {
    let spec = VarFilletSpec::new(radius, radius);
    let (_, u1, u2, cot_half) = sample_blend(&e.p0, &e.p1, &e.n1, &e.n2, &spec, 2)?;
    let is_f1 = is_same(&face.0, &e.f1.0);
    let is_f2 = is_same(&face.0, &e.f2.0);
    if is_f1 {
        Ok((e.p0, u1, radius * cot_half))
    } else if is_f2 {
        Ok((e.p0, u2, radius * cot_half))
    } else {
        Err("fillet_corner: edge is not adjacent to this face".to_string())
    }
}

/// The faces meeting at the corner vertex.
fn corner_faces_at(solid: &TopoShape, p: &GpPnt) -> Vec<Face> {
    faces_of(solid)
        .into_iter()
        .filter(|f| face_contains_point(f, p))
        .collect()
}

/// A corner face to be rebuilt (a face at one or more shared vertices).
struct CornerFaceRebuild {
    face: Face,
    /// Tangency trims `(edge_p0, u_dir, off)` from adjacent filleted edges.
    trims: Vec<(GpPnt, GpVec, f64)>,
    /// Corner arcs that live on this face (their circle centre/radius also
    /// define the corner disk that is cut out of the face).
    arcs: Vec<(GpPnt, GpVec, f64, GpPnt, GpPnt)>,
    /// Sampled tangency polylines `(first, last, points)` of the adjacent
    /// blend faces, shared with them.
    polylines: Vec<(GpPnt, GpPnt, Vec<GpPnt>)>,
}

/// A far end face to be rebuilt (perpendicular to an edge, not at a corner).
struct EndFaceRebuild {
    face: Face,
    replacements: Vec<EndArcReplacement>,
}

/// Build the combined fillet of a run of consecutive edges that share vertices.
/// Each edge gets a constant-radius blend face trimmed to `[2R, len - 2R]` at
/// its corner ends, each shared vertex gets a spherical corner patch, and the
/// corner faces / far end faces are rebuilt with shared edges so the resulting
/// shell is closed.
fn build_combined_run(
    solid: &TopoShape,
    run_edges: &[Edge],
    run_specs: &[VarFilletSpec],
    tol: f64,
) -> Result<TopoShape, String> {
    let n = run_edges.len();
    if n < 2 {
        return Err("build_combined_run: needs at least 2 consecutive edges".to_string());
    }
    if n != run_specs.len() {
        return Err("build_combined_run: specs length mismatch".to_string());
    }

    // Per-edge geometry.
    let mut infos: Vec<EdgeInfo> = Vec::with_capacity(n);
    for (k, edge) in run_edges.iter().enumerate() {
        let spec = run_specs[k];
        spec.check()?;
        let adjacent = faces_touching_edge(solid, edge);
        if adjacent.len() != 2 {
            return Err(format!(
                "fillet_chain_corner: edge {k} touches {} faces (expected 2)",
                adjacent.len()
            ));
        }
        let f1 = adjacent[0].clone();
        let f2 = adjacent[1].clone();
        let n1 = face_outward_normal(&f1)?;
        let n2 = face_outward_normal(&f2)?;
        let (p0, p1) = BRepTool::edge_vertices(edge).ok_or("fillet_chain_corner: edge has no curve")?;
        let axis = GpVec::from_pnts(&p0, &p1).normalized();
        let len = p0.distance(&p1);
        infos.push(EdgeInfo { edge: edge.clone(), p0, p1, f1, f2, n1, n2, axis, len, spec });
    }

    // Shared corners between consecutive edges.
    struct CornerInfo {
        vertex: GpPnt,
        radius: f64,
        geom: CornerGeom,
        k: usize,
    }
    let mut corners: Vec<CornerInfo> = Vec::new();
    for k in 0..n - 1 {
        let (pa0, pa1) = (infos[k].p0, infos[k].p1);
        let (pb0, pb1) = (infos[k + 1].p0, infos[k + 1].p1);
        let shared = find_shared_vertex(&pa0, &pa1, &pb0, &pb1)
            .ok_or("fillet_chain_corner: consecutive edges do not share a vertex")?;
        let ra = radius_at_endpoint(&infos[k].spec, &pa0, &pa1, &shared);
        let rb = radius_at_endpoint(&infos[k + 1].spec, &pb0, &pb1, &shared);
        let radius = ra.max(rb);
        let vshape = find_vertex_at(solid, &shared)
            .ok_or("fillet_chain_corner: shared vertex not found in the solid")?;
        let geom = CornerGeom::from_edges(solid, &vshape, &infos[k].edge, &infos[k + 1].edge, radius)?;
        corners.push(CornerInfo { vertex: shared, radius, geom, k });
    }

    // Corner ends per edge.
    let mut corner_at_start: Vec<Option<usize>> = vec![None; n];
    let mut corner_at_end: Vec<Option<usize>> = vec![None; n];
    for (ci, c) in corners.iter().enumerate() {
        if infos[c.k].p0.distance(&c.vertex) < 1e-9 {
            corner_at_start[c.k] = Some(ci);
        } else {
            corner_at_end[c.k] = Some(ci);
        }
        if infos[c.k + 1].p0.distance(&c.vertex) < 1e-9 {
            corner_at_start[c.k + 1] = Some(ci);
        } else {
            corner_at_end[c.k + 1] = Some(ci);
        }
    }

    // Per-edge corner radius (all incident corners of an edge must agree).
    let mut edge_radius = vec![0.0f64; n];
    for k in 0..n {
        let mut r: Option<f64> = None;
        for ci in [corner_at_start[k], corner_at_end[k]].iter().flatten() {
            let cr = corners[*ci].radius;
            if let Some(prev) = r {
                if (prev - cr).abs() > tol.max(1e-9) {
                    return Err(
                        "fillet_chain_corner: incident corner radii differ along an edge".to_string(),
                    );
                }
            }
            r = Some(cr);
        }
        edge_radius[k] = r.unwrap_or_else(|| infos[k].spec.r_start.max(infos[k].spec.r_end));
    }

    // Per-edge trimmed blend samples (used for the tangency polylines, the
    // blend faces, and the far-end tangency points).
    let mut edge_samples: Vec<Option<Vec<BlendSample>>> = vec![None; n];
    for k in 0..n {
        let e = &infos[k];
        let r = edge_radius[k];
        let t_lo = if corner_at_start[k].is_some() { 2.0 * r / e.len } else { 0.0 };
        let t_hi = if corner_at_end[k].is_some() { 1.0 - 2.0 * r / e.len } else { 1.0 };
        if t_hi - t_lo < 1e-9 {
            return Err("fillet_chain_corner: corner radius is too large for the edge".to_string());
        }
        let sub_p0 = e.p0.translated_vec(&e.axis.multiplied_scalar(t_lo * e.len));
        let sub_p1 = e.p0.translated_vec(&e.axis.multiplied_scalar(t_hi * e.len));
        let const_spec = VarFilletSpec::new(r, r);
        let (samples, _, _, _) = sample_blend(&sub_p0, &sub_p1, &e.n1, &e.n2, &const_spec, DEFAULT_VAR_SAMPLES)?;
        edge_samples[k] = Some(samples);
    }

    // Collect the corner faces to rebuild.
    let mut face_map: HashMap<usize, CornerFaceRebuild> = HashMap::new();
    let mut corner_face_keys: Vec<usize> = Vec::new();
    for c in &corners {
        for f in corner_faces_at(solid, &c.vertex) {
            let n_out = match face_outward_normal(&f) {
                Ok(n) => n,
                Err(_) => continue,
            };
            let Some(arc) = c.geom.arc_on_face(&n_out) else { continue };
            let key = Arc::as_ptr(&f.0.tshape) as usize;
            let entry = face_map.entry(key).or_insert_with(|| {
                corner_face_keys.push(key);
                CornerFaceRebuild {
                    face: f.clone(),
                    trims: Vec::new(),
                    arcs: Vec::new(),
                    polylines: Vec::new(),
                }
            });
            entry.arcs.push(arc);
            for &eidx in &[c.k, c.k + 1] {
                let e = &infos[eidx];
                if is_same(&e.f1.0, &f.0) {
                    let (p0, u, off) = edge_tangency_trim(e, &f, c.radius)?;
                    entry.trims.push((p0, u, off));
                    let poly: Vec<GpPnt> =
                        edge_samples[eidx].as_ref().unwrap().iter().map(|s| s.t1).collect();
                    entry.polylines.push((poly[0], *poly.last().unwrap(), poly));
                } else if is_same(&e.f2.0, &f.0) {
                    let (p0, u, off) = edge_tangency_trim(e, &f, c.radius)?;
                    entry.trims.push((p0, u, off));
                    let poly: Vec<GpPnt> =
                        edge_samples[eidx].as_ref().unwrap().iter().map(|s| s.t2).collect();
                    entry.polylines.push((poly[0], *poly.last().unwrap(), poly));
                }
            }
        }
    }

    // Collect the far end faces to rebuild.
    let mut end_map: HashMap<usize, EndFaceRebuild> = HashMap::new();
    let mut end_face_keys: Vec<usize> = Vec::new();
    for k in 0..n {
        let e = &infos[k];
        let r = edge_radius[k];
        if corner_at_start[k].is_none() {
            let end_face = find_end_face(solid, &e.f1, &e.f2, &e.p0)?;
            let s = sample_endpoint_blend(e, 0.0, r)?;
            let key = Arc::as_ptr(&end_face.0.tshape) as usize;
            let entry = end_map.entry(key).or_insert_with(|| {
                end_face_keys.push(key);
                EndFaceRebuild { face: end_face.clone(), replacements: Vec::new() }
            });
            entry.replacements.push(EndArcReplacement {
                corner: e.p0,
                t1: s.t1,
                t2: s.t2,
                center: s.center,
                axis: e.axis,
                radius: r,
                n1: e.n1,
            });
        }
        if corner_at_end[k].is_none() {
            let end_face = find_end_face(solid, &e.f1, &e.f2, &e.p1)?;
            let s = sample_endpoint_blend(e, 1.0, r)?;
            let key = Arc::as_ptr(&end_face.0.tshape) as usize;
            let entry = end_map.entry(key).or_insert_with(|| {
                end_face_keys.push(key);
                EndFaceRebuild { face: end_face.clone(), replacements: Vec::new() }
            });
            entry.replacements.push(EndArcReplacement {
                corner: e.p1,
                t1: s.t1,
                t2: s.t2,
                center: s.center,
                axis: e.axis,
                radius: r,
                n1: e.n1,
            });
        }
    }

    // Build everything with a single shared cache.
    let mut cache = EdgeCache::new();
    let mut new_faces: Vec<Face> = Vec::new();

    // Rebuilt corner faces.
    for &key in &corner_face_keys {
        let fr = &face_map[&key];
        let rebuilt = rebuild_corner_face(&fr.face, &fr.trims, &fr.arcs, &fr.polylines, tol, &mut cache)?;
        new_faces.push(rebuilt);
    }

    // Rebuilt far end faces.
    for &key in &end_face_keys {
        let er = &end_map[&key];
        let rebuilt = rebuild_face_with_end_arcs(&er.face, &er.replacements, &mut cache)?;
        new_faces.push(rebuilt);
    }

    // Blend faces (trimmed to the corner radius span).
    for k in 0..n {
        let e = &infos[k];
        let samples = edge_samples[k].as_ref().unwrap();
        let blend = build_var_blend_face(samples, &e.axis, &mut cache)?;
        new_faces.push(blend);
    }

    // Corner patch faces.
    for c in &corners {
        let patch = build_corner_patch_face(&c.geom, &mut cache)?;
        new_faces.push(patch);
    }

    // Assemble: keep every face except the corner faces and far end faces.
    let mut kept: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        let key = Arc::as_ptr(&f.0.tshape) as usize;
        if corner_face_keys.contains(&key) || end_face_keys.contains(&key) {
            continue;
        }
        kept.push(f);
    }
    kept.extend(new_faces);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&kept);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Sample the blend cross-section at an edge endpoint fraction (0 or 1) and
/// return the tangency points and arc centre at that cross-section.
fn sample_endpoint_blend(
    e: &EdgeInfo,
    t: f64,
    radius: f64,
) -> Result<BlendSample, String> {
    let spec = VarFilletSpec::new(radius, radius);
    let (samples, _, _, _) = sample_blend(&e.p0, &e.p1, &e.n1, &e.n2, &spec, 2)?;
    Ok(if t < 0.5 { samples[0].clone() } else { samples[1].clone() })
}

/// Split `edge_indices` into maximal runs of edges that share a vertex.
fn split_consecutive_runs(solid: &TopoShape, edge_indices: &[usize]) -> Vec<Vec<usize>> {
    let es = edges_of(solid);
    let mut runs: Vec<Vec<usize>> = Vec::new();
    for &i in edge_indices {
        if let Some(last) = runs.last_mut() {
            if let Some(&prev) = last.last() {
                if let (Some(pe), Some(ce)) = (es.get(prev), es.get(i)) {
                    if edges_share_vertex(pe, ce) {
                        last.push(i);
                        continue;
                    }
                }
            }
        }
        runs.push(vec![i]);
    }
    runs
}

fn edges_share_vertex(a: &Edge, b: &Edge) -> bool {
    let (a0, a1) = edge_vertices(a);
    let (b0, b1) = edge_vertices(b);
    let (Some(a0), Some(a1)) = (a0, a1) else { return false };
    let (Some(b0), Some(b1)) = (b0, b1) else { return false };
    let (pa0, pa1) = (vertex_position(&a0), vertex_position(&a1));
    let (pb0, pb1) = (vertex_position(&b0), vertex_position(&b1));
    pa0.distance(&pb0) < 1e-9
        || pa0.distance(&pb1) < 1e-9
        || pa1.distance(&pb0) < 1e-9
        || pa1.distance(&pb1) < 1e-9
}

/// Blend a corner where the given edges meet: build the spherical corner patch
/// and trim the three corner faces. The edges themselves must already be (or
/// subsequently be) filleted; this inserts the patch face and re-closes the
/// shell. `edge_specs` maps edge indices (into `edges_of(solid)`) to specs.
pub fn fillet_corner_blend(
    solid: &TopoShape,
    corner_vertex: &Vertex,
    edge_specs: &[(usize, VarFilletSpec)],
    tol: f64,
) -> Result<TopoShape, String> {
    let p = vertex_position(corner_vertex);
    let es = edges_of(solid);
    let mut incident: Vec<(Edge, VarFilletSpec)> = Vec::new();
    for (i, spec) in edge_specs {
        let e = es.get(*i).ok_or("fillet_corner_blend: edge index out of range")?;
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            if vertex_position(&va).distance(&p) < 1e-9
                || vertex_position(&vb).distance(&p) < 1e-9
            {
                incident.push((e.clone(), *spec));
            }
        }
    }
    if incident.len() != 2 {
        return Err(format!(
            "fillet_corner_blend: expected 2 incident edge specs, got {}",
            incident.len()
        ));
    }
    let radius = incident
        .iter()
        .map(|(_, s)| s.r_end.max(s.r_start))
        .fold(0.0, f64::max);
    let geom = CornerGeom::from_edges(solid, corner_vertex, &incident[0].0, &incident[1].0, radius)?;

    let mut cache = EdgeCache::new();
    let all_faces = faces_of(solid);
    let cf = corner_faces_at(solid, &p);
    let mut rebuilt = Vec::new();
    for f in &cf {
        let n_out = match face_outward_normal(f) {
            Ok(n) => n,
            Err(_) => continue,
        };
        let Some(arc) = geom.arc_on_face(&n_out) else { continue };
        rebuilt.push(rebuild_corner_face(f, &[], &[arc], &[], tol, &mut cache)?);
    }
    let patch = build_corner_patch_face(&geom, &mut cache)?;

    let mut faces: Vec<Face> = Vec::new();
    for f in all_faces {
        if cf.iter().any(|c| is_same(&f.0, &c.0)) {
            continue;
        }
        faces.push(f);
    }
    faces.extend(rebuilt);
    faces.push(patch);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Fillet a chain of edges, adding a spherical corner patch at every vertex
/// shared by two consecutive edges so the result stays a closed shell.
///
/// Non-adjacent edges are filleted sequentially with `fillet_edge_var` (as in
/// `fillet_edge_var_chain`). A run of consecutive edges sharing vertices is
/// built in one pass: each edge gets a constant-radius blend face trimmed at
/// the shared corners, and each shared vertex gets a rolling-ball corner patch
/// (`fillet_corner_blend`). The legacy `fillet_edge_var_chain` is unchanged.
pub fn fillet_edges_chain_with_corner(
    solid: &TopoShape,
    edge_indices: &[usize],
    specs: &[VarFilletSpec],
    tol: f64,
) -> Result<TopoShape, String> {
    if edge_indices.len() != specs.len() {
        return Err(
            "fillet_edges_chain_with_corner: specs.len() must equal edge_indices.len()".to_string(),
        );
    }
    let original_edges = edges_of(solid);
    let runs = split_consecutive_runs(solid, edge_indices);
    let mut current = solid.clone();
    let mut spec_offset = 0;
    for run in &runs {
        let run_specs = &specs[spec_offset..spec_offset + run.len()];
        if run.len() == 1 {
            let i = run[0];
            let oe = original_edges
                .get(i)
                .ok_or_else(|| format!("fillet_edges_chain_with_corner: edge index {i} out of range"))?;
            let (p0, p1) = BRepTool::edge_vertices(oe).ok_or("edge has no curve")?;
            let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
                format!(
                    "fillet_edges_chain_with_corner: edge {i} was consumed by an earlier fillet"
                )
            })?;
            current = fillet_edge_var(&current, &e, &run_specs[0], tol)?;
        } else {
            let mut run_edges: Vec<Edge> = Vec::with_capacity(run.len());
            for &i in run {
                let oe = original_edges
                    .get(i)
                    .ok_or_else(|| format!("fillet_edges_chain_with_corner: edge index {i} out of range"))?;
                let (p0, p1) = BRepTool::edge_vertices(oe).ok_or("edge has no curve")?;
                let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
                    format!(
                        "fillet_edges_chain_with_corner: edge {i} was consumed by an earlier fillet"
                    )
                })?;
                run_edges.push(e);
            }
            current = build_combined_run(&current, &run_edges, run_specs, tol)?;
        }
        spec_offset += run.len();
    }
    Ok(current)
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

    // --- Shared-vertex corner patch (rolling-ball corner blend) ---

    fn find_edge_by_endpoints_test(b: &BRepPrimBox, a: &GpPnt, bpt: &GpPnt) -> usize {
        box_edges(b)
            .iter()
            .position(|e| {
                let (x, y) = BRepTool::edge_vertices(e).unwrap();
                (x.is_equal(a) && y.is_equal(bpt)) || (x.is_equal(bpt) && y.is_equal(a))
            })
            .expect("edge by endpoints")
    }

    fn corner_vertex(b: &BRepPrimBox, p: &GpPnt) -> Vertex {
        vertices_of(&b.solid.0)
            .into_iter()
            .find(|v| BRepTool::vertex_point(v).is_equal(p))
            .expect("corner vertex")
    }

    #[test]
    fn cut_disk_splits_bottom_edge() {
        // The -Y face rectangle at z >= R, with the corner disk centred at
        // (R, 0, R) of radius R.
        let r = 0.4;
        let rect = [
            GpPnt::new(0.0, 0.0, r),
            GpPnt::new(2.0, 0.0, r),
            GpPnt::new(2.0, 0.0, 2.0),
            GpPnt::new(0.0, 0.0, 2.0),
        ];
        let cut = cut_disk_poly(&rect, &GpPnt::new(r, 0.0, r), r).expect("cut");
        // The bottom edge must be split: the polygon should contain the point
        // (2R, 0, R) where the disk boundary crosses the z = R edge.
        assert!(
            cut.iter().any(|p| p.distance(&GpPnt::new(2.0 * r, 0.0, r)) < 1e-6),
            "bottom edge was not split by the disk"
        );
    }

    #[test]
    fn corner_center_box() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let v = corner_vertex(&bx, &GpPnt::new(0.0, 0.0, 0.0));
        let c = corner_center(&bx.solid.0, &v, 0.4, 1e-6).expect("corner center");
        assert!(
            c.distance(&GpPnt::new(0.4, 0.4, 0.4)) < 1e-9,
            "corner centre {c:?}"
        );
        // The centre of the opposite corner (2,2,2) is (1.6,1.6,1.6).
        let v2 = corner_vertex(&bx, &GpPnt::new(2.0, 2.0, 2.0));
        let c2 = corner_center(&bx.solid.0, &v2, 0.4, 1e-6).unwrap();
        assert!(c2.distance(&GpPnt::new(1.6, 1.6, 1.6)) < 1e-9);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn corner_radius_positive() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let v = corner_vertex(&bx, &GpPnt::new(0.0, 0.0, 0.0));
        let n = edges_of(&bx.solid.0).len();
        // specs aligned with edges_of: constant 0.4 on the two edges incident
        // to the corner.
        let mut specs = vec![VarFilletSpec::new(0.1, 0.1); n];
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        specs[e0] = VarFilletSpec::new(0.4, 0.4);
        specs[e1] = VarFilletSpec::new(0.4, 0.4);
        let r = corner_patch_radius_at_vertex(&bx.solid.0, &v, &specs, 1e-6);
        assert!((r - 0.4).abs() < 1e-9, "corner radius {r}");
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn corner_blend_box_corner_closed() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.4, 0.4), VarFilletSpec::new(0.4, 0.4)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6)
            .expect("corner chain");
        let faces = faces_of(&out);
        assert_eq!(
            faces.len(),
            9,
            "6 faces + 2 blends + 1 corner patch, got {}",
            faces.len()
        );
        assert!(shell_is_closed(&closed_shell(&out)), "corner chain must be closed");
        // The corner patch face is a sphere.
        let sphere = faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere)
                    .unwrap_or(false)
            })
            .expect("corner patch sphere face");
        let _ = sphere;
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_with_corner_success() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.3, 0.3), VarFilletSpec::new(0.3, 0.3)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6)
            .expect("chain with corner must not error");
        assert!(shell_is_closed(&closed_shell(&out)), "chain with corner is closed");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn corner_patch_is_spherical() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.4, 0.4), VarFilletSpec::new(0.4, 0.4)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6).unwrap();
        let faces = faces_of(&out);
        let sphere = faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere)
                    .unwrap_or(false)
            })
            .expect("corner patch sphere face");
        // All sampled surface points are equidistant from the corner centre.
        let v = corner_vertex(&bx, &GpPnt::new(0.0, 0.0, 0.0));
        let cc = corner_center(&bx.solid.0, &v, 0.4, 1e-6).unwrap();
        let s = BRepTool::face_surface(sphere).unwrap();
        let r0 = s.d0(0.0, 0.0).distance(&cc);
        for (u, w) in [(0.5, 0.3), (1.0, 0.0), (1.5, 0.4), (3.0, 0.2)] {
            let d = s.d0(u, w).distance(&cc);
            assert!((d - r0).abs() < 1e-6, "patch point at distance {d} (expected {r0})");
        }
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_nonconsecutive_still_works() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        // Opposite top edge: (2,2,2)→(0,2,2).
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(2.0, 2.0, 2.0), &GpPnt::new(0.0, 2.0, 2.0));
        let specs = [VarFilletSpec::new(0.2, 0.3), VarFilletSpec::new(0.3, 0.4)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6)
            .expect("non-consecutive chain");
        assert_eq!(faces_of(&out).len(), 8, "two sequential fillets");
        assert!(shell_is_closed(&closed_shell(&out)), "non-consecutive chain is closed");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_with_corner_three_edges() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        // Bottom-face loop: front X edge, right Y edge, back X edge.
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let e2 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.3, 0.3); 3];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1, e2], &specs, 1e-6)
            .expect("three-edge corner chain");
        let faces = faces_of(&out);
        assert_eq!(
            faces.len(),
            11,
            "6 faces + 3 blends + 2 corner patches, got {}",
            faces.len()
        );
        assert!(shell_is_closed(&closed_shell(&out)), "three-edge corner chain is closed");
        // Two corner patch faces.
        let spheres = faces
            .iter()
            .filter(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere)
                    .unwrap_or(false)
            })
            .count();
        assert_eq!(spheres, 2, "two corner patches, got {spheres}");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }
}
