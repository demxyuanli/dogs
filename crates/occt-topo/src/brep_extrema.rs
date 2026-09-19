//! Shape-to-shape distance and extrema queries.
//!
//! Port of `BRepExtrema_DistShapeShape` (TKBRep) plus a few lightweight
//! helpers (closest vertex, point-in-solid, bbox separation).
//!
//! `BRepExtrema_DistShapeShape` computes exact distances between shapes via
//! extremum search on the underlying geometry. This port delegates the
//! point–curve and point–surface minima to the Phase 13 analytic/Newton
//! solvers in `occt_geom` (`extrema_pc` / `extrema_surf`), then clamps the
//! parameters to the bounded edge / face range and re-checks the face
//! boundary edges (a bounded face's closest point can lie on an edge).

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpXyz};
use occt_geom::extrema::point_surface_extrema;
use occt_geom::extrema_pc::point_curve_extrema_all;
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topexp::Explorer;

/// Minimum distance from a point to any part of a shape.
///
/// Exact for vertices; exact analytic / Newton for edges and faces (via
/// [`closest_point_on_edge`] / [`closest_point_on_face`]).
pub fn point_shape_distance(p: &GpPnt, shape: &TopoShape) -> f64 {
    let mut best = f64::INFINITY;

    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let vs = ex.current().clone();
        ex.next();
        if let Some(v) = Vertex::wrap(vs) {
            best = best.min(BRepTool::vertex_point(&v).distance(p));
        }
    }

    let mut ex = Explorer::new(shape, ShapeType::Edge);
    while ex.more() {
        let es = ex.current().clone();
        ex.next();
        if let Some(e) = Edge::wrap(es) {
            let (_, q) = closest_point_on_edge(&e, p, 32);
            best = best.min(q.distance(p));
        }
    }

    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let fs = ex.current().clone();
        ex.next();
        if let Some(f) = Face::wrap(fs) {
            let (_, q) = closest_point_on_face(&f, p, 16, 16);
            best = best.min(q.distance(p));
        }
    }

    best
}

/// Approximate minimum distance between two shapes.
///
/// The distance is minimized over vertex–vertex, vertex–edge, vertex–face,
/// edge–vertex and edge–face sample pairs. Edge–edge pairs are covered
/// indirectly through the sampled edge points vs. faces. The result is a
/// coarse upper bound on the true distance; it is not an exact extremum
/// computation.
pub fn shape_distance(a: &TopoShape, b: &TopoShape) -> f64 {
    let va = collect_vertices(a);
    let vb = collect_vertices(b);
    let ea = collect_edges(a);
    let eb = collect_edges(b);
    let fa = collect_faces(a);
    let fb = collect_faces(b);

    let mut best = f64::INFINITY;

    for av in &va {
        for bv in &vb {
            best = best.min(BRepTool::vertex_point(av).distance(&BRepTool::vertex_point(bv)));
        }
    }
    for av in &va {
        let p = BRepTool::vertex_point(av);
        for be in &eb {
            let (_, q) = closest_point_on_edge(be, &p, 16);
            best = best.min(q.distance(&p));
        }
    }
    for av in &va {
        let p = BRepTool::vertex_point(av);
        for bf in &fb {
            let (_, q) = closest_point_on_face(bf, &p, 8, 8);
            best = best.min(q.distance(&p));
        }
    }
    for bv in &vb {
        let p = BRepTool::vertex_point(bv);
        for ae in &ea {
            let (_, q) = closest_point_on_edge(ae, &p, 16);
            best = best.min(q.distance(&p));
        }
    }
    for ae in &ea {
        let (f, l) = BRepTool::edge_parameters(ae);
        if !f.is_finite() || !l.is_finite() {
            continue;
        }
        let curve = BRepTool::edge_curve(ae);
        if curve.is_none() {
            continue;
        }
        let curve = curve.unwrap();
        for i in 0..=8 {
            let u = f + (l - f) * (i as f64 / 8.0);
            let q = curve.d0(u);
            for bf in &fb {
                let (_, fq) = closest_point_on_face(bf, &q, 8, 8);
                best = best.min(fq.distance(&q));
            }
        }
    }

    best
}

/// The nearest vertex of `shape` to `p`.
pub fn closest_vertex(shape: &TopoShape, p: &GpPnt) -> Option<Vertex> {
    let mut best: Option<Vertex> = None;
    let mut best_d = f64::INFINITY;
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let vs = ex.current().clone();
        ex.next();
        if !seen.insert(Arc::as_ptr(&vs.tshape) as usize) {
            continue;
        }
        if let Some(v) = Vertex::wrap(vs) {
            let d = BRepTool::vertex_point(&v).distance(p);
            if d < best_d {
                best_d = d;
                best = Some(v);
            }
        }
    }
    best
}

/// Closest point on an edge to `p`, returned as `(parameter, point)`.
///
/// Solves the point–curve extrema with `occt_geom`'s analytic/Newton solver
/// ([`point_curve_extrema_all`]), then clamps the winning parameter to the
/// edge's bounded parameter range (matching the previous sampling behaviour for
/// unbounded curves: a generous `±1e6` window). Every returned extremum is
/// clamped and evaluated so the global minimum of the *bounded* edge is found
/// even when the unbounded curve's minimum lies outside the edge. For an edge
/// without a curve the parameter is `NaN`.
pub fn closest_point_on_edge(edge: &Edge, p: &GpPnt, _samples: usize) -> (f64, GpPnt) {
    let curve = match BRepTool::edge_curve(edge) {
        Some(c) => c,
        None => return (f64::NAN, GpPnt::zero()),
    };
    let (mut f, mut l) = BRepTool::edge_parameters(edge);
    if !f.is_finite() || !l.is_finite() {
        // Unbounded curve (infinite line): clamp to a generous window.
        f = -1e6;
        l = 1e6;
    }
    let (f, l) = if f <= l { (f, l) } else { (l, f) };
    let extrema = point_curve_extrema_all(&*curve, p);
    if extrema.is_empty() {
        // Degenerate curve (e.g. zero-length): golden-section fallback.
        let (u, q) = occt_geom::extrema::refine_curve_point(&*curve, p, f, l);
        return (u, q);
    }
    let mut best_u = f;
    let mut best_d = f64::INFINITY;
    let mut consider = |u: f64| {
        let u = u.clamp(f, l);
        let d = curve.d0(u).square_distance(p);
        if d < best_d {
            best_d = d;
            best_u = u;
        }
    };
    for e in &extrema {
        consider(e.u1);
    }
    // The bounded-edge candidate must also include both endpoints; the solver
    // reports them only when the underlying curve's range is finite.
    consider(f);
    consider(l);
    (best_u, curve.d0(best_u))
}

/// Closest point on a face to `p`, returned as `((u, v), point)`.
///
/// Solves the point–surface extrema with `occt_geom`'s analytic/Newton solver
/// ([`point_surface_extrema`]), clamps the parameter to the face's UV window,
/// and then also checks each boundary edge (a bounded face's closest point can
/// lie on an edge; the previous UV-grid sampler found those implicitly). When
/// an edge wins, its `(u, v)` is recovered by projecting the 3D point back onto
/// the surface.
pub fn closest_point_on_face(face: &Face, p: &GpPnt, _nu: usize, _nv: usize) -> ((f64, f64), GpPnt) {
    let surface = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return ((0.0, 0.0), GpPnt::zero()),
    };
    let (u1, u2, v1, v2) = match face_uv_window(face) {
        Some(w) => w,
        None => return ((0.0, 0.0), GpPnt::zero()),
    };
    let e = point_surface_extrema(&*surface, p);
    // The analytic solver reports its own parameter frame (a plane is
    // reconstructed with arbitrary X/Y axes), so take its frame-independent
    // closest POINT and re-project it into the face's own (u, v) before
    // clamping to the face's UV window.
    let mut best_u = u1;
    let mut best_v = v1;
    let mut best_p = surface.d0(best_u, best_v);
    let mut best_d = best_p.square_distance(p);
    let (pu, pv) = project_to_surface(&surface, &e.p2);
    let pu = pu.clamp(u1, u2);
    let pv = pv.clamp(v1, v2);
    let q = surface.d0(pu, pv);
    let d = q.square_distance(p);
    if d < best_d {
        best_d = d;
        best_u = pu;
        best_v = pv;
        best_p = q;
    }
    for edge in face_boundary_edges(face) {
        let (_, q) = closest_point_on_edge(&edge, p, 0);
        let d = q.square_distance(p);
        if d < best_d {
            let (pu, pv) = project_to_surface(&surface, &q);
            best_d = d;
            best_u = pu;
            best_v = pv;
            best_p = q;
        }
    }
    ((best_u, best_v), best_p)
}

/// Even-odd point-in-solid test.
///
/// Casts a ray in +X and counts crossings against a coarse surface mesh of the
/// faces (grid-triangulated). Points strictly outside the bounding box are
/// rejected first. If the shape cannot be meshed (no faces with surfaces), the
/// bounding-box containment is returned as a fallback. Best suited to closed
/// solids; results on open surfaces are not meaningful.
pub fn is_inside(shape: &TopoShape, p: &GpPnt) -> bool {
    let bbox = shape_bbox(shape);
    if let Some(b) = &bbox {
        if b.is_out(p) {
            return false;
        }
    }

    let meshes = mesh_faces(shape, 7, 7);
    if meshes.is_empty() {
        // Unmeshed shape: fall back to bounding-box containment.
        return match bbox {
            Some(b) => !b.is_out(p),
            None => false,
        };
    }

    // Jitter the ray origin by a tiny amount perpendicular to the +X ray so it
    // does not land exactly on a grid cell diagonal (a degenerate ray that can
    // be double-counted or skipped by both triangles of a cell).
    let jittered = GpPnt::new(p.x(), p.y() + 1e-7, p.z() + 1e-7);
    let mut hits: Vec<f64> = Vec::new();
    for (verts, tris) in &meshes {
        for (i, j, k) in tris {
            if let Some(t) = ray_plus_x_hits(&verts[*i], &verts[*j], &verts[*k], &jittered) {
                hits.push(t);
            }
        }
    }
    // A ray passing exactly through a shared triangle edge reports two hits at
    // the same parameter; deduplicate coincident hits so the parity stays
    // correct for closed meshes.
    hits.sort_by(f64::total_cmp);
    let mut unique = 0usize;
    let mut prev: Option<f64> = None;
    for t in hits {
        if prev.map_or(true, |q| (t - q).abs() > 1e-9) {
            unique += 1;
            prev = Some(t);
        }
    }
    unique % 2 == 1
}

/// Minimum separation between two axis-aligned bounding boxes; 0 if they
/// overlap.
pub fn shape_bbox_distance(a: &BndBox, b: &BndBox) -> f64 {
    let (a0, a1, ay0, ay1, az0, az1) = match a.get() {
        Some(v) => v,
        None => return f64::INFINITY,
    };
    let (b0, b1, by0, by1, bz0, bz1) = match b.get() {
        Some(v) => v,
        None => return f64::INFINITY,
    };
    let gx = axis_gap(a0, a1, b0, b1);
    let gy = axis_gap(ay0, ay1, by0, by1);
    let gz = axis_gap(az0, az1, bz0, bz1);
    (gx * gx + gy * gy + gz * gz).sqrt()
}

fn axis_gap(amin: f64, amax: f64, bmin: f64, bmax: f64) -> f64 {
    if amax < bmin {
        return bmin - amax;
    }
    if bmax < amin {
        return amin - bmax;
    }
    0.0
}

// ---------------------------------------------------------------------------
// Fallback shape → face-mesh sampling.
//
// `crate::shape_mesh::mesh_shape` is being written concurrently; until it
// lands, brep_gprop and is_inside use this local grid triangulator. Each face
// is sampled on a regular `(nu+1)×(nv+1)` grid in parametric space. Faces
// with an unbounded parametric range (e.g. planes) get a window inferred from
// their boundary vertices projected onto the surface.
// ---------------------------------------------------------------------------

/// Triangulate every face of `shape` into `(vertices, triangles)` pairs.
///
/// Vertices are `GpPnt`s; triangles are index triples into the corresponding
/// vertex list. Faces without a registered surface (or without any boundary
/// from which to infer an unbounded UV window) are skipped.
pub fn mesh_faces(
    shape: &TopoShape,
    nu: usize,
    nv: usize,
) -> Vec<(Vec<GpPnt>, Vec<(usize, usize, usize)>)> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let fs = ex.current().clone();
        ex.next();
        if let Some(face) = Face::wrap(fs) {
            if let Some(surface) = BRepTool::face_surface(&face) {
                if let Some((u1, u2, v1, v2)) = face_uv_window(&face) {
                    let (verts, mut tris) = sample_grid(&surface, u1, u2, v1, v2, nu, nv);
                    // The grid follows the surface parameterization; a REVERSED
                    // face has its material on the opposite side, so its
                    // triangles must wind the other way (OCCT `BRep_Tool::
                    // Triangulation` + the face orientation the caller applies).
                    if face.orientation().is_reversed() {
                        for t in &mut tris {
                            std::mem::swap(&mut t.1, &mut t.2);
                        }
                    }
                    out.push((verts, tris));
                }
            }
        }
    }
    out
}

/// Sample the surface on a regular parametric grid.
fn sample_grid(
    surface: &Arc<dyn Surface>,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    nu: usize,
    nv: usize,
) -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
    let (nu, nv) = (nu.max(1), nv.max(1));
    let mut verts = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        let u = u1 + (u2 - u1) * (i as f64 / nu as f64);
        for j in 0..=nv {
            let v = v1 + (v2 - v1) * (j as f64 / nv as f64);
            verts.push(surface.d0(u, v));
        }
    }
    let stride = nv + 1;
    let mut tris = Vec::with_capacity(2 * nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            tris.push((a, b, c));
            tris.push((b, d, c));
        }
    }
    (verts, tris)
}

/// A finite parametric window `(u_min, u_max, v_min, v_max)` for a face.
///
/// Uses the surface's own range when finite; otherwise projects the boundary
/// vertices of the face's wires onto the surface (Newton solve) and takes the
/// bounding box of the resulting `(u, v)` points.
fn face_uv_window(face: &Face) -> Option<(f64, f64, f64, f64)> {
    let surface = BRepTool::face_surface(face)?;
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(face);
    if u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite() {
        return Some((u1, u2, v1, v2));
    }
    let pts = face_boundary_points(face);
    if pts.is_empty() {
        return None;
    }
    let mut us: Vec<f64> = Vec::new();
    let mut vs: Vec<f64> = Vec::new();
    for p in &pts {
        let (u, v) = project_to_surface(&surface, p);
        us.push(u);
        vs.push(v);
    }
    let ua = us.iter().copied().fold(f64::INFINITY, f64::min);
    let ub = us.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let va = vs.iter().copied().fold(f64::INFINITY, f64::min);
    let vb = vs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (ua, ub) = if ub - ua < 1e-12 { (ua - 1.0, ub + 1.0) } else { (ua, ub) };
    let (va, vb) = if vb - va < 1e-12 { (va - 1.0, vb + 1.0) } else { (va, vb) };
    Some((ua, ub, va, vb))
}

/// Project a 3D point onto a surface by Newton iteration on the `d1` partials.
///
/// Converges in a single step for planes (affine mapping); a few steps for
/// other analytic surfaces.
fn project_to_surface(surface: &Arc<dyn Surface>, p: &GpPnt) -> (f64, f64) {
    let mut u = 0.0;
    let mut v = 0.0;
    for _ in 0..20 {
        let (pt, du, dv) = surface.d1(u, v);
        let rx = p.x() - pt.x();
        let ry = p.y() - pt.y();
        let rz = p.z() - pt.z();
        let duu = du.coord.square_modulus();
        let dvv = dv.coord.square_modulus();
        let duv = du.coord.dot(&dv.coord);
        let b1 = rx * du.x() + ry * du.y() + rz * du.z();
        let b2 = rx * dv.x() + ry * dv.y() + rz * dv.z();
        let det = duu * dvv - duv * duv;
        if det.abs() < 1e-30 {
            break;
        }
        let d_u = (b1 * dvv - b2 * duv) / det;
        let d_v = (b2 * duu - b1 * duv) / det;
        u += d_u;
        v += d_v;
        if (d_u * d_u + d_v * d_v).sqrt() < 1e-12 {
            break;
        }
    }
    (u, v)
}

/// All edges of a face's wires (boundary edges).
fn face_boundary_edges(face: &Face) -> Vec<Edge> {
    let mut out: Vec<Edge> = Vec::new();
    let face_kids = face.0.tshape.read().unwrap().children.clone();
    for wire in face_kids {
        if wire.shape_type() != ShapeType::Wire {
            continue;
        }
        let wire_kids = wire.tshape.read().unwrap().children.clone();
        for e in wire_kids {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            if let Some(edge) = Edge::wrap(e) {
                out.push(edge);
            }
        }
    }
    out
}

/// Distinct 3D endpoint points of all edges in a face's wires.
fn face_boundary_points(face: &Face) -> Vec<GpPnt> {
    let mut out: Vec<GpPnt> = Vec::new();
    let face_kids = face.0.tshape.read().unwrap().children.clone();
    for wire in face_kids {
        if wire.shape_type() != ShapeType::Wire {
            continue;
        }
        let wire_kids = wire.tshape.read().unwrap().children.clone();
        for e in wire_kids {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            if let Some(edge) = Edge::wrap(e) {
                if let Some((a, b)) = BRepTool::edge_vertices(&edge) {
                    if !out.iter().any(|q| q.distance(&a) < 1e-9) {
                        out.push(a);
                    }
                    if !out.iter().any(|q| q.distance(&b) < 1e-9) {
                        out.push(b);
                    }
                }
            }
        }
    }
    out
}

/// Axis-aligned bounding box of a shape's vertices.
fn shape_bbox(shape: &TopoShape) -> Option<BndBox> {
    let mut b = BndBox::new();
    let mut any = false;
    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let vs = ex.current().clone();
        ex.next();
        if let Some(v) = Vertex::wrap(vs) {
            b.add_point(&BRepTool::vertex_point(&v));
            any = true;
        }
    }
    if any {
        Some(b)
    } else {
        None
    }
}

/// Möller–Trumbore ray/triangle test for a ray from `origin` in the +X
/// direction.
/// Möller–Trumbore ray (+X) / triangle intersection; returns the hit parameter
/// `t` when the ray crosses the triangle.
fn ray_plus_x_hits(a: &GpPnt, b: &GpPnt, c: &GpPnt, origin: &GpPnt) -> Option<f64> {
    const EPS: f64 = 1e-9;
    let dir = GpXyz::new(1.0, 0.0, 0.0);
    let edge1 = b.coord.subtracted(&a.coord);
    let edge2 = c.coord.subtracted(&a.coord);
    let h = dir.crossed(&edge2);
    let det = edge1.dot(&h);
    if det.abs() < EPS {
        return None;
    }
    let inv = 1.0 / det;
    let s = origin.coord.subtracted(&a.coord);
    let u = s.dot(&h) * inv;
    if u < -EPS || u > 1.0 + EPS {
        return None;
    }
    let q = s.crossed(&edge1);
    let v = dir.dot(&q) * inv;
    if v < -EPS || u + v > 1.0 + EPS {
        return None;
    }
    let t = edge2.dot(&q) * inv;
    if t > EPS { Some(t) } else { None }
}

fn collect_vertices(shape: &TopoShape) -> Vec<Vertex> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let s = ex.current().clone();
        ex.next();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            if let Some(v) = Vertex::wrap(s) {
                out.push(v);
            }
        }
    }
    out
}

fn collect_edges(shape: &TopoShape) -> Vec<Edge> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Edge);
    while ex.more() {
        let s = ex.current().clone();
        ex.next();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            if let Some(e) = Edge::wrap(s) {
                out.push(e);
            }
        }
    }
    out
}

fn collect_faces(shape: &TopoShape) -> Vec<Face> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let s = ex.current().clone();
        ex.next();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            if let Some(f) = Face::wrap(s) {
                out.push(f);
            }
        }
    }
    out
}

#[cfg(test)]
pub(crate) mod test_box {
    //! Shared unit-box fixture for the Phase 3 test modules.
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
    use occt_geom::GeomLine;

    use crate::builder::TopoBuilder;
    use crate::abs::Orientation;
    use crate::shape::{Edge, Face, Solid, Vertex};

    pub(crate) struct UnitBox {
        pub solid: Solid,
        pub corners: [GpPnt; 8],
        pub vertices: Vec<Vertex>,
        pub edges: Vec<Edge>,
        pub faces: Vec<Face>,
    }

    fn segment_edge(
        b: &TopoBuilder,
        p1: &GpPnt,
        p2: &GpPnt,
        v1: &Vertex,
        v2: &Vertex,
    ) -> Edge {
        let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).expect("distinct points");
        let lin = GpLin::from_pnt_dir(*p1, dir);
        let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
        b.add_edge_vertices(&mut e, v1, v2);
        e
    }

    /// Build a closed unit box solid spanning `[0,1]³`.
    ///
    /// 8 shared vertices, 12 segment edges (each referencing two shared
    /// vertices), 6 planar faces (each with a closed wire of 4 edges), one
    /// shell, one solid. The face planes are oriented so their `(u,v)` grid
    /// maps `[0,1]²` onto each face quad with a consistent (inward) normal,
    /// which keeps the mesh signed-volume and centroid computations coherent.
    pub(crate) fn unit_box() -> UnitBox {
        let b = TopoBuilder::new();
        let corners = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(0.0, 0.0, 1.0),
            GpPnt::new(1.0, 0.0, 1.0),
            GpPnt::new(1.0, 1.0, 1.0),
            GpPnt::new(0.0, 1.0, 1.0),
        ];
        let vertices: Vec<Vertex> = corners.iter().map(|p| b.make_vertex(*p, 1e-7)).collect();

        let edge_idx: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges = Vec::new();
        for &(i, j) in &edge_idx {
            edges.push(segment_edge(&b, &corners[i], &corners[j], &vertices[i], &vertices[j]));
        }

        // Each face's 4 edges (indices into `edges`), sharing edges across faces.
        let face_edge_sets: [[usize; 4]; 6] = [
            [0, 1, 2, 3],     // bottom (z = 0)
            [4, 5, 6, 7],     // top (z = 1)
            [0, 9, 4, 8],     // front (y = 0)
            [2, 10, 6, 11],   // back (y = 1)
            [3, 11, 7, 8],    // left (x = 0)
            [1, 10, 5, 9],    // right (x = 1)
        ];
        // (origin, outward normal, u-direction) per face.
        let face_planes: [(GpPnt, GpDir, GpDir); 6] = [
            (GpPnt::new(0.0, 0.0, 0.0), GpDir::new(0.0, 0.0, -1.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
            (GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(0.0, 0.0, 0.0), GpDir::new(0.0, -1.0, 0.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(0.0, 1.0, 0.0), GpDir::new(0.0, 1.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(0.0, 0.0, 0.0), GpDir::new(-1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(1.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
        ];
        // Outward corner cycle of every face, in the same order as
        // `BRepPrim_GWedge` / `BRepPrimBox::make_box` (`primitives.rs`). An
        // edge that the face closes against the `edge_idx` chord is stored
        // Reversed, so two faces sharing an edge always see opposite
        // orientations — the invariant `BOPTools_AlgoTools::GetEdgeOff`
        // (`BOPTools_AlgoTools.cxx:1099-1127`) and
        // `BOPAlgo_ShellSplitter::SplitBlock` (`BOPAlgo_ShellSplitter.cxx:319`)
        // rely on.
        let face_cycles: [[usize; 4]; 6] = [
            [0, 1, 2, 3], // bottom (z = 0), reversed below to face -Z
            [4, 5, 6, 7], // top (z = 1)
            [0, 1, 5, 4], // front (y = 0)
            [3, 7, 6, 2], // back (y = 1)
            [0, 4, 7, 3], // left (x = 0)
            [1, 2, 6, 5], // right (x = 1)
        ];
        let mut faces = Vec::new();
        for fi in 0..6 {
            let (origin, normal, u_dir) = face_planes[fi];
            let ax3 = GpAx3::new(origin, normal, &u_dir).expect("perpendicular axes");
            let mut face = b.make_face_plane(&GpPln::new(ax3));
            let cycle = face_cycles[fi];
            let quad: Vec<Edge> = face_edge_sets[fi]
                .iter()
                .map(|&ei| {
                    let (i, j) = edge_idx[ei];
                    // Directed chord of the outward cycle closed by this edge.
                    let k = (0..4)
                        .find(|&k| {
                            let a = cycle[k];
                            let b = cycle[(k + 1) % 4];
                            (a == i && b == j) || (a == j && b == i)
                        })
                        .expect("edge belongs to the face cycle");
                    let (a, b) = (cycle[k], cycle[(k + 1) % 4]);
                    let e = edges[ei].clone();
                    if (a, b) == (i, j) {
                        e
                    } else {
                        Edge(e.0.oriented(Orientation::Reversed))
                    }
                })
                .collect();
            let mut wire = b.make_wire(&quad);
            // Bottom plane: origin (0,0,0), normal -Z, u = +Y so v = +X.
            // The +X then +Y edge loop is clockwise in that UV frame; reverse
            // the wire so the outer ring is CCW (material on the left), matching
            // BRepPrimAPI_MakeBox / IntTools_FClass2d::IsHole.
            if fi == 0 {
                wire.0.reverse();
            }
            b.add_wire(&mut face, &wire);
            faces.push(face);
        }
        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        UnitBox { solid, corners, vertices, edges, faces }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::bnd::BndBox;

    fn box_solid() -> test_box::UnitBox {
        test_box::unit_box()
    }

    #[test]
    fn point_to_box_distance() {
        let b = box_solid();
        let d = point_shape_distance(&GpPnt::new(2.0, 0.0, 0.0), &b.solid.0);
        assert!((d - 1.0).abs() < 1e-6, "distance {d}");
    }

    #[test]
    fn closest_vertex_is_corner() {
        let b = box_solid();
        let v = closest_vertex(&b.solid.0, &GpPnt::new(2.0, 0.0, 0.0)).expect("a vertex");
        let p = BRepTool::vertex_point(&v);
        assert!(p.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9, "vertex at {p:?}");
    }

    #[test]
    fn point_in_box() {
        let b = box_solid();
        assert!(is_inside(&b.solid.0, &GpPnt::new(0.5, 0.5, 0.5)));
        assert!(!is_inside(&b.solid.0, &GpPnt::new(2.0, 0.0, 0.0)));
        assert!(!is_inside(&b.solid.0, &GpPnt::new(-0.5, 0.5, 0.5)));
    }

    #[test]
    fn closest_on_edge_and_face() {
        let b = box_solid();
        let edge = &b.edges[0]; // (0,0,0)-(1,0,0)
        let (u, q) = closest_point_on_edge(edge, &GpPnt::new(0.5, 1.0, 0.0), 32);
        assert!((u - 0.5).abs() < 1e-3, "u {u}");
        assert!(q.distance(&GpPnt::new(0.5, 0.0, 0.0)) < 1e-6, "q {q:?}");

        let face = &b.faces[2]; // front (y = 0)
        let ((u, v), q) = closest_point_on_face(face, &GpPnt::new(0.5, 1.0, 0.5), 8, 8);
        assert!((u - 0.5).abs() < 1e-2, "u {u}");
        assert!((v - 0.5).abs() < 1e-2, "v {v}");
        assert!(q.distance(&GpPnt::new(0.5, 0.0, 0.5)) < 1e-6, "q {q:?}");
    }

    #[test]
    fn shape_self_distance_zero() {
        let b = box_solid();
        let d = shape_distance(&b.solid.0, &b.solid.0);
        assert!(d < 1e-6, "distance {d}");
    }

    #[test]
    fn bbox_separation() {
        let a = BndBox::from_corners(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b = BndBox::from_corners(&GpPnt::new(3.0, 0.0, 0.0), &GpPnt::new(4.0, 1.0, 1.0));
        assert!((shape_bbox_distance(&a, &b) - 2.0).abs() < 1e-12);
        let c = BndBox::from_corners(&GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(2.0, 2.0, 2.0));
        assert_eq!(shape_bbox_distance(&a, &c), 0.0);
    }
}
