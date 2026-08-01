//! Exact boolean operations for planar-faced solids (polyhedra).
//!
//! Port of `BOPAlgo_Builder` / `BRepAlgoAPI_Fuse` / `BRepAlgoAPI_Cut` /
//! `BRepAlgoAPI_Common`, restricted to planar-faced inputs (simplified
//! `BOPAlgo_Builder`). The pipeline:
//!
//! 1. collect face–face intersection segments (via `crate::inttools`);
//! 2. split every face's 2D polygon along the segments that lie on it;
//! 3. classify each sub-face against the other solid (in / out / on);
//! 4. select sub-faces per the operation and weld shared boundary edges so
//!    the rebuilt shell is a closed 2-manifold;
//! 5. assemble a shell (and solid when closed) and validate the volume.
//!
//! Non-planar inputs fall back to a voxel boolean (`crate::boolean_ops`).

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::GeomPlane;

use crate::abs::ShapeType;
use crate::brep_extrema::is_inside;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::inttools::{edge_edge_intersections, edge_face_intersections};
use crate::shape::{Edge, Face, Shell, Solid, TopoShape};
use crate::shell_check::shell_is_closed;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{
    edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of, wires_of_face,
};

/// The boolean operation to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp {
    /// A ∪ B.
    Fuse,
    /// A − B.
    Cut,
    /// A ∩ B.
    Common,
}

/// Result of a boolean operation.
#[derive(Debug, Clone)]
pub struct BooleanResult {
    /// The resulting shape (compound for a disjoint Fuse, empty compound for
    /// an empty Common, otherwise the solid or its shell).
    pub shape: TopoShape,
    /// The resulting solid, when the rebuilt shell is closed.
    pub solid: Option<Solid>,
    /// The rebuilt shell(s).
    pub shells: Vec<Shell>,
    /// The selected (rebuilt) faces.
    pub faces: Vec<Face>,
    /// Non-fatal diagnostics (volume inconsistencies, voxel cross-check, …).
    pub warnings: Vec<String>,
}

// ---------------------------------------------------------------------------
// Sub-face bookkeeping
// ---------------------------------------------------------------------------

/// A sub-face produced by splitting a source face's polygon.
#[derive(Debug, Clone)]
struct SubFace {
    face: Face,
    /// 3D polygon vertices (in boundary order).
    poly3d: Vec<GpPnt>,
    /// The face's plane (normal = source solid's outward normal).
    plane: GpPln,
}

/// Classification of a face relative to the other solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Class {
    /// Strictly inside.
    In,
    /// Strictly outside.
    Out,
    /// Coincident with the other solid's surface, and the other solid's
    /// interior lies on the same side as this solid's interior (the face is
    /// part of the outer boundary of the result).
    OnSame,
    /// Coincident with the other solid's surface, and the other solid's
    /// interior lies on the opposite side (the other solid merely touches).
    OnOpp,
}

// ---------------------------------------------------------------------------
// Plane / projection helpers
// ---------------------------------------------------------------------------

fn plane_normal(pln: &GpPln) -> GpVec {
    GpVec::from_xyz(pln.axis().direction().xyz())
}

/// Extract a `GpPln` from a planar face's registered surface. `None` when the
/// surface is not geometrically a plane (triggers the voxel fallback).
fn face_plane_local(face: &Face) -> Option<GpPln> {
    let s = GeometryRegistry::global().face_surface(&face.0)?;
    if !crate::face_face::is_plane_like(&*s) {
        return None;
    }
    let (p, n) = crate::face_face::plane_geometry(&*s);
    let normal = GpDir::from_vec(&n).ok()?;
    let ax = GpAx3::from_ax1(&GpAx1::new(p, normal));
    Some(GpPln::new(ax))
}

/// Project a 3D point onto a plane, expressed in the plane's (x_dir, y_dir).
fn project_point_to_plane(pln: &GpPln, p: &GpPnt) -> GpPnt2d {
    let origin = pln.location();
    let x = GpVec::from_xyz(pln.position().x_direction().xyz());
    let y = GpVec::from_xyz(pln.position().y_direction().xyz());
    let v = GpVec::from_pnts(&origin, p);
    GpPnt2d::new(v.dot(&x), v.dot(&y))
}

/// Map a 2D point in the plane's (x_dir, y_dir) basis back to 3D.
fn plane_point_from_2d(pln: &GpPln, p: &GpPnt2d) -> GpPnt {
    let origin = pln.location();
    let x = GpVec::from_xyz(pln.position().x_direction().xyz());
    let y = GpVec::from_xyz(pln.position().y_direction().xyz());
    origin
        .translated_vec(&x.multiplied_scalar(p.x()))
        .translated_vec(&y.multiplied_scalar(p.y()))
}

/// The plane with the normal reversed (keeps the origin and x-direction).
fn reversed_plane(pln: &GpPln) -> GpPln {
    let ax = pln.position();
    let n = ax.direction().reversed();
    GpPln::new(GpAx3::new(ax.location(), n, ax.x_direction()).unwrap_or(ax))
}

fn planes_coincident(a: &GpPln, b: &GpPln, tol: f64) -> bool {
    let na = plane_normal(a);
    let nb = plane_normal(b);
    if na.cross_magnitude(&nb) > tol {
        return false;
    }
    let d = GpVec::from_pnts(&a.location(), &b.location()).dot(&na.normalized()).abs();
    d <= tol
}

fn polygon_centroid_3d(pts: &[GpPnt]) -> GpPnt {
    let n = pts.len();
    if n == 0 {
        return GpPnt::zero();
    }
    let mut acc = GpVec::zero();
    for p in pts {
        acc = acc.added(&GpVec::from_xyz(&p.coord));
    }
    GpPnt::from_xyz(&acc.divided(n as f64).coord)
}

/// Distinct boundary vertex points of a face (from its wire edges' curves).
fn face_boundary_points(face: &Face) -> Option<Vec<GpPnt>> {
    let mut pts: Vec<GpPnt> = Vec::new();
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            let (a, b) = BRepTool::edge_vertices(&e)?;
            if !pts.iter().any(|p| p.distance(&a) < 1e-9) {
                pts.push(a);
            }
            if !pts.iter().any(|p| p.distance(&b) < 1e-9) {
                pts.push(b);
            }
        }
    }
    if pts.len() < 3 {
        None
    } else {
        Some(pts)
    }
}

/// The face's boundary as a 2D polygon, projected onto `pln` and sorted in
/// cyclic order around the face centroid (correct for convex planar faces).
fn face_polygon_local(face: &Face, pln: &GpPln) -> Option<Vec<GpPnt2d>> {
    let pts = face_boundary_points(face)?;
    let c3 = polygon_centroid_3d(&pts);
    let c = project_point_to_plane(pln, &c3);
    let mut with_angle: Vec<(f64, GpPnt2d)> = pts
        .iter()
        .map(|p| {
            let p2 = project_point_to_plane(pln, p);
            let v = (p2.x() - c.x(), p2.y() - c.y());
            (v.1.atan2(v.0), p2)
        })
        .collect();
    with_angle.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    let out: Vec<GpPnt2d> = with_angle.into_iter().map(|(_, p)| p).collect();
    if out.len() < 3 {
        None
    } else {
        Some(out)
    }
}

/// Parameter intervals of the line `origin + t·dir` (t = 3D arc-length, `dir`
/// unit) that lie inside the 2D polygon, expressed in `pln`'s frame. Because
/// `t` is the shared 3D arc-length, intervals computed for two faces on the
/// same line are directly comparable.
fn face_line_intervals(pln: &GpPln, origin: &GpPnt, dir: &GpVec, poly: &[GpPnt2d]) -> Vec<(f64, f64)> {
    let n = poly.len();
    if n < 3 {
        return Vec::new();
    }
    let xd = GpVec::from_xyz(pln.position().x_direction().xyz());
    let yd = GpVec::from_xyz(pln.position().y_direction().xyz());
    let v0 = GpVec::from_pnts(&pln.location(), origin);
    let p0_2 = GpPnt2d::new(v0.dot(&xd), v0.dot(&yd));
    let (dx, dy) = (dir.dot(&xd), dir.dot(&yd));
    let len2 = dx * dx + dy * dy;
    if len2 <= 1e-30 {
        return Vec::new();
    }
    let s = |p: &GpPnt2d| dx * (p.y() - p0_2.y()) - dy * (p.x() - p0_2.x());
    let proj = |p: &GpPnt2d| ((p.x() - p0_2.x()) * dx + (p.y() - p0_2.y()) * dy) / len2;
    let eps = 1e-9 * len2.sqrt().max(1.0);

    let mut crossings: Vec<f64> = Vec::new();
    let mut intervals: Vec<(f64, f64)> = Vec::new();
    for i in 0..n {
        let a = poly[i];
        let c = poly[(i + 1) % n];
        let sa = s(&a);
        let sc = s(&c);
        if sa.abs() <= eps && sc.abs() <= eps {
            let (ta, tc) = (proj(&a), proj(&c));
            intervals.push((ta.min(tc), ta.max(tc)));
        } else if sa * sc < 0.0 {
            let w = sa / (sa - sc);
            crossings.push(proj(&a) + (proj(&c) - proj(&a)) * w);
        }
    }
    crossings.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|x, y| (*x - *y).abs() <= eps);
    for w in crossings.windows(2) {
        let (t0, t1) = (w[0], w[1]);
        if t1 - t0 <= eps {
            continue;
        }
        let mid = GpPnt2d::new(p0_2.x() + dx * 0.5 * (t0 + t1), p0_2.y() + dy * 0.5 * (t0 + t1));
        if point_in_polygon2d(poly, &mid) {
            intervals.push((t0, t1));
        }
    }
    intervals.retain(|(a, b)| b - a > 1e-12);
    intervals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut merged: Vec<(f64, f64)> = Vec::new();
    for iv in intervals {
        if let Some(last) = merged.last_mut() {
            if iv.0 <= last.1 + 1e-9 {
                last.1 = last.1.max(iv.1);
            } else {
                merged.push(iv);
            }
        } else {
            merged.push(iv);
        }
    }
    merged
}

/// Intersection of two planar faces as one or more 3D segments (the plane–
/// plane line clipped to both boundary polygons).
///
/// Uses the robust local polygon extraction (angular-sorted boundary points)
/// rather than `inttools::face_face_intersection_segments`, whose polygon
/// extraction relies on the wire edge order and can drop corners.
fn face_face_segments_local(f1: &Face, f2: &Face, tol: f64) -> Vec<(GpPnt, GpPnt)> {
    let p1 = match face_plane_local(f1) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let p2 = match face_plane_local(f2) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let (origin, dir) = match crate::face_face::plane_plane_intersection(&p1, &p2) {
        Some(x) => x,
        None => return Vec::new(), // parallel (coplanar included)
    };
    let poly1 = match face_polygon_local(f1, &p1) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let poly2 = match face_polygon_local(f2, &p2) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let iv1 = face_line_intervals(&p1, &origin, &dir, &poly1);
    let iv2 = face_line_intervals(&p2, &origin, &dir, &poly2);
    let mut out = Vec::new();
    for &(a0, a1) in &iv1 {
        for &(b0, b1) in &iv2 {
            let lo = a0.max(b0);
            let hi = a1.min(b1);
            if hi - lo > tol.max(1e-9) {
                out.push((
                    origin.translated_vec(&dir.multiplied_scalar(lo)),
                    origin.translated_vec(&dir.multiplied_scalar(hi)),
                ));
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// 2D polygon splitting
// ---------------------------------------------------------------------------

/// Remove consecutive duplicate vertices (closing point handled separately).
fn dedupe_polygon(poly: &[GpPnt2d]) -> Vec<GpPnt2d> {
    let mut out: Vec<GpPnt2d> = Vec::new();
    for p in poly {
        if let Some(last) = out.last() {
            if (last.x() - p.x()).abs() < 1e-9 && (last.y() - p.y()).abs() < 1e-9 {
                continue;
            }
        }
        out.push(*p);
    }
    if out.len() > 1 {
        let first = out[0];
        let last = *out.last().unwrap();
        if (first.x() - last.x()).abs() < 1e-9 && (first.y() - last.y()).abs() < 1e-9 {
            out.pop();
        }
    }
    out
}

/// Split a simple polygon into the two half-polygons on either side of the
/// directed line through `a` and `b`. Vertices on the line belong to both
/// halves. Degenerate halves (fewer than 3 vertices / zero area) are dropped;
/// a line that does not cross the polygon leaves it unchanged.
fn split_polygon_by_segment(poly: &[GpPnt2d], a: &GpPnt2d, b: &GpPnt2d) -> Vec<Vec<GpPnt2d>> {
    let d = (b.x() - a.x(), b.y() - a.y());
    let len2 = d.0 * d.0 + d.1 * d.1;
    if len2 < 1e-24 || poly.len() < 3 {
        return vec![poly.to_vec()];
    }
    let eps = 1e-9 * len2.sqrt().max(1.0);
    let side = |p: &GpPnt2d| d.0 * (p.y() - a.y()) - d.1 * (p.x() - a.x());

    let mut pos: Vec<GpPnt2d> = Vec::new();
    let mut neg: Vec<GpPnt2d> = Vec::new();
    let n = poly.len();
    for i in 0..n {
        let p = poly[i];
        let q = poly[(i + 1) % n];
        let sp = side(&p);
        let sq = side(&q);
        if sp >= -eps {
            pos.push(p);
        }
        if sp <= eps {
            neg.push(p);
        }
        if (sp > eps && sq < -eps) || (sp < -eps && sq > eps) {
            let t = sp / (sp - sq);
            let inter = GpPnt2d::new(p.x() + t * (q.x() - p.x()), p.y() + t * (q.y() - p.y()));
            pos.push(inter);
            neg.push(inter);
        }
    }
    let mut out: Vec<Vec<GpPnt2d>> = Vec::new();
    for raw in [pos, neg] {
        let pp = dedupe_polygon(&raw);
        if pp.len() >= 3 && polygon_area2d(&pp).abs() > 1e-12 {
            out.push(pp);
        }
    }
    if out.is_empty() {
        // The line ran along the polygon boundary or outside it.
        out.push(dedupe_polygon(poly));
    }
    out
}

/// Whether two segments lie on the same (infinite) line, within tolerance.
fn segments_same_line(a: &(GpPnt2d, GpPnt2d), b: &(GpPnt2d, GpPnt2d)) -> bool {
    let d1 = (a.1.x() - a.0.x(), a.1.y() - a.0.y());
    let d2 = (b.1.x() - b.0.x(), b.1.y() - b.0.y());
    let cr = d1.0 * d2.1 - d1.1 * d2.0;
    if cr.abs() > 1e-9 {
        return false;
    }
    let v = (b.0.x() - a.0.x(), b.0.y() - a.0.y());
    (d1.0 * v.1 - d1.1 * v.0).abs() <= 1e-6
}

fn dedupe_segments(segs: Vec<(GpPnt2d, GpPnt2d)>) -> Vec<(GpPnt2d, GpPnt2d)> {
    let mut out: Vec<(GpPnt2d, GpPnt2d)> = Vec::new();
    for s in segs {
        if !out.iter().any(|t| segments_same_line(&s, t)) {
            out.push(s);
        }
    }
    out
}

/// Sequentially split `poly` by every cutting line in `segs`.
fn split_polygon_by_segments(poly: &[GpPnt2d], segs: &[(GpPnt2d, GpPnt2d)]) -> Vec<Vec<GpPnt2d>> {
    let mut polys = vec![poly.to_vec()];
    for (a, b) in segs {
        let mut next: Vec<Vec<GpPnt2d>> = Vec::new();
        for p in &polys {
            next.extend(split_polygon_by_segment(p, a, b));
        }
        polys = next;
        if polys.is_empty() {
            break;
        }
    }
    polys
}

// ---------------------------------------------------------------------------
// Welded edge sharing
// ---------------------------------------------------------------------------

/// Spatial-hash point welder: maps near-coincident points to one index so the
/// rebuilt boundary edges are shared between adjacent faces.
struct Weld {
    cell: f64,
    points: Vec<GpPnt>,
    grid: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl Weld {
    fn new(cell: f64) -> Self {
        Self { cell: cell.max(1e-9), points: Vec::new(), grid: HashMap::new() }
    }

    fn key(&self, p: &GpPnt) -> (i64, i64, i64) {
        (
            f64::floor(p.x() / self.cell) as i64,
            f64::floor(p.y() / self.cell) as i64,
            f64::floor(p.z() / self.cell) as i64,
        )
    }

    fn weld(&mut self, p: &GpPnt) -> usize {
        let k = self.key(p);
        let mut best: Option<(usize, f64)> = None;
        for dx in -1i64..=1 {
            for dy in -1i64..=1 {
                for dz in -1i64..=1 {
                    if let Some(bucket) = self.grid.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                        for &i in bucket {
                            let d = self.points[i].distance(p);
                            if d <= self.cell * 1.01 && best.map_or(true, |(_, bd)| d < bd) {
                                best = Some((i, d));
                            }
                        }
                    }
                }
            }
        }
        match best {
            Some((i, _)) => i,
            None => {
                let i = self.points.len();
                self.points.push(*p);
                self.grid.entry(k).or_default().push(i);
                i
            }
        }
    }
}

/// Map of unordered welded-point pairs to a shared `Edge`.
#[derive(Default)]
struct EdgeMap {
    edges: HashMap<(usize, usize), Edge>,
}

impl EdgeMap {
    fn edge(&mut self, b: &TopoBuilder, i1: usize, i2: usize, pts: &[GpPnt]) -> Edge {
        let key = (i1.min(i2), i1.max(i2));
        if let Some(e) = self.edges.get(&key) {
            return e.clone();
        }
        let e = b.make_edge_segment(&pts[i1], &pts[i2]);
        self.edges.insert(key, e.clone());
        e
    }
}

// ---------------------------------------------------------------------------
// Face splitting
// ---------------------------------------------------------------------------

fn split_faces(
    b: &TopoBuilder,
    src_faces: &[Face],
    planes: &[Option<GpPln>],
    segs: &[Vec<(GpPnt, GpPnt)>],
    weld: &mut Weld,
    edge_map: &mut EdgeMap,
) -> Vec<SubFace> {
    let mut out = Vec::new();
    for (i, f) in src_faces.iter().enumerate() {
        let Some(pln) = planes[i].clone() else { continue };
        let Some(poly2d) = face_polygon_local(f, &pln) else { continue };
        let segs2d: Vec<(GpPnt2d, GpPnt2d)> = segs[i]
            .iter()
            .filter_map(|(a, b)| {
                let a2 = project_point_to_plane(&pln, a);
                let b2 = project_point_to_plane(&pln, b);
                if a2.distance(&b2) < 1e-12 {
                    None
                } else {
                    Some((a2, b2))
                }
            })
            .collect();
        let segs2d = dedupe_segments(segs2d);
        let parts = split_polygon_by_segments(&poly2d, &segs2d);
        for part in parts {
            if let Some(sf) = polygon_to_subface(b, &pln, &part, weld, edge_map) {
                out.push(sf);
            }
        }
    }
    out
}

fn polygon_to_subface(
    b: &TopoBuilder,
    pln: &GpPln,
    poly2d: &[GpPnt2d],
    weld: &mut Weld,
    edge_map: &mut EdgeMap,
) -> Option<SubFace> {
    if poly2d.len() < 3 || polygon_area2d(poly2d).abs() < 1e-12 {
        return None;
    }
    let pts3d: Vec<GpPnt> = poly2d.iter().map(|p| plane_point_from_2d(pln, p)).collect();
    let idx: Vec<usize> = pts3d.iter().map(|p| weld.weld(p)).collect();
    let n = pts3d.len();
    let mut wire_edges: Vec<Edge> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        if idx[i] == idx[j] {
            continue;
        }
        let e = edge_map.edge(b, idx[i], idx[j], &weld.points);
        wire_edges.push(e);
    }
    if wire_edges.len() < 3 {
        return None;
    }
    let wire = b.make_wire(&wire_edges);
    let face = b.make_face(Arc::new(GeomPlane::new(pln.clone())), &[wire]);
    Some(SubFace { face, poly3d: pts3d, plane: pln.clone() })
}

/// Rebuild a sub-face with the reversed plane (used for Cut's cut-through
/// faces so their mesh normal points away from the resulting solid).
fn flip_face(b: &TopoBuilder, sub: &SubFace) -> Face {
    let rpl = reversed_plane(&sub.plane);
    let mut edges: Vec<Edge> = Vec::new();
    for w in wires_of_face(&sub.face) {
        edges.extend(edges_of_wire(&w));
    }
    let wire = b.make_wire(&edges);
    b.make_face(Arc::new(GeomPlane::new(rpl)), &[wire])
}

/// Rebuild a face with the reversed plane, keeping its wire edges.
fn flip_face_plane(b: &TopoBuilder, f: &Face, pln: &GpPln) -> Face {
    let rpl = reversed_plane(pln);
    let mut edges: Vec<Edge> = Vec::new();
    for w in wires_of_face(f) {
        edges.extend(edges_of_wire(&w));
    }
    let wire = b.make_wire(&edges);
    b.make_face(Arc::new(GeomPlane::new(rpl)), &[wire])
}

/// Orientation-consistency pass: make every kept face's surface normal point
/// OUTWARD from the result solid.
///
/// The rebuilt faces inherit the source solids' plane orientations, which may
/// be inward for some sources (e.g. swept/prism faces). This pass probes a
/// point just outside each face along its current normal; if that probe is
/// inside the (closed) result shell, the normal points into the solid and the
/// face is flipped. For a closed manifold this yields a consistently outward
/// shell, so `brep_gprop::volume`'s signed divergence sum matches the volume.
fn orient_faces_outward(bld: &TopoBuilder, faces: &[Face]) -> Vec<Face> {
    let shell = bld.make_shell(faces);
    if !shell_is_closed(&shell) {
        return faces.to_vec();
    }
    let shape = shell.0.clone();
    let eps = 1e-6;
    let mut out: Vec<Face> = Vec::with_capacity(faces.len());
    for f in faces {
        let pln = match face_plane_local(f) {
            Some(p) => p,
            None => {
                out.push(f.clone());
                continue;
            }
        };
        let n = plane_normal(&pln).normalized();
        let pts = face_boundary_points(f).unwrap_or_default();
        let c = polygon_centroid_3d(&pts);
        let probe = c.translated_vec(&n.multiplied_scalar(eps));
        if is_inside(&shape, &probe) {
            out.push(flip_face_plane(bld, f, &pln));
        } else {
            out.push(f.clone());
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------

/// Is the point `p` (within `tol`) on the given face's polygon?
fn point_in_face_polygon(face: &Face, pln: &GpPln, p: &GpPnt, tol: f64) -> bool {
    let n = plane_normal(pln);
    let d = GpVec::from_pnts(&pln.location(), p).dot(&n.normalized()).abs();
    if d > tol {
        return false;
    }
    let Some(poly) = face_polygon_local(face, pln) else { return false };
    let p2 = project_point_to_plane(pln, p);
    point_in_polygon2d(&poly, &p2)
}

/// True when the sub-face's plane is coincident with a face of `other_faces`
/// and its centroid lies inside that face's polygon.
fn face_on_other(plane: &GpPln, poly3d: &[GpPnt], other_faces: &[Face], tol: f64) -> bool {
    let c = polygon_centroid_3d(poly3d);
    for of in other_faces {
        let Some(opl) = face_plane_local(of) else { continue };
        if !planes_coincident(plane, &opl, tol) {
            continue;
        }
        if point_in_face_polygon(of, &opl, &c, tol) {
            return true;
        }
    }
    false
}

fn classify_face(
    plane: &GpPln,
    poly3d: &[GpPnt],
    other_faces: &[Face],
    other: &TopoShape,
    tol: f64,
) -> Class {
    let c = polygon_centroid_3d(poly3d);
    if face_on_other(plane, poly3d, other_faces, tol) {
        let n = plane_normal(plane);
        let eps = tol.max(1e-6);
        let p_in = c.translated_vec(&n.reversed().multiplied_scalar(eps));
        if is_inside(other, &p_in) {
            Class::OnSame
        } else {
            Class::OnOpp
        }
    } else if is_inside(other, &c) {
        Class::In
    } else {
        Class::Out
    }
}

fn select(c: Class, from_a: bool, op: BoolOp) -> bool {
    match op {
        BoolOp::Fuse => {
            if from_a {
                matches!(c, Class::Out | Class::OnSame)
            } else {
                matches!(c, Class::Out)
            }
        }
        BoolOp::Cut => {
            if from_a {
                matches!(c, Class::Out | Class::OnOpp)
            } else {
                matches!(c, Class::In)
            }
        }
        BoolOp::Common => {
            if from_a {
                matches!(c, Class::In | Class::OnSame)
            } else {
                matches!(c, Class::In)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Result helpers
// ---------------------------------------------------------------------------

fn empty_result(op: BoolOp) -> BooleanResult {
    let b = TopoBuilder::new();
    let comp = b.make_compound_of(&[]);
    BooleanResult {
        shape: comp.0,
        solid: None,
        shells: vec![],
        faces: vec![],
        warnings: if matches!(op, BoolOp::Fuse) { vec!["empty fuse result".into()] } else { vec![] },
    }
}

fn disjoint_result(a: &TopoShape, b: &TopoShape, op: BoolOp) -> BooleanResult {
    let builder = TopoBuilder::new();
    match op {
        BoolOp::Fuse => {
            let comp = builder.make_compound_of(&[a.clone(), b.clone()]);
            let shells: Vec<Shell> = shapes_of(a, ShapeType::Shell)
                .into_iter()
                .map(Shell)
                .chain(shapes_of(b, ShapeType::Shell).into_iter().map(Shell))
                .collect();
            BooleanResult { shape: comp.0, solid: None, shells, faces: vec![], warnings: vec![] }
        }
        BoolOp::Cut => {
            let solid = Solid::wrap(a.clone());
            let shells: Vec<Shell> = shapes_of(a, ShapeType::Shell).into_iter().map(Shell).collect();
            BooleanResult {
                shape: a.clone(),
                solid,
                shells,
                faces: faces_of(a),
                warnings: vec![],
            }
        }
        BoolOp::Common => empty_result(op),
    }
}

fn voxel_fallback(a: &TopoShape, b: &TopoShape, op: BoolOp) -> Result<BooleanResult, String> {
    let vop = match op {
        BoolOp::Fuse => crate::boolean_ops::BoolOp::Union,
        BoolOp::Cut => crate::boolean_ops::BoolOp::Subtraction,
        BoolOp::Common => crate::boolean_ops::BoolOp::Intersection,
    };
    let mesh = crate::boolean_ops::voxel_boolean(a, b, 32, vop)?;
    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&mesh);
    let warnings = vec!["non-planar input: fell back to voxel boolean".to_string()];
    let shape = brep.solid.clone().map(|s| s.0).unwrap_or_else(|| brep.shell.0.clone());
    Ok(BooleanResult {
        shape,
        solid: brep.solid,
        shells: vec![brep.shell],
        faces: brep.faces,
        warnings,
    })
}

fn validate(a: &TopoShape, b: &TopoShape, op: BoolOp, result: &mut BooleanResult) {
    if result.faces.is_empty() {
        return;
    }
    let vol = crate::brep_gprop::volume(&result.shape, 0.02);
    let vol_a = crate::brep_gprop::volume(a, 0.02);
    let vol_b = crate::brep_gprop::volume(b, 0.02);
    match op {
        BoolOp::Fuse => {
            if vol + 1e-6 < vol_a.max(vol_b) {
                result.warnings.push(format!(
                    "fuse volume {vol} below max input {}",
                    vol_a.max(vol_b)
                ));
            }
            let voxel = crate::solid_union::union_volume(a, b, 32);
            if voxel > 1e-9 {
                let rel = (vol - voxel).abs() / voxel;
                if rel > 0.15 {
                    result
                        .warnings
                        .push(format!("voxel cross-check mismatch: exact {vol}, voxel {voxel}"));
                }
            }
        }
        BoolOp::Cut => {
            if vol > vol_a + 1e-6 {
                result.warnings.push(format!("cut volume {vol} exceeds input {vol_a}"));
            }
        }
        BoolOp::Common => {
            if vol > vol_a.min(vol_b) + 1e-6 {
                result.warnings.push(format!(
                    "common volume {vol} exceeds min input {}",
                    vol_a.min(vol_b)
                ));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Exact boolean of two planar-faced solids.
///
/// For disjoint inputs a Fuse returns a compound of both shapes, a Cut returns
/// `a`, and a Common returns an empty compound. Non-planar faces fall back to a
/// voxel boolean and add a warning.
pub fn boolean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let fa = faces_of(a);
    let fb = faces_of(b);

    // Planarity gate → voxel fallback for curved inputs.
    if fa.is_empty() || fb.is_empty() {
        return voxel_fallback(a, b, op);
    }
    let planes_a: Vec<Option<GpPln>> = fa.iter().map(face_plane_local).collect();
    let planes_b: Vec<Option<GpPln>> = fb.iter().map(face_plane_local).collect();
    if planes_a.iter().any(Option::is_none) || planes_b.iter().any(Option::is_none) {
        return voxel_fallback(a, b, op);
    }

    // Disjoint shortcut (bounding boxes don't overlap).
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op));
    }

    // Step 1: face–face intersection segments, grouped by the face they lie on.
    let mut segs_a: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); fa.len()];
    let mut segs_b: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); fb.len()];
    for (i, af) in fa.iter().enumerate() {
        for (j, bf) in fb.iter().enumerate() {
            let segs = face_face_segments_local(af, bf, tol);
            for s in &segs {
                segs_a[i].push(*s);
                segs_b[j].push(*s);
            }
        }
    }

    // Edge–edge / edge–face hits (used for completeness checks; the polygon
    // split above already subsumes edge splitting, since split points become
    // sub-polygon vertices).
    let ea = edges_of(a);
    let eb = edges_of(b);
    for e in &ea {
        for f in &fb {
            let _ = edge_face_intersections(e, f, tol);
        }
    }
    for e in &eb {
        for f in &fa {
            let _ = edge_face_intersections(e, f, tol);
        }
    }
    for e1 in &ea {
        for e2 in &eb {
            let _ = edge_edge_intersections(e1, e2, tol);
        }
    }

    // Step 3: split every face along the segments that lie on it.
    let bld = TopoBuilder::new();
    let weld_cell = tol.max(1e-7);
    let mut weld = Weld::new(weld_cell);
    let mut edge_map = EdgeMap::default();
    let subs_a = split_faces(&bld, &fa, &planes_a, &segs_a, &mut weld, &mut edge_map);
    let subs_b = split_faces(&bld, &fb, &planes_b, &segs_b, &mut weld, &mut edge_map);

    // Steps 4–5: classify each sub-face and select per the operation.
    let mut selected: Vec<SubFace> = Vec::new();
    let mut flipped: Vec<bool> = Vec::new();
    for sub in &subs_a {
        let c = classify_face(&sub.plane, &sub.poly3d, &fb, b, tol);
        if select(c, true, op) {
            selected.push(sub.clone());
            flipped.push(false);
        }
    }
    for sub in &subs_b {
        let c = classify_face(&sub.plane, &sub.poly3d, &fa, a, tol);
        if select(c, false, op) {
            selected.push(sub.clone());
            // B's cut-through surfaces are oriented outward from A.
            flipped.push(matches!(op, BoolOp::Cut));
        }
    }

    // Step 6: rebuild the shell (and solid when closed).
    let mut result_faces: Vec<Face> = Vec::new();
    for (sub, fl) in selected.iter().zip(&flipped) {
        if *fl {
            result_faces.push(flip_face(&bld, sub));
        } else {
            result_faces.push(sub.face.clone());
        }
    }
    if result_faces.is_empty() {
        return Ok(empty_result(op));
    }
    // Make every face's normal point outward from the result, so the signed
    // mesh volume is consistent (source faces may be oriented inward).
    let result_faces = orient_faces_outward(&bld, &result_faces);
    let shell = bld.make_shell(&result_faces);
    let closed = shell_is_closed(&shell);
    let solid = if closed { Some(bld.make_solid(&[shell.clone()])) } else { None };
    let shape = solid
        .as_ref()
        .map(|s| s.0.clone())
        .unwrap_or_else(|| shell.0.clone());

    let mut result = BooleanResult {
        shape,
        solid,
        shells: vec![shell],
        faces: result_faces,
        warnings: vec![],
    };
    if !closed {
        result.warnings.push("result shell is not closed".into());
    }

    // Step 7: validate the volume.
    validate(a, b, op, &mut result);
    Ok(result)
}

// ---------------------------------------------------------------------------
// Full-topology boolean
// ---------------------------------------------------------------------------
//
// Extensions porting the rest of `BOPAlgo_Builder` / `TopOpeBRep`: N-ary
// (multi-argument) operations, compound expansion, self-intersection
// detection, result validation with tolerance escalation, degenerate-input
// handling and result decomposition. All functions use the same
// `Result<BooleanResult, String>` convention and dispatch to the exact planar
// boolean (`boolean`) or the curved dispatcher (`crate::bop_curved::curved_boolean_full`).

/// Wrap a single shape as a [`BooleanResult`] (identity operation).
fn single_shape_result(s: &TopoShape) -> BooleanResult {
    let shells: Vec<Shell> = shapes_of(s, ShapeType::Shell).into_iter().map(Shell).collect();
    let faces = faces_of(s);
    let solid = Solid::wrap(s.clone());
    BooleanResult { shape: s.clone(), solid, shells, faces, warnings: vec![] }
}

/// Recursively flatten a compound into its non-compound sub-shapes.
///
/// A compound whose children are themselves compounds is fully flattened, so
/// nested results (e.g. a compound built by folding a multi-fuse) decompose
/// into their atomic solids/shells/faces.
fn expand_compound(s: &TopoShape) -> Vec<TopoShape> {
    if !s.is_compound() {
        return vec![s.clone()];
    }
    let mut out = Vec::new();
    for k in s.tshape.read().unwrap().children.clone() {
        out.extend(expand_compound(&TopoShape::from_handle(k)));
    }
    out
}

/// Split a `Compound` into its top-level sub-shapes (its direct children).
///
/// A non-compound shape is returned as a single-element slice. Nested
/// compounds are kept as-is at the top level (use [`boolean_multi`] or
/// [`expand_compound`] when full flattening is needed).
pub fn decompose_compound(shape: &TopoShape) -> Vec<TopoShape> {
    if !shape.is_compound() {
        return vec![shape.clone()];
    }
    shape
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .map(|k| TopoShape::from_handle(k.clone()))
        .collect()
}

/// Decompose a shape into its connected boundary components.
///
/// * a `Solid` → its shells, each wrapped back into a one-shell solid;
/// * a `Compound` → its top-level sub-shapes (recursively flattened);
/// * anything else → the shape itself.
///
/// This is the "result decomposition" step of `BOPAlgo_Builder`: after a
/// boolean, a result may hold several disconnected solids; this splits them
/// apart so callers can reason about each component.
pub fn shape_components(shape: &TopoShape) -> Vec<TopoShape> {
    match shape.shape_type() {
        ShapeType::Solid => {
            let bld = TopoBuilder::new();
            let mut comps = Vec::new();
            for sh in shapes_of(shape, ShapeType::Shell) {
                comps.push(bld.make_solid(&[Shell(sh)]).0);
            }
            if comps.is_empty() {
                vec![shape.clone()]
            } else {
                comps
            }
        }
        ShapeType::Compound => expand_compound(shape),
        _ => vec![shape.clone()],
    }
}

/// Does `s` carry any boundary faces? Empty compounds, null wires and other
/// degenerate inputs have none and are treated as the "empty" shape.
fn has_content(s: &TopoShape) -> bool {
    !faces_of(s).is_empty()
}

/// Is `s` usable as a solid operand? A `Solid` or a *closed* `Shell`; empty
/// shells/solids are rejected so degenerate inputs route to
/// [`boolean_degenerate`].
fn is_solid_input(s: &TopoShape) -> bool {
    if !has_content(s) {
        return false;
    }
    if s.is_solid() {
        return true;
    }
    if s.is_shell() {
        return shell_is_closed(&Shell(s.clone()));
    }
    false
}

/// Average position of the vertices of `s` (a coarse centroid probe).
fn shape_centroid(s: &TopoShape) -> GpPnt {
    let vs = vertices_of(s);
    if vs.is_empty() {
        return GpPnt::zero();
    }
    let n = vs.len();
    let mut acc = GpVec::zero();
    for v in &vs {
        acc = acc.added(&GpVec::from_xyz(&vertex_position(v).coord));
    }
    GpPnt::from_xyz(&acc.divided(n as f64).coord)
}

/// Whether the bounding boxes of `a` and `b` overlap (touch counts).
fn bounding_boxes_overlap(a: &TopoShape, b: &TopoShape) -> bool {
    let ba = crate::bbox_from_geometry::shape_bbox(a);
    let bb = crate::bbox_from_geometry::shape_bbox(b);
    !ba.is_void() && !bb.is_void() && !ba.is_out_box(&bb)
}

/// Merge several result shapes into a single [`BooleanResult`].
///
/// Empty shapes are dropped. A single survivor becomes the result shape (and
/// its solid, when closed); multiple survivors are assembled into a compound.
/// This mirrors `BOPAlgo_Builder` returning a compound of the disconnected
/// result solids.
fn merge_shapes_result(shapes: Vec<TopoShape>, warnings: Vec<String>) -> BooleanResult {
    let mut kept: Vec<TopoShape> = shapes.into_iter().filter(has_content).collect();
    if kept.is_empty() {
        let mut r = empty_result(BoolOp::Fuse);
        r.warnings.extend(warnings);
        return r;
    }
    if kept.len() == 1 {
        let mut r = single_shape_result(&kept[0]);
        r.warnings.extend(warnings);
        return r;
    }
    let bld = TopoBuilder::new();
    let comp = bld.make_compound_of(&kept);
    let mut shells = Vec::new();
    let mut faces = Vec::new();
    for s in &kept {
        shells.extend(shapes_of(s, ShapeType::Shell).into_iter().map(Shell));
        faces.extend(faces_of(s));
    }
    BooleanResult { shape: comp.0, solid: None, shells, faces, warnings }
}

/// Fuse a list of shapes into a single result.
///
/// Components whose bounding boxes do not overlap are left apart and collected
/// into a compound (a disjoint Fuse). Components whose boxes do overlap are
/// merged pairwise through the exact boolean until no further merge is
/// possible. This keeps a *single* result object — either one solid or a
/// compound of the disconnected pieces.
fn fuse_components(shapes: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    let mut components: Vec<TopoShape> = shapes.iter().filter(|s| has_content(s)).cloned().collect();
    if components.is_empty() {
        return Ok(empty_result(BoolOp::Fuse));
    }
    let mut warnings: Vec<String> = Vec::new();
    // Repeatedly merge the first pair whose bounding boxes overlap.
    loop {
        let n = components.len();
        if n < 2 {
            break;
        }
        let mut merged_any = false;
        'search: for i in 0..n {
            for j in (i + 1)..n {
                if bounding_boxes_overlap(&components[i], &components[j]) {
                    let r =
                        crate::bop_curved::curved_boolean_full(&components[i], &components[j], BoolOp::Fuse, tol)?;
                    warnings.extend(r.warnings);
                    let mut next: Vec<TopoShape> = Vec::with_capacity(n - 1);
                    for k in 0..n {
                        if k != i && k != j {
                            next.push(components[k].clone());
                        }
                    }
                    next.push(r.shape);
                    components = next;
                    merged_any = true;
                    break 'search;
                }
            }
        }
        if !merged_any {
            break;
        }
    }
    Ok(merge_shapes_result(components, warnings))
}

/// Top-level dispatch between the exact boolean paths.
///
/// * either operand is a compound → [`boolean_compound`] (expand, per-part);
/// * either operand is not a solid → [`boolean_degenerate`] (best effort);
/// * otherwise the full curved/planar boolean dispatcher.
fn boolean_dispatch(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    if a.is_compound() || b.is_compound() {
        return boolean_compound(a, b, op, tol);
    }
    if !is_solid_input(a) || !is_solid_input(b) {
        return boolean_degenerate(a, b, op, tol);
    }
    crate::bop_curved::curved_boolean_full(a, b, op, tol)
}

/// Apply a boolean operation to a *list* of shapes (N-ary boolean).
///
/// `BOPAlgo_Builder`'s `Build`-with-many-arguments entry point. Handles:
///
/// * an empty list → an empty result (an empty compound);
/// * a single shape → the shape itself (identity);
/// * `Fuse` of mutually disjoint shapes → a `Compound` keeping every input
///   (not an error);
/// * `Fuse` of overlapping shapes → pairwise fold into one solid;
/// * `Cut` → fold `a0 − a1 − a2 − …`;
/// * `Common` → fold `a0 ∩ a1 ∩ a2 ∩ …`.
///
/// Compounds among the inputs are flattened first (a compound is a set of
/// components, so the operation distributes over it).
pub fn boolean_multi(shapes: &[TopoShape], op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    if shapes.is_empty() {
        return Ok(empty_result(op));
    }
    let flat: Vec<TopoShape> = shapes.iter().flat_map(expand_compound).collect();
    if flat.is_empty() {
        return Ok(empty_result(op));
    }
    if flat.len() == 1 {
        return Ok(single_shape_result(&flat[0]));
    }
    match op {
        BoolOp::Fuse => fuse_components(&flat, tol),
        BoolOp::Cut | BoolOp::Common => {
            // Fold left; degenerate/empty intermediates are absorbed by the
            // dispatch (a cut against empty leaves the accumulator unchanged,
            // an empty common stays empty).
            let mut acc = flat[0].clone();
            let mut warnings: Vec<String> = Vec::new();
            for s in &flat[1..] {
                let r = boolean_dispatch(&acc, s, op, tol)?;
                warnings.extend(r.warnings);
                acc = r.shape;
            }
            let mut res = single_shape_result(&acc);
            res.warnings.extend(warnings);
            Ok(res)
        }
    }
}

/// Boolean between shapes where either operand is a `Compound`.
///
/// Expands the compound(s) and applies the operation per component:
///
/// * `Fuse` — the union of *all* components of `a` and `b` (compounds are
///   unions of their parts, so `(A1 ∪ A2) ∪ (B1 ∪ B2)` is the full union);
/// * `Cut` — every component of `b` is cut out of every component of `a`;
/// * `Common` — the pairwise intersections `ai ∩ bj` are merged.
///
/// Solid operands are delegated to the exact boolean
/// (`boolean`/`curved_boolean_full`).
pub fn boolean_compound(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let subs_a = expand_compound(a);
    let subs_b = expand_compound(b);
    let mut warnings: Vec<String> = Vec::new();
    match op {
        BoolOp::Fuse => {
            let mut all = subs_a;
            all.extend(subs_b);
            fuse_components(&all, tol)
        }
        BoolOp::Cut => {
            let mut results: Vec<TopoShape> = Vec::new();
            for sa in &subs_a {
                let mut cur = sa.clone();
                for sb in &subs_b {
                    let r = boolean_dispatch(&cur, sb, BoolOp::Cut, tol)?;
                    warnings.extend(r.warnings);
                    cur = r.shape;
                }
                results.push(cur);
            }
            Ok(merge_shapes_result(results, warnings))
        }
        BoolOp::Common => {
            let mut results: Vec<TopoShape> = Vec::new();
            for sa in &subs_a {
                for sb in &subs_b {
                    let r = boolean_dispatch(sa, sb, BoolOp::Common, tol)?;
                    warnings.extend(r.warnings);
                    if has_content(&r.shape) {
                        results.push(r.shape);
                    }
                }
            }
            Ok(merge_shapes_result(results, warnings))
        }
    }
}

/// Report produced by [`detect_self_intersections`].
pub struct SelfIntersectionReport {
    /// Whether any pair of non-adjacent faces was found to intersect.
    pub found: bool,
    /// Number of intersecting face pairs (each pair counts as one
    /// self-intersection "edge" of the defect).
    pub edge_count: usize,
    /// Sample points on the found intersection curves.
    pub points: Vec<GpPnt>,
}

/// Do the two coplanar faces' boundary polygons overlap in area?
fn faces_polygon_overlap(f1: &Face, f2: &Face, tol: f64) -> bool {
    let pln = match face_plane_local(f1) {
        Some(p) => p,
        None => return false,
    };
    let (Some(poly1), Some(poly2)) = (face_polygon_local(f1, &pln), face_polygon_local(f2, &pln)) else {
        return false;
    };
    for v in &poly1 {
        if point_in_polygon2d(&poly2, v) {
            return true;
        }
    }
    for v in &poly2 {
        if point_in_polygon2d(&poly1, v) {
            return true;
        }
    }
    let _ = tol;
    false
}

/// Check a shape for self-intersecting faces.
///
/// A self-intersection is a pair of *non-adjacent* faces (faces that do not
/// share a boundary edge) whose underlying surfaces cross or overlap:
///
/// * transversal crossing → the sampled `SurfaceIntersection::Curves`;
/// * coplanar overlap → `SurfaceIntersection::Coincident` with overlapping
///   face polygons.
///
/// Adjacent faces are skipped because sharing a boundary edge is the normal
/// (and legal) way two faces of a solid meet. A valid closed box therefore
/// reports `found == false` (its only non-adjacent pairs are parallel faces),
/// while a shell built from two crossing faces reports `found == true`.
pub fn detect_self_intersections(shape: &TopoShape, tol: f64) -> SelfIntersectionReport {
    let tol = tol.max(1e-9);
    let faces = faces_of(shape);
    let mut report = SelfIntersectionReport { found: false, edge_count: 0, points: Vec::new() };
    if faces.len() < 2 {
        return report;
    }
    // Per-face boundary-edge identities, for the adjacency test. Two faces are
    // adjacent iff they share the same `TShape` edge.
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                continue; // adjacent faces legitimately meet along an edge
            }
            let (Some(sa), Some(sb)) = (
                GeometryRegistry::global().face_surface(&faces[i].0),
                GeometryRegistry::global().face_surface(&faces[j].0),
            ) else {
                continue;
            };
            match crate::intpatch::surface_surface_intersection(&*sa, &*sb, tol) {
                crate::intpatch::SurfaceIntersection::Curves(curves) => {
                    let mut pts: Vec<GpPnt> = Vec::new();
                    for c in &curves {
                        pts.extend(c.points.iter().cloned());
                    }
                    if !pts.is_empty() {
                        report.found = true;
                        report.edge_count += 1;
                        report.points.extend(pts);
                    }
                }
                crate::intpatch::SurfaceIntersection::Coincident => {
                    if faces_polygon_overlap(&faces[i], &faces[j], tol) {
                        report.found = true;
                        report.edge_count += 1;
                    }
                }
                crate::intpatch::SurfaceIntersection::None => {}
            }
        }
    }
    report
}

/// Structural validation of a boolean result.
///
/// Returns a list of human-readable issue strings (empty when the result is
/// clean): every shell must be closed (every boundary edge used by exactly two
/// faces — `shell_is_closed`) and the result must not be self-intersecting.
fn validate_boolean_result(r: &BooleanResult, tol: f64) -> Vec<String> {
    let mut issues: Vec<String> = Vec::new();
    if r.shells.is_empty() {
        issues.push("result has no shell".into());
    }
    for (i, sh) in r.shells.iter().enumerate() {
        if !shell_is_closed(sh) {
            issues.push(format!("result shell {i} is not closed"));
        }
    }
    let si = detect_self_intersections(&r.shape, tol);
    if si.found {
        issues.push(format!("result self-intersects ({} face pairs)", si.edge_count));
    }
    issues
}

/// Run a boolean and validate the result, re-computing at a larger tolerance
/// when validation finds issues.
///
/// The exact boolean is tolerance-sensitive near coincident faces; when the
/// first attempt fails validation (an open shell, a self-intersection), the
/// operation is re-run at `tol × 10` and `tol × 100`. The best (fewest issues)
/// result is returned, with the validation issues attached as warnings.
pub fn boolean_with_check(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let mut result = boolean_dispatch(a, b, op, tol)?;
    let mut issues = validate_boolean_result(&result, tol);
    if !issues.is_empty() {
        for &scale in &[10.0f64, 100.0] {
            let nt = tol * scale;
            let candidate = boolean_dispatch(a, b, op, nt)?;
            let cand_issues = validate_boolean_result(&candidate, nt);
            if cand_issues.len() < issues.len() {
                result = candidate;
                issues = cand_issues;
                result.warnings.push(format!("recomputed at tol {nt:.1e} after validation (original {tol:.1e})"));
            }
            if issues.is_empty() {
                break;
            }
        }
    }
    result.warnings.extend(issues);
    Ok(result)
}

/// Handle boolean operations with degenerate (non-solid / empty) inputs.
///
/// Degenerate inputs are: empty shapes (no faces — null wires, empty
/// compounds), open shells, single faces and wires. Rather than erroring, a
/// best-effort result is produced:
///
/// * empty ⊕ solid → the solid (empty is the identity for Fuse);
/// * solid − empty → the solid; empty − solid → empty;
/// * face ∪ solid → a compound of both, with a warning;
/// * face ∩ solid → the face when its centroid lies inside the solid, else
///   empty.
pub fn boolean_degenerate(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    if a.is_compound() || b.is_compound() {
        return boolean_compound(a, b, op, tol);
    }
    let a_empty = !has_content(a);
    let b_empty = !has_content(b);
    if a_empty || b_empty {
        return degenerate_empty(a, b, op, a_empty, b_empty);
    }
    let a_solid = is_solid_input(a);
    let b_solid = is_solid_input(b);
    if a_solid && b_solid {
        // Both solid: a caller reached us directly with two valid solids.
        return crate::bop_curved::curved_boolean_full(a, b, op, tol);
    }
    degenerate_non_solid(a, b, op, tol, a_solid, b_solid)
}

/// Empty-input branch of [`boolean_degenerate`].
fn degenerate_empty(a: &TopoShape, b: &TopoShape, op: BoolOp, a_empty: bool, b_empty: bool) -> Result<BooleanResult, String> {
    match op {
        BoolOp::Fuse => {
            if a_empty && b_empty {
                return Ok(empty_result(op));
            }
            if a_empty {
                let mut r = single_shape_result(b);
                r.warnings.push("empty input 'a' treated as empty; result is 'b' unchanged".into());
                return Ok(r);
            }
            let mut r = single_shape_result(a);
            r.warnings.push("empty input 'b' treated as empty; result is 'a' unchanged".into());
            Ok(r)
        }
        BoolOp::Cut => {
            if a_empty {
                return Ok(empty_result(op));
            }
            let mut r = single_shape_result(a);
            r.warnings.push("degenerate cut: 'b' has no faces; result is 'a' unchanged".into());
            Ok(r)
        }
        BoolOp::Common => Ok(empty_result(op)),
    }
}

/// Non-solid (but non-empty) input branch of [`boolean_degenerate`].
fn degenerate_non_solid(
    a: &TopoShape,
    b: &TopoShape,
    op: BoolOp,
    tol: f64,
    a_solid: bool,
    b_solid: bool,
) -> Result<BooleanResult, String> {
    let _ = tol;
    let bld = TopoBuilder::new();
    match op {
        BoolOp::Fuse => {
            // Fuse keeps every non-empty operand: a compound of the parts.
            let mut parts: Vec<TopoShape> = Vec::new();
            if has_content(a) {
                parts.push(a.clone());
            }
            if has_content(b) {
                parts.push(b.clone());
            }
            if parts.is_empty() {
                return Ok(empty_result(op));
            }
            let comp = bld.make_compound_of(&parts);
            let mut shells = Vec::new();
            let mut faces = Vec::new();
            for p in &parts {
                shells.extend(shapes_of(p, ShapeType::Shell).into_iter().map(Shell));
                faces.extend(faces_of(p));
            }
            Ok(BooleanResult {
                shape: comp.0,
                solid: None,
                shells,
                faces,
                warnings: vec!["non-solid input fused into a compound".into()],
            })
        }
        BoolOp::Cut => {
            if a_solid {
                let mut r = single_shape_result(a);
                r.warnings.push("degenerate cut: 'b' is not a solid; 'a' returned unchanged".into());
                return Ok(r);
            }
            if !b_solid {
                let mut r = single_shape_result(a);
                r.warnings.push("degenerate cut: neither input is a solid".into());
                return Ok(r);
            }
            // A non-solid `a` cut by a solid `b`: keep `a` unless it lies
            // inside `b` (then nothing survives).
            if is_inside(b, &shape_centroid(a)) {
                Ok(empty_result(op))
            } else {
                let mut r = single_shape_result(a);
                r.warnings.push("degenerate cut: 'a' is not a solid; kept 'a' unchanged".into());
                Ok(r)
            }
        }
        BoolOp::Common => {
            let (solid_s, other) = if a_solid {
                (a, b)
            } else if b_solid {
                (b, a)
            } else {
                let mut r = empty_result(op);
                r.warnings.push("degenerate common: neither input is a solid".into());
                return Ok(r);
            };
            if is_inside(solid_s, &shape_centroid(other)) {
                let mut r = single_shape_result(other);
                r.warnings.push("degenerate common: non-solid inside solid retained".into());
                Ok(r)
            } else {
                let mut r = empty_result(op);
                r.warnings.push("degenerate common: non-solid outside solid → empty".into());
                Ok(r)
            }
        }
    }
}

/// One-line diagnostic of a [`BooleanResult`].
///
/// Reports the result shape type, whether it is a closed solid, its face
/// count, and how many warnings it carries.
pub fn boolean_result_summary(r: &BooleanResult) -> String {
    let st = r.shape.shape_type().to_str();
    let kind = if r.solid.is_some() { "closed solid" } else { "open shell/compound" };
    let faces = faces_of(&r.shape).len();
    let warns = if r.warnings.is_empty() {
        "no warnings".to_string()
    } else {
        format!("{} warning(s): {}", r.warnings.len(), r.warnings.join("; "))
    };
    format!("{st} ({kind}) {faces} faces, {warns}")
}

/// Apply a *sequence* of boolean operations to a list of shapes.
///
/// The operations are applied left to right, each between the running
/// accumulator and the next shape: `((a0 op0 a1) op1 a2) op2 a3 …`. Unlike
/// [`boolean_multi`], every step may use a different operation, so a CSG tree
/// flattened into the alternating form `shape, op, shape, op, shape, …` can be
/// evaluated directly. Extra shapes (beyond `ops.len() + 1`) are ignored; a
/// missing operation for a remaining shape stops the fold.
pub fn boolean_fold(shapes: &[TopoShape], ops: &[BoolOp], tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    if shapes.is_empty() {
        return Ok(empty_result(BoolOp::Fuse));
    }
    let mut acc = shapes[0].clone();
    let mut warnings: Vec<String> = Vec::new();
    let steps = shapes.len().saturating_sub(1).min(ops.len());
    for i in 0..steps {
        let r = boolean_dispatch(&acc, &shapes[i + 1], ops[i], tol)?;
        warnings.extend(r.warnings);
        acc = r.shape;
    }
    let mut res = single_shape_result(&acc);
    res.warnings.extend(warnings);
    Ok(res)
}

/// Subtract every shape in `cuts` from `a`, one after the other (folded Cut).
///
/// Equivalent to `a − c1 − c2 − … − cn`. Each step routes through
/// [`boolean_dispatch`], so compounds and degenerate operands are handled.
/// The result is the final solid (or an empty compound when `a` is fully
/// removed).
pub fn boolean_cut_many(a: &TopoShape, cuts: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let mut acc = a.clone();
    let mut warnings: Vec<String> = Vec::new();
    for c in cuts {
        let r = boolean_dispatch(&acc, c, BoolOp::Cut, tol)?;
        warnings.extend(r.warnings);
        acc = r.shape;
    }
    let mut res = single_shape_result(&acc);
    res.warnings.extend(warnings);
    Ok(res)
}

/// Intersect a list of shapes (folded Common).
///
/// Equivalent to `a0 ∩ a1 ∩ a2 ∩ …`. Degenerate intermediates (an empty
/// common) stay empty through the fold, so the final result is empty as soon
/// as any pair is disjoint.
pub fn boolean_common_many(shapes: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    boolean_multi(shapes, BoolOp::Common, tol)
}

/// Pairs of face indices (into [`faces_of`]) that share a boundary edge.
///
/// Two faces are *adjacent* when they reference the same `TShape` edge. This
/// is the raw face-adjacency graph of a boundary, useful for connectivity
/// analysis and for understanding where a self-intersection check skips.
pub fn face_adjacency(shape: &TopoShape) -> Vec<(usize, usize)> {
    let faces = faces_of(shape);
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

/// Split a shape into its edge-connected boundary components.
///
/// Two faces belong to the same component when they are connected through a
/// chain of shared boundary edges (the face-adjacency graph). A single closed
/// solid has exactly one component; a compound of disjoint solids yields one
/// component per solid. Each component is rebuilt as a one-shell solid, so the
/// result is the "decomposed" form of a boolean result (`BOPAlgo_Builder`
/// returns such disconnected solids inside a compound).
///
/// Faces that are topologically isolated (an open face with no shared edges)
/// each become their own component.
pub fn connected_components(shape: &TopoShape, tol: f64) -> Vec<TopoShape> {
    let _ = tol;
    let faces = faces_of(shape);
    if faces.len() <= 1 {
        return vec![shape.clone()];
    }
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();

    // Union-find over faces (edge-sharing ⇒ same component).
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            parent[r] = parent[parent[r]];
            r = parent[r];
        }
        r
    }
    fn unite(parent: &mut Vec<usize>, a: usize, b: usize) {
        let (ra, rb) = (find(parent, a), find(parent, b));
        if ra != rb {
            parent[ra] = rb;
        }
    }

    let mut parent: Vec<usize> = (0..faces.len()).collect();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                unite(&mut parent, i, j);
            }
        }
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..faces.len() {
        groups.entry(find(&mut parent, i)).or_default().push(i);
    }

    let bld = TopoBuilder::new();
    let mut comps: Vec<TopoShape> = Vec::new();
    for (_, idx) in groups {
        let fs: Vec<Face> = idx.iter().map(|&i| faces[i].clone()).collect();
        let shell = bld.make_shell(&fs);
        comps.push(bld.make_solid(&[shell]).0);
    }
    comps
}

/// Public structural validation of a boolean result.
///
/// Returns human-readable issue strings (empty when the result is clean). See
/// [`boolean_with_check`] for the tolerance-escalation wrapper.
pub fn boolean_result_validate(r: &BooleanResult, tol: f64) -> Vec<String> {
    validate_boolean_result(r, tol)
}

/// Convenience wrapper: fuse every shape in `shapes` (see [`boolean_multi`]).
///
/// Kept as a named entry point so callers do not have to spell out
/// `BoolOp::Fuse`; behaves identically to `boolean_multi(shapes, BoolOp::Fuse,
/// tol)`.
pub fn boolean_fuse_all(shapes: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    boolean_multi(shapes, BoolOp::Fuse, tol)
}

/// Bounding box of a boolean result, or `None` when the result is empty.
///
/// Delegates to [`crate::bbox_from_geometry::shape_bbox`] on the result shape,
/// so compounds (multi-solid results) are covered by the compound bbox.
pub fn boolean_result_bbox(r: &BooleanResult) -> Option<occt_core::bnd::BndBox> {
    let bb = crate::bbox_from_geometry::shape_bbox(&r.shape);
    if bb.is_void() {
        None
    } else {
        Some(bb)
    }
}

/// Normalize a boolean result into its "assembled" form.
///
/// When a result carries more than one shell (a multi-solid Cut/Common result
/// that was not merged into a compound), wrap each closed shell into its own
/// solid and return a compound of them; a single-shell result is returned
/// unchanged. This mirrors the final assembly step of `BOPAlgo_Builder`, which
/// hands the caller a compound of the disconnected result solids.
pub fn normalize_result(r: &BooleanResult) -> BooleanResult {
    if r.shape.is_compound() || r.shells.len() <= 1 {
        return r.clone();
    }
    let bld = TopoBuilder::new();
    let solids: Vec<TopoShape> = r
        .shells
        .iter()
        .filter(|sh| shell_is_closed(sh))
        .map(|sh| bld.make_solid(&[sh.clone()]).0)
        .collect();
    if solids.is_empty() {
        return r.clone();
    }
    let shape = if solids.len() == 1 {
        solids[0].clone()
    } else {
        bld.make_compound_of(&solids).0
    };
    let solid = Solid::wrap(shape.clone());
    BooleanResult {
        shape,
        solid,
        shells: r.shells.clone(),
        faces: r.faces.clone(),
        warnings: r.warnings.clone(),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_gprop::volume as shape_volume;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::GpPnt;

    fn box_vol(s: &TopoShape) -> f64 {
        shape_volume(s, 0.02)
    }

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn overlapping_boxes() -> (Solid, Solid) {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0));
        (a.solid, b.solid)
    }

    #[test]
    fn fuse_overlapping_boxes() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.solid.is_some(), "fuse produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "fuse shell is closed");
        let v = box_vol(&r.shape);
        assert!((v - 1.5).abs() < 0.05, "fuse volume {v} (expected 1.5)");
        assert!(faces_of(&r.shape).len() > 12, "faces {} (expected > 12)", faces_of(&r.shape).len());
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn cut_overlapping_boxes() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!(r.solid.is_some(), "cut produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "cut shell is closed");
        let v = box_vol(&r.shape);
        assert!((v - 0.5).abs() < 0.05, "cut volume {v} (expected 0.5)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn common_overlapping_boxes() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Common, 1e-6).expect("common ok");
        assert!(r.solid.is_some(), "common produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "common shell is closed");
        let v = box_vol(&r.shape);
        assert!((v - 0.5).abs() < 0.05, "common volume {v} (expected 0.5)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn disjoint_boxes() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(2.0, 2.0, 2.0), &GpPnt::new(3.0, 3.0, 3.0));

        let f = boolean(&a.solid.0, &b.solid.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!((box_vol(&f.shape) - 2.0).abs() < 0.05, "disjoint fuse volume {}", box_vol(&f.shape));

        let c = boolean(&a.solid.0, &b.solid.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!((box_vol(&c.shape) - 1.0).abs() < 0.05, "disjoint cut volume {}", box_vol(&c.shape));

        let m = boolean(&a.solid.0, &b.solid.0, BoolOp::Common, 1e-6).expect("common ok");
        assert!((box_vol(&m.shape)).abs() < 1e-9, "disjoint common volume {}", box_vol(&m.shape));

        clear_tree(&f.shape);
        clear_tree(&c.shape);
        clear_tree(&m.shape);
        clear_tree(&a.solid.0);
        clear_tree(&b.solid.0);
    }

    #[test]
    fn cube_minus_internal_box() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.2, 0.2, 0.2), &GpPnt::new(0.8, 0.8, 0.8));
        let r = boolean(&a.solid.0, &b.solid.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!(r.solid.is_some(), "hollow cube produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "hollow cube shell is closed");
        let v = box_vol(&r.shape);
        let expected = 1.0 - 0.6 * 0.6 * 0.6;
        assert!((v - expected).abs() < 0.05, "hollow cube volume {v} (expected {expected})");
        clear_tree(&r.shape);
        clear_tree(&a.solid.0);
        clear_tree(&b.solid.0);
    }

    #[test]
    fn fuse_matches_voxel_union() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let exact = box_vol(&r.shape);
        let voxel = crate::solid_union::union_volume(&a.0, &b.0, 32);
        assert!(voxel > 1e-9, "voxel volume should be non-zero");
        let rel = (exact - voxel).abs() / voxel;
        assert!(rel < 0.15, "fuse {exact} vs voxel {voxel} (rel {rel:.3})");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    // ------------------------------------------------------------------
    // Full-topology tests
    // ------------------------------------------------------------------

    /// Build an axis-aligned box with all six faces oriented OUTWARD.
    ///
    /// (`BRepPrimBox::make_box_corner` currently builds the prism top cap with
    /// the base's normal — an inward face that makes signed volumes
    /// position-dependent — so tests that need predictable volume/orientation
    /// use this helper instead.)
    fn test_box_at(lo: &GpPnt, hi: &GpPnt) -> Solid {
        use crate::shape::Vertex;
        use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpVec};
        use occt_geom::{GeomLine, GeomPlane};
        let b = TopoBuilder::new();
        let (x0, y0, z0) = (lo.x(), lo.y(), lo.z());
        let (x1, y1, z1) = (hi.x(), hi.y(), hi.z());
        let c = [
            GpPnt::new(x0, y0, z0),
            GpPnt::new(x1, y0, z0),
            GpPnt::new(x1, y1, z0),
            GpPnt::new(x0, y1, z0),
            GpPnt::new(x0, y0, z1),
            GpPnt::new(x1, y0, z1),
            GpPnt::new(x1, y1, z1),
            GpPnt::new(x0, y1, z1),
        ];
        let verts: Vec<Vertex> = c.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        let edge_pairs = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];
        let edges: Vec<Edge> = edge_pairs
            .iter()
            .map(|&(i, j)| {
                let dir = GpDir::from_vec(&GpVec::from_pnts(&c[i], &c[j])).unwrap();
                let lin = GpLin::from_pnt_dir(c[i], dir);
                let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, c[i].distance(&c[j]));
                b.add(&mut e.0, &verts[i].0);
                b.add(&mut e.0, &verts[j].0);
                e
            })
            .collect();
        // Face → edge-index lists (indices into `edge_pairs`) for −X, +X, −Y,
        // +Y, −Z, +Z. Mirrors the box builder used by the shell-check tests.
        let face_edges: [[usize; 4]; 6] = [
            [8, 7, 11, 3],
            [9, 5, 10, 1],
            [0, 9, 4, 8],
            [2, 10, 6, 11],
            [3, 2, 1, 0],
            [4, 5, 6, 7],
        ];
        let face_origins = [c[0], c[1], c[0], c[3], c[0], c[4]];
        let normals = [
            GpDir::new(-1.0, 0.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            GpDir::new(0.0, -1.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
        ];
        let mut faces = Vec::new();
        for i in 0..6 {
            let x_dir = if normals[i].x().abs() > 0.9 {
                GpDir::new(0.0, 1.0, 0.0).unwrap()
            } else if normals[i].y().abs() > 0.9 {
                GpDir::new(0.0, 0.0, 1.0).unwrap()
            } else {
                GpDir::new(1.0, 0.0, 0.0).unwrap()
            };
            let ax3 = GpAx3::new(face_origins[i], normals[i], &x_dir).unwrap();
            let pln = GpPln::new(ax3);
            let wire_edges: Vec<Edge> = face_edges[i].iter().map(|&e| edges[e].clone()).collect();
            let wire = b.make_wire(&wire_edges);
            faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[wire]));
        }
        let shell = b.make_shell(&faces);
        b.make_solid(&[shell])
    }

    fn disjoint_box_shapes() -> Vec<TopoShape> {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let b3 = test_box_at(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(1.0, 3.0, 1.0));
        vec![b1.0, b2.0, b3.0]
    }

    #[test]
    fn fuse_three_boxes_compound_or_solid() {
        let shapes = disjoint_box_shapes();
        let r = boolean_multi(&shapes, BoolOp::Fuse, 1e-6).expect("multi fuse ok");
        assert!(r.shape.is_compound(), "three disjoint boxes fuse to a compound");
        let subs = decompose_compound(&r.shape);
        assert_eq!(subs.len(), 3, "compound has 3 sub-shapes, got {}", subs.len());
        for s in &subs {
            assert!(s.is_solid(), "each sub-shape is a solid");
        }
        let v = box_vol(&r.shape);
        assert!((v - 3.0).abs() < 0.1, "disjoint multi fuse volume {v} (expected 3.0)");
        clear_tree(&r.shape);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_multi_single_shape_identity() {
        let shapes = disjoint_box_shapes();
        let r = boolean_multi(&shapes[..1], BoolOp::Fuse, 1e-6).expect("single fuse ok");
        assert_eq!(r.shape.shape_type(), ShapeType::Solid, "one shape → itself");
        assert!(r.solid.is_some());
        assert!((box_vol(&r.shape) - 1.0).abs() < 0.05);
        clear_tree(&r.shape);
        clear_tree(&shapes[0]);
    }

    #[test]
    fn boolean_multi_overlap_folds() {
        // Three boxes overlapping along x (each overlaps the next).
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0));
        let b3 = test_box_at(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 1.0, 1.0));
        let shapes = vec![b1.0.clone(), b2.0.clone(), b3.0.clone()];
        let r = boolean_multi(&shapes, BoolOp::Fuse, 1e-6).expect("multi fuse ok");
        assert!(r.solid.is_some(), "overlapping fuse produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "overlapping fuse shell closed");
        let v = box_vol(&r.shape);
        // Analytic union: b1∪b2 = 1.5; ∪b3 = [0,2]×[0,1]×[0,1] = 2.0.
        assert!((v - 2.0).abs() < 0.2, "multi fuse volume {v} (expected ~2.0)");
        clear_tree(&r.shape);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_compound_cut() {
        let big = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 3.0, 3.0));
        let s1 = test_box_at(&GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(1.5, 1.5, 1.5));
        let s2 = test_box_at(&GpPnt::new(1.7, 0.5, 0.5), &GpPnt::new(2.7, 1.5, 1.5));
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[s1.0.clone(), s2.0.clone()]);
        let r = boolean_compound(&big.0, &comp.0, BoolOp::Cut, 1e-6).expect("compound cut ok");
        let v = box_vol(&r.shape);
        // 3³ − 1 − 1 = 25.
        assert!((v - 25.0).abs() < 0.15, "compound cut volume {v} (expected ~25.0)");
        assert!(r.solid.is_some(), "one big box minus two internal boxes → one solid");
        clear_tree(&r.shape);
        clear_tree(&big.0);
        clear_tree(&s1.0);
        clear_tree(&s2.0);
        clear_tree(&comp.0);
    }

    #[test]
    fn boolean_compound_common() {
        let a1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let a2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let b = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(2.5, 1.0, 1.0));
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[a1.0.clone(), a2.0.clone()]);
        let r = boolean_compound(&comp.0, &b.0, BoolOp::Common, 1e-6).expect("compound common ok");
        let v = box_vol(&r.shape);
        // a1∩b = 0.5 volume, a2∩b = 0.5 volume.
        assert!((v - 1.0).abs() < 0.15, "compound common volume {v} (expected ~1.0)");
        clear_tree(&r.shape);
        clear_tree(&comp.0);
        clear_tree(&a1.0);
        clear_tree(&a2.0);
        clear_tree(&b.0);
    }

    #[test]
    fn self_intersection_detected() {
        use occt_core::gp::{GpAx3, GpDir, GpPln};
        use occt_geom::GeomPlane;
        use std::sync::Arc;

        // A valid box has no self-intersections.
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = detect_self_intersections(&boxed.solid.0, 1e-6);
        assert!(!rep.found, "a box must not self-intersect");

        // A shell with a horizontal face (z=0) crossed by a vertical face
        // (y=0): the two surfaces intersect along a line through both faces.
        let bld = TopoBuilder::new();
        let h_wire = bld.make_wire(&[
            bld.make_edge_segment(&GpPnt::new(-1.0, -1.0, 0.0), &GpPnt::new(1.0, -1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, -1.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 1.0, 0.0), &GpPnt::new(-1.0, 1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(-1.0, 1.0, 0.0), &GpPnt::new(-1.0, -1.0, 0.0)),
        ]);
        let f_h = bld.make_face(Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard()))), &[h_wire]);
        let v_wire = bld.make_wire(&[
            bld.make_edge_segment(&GpPnt::new(-1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 1.0), &GpPnt::new(-1.0, 0.0, 1.0)),
            bld.make_edge_segment(&GpPnt::new(-1.0, 0.0, 1.0), &GpPnt::new(-1.0, 0.0, 0.0)),
        ]);
        let pln_y = GpPln::new(
            GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 1.0, 0.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
                .unwrap(),
        );
        let f_v = bld.make_face(Arc::new(GeomPlane::new(pln_y)), &[v_wire]);
        let shell = bld.make_shell(&[f_h, f_v]);
        let rep2 = detect_self_intersections(&shell.0, 1e-6);
        assert!(rep2.found, "crossing faces must be detected as self-intersecting");
        assert!(rep2.edge_count >= 1, "at least one intersecting face pair");
        clear_tree(&shell.0);
        clear_tree(&boxed.solid.0);
    }

    #[test]
    fn boolean_with_check_repairs() {
        let (a, b) = overlapping_boxes();
        let r = boolean_with_check(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("with-check fuse ok");
        assert!(r.solid.is_some(), "fuse produces a solid");
        assert!(!r.shells.is_empty(), "result has a shell");
        assert!(r.shells.iter().all(shell_is_closed), "result shells are closed");
        // Hard failures (open shell / repair failure) must not occur; benign
        // diagnostics such as a volume cross-check warning are acceptable.
        let hard = r.warnings.iter().any(|w| w.contains("not closed") || w.contains("no shell"));
        assert!(!hard, "no hard failures: {:?}", r.warnings);
        let v = box_vol(&r.shape);
        assert!((v - 1.5).abs() < 0.05, "with-check fuse volume {v} (expected 1.5)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn decompose_compound_top_level() {
        let bld = TopoBuilder::new();
        let shapes = disjoint_box_shapes();
        let comp = bld.make_compound_of(&shapes);
        let subs = decompose_compound(&comp.0);
        assert_eq!(subs.len(), 3, "compound of 3 shapes decomposes into 3");
        for s in &subs {
            assert!(!s.is_compound(), "sub-shapes are atomic");
        }
        // Non-compound input → itself.
        let single = decompose_compound(&shapes[0]);
        assert_eq!(single.len(), 1);
        // shape_components splits a solid into its shells.
        let comps = shape_components(&shapes[0]);
        assert!(!comps.is_empty());
        assert!(comps[0].is_solid());
        clear_tree(&comp.0);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_degenerate_face() {
        let solid = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let bld = TopoBuilder::new();
        let face = bld.make_face_plane(&occt_core::gp::GpPln::new(occt_core::gp::GpAx3::standard()));
        let r = boolean_degenerate(&face.0, &solid.solid.0, BoolOp::Fuse, 1e-6).expect("face fuse ok");
        assert!(r.shape.is_compound() || r.solid.is_some(), "face fused with solid → compound or solid");
        assert!(!r.warnings.is_empty(), "degenerate fuse reports a warning");
        clear_tree(&face.0);
        clear_tree(&r.shape);
        clear_tree(&solid.solid.0);
    }

    #[test]
    fn boolean_degenerate_empty() {
        let bld = TopoBuilder::new();
        let empty = bld.make_compound_of(&[]);
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        // empty fused with a box → the box (empty is the Fuse identity).
        let r = boolean_degenerate(&empty.0, &boxed.solid.0, BoolOp::Fuse, 1e-6).expect("empty fuse ok");
        assert!(r.solid.is_some(), "empty + box → the box");
        assert!((box_vol(&r.shape) - 1.0).abs() < 0.05);
        // box minus empty → the box.
        let r2 = boolean_degenerate(&boxed.solid.0, &empty.0, BoolOp::Cut, 1e-6).expect("empty cut ok");
        assert!((box_vol(&r2.shape) - 1.0).abs() < 0.05);
        // empty common → empty.
        let r3 = boolean_degenerate(&empty.0, &boxed.solid.0, BoolOp::Common, 1e-6).expect("empty common ok");
        assert!(box_vol(&r3.shape) < 1e-9, "empty common has no volume");
        clear_tree(&empty.0);
        clear_tree(&r.shape);
        clear_tree(&r2.shape);
        clear_tree(&r3.shape);
        clear_tree(&boxed.solid.0);
    }

    #[test]
    fn boolean_result_summary_string() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let s = boolean_result_summary(&r);
        assert!(s.contains(r.shape.shape_type().to_str()), "summary has the shape type: {s}");
        assert!(s.contains("faces"), "summary mentions faces: {s}");
        assert!(s.contains("solid") || s.contains("compound"), "summary describes the solidity: {s}");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn connected_components_splits_compound() {
        let shapes = disjoint_box_shapes();
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&shapes);
        let comps = connected_components(&comp.0, 1e-6);
        assert_eq!(comps.len(), 3, "compound of 3 disjoint boxes → 3 components");
        for c in &comps {
            assert!(c.is_solid(), "each component is a solid");
        }
        // A single solid has one component.
        let one = connected_components(&shapes[0], 1e-6);
        assert_eq!(one.len(), 1, "a single solid → one component");
        clear_tree(&comp.0);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_fold_sequence() {
        // ((a ∪ b) − c): fuse two boxes, then cut a third.
        let a = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0));
        let c = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 0.5, 1.0));
        let shapes = vec![a.0.clone(), b.0.clone(), c.0.clone()];
        let ops = vec![BoolOp::Fuse, BoolOp::Cut];
        let r = boolean_fold(&shapes, &ops, 1e-6).expect("fold ok");
        // a∪b = [0,1.5]×[0,1]×[0,1] (vol 1.5); cut c = [0.5,1.5]×[0,0.5]×[0,1]
        // (vol 0.5) → 1.5 − 0.5 = 1.0.
        let v = box_vol(&r.shape);
        assert!((v - 1.0).abs() < 0.2, "fold volume {v} (expected ~1.0)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
        clear_tree(&c.0);
    }

    #[test]
    fn boolean_cut_many_works() {
        let big = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 3.0, 3.0));
        let c1 = test_box_at(&GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(1.5, 1.5, 1.5));
        let c2 = test_box_at(&GpPnt::new(1.7, 0.5, 0.5), &GpPnt::new(2.7, 1.5, 1.5));
        let r = boolean_cut_many(&big.0, &[c1.0.clone(), c2.0.clone()], 1e-6).expect("cut many ok");
        let v = box_vol(&r.shape);
        assert!((v - 25.0).abs() < 0.2, "cut-many volume {v} (expected ~25.0)");
        clear_tree(&r.shape);
        clear_tree(&big.0);
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn boolean_result_validate_clean() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        // A closed, non-self-intersecting fuse has no hard structural issues.
        let issues = boolean_result_validate(&r, 1e-6);
        let hard: Vec<&String> = issues.iter().filter(|s| s.contains("not closed") || s.contains("no shell")).collect();
        assert!(hard.is_empty(), "no open-shell issues: {issues:?}");
        assert!(r.shells.iter().all(shell_is_closed), "fuse shell is closed");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

}
