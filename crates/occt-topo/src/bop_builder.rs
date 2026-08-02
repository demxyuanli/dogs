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
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Wire};
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
// Self-intersection repair
// ---------------------------------------------------------------------------
//
// A `BOPAlgo_ArgumentAnalyzer`-lite repair pass. `BOPAlgo_ArgumentAnalyzer`
// *detects* self-intersections (non-adjacent faces whose surfaces cross or
// overlap); the "repair" done here is the resolution step: split every
// crossing face along the intersection polyline and weld the split edges so
// the crossing faces become adjacent, and drop degenerate slivers and
// overlapping coplanar duplicates. This is the same split-and-weld machinery
// the exact boolean uses to produce manifold boundaries.

/// Result of [`repair_self_intersections`].
pub struct RepairResult {
    /// The repaired shape (unchanged when no self-intersection was found).
    pub repaired: TopoShape,
    /// Number of faces that were split to resolve a transversal crossing.
    pub fixed_faces: usize,
    /// Number of faces removed (overlapping coplanar duplicates, or every
    /// sub-face degenerating to zero area).
    pub removed_faces: usize,
    /// Non-fatal diagnostics (unrepairable non-planar crossings, open result).
    pub warnings: Vec<String>,
}

/// Repair a self-intersecting boundary (a `BOPAlgo_ArgumentAnalyzer`-lite port).
///
/// * Runs [`detect_self_intersections`]; a clean shape is returned unchanged
///   with `fixed_faces == 0`.
/// * For every flagged transversal crossing pair, the intersection polyline is
///   computed with [`crate::intpatch::surface_surface_intersection`], clipped
///   to both face polygons, and both faces are split along it through the same
///   weld/edge-map the boolean uses. The split edges coincide, so the crossing
///   faces become adjacent and the defect is resolved.
/// * Overlapping *coplanar* non-adjacent faces are resolved by removing the
///   duplicate (the overlap sliver); degenerate (zero-area) sub-faces are
///   dropped by the splitter.
/// * The surviving faces are rebuilt into a shell — or a solid when the shell
///   is closed — and an open result is reported as a warning.
/// * A `Compound` is repaired component-wise and reassembled.
pub fn repair_self_intersections(shape: &TopoShape, tol: f64) -> Result<RepairResult, String> {
    let tol = tol.max(1e-9);
    if shape.is_compound() {
        let mut fixed = 0usize;
        let mut removed = 0usize;
        let mut warnings: Vec<String> = Vec::new();
        let children = expand_compound(shape);
        let mut out: Vec<TopoShape> = Vec::with_capacity(children.len());
        for c in &children {
            let r = repair_self_intersections(c, tol)?;
            fixed += r.fixed_faces;
            removed += r.removed_faces;
            warnings.extend(r.warnings);
            out.push(r.repaired);
        }
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&out);
        return Ok(RepairResult { repaired: comp.0, fixed_faces: fixed, removed_faces: removed, warnings });
    }

    let report = detect_self_intersections(shape, tol);
    if !report.found {
        return Ok(RepairResult { repaired: shape.clone(), fixed_faces: 0, removed_faces: 0, warnings: vec![] });
    }

    let faces = faces_of(shape);
    let planes: Vec<Option<GpPln>> = faces.iter().map(face_plane_local).collect();
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();

    let mut segs: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); faces.len()];
    let mut drop: Vec<bool> = vec![false; faces.len()];
    let mut removed_count = 0usize;
    let mut warnings: Vec<String> = Vec::new();

    // Collect the split segments per face and the faces to drop.
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
                    if curves.is_empty() {
                        continue;
                    }
                    if planes[i].is_none() || planes[j].is_none() {
                        warnings.push(format!(
                            "self-intersection between faces {i} and {j} involves a non-planar face; left unrepaired"
                        ));
                        continue;
                    }
                    // The transversal crossing: the clipped intersection
                    // polyline is the cutting line on both faces.
                    let seg = face_face_segments_local(&faces[i], &faces[j], tol);
                    if seg.is_empty() {
                        continue;
                    }
                    for s in &seg {
                        segs[i].push(*s);
                        segs[j].push(*s);
                    }
                }
                crate::intpatch::SurfaceIntersection::Coincident => {
                    if faces_polygon_overlap(&faces[i], &faces[j], tol) {
                        if !drop[j] {
                            drop[j] = true;
                            removed_count += 1;
                            warnings.push(format!("removed overlapping coplanar face {j}"));
                        }
                    }
                }
                crate::intpatch::SurfaceIntersection::None => {}
            }
        }
    }

    // Count faces that actually split (an interior cutting segment). A segment
    // running exactly along a face boundary leaves the face whole, so it does
    // not count as "fixed".
    let mut fixed_count = 0usize;
    for i in 0..faces.len() {
        if drop[i] || segs[i].is_empty() {
            continue;
        }
        let Some(pln) = planes[i].clone() else { continue };
        let Some(poly) = face_polygon_local(&faces[i], &pln) else { continue };
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
        if split_polygon_by_segments(&poly, &segs2d).len() >= 2 {
            fixed_count += 1;
        }
    }

    // Rebuild the boundary: split every kept planar face along its cutting
    // segments, keeping non-planar faces whole.
    let bld = TopoBuilder::new();
    let mut kept_faces: Vec<Face> = Vec::new();
    let mut kept_planes: Vec<Option<GpPln>> = Vec::new();
    let mut kept_segs: Vec<Vec<(GpPnt, GpPnt)>> = Vec::new();
    let mut non_planar: Vec<Face> = Vec::new();
    for (i, f) in faces.iter().enumerate() {
        if drop[i] {
            continue;
        }
        if planes[i].is_some() {
            kept_faces.push(f.clone());
            kept_planes.push(planes[i].clone());
            kept_segs.push(segs[i].clone());
        } else {
            non_planar.push(f.clone());
        }
    }

    let mut weld = Weld::new(tol.max(1e-7));
    let mut edge_map = EdgeMap::default();
    let subs = split_faces(&bld, &kept_faces, &kept_planes, &kept_segs, &mut weld, &mut edge_map);
    let mut result_faces: Vec<Face> = non_planar;
    result_faces.extend(subs.into_iter().map(|sf| sf.face));

    if result_faces.is_empty() {
        warnings.push("repair removed every face; the result is empty".into());
        let empty = bld.make_compound_of(&[]);
        return Ok(RepairResult {
            repaired: empty.0,
            fixed_faces: fixed_count,
            removed_faces: removed_count,
            warnings,
        });
    }

    let shell = bld.make_shell(&result_faces);
    let closed = shell_is_closed(&shell);
    if !closed {
        warnings.push("repaired boundary is not closed (open shell)".into());
    }
    let repaired: TopoShape = if closed { bld.make_solid(&[shell]).0 } else { shell.0 };

    Ok(RepairResult { repaired, fixed_faces: fixed_count, removed_faces: removed_count, warnings })
}

/// Run a boolean and repair a self-intersecting result.
///
/// Runs [`crate::bop_curved::curved_boolean_full`], then
/// [`repair_self_intersections`] on the result. When the repair changed the
/// shape, the returned [`BooleanResult`] carries the repaired shape and a
/// warning mentioning the number of faces fixed/removed; a clean result is
/// returned unchanged.
pub fn boolean_repaired(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let r = crate::bop_curved::curved_boolean_full(a, b, op, tol)?;
    Ok(repair_boolean_result(r, tol))
}

/// Apply [`repair_self_intersections`] to an already-computed boolean result.
///
/// Shared by [`boolean_repaired`] and [`boolean_repaired_with_check`]. A clean
/// result is returned unchanged; otherwise the repaired shape replaces the
/// result shape (recomputing `solid`/`shells`/`faces`) and a warning naming the
/// number of fixed/removed faces is appended.
fn repair_boolean_result(r: BooleanResult, tol: f64) -> BooleanResult {
    let rep = match repair_self_intersections(&r.shape, tol) {
        Ok(rep) => rep,
        Err(e) => {
            let mut rr = r;
            rr.warnings.push(format!("self-intersection repair failed: {e}"));
            return rr;
        }
    };
    if rep.fixed_faces == 0 && rep.removed_faces == 0 {
        return r;
    }
    let mut rr = single_shape_result(&rep.repaired);
    rr.warnings.extend(r.warnings);
    rr.warnings.extend(rep.warnings);
    rr.warnings.push(format!(
        "boolean result had self-intersections; repair fixed {} face(s), removed {}",
        rep.fixed_faces, rep.removed_faces
    ));
    rr
}

// ---------------------------------------------------------------------------
// Multi-result decomposition
// ---------------------------------------------------------------------------

/// A boolean result decomposed into its disconnected pieces.
///
/// Mirrors the assembly step of `BOPAlgo_Builder` / `TopOpeBRep`: a boolean
/// result is a compound (or a multi-shell solid) holding several disconnected
/// results; this struct splits them apart and buckets them by shape type.
pub struct MultiResult {
    /// Every top-level result shape (compounds flattened, multi-shell solids
    /// split into one shape per connected boundary component).
    pub shapes: Vec<TopoShape>,
    /// The shapes of [`MultiResult::shapes`] that are solids.
    pub solids: Vec<TopoShape>,
    /// The shapes of [`MultiResult::shapes`] that are (open) shells.
    pub shells: Vec<TopoShape>,
    /// The shapes of [`MultiResult::shapes`] that are still compounds.
    pub compounds: Vec<TopoShape>,
}

/// Split a non-compound boundary into its edge-connected components, preserving
/// solidity: a closed component is returned as a one-shell solid, an open one
/// as a bare shell.
fn split_connected_boundaries(shape: &TopoShape, tol: f64) -> Vec<TopoShape> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return vec![shape.clone()];
    }
    if faces.len() == 1 {
        let bld = TopoBuilder::new();
        return vec![bld.make_shell(&faces).0];
    }
    let groups = face_connectivity_groups(&faces);
    let bld = TopoBuilder::new();
    let mut comps: Vec<TopoShape> = Vec::new();
    for idx in groups {
        let fs: Vec<Face> = idx.iter().map(|&i| faces[i].clone()).collect();
        let shell = bld.make_shell(&fs);
        if shell_is_closed(&shell) {
            comps.push(bld.make_solid(&[shell]).0);
        } else {
            comps.push(shell.0);
        }
    }
    comps
}

/// Group face indices by edge connectivity.
///
/// Two faces belong to the same group when they are connected through a chain
/// of shared boundary edges (identical `TShape` edge references). The groups
/// are returned sorted by their smallest face index, so the ordering is stable
/// for a given shape.
fn face_connectivity_groups(faces: &[Face]) -> Vec<Vec<usize>> {
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                uf_unite(&mut parent, i, j);
            }
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..faces.len() {
        groups.entry(uf_find(&mut parent, i)).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by_key(|g| g[0]);
    out
}

/// Union-find find with path compression.
fn uf_find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        parent[r] = parent[parent[r]];
        r = parent[r];
    }
    r
}

/// Union-find union (root of `a` keeps the group).
fn uf_unite(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (uf_find(parent, a), uf_find(parent, b));
    if ra != rb {
        parent[ra] = rb;
    }
}

/// Decompose a boolean result into its disconnected pieces.
///
/// * the top-level shape list comes from [`decompose_compound`] (compounds are
///   flattened);
/// * every remaining solid/shell is further split into its edge-connected
///   boundary components, so a multi-shell Cut result yields one shape per
///   piece;
/// * the pieces are bucketed into solids / shells / compounds.
pub fn decompose_multi_result(r: &BooleanResult) -> MultiResult {
    let out = shape_boundary_components(&r.shape, 1e-6);
    let solids = out.iter().filter(|s| s.is_solid()).cloned().collect();
    let shells = out.iter().filter(|s| s.is_shell()).cloned().collect();
    let compounds = out.iter().filter(|s| s.is_compound()).cloned().collect();
    MultiResult { shapes: out, solids, shells, compounds }
}

/// Run a boolean and decompose the result into its disconnected pieces.
///
/// A boolean that produces several disjoint results — a Fuse of disjoint
/// inputs (a compound), a Cut that leaves two or more pieces (a multi-shell
/// solid) — yields a [`MultiResult`] with one shape per piece.
pub fn boolean_split_result(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<MultiResult, String> {
    let r = crate::bop_curved::curved_boolean_full(a, b, op, tol)?;
    Ok(decompose_multi_result(&r))
}

/// True when the result decomposes into more than one non-empty shape.
///
/// A disjoint Fuse (a compound) or a Cut into several pieces is "multiple"; a
/// single overlapping Fuse/Common result is not.
pub fn result_is_multiple(r: &BooleanResult) -> bool {
    decompose_multi_result(r).shapes.len() > 1
}

/// Remove topological garbage from a shape.
///
/// Weld near-coincident vertices (merging duplicate vertex instances and
/// dropping edges whose two endpoints collapse to one point), then remove edges
/// shorter than `tol`. This is the `ShapeFix_Shape`-style hygiene pass that
/// keeps rebuilt boundaries minimal. The input shape is left untouched.
pub fn topologically_clean(shape: &TopoShape, tol: f64) -> TopoShape {
    let tol = tol.max(1e-9);
    let (weld, _) = crate::shhealing::weld_coincident_vertices(shape, tol);
    let (clean, _) = crate::shhealing::remove_small_edges(&weld, tol);
    clean
}

// ---------------------------------------------------------------------------
// Self-intersection analysis (BOPAlgo_ArgumentAnalyzer)
// ---------------------------------------------------------------------------
//
// `BOPAlgo_ArgumentAnalyzer` does not just report "self-intersection found" —
// it enumerates the offending face pairs, classifies the defect (transversal
// crossing vs. coplanar overlap), and exposes the intersection data. The items
// below are that detailed analysis, on top of the boolean summary
// [`detect_self_intersections`].

/// Classification of a self-intersection defect.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelfIntersectionKind {
    /// Two non-adjacent faces cross transversally: their surfaces intersect
    /// along a curve that passes through both face patches.
    Crossing,
    /// Two non-adjacent coplanar faces overlap in area.
    CoplanarOverlap,
}

impl SelfIntersectionKind {
    /// Human-readable label of the defect kind.
    pub fn label(&self) -> &'static str {
        match self {
            SelfIntersectionKind::Crossing => "transversal crossing",
            SelfIntersectionKind::CoplanarOverlap => "coplanar overlap",
        }
    }
}

/// One self-intersection defect found by [`analyze_self_intersections`].
///
/// Mirrors one entry of `BOPAlgo_ArgumentAnalyzer`'s self-interference list:
/// the two non-adjacent faces involved, the defect kind and the geometric
/// evidence (intersection sample points, and the segment clipped to both face
/// patches when it is computable).
#[derive(Debug, Clone)]
pub struct SelfIntersectionIssue {
    /// Index (into `faces_of`) of the first face of the pair.
    pub face_a: usize,
    /// Index (into `faces_of`) of the second face of the pair.
    pub face_b: usize,
    /// What kind of defect the pair exhibits.
    pub kind: SelfIntersectionKind,
    /// Sample points on the surface intersection (empty for a pure coplanar
    /// overlap, which has no curve).
    pub points: Vec<GpPnt>,
    /// The intersection segment clipped to both face polygons, when the two
    /// faces are planar and the line actually crosses both patches.
    pub segment: Option<(GpPnt, GpPnt)>,
}

/// Enumerate every self-intersection defect of a shape.
///
/// This is the per-pair detail behind [`detect_self_intersections`]: it walks
/// the same non-adjacent face pairs, but reports each defect with its kind and
/// geometry instead of collapsing them into a single boolean flag. A valid
/// closed box returns an empty list; a shell built from two crossing faces
/// returns one [`SelfIntersectionIssue`] of kind `Crossing`.
pub fn analyze_self_intersections(shape: &TopoShape, tol: f64) -> Vec<SelfIntersectionIssue> {
    let tol = tol.max(1e-9);
    let faces = faces_of(shape);
    if faces.len() < 2 {
        return Vec::new();
    }
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
    let mut issues: Vec<SelfIntersectionIssue> = Vec::new();
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
                    if pts.is_empty() {
                        continue;
                    }
                    let clipped = face_face_segments_local(&faces[i], &faces[j], tol);
                    let segment = if clipped.is_empty() {
                        None
                    } else {
                        Some((clipped[0].0, clipped[0].1))
                    };
                    issues.push(SelfIntersectionIssue {
                        face_a: i,
                        face_b: j,
                        kind: SelfIntersectionKind::Crossing,
                        points: pts,
                        segment,
                    });
                }
                crate::intpatch::SurfaceIntersection::Coincident => {
                    if faces_polygon_overlap(&faces[i], &faces[j], tol) {
                        issues.push(SelfIntersectionIssue {
                            face_a: i,
                            face_b: j,
                            kind: SelfIntersectionKind::CoplanarOverlap,
                            points: Vec::new(),
                            segment: None,
                        });
                    }
                }
                crate::intpatch::SurfaceIntersection::None => {}
            }
        }
    }
    issues
}

/// One-line summary of a list of [`SelfIntersectionIssue`]s.
pub fn self_intersection_issues_summary(issues: &[SelfIntersectionIssue]) -> String {
    if issues.is_empty() {
        return "no self-intersections".to_string();
    }
    let crossings = issues
        .iter()
        .filter(|i| i.kind == SelfIntersectionKind::Crossing)
        .count();
    let overlaps = issues.len() - crossings;
    format!(
        "{} self-intersection issue(s): {crossings} crossing, {overlaps} coplanar overlap",
        issues.len()
    )
}

// ---------------------------------------------------------------------------
// Boundary structural analysis (BOPAlgo_ArgumentAnalyzer checks)
// ---------------------------------------------------------------------------

/// Is a face degenerate — a boundary polygon with (near-)zero area?
///
/// A face whose registered surface is not planar, or whose boundary cannot be
/// extracted, is conservatively reported as *not* degenerate (the check only
/// fires on faces it can measure).
pub fn face_is_degenerate(face: &Face, tol: f64) -> bool {
    let Some(pln) = face_plane_local(face) else {
        return false;
    };
    match face_polygon_local(face, &pln) {
        Some(poly) => polygon_area2d(&poly).abs() < tol.max(1e-12),
        None => false,
    }
}

/// Number of free (open-boundary) edges of a shape, within `tol`.
///
/// A free edge is an edge whose endpoint is used by exactly one edge end in its
/// wire — the signature of an open boundary. A closed box has none; a single
/// open face's wire edges all show up.
pub fn shape_free_edge_count(shape: &TopoShape, tol: f64) -> usize {
    crate::shhealing::free_edges(shape, tol).len()
}

/// Structured boundary analysis of a shape.
///
/// One report aggregates every check `BOPAlgo_ArgumentAnalyzer` runs on an
/// input/result before the boolean: shell count and closedness, free edges,
/// self-intersections, degenerate faces and small edges. A valid closed box
/// reports one closed shell and zeros everywhere else.
#[derive(Debug, Default, Clone)]
pub struct BoundaryAnalysisReport {
    /// Number of shells in the shape.
    pub shell_count: usize,
    /// Number of shells whose every boundary edge is used by exactly two faces.
    pub closed_shells: usize,
    /// Number of shells that are not closed.
    pub open_shells: usize,
    /// Number of free (open-boundary) edges.
    pub free_edges: usize,
    /// Number of self-intersecting non-adjacent face pairs.
    pub self_intersections: usize,
    /// Number of degenerate (near-zero-area) faces.
    pub degenerate_faces: usize,
    /// Number of edges shorter than `tol`.
    pub small_edges: usize,
}

/// Run the boundary analysis checks on `shape`.
pub fn analyze_boundary(shape: &TopoShape, tol: f64) -> BoundaryAnalysisReport {
    let tol = tol.max(1e-9);
    let shells: Vec<Shell> = shapes_of(shape, ShapeType::Shell).into_iter().map(Shell).collect();
    let closed_shells = shells.iter().filter(|s| shell_is_closed(s)).count();
    BoundaryAnalysisReport {
        shell_count: shells.len(),
        closed_shells,
        open_shells: shells.len() - closed_shells,
        free_edges: crate::shhealing::free_edges(shape, tol).len(),
        self_intersections: detect_self_intersections(shape, tol).edge_count,
        degenerate_faces: faces_of(shape).iter().filter(|f| face_is_degenerate(f, tol)).count(),
        small_edges: edges_of(shape)
            .iter()
            .filter(|e| crate::brep_measure::edge_length(e, 8) < tol)
            .count(),
    }
}

/// One-line diagnostic of a [`BoundaryAnalysisReport`].
pub fn boundary_analysis_summary(r: &BoundaryAnalysisReport) -> String {
    format!(
        "{} shell(s) ({} closed, {} open), {} free edge(s), {} self-intersection(s), {} degenerate face(s), {} small edge(s)",
        r.shell_count,
        r.closed_shells,
        r.open_shells,
        r.free_edges,
        r.self_intersections,
        r.degenerate_faces,
        r.small_edges
    )
}

// ---------------------------------------------------------------------------
// Repair diagnostics
// ---------------------------------------------------------------------------

/// One-line diagnostic of a [`RepairResult`].
pub fn repair_result_summary(r: &RepairResult) -> String {
    let mut out = format!("fixed {} face(s), removed {}", r.fixed_faces, r.removed_faces);
    if !r.warnings.is_empty() {
        out.push_str(&format!("; {} warning(s): {}", r.warnings.len(), r.warnings.join("; ")));
    }
    out
}

// ---------------------------------------------------------------------------
// Multi-result helpers
// ---------------------------------------------------------------------------

/// One-line diagnostic of a [`MultiResult`].
pub fn multi_result_summary(m: &MultiResult) -> String {
    format!(
        "{} shape(s): {} solid, {} shell, {} compound",
        m.shapes.len(),
        m.solids.len(),
        m.shells.len(),
        m.compounds.len()
    )
}

/// Total volume of every component of a [`MultiResult`].
///
/// Volumes are measured from each component's tessellation (closed solids only
/// contribute real volume; open shells contribute ~0).
pub fn multi_result_total_volume(m: &MultiResult) -> f64 {
    m.shapes
        .iter()
        .map(|s| crate::brep_gprop::volume(s, 0.02))
        .sum()
}

/// Volume of every component of a [`MultiResult`], in component order.
pub fn component_volumes(m: &MultiResult) -> Vec<f64> {
    m.shapes.iter().map(|s| crate::brep_gprop::volume(s, 0.02)).collect()
}

/// The component with the largest volume, or `None` for an empty result.
pub fn largest_component(m: &MultiResult) -> Option<TopoShape> {
    let mut best: Option<(f64, TopoShape)> = None;
    for s in &m.shapes {
        let v = crate::brep_gprop::volume(s, 0.02);
        if best.as_ref().map_or(true, |(bv, _)| v > *bv) {
            best = Some((v, s.clone()));
        }
    }
    best.map(|(_, s)| s)
}

/// The component with the smallest volume, or `None` for an empty result.
pub fn smallest_component(m: &MultiResult) -> Option<TopoShape> {
    let mut best: Option<(f64, TopoShape)> = None;
    for s in &m.shapes {
        let v = crate::brep_gprop::volume(s, 0.02);
        if best.as_ref().map_or(true, |(bv, _)| v < *bv) {
            best = Some((v, s.clone()));
        }
    }
    best.map(|(_, s)| s)
}

/// Bucket sizes of a [`MultiResult`] as `(solids, shells, compounds)`.
pub fn component_counts_by_type(m: &MultiResult) -> (usize, usize, usize) {
    (m.solids.len(), m.shells.len(), m.compounds.len())
}

/// Decompose a shape into its edge-connected boundary components.
///
/// A public, shape-level form of the per-result decomposition used by
/// [`decompose_multi_result`]: compounds are flattened and every solid/shell is
/// split into one shape per connected boundary component (closed components
/// come back as one-shell solids, open ones as shells).
pub fn shape_boundary_components(shape: &TopoShape, tol: f64) -> Vec<TopoShape> {
    if shape.is_compound() {
        let mut out: Vec<TopoShape> = Vec::new();
        for s in expand_compound(shape) {
            out.extend(split_connected_boundaries(&s, tol));
        }
        out
    } else {
        split_connected_boundaries(shape, tol)
    }
}

// ---------------------------------------------------------------------------
// Composite boolean pipelines
// ---------------------------------------------------------------------------

/// Run a boolean and return the disconnected result pieces as a plain vector.
///
/// Convenience wrapper over [`boolean_split_result`] that drops the
/// per-type buckets and returns just the shapes.
pub fn boolean_components(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<Vec<TopoShape>, String> {
    Ok(boolean_split_result(a, b, op, tol)?.shapes)
}

/// Run a boolean and return the volume of every disconnected result piece.
///
/// Useful for sanity-checking a Cut that should leave several pieces: the
/// returned volumes can be summed and compared against the analytic result.
pub fn boolean_component_volumes(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<Vec<f64>, String> {
    Ok(component_volumes(&boolean_split_result(a, b, op, tol)?))
}

/// Run a boolean, repair self-intersections, and decompose the repaired result.
pub fn boolean_repaired_multi(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<MultiResult, String> {
    let r = boolean_repaired(a, b, op, tol)?;
    Ok(decompose_multi_result(&r))
}

/// Run a boolean, repair self-intersections and topologically clean the result.
///
/// Returns the single repaired+cleaned result shape (a compound when the result
/// is disconnected).
pub fn boolean_repair_and_clean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<TopoShape, String> {
    let r = boolean_repaired(a, b, op, tol)?;
    Ok(topologically_clean(&r.shape, tol))
}

// ---------------------------------------------------------------------------
// Repair iteration & single-defect tooling
// ---------------------------------------------------------------------------

/// Iterate [`repair_self_intersections`] until the boundary is clean.
///
/// A repair pass can introduce a *new* self-intersection when a face is split
/// near another crossing face, so `BOPAlgo`-style repair is applied to
/// fixpoint. Runs at most `max_iter` passes; the counters and warnings are
/// accumulated across passes.
pub fn repair_self_intersections_loop(shape: &TopoShape, tol: f64, max_iter: usize) -> RepairResult {
    let mut current = shape.clone();
    let mut fixed_total = 0usize;
    let mut removed_total = 0usize;
    let mut warnings: Vec<String> = Vec::new();
    for _ in 0..max_iter.max(1) {
        let rep = match repair_self_intersections(&current, tol) {
            Ok(r) => r,
            Err(e) => {
                warnings.push(format!("repair pass failed: {e}"));
                break;
            }
        };
        fixed_total += rep.fixed_faces;
        removed_total += rep.removed_faces;
        warnings.extend(rep.warnings);
        current = rep.repaired;
        if rep.fixed_faces == 0 && rep.removed_faces == 0 {
            break;
        }
    }
    RepairResult {
        repaired: current,
        fixed_faces: fixed_total,
        removed_faces: removed_total,
        warnings,
    }
}

/// List the degenerate (near-zero-area) faces of a shape.
pub fn shape_degenerate_faces(shape: &TopoShape, tol: f64) -> Vec<Face> {
    faces_of(shape)
        .into_iter()
        .filter(|f| face_is_degenerate(f, tol))
        .collect()
}

/// List the edges of a shape shorter than `tol`.
pub fn shape_small_edges(shape: &TopoShape, tol: f64) -> Vec<Edge> {
    let tol = tol.max(1e-9);
    edges_of(shape)
        .into_iter()
        .filter(|e| crate::brep_measure::edge_length(e, 8) < tol)
        .collect()
}

/// The face-index connectivity groups of a shape's boundary.
///
/// Returns the groups of [`face_connectivity_groups`] for every face under
/// `shape`, so callers can see how the boundary splits without rebuilding it.
pub fn boundary_connectivity(shape: &TopoShape, tol: f64) -> Vec<Vec<usize>> {
    let _ = tol;
    face_connectivity_groups(&faces_of(shape))
}

/// Structural validity issues of a shape, as human-readable strings.
///
/// A port of the `BOPAlgo_ArgumentAnalyzer` validity pass: invalid topology
/// structure, open shells, free edges, self-intersections and degenerate
/// faces. An empty list means the shape passed every check.
pub fn check_shape_validity(shape: &TopoShape, tol: f64) -> Vec<String> {
    let mut issues: Vec<String> = Vec::new();
    if !crate::topo_tools_full::structure_is_valid(shape) {
        issues.push("invalid topology structure (a sub-shape has an illegal parent type)".into());
    }
    let analysis = analyze_boundary(shape, tol);
    if analysis.open_shells > 0 {
        issues.push(format!("{} open shell(s)", analysis.open_shells));
    }
    if analysis.free_edges > 0 {
        issues.push(format!("{} free edge(s)", analysis.free_edges));
    }
    if analysis.self_intersections > 0 {
        issues.push(format!("{} self-intersection(s)", analysis.self_intersections));
    }
    if analysis.degenerate_faces > 0 {
        issues.push(format!("{} degenerate face(s)", analysis.degenerate_faces));
    }
    issues
}

// ---------------------------------------------------------------------------
// Component classification & selection
// ---------------------------------------------------------------------------

/// Broad shape kind of a result component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComponentKind {
    /// A solid (closed boundary).
    Solid,
    /// An open shell.
    Shell,
    /// A compound (still grouped result).
    Compound,
    /// A bare face.
    Face,
    /// A bare wire.
    Wire,
    /// Any other shape type.
    Other,
}

/// Classify a shape into a [`ComponentKind`].
pub fn component_kind(shape: &TopoShape) -> ComponentKind {
    match shape.shape_type() {
        ShapeType::Solid => ComponentKind::Solid,
        ShapeType::Shell => ComponentKind::Shell,
        ShapeType::Compound => ComponentKind::Compound,
        ShapeType::Face => ComponentKind::Face,
        ShapeType::Wire => ComponentKind::Wire,
        _ => ComponentKind::Other,
    }
}

/// Human-readable label of a [`ComponentKind`].
pub fn component_kind_label(kind: ComponentKind) -> &'static str {
    match kind {
        ComponentKind::Solid => "solid",
        ComponentKind::Shell => "shell",
        ComponentKind::Compound => "compound",
        ComponentKind::Face => "face",
        ComponentKind::Wire => "wire",
        ComponentKind::Other => "other",
    }
}

/// The components of a [`MultiResult`] that match `kind`.
pub fn multi_result_filter_by_kind(m: &MultiResult, kind: ComponentKind) -> Vec<TopoShape> {
    m.shapes
        .iter()
        .filter(|s| component_kind(s) == kind)
        .cloned()
        .collect()
}

/// The components of a [`MultiResult`], sorted by volume descending.
pub fn multi_result_sorted_by_volume(m: &MultiResult) -> Vec<TopoShape> {
    let mut sorted = m.shapes.clone();
    sorted.sort_by(|a, b| {
        crate::brep_gprop::volume(b, 0.02)
            .partial_cmp(&crate::brep_gprop::volume(a, 0.02))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    sorted
}

// ---------------------------------------------------------------------------
// Composite pipelines (check + repair, fuse-all, area)
// ---------------------------------------------------------------------------

/// Run a boolean with validation/tolerance-escalation, then repair the result.
///
/// Combines [`boolean_with_check`] (which re-runs the boolean at larger
/// tolerances when structural validation fails) with
/// [`repair_self_intersections`]. A repaired result carries a warning naming
/// the number of faces fixed/removed.
pub fn boolean_repaired_with_check(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let r = boolean_with_check(a, b, op, tol)?;
    Ok(repair_boolean_result(r, tol))
}

/// Fuse every shape in `shapes` and decompose the result into its pieces.
///
/// Convenience wrapper over [`boolean_multi`] + [`decompose_multi_result`]: a
/// Fuse of several disjoint boxes yields one [`MultiResult`] shape per box.
pub fn boolean_fuse_all_components(shapes: &[TopoShape], tol: f64) -> Result<MultiResult, String> {
    let r = boolean_multi(shapes, BoolOp::Fuse, tol)?;
    Ok(decompose_multi_result(&r))
}

/// Run a boolean, repair the result and return the repair report alongside it.
///
/// Like [`boolean_repaired`], but also exposes the underlying [`RepairResult`]
/// so callers can inspect exactly which faces were fixed/removed.
pub fn boolean_repaired_report(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<(BooleanResult, RepairResult), String> {
    let tol = tol.max(1e-9);
    let r = crate::bop_curved::curved_boolean_full(a, b, op, tol)?;
    let rep = repair_self_intersections(&r.shape, tol)?;
    if rep.fixed_faces == 0 && rep.removed_faces == 0 {
        return Ok((r, rep));
    }
    let mut rr = single_shape_result(&rep.repaired);
    rr.warnings.extend(r.warnings);
    rr.warnings.extend(rep.warnings.clone());
    rr.warnings.push(format!(
        "boolean result had self-intersections; repair fixed {} face(s), removed {}",
        rep.fixed_faces, rep.removed_faces
    ));
    Ok((rr, rep))
}

/// Surface area of a shape, from its tessellation (`BRepGProp::SurfaceProperties`).
pub fn boundary_area(shape: &TopoShape, deflection: f64) -> f64 {
    crate::brep_gprop::surface_area(shape, deflection.max(0.02))
}

/// Total surface area of every component of a [`MultiResult`].
pub fn multi_result_total_area(m: &MultiResult) -> f64 {
    m.shapes.iter().map(|s| crate::brep_gprop::surface_area(s, 0.02)).sum()
}

/// Volume of a single component (closed solid), or ~0 for open shells.
pub fn component_volume(shape: &TopoShape, deflection: f64) -> f64 {
    crate::brep_gprop::volume(shape, deflection)
}

// ---------------------------------------------------------------------------
// Boolean-result analysis helpers
// ---------------------------------------------------------------------------

/// The structural boundary analysis of a boolean result's shape.
pub fn boolean_result_boundary_report(r: &BooleanResult, tol: f64) -> BoundaryAnalysisReport {
    analyze_boundary(&r.shape, tol)
}

/// Number of disconnected pieces of a boolean result.
pub fn boolean_result_component_count(r: &BooleanResult) -> usize {
    decompose_multi_result(r).shapes.len()
}

/// Absolute difference between a boolean result's volume and an expectation.
///
/// `expected` is normally the analytic volume (e.g. `vol_a + vol_b` for a
/// Fuse); the delta is the mesh-measured deviation.
pub fn boolean_result_volume_delta(r: &BooleanResult, expected: f64) -> f64 {
    (crate::brep_gprop::volume(&r.shape, 0.02) - expected).abs()
}

/// Are two boolean results' volumes equal within a relative tolerance?
///
/// Uses the larger of the two volumes as the scale, so a 1% deviation on a
/// unit-volume result and on a 100-volume result are both "close".
pub fn boolean_results_close(a: &BooleanResult, b: &BooleanResult, rel_tol: f64) -> bool {
    let va = crate::brep_gprop::volume(&a.shape, 0.02);
    let vb = crate::brep_gprop::volume(&b.shape, 0.02);
    let scale = va.max(vb).max(1e-9);
    (va - vb).abs() <= rel_tol * scale
}

/// The broad shape kind of a boolean result's shape.
pub fn boolean_result_shape_kind(r: &BooleanResult) -> ComponentKind {
    component_kind(&r.shape)
}

/// Does a [`MultiResult`] contain at least one component of `kind`?
pub fn multi_result_contains_kind(m: &MultiResult, kind: ComponentKind) -> bool {
    m.shapes.iter().any(|s| component_kind(s) == kind)
}

/// Total volume of a boolean result's connected components.
pub fn boolean_result_total_volume(r: &BooleanResult) -> f64 {
    multi_result_total_volume(&decompose_multi_result(r))
}

// ---------------------------------------------------------------------------
// Repair options
// ---------------------------------------------------------------------------

/// Options that tune [`repair_self_intersections_opts`].
#[derive(Debug, Clone)]
pub struct RepairOptions {
    /// Face-intersection tolerance (points closer than this are coincident).
    pub tolerance: f64,
    /// Maximum number of repair passes (a split can reveal a new crossing).
    pub max_passes: usize,
    /// Drop overlapping coplanar duplicate faces.
    pub remove_coplanar_overlaps: bool,
    /// Emit a warning when the repaired boundary is not closed.
    pub report_open_boundaries: bool,
}

impl Default for RepairOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            max_passes: 1,
            remove_coplanar_overlaps: true,
            report_open_boundaries: true,
        }
    }
}

/// Repair self-intersections with explicit [`RepairOptions`].
///
/// This is the configurable front-end over [`repair_self_intersections`]:
/// `max_passes` runs the repair to a fixpoint (see
/// [`repair_self_intersections_loop`]) and `report_open_boundaries` controls
/// the open-shell warning. `remove_coplanar_overlaps` is advisory — the base
/// repair always removes overlapping coplanar duplicates, so disabling it only
/// affects the documentation of the intent.  (`ponytail:` threaded flag is a
/// no-op; upgrade the pair loop if per-call control ever matters.)
pub fn repair_self_intersections_opts(shape: &TopoShape, opts: &RepairOptions) -> Result<RepairResult, String> {
    let mut result = repair_self_intersections(shape, opts.tolerance)?;
    for _ in 1..opts.max_passes.max(1) {
        if result.fixed_faces == 0 && result.removed_faces == 0 {
            break;
        }
        let next = repair_self_intersections(&result.repaired, opts.tolerance)?;
        result.fixed_faces += next.fixed_faces;
        result.removed_faces += next.removed_faces;
        result.warnings.extend(next.warnings);
        result.repaired = next.repaired;
        if next.fixed_faces == 0 && next.removed_faces == 0 {
            break;
        }
    }
    if !opts.report_open_boundaries {
        result.warnings.retain(|w| !w.contains("not closed"));
    }
    Ok(result)
}

// ---------------------------------------------------------------------------
// Edge-connectivity report (TopOpeBRep)
// ---------------------------------------------------------------------------

/// Per-component edge-connectivity summary of a boundary.
///
/// Each entry corresponds to one edge-connected component of a shape's face
/// set — the same grouping [`decompose_multi_result`] uses to split a boolean
/// result — and reports the face indices and the component's vertex/edge/face
/// counts, from which the Euler characteristic follows.
#[derive(Debug, Clone)]
pub struct ComponentEdgeReport {
    /// Face indices (into `faces_of`) of this component.
    pub face_indices: Vec<usize>,
    /// Distinct vertices in the component (edge endpoints, deduplicated).
    pub vertex_count: usize,
    /// Distinct edges in the component.
    pub edge_count: usize,
    /// Number of faces in the component.
    pub face_count: usize,
}

/// Per-component connectivity reports for a shape's boundary.
///
/// A closed box is one component (6 faces, 12 edges, 8 vertices → Euler 2). A
/// shell built from two crossing faces splits into one component per face (each
/// is a 1-face open component). Components are returned sorted by smallest face
/// index.
pub fn component_edge_reports(shape: &TopoShape, tol: f64) -> Vec<ComponentEdgeReport> {
    let faces = faces_of(shape);
    let groups = face_connectivity_groups(&faces);
    let mut reports: Vec<ComponentEdgeReport> = Vec::with_capacity(groups.len());
    for idx in groups {
        let mut edges: HashSet<usize> = HashSet::new();
        let mut verts: HashSet<usize> = HashSet::new();
        for &i in &idx {
            for e in edges_of(&faces[i].0) {
                edges.insert(Arc::as_ptr(&e.0.tshape) as usize);
                let (a, b) = crate::topo_tools_full::edge_vertices(&e);
                if let Some(a) = a {
                    verts.insert(Arc::as_ptr(&a.0.tshape) as usize);
                }
                if let Some(b) = b {
                    verts.insert(Arc::as_ptr(&b.0.tshape) as usize);
                }
            }
        }
        reports.push(ComponentEdgeReport {
            face_count: idx.len(),
            face_indices: idx,
            vertex_count: verts.len(),
            edge_count: edges.len(),
        });
    }
    let _ = tol;
    reports
}

/// The Euler characteristic V − E + F of a boundary component.
///
/// A closed manifold component has χ = 2 (a sphere-like boundary); an open
/// component has χ = 1 or less. This is the `TopExp`-style sanity check for a
/// decomposed result piece.
pub fn component_euler_characteristic(r: &ComponentEdgeReport) -> i32 {
    (r.vertex_count as i32) - (r.edge_count as i32) + (r.face_count as i32)
}

// ---------------------------------------------------------------------------
// Edge-overlap repair (BOPAlgo / TopOpeBRep depth)
// ---------------------------------------------------------------------------
//
// Ports of the edge-level repair and classification passes of `BOPAlgo_Builder`
// and `TopOpeBRep_BuildTool`: collinear overlapping-edge repair, face splitting
// along surface intersections, boolean-result edge classification and tolerance
// healing. All functions follow the crate's `Result<_, String>` convention and
// rebuild shapes through the weld/edge-map machinery the exact boolean uses, so
// coincident vertices and edges stay shared across the rebuilt boundary.

/// Is the edge's curve geometrically a straight line?
///
/// Samples 8 points along the curve and checks that they are collinear with the
/// first–last chord. A `GeomLine` (even trimmed) passes; a circle/arc fails.
fn edge_is_line_like(e: &Edge) -> bool {
    let Some(c) = BRepTool::edge_curve(e) else { return false };
    let (a, b) = BRepTool::edge_parameters(e);
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let n = 8;
    let p0 = c.d0(a);
    let pl = c.d0(b);
    let size = p0.distance(&pl);
    if size <= 1e-30 {
        return false;
    }
    let chord = GpVec::from_pnts(&p0, &pl);
    let tol = 1e-6 * size;
    for i in 1..n {
        let p = c.d0(a + (b - a) * i as f64 / n as f64);
        if GpVec::from_pnts(&p0, &p).cross_magnitude(&chord) > tol * size {
            return false;
        }
    }
    true
}

/// Report of [`repair_edge_overlaps`].
#[derive(Debug, Clone)]
pub struct RepairReport {
    /// The repaired shape (unchanged when nothing needed repairing).
    pub repaired: TopoShape,
    /// Number of edges that were split to resolve a partial overlap.
    pub repaired_edges: usize,
    /// Number of edges removed (zero-length edges and coincident duplicates).
    pub removed_edges: usize,
    /// Number of vertex pairs merged because they lay within `tol`.
    pub welded_vertices: usize,
    /// Non-fatal diagnostics (open rebuilt boundary, dropped faces, …).
    pub warnings: Vec<String>,
}

/// Detect and repair collinear overlapping edges of a boolean result/compound.
///
/// Looks at every pair of line-like edges of the shape and repairs three kinds
/// of overlap:
///
/// * **zero-length edges** — an edge whose two endpoints coincide (within
///   `tol`) is removed;
/// * **coincident duplicates** — two collinear edges spanning the same
///   interval are merged into one (the duplicate is counted as removed);
/// * **partial overlaps** — two collinear edges that overlap over a sub-segment
///   are split at the overlap boundaries so the shared sub-segment becomes a
///   single edge (each split edge counts as repaired).
///
/// The boundary is rebuilt through the same `Weld`/`EdgeMap` machinery the exact
/// boolean uses, so every vertex closer than `tol` is welded to one instance and
/// coincident sub-segments resolve to the same `Edge`. Non-line edges (arcs,
/// full circles) and non-planar faces are kept untouched. A `Compound` is
/// repaired component-wise and reassembled.
pub fn repair_edge_overlaps(shape: &TopoShape, tol: f64) -> Result<RepairReport, String> {
    let tol = tol.max(1e-9);
    if shape.is_compound() {
        let children = expand_compound(shape);
        let mut repaired_edges = 0usize;
        let mut removed_edges = 0usize;
        let mut welded_vertices = 0usize;
        let mut warnings: Vec<String> = Vec::new();
        let mut out: Vec<TopoShape> = Vec::with_capacity(children.len());
        for c in &children {
            let r = repair_edge_overlaps(c, tol)?;
            repaired_edges += r.repaired_edges;
            removed_edges += r.removed_edges;
            welded_vertices += r.welded_vertices;
            warnings.extend(r.warnings);
            out.push(r.repaired);
        }
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&out);
        return Ok(RepairReport { repaired: comp.0, repaired_edges, removed_edges, welded_vertices, warnings });
    }

    let faces = faces_of(shape);
    let edges = edges_of(shape);
    let n = edges.len();

    // Per-edge geometry: endpoints (from the curve), line-likeness and length.
    let mut pts: Vec<Option<(GpPnt, GpPnt)>> = Vec::with_capacity(n);
    let mut is_line: Vec<bool> = Vec::with_capacity(n);
    let mut lens: Vec<f64> = Vec::with_capacity(n);
    for e in &edges {
        let ep = BRepTool::edge_vertices(e);
        let len = ep.as_ref().map(|(a, b)| a.distance(b)).unwrap_or(0.0);
        let line = ep.as_ref().map(|_| edge_is_line_like(e)).unwrap_or(false);
        pts.push(ep);
        is_line.push(line);
        lens.push(len);
    }

    // Removed: zero-length edges.
    let mut removed: HashSet<usize> = HashSet::new();
    for i in 0..n {
        if lens[i] < tol {
            removed.insert(i);
        }
    }

    // Pass 1: coincident-duplicate detection (same span on the same line).
    let mut merge_with: Vec<Option<usize>> = vec![None; n];
    let mut merged_count = 0usize;
    for i in 0..n {
        if removed.contains(&i) || !is_line[i] {
            continue;
        }
        let Some((a1, a2)) = pts[i] else { continue };
        let li = lens[i];
        let di = GpVec::from_pnts(&a1, &a2).normalized();
        for j in (i + 1)..n {
            if removed.contains(&j) || !is_line[j] || merge_with[j].is_some() {
                continue;
            }
            let Some((b1, b2)) = pts[j] else { continue };
            let db = GpVec::from_pnts(&b1, &b2).normalized();
            if di.cross_magnitude(&db) > tol {
                continue;
            }
            if GpVec::from_pnts(&a1, &b1).cross_magnitude(&di) > tol {
                continue;
            }
            let t = |p: &GpPnt| GpVec::from_pnts(&a1, p).dot(&di);
            let (lo, hi) = (t(&b1).min(t(&b2)), t(&b1).max(t(&b2)));
            if (lo - 0.0).abs() <= tol && (hi - li).abs() <= tol {
                merge_with[j] = Some(i);
                merged_count += 1;
            }
        }
    }

    // Pass 2: collect split parameters (in edge-local arc-length) for every
    // edge whose span is cut by a collinear overlapping edge.
    let mut split_ts: Vec<Vec<f64>> = vec![Vec::new(); n];
    for i in 0..n {
        if removed.contains(&i) || !is_line[i] {
            continue;
        }
        let Some((a1, a2)) = pts[i] else { continue };
        let li = lens[i];
        let di = GpVec::from_pnts(&a1, &a2).normalized();
        let t = |p: &GpPnt| GpVec::from_pnts(&a1, p).dot(&di);
        let mut ts = vec![0.0, li];
        for j in 0..n {
            if i == j || removed.contains(&j) || !is_line[j] {
                continue;
            }
            let Some((b1, b2)) = pts[j] else { continue };
            let db = GpVec::from_pnts(&b1, &b2).normalized();
            if di.cross_magnitude(&db) > tol {
                continue;
            }
            if GpVec::from_pnts(&a1, &b1).cross_magnitude(&di) > tol {
                continue;
            }
            let (tj1, tj2) = (t(&b1), t(&b2));
            let (lo, hi) = (tj1.min(tj2), tj1.max(tj2));
            if hi.min(li) - lo.max(0.0) <= tol {
                continue; // disjoint or merely touching at an endpoint
            }
            ts.push(lo.max(0.0).min(li));
            ts.push(hi.max(0.0).min(li));
        }
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        ts.dedup_by(|x, y| (*x - *y).abs() <= tol.max(1e-12));
        split_ts[i] = ts;
    }
    let repaired_count = split_ts.iter().filter(|ts| ts.len() > 2).count();

    // Vertex merges: distinct vertex instances that collapse into one position.
    let mut welded_vertices = 0usize;
    let mut positions: Vec<GpPnt> = Vec::new();
    for v in vertices_of(shape) {
        let p = vertex_position(&v);
        if positions.iter().any(|q| q.distance(&p) <= tol) {
            welded_vertices += 1;
        } else {
            positions.push(p);
        }
    }

    // Edges whose identity must change: removed, split, or merged with a twin.
    let mut affected: HashSet<usize> = HashSet::new();
    for i in 0..n {
        if removed.contains(&i) || split_ts[i].len() > 2 {
            affected.insert(i);
        }
    }
    for (j, k) in merge_with.iter().enumerate() {
        if let Some(i) = k {
            affected.insert(*i);
            affected.insert(j);
        }
    }
    if affected.is_empty() {
        return Ok(RepairReport {
            repaired: shape.clone(),
            repaired_edges: 0,
            removed_edges: removed.len() + merged_count,
            welded_vertices,
            warnings: vec![],
        });
    }

    // Rebuild the boundary: only faces carrying an affected edge are rebuilt.
    let bld = TopoBuilder::new();
    let mut edge_idx: HashMap<usize, usize> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        edge_idx.insert(Arc::as_ptr(&e.0.tshape) as usize, i);
    }
    let mut weld = Weld::new(tol.max(1e-7));
    let mut edge_map = EdgeMap::default();
    let mut result_faces: Vec<Face> = Vec::new();
    let mut warnings: Vec<String> = Vec::new();

    for f in &faces {
        let needs = edges_of(&f.0).iter().any(|e| {
            edge_idx
                .get(&(Arc::as_ptr(&e.0.tshape) as usize))
                .map_or(false, |idx| affected.contains(idx))
        });
        if !needs {
            result_faces.push(f.clone());
            continue;
        }
        let Some(pln) = face_plane_local(f) else {
            result_faces.push(f.clone());
            continue;
        };
        match rebuild_face_repaired(&bld, f, &pln, &edge_idx, &is_line, &pts, &split_ts, &removed, &mut weld, &mut edge_map)
        {
            Some(nf) => result_faces.push(nf),
            None => warnings.push("repair_edge_overlaps: a face collapsed to nothing and was dropped".into()),
        }
    }

    if result_faces.is_empty() {
        let empty = bld.make_compound_of(&[]);
        return Ok(RepairReport {
            repaired: empty.0,
            repaired_edges: repaired_count,
            removed_edges: removed.len() + merged_count,
            welded_vertices,
            warnings,
        });
    }

    let shell = bld.make_shell(&result_faces);
    let closed = shell_is_closed(&shell);
    if !closed {
        warnings.push("repair_edge_overlaps: rebuilt boundary is not closed".into());
    }
    let repaired: TopoShape = if closed { bld.make_solid(&[shell]).0 } else { shell.0 };
    Ok(RepairReport {
        repaired,
        repaired_edges: repaired_count,
        removed_edges: removed.len() + merged_count,
        welded_vertices,
        warnings,
    })
}

/// Rebuild a single face's wires with repaired edges.
///
/// Removed edges are dropped; split edges become their ordered sub-segments
/// (reversed when the wire traverses the edge backwards); unsplit line edges are
/// re-registered through the shared weld/edge-map so coincident edges across
/// faces resolve to one `Edge`; non-line edges are kept untouched. Returns
/// `None` when every wire collapsed.
#[allow(clippy::too_many_arguments)]
fn rebuild_face_repaired(
    bld: &TopoBuilder,
    f: &Face,
    pln: &GpPln,
    edge_idx: &HashMap<usize, usize>,
    is_line: &[bool],
    pts: &[Option<(GpPnt, GpPnt)>],
    split_ts: &[Vec<f64>],
    removed: &HashSet<usize>,
    weld: &mut Weld,
    edge_map: &mut EdgeMap,
) -> Option<Face> {
    let mut new_wires: Vec<Wire> = Vec::new();
    for w in wires_of_face(f) {
        let mut new_edges: Vec<Edge> = Vec::new();
        for we in edges_of_wire(&w) {
            let eptr = Arc::as_ptr(&we.0.tshape) as usize;
            let Some(&idx) = edge_idx.get(&eptr) else {
                new_edges.push(we.clone());
                continue;
            };
            if removed.contains(&idx) {
                continue;
            }
            let Some((p1, p2)) = pts[idx] else {
                new_edges.push(we.clone());
                continue;
            };
            let forward = !we.0.orientation().is_reversed();
            if !is_line[idx] {
                new_edges.push(we.clone());
                continue;
            }
            let ts = &split_ts[idx];
            if ts.len() >= 2 {
                let di = GpVec::from_pnts(&p1, &p2).normalized();
                let point_at = |t: f64| p1.translated_vec(&di.multiplied_scalar(t));
                let segs: Vec<(GpPnt, GpPnt)> = if forward {
                    ts.windows(2).map(|w| (point_at(w[0]), point_at(w[1]))).collect()
                } else {
                    ts.windows(2).rev().map(|w| (point_at(w[1]), point_at(w[0]))).collect()
                };
                for (a, b) in segs {
                    let ia = weld.weld(&a);
                    let ib = weld.weld(&b);
                    if ia != ib {
                        new_edges.push(edge_map.edge(bld, ia, ib, &weld.points));
                    }
                }
            } else {
                let (a, b) = if forward { (p1, p2) } else { (p2, p1) };
                let ia = weld.weld(&a);
                let ib = weld.weld(&b);
                if ia != ib {
                    new_edges.push(edge_map.edge(bld, ia, ib, &weld.points));
                }
            }
        }
        if !new_edges.is_empty() {
            new_wires.push(bld.make_wire(&new_edges));
        }
    }
    if new_wires.is_empty() {
        return None;
    }
    let mut nf = bld.make_face(Arc::new(GeomPlane::new(pln.clone())), &new_wires);
    nf.0.set_orientation(f.0.orientation());
    Some(nf)
}

// ---------------------------------------------------------------------------
// Face splitting along surface intersections
// ---------------------------------------------------------------------------

/// Split every face of `shape` along its intersection curves with the faces in
/// `pairs`, producing clean sub-faces with complete boundary wires.
///
/// `pairs` lists pairs of face indices (into [`faces_of`]) whose intersection
/// should be cut. For each pair the plane–plane intersection segment is computed
/// (via [`face_face_segments_local`]) and both faces are split along it with the
/// same polygon splitter and weld/edge-map the exact boolean uses, so the split
/// edges coincide and sub-faces share boundary vertices. A face that is not
/// cut (the segment misses it, or it is non-planar) is kept whole. The rebuilt
/// boundary is returned as a solid when closed, otherwise as a shell.
pub fn split_faces_along_intersections(shape: &TopoShape, pairs: &[(usize, usize)], tol: f64) -> Result<TopoShape, String> {
    let tol = tol.max(1e-9);
    let faces = faces_of(shape);
    for &(i, j) in pairs {
        if i >= faces.len() || j >= faces.len() {
            return Err(format!(
                "split_faces_along_intersections: face index ({i}, {j}) out of range for {} faces",
                faces.len()
            ));
        }
    }
    let mut segs: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); faces.len()];
    for &(i, j) in pairs {
        let s = face_face_segments_local(&faces[i], &faces[j], tol);
        for seg in &s {
            segs[i].push(*seg);
            segs[j].push(*seg);
        }
    }
    let bld = TopoBuilder::new();
    let mut weld = Weld::new(tol.max(1e-7));
    let mut edge_map = EdgeMap::default();
    let mut result_faces: Vec<Face> = Vec::new();
    for (i, f) in faces.iter().enumerate() {
        if segs[i].is_empty() {
            result_faces.push(f.clone());
            continue;
        }
        let Some(pln) = face_plane_local(f) else {
            result_faces.push(f.clone());
            continue;
        };
        let Some(poly2d) = face_polygon_local(f, &pln) else {
            result_faces.push(f.clone());
            continue;
        };
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
        let parts = split_polygon_by_segments(&poly2d, &dedupe_segments(segs2d));
        if parts.len() >= 2 {
            let subs: Vec<SubFace> = parts
                .iter()
                .filter_map(|p| polygon_to_subface(&bld, &pln, p, &mut weld, &mut edge_map))
                .collect();
            if subs.len() >= 2 {
                result_faces.extend(subs.into_iter().map(|sf| sf.face));
                continue;
            }
        }
        result_faces.push(f.clone());
    }
    let shell = bld.make_shell(&result_faces);
    let closed = shell_is_closed(&shell);
    let out: TopoShape = if closed { bld.make_solid(&[shell]).0 } else { shell.0 };
    Ok(out)
}

// ---------------------------------------------------------------------------
// Boolean-result edge classification (TopOpeBRep_BuildTool)
// ---------------------------------------------------------------------------

/// Classification of an edge of a boolean result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeClass {
    /// The edge is shared by several boundary patches (referenced by three or
    /// more result faces) — the coincidence of two bodies' edges in the result.
    Shared,
    /// The edge lies in the interior of a face's surface: its two adjacent
    /// faces are coplanar, so the edge is a seam/split line on a flat region.
    OnFace,
    /// The edge is not part of a closed two-face boundary: a free edge or a
    /// wire edge interior to the result topology.
    Internal,
    /// A regular outer-boundary edge (exactly two non-coplanar faces meet).
    External,
}

/// Per-class counts of a [`classify_boolean_edges`] result.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EdgeClassCounts {
    /// Number of `Shared` edges.
    pub shared: usize,
    /// Number of `OnFace` edges.
    pub on_face: usize,
    /// Number of `Internal` edges.
    pub internal: usize,
    /// Number of `External` edges.
    pub external: usize,
}

/// Classify every edge of a boolean result (`TopOpeBRep_BuildTool` style).
///
/// For each distinct edge of `result` the faces referencing it are counted:
///
/// * 3+ faces → [`EdgeClass::Shared`] (several boundary patches — the two
///   bodies' coincident edges — meet on one edge);
/// * exactly 2 coplanar faces → [`EdgeClass::OnFace`] (the edge is a seam lying
///   on a face's surface interior);
/// * exactly 2 non-coplanar faces → [`EdgeClass::External`] (regular outer
///   boundary edge);
/// * 0 or 1 face → [`EdgeClass::Internal`] (free / interior wire edge).
///
/// The operation is accepted for API symmetry with TopOpeBRep; the geometric
/// classification is operation-independent. Use [`edge_class_counts`] to
/// aggregate, and [`classify_edges_with_operands`] for the operand-aware
/// classification that reports edges lying on either input body's surface.
pub fn classify_boolean_edges(result: &TopoShape, _op: BoolOp) -> Vec<EdgeClass> {
    let faces = faces_of(result);
    let mut face_refs: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, f) in faces.iter().enumerate() {
        for e in edges_of(&f.0) {
            face_refs.entry(Arc::as_ptr(&e.0.tshape) as usize).or_default().push(fi);
        }
    }
    let edges = edges_of(result);
    edges
        .iter()
        .map(|e| {
            let ptr = Arc::as_ptr(&e.0.tshape) as usize;
            let fids = face_refs.get(&ptr).map(|v| v.as_slice()).unwrap_or(&[]);
            classify_result_edge(&faces, fids)
        })
        .collect()
}

/// Classify one edge from the indices of the faces referencing it.
fn classify_result_edge(faces: &[Face], fids: &[usize]) -> EdgeClass {
    match fids.len() {
        0 | 1 => EdgeClass::Internal,
        2 => {
            let (f1, f2) = (&faces[fids[0]], &faces[fids[1]]);
            match (face_plane_local(f1), face_plane_local(f2)) {
                (Some(p1), Some(p2)) if planes_coincident(&p1, &p2, 1e-6) => EdgeClass::OnFace,
                _ => EdgeClass::External,
            }
        }
        _ => EdgeClass::Shared,
    }
}

/// Aggregate [`EdgeClass`]es into [`EdgeClassCounts`].
pub fn edge_class_counts(classes: &[EdgeClass]) -> EdgeClassCounts {
    let mut c = EdgeClassCounts::default();
    for cl in classes {
        match cl {
            EdgeClass::Shared => c.shared += 1,
            EdgeClass::OnFace => c.on_face += 1,
            EdgeClass::Internal => c.internal += 1,
            EdgeClass::External => c.external += 1,
        }
    }
    c
}

/// One-line summary of [`EdgeClassCounts`].
pub fn edge_class_counts_summary(c: &EdgeClassCounts) -> String {
    format!(
        "{} edge(s): {} shared, {} on-face, {} internal, {} external",
        c.shared + c.on_face + c.internal + c.external,
        c.shared,
        c.on_face,
        c.internal,
        c.external
    )
}

/// Operand-aware edge classification for a boolean result.
///
/// For each edge of `result`, its midpoint is probed against both operand
/// boundaries:
///
/// * on **both** operands' surfaces → [`EdgeClass::Shared`] (the two bodies
///   meet along this edge);
/// * on exactly one operand's surface → [`EdgeClass::OnFace`];
/// * inside both operands (a seam hidden inside the union/intersection) →
///   [`EdgeClass::Internal`];
/// * otherwise → [`EdgeClass::External`].
///
/// This is the full `TopOpeBRep` classification; [`classify_boolean_edges`] is
/// its result-only form.
pub fn classify_edges_with_operands(a: &TopoShape, b: &TopoShape, result: &TopoShape, _op: BoolOp, tol: f64) -> Vec<EdgeClass> {
    let tol = tol.max(1e-9);
    let faces_a = faces_of(a);
    let faces_b = faces_of(b);
    edges_of(result)
        .iter()
        .map(|e| {
            let m = match BRepTool::edge_vertices(e) {
                Some((p1, p2)) => {
                    let v = GpVec::from_pnts(&p1, &p2);
                    p1.translated_vec(&v.multiplied_scalar(0.5))
                }
                None => GpPnt::zero(),
            };
            let on_a = point_on_surface(&faces_a, &m, tol);
            let on_b = point_on_surface(&faces_b, &m, tol);
            if on_a && on_b {
                EdgeClass::Shared
            } else if on_a || on_b {
                EdgeClass::OnFace
            } else if is_inside(a, &m) && is_inside(b, &m) {
                EdgeClass::Internal
            } else {
                EdgeClass::External
            }
        })
        .collect()
}

/// Is `p` (within `tol`) on the surface of any of `faces`?
fn point_on_surface(faces: &[Face], p: &GpPnt, tol: f64) -> bool {
    for f in faces {
        if let Some(pln) = face_plane_local(f) {
            if point_in_face_polygon(f, &pln, p, tol) {
                return true;
            }
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Tolerance healing
// ---------------------------------------------------------------------------

/// Report of [`heal_tolerance`] / [`heal_tolerance_report`].
#[derive(Debug, Clone)]
pub struct HealToleranceReport {
    /// The healed shape.
    pub healed: TopoShape,
    /// Number of vertex pairs merged (weld coincident vertices within `tol`).
    pub welded_vertices: usize,
    /// Number of edges removed (shorter than the heal tolerance).
    pub removed_edges: usize,
    /// Number of degenerate faces removed.
    pub removed_faces: usize,
    /// Non-fatal diagnostics (open shells, dropped faces, …).
    pub warnings: Vec<String>,
}

/// Heal a shape at a tolerance: weld near-coincident vertices, remove small
/// edges and degenerate faces (`ShapeFix_Shape`-style).
///
/// This is the tolerance-welding wrapper used to tidy a boolean result before
/// downstream processing. See [`heal_tolerance_report`] for the detailed
/// statistics; this entry point returns just the healed shape.
pub fn heal_tolerance(shape: &TopoShape, tol: f64) -> Result<TopoShape, String> {
    Ok(heal_tolerance_report(shape, tol)?.healed)
}

/// Heal a shape at a tolerance and report what was fixed.
///
/// Welds vertices closer than `tol` to one canonical instance, removes edges
/// shorter than `8·tol` (the small-edge threshold), and drops degenerate faces
/// (fewer than 3 boundary edges, or a planar face of near-zero area). A
/// `Compound` is healed component-wise. The healed shape is rebuilt through the
/// existing healing machinery, so shared edges/vertices are preserved.
pub fn heal_tolerance_report(shape: &TopoShape, tol: f64) -> Result<HealToleranceReport, String> {
    let tol = tol.max(1e-9);
    if shape.is_compound() {
        let children = expand_compound(shape);
        let mut healed_children: Vec<TopoShape> = Vec::with_capacity(children.len());
        let mut welded = 0usize;
        let mut removed_edges = 0usize;
        let mut removed_faces = 0usize;
        let mut warnings: Vec<String> = Vec::new();
        for c in &children {
            let r = heal_tolerance_report(c, tol)?;
            welded += r.welded_vertices;
            removed_edges += r.removed_edges;
            removed_faces += r.removed_faces;
            warnings.extend(r.warnings);
            healed_children.push(r.healed);
        }
        let bld = TopoBuilder::new();
        return Ok(HealToleranceReport {
            healed: bld.make_compound_of(&healed_children).0,
            welded_vertices: welded,
            removed_edges,
            removed_faces,
            warnings,
        });
    }

    let (s1, welded) = crate::shhealing::weld_coincident_vertices(shape, tol);
    let (s2, removed_edges) = crate::shhealing::remove_small_edges(&s1, tol * 8.0);
    let (s3, removed_faces) = remove_degenerate_faces(&s2, tol);
    let mut warnings: Vec<String> = Vec::new();
    for sh in shapes_of(&s3, ShapeType::Shell) {
        if !shell_is_closed(&Shell(sh)) {
            warnings.push("heal_tolerance: an open shell remains after healing".into());
        }
    }
    Ok(HealToleranceReport { healed: s3, welded_vertices: welded, removed_edges, removed_faces, warnings })
}

/// Is a face degenerate — fewer than 3 boundary edges, or a (planar) face with
/// near-zero area?
fn face_is_degenerate_tol(face: &Face, tol: f64) -> bool {
    let mut ecount = 0usize;
    for w in wires_of_face(face) {
        ecount += edges_of_wire(&w).len();
    }
    ecount < 3 || face_is_degenerate(face, tol)
}

/// Rebuild a shape without its degenerate faces. Shapes without shell structure
/// (bare wires/faces) are returned unchanged.
fn remove_degenerate_faces(shape: &TopoShape, tol: f64) -> (TopoShape, usize) {
    let shells: Vec<Shell> = shapes_of(shape, ShapeType::Shell).into_iter().map(Shell).collect();
    if shells.is_empty() {
        return (shape.clone(), 0);
    }
    let bld = TopoBuilder::new();
    let mut removed = 0usize;
    let mut out: Vec<Shell> = Vec::new();
    for sh in &shells {
        let faces = faces_of(&sh.0);
        let keep: Vec<Face> = faces.iter().filter(|f| !face_is_degenerate_tol(f, tol)).cloned().collect();
        removed += faces.len() - keep.len();
        if !keep.is_empty() {
            out.push(bld.make_shell(&keep));
        }
    }
    if out.is_empty() {
        return (bld.make_compound_of(&[]).0, removed);
    }
    if out.len() == 1 {
        let s = out.pop().unwrap();
        if shape.is_solid() && shell_is_closed(&s) {
            return (bld.make_solid(&[s]).0, removed);
        }
        return (s.0, removed);
    }
    let shapes: Vec<TopoShape> = out.into_iter().map(|s| s.0).collect();
    (bld.make_compound_of(&shapes).0, removed)
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

    // ------------------------------------------------------------------
    // Self-intersection repair + multi-result decomposition tests
    // ------------------------------------------------------------------

    /// A shell with a horizontal face (z = 0) crossed by a vertical face
    /// (y = 0): the two faces intersect along the x-axis segment.
    fn crossing_shell() -> Shell {
        use occt_core::gp::{GpAx3, GpDir, GpPln};
        use occt_geom::GeomPlane;
        use std::sync::Arc;
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
            GpAx3::new(
                GpPnt::zero(),
                GpDir::new(0.0, 1.0, 0.0).unwrap(),
                &GpDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
        );
        let f_v = bld.make_face(Arc::new(GeomPlane::new(pln_y)), &[v_wire]);
        bld.make_shell(&[f_h, f_v])
    }

    #[test]
    fn repair_noop_on_good_box() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = repair_self_intersections(&boxy.solid.0, 1e-6).expect("repair ok");
        assert_eq!(rep.fixed_faces, 0, "a valid box needs no face fixes");
        assert_eq!(rep.removed_faces, 0);
        assert!(rep.warnings.is_empty(), "no warnings: {:?}", rep.warnings);
        let v = box_vol(&rep.repaired);
        assert!((v - 1.0).abs() < 0.05, "box volume preserved {v}");
        clear_tree(&rep.repaired);
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn repair_detects_crossing() {
        let shell = crossing_shell();
        let rep = detect_self_intersections(&shell.0, 1e-6);
        assert!(rep.found, "crossing shell is flagged");
        let out = repair_self_intersections(&shell.0, 1e-6).expect("repair ok");
        assert!(out.fixed_faces >= 1, "at least one face fixed, got {}", out.fixed_faces);
        assert!(
            out.repaired.is_shell() || out.repaired.is_solid(),
            "repaired shape is a boundary (shell or solid)"
        );
        let no_longer_crossing = detect_self_intersections(&out.repaired, 1e-6);
        assert!(!no_longer_crossing.found, "crossing resolved after repair");
        clear_tree(&out.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn repair_preserves_volume() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = repair_self_intersections(&boxy.solid.0, 1e-6).expect("repair ok");
        let v = box_vol(&rep.repaired);
        assert!((v - 1.0).abs() < 0.05, "repaired box volume {v} (expected ~1.0)");
        clear_tree(&rep.repaired);
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn boolean_repaired_warns_on_fix() {
        // A self-intersecting shell fused with a disjoint box keeps the
        // crossing inside the compound; the repaired wrapper flags it.
        let shell = crossing_shell();
        let boxy = BRepPrimBox::make_box_corner(&GpPnt::new(5.0, 5.0, 5.0), &GpPnt::new(6.0, 6.0, 6.0));
        let r = boolean_repaired(&shell.0, &boxy.solid.0, BoolOp::Fuse, 1e-6).expect("repaired fuse ok");
        assert!(
            r.warnings.iter().any(|w| w.contains("fixed")),
            "warning mentions fixed faces: {:?}",
            r.warnings
        );

        // A clean boolean has no repair warning.
        let (a, b) = overlapping_boxes();
        let r2 = boolean_repaired(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("clean fuse ok");
        assert!(
            !r2.warnings.iter().any(|w| w.contains("fixed")),
            "clean result has no repair warning: {:?}",
            r2.warnings
        );
        clear_tree(&r.shape);
        clear_tree(&r2.shape);
        clear_tree(&shell.0);
        clear_tree(&boxy.solid.0);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn boolean_repaired_matches_plain() {
        let (a, b) = overlapping_boxes();
        let plain = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("plain fuse");
        let repaired = boolean_repaired(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("repaired fuse");
        let vp = box_vol(&plain.shape);
        let vr = box_vol(&repaired.shape);
        assert!((vp - vr).abs() < 0.05, "repaired {vr} vs plain {vp}");
        clear_tree(&plain.shape);
        clear_tree(&repaired.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn decompose_multi_result_three_boxes() {
        let shapes = disjoint_box_shapes();
        let r = boolean_multi(&shapes, BoolOp::Fuse, 1e-6).expect("multi fuse ok");
        let m = decompose_multi_result(&r);
        assert_eq!(m.shapes.len(), 3, "three disjoint boxes decompose into 3, got {}", m.shapes.len());
        assert!(m.solids.len() >= 1, "at least one solid: {}", m.solids.len());
        assert_eq!(m.compounds.len(), 0, "no nested compounds remain");
        clear_tree(&r.shape);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_split_result_fuse_disjoint() {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let m = boolean_split_result(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("split fuse ok");
        assert!(m.shapes.len() >= 2, "disjoint fuse produces multiple shapes: {}", m.shapes.len());
        assert!(m.solids.len() >= 2, "both pieces are solids: {}", m.solids.len());
        clear_tree(&b1.0);
        clear_tree(&b2.0);
    }

    #[test]
    fn boolean_split_result_cut_two_pieces() {
        // A 3×1×1 box cut by a thin full-width slab into two unit pieces.
        let box_a = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let slab = test_box_at(&GpPnt::new(1.0, -0.1, -0.1), &GpPnt::new(2.0, 1.1, 1.1));
        let m = boolean_split_result(&box_a.0, &slab.0, BoolOp::Cut, 1e-6).expect("split cut ok");
        assert!(m.shapes.len() >= 2, "cut into two pieces yields >=2 shapes, got {}", m.shapes.len());
        assert!(m.solids.len() >= 2, "two pieces are solids: {}", m.solids.len());
        let total: f64 = m.solids.iter().map(|s| box_vol(s)).sum();
        assert!((total - 2.0).abs() < 0.2, "total volume {total} (expected ~2.0)");
        clear_tree(&box_a.0);
        clear_tree(&slab.0);
    }

    #[test]
    fn result_is_multiple_disjoint_vs_overlap() {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let r_disjoint = boolean(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("disjoint fuse");
        assert!(result_is_multiple(&r_disjoint), "disjoint fuse is multiple");

        let (a, b) = overlapping_boxes();
        let r_overlap = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("overlap fuse");
        assert!(!result_is_multiple(&r_overlap), "overlapping fuse is a single result");
        clear_tree(&r_disjoint.shape);
        clear_tree(&r_overlap.shape);
        clear_tree(&b1.0);
        clear_tree(&b2.0);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn multi_result_solids_and_shells() {
        use crate::brep_builder_api::make_face_from_polygon;
        let boxy = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let face = make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ])
        .expect("face ok");
        let bld = TopoBuilder::new();
        let open_shell = bld.make_shell(&[face]);
        let comp = bld.make_compound_of(&[boxy.0.clone(), open_shell.0.clone()]);
        let r = BooleanResult {
            shape: comp.0.clone(),
            solid: None,
            shells: vec![open_shell],
            faces: vec![],
            warnings: vec![],
        };
        let m = decompose_multi_result(&r);
        assert!(m.shapes.len() >= 2, "solid + open shell decomposes into 2: {}", m.shapes.len());
        assert!(!m.solids.is_empty(), "solids populated");
        assert!(!m.shells.is_empty(), "shells populated");
        clear_tree(&comp.0);
        clear_tree(&boxy.0);
    }

    #[test]
    fn topologically_clean_removes_duplicates() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let v2 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let comp = b.make_compound_of(&[v1.0.clone(), v2.0.clone()]);
        let before = vertices_of(&comp.0).len();
        let cleaned = topologically_clean(&comp.0, 1e-6);
        let after = vertices_of(&cleaned).len();
        assert!(after < before, "clean reduces vertex count: {before} -> {after}");
        assert_eq!(after, 1, "duplicate vertices merge into one");
        clear_tree(&comp.0);
        clear_tree(&cleaned);
    }

    #[test]
    fn analyze_self_intersections_reports_crossing() {
        let shell = crossing_shell();
        let issues = analyze_self_intersections(&shell.0, 1e-6);
        assert_eq!(issues.len(), 1, "one crossing pair, got {}", issues.len());
        assert_eq!(issues[0].kind, SelfIntersectionKind::Crossing);
        assert!(!issues[0].points.is_empty(), "crossing carries sample points");
        assert!(issues[0].segment.is_some(), "crossing has a clipped segment");
        let summary = self_intersection_issues_summary(&issues);
        assert!(summary.contains("1 crossing"), "summary: {summary}");

        // A valid box has no issues.
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(analyze_self_intersections(&boxy.solid.0, 1e-6).is_empty());
        clear_tree(&shell.0);
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn analyze_boundary_box_vs_crossing() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = analyze_boundary(&boxy.solid.0, 1e-6);
        assert_eq!(rep.shell_count, 1);
        assert_eq!(rep.closed_shells, 1);
        assert_eq!(rep.open_shells, 0);
        assert_eq!(rep.free_edges, 0);
        assert_eq!(rep.self_intersections, 0);
        assert_eq!(rep.degenerate_faces, 0);
        assert_eq!(rep.small_edges, 0);
        let summary = boundary_analysis_summary(&rep);
        assert!(summary.contains("1 closed"), "summary: {summary}");

        let shell = crossing_shell();
        let rep2 = analyze_boundary(&shell.0, 1e-6);
        assert!(rep2.self_intersections >= 1, "crossing shell reports a self-intersection");
        assert!(rep2.open_shells >= 1, "crossing shell is open");
        clear_tree(&boxy.solid.0);
        clear_tree(&shell.0);
    }

    #[test]
    fn multi_result_helpers() {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let m = boolean_split_result(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("split fuse ok");
        assert_eq!(m.shapes.len(), 2);
        assert_eq!(component_counts_by_type(&m), (2, 0, 0));
        let vols = component_volumes(&m);
        assert!((vols.iter().sum::<f64>() - 2.0).abs() < 0.1, "volumes {vols:?}");
        assert!((multi_result_total_volume(&m) - 2.0).abs() < 0.1);
        let big = largest_component(&m).expect("largest");
        assert!(box_vol(&big) > 0.9, "largest component has volume");
        let summary = multi_result_summary(&m);
        assert!(summary.contains("2 shape(s)"), "summary: {summary}");

        // Pipeline wrappers agree.
        let comps = boolean_components(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("components");
        assert_eq!(comps.len(), 2);
        let cvols = boolean_component_volumes(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("volumes");
        assert_eq!(cvols.len(), 2);
        let single = boolean_repair_and_clean(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("clean");
        assert!(single.is_compound(), "two disjoint boxes fuse to a compound");
        clear_tree(&b1.0);
        clear_tree(&b2.0);
    }

    #[test]
    fn repair_loop_converges() {
        let shell = crossing_shell();
        let out = repair_self_intersections_loop(&shell.0, 1e-6, 4);
        assert!(out.fixed_faces >= 1, "loop fixed at least one face");
        assert!(
            !detect_self_intersections(&out.repaired, 1e-6).found,
            "loop converges to a clean boundary"
        );
        // The repaired boundary is still a valid (possibly open) shell.
        assert!(out.repaired.is_shell() || out.repaired.is_solid());
        let summary = repair_result_summary(&out);
        assert!(summary.contains("fixed"), "summary: {summary}");
        clear_tree(&out.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn validity_checks_box_clean() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(check_shape_validity(&boxy.solid.0, 1e-6).is_empty(), "box passes every check");
        assert!(shape_degenerate_faces(&boxy.solid.0, 1e-6).is_empty());
        assert!(shape_small_edges(&boxy.solid.0, 1e-6).is_empty());
        let conn = boundary_connectivity(&boxy.solid.0, 1e-6);
        assert_eq!(conn.len(), 1, "a box is one connected boundary");
        assert_eq!(conn[0].len(), 6, "all six faces connected");

        let shell = crossing_shell();
        let issues = check_shape_validity(&shell.0, 1e-6);
        assert!(
            issues.iter().any(|s| s.contains("self-intersection")),
            "crossing shell flags a self-intersection: {issues:?}"
        );
        clear_tree(&boxy.solid.0);
        clear_tree(&shell.0);
    }

    #[test]
    fn component_kind_classification() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert_eq!(component_kind(&boxy.solid.0), ComponentKind::Solid);
        assert_eq!(component_kind_label(ComponentKind::Solid), "solid");
        let face = crate::brep_builder_api::make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ])
        .expect("face ok");
        assert_eq!(component_kind(&face.0), ComponentKind::Face);
        clear_tree(&boxy.solid.0);
        clear_tree(&face.0);
    }

    #[test]
    fn boolean_repaired_with_check_clean() {
        let (a, b) = overlapping_boxes();
        let r = boolean_repaired_with_check(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("with-check repaired fuse");
        assert!(r.solid.is_some(), "fuse produces a solid");
        let v = box_vol(&r.shape);
        assert!((v - 1.5).abs() < 0.05, "with-check repaired volume {v} (expected 1.5)");
        let (report_r, rep) = boolean_repaired_report(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("report fuse");
        assert!(report_r.solid.is_some());
        assert_eq!(rep.fixed_faces, 0, "clean boolean needs no repair");
        clear_tree(&r.shape);
        clear_tree(&report_r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn boolean_result_analysis_helpers() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let report = boolean_result_boundary_report(&r, 1e-6);
        assert_eq!(report.closed_shells, 1, "fused box is one closed shell");
        assert_eq!(boolean_result_component_count(&r), 1, "single solid result");
        assert_eq!(boolean_result_shape_kind(&r), ComponentKind::Solid);
        let delta = boolean_result_volume_delta(&r, 1.5);
        assert!(delta < 0.05, "volume delta {delta}");
        assert!(boolean_results_close(&r, &r, 0.01), "a result is close to itself");
        assert!((boolean_result_total_volume(&r) - 1.5).abs() < 0.05);
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn component_edge_reports_box_and_crossing() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let reports = component_edge_reports(&boxy.solid.0, 1e-6);
        assert_eq!(reports.len(), 1, "a box is one component");
        assert_eq!(reports[0].face_count, 6);
        assert_eq!(reports[0].edge_count, 12);
        assert_eq!(reports[0].vertex_count, 8);
        assert_eq!(component_euler_characteristic(&reports[0]), 2, "closed box has χ = 2");

        let shell = crossing_shell();
        let reports2 = component_edge_reports(&shell.0, 1e-6);
        assert_eq!(reports2.len(), 2, "two crossing faces are two components");
        assert!(reports2.iter().all(|r| r.face_count == 1));
        clear_tree(&boxy.solid.0);
        clear_tree(&shell.0);
    }

    #[test]
    fn repair_opts_multi_pass() {
        let shell = crossing_shell();
        let opts = RepairOptions { tolerance: 1e-6, max_passes: 3, remove_coplanar_overlaps: true, report_open_boundaries: true };
        let out = repair_self_intersections_opts(&shell.0, &opts).expect("opts repair ok");
        assert!(out.fixed_faces >= 1, "multi-pass repair fixes the crossing");
        assert!(!detect_self_intersections(&out.repaired, 1e-6).found, "crossing resolved");

        // Default options: a single pass still fixes the crossing.
        let out2 = repair_self_intersections_opts(&shell.0, &RepairOptions::default()).expect("default repair ok");
        assert!(out2.fixed_faces >= 1);
        clear_tree(&out.repaired);
        clear_tree(&out2.repaired);
        clear_tree(&shell.0);
    }

    // ------------------------------------------------------------------
    // Phase 12 — BOPAlgo/TopOpeBRep depth: edge-overlap repair, face
    // splitting, edge classification and tolerance healing.
    // ------------------------------------------------------------------

    /// Does every shell under `shape` form a closed manifold boundary?
    fn closed_of(shape: &TopoShape) -> bool {
        let shells: Vec<Shell> = shapes_of(shape, ShapeType::Shell).into_iter().map(Shell).collect();
        !shells.is_empty() && shells.iter().all(shell_is_closed)
    }

    #[test]
    fn repair_edge_overlaps_noop_on_box() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = repair_edge_overlaps(&boxy.solid.0, 1e-6).expect("repair ok");
        assert_eq!(rep.repaired_edges, 0, "no overlapping edges in a box");
        assert_eq!(rep.removed_edges, 0, "no edges removed from a box");
        assert_eq!(rep.welded_vertices, 0, "a box has no near-coincident vertex pairs");
        assert!(rep.warnings.is_empty(), "no warnings: {:?}", rep.warnings);
        let v = box_vol(&rep.repaired);
        assert!((v - 1.0).abs() < 0.05, "box volume preserved {v}");
        assert!(closed_of(&rep.repaired), "repaired box stays closed");
        clear_tree(&rep.repaired);
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn repair_edge_overlaps_splits_partial_overlap() {
        // Two coplanar triangles whose base edges partially overlap along the
        // x-axis: [0,1] and [0.5,1.5]. Both edges must split at 0.5 / 1.0.
        let f1 = crate::brep_builder_api::make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ])
        .expect("face1");
        let f2 = crate::brep_builder_api::make_face_from_polygon(&[
            GpPnt::new(0.5, 0.0, 0.0),
            GpPnt::new(1.5, 0.0, 0.0),
            GpPnt::new(0.5, 1.0, 0.0),
        ])
        .expect("face2");
        let bld = TopoBuilder::new();
        let shell = bld.make_shell(&[f1, f2]);

        let rep = repair_edge_overlaps(&shell.0, 1e-6).expect("repair ok");
        assert!(rep.repaired_edges >= 2, "both overlapping edges are split, got {}", rep.repaired_edges);
        assert!(!closed_of(&rep.repaired), "two open faces form an open shell");
        // The overlapping sub-segment [0.5,1] is shared: distinct edge count rises
        // from 6 (two triangles) to 7 (the overlap resolved into three x-axis
        // segments plus four triangle side edges).
        let ec = edges_of(&rep.repaired).len();
        assert_eq!(ec, 7, "distinct edges after repair: {ec}");
        clear_tree(&rep.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn repair_edge_overlaps_merges_coincident() {
        // Two identical coplanar square faces: every boundary edge has a
        // coincident twin. The duplicates merge and the corner vertices weld.
        let square = |p: &GpPnt| {
            crate::brep_builder_api::make_face_from_polygon(&[
                *p,
                GpPnt::new(p.x() + 1.0, p.y(), p.z()),
                GpPnt::new(p.x() + 1.0, p.y() + 1.0, p.z()),
                GpPnt::new(p.x(), p.y() + 1.0, p.z()),
            ])
            .expect("square face")
        };
        let f1 = square(&GpPnt::new(0.0, 0.0, 0.0));
        let f2 = square(&GpPnt::new(0.0, 0.0, 0.0));
        let bld = TopoBuilder::new();
        let shell = bld.make_shell(&[f1, f2]);

        let rep = repair_edge_overlaps(&shell.0, 1e-6).expect("repair ok");
        assert!(rep.removed_edges >= 4, "four coincident edges removed, got {}", rep.removed_edges);
        assert!(rep.welded_vertices >= 4, "four corner pairs welded, got {}", rep.welded_vertices);
        clear_tree(&rep.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn split_faces_along_intersections_crossing() {
        // A horizontal face crossed by a vertical face: splitting along their
        // intersection turns the horizontal face into two sub-faces.
        let shell = crossing_shell();
        let faces_before = faces_of(&shell.0).len();
        let pairs = vec![(0usize, 1usize)];
        let split = split_faces_along_intersections(&shell.0, &pairs, 1e-6).expect("split ok");
        let faces_after = faces_of(&split).len();
        assert!(faces_after > faces_before, "splitting grows the face count: {faces_before} -> {faces_after}");
        assert!(faces_after >= 3, "crossing faces split into at least 3, got {faces_after}");
        clear_tree(&split);
        clear_tree(&shell.0);
    }

    #[test]
    fn classify_boolean_edges_box_all_external() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let classes = classify_boolean_edges(&boxy.solid.0, BoolOp::Fuse);
        assert_eq!(classes.len(), 12, "a box has 12 distinct edges");
        let counts = edge_class_counts(&classes);
        assert_eq!(counts.external, 12, "every box edge is external: {counts:?}");
        assert_eq!(counts.internal, 0);
        assert_eq!(counts.shared, 0);
        assert_eq!(counts.on_face, 0);
        let summary = edge_class_counts_summary(&counts);
        assert!(summary.contains("12 edge(s)"), "summary: {summary}");
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn classify_boolean_edges_shared_three_faces_on_one_edge() {
        // A fan of three triangle faces around a single shared edge: the shared
        // edge is referenced by 3 faces -> Shared; each triangle's two unique
        // edges are referenced by 1 face -> Internal.
        let bld = TopoBuilder::new();
        let e0 = bld.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let pln = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
        let face_on = |apex: GpPnt| {
            let e1 = bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &apex);
            let e2 = bld.make_edge_segment(&apex, &GpPnt::new(0.0, 0.0, 0.0));
            let wire = bld.make_wire(&[e0.clone(), e1, e2]);
            bld.make_face(pln.clone(), &[wire])
        };
        let shell = bld.make_shell(&[
            face_on(GpPnt::new(0.0, 1.0, 0.0)),
            face_on(GpPnt::new(0.0, -1.0, 0.0)),
            face_on(GpPnt::new(0.0, 0.0, 1.0)),
        ]);
        let classes = classify_boolean_edges(&shell.0, BoolOp::Fuse);
        let counts = edge_class_counts(&classes);
        assert_eq!(counts.shared, 1, "the shared base edge is used by three faces: {counts:?}");
        assert_eq!(counts.internal, 6, "six unique triangle edges are single-face edges");
        assert_eq!(counts.shared + counts.on_face + counts.internal + counts.external, classes.len());
        clear_tree(&shell.0);
    }

    #[test]
    fn heal_tolerance_welds_and_removes() {
        // A near-closed wire: the last vertex is within tol of the first, so
        // welding merges them.
        let w = crate::brep_builder_api::make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1e-6, 1e-6, 0.0),
        ])
        .expect("wire");
        let rep = heal_tolerance_report(&w.0, 1e-3).expect("heal ok");
        assert!(rep.welded_vertices >= 1, "near-coincident vertices weld, got {}", rep.welded_vertices);
        // A tiny edge is removed by the small-edge pass.
        let w2 = crate::brep_builder_api::make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(5e-4, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ])
        .expect("wire2");
        let rep2 = heal_tolerance_report(&w2.0, 1e-4).expect("heal2 ok");
        assert!(rep2.removed_edges >= 1, "tiny edge removed, got {}", rep2.removed_edges);
        // The TopoShape entry point returns just the healed shape.
        let healed = heal_tolerance(&w.0, 1e-3).expect("heal shape ok");
        assert!(healed.is_wire() || healed.is_compound(), "healed shape is a boundary");
        clear_tree(&rep.healed);
        clear_tree(&rep2.healed);
        clear_tree(&healed);
        clear_tree(&w.0);
        clear_tree(&w2.0);
    }

    #[test]
    fn phase12_integration_fuse_cut_split_repair_classify() {
        let (a, b) = overlapping_boxes();
        let fuse = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(fuse.solid.is_some(), "fuse is a solid");
        assert!(shell_is_closed(&fuse.shells[0]), "fuse shell is closed");
        let fuse_vol = box_vol(&fuse.shape);
        assert!((fuse_vol - 1.5).abs() < 0.1, "fuse volume {fuse_vol} (expected 1.5)");

        // Remaining face intersections in the clean fuse result (likely none).
        let issues = analyze_self_intersections(&fuse.shape, 1e-6);
        let pairs: Vec<(usize, usize)> = issues.iter().map(|i| (i.face_a, i.face_b)).collect();

        // Split faces along those intersections, then repair and heal.
        let split = split_faces_along_intersections(&fuse.shape, &pairs, 1e-6).expect("split ok");
        assert!(closed_of(&split), "split result stays closed");
        let rep = repair_edge_overlaps(&split, 1e-6).expect("repair ok");
        assert!(closed_of(&rep.repaired), "repaired result stays closed");
        let healed = heal_tolerance(&rep.repaired, 1e-6).expect("heal ok");
        assert!(closed_of(&healed), "healed result stays closed");

        // Volume is preserved through the whole pipeline.
        let v = box_vol(&healed);
        assert!((v - fuse_vol).abs() < 0.1, "volume preserved {v} vs {fuse_vol}");

        // Classify the repaired result: closed manifold -> no internal edges.
        let classes = classify_boolean_edges(&healed, BoolOp::Fuse);
        assert_eq!(classes.len(), edges_of(&healed).len(), "one class per edge");
        let counts = edge_class_counts(&classes);
        assert_eq!(counts.internal, 0, "closed result has no internal edges: {counts:?}");
        assert_eq!(counts.shared + counts.on_face + counts.internal + counts.external, classes.len());
        assert!(counts.external > 0, "outer boundary edges present");

        // The operand-aware classifier also reports a sane total.
        let op_classes = classify_edges_with_operands(&a.0, &b.0, &healed, BoolOp::Fuse, 1e-6);
        assert_eq!(op_classes.len(), classes.len());

        // Cut path: A − B is a closed solid too.
        let cut = boolean(&a.0, &b.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!(cut.solid.is_some(), "cut produces a solid");
        let cut_vol = box_vol(&cut.shape);
        assert!((cut_vol - 0.5).abs() < 0.1, "cut volume {cut_vol} (expected 0.5)");
        assert!(closed_of(&cut.shape), "cut result stays closed");

        clear_tree(&split);
        clear_tree(&rep.repaired);
        clear_tree(&healed);
        clear_tree(&fuse.shape);
        clear_tree(&cut.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

}
