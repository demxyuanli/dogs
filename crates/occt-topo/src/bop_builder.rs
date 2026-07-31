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

use std::collections::HashMap;
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
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, shapes_of, wires_of_face};

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

}
