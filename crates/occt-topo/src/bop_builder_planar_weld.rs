//! Weld / split / classify for the planar boolean arrangement.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::{GeomPlane, Surface};

use crate::abs::ShapeType;
use crate::brep_extrema::is_inside;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::inttools::{edge_edge_intersections, edge_face_intersections};
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::shell_check::{shell_invariants, shell_is_closed};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of,
    wires_of_face,
};

use crate::bop_builder_core::BoolOp;
use crate::bop_builder_planar_geom::*;
use crate::bop_builder_planar_trace::*;

// ---------------------------------------------------------------------------
// Welded edge sharing
// ---------------------------------------------------------------------------

/// Spatial-hash point welder: maps near-coincident points to one index so the
/// rebuilt edges of adjacent faces share vertices.
pub(crate) struct Weld {
    cell: f64,
    pub(crate) points: Vec<GpPnt>,
    grid: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl Weld {
    pub(crate) fn new(cell: f64) -> Self {
        Self { cell: cell.max(1e-9), points: Vec::new(), grid: HashMap::new() }
    }

    fn key(&self, p: &GpPnt) -> (i64, i64, i64) {
        (
            f64::floor(p.x() / self.cell) as i64,
            f64::floor(p.y() / self.cell) as i64,
            f64::floor(p.z() / self.cell) as i64,
        )
    }

    pub(crate) fn weld(&mut self, p: &GpPnt) -> usize {
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
pub(crate) struct EdgeMap {
    edges: HashMap<(usize, usize), Edge>,
}

impl EdgeMap {
    pub(crate) fn edge(&mut self, b: &TopoBuilder, i1: usize, i2: usize, pts: &[GpPnt]) -> Edge {
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

pub(crate) fn split_faces(
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
        // BOPAlgo_BuilderFace-style planar arrangement: the face's region(s)
        // come from tracing the subdivision of the boundary by the on-face
        // section segments (open lines split the face, closed loops carve
        // holes). A face with no sections traces back to its single region.
        let regions = trace_planar_regions(&poly2d, &segs2d);
        for r in regions {
            if let Some(sf) = region_to_subface(b, &pln, &r, weld, edge_map) {
                out.push(sf);
            }
        }
    }
    out
}

/// Build a [`SubFace`] from a traced planar region: the outer loop plus any
/// hole loops become separate wires of the face, and the region's interior
/// point is recorded for classification.
pub(crate) fn region_to_subface(
    b: &TopoBuilder,
    pln: &GpPln,
    r: &Region2d,
    weld: &mut Weld,
    edge_map: &mut EdgeMap,
) -> Option<SubFace> {
    if r.outer.len() < 3 || polygon_area2d(&r.outer).abs() < 1e-12 {
        return None;
    }
    let pts3d: Vec<GpPnt> = r.outer.iter().map(|p| plane_point_from_2d(pln, p)).collect();
    let mut wires: Vec<Wire> = Vec::new();
    // Outer wire.
    if let Some(w) = loop_wire(b, pln, &r.outer, weld, edge_map) {
        wires.push(w);
    } else {
        return None;
    }
    // Hole wires.
    for h in &r.holes {
        if let Some(w) = loop_wire(b, pln, h, weld, edge_map) {
            wires.push(w);
        }
    }
    let face = b.make_face(Arc::new(GeomPlane::new(pln.clone())), &wires);
    let interior = plane_point_from_2d(pln, &r.interior);
    Some(SubFace { face, poly3d: pts3d, plane: pln.clone(), interior })
}

/// Build a `Face` on `surface` from a traced 2D region (outer loop plus hole
/// loops), used by the BOPAlgo face-image stage's planar fast path.
pub(crate) fn region_to_face(r: &Region2d, pln: &GpPln, surface: Arc<dyn Surface>) -> Option<Face> {
    if r.outer.len() < 3 || polygon_area2d(&r.outer).abs() < 1e-12 {
        return None;
    }
    let b = TopoBuilder::new();
    let wire = |poly: &[GpPnt2d]| -> Option<Wire> {
        if poly.len() < 3 {
            return None;
        }
        let pts3d: Vec<GpPnt> = poly.iter().map(|p| plane_point_from_2d(pln, p)).collect();
        let edges: Vec<Edge> = (0..pts3d.len())
            .map(|i| b.make_edge_segment(&pts3d[i], &pts3d[(i + 1) % pts3d.len()]))
            .collect();
        Some(b.make_wire(&edges))
    };
    let outer = wire(&r.outer)?;
    let mut wires = vec![outer];
    for h in &r.holes {
        if let Some(w) = wire(h) {
            wires.push(w);
        }
    }
    Some(b.make_face(surface, &wires))
}

/// Split `face` into its planar regions with the 2-D arrangement
/// ([`trace_planar_regions`]), building one face per region on the face's own
/// surface. `on_edges` are the on-face (section) edges of the face.
///
/// This is the planar fast path used by `BOPAlgo_Builder::BuildSplitFaces`
/// (via `crate::bop_build_faces`): the arrangement splits the boundary edges
/// at section endpoints itself, so it does not depend on the pave-block
/// machinery having pre-split them.
pub(crate) fn split_face_planar_regions(face: &Face, on_edges: &[Edge]) -> Option<Vec<TopoShape>> {
    if on_edges.is_empty() {
        return None;
    }
    let pln = face_plane_local(face)?;
    let poly = face_polygon_local(face, &pln)?;
    let segs2d: Vec<(GpPnt2d, GpPnt2d)> = on_edges
        .iter()
        .filter_map(|e| {
            let (a, b) = edge_vertices(e);
            match (a, b) {
                (Some(av), Some(bv)) => {
                    let a2 = project_point_to_plane(&pln, &vertex_position(&av));
                    let b2 = project_point_to_plane(&pln, &vertex_position(&bv));
                    if a2.distance(&b2) > 1e-9 {
                        Some((a2, b2))
                    } else {
                        None
                    }
                }
                _ => None,
            }
        })
        .collect();
    let regions = trace_planar_regions(&poly, &segs2d);
    if regions.is_empty() {
        return None;
    }
    let surf = BRepTool::face_surface(face)?;
    let mut out = Vec::new();
    for r in &regions {
        if let Some(fc) = region_to_face(r, &pln, surf.clone()) {
            out.push(fc.0);
        }
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

/// Build the closed wire of one 2D loop on `pln`, welding vertices and
/// reusing shared edges through `edge_map`.
pub(crate) fn loop_wire(
    b: &TopoBuilder,
    pln: &GpPln,
    poly2d: &[GpPnt2d],
    weld: &mut Weld,
    edge_map: &mut EdgeMap,
) -> Option<Wire> {
    if poly2d.len() < 3 {
        return None;
    }
    let pts3d: Vec<GpPnt> = poly2d.iter().map(|p| plane_point_from_2d(pln, p)).collect();
    let idx: Vec<usize> = pts3d.iter().map(|p| weld.weld(p)).collect();
    let mut wire_edges: Vec<Edge> = Vec::new();
    for i in 0..pts3d.len() {
        let j = (i + 1) % pts3d.len();
        if idx[i] == idx[j] {
            continue;
        }
        wire_edges.push(edge_map.edge(b, idx[i], idx[j], &weld.points));
    }
    if wire_edges.len() < 3 {
        return None;
    }
    Some(b.make_wire(&wire_edges))
}

/// Rebuild a sub-face with the reversed plane (used for Cut's cut-through
/// faces so their mesh normal points away from the resulting solid). The
/// per-wire structure is preserved (a ring face keeps its hole loops).
pub(crate) fn flip_face(b: &TopoBuilder, sub: &SubFace) -> Face {
    let rpl = reversed_plane(&sub.plane);
    let wires: Vec<Wire> =
        wires_of_face(&sub.face).into_iter().map(|w| b.make_wire(&edges_of_wire(&w))).collect();
    b.make_face(Arc::new(GeomPlane::new(rpl)), &wires)
}

/// Rebuild a face with the reversed plane, keeping its wire edges.
pub(crate) fn flip_face_plane(b: &TopoBuilder, f: &Face, pln: &GpPln) -> Face {
    let rpl = reversed_plane(pln);
    let wires: Vec<Wire> =
        wires_of_face(f).into_iter().map(|w| b.make_wire(&edges_of_wire(&w))).collect();
    b.make_face(Arc::new(GeomPlane::new(rpl)), &wires)
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
/// A point inside the trimmed region of `f` on `pln` — the outer loop's
/// interior (not the centroid of all boundary points, which for a ring face
/// lands in the hole). Used as the probe origin for outward-orientation.
pub(crate) fn face_region_probe(f: &Face, pln: &GpPln) -> Option<GpPnt> {
    let mut loops: Vec<(f64, Vec<GpPnt2d>)> = Vec::new();
    for w in wires_of_face(f) {
        let mut edges: Vec<(GpPnt, GpPnt)> = Vec::new();
        for e in edges_of_wire(&w) {
            let (v1, v2) = edge_vertices(&e);
            if let (Some(a), Some(b)) = (v1, v2) {
                edges.push((vertex_position(&a), vertex_position(&b)));
            }
        }
        if edges.is_empty() {
            continue;
        }
        let mut pts: Vec<GpPnt> = vec![edges[0].0, edges[0].1];
        let mut used = vec![false; edges.len()];
        used[0] = true;
        let mut last = edges[0].1;
        for _ in 0..edges.len() {
            let mut found: Option<(usize, GpPnt)> = None;
            for (i, (a, b)) in edges.iter().enumerate() {
                if used[i] {
                    continue;
                }
                if a.distance(&last) < 1e-9 {
                    found = Some((i, *b));
                    break;
                }
                if b.distance(&last) < 1e-9 {
                    found = Some((i, *a));
                    break;
                }
            }
            match found {
                Some((i, nxt)) => {
                    used[i] = true;
                    pts.push(nxt);
                    last = nxt;
                }
                None => break,
            }
        }
        if pts.len() >= 3 {
            let pts2: Vec<GpPnt2d> = pts.iter().map(|p| project_point_to_plane(pln, p)).collect();
            let a = polygon_area2d(&pts2);
            if a.abs() > 1e-12 {
                loops.push((a, pts2));
            }
        }
    }
    if loops.is_empty() {
        return None;
    }
    let mut outer_i = 0;
    for (i, (a, _)) in loops.iter().enumerate() {
        if a.abs() > loops[outer_i].0.abs() {
            outer_i = i;
        }
    }
    let outer = loops[outer_i].1.clone();
    let holes: Vec<Vec<GpPnt2d>> = loops
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != outer_i)
        .map(|(_, (_, p))| p.clone())
        .collect();
    Some(plane_point_from_2d(pln, &region_interior2(&outer, &holes)))
}

pub(crate) fn orient_faces_outward(bld: &TopoBuilder, faces: &[Face]) -> Vec<Face> {
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
        let c = match face_region_probe(f, &pln) {
            Some(c) => c,
            None => {
                let pts = face_boundary_points(f).unwrap_or_default();
                polygon_centroid_3d(&pts)
            }
        };
        let probe = c.translated_vec(&n.multiplied_scalar(eps));
        if is_inside(&shape, &probe) {
            out.push(flip_face_plane(bld, f, &pln));
        } else {
            out.push(f.clone());
        }
    }
    out
}

/// Unify boundary edges across result faces (BOPDS-style): subdivide every
/// face's boundary at every vertex that lies on it, so faces sharing a
/// geometric boundary use the same edge partition and thus the same `Edge`
/// `TShape`. Without this, a face split by intersection lines (e.g. a box top
/// cut by a boss cylinder's chord lines) keeps subdivided boundary edges while
/// its neighbour keeps the original long edge, and the shell never closes.
///
/// Each face's boundary loop is rebuilt by CHAINING its wire edges at shared
/// endpoints (a shared edge's stored vertex order is the first creator's, so
/// `edge_vertices` order is not the traversal order), then refined at all
/// globally-welded vertices lying strictly inside a segment; every refined
/// adjacent pair maps to one canonical `Edge`.
pub(crate) fn unify_result_edges(bld: &TopoBuilder, faces: &[Face], tol: f64) -> Vec<Face> {
    let mut weld = Weld::new(tol.max(1e-7));
    // Per-face boundary loops (per wire) as welded vertex-index sequences.
    let mut face_loops: Vec<Vec<Vec<usize>>> = Vec::with_capacity(faces.len());
    let mut all_verts: Vec<usize> = Vec::new();
    for f in faces {
        let mut face_l: Vec<Vec<usize>> = Vec::new();
        for w in wires_of_face(f) {
            // Chain wire edges at shared endpoints into one vertex loop.
            let mut edges: Vec<(GpPnt, GpPnt)> = Vec::new();
            for e in edges_of_wire(&w) {
                let (v1, v2) = edge_vertices(&e);
                if let (Some(a), Some(b)) = (v1, v2) {
                    edges.push((vertex_position(&a), vertex_position(&b)));
                }
            }
            if edges.is_empty() {
                continue;
            }
            let mut loop_pts: Vec<GpPnt> = vec![edges[0].0, edges[0].1];
            let mut used = vec![false; edges.len()];
            used[0] = true;
            let mut last = edges[0].1;
            for _ in 0..edges.len() {
                let mut found: Option<(usize, GpPnt)> = None;
                for (i, (a, b)) in edges.iter().enumerate() {
                    if used[i] {
                        continue;
                    }
                    if a.distance(&last) < 1e-9 {
                        found = Some((i, *b));
                        break;
                    }
                    if b.distance(&last) < 1e-9 {
                        found = Some((i, *a));
                        break;
                    }
                }
                match found {
                    Some((i, next)) => {
                        used[i] = true;
                        loop_pts.push(next);
                        last = next;
                    }
                    None => break,
                }
            }
            if loop_pts.len() < 3 {
                continue;
            }
            if loop_pts.first().unwrap().distance(loop_pts.last().unwrap()) < 1e-9 {
                loop_pts.pop();
            }
            let loop_idx: Vec<usize> = loop_pts
                .iter()
                .map(|p| {
                    let i = weld.weld(p);
                    all_verts.push(i);
                    i
                })
                .collect();
            face_l.push(loop_idx);
        }
        face_loops.push(face_l);
    }

    let mut edge_map: HashMap<(usize, usize), Edge> = HashMap::new();
    // One shared vertex TShape per welded point, so edges sharing an endpoint
    // reference the same vertex (meaningful Euler characteristic).
    let mut vertex_map: HashMap<usize, Vertex> = HashMap::new();
    let mut shared_edge = |bld: &TopoBuilder,
                           edge_map: &mut HashMap<(usize, usize), Edge>,
                           vertex_map: &mut HashMap<usize, Vertex>,
                           a: usize,
                           b: usize|
     -> Edge {
        let key = (a.min(b), a.max(b));
        if let Some(e) = edge_map.get(&key) {
            return e.clone();
        }
        let va = vertex_map
            .entry(a)
            .or_insert_with(|| bld.make_vertex(weld.points[a], 0.0))
            .clone();
        let vb = vertex_map
            .entry(b)
            .or_insert_with(|| bld.make_vertex(weld.points[b], 0.0))
            .clone();
        let e = bld.make_edge_segment_with_vertices(&weld.points[a], &weld.points[b], &va, &vb);
        edge_map.insert(key, e.clone());
        e
    };
    let mut out: Vec<Face> = Vec::with_capacity(faces.len());
    for (fi, f) in faces.iter().enumerate() {
        let face_l = &face_loops[fi];
        let nw = wires_of_face(f).len();
        if face_l.is_empty() || face_l.len() != nw {
            out.push(f.clone());
            continue;
        }
        let mut wires: Vec<Wire> = Vec::new();
        let mut ok = true;
        for loop_idx in face_l {
            let n = loop_idx.len();
            if n < 3 {
                ok = false;
                break;
            }
            let mut refined: Vec<usize> = Vec::new();
            for k in 0..n {
                let a = loop_idx[k];
                let b = loop_idx[(k + 1) % n];
                refined.push(a);
                if a == b {
                    continue;
                }
                let pa = weld.points[a];
                let pb = weld.points[b];
                let ab = pb.coord.subtracted(&pa.coord);
                let len2 = ab.square_modulus();
                if len2 < 1e-24 {
                    continue;
                }
                let mut mid: Vec<(f64, usize)> = Vec::new();
                for &m in &all_verts {
                    if m == a || m == b {
                        continue;
                    }
                    let pm = weld.points[m];
                    let am = pm.coord.subtracted(&pa.coord);
                    if ab.crossed(&am).modulus() > 1e-7 * len2.sqrt() {
                        continue;
                    }
                    let t = am.dot(&ab) / len2;
                    if t > 1e-9 && t < 1.0 - 1e-9 {
                        mid.push((t, m));
                    }
                }
                mid.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
                for (_, m) in mid {
                    refined.push(m);
                }
            }
            let mut new_edges: Vec<Edge> = Vec::new();
            for k in 0..refined.len() {
                let (a, b) = (refined[k], refined[(k + 1) % refined.len()]);
                if a == b {
                    continue;
                }
                new_edges.push(shared_edge(bld, &mut edge_map, &mut vertex_map, a, b));
            }
            if new_edges.len() < 3 {
                ok = false;
                break;
            }
            wires.push(bld.make_wire(&new_edges));
        }
        if ok && !wires.is_empty() {
            if let Some(surf) = BRepTool::face_surface(f) {
                out.push(bld.make_face(surf, &wires));
            } else {
                out.push(f.clone());
            }
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
pub(crate) fn point_in_face_polygon(face: &Face, pln: &GpPln, p: &GpPnt, tol: f64) -> bool {
    let n = plane_normal(pln);
    let d = GpVec::from_pnts(&pln.location(), p).dot(&n.normalized()).abs();
    if d > tol {
        return false;
    }
    let Some(poly) = face_polygon_local(face, pln) else { return false };
    let p2 = project_point_to_plane(pln, p);
    point_in_polygon2d(&poly, &p2)
}

/// Hole-aware point-in-face: is `p` (on `pln`, within `tol`) inside the
/// trimmed region of `face`?
///
/// Mirrors `IntTools_FClass2d::Perform`'s region semantics (a point is `On`
/// within `tol` of a boundary ring, `In` when inside the outer ring and
/// outside every hole ring) but builds the rings in the face's own plane frame
/// from the wire boundaries, so the result is self-consistent with
/// `project_point_to_plane` and needs no pcurves. A point in a ring face's
/// hole is therefore `Out` (not "on" the face) while a ring-material point is
/// `In`.
pub(crate) fn point_in_face_holes(face: &Face, pln: &GpPln, p: &GpPnt, tol: f64) -> bool {
    let n = plane_normal(pln);
    let d = GpVec::from_pnts(&pln.location(), p).dot(&n.normalized()).abs();
    if d > tol {
        return false;
    }
    let p2 = project_point_to_plane(pln, p);
    // Per-wire boundary loops, chained into order and projected to `pln`'s
    // frame; the largest-|area| loop is the outer ring, the rest are holes.
    let mut loops: Vec<(f64, Vec<GpPnt2d>)> = Vec::new();
    for w in wires_of_face(face) {
        let mut pts: Vec<GpPnt> = Vec::new();
        let mut edges: Vec<(GpPnt, GpPnt)> = Vec::new();
        for e in edges_of_wire(&w) {
            let (v1, v2) = edge_vertices(&e);
            if let (Some(a), Some(b)) = (v1, v2) {
                edges.push((vertex_position(&a), vertex_position(&b)));
            }
        }
        if edges.is_empty() {
            continue;
        }
        pts.push(edges[0].0);
        pts.push(edges[0].1);
        let mut used = vec![false; edges.len()];
        used[0] = true;
        let mut last = edges[0].1;
        for _ in 0..edges.len() {
            let mut found: Option<(usize, GpPnt)> = None;
            for (i, (a, b)) in edges.iter().enumerate() {
                if used[i] {
                    continue;
                }
                if a.distance(&last) < 1e-9 {
                    found = Some((i, *b));
                    break;
                }
                if b.distance(&last) < 1e-9 {
                    found = Some((i, *a));
                    break;
                }
            }
            match found {
                Some((i, nxt)) => {
                    used[i] = true;
                    pts.push(nxt);
                    last = nxt;
                }
                None => break,
            }
        }
        if pts.len() >= 3 {
            let poly: Vec<GpPnt2d> = pts.iter().map(|q| project_point_to_plane(pln, q)).collect();
            let a = polygon_area2d(&poly);
            if a.abs() > 1e-12 {
                loops.push((a, poly));
            }
        }
    }
    if loops.is_empty() {
        return point_in_face_polygon(face, pln, p, tol);
    }
    let mut outer_i = 0;
    for (i, (a, _)) in loops.iter().enumerate() {
        if a.abs() > loops[outer_i].0.abs() {
            outer_i = i;
        }
    }
    // `On` within `tol` of a boundary ring (both outer and holes).
    let near = |poly: &[GpPnt2d]| {
        let n = poly.len();
        for i in 0..n {
            let a = poly[i];
            let b = poly[(i + 1) % n];
            let (abx, aby) = (b.x() - a.x(), b.y() - a.y());
            let (apx, apy) = (p2.x() - a.x(), p2.y() - a.y());
            let len2 = abx * abx + aby * aby;
            let t = if len2 < 1e-24 { 0.0 } else { ((apx * abx + apy * aby) / len2).clamp(0.0, 1.0) };
            let qx = a.x() + abx * t;
            let qy = a.y() + aby * t;
            if (p2.x() - qx).hypot(p2.y() - qy) < tol {
                return true;
            }
        }
        false
    };
    if near(&loops[outer_i].1) || loops.iter().enumerate().any(|(i, (_, l))| i != outer_i && near(l)) {
        return true;
    }
    if !point_in_polygon2d(&loops[outer_i].1, &p2) {
        return false;
    }
    for (i, (_, l)) in loops.iter().enumerate() {
        if i != outer_i && point_in_polygon2d(l, &p2) {
            return false;
        }
    }
    true
}

/// True when the sub-face's plane is coincident with a face of `other_faces`
/// and the sub-face's region-interior point lies inside that face's trimmed
/// region (hole-aware via [`point_in_face_holes`]).
pub(crate) fn face_on_other(plane: &GpPln, interior: &GpPnt, other_faces: &[Face], tol: f64) -> bool {
    for of in other_faces {
        let Some(opl) = face_plane_local(of) else { continue };
        if !planes_coincident(plane, &opl, tol) {
            continue;
        }
        if point_in_face_holes(of, &opl, interior, tol) {
            return true;
        }
    }
    false
}

pub(crate) fn classify_face(
    plane: &GpPln,
    interior: &GpPnt,
    other_faces: &[Face],
    other: &TopoShape,
    tol: f64,
) -> Class {
    if face_on_other(plane, interior, other_faces, tol) {
        let n = plane_normal(plane);
        let eps = tol.max(1e-6);
        let p_in = interior.translated_vec(&n.reversed().multiplied_scalar(eps));
        if is_inside(other, &p_in) {
            Class::OnSame
        } else {
            Class::OnOpp
        }
    } else if is_inside(other, interior) {
        Class::In
    } else {
        Class::Out
    }
}

pub(crate) fn select(c: Class, from_a: bool, op: BoolOp) -> bool {
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
