//! Plane / polygon helpers for the planar boolean arrangement.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::{GeomPlane, Surface};

use crate::abs::ShapeType;
use crate::bop_builder_core::{
    disjoint_result, empty_result, validate, BoolOp, BooleanResult,
};
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

// ---------------------------------------------------------------------------
// Sub-face bookkeeping
// ---------------------------------------------------------------------------

/// A sub-face produced by splitting a source face's polygon.
#[derive(Debug, Clone)]
pub(crate) struct SubFace {
    pub(crate) face: Face,
    /// 3D polygon vertices (outer loop, in boundary order).
    pub(crate) poly3d: Vec<GpPnt>,
    /// The face's plane (normal = source solid's outward normal).
    pub(crate) plane: GpPln,
    /// A 3D point guaranteed inside the sub-face's region (avoids the outer
    /// centroid landing in a hole for ring faces).
    pub(crate) interior: GpPnt,
}

/// Classification of a face relative to the other solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Class {
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

pub(crate) fn plane_normal(pln: &GpPln) -> GpVec {
    GpVec::from_xyz(pln.axis().direction().xyz())
}

/// Extract a `GpPln` from a planar face's registered surface. `None` when the
/// surface is not geometrically a plane (triggers the voxel fallback).
pub(crate) fn face_plane_local(face: &Face) -> Option<GpPln> {
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
pub(crate) fn project_point_to_plane(pln: &GpPln, p: &GpPnt) -> GpPnt2d {
    let origin = pln.location();
    let x = GpVec::from_xyz(pln.position().x_direction().xyz());
    let y = GpVec::from_xyz(pln.position().y_direction().xyz());
    let v = GpVec::from_pnts(&origin, p);
    GpPnt2d::new(v.dot(&x), v.dot(&y))
}

/// Map a 2D point in the plane's (x_dir, y_dir) basis back to 3D.
pub(crate) fn plane_point_from_2d(pln: &GpPln, p: &GpPnt2d) -> GpPnt {
    let origin = pln.location();
    let x = GpVec::from_xyz(pln.position().x_direction().xyz());
    let y = GpVec::from_xyz(pln.position().y_direction().xyz());
    origin
        .translated_vec(&x.multiplied_scalar(p.x()))
        .translated_vec(&y.multiplied_scalar(p.y()))
}

/// The plane with the normal reversed (keeps the origin and x-direction).
pub(crate) fn reversed_plane(pln: &GpPln) -> GpPln {
    let ax = pln.position();
    let n = ax.direction().reversed();
    GpPln::new(GpAx3::new(ax.location(), n, ax.x_direction()).unwrap_or(ax))
}

pub(crate) fn planes_coincident(a: &GpPln, b: &GpPln, tol: f64) -> bool {
    let na = plane_normal(a);
    let nb = plane_normal(b);
    if na.cross_magnitude(&nb) > tol {
        return false;
    }
    let d = GpVec::from_pnts(&a.location(), &b.location()).dot(&na.normalized()).abs();
    d <= tol
}

pub(crate) fn polygon_centroid_3d(pts: &[GpPnt]) -> GpPnt {
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
pub(crate) fn face_boundary_points(face: &Face) -> Option<Vec<GpPnt>> {
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
pub(crate) fn face_polygon_local(face: &Face, pln: &GpPln) -> Option<Vec<GpPnt2d>> {
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
pub(crate) fn face_line_intervals(pln: &GpPln, origin: &GpPnt, dir: &GpVec, poly: &[GpPnt2d]) -> Vec<(f64, f64)> {
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
pub(crate) fn face_face_segments_local(f1: &Face, f2: &Face, tol: f64) -> Vec<(GpPnt, GpPnt)> {
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
