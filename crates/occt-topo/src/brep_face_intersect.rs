//! Face–shape intersection. Port of `IntCurvesFace_ShapeIntersector` /
//! `IntCurvesFace_Intersector` (TKTopAlgo).
//!
//! Computes the intersection of a FACE with an arbitrary SHAPE (vertex, edge,
//! wire, face, shell/solid): isolated intersection POINTS (a vertex inside the
//! face, an edge crossing the face) and SEGMENTS / CURVES that lie on the face,
//! bounded by the face's edges.
//!
//! Gap analysis vs the existing modules:
//! - `face_face::face_face_intersection` already intersects two faces (exact
//!   plane–plane, sampled general) but does NOT trim the result to the face
//!   boundaries.
//! - `intpatch::surface_surface_intersection` and the `intersect_*` helpers
//!   return analytic + grid-traced intersection CURVES (surface-level,
//!   unbounded).
//! - `inttools::face_face_intersection_segments` already trims the plane–plane
//!   line to both faces' boundary polygons, and `inttools::edge_face_intersections`
//!   already finds an edge's crossing points on a face.
//!
//! This module therefore stays thin: it reuses `inttools` for the bounded
//! segments / points and falls back to `face_face` / `intpatch` for curved
//! intersections. `occt_geom::intana` is used where it is strictly better than
//! the grid tracer — plane∩sphere and sphere∩sphere return exact points /
//! circles instead of a sampled cloud.

use occt_core::elib::clib;
use occt_core::gp::{GpAx3, GpCirc, GpDir, GpPnt, GpSphere};
use occt_geom::intana::{quadric_quadric_plane_sphere, quadric_quadric_sphere_sphere, QuadricIntersection};

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::face_face::FaceIntersect;
use crate::intpatch::{plane_from_surface, sphere_params};
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of};

/// Default geometric tolerance for point/edge membership tests.
const POINT_TOL: f64 = 1e-7;

/// Result of intersecting a face with a shape.
#[derive(Debug, Clone)]
pub enum FaceShapeIntersect {
    /// No intersection.
    None,
    /// Isolated intersection points lying on the face.
    Points(Vec<GpPnt>),
    /// Line segments lying on the face (trimmed to the face's edges).
    Segments(Vec<[GpPnt; 2]>),
    /// Curved intersections, each a sampled polyline of points on the face.
    Curves(Vec<Vec<GpPnt>>),
}

/// Intersect `face` with `shape`, returning the intersection points / segments
/// / curves that lie on `face`, bounded by the face's edges. A non-positive
/// `tol` selects [`POINT_TOL`].
pub fn face_shape_intersection(face: &Face, shape: &TopoShape, tol: f64) -> FaceShapeIntersect {
    let tol = if tol > 0.0 { tol } else { POINT_TOL };
    match shape.shape_type() {
        ShapeType::Vertex => face_vertex(face, shape, tol),
        ShapeType::Edge => face_edge(face, shape, tol),
        ShapeType::Wire => match Wire::wrap(shape.clone()) {
            Some(w) => {
                let mut acc = FaceShapeIntersect::None;
                for e in edges_of_wire(&w) {
                    acc = merge(acc, face_edge(face, &e.0, tol));
                }
                acc
            }
            None => FaceShapeIntersect::None,
        },
        ShapeType::Face => match Face::wrap(shape.clone()) {
            Some(f2) => face_face(face, &f2, tol),
            None => FaceShapeIntersect::None,
        },
        // Shell / Solid / Compound: recurse over all sub-shapes.
        _ => face_container(face, shape, tol),
    }
}

/// A vertex that lies on the face.
fn face_vertex(face: &Face, shape: &TopoShape, tol: f64) -> FaceShapeIntersect {
    match Vertex::wrap(shape.clone()) {
        Some(v) => {
            let p = BRepTool::vertex_point(&v);
            if crate::inttools::point_on_face(face, &p, tol) {
                FaceShapeIntersect::Points(vec![p])
            } else {
                FaceShapeIntersect::None
            }
        }
        None => FaceShapeIntersect::None,
    }
}

/// Points where an edge crosses the face.
fn face_edge(face: &Face, shape: &TopoShape, tol: f64) -> FaceShapeIntersect {
    match Edge::wrap(shape.clone()) {
        Some(e) => {
            let pts: Vec<GpPnt> = crate::inttools::edge_face_intersections(&e, face, tol)
                .into_iter()
                .map(|(_, p)| p)
                .collect();
            if pts.is_empty() {
                FaceShapeIntersect::None
            } else {
                FaceShapeIntersect::Points(pts)
            }
        }
        None => FaceShapeIntersect::None,
    }
}

/// Intersect a face with a container shape (shell / solid / compound): collect
/// the per-sub-shape results.
fn face_container(face: &Face, shape: &TopoShape, tol: f64) -> FaceShapeIntersect {
    let mut acc = FaceShapeIntersect::None;
    for f2 in faces_of(shape) {
        acc = merge(acc, face_face(face, &f2, tol));
    }
    for e in edges_of(shape) {
        acc = merge(acc, face_edge(face, &e.0, tol));
    }
    for v in vertices_of(shape) {
        acc = merge(acc, face_vertex(face, &v.0, tol));
    }
    acc
}

/// Face–face intersection, trimmed to both faces' boundaries.
fn face_face(f1: &Face, f2: &Face, tol: f64) -> FaceShapeIntersect {
    // 1. Analytic quadric pairs via intana — exact where the grid tracer is
    //    approximate (tangent plane–sphere gives a clean point, secant gives a
    //    clean circle).
    if let Some(r) = analytic_quadric_face(f1, f2) {
        return r;
    }
    // 2. Bounded segments: inttools trims the plane–plane line to both face
    //    polygons (exact) and chains sampled points for non-planar faces.
    let segs = crate::inttools::face_face_intersection_segments(f1, f2, tol);
    if !segs.is_empty() {
        return FaceShapeIntersect::Segments(segs.into_iter().map(|(a, b)| [a, b]).collect());
    }
    // 3. Fallback: the sampled intersection point cloud.
    match crate::face_face::face_face_intersection(f1, f2, tol) {
        FaceIntersect::None => FaceShapeIntersect::None,
        // ponytail: coplanar overlapping faces should return the overlap
        // region; the polygon-overlap of two coplanar faces is not yet
        // computed. Return None until it is.
        FaceIntersect::Coplanar => FaceShapeIntersect::None,
        FaceIntersect::Line { .. } => FaceShapeIntersect::None, // trimmed line missed both polygons
        FaceIntersect::Curve { points } => {
            if points.is_empty() {
                FaceShapeIntersect::None
            } else {
                FaceShapeIntersect::Curves(vec![points])
            }
        }
    }
}

/// Analytic quadric intersection for the pairs this module can reconstruct
/// from a face's surface: plane∩sphere and sphere∩sphere.
fn analytic_quadric_face(f1: &Face, f2: &Face) -> Option<FaceShapeIntersect> {
    let s1 = BRepTool::face_surface_world(f1)?;
    let s2 = BRepTool::face_surface_world(f2)?;
    let pl1 = plane_from_surface(&*s1);
    let pl2 = plane_from_surface(&*s2);
    // A plane-like surface is never a sphere (numeric noise in `sphere_center`
    // can otherwise return a far-away circumsphere for a plane's near-coplanar
    // sample points). Only classify as a sphere when the surface is not planar.
    let sp1 = if pl1.is_some() { None } else { sphere_params(&*s1) };
    let sp2 = if pl2.is_some() { None } else { sphere_params(&*s2) };

    // plane ∩ sphere (either argument order).
    if let (Some(pl), Some((c, r))) = (pl1.as_ref(), sp2.as_ref()) {
        if let Some(sph) = build_sphere(*c, *r) {
            return Some(classify_analytic(quadric_quadric_plane_sphere(pl, &sph), f1, f2));
        }
    }
    if let (Some(pl), Some((c, r))) = (pl2.as_ref(), sp1.as_ref()) {
        if let Some(sph) = build_sphere(*c, *r) {
            return Some(classify_analytic(quadric_quadric_plane_sphere(pl, &sph), f1, f2));
        }
    }
    // sphere ∩ sphere.
    if let (Some((c1, r1)), Some((c2, r2))) = (sp1, sp2) {
        if let (Some(s1), Some(s2)) = (build_sphere(c1, r1), build_sphere(c2, r2)) {
            return Some(classify_analytic(
                quadric_quadric_sphere_sphere(&s1, &s2, 1e-7),
                f1,
                f2,
            ));
        }
    }
    None
}

/// Build a `gp_Sphere` from a solved center and radius (the face-surface
/// classifier does not retain the analytic type, so we rebuild it).
fn build_sphere(center: GpPnt, radius: f64) -> Option<GpSphere> {
    let ax3 = GpAx3::new(
        center,
        GpDir::new(0.0, 0.0, 1.0).ok()?,
        &GpDir::new(1.0, 0.0, 0.0).ok()?,
    )
    .ok()?;
    GpSphere::new(ax3, radius).ok()
}

/// Map an `intana::QuadricIntersection` onto the face-bounded result.
fn classify_analytic(inter: QuadricIntersection, f1: &Face, f2: &Face) -> FaceShapeIntersect {
    match inter {
        QuadricIntersection::None => FaceShapeIntersect::None,
        QuadricIntersection::Point(p) => {
            if crate::inttools::point_on_face(f1, &p, POINT_TOL)
                && crate::inttools::point_on_face(f2, &p, POINT_TOL)
            {
                FaceShapeIntersect::Points(vec![p])
            } else {
                FaceShapeIntersect::None
            }
        }
        QuadricIntersection::Circle(c) => {
            // A tangent pair can come back as a numerically-degenerate circle
            // (radius ~1e-8 from reconstruction noise); collapse it to its
            // center point so the API reports the tangency point.
            if c.radius() < 1e-6 {
                let p = c.location();
                if crate::inttools::point_on_face(f1, &p, POINT_TOL)
                    && crate::inttools::point_on_face(f2, &p, POINT_TOL)
                {
                    FaceShapeIntersect::Points(vec![p])
                } else {
                    FaceShapeIntersect::None
                }
            } else {
                circle_to_curves(&c, f1, f2)
            }
        }
        // Line (plane∩plane) is handled by inttools' polygon trimming; the
        // remaining conics (ellipse / parabola / hyperbola / two lines from
        // plane∩cone|cylinder) are not yet trimmed to the face boundary — the
        // grid tracer covers those pairs.
        _ => FaceShapeIntersect::None,
    }
}

/// Sample a `gp_Circ` and keep the points that lie on both faces.
fn circle_to_curves(c: &GpCirc, f1: &Face, f2: &Face) -> FaceShapeIntersect {
    let n = 48usize;
    let mut pts = Vec::new();
    for i in 0..n {
        let t = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let p = clib::circle_value(c, t);
        if crate::inttools::point_on_face(f1, &p, POINT_TOL)
            && crate::inttools::point_on_face(f2, &p, POINT_TOL)
        {
            pts.push(p);
        }
    }
    if pts.is_empty() {
        FaceShapeIntersect::None
    } else {
        FaceShapeIntersect::Curves(vec![pts])
    }
}

/// Combine two sub-results. When the kinds differ, the "richer" one wins
/// (segments > curves > points); a solid/compound intersection may contain all
/// three, which this single-variant enum cannot represent at once.
fn merge(a: FaceShapeIntersect, b: FaceShapeIntersect) -> FaceShapeIntersect {
    use FaceShapeIntersect::*;
    match (a, b) {
        (None, b) => b,
        (a, None) => a,
        (Points(mut v1), Points(v2)) => {
            v1.extend(v2);
            Points(dedupe_points(v1))
        }
        (Segments(mut v1), Segments(v2)) => {
            v1.extend(v2);
            Segments(v1)
        }
        (Curves(mut v1), Curves(v2)) => {
            v1.extend(v2);
            Curves(v1)
        }
        (Segments(v), _) | (_, Segments(v)) => Segments(v),
        (Curves(v), _) | (_, Curves(v)) => Curves(v),
    }
}

fn dedupe_points(pts: Vec<GpPnt>) -> Vec<GpPnt> {
    let mut out: Vec<GpPnt> = Vec::new();
    for p in pts {
        if !out.iter().any(|q| q.distance(&p) < 1e-6) {
            out.push(p);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpDir, GpPln};
    use occt_geom::{GeomSphere, Surface};

    use crate::builder::TopoBuilder;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn box_face_by_normal(bx: &crate::primitives::BRepPrimBox, n: (f64, f64, f64)) -> Face {
        crate::topo_tools_full::faces_of(&bx.solid.0)
            .into_iter()
            .find(|f| {
                crate::inttools::face_plane_from_face(f).map_or(false, |pln| {
                    let ax = pln.axis();
                    let d = ax.direction();
                    (d.x() - n.0).abs() < 1e-9 && (d.y() - n.1).abs() < 1e-9 && (d.z() - n.2).abs() < 1e-9
                })
            })
            .expect("box face with normal")
    }

    fn sphere_face(center: GpPnt, r: f64) -> Face {
        let b = TopoBuilder::new();
        let ax3 = GpAx3::new(
            center,
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let s: Arc<dyn Surface> = Arc::new(GeomSphere::new(
            occt_core::gp::GpSphere::new(ax3, r).unwrap(),
        ));
        b.make_face(s, &[])
    }

    fn plane_face(origin: GpPnt, normal: GpDir) -> Face {
        let b = TopoBuilder::new();
        b.make_face_plane(&GpPln::new(
            GpAx3::new(origin, normal, &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(),
        ))
    }

    #[test]
    fn adjacent_box_faces_share_edge_segment() {
        let bx = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let fx = box_face_by_normal(&bx, (1.0, 0.0, 0.0)); // x = 1
        let fy = box_face_by_normal(&bx, (0.0, 1.0, 0.0)); // y = 1
        match face_shape_intersection(&fx, &fy.0, 1e-9) {
            FaceShapeIntersect::Segments(segs) => {
                assert_eq!(segs.len(), 1, "segs: {segs:?}");
                let (p, q) = (segs[0][0], segs[0][1]);
                assert!((p.x() - 1.0).abs() < 1e-6 && (p.y() - 1.0).abs() < 1e-6, "p {p:?}");
                assert!((q.x() - 1.0).abs() < 1e-6 && (q.y() - 1.0).abs() < 1e-6, "q {q:?}");
                let zs = [p.z(), q.z()];
                assert!(zs.contains(&0.0) && zs.contains(&1.0), "z endpoints {zs:?}");
            }
            other => panic!("expected Segments, got {other:?}"),
        }
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn opposite_box_faces_miss() {
        let bx = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let fx = box_face_by_normal(&bx, (1.0, 0.0, 0.0)); // x = 1
        let fx0 = box_face_by_normal(&bx, (-1.0, 0.0, 0.0)); // x = 0
        assert!(matches!(face_shape_intersection(&fx, &fx0.0, 1e-9), FaceShapeIntersect::None));
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn tangent_plane_sphere_gives_point() {
        let sphere = sphere_face(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        let plane = plane_face(GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap());
        match face_shape_intersection(&sphere, &plane.0, 1e-6) {
            FaceShapeIntersect::Points(pts) => {
                assert_eq!(pts.len(), 1, "pts: {pts:?}");
                assert!(pts[0].distance(&GpPnt::new(0.0, 0.0, 1.0)) < 1e-6, "pt {:?}", pts[0]);
            }
            other => panic!("expected Points, got {other:?}"),
        }
        clear_tree(&sphere.0);
        clear_tree(&plane.0);
    }

    #[test]
    fn secant_plane_sphere_gives_circle() {
        let sphere = sphere_face(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        let plane = plane_face(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap());
        match face_shape_intersection(&sphere, &plane.0, 0.01) {
            FaceShapeIntersect::Curves(curves) => {
                assert_eq!(curves.len(), 1, "curves: {curves:?}");
                let pts = &curves[0];
                assert!(!pts.is_empty());
                for p in pts {
                    let d = p.distance(&GpPnt::new(0.0, 0.0, 0.0));
                    assert!((d - 1.0).abs() < 0.05, "point {p:?} at radius {d}");
                    assert!(p.z().abs() < 1e-9, "point {p:?} off the plane");
                }
            }
            other => panic!("expected Curves, got {other:?}"),
        }
        clear_tree(&sphere.0);
        clear_tree(&plane.0);
    }

    #[test]
    fn edge_crossing_face_gives_point() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard())); // z = 0
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));
        match face_shape_intersection(&face, &e.0, 1e-9) {
            FaceShapeIntersect::Points(pts) => {
                assert_eq!(pts.len(), 1, "pts: {pts:?}");
                assert!(pts[0].distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6, "pt {:?}", pts[0]);
            }
            other => panic!("expected Points, got {other:?}"),
        }
        clear_tree(&face.0);
        clear_tree(&e.0);
    }

    #[test]
    fn vertex_inside_face_and_outside() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard())); // unbounded z = 0
        let inside = b.make_vertex(GpPnt::new(1.0, 2.0, 0.0), 0.0);
        let outside = b.make_vertex(GpPnt::new(1.0, 2.0, 5.0), 0.0);
        assert!(matches!(
            face_shape_intersection(&face, &inside.0, 1e-7),
            FaceShapeIntersect::Points(ref pts) if pts.len() == 1
        ));
        assert!(matches!(face_shape_intersection(&face, &outside.0, 1e-7), FaceShapeIntersect::None));
        clear_tree(&face.0);
        clear_tree(&inside.0);
        clear_tree(&outside.0);
    }
}
