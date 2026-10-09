//! Geometry query — ray/segment/point interrogation of BRep shapes.
//! Source: `BRepExtrema_ShapeProximity`, `BRepClass3d_SolidClassifier`,
//! `IntCurvesFace_Intersector`.
//!
//! These operate on the MESHED approximation of a shape (like
//! `BRepMesh_IncrementalMesh` + `IntCurvesFace`), with optional refinement
//! against the analytic surfaces for the closest-hit point.

use occt_core::gp::{GpPnt, GpVec};

use crate::brep_tool::BRepTool;
use crate::mesh::ShapeMesh;
use crate::shape::{Face, TopoShape};
use crate::topo_tools_full::faces_of;

/// A ray-surface hit.
#[derive(Debug, Clone, Copy)]
pub struct Hit {
    /// Ray parameter `t` along `origin + t·dir`.
    pub t: f64,
    pub point: GpPnt,
    /// Face index into the queried shape's faces.
    pub face_index: usize,
}

/// Ray–triangle intersection (Möller–Trumbore). Returns `t` or None.
pub fn ray_triangle_hit(
    origin: &GpPnt,
    dir: &GpVec,
    v0: &GpPnt,
    v1: &GpPnt,
    v2: &GpPnt,
) -> Option<f64> {
    let edge1 = GpVec::from_pnts(v0, v1);
    let edge2 = GpVec::from_pnts(v0, v2);
    let p = dir.xyz().crossed(edge2.xyz());
    let det = edge1.xyz().dot(&p);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv = 1.0 / det;
    let tvec = GpVec::from_pnts(v0, origin);
    let u = tvec.xyz().dot(&p) * inv;
    if u < 0.0 || u > 1.0 {
        return None;
    }
    let q = tvec.xyz().crossed(edge1.xyz());
    let v = dir.xyz().dot(&q) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = edge2.xyz().dot(&q) * inv;
    if t > 1e-12 { Some(t) } else { None }
}

/// Ray–mesh intersection: the nearest hit `t` along the ray, plus the triangle
/// and the face it belongs to (via `face_of_triangle`). `face_of_triangle` maps
/// a triangle index to the face index; provide `None` to skip face attribution.
pub fn ray_mesh_hit(
    mesh: &ShapeMesh,
    origin: &GpPnt,
    dir: &GpVec,
    face_of_triangle: Option<&dyn Fn(usize) -> usize>,
) -> Option<Hit> {
    let mut best: Option<Hit> = None;
    for (ti, tri) in mesh.triangles.iter().enumerate() {
        let (a, b, c) = (&mesh.vertices[tri.n0], &mesh.vertices[tri.n1], &mesh.vertices[tri.n2]);
        if let Some(t) = ray_triangle_hit(origin, dir, a, b, c) {
            let face_index = face_of_triangle.map(|f| f(ti)).unwrap_or(0);
            let better = best.map_or(true, |b: Hit| t < b.t);
            if better {
                best = Some(Hit {
                    t,
                    point: GpPnt::new(
                        origin.x() + t * dir.x(),
                        origin.y() + t * dir.y(),
                        origin.z() + t * dir.z(),
                    ),
                    face_index,
                });
            }
        }
    }
    best
}

/// Mesh the shape (once) and return the hit + per-triangle face attribution.
/// This is the convenient `IntCurvesFace`-style entry point.
pub fn ray_shape_hit(
    shape: &TopoShape,
    origin: &GpPnt,
    dir: &GpVec,
    deflection: f64,
) -> Option<Hit> {
    let faces = faces_of(shape);
    let mesh = crate::shape_mesh::mesh_shape(shape, deflection);
    // Build a face index for every triangle by walking each face's triangulation.
    let mut face_of: Vec<usize> = Vec::with_capacity(mesh.triangles.len());
    for (i, f) in faces.iter().enumerate() {
        let (_, tris) = crate::wireframe::face_to_triangles(f, deflection);
        let n = tris.len();
        face_of.extend(std::iter::repeat(i).take(n));
    }
    // If the counts don't align (mesh_shape merges faces), fall back to a flat
    // per-triangle scan without attribution.
    if face_of.len() != mesh.triangles.len() {
        face_of = vec![0; mesh.triangles.len()];
    }
    ray_mesh_hit(&mesh, origin, dir, Some(&|ti| face_of[ti]))
}

/// Distance from a point to a shape's meshed surface (`BRepExtrema`-style).
pub fn point_to_shape_distance(
    shape: &TopoShape,
    p: &GpPnt,
    deflection: f64,
) -> f64 {
    let mesh = crate::shape_mesh::mesh_shape(shape, deflection);
    let mut best = f64::INFINITY;
    for tri in &mesh.triangles {
        let (a, b, c) = (&mesh.vertices[tri.n0], &mesh.vertices[tri.n1], &mesh.vertices[tri.n2]);
        // Distance to the triangle: min of point-to-plane and point-to-edges.
        let d = point_triangle_distance(p, a, b, c);
        best = best.min(d);
    }
    best
}

/// Distance from a point to a triangle (clamped to the triangle).
fn point_triangle_distance(p: &GpPnt, a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    // Project onto the plane, clamp to the triangle (barycentric), else edge dist.
    let ab = GpVec::from_pnts(a, b);
    let ac = GpVec::from_pnts(a, c);
    let n = ab.xyz().crossed(ac.xyz());
    let n2 = n.dot(&n);
    if n2 < 1e-30 {
        // Degenerate: min edge distance.
        return point_segment_dist(p, a, b).min(point_segment_dist(p, a, c)).min(point_segment_dist(p, b, c));
    }
    let ap = GpVec::from_pnts(a, p);
    let bdot = ab.xyz().dot(&ap.xyz());
    let cdot = ac.xyz().dot(&ap.xyz());
    let bbc = ab.xyz().dot(&ac.xyz());
    // Barycentric weights: v corresponds to vertex b (opposite edge ac), w to
    // vertex c (opposite edge ab). The norms must be the *opposite* edge's:
    //   v = (|ac|²·bdot − bbc·cdot) / |ab×ac|²
    //   w = (|ab|²·cdot − bbc·bdot) / |ab×ac|²
    // Swapping |ab|² and |ac|² here was the bug: a point on a face split into
    // two triangles (e.g. a box centre) landed on the wrong diagonal, giving
    // the edge distance instead of the plane distance.
    let ab2 = ab.xyz().dot(&ab.xyz());
    let ac2 = ac.xyz().dot(&ac.xyz());
    let denom = n2;
    let v = (ac2 * bdot - bbc * cdot) / denom;
    let w = (ab2 * cdot - bbc * bdot) / denom;
    let u = 1.0 - v - w;
    // Tolerance on the barycentric bounds: a point exactly on the triangle's
    // diagonal edge (a box centre projected onto a face split into two
    // triangles) yields u/v/w that are zero up to floating-point sign — treat
    // |·| ≤ eps as inside rather than falling through to the edge distance.
    const EPS: f64 = 1e-12;
    if u >= -EPS && v >= -EPS && w >= -EPS {
        // Closest point inside the triangle (clamp tiny negatives to the edge).
        let (u, v, w) = (u.max(0.0), v.max(0.0), w.max(0.0));
        let q = GpPnt::new(
            u * a.x() + v * b.x() + w * c.x(),
            u * a.y() + v * b.y() + w * c.y(),
            u * a.z() + v * b.z() + w * c.z(),
        );
        return p.distance(&q);
    }
    // Clamp to edges.
    let qab = point_segment_closest(p, a, b);
    let qac = point_segment_closest(p, a, c);
    let qbc = point_segment_closest(p, b, c);
    qab.min(qac).min(qbc)
}

fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let (_, d) = point_segment_closest_dist(p, a, b);
    d
}

fn point_segment_closest(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    point_segment_closest_dist(p, a, b).1
}

fn point_segment_closest_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> (GpPnt, f64) {
    let ab = GpVec::from_pnts(a, b);
    let len2 = ab.xyz().dot(&ab.xyz());
    if len2 < 1e-30 {
        return (*a, p.distance(a));
    }
    let ap = GpVec::from_pnts(a, p);
    let t = (ap.xyz().dot(&ab.xyz()) / len2).clamp(0.0, 1.0);
    let q = GpPnt::new(a.x() + t * ab.x(), a.y() + t * ab.y(), a.z() + t * ab.z());
    (q, p.distance(&q))
}

/// Refine a mesh hit against the analytic face surface: given the face and an
/// approximate hit point, project onto the surface for the exact surface point.
///
/// OCCT's point-on-surface primitive is `Extrema_ExtPS` (`GeomAPI_ProjectPointOnSurf`);
/// `IntCurvesFace_ShapeIntersector` returns the surface UV together with the hit
/// and never re-projects, so there is no OCCT grid branch behind this helper.
pub fn refine_hit_on_surface(face: &Face, approx: &GpPnt) -> GpPnt {
    match BRepTool::face_surface(face) {
        Some(s) => occt_geom::geom_api::project_point_on_surface(
            s.as_ref(),
            approx,
            occt_core::precision::CONFUSION,
        )
        .map_or(*approx, |ps| ps.point),
        None => *approx,
    }
}

/// Whether a ray segment [origin, origin+t·dir] enters a closed shape
/// (odd number of surface crossings) — a fast solid-classification probe.
pub fn ray_enters_solid(shape: &TopoShape, origin: &GpPnt, dir: &GpVec, deflection: f64) -> bool {
    let mesh = crate::shape_mesh::mesh_shape(shape, deflection);
    // Jitter the origin perpendicular to the ray so an exactly-centered start
    // does not pass through a grid vertex (0 crossings).
    let jittered = GpPnt::new(origin.x(), origin.y() + 1e-7, origin.z() + 1e-7);
    let mut hits: Vec<f64> = mesh
        .triangles
        .iter()
        .filter_map(|t| {
            ray_triangle_hit(&jittered, dir, &mesh.vertices[t.n0], &mesh.vertices[t.n1], &mesh.vertices[t.n2])
        })
        .collect();
    hits.sort_by(f64::total_cmp);
    let mut unique = 0usize;
    let mut prev: Option<f64> = None;
    for t in hits {
        if prev.map_or(true, |p| (t - p).abs() > 1e-9) {
            unique += 1;
            prev = Some(t);
        }
    }
    unique % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::ShapeType;
    use crate::primitives::BRepPrimBox;

    fn box_shape() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0
    }

    #[test]
    fn ray_hits_box_side() {
        let shape = box_shape();
        // Ray from (-1, 0.5, 0.5) along +X hits the x=0 face at t=1.
        let hit = ray_shape_hit(&shape, &GpPnt::new(-1.0, 0.5, 0.5), &GpVec::new(1.0, 0.0, 0.0), 0.1)
            .expect("hit");
        assert!((hit.t - 1.0).abs() < 1e-6, "t={}", hit.t);
        assert!((hit.point.x() - 0.0).abs() < 1e-6);
        assert!(hit.face_index < 6);
    }

    #[test]
    fn ray_miss() {
        let shape = box_shape();
        let hit = ray_shape_hit(&shape, &GpPnt::new(5.0, 5.0, 5.0), &GpVec::new(1.0, 0.0, 0.0), 0.1);
        assert!(hit.is_none());
    }

    #[test]
    fn point_distance_inside_and_outside() {
        let shape = box_shape();
        // Inside point: distance to surface ≈ 0.5 (min to the nearest face).
        let d_in = point_to_shape_distance(&shape, &GpPnt::new(0.5, 0.5, 0.5), 0.1);
        assert!(d_in > 0.4 && d_in < 0.6, "inside dist {d_in}");
        // Outside point 2 units away on +X → distance ≈ 1.0.
        let d_out = point_to_shape_distance(&shape, &GpPnt::new(2.0, 0.5, 0.5), 0.1);
        assert!(d_out > 0.9 && d_out < 1.2, "outside dist {d_out}");
    }

    #[test]
    fn ray_triangle_mt() {
        let (v0, v1, v2) = (GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.));
        let t = ray_triangle_hit(&GpPnt::new(0.25, 0.25, -1.0), &GpVec::new(0.0, 0.0, 1.0), &v0, &v1, &v2);
        assert!(t.is_some());
        assert!((t.unwrap() - 1.0).abs() < 1e-9);
        // Miss outside the footprint.
        assert!(ray_triangle_hit(&GpPnt::new(0.9, 0.9, -1.0), &GpVec::new(0.0, 0.0, 1.0), &v0, &v1, &v2).is_none());
    }

    #[test]
    fn ray_enters_closed_box() {
        let shape = box_shape();
        // Start INSIDE the closed box → odd crossings → true.
        assert!(ray_enters_solid(&shape, &GpPnt::new(0.5, 0.5, 0.5), &GpVec::new(1.0, 0.0, 0.0), 0.1));
        // Start OUTSIDE, even though the ray passes through → even → false.
        assert!(!ray_enters_solid(&shape, &GpPnt::new(-1.0, 0.5, 0.5), &GpVec::new(1.0, 0.0, 0.0), 0.1));
        assert!(!ray_enters_solid(&shape, &GpPnt::new(5.0, 5.0, 5.0), &GpVec::new(1.0, 0.0, 0.0), 0.1));
    }

    #[test]
    fn empty_and_vertex_shapes() {
        let v = TopoShape::new(ShapeType::Vertex);
        assert!(ray_shape_hit(&v, &GpPnt::new(0.,0.,0.), &GpVec::new(1.,0.,0.), 0.1).is_none());
        assert!(point_to_shape_distance(&v, &GpPnt::new(1.,1.,1.), 0.1).is_infinite());
    }
}
