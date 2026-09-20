//! Phase 4 module: hlr — hidden-line removal / projection.
//!
//! **UNPORTED (audit A14)**: not a translation of `HLRBRep_*` (TKHlr). Only the
//! common case of a lightweight projection/wireframe is implemented here; the
//! OCCT hidden-line algorithm (`HLRBRep_Algo`/`HLRBRep_HLRToShape`) has no port.
//!
//! A lightweight port of OCCT's `HLRBRep` for the common case of a
//! triangulated shape under orthographic projection. Vertices are projected
//! onto the plane perpendicular to a view direction, triangles are
//! backface-culled, then depth-sorted (painter's algorithm). `visible_edges`
//! walks the triangles front-to-back and reports each edge that is not covered
//! by a nearer triangle.
//!
//! Approximation: an edge is considered hidden only when a nearer triangle's
//! projected area contains *both* of its endpoints. A full HLR (OCCT's
//! `HLRBRep`) would also handle partial segment/area overlaps, edge
//! intersections and silhouette outlines. This is adequate for clean convex
//! shapes viewed along a fixed axis.

use occt_core::gp::{GpVec, GpXyz};

use crate::mesh::ShapeMesh;

/// A vertex projected onto the view plane. `z` is the depth along the view
/// direction; `x`/`y` are the view-plane coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProjectedPoint {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// A triangle in the projected mesh. `a`/`b`/`c` index into
/// [`ProjectedMesh::points`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProjectedTriangle {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    /// Average depth (view-axis coordinate) of the three vertices.
    pub depth: f64,
    /// View-axis component of the triangle normal (`normal · view_dir`).
    /// Negative ⇒ facing the viewer; positive ⇒ back-facing.
    pub normal_z: f64,
}

/// A mesh projected onto the plane perpendicular to the view direction.
#[derive(Debug, Clone)]
pub struct ProjectedMesh {
    pub points: Vec<ProjectedPoint>,
    pub triangles: Vec<ProjectedTriangle>,
}

/// Project `mesh` onto the plane perpendicular to `view_dir`.
///
/// An orthonormal basis `(u, v, n)` is built with `n = normalize(view_dir)`;
/// each vertex maps to `(x, y, z) = (p·u, p·v, p·n)`. Triangles whose normal
/// points away from the viewer (`normal · view_dir > 0`) are backface-culled;
/// `depth` is the average of the three projected `z` values.
pub fn orthographic_project(mesh: &ShapeMesh, view_dir: &GpVec) -> ProjectedMesh {
    let n = view_dir.normalized().coord;
    let (u, v) = orthonormal_basis(n);

    let points: Vec<ProjectedPoint> = mesh
        .vertices
        .iter()
        .map(|p| ProjectedPoint {
            x: p.coord.dot(&u),
            y: p.coord.dot(&v),
            z: p.coord.dot(&n),
        })
        .collect();

    let triangles = mesh
        .triangles
        .iter()
        .filter_map(|t| {
            let pa = mesh.vertices[t.n0].coord;
            let pb = mesh.vertices[t.n1].coord;
            let pc = mesh.vertices[t.n2].coord;
            let normal = pb.subtracted(&pa).crossed(&pc.subtracted(&pa));
            let normal_z = normal.dot(&n);
            if normal_z > 0.0 {
                return None; // back-facing → cull
            }
            let depth = (points[t.n0].z + points[t.n1].z + points[t.n2].z) / 3.0;
            Some(ProjectedTriangle { a: t.n0, b: t.n1, c: t.n2, depth, normal_z })
        })
        .collect();

    ProjectedMesh { points, triangles }
}

/// Build an orthonormal basis `(u, v, n)` with `u × v = n`.
fn orthonormal_basis(n: GpXyz) -> (GpXyz, GpXyz) {
    if n.square_modulus() < 1e-30 {
        return (GpXyz::new(1.0, 0.0, 0.0), GpXyz::new(0.0, 1.0, 0.0));
    }
    let n = n.normalized();
    // Pick a reference axis not parallel to n, then u = ref × n, v = n × u.
    let reference = if n.z.abs() < 0.9 {
        GpXyz::new(0.0, 0.0, 1.0)
    } else {
        GpXyz::new(1.0, 0.0, 0.0)
    };
    let u = reference.crossed(&n).normalized();
    let v = n.crossed(&u);
    (u, v)
}

/// Triangle indices sorted far-to-near by depth (painter's algorithm).
pub fn painter_sort(mesh: &ProjectedMesh) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..mesh.triangles.len()).collect();
    indices.sort_by(|&a, &b| {
        mesh.triangles[b]
            .depth
            .partial_cmp(&mesh.triangles[a].depth)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    indices
}

/// Visible 2D edge segments of the projected mesh.
///
/// Walks triangles front-to-back (the reverse of [`painter_sort`]). A triangle
/// edge is emitted unless a nearer, already-visited triangle's projected area
/// contains both of its endpoints. `tol` is the tolerance for the
/// point-in-triangle test and for discarding degenerate (zero-length) edges.
///
/// ponytail: conservative segment test — a partial occlusion (an edge crosses
/// a nearer triangle but leaves it) is not detected; upgrade to a
/// segment-clipping pass if partial overlaps matter.
pub fn visible_edges(mesh: &ProjectedMesh, tol: f64) -> Vec<((f64, f64), (f64, f64))> {
    let order = painter_sort(mesh); // far → near
    let mut nearer: Vec<usize> = Vec::new();
    let mut out = Vec::new();

    for &ti in order.iter().rev() {
        let tri = &mesh.triangles[ti];
        for (i, j) in [(tri.a, tri.b), (tri.b, tri.c), (tri.c, tri.a)] {
            let p = &mesh.points[i];
            let q = &mesh.points[j];
            if (p.x - q.x).hypot(p.y - q.y) <= tol {
                continue; // degenerate edge (edge-on triangle)
            }
            let covered = nearer.iter().any(|&nj| {
                let nt = &mesh.triangles[nj];
                point_in_triangle(&mesh.points, nt, p, tol)
                    && point_in_triangle(&mesh.points, nt, q, tol)
            });
            if !covered {
                out.push(((p.x, p.y), (q.x, q.y)));
            }
        }
        nearer.push(ti);
    }
    out
}

/// Point-in-triangle test in view-plane coordinates, inclusive within `tol`.
fn point_in_triangle(
    points: &[ProjectedPoint],
    tri: &ProjectedTriangle,
    p: &ProjectedPoint,
    tol: f64,
) -> bool {
    let a = (points[tri.a].x, points[tri.a].y);
    let b = (points[tri.b].x, points[tri.b].y);
    let c = (points[tri.c].x, points[tri.c].y);
    let d1 = cross2d(a, b, p.x, p.y);
    let d2 = cross2d(b, c, p.x, p.y);
    let d3 = cross2d(c, a, p.x, p.y);
    let has_neg = d1 < -tol || d2 < -tol || d3 < -tol;
    let has_pos = d1 > tol || d2 > tol || d3 > tol;
    !(has_neg && has_pos)
}

/// 2D cross product of `(b − a)` with `(p − a)`.
fn cross2d(a: (f64, f64), b: (f64, f64), px: f64, py: f64) -> f64 {
    (b.0 - a.0) * (py - a.1) - (b.1 - a.1) * (px - a.0)
}

/// Mesh `shape`, project it, and remove hidden edges.
///
/// The edge tolerance is scaled to the projected extent of the mesh.
pub fn wireframe_projection(
    shape: &crate::shape::TopoShape,
    view_dir: &GpVec,
    deflection: f64,
) -> Vec<((f64, f64), (f64, f64))> {
    let mesh = crate::shape_mesh::mesh_shape(shape, deflection);
    let pm = orthographic_project(&mesh, view_dir);
    let scale = pm
        .points
        .iter()
        .fold(1.0f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));
    visible_edges(&pm, 1e-9 * scale)
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpPnt;
    use occt_core::poly::triangulation::Triangle;

    fn box_mesh() -> ShapeMesh {
        crate::mesh::mesh_box((GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 1.0)))
    }

    #[test]
    fn box_along_z_backface_culls_front_face() {
        let pm = orthographic_project(&box_mesh(), &GpVec::new(0.0, 0.0, 1.0));
        // 12 input triangles; the two +Z (front) triangles are culled. The
        // four side faces are edge-on (normal_z == 0) and survive.
        assert_eq!(pm.triangles.len(), 10);
        // No surviving triangle faces the viewer.
        assert!(pm.triangles.iter().all(|t| t.normal_z <= 0.0));
        // The culled triangles are exactly the two +Z ones.
        assert_eq!(
            pm.triangles.iter().filter(|t| t.normal_z < 0.0).count(),
            2,
            "only the -Z face is strictly front-facing"
        );
    }

    #[test]
    fn visible_edges_of_single_triangle() {
        // Triangle in the z=0 plane wound to face -Z, so it survives a +Z view.
        let mesh = ShapeMesh {
            vertices: vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
            ],
            triangles: vec![Triangle::new(0, 1, 2)],
            source_shape: crate::abs::ShapeType::Face,
        };
        let pm = orthographic_project(&mesh, &GpVec::new(0.0, 0.0, 1.0));
        assert_eq!(pm.triangles.len(), 1);
        let edges = visible_edges(&pm, 1e-9);
        assert_eq!(edges.len(), 3);
    }

    #[test]
    fn painter_sort_orders_by_depth() {
        // Two -Z-facing triangles at z=0 and z=2 (both survive a +Z view).
        let mesh = ShapeMesh {
            vertices: vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(0.0, 0.0, 2.0),
                GpPnt::new(0.0, 1.0, 2.0),
                GpPnt::new(1.0, 0.0, 2.0),
            ],
            triangles: vec![
                Triangle::new(0, 1, 2), // depth 0 (near)
                Triangle::new(3, 4, 5), // depth 2 (far)
            ],
            source_shape: crate::abs::ShapeType::Face,
        };
        let pm = orthographic_project(&mesh, &GpVec::new(0.0, 0.0, 1.0));
        assert_eq!(pm.triangles.len(), 2);
        let order = painter_sort(&pm);
        // Far-to-near: triangle 1 (z=2) first, triangle 0 (z=0) last.
        assert!(pm.triangles[order[0]].depth >= pm.triangles[order[1]].depth);
        assert_eq!(order, vec![1, 0]);
    }

    #[test]
    fn wireframe_projection_of_box_nonempty() {
        let box_shape = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let edges = wireframe_projection(&box_shape, &GpVec::new(0.0, 0.0, 1.0), 0.25);
        assert!(!edges.is_empty());
    }
}
