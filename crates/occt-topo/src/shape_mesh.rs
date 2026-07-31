//! Arbitrary BRep meshing: faces and edges of any shape → a triangle soup.
//!
//! Collects every face via `TopExp_Explorer`, tessellates each with
//! `wireframe::face_to_triangles`, and concatenates them into one indexed
//! mesh. Also provides vertex welding for watertight output and analytic
//! surface-area / rough-volume helpers.

use std::collections::HashMap;

use occt_core::gp::GpPnt;
use occt_core::poly::triangulation::Triangle;

use crate::abs::ShapeType;
use crate::mesh::{compute_vertex_normals, mesh_surface_area, ShapeMesh};
use crate::shape::{Edge, TopoShape};
use crate::topexp::Explorer;
use crate::wireframe::{edge_to_polyline, face_to_triangles};

/// Mesh every face of `shape` into a single triangle soup.
///
/// Vertex lists of per-face tessellations are concatenated and triangle
/// indices shifted by the running offset.
pub fn mesh_shape(shape: &TopoShape, deflection: f64) -> ShapeMesh {
    let mut vertices: Vec<GpPnt> = Vec::new();
    let mut triangles: Vec<Triangle> = Vec::new();
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let face = crate::shape::Face(ex.current().clone());
        let (vs, ts) = face_to_triangles(&face, deflection);
        let offset = vertices.len();
        vertices.extend(vs);
        for t in ts {
            triangles.push(Triangle::new(offset + t.n0, offset + t.n1, offset + t.n2));
        }
        ex.next();
    }
    ShapeMesh { vertices, triangles, source_shape: shape.shape_type() }
}

/// Mesh `shape` and compute per-vertex (outward) normals.
pub fn mesh_shape_with_normals(shape: &TopoShape, deflection: f64) -> (ShapeMesh, Vec<GpPnt>) {
    let mesh = mesh_shape(shape, deflection);
    let normals = compute_vertex_normals(&mesh);
    (mesh, normals)
}

/// Discretize every edge of `shape` into a polyline.
pub fn mesh_wireframe(shape: &TopoShape, deflection: f64) -> Vec<Vec<GpPnt>> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(shape, ShapeType::Edge);
    while ex.more() {
        let e = Edge(ex.current().clone());
        out.push(edge_to_polyline(&e, deflection));
        ex.next();
    }
    out
}

/// Merge near-duplicate vertices within `tol` (spatial hash grid), rebuilding
/// triangle indices so the mesh becomes watertight for STL-style export.
pub fn weld_vertices(mesh: &mut ShapeMesh, tol: f64) {
    let tol = tol.max(1e-12);
    let cell = |p: &GpPnt| -> (i64, i64, i64) {
        (
            f64::floor(p.x() / tol) as i64,
            f64::floor(p.y() / tol) as i64,
            f64::floor(p.z() / tol) as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    let mut unique: Vec<GpPnt> = Vec::new();
    let mut remap = vec![0usize; mesh.vertices.len()];
    for (i, p) in mesh.vertices.iter().enumerate() {
        let c = cell(p);
        let mut found = None;
        'search: for dx in -1i64..=1 {
            for dy in -1i64..=1 {
                for dz in -1i64..=1 {
                    if let Some(bucket) = grid.get(&(c.0 + dx, c.1 + dy, c.2 + dz)) {
                        for &j in bucket {
                            if p.distance(&unique[j]) <= tol {
                                found = Some(j);
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
        match found {
            Some(j) => remap[i] = j,
            None => {
                let j = unique.len();
                unique.push(*p);
                grid.entry(c).or_default().push(j);
                remap[i] = j;
            }
        }
    }
    mesh.vertices = unique;
    for t in &mut mesh.triangles {
        t.n0 = remap[t.n0];
        t.n1 = remap[t.n1];
        t.n2 = remap[t.n2];
    }
}

/// Total surface area of `shape`'s mesh.
pub fn shape_surface_area(shape: &TopoShape, deflection: f64) -> f64 {
    mesh_surface_area(&mesh_shape(shape, deflection))
}

/// Rough volume estimate for a closed shell via the divergence theorem:
/// `V = (1/6) |Σ a·(b×c)|` over the triangle mesh. Open or inconsistently
/// oriented meshes give unreliable results — treat this as an estimate, not a
/// precise solid volume.
pub fn shape_volume(shape: &TopoShape, deflection: f64) -> f64 {
    let mesh = mesh_shape(shape, deflection);
    let mut vol = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t.n0].coord;
        let b = mesh.vertices[t.n1].coord;
        let c = mesh.vertices[t.n2].coord;
        vol += a.dot_cross(&b, &c);
    }
    (vol / 6.0).abs()
}

/// Number of faces in `shape` (including itself if it is a face).
pub fn count_faces(shape: &TopoShape) -> usize {
    let mut n = 0;
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        n += 1;
        ex.next();
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::shape::Compound;
    use crate::wireframe::tests::square_face;

    fn square_face_compound() -> TopoShape {
        let f = square_face();
        let b = TopoBuilder::new();
        let mut c = Compound::new();
        b.add_compound(&mut c, &f.0);
        c.0
    }

    #[test]
    fn mesh_shape_produces_correct_area() {
        let shape = square_face_compound();
        let mesh = mesh_shape(&shape, 0.25);
        assert_eq!(mesh.source_shape, ShapeType::Compound);
        assert!(mesh.triangles.len() >= 6, "tris {}", mesh.triangles.len());
        assert!((mesh_surface_area(&mesh) - 1.0).abs() < 1e-6, "area {}", mesh_surface_area(&mesh));
    }

    #[test]
    fn mesh_shape_with_normals_is_unit() {
        let shape = square_face_compound();
        let (mesh, normals) = mesh_shape_with_normals(&shape, 0.25);
        assert_eq!(mesh.vertices.len(), normals.len());
        for n in &normals {
            assert!((n.coord.modulus() - 1.0).abs() < 1e-6, "normal modulus {}", n.coord.modulus());
        }
    }

    #[test]
    fn mesh_wireframe_returns_per_edge_polylines() {
        let shape = square_face_compound();
        let lines = mesh_wireframe(&shape, 0.01);
        assert_eq!(lines.len(), 4, "expected 4 edges, got {}", lines.len());
        for l in &lines {
            assert!(l.len() >= 2);
        }
    }

    #[test]
    fn weld_vertices_merges_duplicates() {
        let mut mesh = crate::mesh::mesh_box((GpPnt::zero(), GpPnt::new(1.0, 1.0, 1.0)));
        let n = mesh.vertices.len();
        let dup = mesh.vertices[0];
        mesh.vertices.push(dup);
        mesh.triangles.push(Triangle::new(0, 0, n));
        weld_vertices(&mut mesh, 1e-9);
        assert_eq!(mesh.vertices.len(), n, "duplicate not merged");
    }

    #[test]
    fn count_faces_in_compound() {
        let shape = square_face_compound();
        assert_eq!(count_faces(&shape), 1);
    }

    #[test]
    fn shape_surface_area_of_square() {
        let shape = square_face_compound();
        assert!((shape_surface_area(&shape, 0.25) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn shape_volume_flat_face_is_zero() {
        // An open planar mesh encloses no volume.
        let shape = square_face_compound();
        assert!(shape_volume(&shape, 0.25) < 1e-12, "vol {}", shape_volume(&shape, 0.25));
    }

    #[test]
    fn volume_formula_holds_for_box_mesh() {
        let m = crate::mesh::mesh_box((GpPnt::zero(), GpPnt::new(1.0, 1.0, 1.0)));
        let mut vol = 0.0;
        for t in &m.triangles {
            let a = m.vertices[t.n0].coord;
            let b = m.vertices[t.n1].coord;
            let c = m.vertices[t.n2].coord;
            vol += a.dot_cross(&b, &c);
        }
        assert!((vol / 6.0 - 1.0).abs() < 1e-9, "vol {}", vol / 6.0);
    }
}
