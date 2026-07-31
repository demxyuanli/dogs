//! Module 1: mesh_to_brep — convert a `Poly_Triangulation` back into BRep shapes.
//!
//! Inverse of BRepMesh: take a triangulated surface (nodes + triangles) and
//! rebuild vertices, edges, wires, faces, a shell, and (when the mesh is a
//! closed manifold) a solid. Each triangle becomes a planar face whose
//! surface normal matches the triangle winding.
//! Source: `BRep_Builder` / `BRepBuilderAPI_MakeSolid` + `Poly_Triangulation`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt};
use occt_core::poly::Triangulation;
use occt_geom::GeomPlane;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::mesh::ShapeMesh;
use crate::shape::{Edge, Face, Shell, Solid, Vertex};

/// Result of converting a triangle mesh into BRep topology.
pub struct BrepFromMesh {
    pub solid: Option<Solid>,
    pub shell: Shell,
    pub faces: Vec<Face>,
    pub edges: Vec<Edge>,
    pub vertices: Vec<Vertex>,
}

/// Spatial-hash key for a point: coordinates rounded to 1e-9.
type PointKey = (i64, i64, i64);
fn point_key(p: &GpPnt) -> PointKey {
    (
        (p.x() * 1e9).round() as i64,
        (p.y() * 1e9).round() as i64,
        (p.z() * 1e9).round() as i64,
    )
}

/// Get (or create) the shared edge between two deduplicated vertex indices.
fn ensure_edge(
    builder: &TopoBuilder,
    edge_of_pair: &mut HashMap<(usize, usize), Edge>,
    edges: &mut Vec<Edge>,
    vertices: &[Vertex],
    a: usize,
    b: usize,
) -> Edge {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(e) = edge_of_pair.get(&key) {
        return e.clone();
    }
    let pa = BRepTool::vertex_point(&vertices[a]);
    let pb = BRepTool::vertex_point(&vertices[b]);
    let e = builder.make_edge_segment(&pa, &pb);
    edge_of_pair.insert(key, e.clone());
    edges.push(e.clone());
    e
}

/// Convert a `Poly_Triangulation` into a BRep shell (and a solid when closed).
pub fn triangulation_to_brep(tri: &Triangulation) -> BrepFromMesh {
    let builder = TopoBuilder::new();

    // 1. Vertices — one per distinct position, so shared triangle corners map
    //    to a single Vertex (dedupe by position via spatial hash).
    let mut vertices: Vec<Vertex> = Vec::new();
    let mut vertex_of_key: HashMap<PointKey, usize> = HashMap::new();
    let mut node_to_vidx: Vec<usize> = vec![usize::MAX; tri.nodes.len()];
    for (node, p) in tri.nodes.iter().enumerate() {
        let idx = *vertex_of_key.entry(point_key(p)).or_insert_with(|| {
            vertices.push(builder.make_vertex(*p, 0.0));
            vertices.len() - 1
        });
        node_to_vidx[node] = idx;
    }

    // 2. Unique edges (shared between adjacent triangles) and
    // 3. one planar face per triangle, wound to match the input triangle.
    let mut edge_of_pair: HashMap<(usize, usize), Edge> = HashMap::new();
    let mut edges: Vec<Edge> = Vec::new();
    let mut faces: Vec<Face> = Vec::new();
    let mut edge_use: HashMap<(usize, usize), usize> = HashMap::new();

    for t in &tri.triangles {
        let a = node_to_vidx[t.n0];
        let b = node_to_vidx[t.n1];
        let c = node_to_vidx[t.n2];
        let e_ab = ensure_edge(&builder, &mut edge_of_pair, &mut edges, &vertices, a, b);
        let e_bc = ensure_edge(&builder, &mut edge_of_pair, &mut edges, &vertices, b, c);
        let e_ca = ensure_edge(&builder, &mut edge_of_pair, &mut edges, &vertices, c, a);

        let pa = BRepTool::vertex_point(&vertices[a]);
        let pb = BRepTool::vertex_point(&vertices[b]);
        let pc = BRepTool::vertex_point(&vertices[c]);

        // Plane through the triangle; normal from winding (right-hand rule).
        let n = pb.coord.subtracted(&pa.coord).crossed(&pc.coord.subtracted(&pa.coord));
        let normal = GpDir::from_xyz(&n).unwrap_or(GpDir::new(0.0, 0.0, 1.0).unwrap());
        let pln = GpPln::new(GpAx3::from_ax1(&GpAx1::new(pa, normal)));

        let wire = builder.make_wire(&[e_ab, e_bc, e_ca]);
        faces.push(builder.make_face(Arc::new(GeomPlane::new(pln)), &[wire]));

        for &(x, y) in &[(a, b), (b, c), (c, a)] {
            let key = if x < y { (x, y) } else { (y, x) };
            *edge_use.entry(key).or_insert(0) += 1;
        }
    }

    // 4. Shell; solid iff every edge borders exactly two faces (closed
    //    2-manifold mesh).
    let shell = builder.make_shell(&faces);
    let closed = !edge_use.is_empty() && edge_use.values().all(|&n| n == 2);
    let solid = if closed { Some(builder.make_solid(&[shell.clone()])) } else { None };

    BrepFromMesh { solid, shell, faces, edges, vertices }
}

/// Convert a `ShapeMesh` into BRep topology (via a `Triangulation`).
pub fn shape_mesh_to_brep(mesh: &ShapeMesh) -> BrepFromMesh {
    let tri = crate::mesh::to_triangulation(mesh);
    triangulation_to_brep(&tri)
}

/// Number of vertices in the result.
pub fn brep_vertex_count(b: &BrepFromMesh) -> usize { b.vertices.len() }
/// Number of edges in the result.
pub fn brep_edge_count(b: &BrepFromMesh) -> usize { b.edges.len() }
/// Number of faces in the result.
pub fn brep_face_count(b: &BrepFromMesh) -> usize { b.faces.len() }

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::poly::triangulation::Triangle;

    #[test]
    fn unit_triangle_to_brep() {
        let p0 = GpPnt::new(0.0, 0.0, 0.0);
        let p1 = GpPnt::new(1.0, 0.0, 0.0);
        let p2 = GpPnt::new(0.0, 1.0, 0.0);
        let tri = Triangulation::new(vec![p0, p1, p2], vec![Triangle::new(0, 1, 2)]);
        let b = triangulation_to_brep(&tri);

        assert_eq!(brep_vertex_count(&b), 3);
        assert_eq!(brep_edge_count(&b), 3);
        assert_eq!(brep_face_count(&b), 1);
        assert!(b.solid.is_none(), "an open triangle is not a closed solid");

        // Vertex points match the input nodes.
        for p in [p0, p1, p2] {
            assert!(b.vertices.iter().any(|v| BRepTool::vertex_point(v).is_equal(&p)));
        }

        // The face's registered plane surface contains all three points.
        let surf = BRepTool::face_surface(&b.faces[0]).expect("face surface registered");
        let origin = surf.d0(0.0, 0.0);
        assert!(origin.is_equal(&p0), "plane passes through the first vertex");
        let xdir = surf.d0(1.0, 0.0).coord.subtracted(&origin.coord);
        let ydir = surf.d0(0.0, 1.0).coord.subtracted(&origin.coord);
        let nrm = xdir.crossed(&ydir);
        for p in [p0, p1, p2] {
            let off = p.coord.subtracted(&origin.coord).dot(&nrm);
            assert!(off.abs() < 1e-9, "point lies on the surface plane, offset {off}");
        }

        // Projecting the centroid onto the plane basis evaluates to the centroid.
        let centroid = GpPnt::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
        let v = centroid.coord.subtracted(&origin.coord);
        let at = surf.d0(v.dot(&xdir), v.dot(&ydir));
        assert!(at.distance(&centroid) < 1e-9);
    }
}
