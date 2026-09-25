//! Module 1: mesh_to_brep — convert a `Poly_Triangulation` back into BRep shapes.
//!
//! Inverse of BRepMesh: take a triangulated surface (nodes + triangles) and
//! rebuild vertices, edges, wires, faces, a shell, and (when the mesh is a
//! closed manifold) a solid. Each triangle becomes a planar face whose
//! surface normal matches the triangle winding.
//!
//! Source: `BRep_Builder` / `BRepBuilderAPI_MakeSolid` + `Poly_Triangulation`,
//! plus `BRepLib_MakeWire` for the wires (`BRepBuilderAPI_MakeWire`): the
//! triangle's three edges are **shared** with the neighbouring triangles (one
//! `TShape` per mesh edge) and each face adds them with its own traversal
//! orientation, then they are chained exactly as `BRepLib_MakeWire::Add`
//! (`BRepLib_MakeWire.cxx:123-453`) does — see [`crate::brep_lib_make_wire`]
//! (board task T-32).

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt};
use occt_core::poly::Triangulation;
use occt_geom::GeomPlane;

use crate::abs::Orientation;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::mesh::ShapeMesh;
use crate::shape::{Edge, Face, Shell, Solid, Vertex, Wire};

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

/// Get (or create) the edge between two deduplicated vertex indices, **oriented
/// along `a -> b`**.
///
/// The edge `TShape` is shared between the two triangles that use it (the
/// canonical copy always runs `lo -> hi`), while the orientation is per use:
/// `BRep_Builder`/`BRepBuilderAPI_MakeFace` add the *same* `TopoDS_Edge` to each
/// face with that face's orientation (`e.Oriented(FORWARD|REVERSED)`), which is
/// what keeps every face's wire wound with its own triangle. Handing the cached
/// edge to the wire without re-orienting it would instead let whichever triangle
/// created the edge first dictate the wire direction, and a wire wound against
/// its plane normal reads as clockwise in UV — breaking the `FORWARD` face
/// invariant `IntTools_FClass2d` relies on (`IsHole`).
fn ensure_edge(
    builder: &TopoBuilder,
    edge_of_pair: &mut HashMap<(usize, usize), Edge>,
    edges: &mut Vec<Edge>,
    vertices: &[Vertex],
    a: usize,
    b: usize,
) -> Edge {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    let e = if let Some(e) = edge_of_pair.get(&(lo, hi)) {
        e.clone()
    } else {
        let pa = BRepTool::vertex_point(&vertices[lo]);
        let pb = BRepTool::vertex_point(&vertices[hi]);
        let e = builder.make_edge_segment_with_vertices(&pa, &pb, &vertices[lo], &vertices[hi]);
        edge_of_pair.insert((lo, hi), e.clone());
        edges.push(e.clone());
        e
    };
    let ori = if a < b {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    Edge(e.0.oriented(ori))
}

/// Build a face's wire through `BRepLib_MakeWire` (`BRepBuilderAPI_MakeWire`),
/// which **decides the orientation of every edge** so that the wire chains
/// head-to-tail (`brep_lib_make_wire`, `BRepLib_MakeWire.cxx:123-453`).
///
/// OCCT's caller raises `StdFail_NotDone` when the builder is not done. This
/// function has no failure channel (it feeds the port-local mesh→BRep
/// conversion), so a non-done result keeps the `BRep_Builder` level append as a
/// documented fallback — unreachable for mesh triangles, whose three edges share
/// deduplicated vertices by identity.
fn build_wire(builder: &TopoBuilder, edges: &[Edge]) -> Wire {
    let mut mw = MakeWire::new();
    for e in edges {
        if mw.add(e).is_err() {
            return builder.make_wire(edges);
        }
    }
    mw.wire()
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

        let wire = build_wire(&builder, &[e_ab, e_bc, e_ca]);
        faces.push(builder.make_face(Arc::new(GeomPlane::new(pln)), &[wire]));

        for &(x, y) in &[(a, b), (b, c), (c, a)] {
            let key = if x < y { (x, y) } else { (y, x) };
            *edge_use.entry(key).or_insert(0) += 1;
        }
    }

    // 4. Shell; solid iff every edge borders exactly two faces (closed
    //    2-manifold mesh). A closed solid must have its material **inside**:
    //    `BRepTools::OrientClosedSolid` (`brep_class3d::orient_closed_solid`)
    //    reverses it when the infinite point classifies as IN, which is the
    //    case when the input triangles are wound the other way round.
    //
    //    The winding of the input mesh is read once, deterministically, from the
    //    divergence-theorem volume (`Σ a·(b×c)`: negative for a mesh whose
    //    normals point into the material). A face is FORWARD when its plane
    //    normal (the triangle normal) already points away from the material, and
    //    REVERSED otherwise; the wires stay wound with their triangle, so the
    //    outer wire is always CCW around the plane normal, which is the
    //    invariant `IntTools_FClass2d` (and every UV area) relies on. Relying on
    //    the infinite-point classifier alone would leave the material side of
    //    the faces dependent on the (unordered) result of that classification.
    let mut signed6 = 0.0;
    for t in &tri.triangles {
        let a = tri.nodes[t.n0].coord;
        let b = tri.nodes[t.n1].coord;
        let c = tri.nodes[t.n2].coord;
        signed6 += a.dot_cross(&b, &c);
    }
    if signed6 < 0.0 {
        for f in faces.iter_mut() {
            f.0.set_orientation(Orientation::Reversed);
        }
    }
    let shell = builder.make_shell(&faces);
    let closed = !edge_use.is_empty() && edge_use.values().all(|&n| n == 2);
    let _solid = if closed {
        let mut s = builder.make_solid(&[shell.clone()]).0;
        // Verification step (`BRepTools::OrientClosedSolid`): with the winding
        // above it must report the material already inside and leave the solid
        // alone.
        crate::brep_class3d::orient_closed_solid(&mut s);
        Solid::wrap(s)
    } else {
        None
    };
    let solid = if closed {
        let mut s = builder.make_solid(&[shell.clone()]).0;
        crate::brep_class3d::orient_closed_solid(&mut s);
        Solid::wrap(s)
    } else {
        None
    };

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
