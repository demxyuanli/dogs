//! Swept shapes — extrusion (prism) of a planar face or polygon.
//! Source: `BRepPrimAPI_MakePrism.hxx`.
//!
//! Extrudes a planar boundary along a direction vector by a height, building
//! a closed solid: base face (the input), top face (translated), and one
//! planar lateral face per boundary edge. The sweep direction must have a
//! non-zero component along the face normal (a pure in-plane direction would
//! produce a degenerate solid).

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_geom::{Curve, GeomLine, GeomPlane, Surface};

use crate::brep_surface;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Solid, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full;

/// Result of a prism sweep: the solid plus its distinguishing faces.
#[derive(Debug, Clone)]
pub struct Prism {
    pub solid: Solid,
    pub base_face: Face,
    pub top_face: Face,
    pub lateral_faces: Vec<Face>,
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
}

/// Shift a point by a vector.
fn shift(p: &GpPnt, d: &GpVec) -> GpPnt {
    GpPnt::new(p.x() + d.x(), p.y() + d.y(), p.z() + d.z())
}

/// Plane through three non-collinear points.
fn plane_through3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> GpPln {
    let ab = GpVec::from_pnts(a, b);
    let ac = GpVec::from_pnts(a, c);
    let n = ab.xyz().crossed(ac.xyz());
    let d = GpDir::from_vec(&GpVec::new(n.x, n.y, n.z)).expect("sweep: degenerate lateral face");
    let z_axis = GpDir::new(0.0, 0.0, 1.0).unwrap();
    let x_dir = if d.is_normal(&z_axis) { z_axis } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
    GpPln::new(GpAx3::new(*a, d, &x_dir).expect("sweep: lateral plane frame"))
}

/// A line edge through two existing vertices, registering geometry and
/// attaching the shared vertices as children.
fn edge_through(b: &TopoBuilder, v1: &Vertex, v2: &Vertex) -> Edge {
    let p1 = GeometryRegistry::global().vertex_point(&v1.0);
    let p2 = GeometryRegistry::global().vertex_point(&v2.0);
    let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2))
        .unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
    let mut e = b.make_edge(Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir))), 0.0, p1.distance(&p2));
    b.add(&mut e.0, &v1.0);
    b.add(&mut e.0, &v2.0);
    e
}

/// Boundary vertices of a planar face's outer wire, in traversal order.
///
/// The reconstruction is direction-robust: it walks the edges finding the next
/// vertex by matching endpoints (edges may be stored in either orientation
/// relative to the wire), and stops when it returns to the first vertex.
fn base_ring_vertices(face: &Face) -> Vec<Vertex> {
    // ponytail: only the outer wire is swept; inner (hole) wires are ignored.
    let Some(w) = topo_tools_full::wires_of_face(face).into_iter().next() else {
        return Vec::new();
    };
    let edges = topo_tools_full::edges_of_wire(&w);
    if edges.is_empty() {
        return Vec::new();
    }
    let pos = |v: &Vertex| GeometryRegistry::global().vertex_point(&v.0);
    let start = pos(&edge_first(&edges[0]));
    let mut ring: Vec<Vertex> = Vec::new();
    let mut cur = edge_first(&edges[0]);
    ring.push(cur.clone());
    let mut used = vec![false; edges.len()];
    used[0] = true;
    loop {
        // Find the next edge sharing `cur` as one endpoint.
        let mut found = false;
        for (i, e) in edges.iter().enumerate() {
            if used[i] {
                continue;
            }
            let (a, b) = topo_tools_full::edge_vertices(e);
            let (Some(a), Some(b)) = (a, b) else { continue };
            let next;
            if pos(&a).distance(&pos(&cur)) < 1e-12 {
                next = b;
            } else if pos(&b).distance(&pos(&cur)) < 1e-12 {
                next = a;
            } else {
                continue;
            }
            if pos(&next).distance(&start) < 1e-12 {
                // Ring closes; do not push the repeated start.
                return ring;
            }
            ring.push(next.clone());
            cur = next;
            used[i] = true;
            found = true;
            break;
        }
        if !found {
            return ring;
        }
    }
}

fn edge_first(e: &Edge) -> Vertex {
    let (a, _) = topo_tools_full::edge_vertices(e);
    a.expect("sweep: edge has no vertices")
}

/// Extrude a planar face by the full displacement vector `d`.
pub fn prism_from_face(face: &Face, d: &GpVec) -> Prism {
    assert!(d.xyz().square_modulus() > 0.0, "sweep: zero sweep displacement");
    let b = TopoBuilder::new();

    let base_verts = base_ring_vertices(face);
    assert!(base_verts.len() >= 3, "sweep: base face needs >= 3 boundary vertices");
    let n = base_verts.len();
    let base_pts: Vec<GpPnt> = base_verts.iter().map(|v| GeometryRegistry::global().vertex_point(&v.0)).collect();

    // The displacement vector is passed through unchanged (it already encodes
    // direction × magnitude).
    let d = d.clone();

    // Top ring: one new vertex per base corner, displaced by `d`.
    let top_verts: Vec<Vertex> = base_pts.iter().map(|p| b.make_vertex(shift(p, &d), 0.0)).collect();

    // Top edges and vertical edges.
    let mut top_edges = Vec::with_capacity(n);
    let mut vert_edges = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        top_edges.push(edge_through(&b, &top_verts[i], &top_verts[j]));
        vert_edges.push(edge_through(&b, &base_verts[i], &top_verts[i]));
    }

    // Base edges: the input face's boundary edges (already registered).
    let base_edges: Vec<Edge> = topo_tools_full::wires_of_face(face)
        .first()
        .map(|w| topo_tools_full::edges_of_wire(w))
        .unwrap_or_default();
    assert_eq!(base_edges.len(), n, "sweep: base wire edge count mismatch");

    // Top face: base plane translated by `d`.
    let base_plane = brep_surface::face_plane(face).expect("sweep: base face must be planar");
    let mut top_plane = base_plane.clone();
    top_plane.set_location(&shift(&base_plane.location(), &d));
    let top_wire = b.make_wire(&top_edges);
    let top_face = b.make_face(Arc::new(GeomPlane::new(top_plane)), &[top_wire]);

    // Lateral faces: one planar quad per base edge.
    let mut lateral_faces = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let wire = b.make_wire(&[base_edges[i].clone(), vert_edges[j].clone(), top_edges[i].clone(), vert_edges[i].clone()]);
        let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_through3(
            &base_pts[i], &base_pts[j], &shift(&base_pts[j], &d),
        )));
        lateral_faces.push(b.make_face(surface, &[wire]));
    }

    let mut all_faces = vec![face.clone(), top_face.clone()];
    all_faces.extend(lateral_faces.iter().cloned());
    let shell = b.make_shell(&all_faces);
    let solid = b.make_solid(&[shell]);

    let mut vertices = base_verts;
    vertices.extend(top_verts);
    let mut edges = base_edges;
    edges.extend(top_edges);
    edges.extend(vert_edges);

    Prism { solid, base_face: face.clone(), top_face, lateral_faces, vertices, edges }
}

/// Build a prism by first constructing a planar base face from a polygon of
/// points (coplanar, listed in boundary order), then extruding it.
pub fn prism_from_polygon(points: &[GpPnt], direction: &GpVec, height: f64) -> Prism {
    assert!(points.len() >= 3, "sweep: polygon needs >= 3 points");
    let b = TopoBuilder::new();

    // Base ring vertices.
    let mut verts: Vec<Vertex> = points.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
    // Top ring.
    let mag = direction.xyz().modulus();
    assert!(mag > 1e-30, "sweep: zero sweep direction");
    let d = GpVec::new(
        direction.x() / mag * height,
        direction.y() / mag * height,
        direction.z() / mag * height,
    );
    let top_pts: Vec<GpPnt> = points.iter().map(|p| shift(p, &d)).collect();

    // Base wire.
    let mut base_edges = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let mut e = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(points[i], {
                let dir = GpDir::from_vec(&GpVec::from_pnts(&points[i], &points[j])).unwrap();
                dir
            }))),
            0.0, points[i].distance(&points[j]),
        );
        b.add(&mut e.0, &verts[i].0);
        b.add(&mut e.0, &verts[j].0);
        base_edges.push(e);
    }
    let base_wire = b.make_wire(&base_edges);
    let base_plane = plane_through3(&points[0], &points[1], &points[2]);
    let base_face = b.make_face(Arc::new(GeomPlane::new(base_plane)), &[base_wire]);

    // Build the prism using the shared path; `d` is the full displacement.
    prism_from_face(&base_face, &d)
}

/// Volume of a straight prism: base polygon area × sweep height.
pub fn prism_volume(prism: &Prism) -> f64 {
    let n = prism.lateral_faces.len();
    if n < 3 {
        return 0.0;
    }
    let base_pts: Vec<GpPnt> = prism.vertices[..n]
        .iter()
        .map(|v| GeometryRegistry::global().vertex_point(&v.0))
        .collect();
    let top_pts: Vec<GpPnt> = prism.vertices[n..]
        .iter()
        .map(|v| GeometryRegistry::global().vertex_point(&v.0))
        .collect();
    // Base area of the planar polygon (3D cross-product shoelace) / 2.
    let mut acc = occt_core::gp::GpXyz::zero();
    for i in 0..n {
        let j = (i + 1) % n;
        acc = acc.added(&base_pts[i].coord.crossed(&base_pts[j].coord));
    }
    let area = 0.5 * acc.modulus();
    // Height = mean displacement magnitude of corresponding base/top corners.
    let mut h = 0.0;
    for i in 0..n {
        h += base_pts[i].distance(&top_pts[i]);
    }
    area * h / n as f64
}

/// Boundary counts of the prism solid: (vertices, edges, faces).
pub fn prism_counts(prism: &Prism) -> (usize, usize, usize) {
    (
        prism.vertices.len(),
        prism.edges.len(),
        2 + prism.lateral_faces.len(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::primitives::BRepPrimBox;

    fn approx(a: f64, b: f64) -> bool { (a - b).abs() < 1e-9 * b.abs().max(1.0) }

    #[test]
    fn unit_square_prism_is_cuboid() {
        // Base square in z=0 plane, extrude +Z by 2 → 2×2×2 cuboid.
        let sq = [GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(2.,2.,0.), GpPnt::new(0.,2.,0.)];
        let prism = prism_from_polygon(&sq, &GpVec::new(0.0, 0.0, 1.0), 2.0);
        let (nv, ne, nf) = prism_counts(&prism);
        assert_eq!((nv, ne, nf), (8, 12, 6));
        assert!(approx(prism_volume(&prism), 8.0));
        // Top face sits at z=2.
        let top_plane = brep_surface::face_plane(&prism.top_face).unwrap();
        assert!(approx(top_plane.location().z(), 2.0));
    }

    #[test]
    fn prism_from_box_face() {
        let box_ = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let face = crate::topo_tools_full::faces_of(&box_.solid.0)
            .into_iter()
            .find(|f| BRepTool::face_surface(f).map(|s| (s.d0(0.0, 0.0).z() - 0.0).abs() < 1e-9).unwrap_or(false))
            .expect("bottom face");
        // Sweep the unit bottom face upward by 3 → 1×1×3 box.
        let prism = prism_from_face(&face, &GpVec::new(0.0, 0.0, 3.0));
        assert_eq!(prism_counts(&prism), (8, 12, 6));
        assert!(approx(prism_volume(&prism), 3.0));
    }

    #[test]
    fn triangle_prism() {
        let tri = [GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.)];
        let prism = prism_from_polygon(&tri, &GpVec::new(0.0, 0.0, 5.0), 5.0);
        // 6 vertices, 9 edges, 5 faces.
        assert_eq!(prism_counts(&prism), (6, 9, 5));
        // Volume = base area (1/2) × height (5) = 2.5.
        assert!(approx(prism_volume(&prism), 2.5));
    }

    #[test]
    #[should_panic]
    fn zero_height_rejected() {
        let sq = [GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.)];
        prism_from_polygon(&sq, &GpVec::new(0.0, 0.0, 1.0), 0.0);
    }
}
