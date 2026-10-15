//! Module 1: mesh_to_brep — convert a `Poly_Triangulation` back into BRep shapes.
//!
//! Inverse of BRepMesh: take a triangulated surface (nodes + triangles) and
//! rebuild vertices, edges, wires and faces, plus (as a port-local extension
//! described below) a shell and a solid.
//!
//! Source: `BRepBuilderAPI_MakeShapeOnMesh` (`BRepBuilderAPI_MakeShapeOnMesh.cxx:74-263`).
//! Its control flow is the one implemented by [`MakeShapeOnMesh`]:
//!
//! * mesh nodes map 1:1 to vertices (`cxx:96-111`, `aPnt2VertexMap.Add(i, aV)`
//!   per node index, no welding); the vertex is a `BRepBuilderAPI_MakeVertex`,
//!   i.e. `BRepLib_MakeVertex` -> `BRep_Builder::MakeVertex(V, P,
//!   BRepLib::Precision())` (`BRepLib_MakeVertex.cxx:26-31`), tolerance
//!   `Precision::Confusion()` (`BRepLib.cxx:84`), not 0;
//! * mesh edges become shared topological edges keyed on the node-index pair
//!   `(min, max)` (`cxx:31-55`, `cxx:150-179`); `NCollection_IndexedDataMap::Add`
//!   keeps the value first bound to a key, and in the cxx an instance's curve
//!   direction and its orientation flag are complementary, so the stored edge's
//!   *effective* direction is always `lo -> hi` -- this port stores that
//!   canonical geometry and puts the per-face orientation on the instance;
//! * each triangle becomes a planar face (`cxx:181-261`): the three edges are
//!   re-fetched from that map and re-oriented along the triangle's own traversal
//!   (`cxx:213-228`), chained by `BRepBuilderAPI_MakeWire`
//!   (`BRepLib_MakeWire::Add`, `BRepLib_MakeWire.cxx:123-453` — see
//!   [`crate::brep_lib_make_wire`], board task T-32), the plane normal is the
//!   cross product of the first two edges' line directions (`cxx:238-243`) and
//!   the plane goes through `Node(anIdx[0])` (`cxx:255`);
//! * the result is a `TopoDS_Compound` holding those faces (`cxx:182-184`,
//!   `cxx:257-263`), and a null or empty mesh leaves the builder not-done
//!   (`cxx:80-85`).
//!
//! Port-local extensions, **not** in the OCCT class and marked `LOCAL:` at each
//! site: the node welding of [`weld_nodes`] (so that mesh sources repeating a
//! position per incident triangle still yield a closed shell — OCCT keeps the
//! nodes 1:1 and answers a compound of disconnected faces for such input) and
//! the shell/solid part of [`BrepFromMesh`] (closed-manifold test,
//! `OrientClosedSolid`, divergence-theorem winding pass).

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_core::poly::triangulation::Triangle;
use occt_core::poly::Triangulation;
use occt_core::precision::{CONFUSION, REAL_SMALL, SQUARE_CONFUSION};
use occt_geom::GeomPlane;

use crate::abs::Orientation;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::mesh::ShapeMesh;
use crate::shape::{Compound, Edge, Face, Shell, Solid, Vertex, Wire};

/// Result of converting a triangle mesh into BRep topology.
///
/// `faces` / `edges` / `vertices` are the faithful `MakeShapeOnMesh`
/// counterparts; `shell` and `solid` are the port-local assembly on top.
pub struct BrepFromMesh {
    pub solid: Option<Solid>,
    pub shell: Shell,
    pub faces: Vec<Face>,
    pub edges: Vec<Edge>,
    pub vertices: Vec<Vertex>,
}

/// `BRepBuilderAPI_MakeShapeOnMesh` (`BRepBuilderAPI_MakeShapeOnMesh.cxx:26-263`).
pub struct MakeShapeOnMesh {
    done: bool,
    compound: Compound,
    faces: Vec<Face>,
    edges: Vec<Edge>,
    vertices: Vec<Vertex>,
}

impl MakeShapeOnMesh {
    /// `BRepBuilderAPI_MakeShapeOnMesh::Build` (`cxx:74-263`).
    pub fn new(mesh: &Triangulation) -> Self {
        let mut m = Self {
            done: false,
            compound: Compound::new(),
            faces: Vec::new(),
            edges: Vec::new(),
            vertices: Vec::new(),
        };
        m.build(mesh);
        m
    }

    fn build(&mut self, mesh: &Triangulation) {
        let builder = TopoBuilder::new();

        // `cxx:80-85`: null / empty mesh leaves the builder not-done.
        if mesh.nodes.is_empty() || mesh.triangles.is_empty() {
            return;
        }

        // `cxx:96-111`: one `BRepBuilderAPI_MakeVertex` per node, node index
        // 1:1 (`aPnt2VertexMap.Add(i, aV)`), tolerance `BRepLib::Precision()`
        // = `Precision::Confusion()` (`BRepLib.cxx:84`), not 0.
        for p in mesh.nodes.iter() {
            self.vertices.push(builder.make_vertex(*p, CONFUSION));
        }

        // `cxx:113-179`: shared edges. The first triangle that uses a pair
        // decides its curve direction; the stored edge is canonical (lo -> hi).
        let mut edge_of_pair: HashMap<(usize, usize), Edge> = HashMap::new();
        for t in &mesh.triangles {
            let idx = [t.n0, t.n1, t.n2];
            // `cxx:126-130`: index-degenerated triangles are skipped.
            if idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2] {
                continue;
            }
            // `cxx:133-141`: geometrically degenerated triangles are skipped
            // (`aD1/aD2/aD3 < gp::Resolution()` — `gp::Resolution()` is
            // `RealSmall()`, i.e. `REAL_SMALL`, so only exactly coincident nodes
            // are caught; see `occt_core::precision::RESOLUTION`).
            let p = [mesh.nodes[idx[0]], mesh.nodes[idx[1]], mesh.nodes[idx[2]]];
            if p[0].square_distance(&p[1]) < REAL_SMALL
                || p[0].square_distance(&p[2]) < REAL_SMALL
                || p[1].square_distance(&p[2]) < REAL_SMALL
            {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (idx[k], idx[(k + 1) % 3]);
                let key = if a < b { (a, b) } else { (b, a) };
                // `NCollection_IndexedDataMap::Add` keeps the first value.
                if edge_of_pair.contains_key(&key) {
                    continue;
                }
                let (lo, hi) = (p[k], p[(k + 1) % 3]);
                let (lo, hi) = if a < b { (lo, hi) } else { (hi, lo) };
                let (vlo, vhi) = if a < b { (a, b) } else { (b, a) };
                let e = builder.make_edge_segment_with_vertices(
                    &lo,
                    &hi,
                    &self.vertices[vlo],
                    &self.vertices[vhi],
                );
                edge_of_pair.insert(key, e.clone());
                self.edges.push(e);
            }
        }

        // `cxx:181-261`: one planar face per triangle, added to a compound.
        for t in &mesh.triangles {
            let idx = [t.n0, t.n1, t.n2];
            let keys = [
                if idx[0] < idx[1] { (idx[0], idx[1]) } else { (idx[1], idx[0]) },
                if idx[1] < idx[2] { (idx[1], idx[2]) } else { (idx[2], idx[1]) },
                if idx[2] < idx[0] { (idx[2], idx[0]) } else { (idx[0], idx[2]) },
            ];
            // `cxx:204-211`: a triangle whose edges were skipped by the
            // degeneracy tests is skipped here too.
            let Some(stored) = keys.iter().map(|k| edge_of_pair.get(k)).collect::<Option<Vec<_>>>()
            else {
                continue;
            };
            // `cxx:213-228`: re-orient each edge along the triangle's traversal.
            let ori = [
                if idx[0] < idx[1] { Orientation::Forward } else { Orientation::Reversed },
                if idx[1] < idx[2] { Orientation::Forward } else { Orientation::Reversed },
                if idx[2] < idx[0] { Orientation::Forward } else { Orientation::Reversed },
            ];
            let tri_edges: Vec<Edge> =
                stored.iter().zip(ori).map(|(e, o)| Edge(e.0.oriented(o))).collect();

            // `cxx:230-234`: `BRepBuilderAPI_MakeWire` chains the three edges.
            let wire = build_wire(&builder, &tri_edges);

            // `cxx:238-247`: the plane normal is `d1 x d2` of the first two
            // edges' line directions. `BRepAdaptor_Curve::Line()` answers the
            // edge curve's `gp_Lin` (`GeomAdaptor_Curve.cxx:1158-1163`) and
            // `gp_Dir` is normalised, so the cross product is taken on **unit**
            // directions and the `Precision::SquareConfusion()` test rejects a
            // near-parallel (sliver) triangle, not a small one.
            //
            // `cxx:248-253` reverses that normal when one of the two edges is
            // `TopAbs_REVERSED`, turning the raw curve direction into the
            // edge's *effective* direction, i.e. the triangle's traversal.
            // Every shared edge here is stored canonically (`lo -> hi`) with
            // the per-face orientation set by `cxx:213-228`, so the two
            // directions below are taken straight from the traversal
            // (`Node(idx[k]) -> Node(idx[k+1])`) and both `aN.Reverse()` tests
            // are already folded in.
            let p0 = mesh.nodes[idx[0]];
            let (Ok(d1), Ok(d2)) = (
                GpDir::from_vec(&GpVec::from_pnts(&p0, &mesh.nodes[idx[1]])),
                GpDir::from_vec(&GpVec::from_pnts(&mesh.nodes[idx[1]], &mesh.nodes[idx[2]])),
            ) else {
                continue;
            };
            let n = d1.xyz().crossed(d2.xyz());
            if n.square_modulus() < SQUARE_CONFUSION {
                continue;
            }
            // `cxx:255`: the plane goes through `Node(anIdx[0])`.
            let normal = GpDir::from_xyz(&n)
                .unwrap_or_else(|_| GpDir::new(0.0, 0.0, 1.0).expect("canonical +Z"));
            let pln = GpPln::new(GpAx3::from_ax1(&GpAx1::new(p0, normal)));
            let face = builder.make_face(Arc::new(GeomPlane::new(pln)), &[wire]);
            builder.add_compound(&mut self.compound, &face.0);
            self.faces.push(face);
        }

        self.done = true;
    }

    /// `BRepBuilderAPI_MakeShape::IsDone()`: true once `Build` reached
    /// `Done()` (`cxx:263`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `myShape` — the `TopoDS_Compound` of `cxx:182-184`, `cxx:257-261`.
    pub fn shape(&self) -> &Compound {
        &self.compound
    }

    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
}

/// LOCAL: weld nodes that repeat a position, so a mesh whose triangles each
/// carry their own copy of a corner (STL/OBJ imports, `Poly_Triangulation` from
/// a faceted writer) still yields one shared `TShape` per geometric edge.
/// OCCT has no such step in `BRepBuilderAPI_MakeShapeOnMesh`: it keeps the nodes
/// 1:1 (`cxx:96-111`) and returns a compound of disconnected faces for that
/// input. The grid is the port's own (`1e-9`), not an OCCT tolerance.
pub(crate) fn weld_nodes(tri: &Triangulation) -> (Vec<GpPnt>, Vec<Triangle>) {
    type PointKey = (i64, i64, i64);
    let point_key = |p: &GpPnt| -> PointKey {
        (
            (p.x() * 1e9).round() as i64,
            (p.y() * 1e9).round() as i64,
            (p.z() * 1e9).round() as i64,
        )
    };
    let mut nodes: Vec<GpPnt> = Vec::new();
    let mut vertex_of_key: HashMap<PointKey, usize> = HashMap::new();
    let mut node_to_welded: Vec<usize> = vec![usize::MAX; tri.nodes.len()];
    for (node, p) in tri.nodes.iter().enumerate() {
        let idx = *vertex_of_key.entry(point_key(p)).or_insert_with(|| {
            nodes.push(*p);
            nodes.len() - 1
        });
        node_to_welded[node] = idx;
    }
    let triangles = tri
        .triangles
        .iter()
        .map(|t| Triangle::new(node_to_welded[t.n0], node_to_welded[t.n1], node_to_welded[t.n2]))
        .collect();
    (nodes, triangles)
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
///
/// The conversion itself is [`MakeShapeOnMesh`] (`BRepBuilderAPI_MakeShapeOnMesh.cxx:74-263`);
/// the welding and the shell/solid assembly are the port-local extension
/// described in the module header.
pub fn triangulation_to_brep(tri: &Triangulation) -> BrepFromMesh {
    let builder = TopoBuilder::new();

    // LOCAL: weld repeated node positions before the faithful conversion, so a
    // mesh that carries one node copy per incident triangle still closes.
    let (nodes, triangles) = weld_nodes(tri);
    let welded = Triangulation::new(nodes, triangles);
    let conv = MakeShapeOnMesh::new(&welded);

    let vertices = conv.vertices().to_vec();
    let edges = conv.edges().to_vec();
    let mut faces = conv.faces().to_vec();

    // LOCAL: edge use counts, for the closed 2-manifold test below.
    let mut edge_use: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &welded.triangles {
        // A pole triangle has UV-distinct corners that weld to one 3-D point,
        // so two of its indices coincide; `BRepBuilderAPI_MakeEdge` returns
        // `NotDone` for a zero-length edge, so such a triangle has no face
        // (`BRepBuilderAPI_MakeShapeOnMesh.cxx:126-141`, `:204-211`).
        if t.n0 == t.n1 || t.n1 == t.n2 || t.n2 == t.n0 {
            continue;
        }
        for &(x, y) in &[(t.n0, t.n1), (t.n1, t.n2), (t.n2, t.n0)] {
            let key = if x < y { (x, y) } else { (y, x) };
            *edge_use.entry(key).or_insert(0) += 1;
        }
    }

    // LOCAL: shell; solid iff every edge borders exactly two faces (closed
    // 2-manifold mesh). A closed solid must have its material **inside**:
    // `BRepTools::OrientClosedSolid` (`brep_class3d::orient_closed_solid`)
    // reverses it when the infinite point classifies as IN, which is the
    // case when the input triangles are wound the other way round.
    //
    // The winding of the input mesh is read once, deterministically, from the
    // divergence-theorem volume (`sum a.(b x c)`: negative for a mesh whose
    // normals point into the material). A face is FORWARD when its plane
    // normal (the triangle normal, `BRepBuilderAPI_MakeShapeOnMesh.cxx:238-255`)
    // already points away from the material, and REVERSED otherwise; the wires
    // stay wound with their triangle, so the outer wire is always CCW around the
    // plane normal, which is the invariant `IntTools_FClass2d` (and every UV
    // area) relies on. Relying on the infinite-point classifier alone would
    // leave the material side of the faces dependent on the (unordered) result
    // of that classification. OCCT's class stops at the compound of faces
    // (`cxx:257-261`) and performs none of this.
    let mut signed6 = 0.0;
    for t in &welded.triangles {
        let a = welded.nodes[t.n0].coord;
        let b = welded.nodes[t.n1].coord;
        let c = welded.nodes[t.n2].coord;
        signed6 += a.dot_cross(&b, &c);
    }
    if signed6 < 0.0 {
        for f in faces.iter_mut() {
            f.0.set_orientation(Orientation::Reversed);
        }
    }
    let shell = builder.make_shell(&faces);
    let closed = !edge_use.is_empty() && edge_use.values().all(|&n| n == 2);
    let solid = if closed {
        let mut s = builder.make_solid(&[shell.clone()]).0;
        // Verification step (`BRepTools::OrientClosedSolid`): with the winding
        // above it must report the material already inside and leave the solid
        // alone.
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
