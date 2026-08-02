//! Port of OCCT model healer + post-processor — Wave 4 BRepMesh.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_ModelHealer.{hxx,cxx}` — heals the discretized model.
//! - `BRepMesh_ModelPostProcessor.hxx` + `.cxx` — writes Delaunay triangles back
//!   as model-consumable `Poly_Triangulation` / `Poly_PolygonOnTriangulation`.
//!
//! The OCCT classes operate on the `IMeshData` model's wires and pcurves. In
//! this port the `BRepMeshData` layer (`data_model`) carries no face
//! triangulation store yet, so the healing is applied at the Delaunay level
//! (where the triangles actually live):
//!
//! * [`ModelHealer::heal`] deletes degenerate (zero-area) triangles, unifies
//!   triangle orientation so adjacent triangles share edges in opposite
//!   directions (=> consistent surface normals) and removes orphan nodes.
//! * [`ModelHealer::heal_model`] normalizes the discrete model's edge curves by
//!   dropping consecutive coincident points.
//! * [`ModelPostProcessor::process`] converts a `DelaunDataStructure` into a
//!   [`FaceTriangulation`] (3D nodes + per-node UV + index triples) — the
//!   structure a `MeshModel` face is meshed from, mirroring `Poly_Triangulation`.

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{PCONFUSION, RESOLUTION};

use super::data_model::{MeshCurve, MeshModel, MeshStatus};
use super::delaun_data::DelaunDataStructure;
use super::delaun_types::{DelaunLink, DelaunTriangle, VertexState};

/// A triangle whose UV area is at or below this magnitude is treated as
/// degenerate and removed by the healer. `Resolution` (1e-12) is far below any
/// meaningful feature size while still catching exactly-collinear slivers.
const DEGENERATE_AREA_EPS: f64 = RESOLUTION;

/// Outcome counters of [`ModelHealer::heal`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct HealStats {
    /// Triangles removed for being degenerate (zero-area / repeated vertex).
    pub degenerate_triangles: usize,
    /// Triangles whose winding was flipped to a consistent (CCW) orientation.
    pub flipped_triangles: usize,
    /// Nodes that were not referenced by any live triangle and got removed.
    pub orphan_nodes: usize,
}

/// Heals a triangulated mesh.
///
/// The OCCT `BRepMesh_ModelHealer` fixes wire boundaries and self-intersections
/// on the `IMeshData` model; here the same responsibility is applied to the
/// Delaunay structure where the triangles live. Vertex "welding" is handled at
/// insertion time by the data structure's cell filter, so it is not repeated.
pub struct ModelHealer;

impl ModelHealer {
    /// Deletes degenerate (zero-area) triangles, unifies orientation so every
    /// live triangle is counter-clockwise in UV space (adjacent triangles then
    /// traverse shared edges in opposite directions, giving consistent surface
    /// normals), and removes nodes not referenced by any live triangle.
    pub fn heal(mesh_data: &mut DelaunDataStructure) -> Result<HealStats, String> {
        let mut stats = HealStats::default();

        // 1. Delete degenerate (zero-area / repeated-vertex) triangles.
        let ids: Vec<i32> = mesh_data.elements_of_domain().iter().copied().collect();
        for id in ids {
            let tri = mesh_data.get_element(id);
            if Self::is_degenerate(mesh_data, &tri) {
                mesh_data.remove_element(id);
                stats.degenerate_triangles += 1;
            }
        }

        // 2. Drop links orphaned by the deletion above. Constraint edges
        //    (`Frontier`/`Fixed`) are preserved as boundary of the domain.
        let link_ids: Vec<i32> = mesh_data.links_of_domain().iter().copied().collect();
        for id in link_ids {
            if mesh_data.elements_connected_to(id).is_empty()
                && mesh_data.link_movability(id) == VertexState::Free
            {
                mesh_data.remove_link(id, false);
            }
        }

        // 3. Unify orientation: flip clockwise triangles to counter-clockwise.
        let ids: Vec<i32> = mesh_data.elements_of_domain().iter().copied().collect();
        for id in ids {
            let tri = mesh_data.get_element(id);
            if Self::signed_area(mesh_data, &tri) < 0.0 {
                Self::flip_triangle(mesh_data, id);
                stats.flipped_triangles += 1;
            }
        }

        // 4. Remove nodes not referenced by any live triangle (only fully
        //    isolated ones; `remove_node` refuses nodes with live links, so
        //    boundary nodes stay untouched).
        let live: HashSet<i32> = {
            let mut set = HashSet::new();
            let ids: Vec<i32> = mesh_data.elements_of_domain().iter().copied().collect();
            for id in ids {
                set.extend(mesh_data.get_element(id).vertex_indices.iter().copied());
            }
            set
        };
        let orphans: Vec<i32> = (1..=mesh_data.nb_nodes() as i32)
            .filter(|&n| mesh_data.get_node(n).state != VertexState::Deleted && !live.contains(&n))
            .collect();
        for n in orphans {
            mesh_data.remove_node(n, true);
            if mesh_data.get_node(n).state == VertexState::Deleted {
                stats.orphan_nodes += 1;
            }
        }

        Ok(stats)
    }

    /// Heals a discrete model by dropping consecutive coincident points from
    /// every edge's 3D curve discretization.
    ///
    /// ponytail: the `BRepMeshData` layer exposes no mutable accessor for edge
    /// pcurves (and no face triangulation store), so only the 3D curve is
    /// deduplicated here. Extend when `data_model` grows a pcurve mutator.
    pub fn heal_model(model: &mut MeshModel) -> Result<usize, String> {
        let mut removed = 0usize;
        for edge_index in 0..model.edges_nb() {
            let edge = model.edge_mut(edge_index)?;
            removed += Self::dedup_consecutive_points(edge.discretization_mut())?;
        }
        Ok(removed)
    }

    /// Drops consecutive curve points closer than `Precision::PConfusion`.
    fn dedup_consecutive_points(curve: &mut MeshCurve) -> Result<usize, String> {
        let mut removed = 0usize;
        let mut i = 1usize;
        while i < curve.parameters_nb() {
            let a = curve.get_point(i - 1)?;
            let b = curve.get_point(i)?;
            if a.distance(&b) <= PCONFUSION {
                curve.remove_point(i)?;
                removed += 1;
            } else {
                i += 1;
            }
        }
        Ok(removed)
    }

    /// True when the triangle is degenerate: repeated vertex indices or a
    /// (near-)zero signed area in UV space.
    fn is_degenerate(mesh_data: &DelaunDataStructure, tri: &DelaunTriangle) -> bool {
        let v = tri.vertex_indices;
        if v[0] == v[1] || v[1] == v[2] || v[0] == v[2] {
            return true;
        }
        Self::signed_area(mesh_data, tri).abs() <= DEGENERATE_AREA_EPS
    }

    /// Twice the signed area of the triangle in UV space (positive = CCW).
    fn signed_area(mesh_data: &DelaunDataStructure, tri: &DelaunTriangle) -> f64 {
        let a = mesh_data.get_node(tri.vertex_indices[0]).location;
        let b = mesh_data.get_node(tri.vertex_indices[1]).location;
        let c = mesh_data.get_node(tri.vertex_indices[2]).location;
        (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x())
    }

    /// Rebuilds a triangle with the winding flipped: the last two corners are
    /// swapped and each link is retargeted (with orientation) so link `k` still
    /// joins vertex `k` to vertex `(k+1) % 3` of the new order.
    fn flip_triangle(mesh_data: &mut DelaunDataStructure, element_id: i32) {
        let tri = mesh_data.get_element(element_id);
        let v = tri.vertex_indices;
        let new_verts = [v[0], v[2], v[1]];
        let mut new_links = [0i32; 3];
        for k in 0..3 {
            let (p, q) = (new_verts[k], new_verts[(k + 1) % 3]);
            let link_id = mesh_data.index_of_link(&DelaunLink::new_unparameterized(p, q, 0));
            let link = mesh_data.get_link(link_id);
            // `index_of_link` reports only the (positive) id; re-derive the sign:
            // the stored link is traversed forward exactly when it runs p -> q.
            let sign = if link.first_node() == p && link.last_node() == q { 1 } else { -1 };
            new_links[k] = sign * link_id;
        }
        mesh_data.substitute_element(element_id, DelaunTriangle::new(new_links, new_verts));
    }
}

/// Triangulation of one face produced by [`ModelPostProcessor::process`].
///
/// `nodes[i]` (3D) and `uv_nodes[i]` (UV) are the same mesh vertex; `triangles`
/// reference positions in `nodes`. This is the Rust counterpart of OCCT's
/// `Poly_Triangulation` and is what a `MeshModel` face is meshed from.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceTriangulation {
    /// Model index of the face the triangulation belongs to.
    pub face_index: usize,
    /// 3D positions of the mesh vertices.
    pub nodes: Vec<GpPnt>,
    /// Parametric (UV) position of each mesh vertex, aligned with `nodes`.
    pub uv_nodes: Vec<GpPnt2d>,
    /// Triangle corners as 0-based indices into `nodes`.
    pub triangles: Vec<[i32; 3]>,
}

impl FaceTriangulation {
    /// Empty triangulation (used for faces that were skipped).
    pub fn empty(face_index: usize) -> Self {
        Self { face_index, nodes: Vec::new(), uv_nodes: Vec::new(), triangles: Vec::new() }
    }

    /// Number of mesh vertices.
    pub fn nodes_nb(&self) -> usize {
        self.nodes.len()
    }

    /// Number of triangles.
    pub fn triangles_nb(&self) -> usize {
        self.triangles.len()
    }
}

/// Writes Delaunay triangles back as a model-consumable face triangulation.
pub struct ModelPostProcessor;

impl ModelPostProcessor {
    /// Converts the live Delaunay triangles into a [`FaceTriangulation`] for the
    /// given model face: 3D nodes + per-node UV + index triples.
    ///
    /// Faces marked `Failure`/`Reused` are skipped (empty result), matching the
    /// OCCT post-processor. The Delaunay structure is not modified.
    ///
    /// ponytail: the `BRepMeshData` layer has no face triangulation store to
    /// write into, so the triangulation is returned for the caller to attach
    /// (as OCCT attaches `Poly_Triangulation`). No pcurve index mapping is
    /// written back; `data_model` lacks a mutable pcurve accessor.
    pub fn process(
        delaun: &DelaunDataStructure,
        model: &MeshModel,
        face_index: usize,
    ) -> Result<FaceTriangulation, String> {
        let face = model.face(face_index)?;
        if face.is_status(MeshStatus::FAILURE) || face.is_status(MeshStatus::REUSED) {
            return Ok(FaceTriangulation::empty(face_index));
        }

        let element_ids: Vec<i32> = delaun.elements_of_domain().iter().copied().collect();
        let mut remap: HashMap<i32, usize> = HashMap::with_capacity(element_ids.len() * 2);
        let mut nodes: Vec<GpPnt> = Vec::with_capacity(element_ids.len() * 2);
        let mut uv_nodes: Vec<GpPnt2d> = Vec::with_capacity(element_ids.len() * 2);
        let mut triangles: Vec<[i32; 3]> = Vec::with_capacity(element_ids.len());

        for id in element_ids {
            let tri = delaun.get_element(id);
            let mut corner = [0i32; 3];
            for (k, &node_id) in tri.vertex_indices.iter().enumerate() {
                let pos = *remap.entry(node_id).or_insert_with(|| {
                    let vertex = delaun.get_node(node_id);
                    nodes.push(vertex.p3d);
                    uv_nodes.push(vertex.location);
                    nodes.len() - 1
                });
                corner[k] = pos as i32;
            }
            triangles.push(corner);
        }

        Ok(FaceTriangulation { face_index, nodes, uv_nodes, triangles })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::ShapeType;
    use crate::builder::TopoBuilder;
    use crate::meshing::delaun::Delaun;
    use crate::meshing::delaun_types::DelaunVertex;
    use crate::tgeometry::GeometryRegistry;
    use occt_core::gp::{GpAx3, GpPln};

    fn v(u: f64, w: f64) -> DelaunVertex {
        DelaunVertex::new_parametric(u, w, VertexState::Free)
    }

    fn grid_vertex(u: f64, w: f64) -> DelaunVertex {
        DelaunVertex::new(GpPnt2d::new(u, w), GpPnt::new(u, w, 0.0), 0, VertexState::Free)
    }

    /// Adds the three links of a triangle and registers it, returning its id.
    fn add_triangle(ds: &mut DelaunDataStructure, corners: [i32; 3]) -> i32 {
        let [a, b, c] = corners;
        let ab = ds.add_link(a, b, VertexState::Free);
        let bc = ds.add_link(b, c, VertexState::Free);
        let ca = ds.add_link(c, a, VertexState::Free);
        ds.add_element(DelaunTriangle::new([ab, bc, ca], corners))
    }

    fn signed_area_of(ds: &DelaunDataStructure, id: i32) -> f64 {
        let tri = ds.get_element(id);
        ModelHealer::signed_area(ds, &tri)
    }

    fn clear_tree(s: &crate::shape::TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&crate::shape::TopoShape::from_handle(c));
        }
    }

    #[test]
    fn heal_deletes_degenerate_triangle_and_orphans() {
        let mut ds = DelaunDataStructure::new(16);
        let a = ds.add_node(v(0.0, 0.0));
        let b = ds.add_node(v(1.0, 0.0));
        let c = ds.add_node(v(0.0, 1.0));
        // Degenerate: three collinear points on the U axis => zero area.
        let d = ds.add_node(v(2.0, 0.0));
        let e = ds.add_node(v(3.0, 0.0));
        let f = ds.add_node(v(4.0, 0.0));
        let valid = add_triangle(&mut ds, [a, b, c]);
        let degen = add_triangle(&mut ds, [d, e, f]);
        assert_eq!(ds.elements_of_domain().len(), 2);

        let stats = ModelHealer::heal(&mut ds).expect("heal");

        assert_eq!(stats.degenerate_triangles, 1, "collinear triangle removed");
        assert_eq!(stats.orphan_nodes, 3, "the three collinear nodes become isolated");
        assert_eq!(ds.elements_of_domain().len(), 1);
        assert!(ds.elements_of_domain().contains(&valid));
        assert!(!ds.elements_of_domain().contains(&degen));
        // The valid triangle keeps a positive (CCW) orientation.
        assert!(signed_area_of(&ds, valid) > 0.0);
        // Orphan nodes are marked deleted.
        for n in [d, e, f] {
            assert_eq!(ds.get_node(n).state, VertexState::Deleted);
        }
    }

    #[test]
    fn heal_unifies_orientation_of_adjacent_triangles() {
        let mut ds = DelaunDataStructure::new(16);
        let n0 = ds.add_node(v(0.0, 0.0));
        let n1 = ds.add_node(v(1.0, 0.0));
        let n2 = ds.add_node(v(0.0, 1.0));
        let n3 = ds.add_node(v(1.0, 1.0));
        // Triangle 1 is CCW; triangle 2 (n1,n2,n3) is deliberately CW.
        let t1 = add_triangle(&mut ds, [n0, n1, n2]);
        let t2 = add_triangle(&mut ds, [n1, n2, n3]);
        let shared = ds.index_of_link(&DelaunLink::new_unparameterized(n1, n2, 0));
        assert!(shared > 0);
        assert_eq!(ds.elements_connected_to(shared).extent(), 2, "triangles are adjacent");
        assert!(signed_area_of(&ds, t2) < 0.0, "t2 starts clockwise");

        let stats = ModelHealer::heal(&mut ds).expect("heal");

        assert_eq!(stats.flipped_triangles, 1);
        assert_eq!(ds.elements_of_domain().len(), 2);
        // Both triangles now share the edge in opposite directions => consistent normals.
        let te1 = ds.get_element(t1);
        let te2 = ds.get_element(t2);
        assert!(te1.link_indices.contains(&shared), "t1 traverses shared link forward");
        assert!(te2.link_indices.contains(&-shared), "t2 traverses shared link backward");
        for id in [t1, t2] {
            let tri = ds.get_element(id);
            let area = ModelHealer::signed_area(&ds, &tri);
            assert!(area > 0.0, "every live triangle is counter-clockwise, got {area}");
        }
    }

    #[test]
    fn post_processor_writes_back_3x3_triangulation() {
        // 3x3 lattice: 8 hull vertices (octagon) + 1 interior => 8 triangles.
        let pts: Vec<DelaunVertex> = (0..3)
            .flat_map(|i| (0..3).map(move |j| grid_vertex(i as f64, j as f64)))
            .collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.into_result();
        assert_eq!(ds.elements_of_domain().len(), 8);

        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let model_face = face.clone();
        let mut model = MeshModel::new(crate::shape::TopoShape::new(ShapeType::Compound));
        let face_index = model.add_face(model_face);

        let out = ModelPostProcessor::process(&ds, &model, face_index).expect("process");

        assert_eq!(out.nodes.len(), 9, "all 9 grid vertices written back");
        assert_eq!(out.uv_nodes.len(), 9, "per-node UV aligned with 3D nodes");
        assert_eq!(out.triangles.len(), 8, "8 triangles for the 3x3 grid");
        assert_eq!(out.nodes_nb(), 9);
        assert_eq!(out.triangles_nb(), 8);
        assert_eq!(out.face_index, face_index);
        // Triangle corners index into the compact node array.
        for tri in &out.triangles {
            for &corner in tri {
                assert!((corner as usize) < out.nodes.len());
            }
        }

        clear_tree(&face.0);
        clear_tree(&model.shape().unwrap().clone());
    }

    #[test]
    fn post_processor_skips_failed_face() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let model_face = face.clone();
        let mut model = MeshModel::new(crate::shape::TopoShape::new(ShapeType::Compound));
        let face_index = model.add_face(model_face);
        model.face_mut(face_index).unwrap().set_status(MeshStatus::FAILURE);

        let ds = DelaunDataStructure::new(4);
        let out = ModelPostProcessor::process(&ds, &model, face_index).expect("process");
        assert_eq!(out.nodes_nb(), 0);
        assert_eq!(out.triangles_nb(), 0);

        clear_tree(&face.0);
        clear_tree(&model.shape().unwrap().clone());
    }

    #[test]
    fn heal_model_dedups_consecutive_coincident_curve_points() {
        let b = TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let mut model = MeshModel::new(crate::shape::TopoShape::new(ShapeType::Compound));
        let edge_index = model.add_edge(edge);
        {
            let curve = model.edge_mut(edge_index).unwrap().discretization_mut();
            curve.add_point(GpPnt::new(0.0, 0.0, 0.0), 0.0);
            curve.add_point(GpPnt::new(0.5, 0.0, 0.0), 0.5);
            // Same 3D point twice => consecutive coincident points to drop.
            curve.add_point(GpPnt::new(0.5, 0.0, 0.0), 0.500001);
            curve.add_point(GpPnt::new(1.0, 0.0, 0.0), 1.0);
        }
        assert_eq!(model.edge(edge_index).unwrap().discretization().parameters_nb(), 4);

        let removed = ModelHealer::heal_model(&mut model).expect("heal_model");

        assert_eq!(removed, 1);
        let curve = model.edge(edge_index).unwrap().discretization();
        assert_eq!(curve.parameters_nb(), 3);
        assert_eq!(curve.points().len(), 3);

        clear_tree(&model.shape().unwrap().clone());
    }
}
