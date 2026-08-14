//! Port of OCCT BRepMesh node insertion — Wave 3 mesh algorithms.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_NodeInsertionMeshAlgo.hxx` — template mixin enabling insertion of
//!   free vertices into the mesh: boundary nodes from discretized wires/pcurves
//!   plus internal face vertices.
//! - `BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx` — Delaunay node insertion:
//!   feeds the face's discretized boundary UV points (and, optionally, generated
//!   surface nodes) into a [`Delaun`] triangulator, marks the constraint edges,
//!   and collects the resulting triangles.
//!
//! The OCCT classes are templates over a `RangeSplitter` + `BaseAlgo` (the factory
//! instantiates `BRepMesh_DelaunayNodeInsertionMeshAlgo<RangeSplitter,
//! BRepMesh_DelaunayBaseMeshAlgo>`). The Rust port collapses the mixin chain into
//! a [`NodeInsertionMeshAlgo`] trait + a concrete [`DelaunayNodeInsertionMeshAlgo`]
//! struct.
//!
//! ponytail: the `RangeSplitter` UV normalization/scaling and the deflection-driven
//! `GenerateSurfaceNodes` interior-point generation belong to the sibling
//! `mesh_algo.rs` range-splitter port (still a stub), so both are omitted here:
//! nodes are stored at their raw UV coordinates and surface-node generation
//! currently yields no interior points.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::geom_api::project_point_on_surface;

use super::data_model::*;
use super::delaun::Delaun;
use super::delaun_data::DelaunDataStructure;
use super::delaun_types::{DelaunTriangle, DelaunVertex, VertexState};
// ponytail: `mesh_algo.rs` (BRepMesh_{BaseMeshAlgo,ConstrainedBaseMeshAlgo}) is a
// stub; the base-algo contract it will provide is modelled here by
// `NodeInsertionMeshAlgo` + `DelaunayNodeInsertionMeshAlgo` until it lands.
#[allow(unused_imports)]
use super::mesh_algo::*;
use super::parameters::MeshParameters;
use super::range_splitter::{create_range_splitter, RangeSplitter};

/// Abstract insertion of free vertices into a Delaunay mesh.
/// Source: `BRepMesh_NodeInsertionMeshAlgo`.
pub trait NodeInsertionMeshAlgo {
    /// Adds the given 2D point to the mesh data structure and returns its
    /// (1-based) node index. `p3d` is the associated 3D point (inlined on the
    /// vertex, matching this port's `DelaunVertex`); coincident points collapse
    /// to the existing node unless `is_force_add` is set. Source:
    /// `addNodeToStructure`.
    fn add_node_to_structure(
        &mut self,
        point: GpPnt2d,
        p3d: GpPnt,
        movability: VertexState,
        is_force_add: bool,
    ) -> i32;

    /// Returns the 2D point associated with the given vertex.
    /// Source: `getNodePoint2d`.
    fn get_node_point_2d(&self, vertex: &DelaunVertex) -> GpPnt2d;

    /// Registers the boundary UV points of a face as `Frontier` nodes and adds
    /// the given constraint edges — `(first, last)` positions into `uv_points` —
    /// as `Frontier` links. Returns the number of registered boundary points.
    /// Source: `BRepMesh_BaseMeshAlgo::initDataStructure` (wire/pcurve node
    /// registration) + `addLinkToMesh`.
    fn insert_boundary_nodes(
        &mut self,
        uv_points: &[GpPnt2d],
        p3d_points: &[GpPnt],
        constraint_edges: &[(i32, i32)],
    ) -> Result<usize, String>;

    /// Registers internal (in-face) UV points as `Fixed` nodes.
    /// Returns the number of registered points.
    /// Source: `insertInternalVertices` / `insertInternalVertex`.
    fn insert_internal_nodes(&mut self, uv_points: &[GpPnt2d]) -> Result<usize, String>;
}

/// Result of a node-insertion triangulation run.
#[derive(Debug, Clone, PartialEq)]
pub struct TriangulationResult {
    /// Triangles of the triangulated domain, referencing 1-based structure node
    /// indices (as in the OCCT `BRepMesh_DataStructureOfDelaun`).
    pub triangles: Vec<DelaunTriangle>,
    /// Compact list of live nodes: `nodes[i - 1]` is the node with structure
    /// index `i`, so a triangle's vertex indices index directly into this list.
    pub nodes: Vec<DelaunVertex>,
}

/// Delaunay node insertion: builds a 2D UV triangulation of a face's discretized
/// boundary (plus optional interior nodes) with [`Delaun`] and collects the result.
/// Source: `BRepMesh_DelaunayNodeInsertionMeshAlgo`.
pub struct DelaunayNodeInsertionMeshAlgo {
    structure: DelaunDataStructure,
    nodes_map: Vec<GpPnt>,
    used_nodes: HashMap<i32, i32>,
    boundary_indices: Vec<i32>,
    boundary_uv: Vec<GpPnt2d>,
    pre_process_surface_nodes: bool,
    last_result: Option<TriangulationResult>,
}

impl DelaunayNodeInsertionMeshAlgo {
    /// Creates an empty insertion algorithm (surface nodes pre-processed).
    pub fn new() -> Self {
        Self::with_options(true)
    }

    /// Creates an insertion algorithm controlling whether generated surface
    /// (interior) nodes are registered before triangulation.
    pub fn with_options(pre_process_surface_nodes: bool) -> Self {
        Self {
            structure: DelaunDataStructure::new(16),
            nodes_map: Vec::new(),
            used_nodes: HashMap::new(),
            boundary_indices: Vec::new(),
            boundary_uv: Vec::new(),
            pre_process_surface_nodes,
            last_result: None,
        }
    }

    /// `PreProcessSurfaceNodes` flag.
    pub fn is_pre_process_surface_nodes(&self) -> bool {
        self.pre_process_surface_nodes
    }

    /// Sets the `PreProcessSurfaceNodes` flag.
    pub fn set_pre_process_surface_nodes(&mut self, value: bool) {
        self.pre_process_surface_nodes = value;
    }

    /// The mesh data structure being built (nodes + constraint links).
    pub fn structure(&self) -> &DelaunDataStructure {
        &self.structure
    }

    /// The 3D nodes map, parallel to the structure's node indices.
    pub fn nodes_map(&self) -> &[GpPnt] {
        &self.nodes_map
    }

    /// Identity map of structure node indices that are used (non-`Free`), i.e.
    /// nodes that must survive into the output triangulation. Source:
    /// `BRepMesh_BaseMeshAlgo::myUsedNodes`.
    pub fn used_nodes(&self) -> &HashMap<i32, i32> {
        &self.used_nodes
    }

    /// Structure node index registered for each boundary UV point, in input order.
    pub fn boundary_indices(&self) -> &[i32] {
        &self.boundary_indices
    }

    /// Result of the last triangulation run, if any.
    pub fn last_result(&self) -> Option<&TriangulationResult> {
        self.last_result.as_ref()
    }

    /// Registers a point in the 3D nodes map and adds its 2D projection to the
    /// structure. Returns the structure node index. Source:
    /// `BRepMesh_BaseMeshAlgo::registerNode`.
    pub fn register_node(
        &mut self,
        point: GpPnt,
        point2d: GpPnt2d,
        movability: VertexState,
        is_force_add: bool,
    ) -> i32 {
        let node_index =
            self.add_node_to_structure(point2d, point, movability, is_force_add);
        if node_index as usize > self.nodes_map.len() {
            self.nodes_map.push(point);
            // Pre-bind frontier/fixed nodes with identity mapping so boundary
            // node indices stay stable; Free (internal) nodes enter only when
            // actually referenced by a triangle.
            if movability != VertexState::Free {
                self.used_nodes.insert(node_index, node_index);
            }
        }
        node_index
    }

    /// Generates surface (interior) nodes for the face and registers them.
    /// Source: `DelaunayNodeInsertionMeshAlgo::registerSurfaceNodes`, driving the
    /// analytical `RangeSplitter` (cylinder/cone/sphere/torus/NURBS) to place the
    /// interior UV nodes that a periodic surface's degenerate seam boundary alone
    /// cannot cover.
    pub fn generate_surface_nodes(
        &mut self,
        model: &MeshModel,
        face_index: usize,
        params: &MeshParameters,
    ) -> Result<usize, String> {
        let face = model.face(face_index)?;
        let Some(surface) = face.surface() else {
            return Ok(0);
        };
        let mut splitter = create_range_splitter(surface.as_ref());
        splitter.reset(face, params);
        for &uv in &self.boundary_uv {
            splitter.add_point(uv);
        }
        splitter.adjust_range();
        let Some(nodes) = splitter.generate_surface_nodes(params) else {
            return Ok(0);
        };
        let count = nodes.len();
        self.insert_internal_nodes(&nodes)?;
        Ok(count)
    }

    /// Triangulates the registered nodes (boundary + internal + generated
    /// surface) and returns the triangles. The structure is consumed by the
    /// [`Delaun`]; call this once after inserting all nodes.
    pub fn triangulate(&mut self) -> Result<TriangulationResult, String> {
        if self.structure.nb_nodes() == 0 {
            return Err("DelaunayNodeInsertionMeshAlgo::triangulate: no nodes registered".to_string());
        }

        let structure = std::mem::replace(&mut self.structure, DelaunDataStructure::new(0));
        let mut indices: Vec<i32> = (1..=structure.nb_nodes() as i32)
            .filter(|&i| structure.get_node(i).state != VertexState::Deleted)
            .collect();
        let mut delaun = Delaun::new_with_data(structure, &mut indices);

        // Constraint edges. OCCT's `BRepMesh_Delaun::UseEdge` is an empty stub
        // (the Rust port keeps that faithfully) — the constraint edges are
        // enforced through the Frontier/Fixed links pre-registered in the
        // structure, which `process_constraints` honours during triangulation.
        let frontier: Vec<i32> = delaun.frontier().into_iter().collect();
        for &e in &frontier {
            let _ = delaun.use_edge(e);
        }

        let ds = delaun.result();
        let triangles: Vec<DelaunTriangle> = ds
            .elements_of_domain()
            .iter()
            .map(|&id| ds.get_element(id))
            .collect();
        let nodes: Vec<DelaunVertex> = (1..=ds.nb_nodes() as i32)
            .filter(|&i| ds.get_node(i).state != VertexState::Deleted)
            .map(|i| *ds.get_node(i))
            .collect();

        let result = TriangulationResult { triangles, nodes };
        self.last_result = Some(result.clone());
        Ok(result)
    }

    /// Performs node insertion for a model face: extracts the boundary UV points
    /// and constraint edges from the face's wires → edges → pcurves (populated by
    /// `edge_discret`/`face_discret`), registers them together with any internal
    /// face points, triangulates, and returns the triangles.
    ///
    /// ponytail: `MeshModel` has no triangle slot yet — the committing
    /// `collectTriangles`/`commitSurfaceTriangulation` belong to the stub
    /// `mesh_algo.rs` module — so the triangulation is returned and stored in
    /// [`Self::last_result`] for the caller to commit.
    pub fn perform(
        &mut self,
        model: &mut MeshModel,
        face_index: usize,
        params: &MeshParameters,
    ) -> Result<TriangulationResult, String> {
        let (uv, p3d, constraints) = Self::collect_boundary_uv(model, face_index)?;
        if uv.is_empty() {
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index} has no boundary UV points"
            ));
        }
        self.boundary_uv = uv.clone();
        self.insert_boundary_nodes(&uv, &p3d, &constraints)?;

        // Internal (in-face) 3D points -> UV via the face surface.
        let mut internal_uv: Vec<GpPnt2d> = Vec::new();
        let face = model.face(face_index)?;
        if let Some(surface) = face.surface() {
            for p in face.points() {
                if let Some(proj) = project_point_on_surface(surface.as_ref(), p, 1e-7) {
                    internal_uv.push(GpPnt2d::new(proj.u, proj.v));
                }
            }
        }
        self.insert_internal_nodes(&internal_uv)?;

        if self.pre_process_surface_nodes {
            self.generate_surface_nodes(model, face_index, params)?;
        }

        self.triangulate()
    }

    /// Collects the boundary UV points and constraint edges of a model face from
    /// its wires → edges → pcurves (in wire order, honoring each pcurve
    /// orientation). Each wire's concatenated pcurve points form a closed chain;
    /// consecutive points (plus the closing edge) become the constraint edges.
    fn collect_boundary_uv(
        model: &MeshModel,
        face_index: usize,
    ) -> Result<(Vec<GpPnt2d>, Vec<GpPnt>, Vec<(i32, i32)>), String> {
        let face = model.face(face_index)?;
        let mut uv: Vec<GpPnt2d> = Vec::new();
        let mut p3d: Vec<GpPnt> = Vec::new();
        let mut constraints: Vec<(i32, i32)> = Vec::new();

        for &wire_index in face.wires() {
            let wire = model.wire(wire_index)?;
            let chain_start = uv.len();
            for j in 0..wire.edges_nb() {
                let edge_index = wire.edge(j)?;
                let orientation = wire.edge_orientation(j)?;
                let edge = model.edge(edge_index)?;
                let Some(pcurve) = edge.pcurve_for(face_index, orientation) else {
                    continue;
                };
                // The pcurve (2D) and the edge's shared 3D polyline are discretized
                // at the same parameters, so their i-th points coincide; the wire
                // orientation reverses both together (`BRep_Tool::CurveOnSurface`
                // + the shared `BRepMeshData_Edge` polyline).
                let pts = pcurve.points();
                let pts3d = edge.discretization().points();
                if pcurve.is_forward() {
                    uv.extend_from_slice(pts);
                    p3d.extend_from_slice(pts3d);
                } else {
                    uv.extend(pts.iter().rev());
                    p3d.extend(pts3d.iter().rev());
                }
            }

            let n = uv.len();
            if n - chain_start >= 2 {
                for t in chain_start..n - 1 {
                    if uv[t] != uv[t + 1] {
                        constraints.push((t as i32, (t + 1) as i32));
                    }
                }
                if uv[chain_start] != uv[n - 1] {
                    constraints.push(((n - 1) as i32, chain_start as i32));
                }
            }
        }

        Ok((uv, p3d, constraints))
    }
}

impl NodeInsertionMeshAlgo for DelaunayNodeInsertionMeshAlgo {
    fn add_node_to_structure(
        &mut self,
        point: GpPnt2d,
        p3d: GpPnt,
        movability: VertexState,
        is_force_add: bool,
    ) -> i32 {
        let vertex = DelaunVertex::new(point, p3d, 0, movability);
        if is_force_add {
            self.structure.add_node_force(vertex)
        } else {
            self.structure.add_node(vertex)
        }
    }

    fn get_node_point_2d(&self, vertex: &DelaunVertex) -> GpPnt2d {
        // ponytail: no range splitter, so the stored parametric location is
        // returned verbatim (OCCT un-scales it via the splitter).
        vertex.location
    }

    fn insert_boundary_nodes(
        &mut self,
        uv_points: &[GpPnt2d],
        p3d_points: &[GpPnt],
        constraint_edges: &[(i32, i32)],
    ) -> Result<usize, String> {
        self.boundary_indices.clear();
        self.boundary_indices.reserve(uv_points.len());
        for (i, &p) in uv_points.iter().enumerate() {
            // The boundary 3D point is the shared edge polyline vertex — not a
            // placeholder — so adjacent faces meet exactly on the shared edge.
            let node = self.register_node(p3d_points[i], p, VertexState::Frontier, false);
            self.boundary_indices.push(node);
        }

        for &(a, b) in constraint_edges {
            let na = *self.boundary_indices.get(a as usize).ok_or_else(|| {
                format!(
                    "DelaunayNodeInsertionMeshAlgo::insert_boundary_nodes: constraint start {a} out of range ({} points)",
                    uv_points.len()
                )
            })?;
            let nb = *self.boundary_indices.get(b as usize).ok_or_else(|| {
                format!(
                    "DelaunayNodeInsertionMeshAlgo::insert_boundary_nodes: constraint end {b} out of range ({} points)",
                    uv_points.len()
                )
            })?;
            if na != nb {
                self.structure.add_link(na, nb, VertexState::Frontier);
            }
        }

        Ok(uv_points.len())
    }

    fn insert_internal_nodes(&mut self, uv_points: &[GpPnt2d]) -> Result<usize, String> {
        for &p in uv_points {
            self.register_node(GpPnt::zero(), p, VertexState::Fixed, false);
        }
        Ok(uv_points.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshing::model_builder::ModelBuilder;
    use crate::primitives::BRepPrimBox;
    use crate::tgeometry::GeometryRegistry;

    /// The 3x3 lattice of UV points.
    fn grid_3x3() -> Vec<GpPnt2d> {
        (0..3)
            .flat_map(|i| (0..3).map(move |j| GpPnt2d::new(i as f64, j as f64)))
            .collect()
    }

    /// Octagonal boundary (counter-clockwise, positive area) of the 3x3 grid's
    /// convex hull plus its closing loop.
    fn octagon() -> (Vec<GpPnt2d>, Vec<(i32, i32)>) {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(2.0, 0.0),
            GpPnt2d::new(2.0, 1.0),
            GpPnt2d::new(2.0, 2.0),
            GpPnt2d::new(1.0, 2.0),
            GpPnt2d::new(0.0, 2.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let edges = (0..8).map(|i| (i as i32, ((i + 1) % 8) as i32)).collect();
        (pts, edges)
    }

    #[test]
    fn grid_3x3_uv_points_inserted_into_delaun_produce_triangles() {
        let pts = grid_3x3();
        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        assert_eq!(algo.insert_internal_nodes(&pts).unwrap(), 9);

        let res = algo.triangulate().expect("triangulate");
        assert!(!res.triangles.is_empty());
        // 9 lattice points, octagonal hull (h = 8): 2N - 2 - h = 8 triangles.
        assert_eq!(res.triangles.len(), 8);
        assert_eq!(res.nodes.len(), 9);
        // Every triangle references live node indices only.
        for t in &res.triangles {
            for &v in &t.vertex_indices {
                assert!((1..=9).contains(&v), "triangle references node {v}");
            }
        }
    }

    #[test]
    fn boundary_nodes_registered_as_frontier_with_constraint_links() {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let edges = vec![(0i32, 1), (1, 2), (2, 3), (3, 0)];

        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        assert_eq!(algo.insert_boundary_nodes(&pts, &vec![GpPnt::zero(); 4], &edges).unwrap(), 4);
        assert_eq!(algo.boundary_indices(), &[1, 2, 3, 4][..]);

        let ds = algo.structure();
        assert_eq!(ds.nb_nodes(), 4);
        for i in 1..=4 {
            assert_eq!(ds.get_node(i).state, VertexState::Frontier);
        }
        // The four constraint edges are registered as Frontier links.
        assert_eq!(ds.links_of_domain().len(), 4);
        for &l in ds.links_of_domain() {
            assert_eq!(ds.link_movability(l), VertexState::Frontier);
        }
    }

    #[test]
    fn internal_insertion_yields_2n_minus_2_minus_h() {
        let (boundary, edges) = octagon();
        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        assert_eq!(algo.insert_boundary_nodes(&boundary, &vec![GpPnt::zero(); 8], &edges).unwrap(), 8);
        assert_eq!(algo.insert_internal_nodes(&[GpPnt2d::new(1.0, 1.0)]).unwrap(), 1);

        let res = algo.triangulate().expect("triangulate");
        let n = 9; // 8 boundary + 1 internal point
        let h = 8; // octagonal hull
        assert_eq!(res.triangles.len() as i32, 2 * n - 2 - h);
        assert_eq!(res.triangles.len(), 8);
    }

    #[test]
    fn perform_errors_on_face_without_boundary_pcurves() {
        let shape = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone();
        let params = MeshParameters::default();
        let mut model = ModelBuilder::build_model(&shape, &params).expect("model built");
        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        // Pcurves are empty until the edge discretizers run, so no boundary UV
        // points can be collected and the face must be rejected.
        let err = algo.perform(&mut model, 0, &params);
        assert!(err.is_err(), "unpopulated pcurves yield no boundary UV points");
        GeometryRegistry::global().clear_shape(&shape);
    }
}
