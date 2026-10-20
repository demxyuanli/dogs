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
//! ponytail: interior-node generation used to live only on the range-splitter
//! sibling; this algo now owns the splitter for `Scale` / `AddPoint` /
//! `GenerateSurfaceNodes` (`BRepMesh_NodeInsertionMeshAlgo.hxx`).

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d, GpXY};
use occt_core::precision::SQUARE_CONFUSION;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::abs::Orientation;

use super::data_model::*;
use super::delaun::Delaun;
use super::delaun_data::DelaunDataStructure;
use super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::face_discret::Classifier;
// ponytail: `mesh_algo.rs` (BRepMesh_{BaseMeshAlgo,ConstrainedBaseMeshAlgo}) is a
// stub; the base-algo contract it will provide is modelled here by
// `NodeInsertionMeshAlgo` + `DelaunayNodeInsertionMeshAlgo` until it lands.
#[allow(unused_imports)]
use super::mesh_algo::*;
use super::face_discret::PointState;
use super::deflection_control::{
    DelaunayDeflectionControlMeshAlgo,
};
use super::parameters::MeshParameters;
use super::range_splitter::{
    mesh_surface_type, create_range_splitter, factory_uses_deflection_control,
    grab_params_of_internal_edges, RangeSplitter, MeshSurfaceType,
};

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
    fn node_point_2d(&self, vertex: &DelaunVertex) -> GpPnt2d;

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

/// TEMPORARY diagnostic (DELDIAG): dumps the result triangles of the face whose
/// structure contains the node at (13.0176, -24.7899, 290).
pub fn diag_dump_result(result: &TriangulationResult, face_index: usize) {
    if std::env::var("DELDIAG").is_err() {
        return;
    }
    let is_target = result.nodes.iter().any(|n| {
        (n.p3d.x() - 13.017621944).abs() < 1e-3
            && (n.p3d.y() + 24.789947988).abs() < 1e-3
            && (n.p3d.z() - 290.0).abs() < 1e-3
    });
    if !is_target {
        return;
    }
    eprintln!(
        "DELDIAGRES face={face_index} nodes={} tris={}",
        result.nodes.len(),
        result.triangles.len()
    );
    for (i, n) in result.nodes.iter().enumerate() {
        eprintln!(
            "DELDIAGRES n i={} state={} uv=({:.6},{:.6}) p=({:.6},{:.6},{:.6})",
            i + 1,
            n.state.to_str(),
            n.location.x(),
            n.location.y(),
            n.p3d.x(),
            n.p3d.y(),
            n.p3d.z()
        );
    }
    for (i, t) in result.triangles.iter().enumerate() {
        let v = t.vertex_indices;
        let p = |k: i32| -> &DelaunVertex { &result.nodes[(k - 1) as usize] };
        let (a, b, c) = (p(v[0]), p(v[1]), p(v[2]));
        eprintln!(
            "DELDIAGRES t i={i} v=({},{},{}) uv=[({:.6},{:.6}),({:.6},{:.6}),({:.6},{:.6})]",
            v[0],
            v[1],
            v[2],
            a.location.x(),
            a.location.y(),
            b.location.x(),
            b.location.y(),
            c.location.x(),
            c.location.y()
        );
    }
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
    /// Point-in-face classifier over the face's registered wires
    /// (`BRepMesh_Classifier`), used to drop surface/interior nodes that fall
    /// outside the face's real region.
    classifier: Classifier,
    pre_process_surface_nodes: bool,
    last_result: Option<TriangulationResult>,
    /// Set when [`Delaun::failed`] aborted `generateMesh`: the face must keep
    /// `IMeshData_Failure` and no triangulation (`BRepMesh_BaseMeshAlgo.cxx:59-62`).
    delaun_failed: bool,
    /// Face range splitter: `Scale` to face basis, `AddPoint`, interior grid.
    /// Source: `BRepMesh_NodeInsertionMeshAlgo::myRangeSplitter`.
    splitter: Option<Box<dyn RangeSplitter>>,
    /// Face deflection written back after `optimizeMesh`
    /// (`BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:117-120`).
    optimized_face_deflection: Option<f64>,
}

impl DelaunayNodeInsertionMeshAlgo {
    /// `BRepMesh_DelaunayNodeInsertionMeshAlgo` constructor: `PreProcessSurfaceNodes` is false.
    /// `BRepMesh_MeshAlgoFactory` never sets the flag; only `DelabellaMeshAlgoFactory` does.
    pub fn new() -> Self {
        Self::with_options(false)
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
            classifier: Classifier::new(),
            pre_process_surface_nodes,
            last_result: None,
            delaun_failed: false,
            splitter: None,
            optimized_face_deflection: None,
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
}

mod algo_impl;
mod build;
mod debug_dump;
mod links;
mod mesh_finish;
mod parametric;
mod surface_nodes;
#[cfg(test)]
mod tests;

