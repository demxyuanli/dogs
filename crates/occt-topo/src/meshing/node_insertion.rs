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
use super::deflection_control::{
    DelaunayDeflectionControlMeshAlgo, DeflectionMesh,
};
use super::parameters::MeshParameters;
use super::range_splitter::{
    classify_surface, create_range_splitter, factory_uses_deflection_control,
    grab_params_of_internal_edges, RangeSplitter, SurfaceType,
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

    /// Maps a UV point into the Delaunay face basis. Source:
    /// `BRepMesh_NodeInsertionMeshAlgo::addNodeToStructure` /
    /// `BRepMesh_DefaultRangeSplitter::Scale(thePoint, true)`.
    fn scale_to_face(&self, point: GpPnt2d) -> GpPnt2d {
        let Some(s) = self.splitter.as_ref() else {
            return point;
        };
        let (du, dv) = s.delta();
        if du.abs() <= 1e-16 || dv.abs() <= 1e-16 {
            return point;
        }
        s.scale(point, true)
    }

    /// Maps a face-basis point back to surface UV. Source:
    /// `BRepMesh_NodeInsertionMeshAlgo::getNodePoint2d` /
    /// `Scale(thePoint, false)`.
    fn scale_from_face(&self, point: GpPnt2d) -> GpPnt2d {
        let Some(s) = self.splitter.as_ref() else {
            return point;
        };
        s.scale(point, false)
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
        let nodes = self.list_surface_nodes(model, face_index, params)?;
        let count = nodes.len();
        for (uv, p3d) in nodes {
            self.register_node(p3d, uv, VertexState::Free, false);
        }
        Ok(count)
    }

    /// Surface nodes classified `IN` (`registerSurfaceNodes` / `insertNodes` loop).
    /// Source: `BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:115-123, 144-152`.
    fn list_surface_nodes(
        &mut self,
        model: &MeshModel,
        face_index: usize,
        params: &MeshParameters,
    ) -> Result<Vec<(GpPnt2d, GpPnt)>, String> {
        let face = model.face(face_index)?;
        let Some(surface) = face.surface() else {
            return Ok(Vec::new());
        };
        let nodes = if let Some(splitter) = self.splitter.as_mut() {
            grab_params_of_internal_edges(model, face_index, splitter.as_mut());
            splitter.generate_surface_nodes(params)
        } else {
            let mut splitter = create_range_splitter(surface.as_ref());
            splitter.reset(face, params);
            for &uv in &self.boundary_uv {
                splitter.add_point(uv);
            }
            splitter.adjust_range();
            grab_params_of_internal_edges(model, face_index, splitter.as_mut());
            splitter.generate_surface_nodes(params)
        };
        let Some(nodes) = nodes else {
            return Ok(Vec::new());
        };
        let mut out = Vec::new();
        for n in nodes {
            if self.classifier.is_inside(&n) {
                out.push((n, surface.d0(n.x(), n.y())));
            }
        }
        Ok(out)
    }

    /// Triangulates the registered nodes (boundary + internal + generated
    /// surface) and returns the triangles. The structure is consumed by the
    /// [`Delaun`]; call this once after inserting all nodes.
    pub fn triangulate(&mut self) -> Result<TriangulationResult, String> {
        self.finish_mesh(&[], &MeshParameters::default(), 0.0, usize::MAX)
    }

    /// Base mesh, then optional `insertNodes` (`AddVertices`), then `collectTriangles`.
    /// Source: `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` +
    /// `DelaunayNodeInsertionMeshAlgo::postProcessMesh` +
    /// `BRepMesh_BaseMeshAlgo::collectTriangles`.
    fn finish_mesh(
        &mut self,
        insert: &[(GpPnt2d, GpPnt)],
        params: &MeshParameters,
        face_deflection: f64,
        face_index: usize,
    ) -> Result<TriangulationResult, String> {
        if self.structure.nb_nodes() == 0 {
            return Err("DelaunayNodeInsertionMeshAlgo::triangulate: no nodes registered".to_string());
        }

        let structure = std::mem::replace(&mut self.structure, DelaunDataStructure::new(0));
        let mut indices: Vec<i32> = (1..=structure.nb_nodes() as i32)
            .filter(|&i| structure.get_node(i).state != VertexState::Deleted)
            .collect();
        let (cells_u, cells_v) = self.cells_count(indices.len() as i32);
        let mut delaun = Delaun::new_with_data_cells(structure, &mut indices, cells_u, cells_v);
        // `BRepMesh_Delaun::addTriangle` overflowed a link's triangle pair:
        // OCCT's `Standard_OutOfRange` (`BRepMesh_PairOfIndex.hxx:41`) unwinds out
        // of `generateMesh`, `BRepMesh_BaseMeshAlgo::Perform` swallows it
        // (`BRepMesh_BaseMeshAlgo.cxx:59-62`) and `commitSurfaceTriangulation` is
        // never reached, so the face ends up with `IMeshData_Failure` and no
        // triangulation at all — not with the partial mesh built so far. Record
        // the flag and stop here; `perform` turns it into the face's FAILURE
        // status (audit A25 / task T-61).
        if delaun.failed() {
            self.delaun_failed = true;
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::finish_mesh: BRepMesh_Delaun::addTriangle failed (link pair overflow)"
            ));
        }
        // `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` (`cxx:45-46`).
        delaun.erase_free_links();

        // Constraint edges. OCCT's `BRepMesh_Delaun::UseEdge` is an empty stub
        // (the Rust port keeps that faithfully) — the constraint edges are
        // enforced through the Frontier/Fixed links pre-registered in the
        // structure, which `process_constraints` honours during triangulation.
        let frontier: Vec<i32> = delaun.frontier().into_iter().collect();
        for &e in &frontier {
            let _ = delaun.use_edge(e);
        }
        let bnd_nodes = delaun.result().nb_nodes();
        let bnd_tris = delaun.result().elements_of_domain().len();

        if !insert.is_empty() {
            let mut idxs: Vec<i32> = Vec::with_capacity(insert.len());
            for &(uv, p3d) in insert {
                let vertex = DelaunVertex::new(self.scale_to_face(uv), p3d, 0, VertexState::Free);
                let idx = delaun.add_node(vertex);
                if idx as usize > self.nodes_map.len() {
                    self.nodes_map.push(p3d);
                }
                idxs.push(idx);
            }
            delaun.add_vertices(&mut idxs);
        }

        // `DelaunayDeflectionControlMeshAlgo::postProcessMesh` (`hxx:50-68`):
        // after NodeInsertion surface nodes, `optimizeMesh` when the factory
        // selected DeflectionControl and `ControlSurfaceDeflection` is on.
        if params.control_surface_deflection {
            if let Some(surf) = self.splitter.as_ref().and_then(|s| s.surface().cloned()) {
                if factory_uses_deflection_control(surf.as_ref(), params)
                    && !delaun.result().elements_of_domain().is_empty()
                {
                    self.optimize_mesh(&mut delaun, surf.as_ref(), params, face_deflection);
                }
            }
        }

        self.collect_triangles(delaun)
    }

    /// `BRepMesh_DelaunayDeflectionControlMeshAlgo::optimizeMesh` + `insertNodes`.
    fn optimize_mesh(
        &mut self,
        delaun: &mut Delaun,
        surface: &dyn Surface,
        params: &MeshParameters,
        face_deflection: f64,
    ) {
        let (ru, rv, delta) = match self.splitter.as_ref() {
            Some(sp) => (sp.range_u(), sp.range_v(), sp.delta()),
            None => return,
        };
        let mut ctrl = DelaunayDeflectionControlMeshAlgo::new(params);
        let mut view = ParametricDelaun { delaun, ru, rv, delta };
        ctrl.post_process(&mut view, surface, face_deflection, &mut |mesh, nodes| {
            mesh.insert_parametric(self, surface, nodes)
        });
    }

    /// `BRepMesh_DelaunayNodeInsertionMeshAlgo::getCellsCount` /
    /// `BRepMesh_GeomTool::CellsCount` (`BRepMesh_GeomTool.cxx:465-511`).
    fn cells_count(&self, vertices_nb: i32) -> (i32, i32) {
        let Some(sp) = self.splitter.as_ref() else {
            return (-1, -1);
        };
        let Some(surf) = sp.surface() else {
            return (-1, -1);
        };
        geom_tool_cells_count(
            surf.as_ref(),
            vertices_nb,
            sp.deflection(),
            sp.range_u(),
            sp.range_v(),
            sp.delta(),
        )
    }

    /// `BRepMesh_BaseMeshAlgo::collectTriangles` (`BRepMesh_BaseMeshAlgo.cxx:245-277`).
    /// Frontier/Fixed are pre-bound to their structure index; Free enter only
    /// when a domain triangle references them (`ElementNodes` from links, not
    /// unused Free and not cached `vertex_indices`).
    fn collect_triangles(&mut self, delaun: Delaun) -> Result<TriangulationResult, String> {
        let ds = delaun.result();
        let mut used = self.used_nodes.clone();
        let mut triangles: Vec<DelaunTriangle> = Vec::new();
        for &id in ds.elements_of_domain() {
            let mut elem = ds.get_element(id);
            let mut a = ds.element_nodes(&elem);
            for n in &mut a {
                if !used.contains_key(n) {
                    let next = used.len() as i32 + 1;
                    used.insert(*n, next);
                }
                *n = *used.get(n).expect("collectTriangles: bound node");
            }
            elem.vertex_indices = a;
            triangles.push(elem);
        }

        let n_out = used.values().copied().max().unwrap_or(0) as usize;
        let mut slots: Vec<Option<DelaunVertex>> = vec![None; n_out];
        for i in 1..=ds.nb_nodes() as i32 {
            if ds.get_node(i).state == VertexState::Deleted {
                continue;
            }
            if let Some(&out) = used.get(&i) {
                let mut n = *ds.get_node(i);
                n.location = self.scale_from_face(n.location);
                let slot = (out as usize).saturating_sub(1);
                if slot < slots.len() {
                    slots[slot] = Some(n);
                }
            }
        }
        let nodes: Vec<DelaunVertex> = slots
            .into_iter()
            .map(|n| {
                n.unwrap_or_else(|| {
                    DelaunVertex::new_parametric(0.0, 0.0, VertexState::Deleted)
                })
            })
            .collect();
        self.used_nodes = used;

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
        let (uv, _p3d, _constraints, wires_uv) = Self::collect_boundary_uv(model, face_index)?;
        if uv.is_empty() {
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index} has no boundary UV points"
            ));
        }
        self.boundary_uv = uv.clone();

        // `BRepMesh_NodeInsertionMeshAlgo::initDataStructure`: Reset splitter,
        // AddPoint from collectWirePoints, AdjustRange, SetCellSize/SetTolerance,
        // then RegisterWire with the splitter range.
        self.splitter = None;
        self.classifier = Classifier::new();
        let mut range_invalid = false;
        {
            let face = model.face(face_index)?;
            if let Some(surface) = face.surface() {
                let mut sp = create_range_splitter(surface.as_ref());
                sp.reset(face, params);
                for w in &wires_uv {
                    for &p in w {
                        sp.add_point(p);
                    }
                }
                sp.adjust_range();
                if sp.is_valid() {
                    let (du, dv) = sp.delta();
                    let (tu, tv) = sp.tolerance_uv();
                    if du.abs() > 1e-16 && dv.abs() > 1e-16 {
                        self.structure
                            .set_cell_size(14.0 * tu / du, 14.0 * tv / dv);
                        self.structure.set_tolerance(tu / du, tv / dv);
                    }
                    for w in &wires_uv {
                        self.classifier
                            .register_wire(w, sp.tolerance_uv(), sp.range_u(), sp.range_v());
                    }
                } else {
                    range_invalid = true;
                }
                if !range_invalid {
                    self.splitter = Some(sp);
                }
            }
        }
        // `BRepMesh_NodeInsertionMeshAlgo::initDataStructure`
        // (`BRepMesh_NodeInsertionMeshAlgo.hxx:79-83`): when `AdjustRange` leaves
        // the range splitter invalid, OCCT marks the face failed and returns
        // false, so `BRepMesh_BaseMeshAlgo::process` (`BRepMesh_BaseMeshAlgo.cxx:52-59`)
        // skips `generateMesh`/`commitSurfaceTriangulation` for it. That is what
        // keeps `BRepMesh_GeomTool::CellsCount` (`BRepMesh_GeomTool.cxx:465-511`)
        // from ever seeing a degenerate discrete range: the cylinder branch
        // divides by the V extent (`BRepMesh_GeomTool.cxx:495-501`), and an
        // invalid splitter implies that extent is zero (a zero-extent range
        // samples to a zero-length `computeLengthV`, hence `IsValid` == false).
        // Without this guard a zero V extent reaches `initCirclesTool`
        // (`BRepMesh_Delaun.cxx:280-304`) as a huge cell count and a microscopic
        // cell size, and binding the circumcircles grows the cell grid without
        // bound.
        if range_invalid {
            model
                .face_mut(face_index)
                .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                .set_status(MeshStatus::FAILURE);
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index} has an invalid discrete range"
            ));
        }
        if self.classifier.wires_nb() == 0 {
            let (umin, umax, vmin, vmax) = uv_bounds(&uv);
            for w in &wires_uv {
                self.classifier
                    .register_wire(w, (1e-9, 1e-9), (umin, umax), (vmin, vmax));
            }
        }

        self.init_data_structure(model, face_index)?;

        // Internal (in-face) 3D points -> UV via the face surface, dropped when
        // outside the face (`insertInternalVertex` classifier check).
        let mut internal_uv: Vec<GpPnt2d> = Vec::new();
        let face = model.face(face_index)?;
        if let Some(surface) = face.surface() {
            for p in face.points() {
                if let Some(proj) = project_point_on_surface(surface.as_ref(), p, 1e-7) {
                    let uv_p = GpPnt2d::new(proj.u, proj.v);
                    if self.classifier.is_inside(&uv_p) {
                        internal_uv.push(uv_p);
                    }
                }
            }
        }
        self.insert_internal_nodes(&internal_uv)?;

        let face_deflection = model.face(face_index)?.deflection();
        if self.pre_process_surface_nodes {
            self.generate_surface_nodes(model, face_index, params)?;
            let result = self.finish_mesh(&[], params, face_deflection, face_index);
            if self.delaun_failed {
                model
                    .face_mut(face_index)
                    .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                    .set_status(MeshStatus::FAILURE);
                return Err(format!(
                    "DelaunayNodeInsertionMeshAlgo::perform: face {face_index}: BRepMesh_Delaun::addTriangle failed, no triangulation committed"
                ));
            }
            return result;
        }
        let insert = self.list_surface_nodes(model, face_index, params)?;
        let result = self.finish_mesh(&insert, params, face_deflection, face_index);
        if self.delaun_failed {
            // `BRepMesh_BaseMeshAlgo::Perform` (`BRepMesh_BaseMeshAlgo.cxx:40-62`):
            // the swallowed `Standard_OutOfRange` leaves the face with
            // `IMeshData_Failure` and no triangulation, and the remaining faces of
            // the shape are still meshed.
            model
                .face_mut(face_index)
                .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                .set_status(MeshStatus::FAILURE);
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index}: BRepMesh_Delaun::addTriangle failed, no triangulation committed"
            ));
        }
        result
    }

    /// Collects boundary UV for two OCCT paths that must not be mixed:
    ///
    /// * Delaunay frontier: every pcurve sample `0..=last` with consecutive
    ///   links only (`BRepMesh_BaseMeshAlgo.cxx:99-123`). No last-to-first
    ///   wrap of the concatenated wire — a periodic seam's ends differ in UV
    ///   and are the same 3D vertex; OCCT never adds that diagonal.
    /// * Classifier wire: skip last sample in traversal order
    ///   (`BRepMesh_NodeInsertionMeshAlgo.hxx:152-180` `collectWirePoints`).
    fn collect_boundary_uv(
        model: &MeshModel,
        face_index: usize,
    ) -> Result<(Vec<GpPnt2d>, Vec<GpPnt>, Vec<(i32, i32)>, Vec<Vec<GpPnt2d>>), String> {
        let face = model.face(face_index)?;
        let mut uv: Vec<GpPnt2d> = Vec::new();
        let mut p3d: Vec<GpPnt> = Vec::new();
        let mut constraints: Vec<(i32, i32)> = Vec::new();
        let mut wires_uv: Vec<Vec<GpPnt2d>> = Vec::new();

        for (wire_it, &wire_index) in face.wires().iter().enumerate() {
            let wire = model.wire(wire_index)?;
            // `NodeInsertionMeshAlgo::initDataStructure` (`hxx:70-74`):
            // skip a self-intersecting wire, and skip an OpenWire unless it is
            // the outer wire (`aWireIt == 0`). Those wires are not AddPoint'd
            // and not RegisterWire'd, so they cannot flip registerSurfaceNodes
            // IN to OUT (`BRepMesh_Classifier.cxx:35-59`).
            let skip_classifier = wire.is_status(MeshStatus::SELF_INTERSECTING_WIRE)
                || (wire.is_status(MeshStatus::OPEN_WIRE) && wire_it != 0);
            let mut classifier_wire: Vec<GpPnt2d> = Vec::new();
            for j in 0..wire.edges_nb() {
                let edge_index = wire.edge(j)?;
                let orientation = wire.edge_orientation(j)?;
                let reverse_walk = wire.edge_reverse_walk(j);
                let edge = model.edge(edge_index)?;
                let Some(pcurve) = edge.pcurve_for(face_index, orientation) else {
                    continue;
                };
                let pts = pcurve.points();
                let pts3d = edge.discretization().points();
                let n = pts.len();
                if n == 0 {
                    continue;
                }
                // Frontier links: every sample in the same traversal order as
                // `collectWirePoints`, including the last (`initDataStructure`
                // registers 0..=last). `Ordered < 0` walks the CheckOrder
                // pcurve backwards (same iso, not PCurve2).
                let start = uv.len();
                let reverse = !pcurve.is_forward() ^ reverse_walk;
                if !reverse {
                    for i in 0..n {
                        uv.push(pts[i]);
                        p3d.push(pts3d[i]);
                    }
                } else {
                    for i in (0..n).rev() {
                        uv.push(pts[i]);
                        p3d.push(pts3d[i]);
                    }
                }
                for t in start..uv.len() - 1 {
                    if uv[t] != uv[t + 1] {
                        constraints.push((t as i32, (t + 1) as i32));
                    }
                }

                // `collectWirePoints` (`hxx:153-180`): walk `IsForward` only.
                // `reverse_walk` is ShapeExtend `Ordered<0` / seam PCurve2; cxx
                // does not XOR it into the classifier polygon.
                if pcurve.is_forward() {
                    classifier_wire.extend_from_slice(&pts[..n.saturating_sub(1)]);
                } else {
                    for i in (1..n).rev() {
                        classifier_wire.push(pts[i]);
                    }
                }
            }
            if !skip_classifier && classifier_wire.len() >= 2 {
                wires_uv.push(classifier_wire);
            }
        }

        Ok((uv, p3d, constraints, wires_uv))
    }

    /// Frontier/Fixed links from each pcurve in parameter order.
    /// Source: `BRepMesh_BaseMeshAlgo::initDataStructure` (`BRepMesh_BaseMeshAlgo.cxx:75-128`)
    /// via `NodeInsertionMeshAlgo::initDataStructure` returning `BaseAlgo::initDataStructure()`.
    fn init_data_structure(
        &mut self,
        model: &MeshModel,
        face_index: usize,
    ) -> Result<(), String> {
        let wire_indices: Vec<usize> = model.face(face_index)?.wires().to_vec();
        for (wire_it, &wire_index) in wire_indices.iter().enumerate() {
            let (skip, slots) = {
                let wire = model.wire(wire_index)?;
                // `NodeInsertionMeshAlgo::initDataStructure` (`hxx:70-74`):
                // skip a self-intersecting wire, and skip an OpenWire unless
                // it is the outer wire (`aWireIt == 0`).
                if wire.is_status(MeshStatus::SELF_INTERSECTING_WIRE)
                    || (wire.is_status(MeshStatus::OPEN_WIRE) && wire_it != 0)
                {
                    (true, Vec::new())
                } else {
                    let mut slots = Vec::new();
                    for j in 0..wire.edges_nb() {
                        slots.push((
                            wire.edge(j)?,
                            wire.edge_orientation(j)?,
                            wire.edge_reverse_walk(j),
                        ));
                    }
                    (false, slots)
                }
            };
            if skip {
                continue;
            }
            // One CheckOrder slot → one pcurve. Walk points in the same order
            // as `collectWirePoints` so `addLinkToMesh` left-of-frontier is
            // the classified interior (`meshLeftPolygonOf`).
            for (edge_index, slot_ori, reverse_walk) in slots {
                let job = {
                    let edge = model.edge(edge_index)?;
                    let pts3d = edge.discretization().points().to_vec();
                    let mut found: Option<(Orientation, Vec<GpPnt2d>, Vec<GpPnt>)> = None;
                    for pc_index in edge.pcurves_for(face_index) {
                        let pc = edge.pcurve(pc_index)?;
                        if pc.orientation() != slot_ori {
                            continue;
                        }
                        let mut ori = Self::fix_seam_edge_orientation(&edge, pc_index, face_index)?;
                        if reverse_walk && ori != Orientation::Internal {
                            ori = ori.reversed();
                        }
                        found = Some((ori, pc.points().to_vec(), pts3d.clone()));
                        break;
                    }
                    found
                };
                let Some((ori, pts, pts3d)) = job else {
                    continue;
                };
                let mut prev = -1i32;
                for (i, &uv) in pts.iter().enumerate() {
                    let p3d = pts3d.get(i).copied().unwrap_or_else(GpPnt::zero);
                    let node = self.register_node(p3d, uv, VertexState::Frontier, false);
                    if prev != -1 && prev != node {
                        let links_nb = self.structure.nb_links() as i32;
                        let link_index = self.add_link_to_mesh(prev, node, ori);
                        if wire_it != 0 && link_index <= links_nb {
                            self.structure
                                .set_link_movability(link_index, VertexState::Fixed);
                        }
                    }
                    prev = node;
                }
            }
        }
        Ok(())
    }

    /// Source: `BRepMesh_BaseMeshAlgo::addLinkToMesh`.
    fn add_link_to_mesh(&mut self, first: i32, last: i32, orientation: Orientation) -> i32 {
        let id = if orientation == Orientation::Reversed {
            self.structure
                .add_link(last, first, VertexState::Frontier)
        } else if orientation == Orientation::Internal {
            self.structure.add_link(first, last, VertexState::Fixed)
        } else {
            self.structure
                .add_link(first, last, VertexState::Frontier)
        };
        id.abs()
    }

    /// If another pcurve of the same edge on this face coincides in UV, the
    /// seam is `INTERNAL` (`Fixed`). Else the pcurve's own orientation.
    /// Source: `BRepMesh_BaseMeshAlgo::fixSeamEdgeOrientation`.
    fn fix_seam_edge_orientation(
        edge: &MeshEdge,
        pcurve_index: usize,
        face_index: usize,
    ) -> Result<Orientation, String> {
        let the_pc = edge.pcurve(pcurve_index)?;
        let pts = the_pc.points();
        if pts.is_empty() {
            return Ok(the_pc.orientation());
        }
        let pnt1_1 = pts[0];
        let pnt2_1 = pts[pts.len() - 1];
        for other_index in edge.pcurves_for(face_index) {
            if other_index == pcurve_index {
                continue;
            }
            let other = edge.pcurve(other_index)?;
            let other_pts = other.points();
            if other_pts.is_empty() {
                continue;
            }
            let pnt1_2 = other_pts[0];
            let pnt2_2 = other_pts[other_pts.len() - 1];
            let sq1 = pnt1_1
                .square_distance(&pnt1_2)
                .min(pnt1_1.square_distance(&pnt2_2));
            let sq2 = pnt2_1
                .square_distance(&pnt1_2)
                .min(pnt2_1.square_distance(&pnt2_2));
            if sq1 < SQUARE_CONFUSION && sq2 < SQUARE_CONFUSION {
                return Ok(Orientation::Internal);
            }
        }
        Ok(the_pc.orientation())
    }
}

/// Delaun in face-basis coords, exposing parametric UV to deflection control
/// (`getNodePoint2d` = `Scale(coord, false)`).
struct ParametricDelaun<'a> {
    delaun: &'a mut Delaun,
    ru: (f64, f64),
    rv: (f64, f64),
    delta: (f64, f64),
}

impl ParametricDelaun<'_> {
    fn to_param(&self, p: GpPnt2d) -> GpPnt2d {
        GpPnt2d::new(
            p.x() * self.delta.0 + self.ru.0,
            p.y() * self.delta.1 + self.rv.0,
        )
    }

    fn to_face(&self, p: GpPnt2d) -> GpPnt2d {
        if self.delta.0.abs() <= 1e-16 || self.delta.1.abs() <= 1e-16 {
            return p;
        }
        GpPnt2d::new(
            (p.x() - self.ru.0) / self.delta.0,
            (p.y() - self.rv.0) / self.delta.1,
        )
    }

    /// `DelaunayNodeInsertionMeshAlgo::insertNodes` (`hxx:104-131`).
    fn insert_parametric(
        &mut self,
        algo: &mut DelaunayNodeInsertionMeshAlgo,
        surface: &dyn Surface,
        nodes: &[GpPnt2d],
    ) -> bool {
        // `hxx:108-111`: the early return tests the INPUT list, not the
        // classifier's output.
        if nodes.is_empty() {
            return false;
        }
        let mut idxs = Vec::new();
        for &uv in nodes {
            if !algo.classifier.is_inside(&uv) {
                continue;
            }
            let p3d = surface.d0(uv.x(), uv.y());
            let vertex = DelaunVertex::new(self.to_face(uv), p3d, 0, VertexState::Free);
            let idx = self.delaun.add_node(vertex);
            if idx as usize > algo.nodes_map.len() {
                algo.nodes_map.push(p3d);
            }
            idxs.push(idx);
        }
        // `hxx:125`: `AddVertices` is called even when every node was
        // classified out, so `ProcessConstraints()` still runs.
        self.delaun.add_vertices(&mut idxs);
        // `hxx:130`: only the RETURN VALUE reflects the classifier's output.
        !idxs.is_empty()
    }
}

impl DeflectionMesh for ParametricDelaun<'_> {
    fn elements_of_domain(&self) -> Vec<i32> {
        self.delaun.elements_of_domain()
    }
    fn triangle(&self, id: i32) -> DelaunTriangle {
        self.delaun.triangle(id)
    }
    fn link(&self, id: i32) -> DelaunLink {
        self.delaun.link(id)
    }
    fn link_movability(&self, id: i32) -> VertexState {
        self.delaun.link_movability(id)
    }
    fn vertex(&self, id: i32) -> DelaunVertex {
        let mut v = self.delaun.get_vertex(id);
        v.location = self.to_param(v.location);
        v
    }
    fn shot_triangles(&self, p: GpXY) -> Option<Vec<i32>> {
        // `rejectByMinSize` (`hxx:425-427`): `CircleTool::Select(Scale(p, true))`.
        let face = self.to_face(GpPnt2d::from_xy(p));
        Some(self.delaun.circles().select(face.coord))
    }
}

/// Axis-aligned UV bounds of a point set.
fn uv_bounds(points: &[GpPnt2d]) -> (f64, f64, f64, f64) {
    let (mut umin, mut umax, mut vmin, mut vmax) =
        (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for p in points {
        umin = umin.min(p.x());
        umax = umax.max(p.x());
        vmin = vmin.min(p.y());
        vmax = vmax.max(p.y());
    }
    if !umin.is_finite() {
        return (0.0, 1.0, 0.0, 1.0);
    }
    (umin, umax, vmin, vmax)
}

/// `BRepMesh_GeomTool::CellsCount` + `ComputeErrFactors` + `AdjustCellsCounts`
/// (`BRepMesh_GeomTool.cxx:32-84, 86-144, 465-511`).
fn geom_tool_cells_count(
    surface: &dyn Surface,
    vertices_nb: i32,
    deflection: f64,
    range_u: (f64, f64),
    range_v: (f64, f64),
    delta: (f64, f64),
) -> (i32, i32) {
    let kind = classify_surface(surface);
    let (err_u, err_v) = compute_err_factors(surface, kind, deflection);
    let du = range_u.1 - range_u.0;
    let dv = range_v.1 - range_v.0;
    let mut cells_u;
    let mut cells_v;
    if kind == SurfaceType::Torus {
        cells_u = 2.0_f64.powf((du / delta.0).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / delta.1).log10()).ceil() as i32;
    } else if kind == SurfaceType::Cylinder {
        cells_u = 2.0_f64.powf((du / delta.0 / dv).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / err_v).log10()).ceil() as i32;
    } else {
        cells_u = 2.0_f64.powf((du / delta.0 / err_u).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / delta.1 / err_v).log10()).ceil() as i32;
    }
    if kind == SurfaceType::OtherSurface {
        return (-1, -1);
    }
    let logn = 2.0_f64.powf((vertices_nb as f64).log10()).ceil() as i32;
    match kind {
        SurfaceType::Plane => {
            cells_u = logn;
            cells_v = logn;
        }
        SurfaceType::Cylinder | SurfaceType::Cone => {
            cells_v = logn;
        }
        SurfaceType::SurfaceOfExtrusion | SurfaceType::SurfaceOfRevolution => {
            if let Some(c) = basis_curve_for_err(surface, kind) {
                let planar_like = c.is_line()
                    || (c.bspline_poles().is_some()
                        && c.nurbs_degree().map(|d| d < 2).unwrap_or(false));
                if planar_like {
                    if kind == SurfaceType::SurfaceOfExtrusion {
                        cells_u = logn;
                    } else {
                        cells_v = logn;
                    }
                }
            }
            if kind == SurfaceType::SurfaceOfExtrusion {
                cells_v = logn;
            }
        }
        SurfaceType::BezierSurface | SurfaceType::BSplineSurface => {
            if surface.u_degree() < 2 {
                cells_u = logn;
            }
            if surface.v_degree() < 2 {
                cells_v = logn;
            }
        }
        _ => {}
    }
    (cells_u.max(2), cells_v.max(2))
}

/// `ComputeErrFactors` (`BRepMesh_GeomTool.cxx:32-84`).
fn compute_err_factors(surface: &dyn Surface, kind: SurfaceType, deflection: f64) -> (f64, f64) {
    let mut err_u = deflection * 10.0;
    let mut err_v = deflection * 10.0;
    match kind {
        SurfaceType::Cylinder | SurfaceType::Cone | SurfaceType::Sphere | SurfaceType::Torus => {}
        SurfaceType::SurfaceOfExtrusion | SurfaceType::SurfaceOfRevolution => {
            if let Some(c) = basis_curve_for_err(surface, kind) {
                if c.bspline_poles().is_some() {
                    if let Some(deg) = c.nurbs_degree() {
                        if deg > 2 {
                            let nk = (c.nb_intervals(0) + 1) as f64;
                            err_v /= deg as f64 * nk;
                        }
                    }
                }
            }
        }
        SurfaceType::BezierSurface => {
            let du = surface.u_degree();
            let dv = surface.v_degree();
            if du > 2 {
                err_u /= du as f64;
            }
            if dv > 2 {
                err_v /= dv as f64;
            }
        }
        SurfaceType::BSplineSurface => {
            let du = surface.u_degree();
            let dv = surface.v_degree();
            if du > 2 {
                err_u /= du as f64 * (surface.nb_u_intervals(0) + 1) as f64;
            }
            if dv > 2 {
                err_v /= dv as f64 * (surface.nb_v_intervals(0) + 1) as f64;
            }
        }
        SurfaceType::Plane | SurfaceType::OffsetSurface | SurfaceType::OtherSurface => {
            err_u = 1.0;
            err_v = 1.0;
        }
    }
    (err_u, err_v)
}

fn basis_curve_for_err(
    surface: &dyn Surface,
    kind: SurfaceType,
) -> Option<std::sync::Arc<dyn occt_geom::Curve>> {
    match kind {
        SurfaceType::SurfaceOfExtrusion => surface.extrusion_basis_curve(),
        SurfaceType::SurfaceOfRevolution => surface.revolution_basis_curve(),
        _ => None,
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
        let vertex = DelaunVertex::new(self.scale_to_face(point), p3d, 0, movability);
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
            self.register_node(GpPnt::zero(), p, VertexState::Free, false);
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
        // OCCT BRepMesh_Delaun.cxx:703 calls ProcessConstraints() unconditionally;
        // frontierAdjust() ends with cleanupMesh() (cxx:1028) which prunes boundary
        // triangles whose neighbour touches the super-triangle. These nodes carry only
        // Free links, so the mesh keeps 6 triangles over 7 nodes instead of 8 over 9.
        assert_eq!(res.triangles.len(), 6);
        assert_eq!(res.nodes.len(), 7);
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
