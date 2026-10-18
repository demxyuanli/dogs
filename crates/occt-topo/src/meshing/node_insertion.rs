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
    DelaunayDeflectionControlMeshAlgo, DeflectionMesh,
};
use super::parameters::MeshParameters;
use super::range_splitter::{
    classify_surface, create_range_splitter, factory_uses_deflection_control,
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
        let cls_dbg = std::env::var("CLSNODE").is_ok();
        let (mut n_in, mut n_out, mut n_on) = (0usize, 0usize, 0usize);
        for n in nodes {
            let st = self.classifier.perform(&n);
            match st {
                PointState::In => n_in += 1,
                PointState::Out => n_out += 1,
                PointState::On => n_on += 1,
            }
            if cls_dbg && st != PointState::In {
                eprintln!(
                    "CLSNODE face={face_index} st={:?} uv=({:.12e},{:.12e})",
                    st,
                    n.x(),
                    n.y()
                );
            }
            if st == PointState::In {
                out.push((n, surface.d0(n.x(), n.y())));
            }
        }
        if cls_dbg {
            eprintln!(
                "CLSNODE face={face_index} wires={} total={} in={n_in} out={n_out} on={n_on}",
                self.classifier.wires_nb(),
                n_in + n_out + n_on
            );
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
        _face_index: usize,
    ) -> Result<TriangulationResult, String> {
        if self.structure.nb_nodes() == 0 {
            return Err("DelaunayNodeInsertionMeshAlgo::triangulate: no nodes registered".to_string());
        }

        let structure = std::mem::replace(&mut self.structure, DelaunDataStructure::new(0));
        let mut indices: Vec<i32> = (1..=structure.nb_nodes() as i32)
            .filter(|&i| structure.get_node(i).state != VertexState::Deleted)
            .collect();
        let (cells_u, cells_v) = self.cells_count(indices.len() as i32);
        // TEMPORARY diagnostic (DELAUNDUMP): writes the exact pre-triangulation
        // structure (nodes, links, insertion order) of the `DELSTAGE` face so
        // `occt_probe --delaunstruct` can run OCCT's own `BRepMesh_Delaun` over
        // the identical input.
        if let Ok(dir) = std::env::var("DELAUNDUMP") {
            let tagged = std::env::var("DELSTAGE")
                .ok()
                .and_then(|v| v.parse::<i32>().ok())
                == Some(_face_index as i32);
            if tagged {
                let free_xy: Vec<(f64, f64)> = insert
                    .iter()
                    .map(|(uv, _p3d)| {
                        let c = self.scale_to_face(*uv);
                        (c.x(), c.y())
                    })
                    .collect();
                Self::dump_delaun_structure(&structure, &dir, _face_index, cells_u, cells_v, &indices, &free_xy);
            }
        }
        let mut delaun = Delaun::new_with_data_cells_diag(structure, &mut indices, cells_u, cells_v, _face_index as i32);
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
        let _bnd_nodes = delaun.result().nb_nodes();
        let _bnd_tris = delaun.result().elements_of_domain().len();

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
            // TEMPORARY diagnostic (DELAUNDUMP_GRID): two-phase replay dump of
            // the state right before `AddVertices` (base mesh + registered
            // surface nodes), so `occt_probe --delaunstruct` can replay OCCT's
            // own `BRepMesh_Delaun` over the identical structure and order.
            let grid_dump_dir = std::env::var("DELAUNDUMP_GRID").ok().filter(|_| {
                std::env::var("DELSTAGE").ok().and_then(|v| v.parse::<i32>().ok())
                    == Some(_face_index as i32)
            });
            if let Some(dir) = grid_dump_dir.as_deref() {
                Self::dump_delaun_grid(
                    delaun.result(),
                    dir,
                    _face_index,
                    cells_u,
                    cells_v,
                    &indices,
                    &idxs,
                );
            }
            if let Some(dir) = grid_dump_dir.as_deref() {
                Self::dump_domain_triangles(delaun.result(), dir, _face_index, "base");
            }
            // TEMPORARY diagnostic (ADDDUM): shifts the element ids of the
            // insertions below by `ADDDUM` dead slots, to test whether OCCT's
            // extra base-mesh element creation is what aligns the two meshes.
            if let Ok(n) = std::env::var("ADDDUM") {
                if let Ok(n) = n.parse::<usize>() {
                    delaun.debug_bump_elements(n);
                }
            }
            // TEMPORARY diagnostic (DELONESTEP): instead of one batched
            // `AddVertices`, feed the surface nodes one at a time (the OCCT
            // `AddVertices(single)` call) and dump the domain triangle set after
            // every step, so `occt_probe --delaunsteps` can be diffed step by
            // step and the first divergent insertion identified.
            //
            // The batched path heap sorts first (`AddVertices`, `cxx:2415-2417`),
            // so the step dump applies the same sort: otherwise the steps walk an
            // order neither OCCT run uses. `DELONESTEP_RAW` keeps the raw
            // registration order.
            if let Ok(dir) = std::env::var("DELONESTEP") {
                let mut order: Vec<i32> = idxs.clone();
                if std::env::var("DELONESTEP_RAW").is_err() {
                    let ds = delaun.result();
                    crate::meshing::delaun::heap_sort_by_xy_sum(&mut order, |i| {
                        let p = ds.get_node(*i).location;
                        (p.x(), p.y())
                    });
                }
                let mut out = String::new();
                for (k, &idx) in order.iter().enumerate() {
                    let mut one = vec![idx];
                    delaun.add_vertices(&mut one);
                    Self::append_step_triangles(&mut out, k, delaun.result());
                }
                let path = std::path::Path::new(&dir)
                    .join(format!("delaun_steps_port_f{_face_index}.txt"));
                if std::fs::write(&path, out).is_err() {
                    eprintln!("DELONESTEP: cannot write {}", path.display());
                }
            } else {
                delaun.add_vertices(&mut idxs);
            }
            if let Some(dir) = grid_dump_dir.as_deref() {
                Self::dump_domain_triangles(delaun.result(), dir, _face_index, "grid");
            }
        }

        // `DelaunayDeflectionControlMeshAlgo::postProcessMesh` (`hxx:50-68`):
        // after NodeInsertion surface nodes, `optimizeMesh` when the factory
        // selected DeflectionControl and `ControlSurfaceDeflection` is on.
        if params.control_surface_deflection {
            if let Some(surf) = self.splitter.as_ref().and_then(|s| s.surface().cloned()) {
                if factory_uses_deflection_control(surf.as_ref(), params)
                    && !delaun.result().elements_of_domain().is_empty()
                {
                    if std::env::var("DCDBG").is_ok() {
                        let ds = delaun.result();
                        let mut edge_fp: u64 = 0;
                        for &tid in ds.elements_of_domain() {
                            let e = ds.get_element(tid);
                            let n = ds.element_nodes(&e);
                            for k in 0..3 {
                                let a = n[k].min(n[(k + 1) % 3]) as u64;
                                let b = n[k].max(n[(k + 1) % 3]) as u64;
                                edge_fp = edge_fp
                                    .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                                    .wrapping_add(a.wrapping_mul(0x10003) ^ b);
                            }
                        }
                        eprintln!(
                            "DCDBG pre_dc f={} nodes={} tris={} face_defl={:.6} cells=({},{}) delta=({:.6},{:.6}) insert={} edge_fp={edge_fp:016x}",
                            _face_index,
                            ds.nb_nodes(),
                            ds.elements_of_domain().len(),
                            face_deflection,
                            cells_u,
                            cells_v,
                            self.splitter.as_ref().map(|s| s.delta().0).unwrap_or(0.0),
                            self.splitter.as_ref().map(|s| s.delta().1).unwrap_or(0.0),
                            insert.len(),
                        );
                    }
                    // TEMP STRUCTDUMP: full pre-DC node/element dump, so it can
                    // be diffed entry by entry against the OCCT probe's
                    // `OCCTDC PRE` / `OCCTDCN PRE` / `OCCTEL PRE` lines.
                    if std::env::var("PORTDUMP").ok().filter(|_| {
                        std::env::var("DELSTAGE").ok().and_then(|v| v.parse::<i32>().ok())
                            == Some(_face_index as i32)
                    }).is_some() {
                        let ds = delaun.result();
                        eprintln!(
                            "PORTDC PRE f={} struct_nodes={} nodes={} domain={}",
                            _face_index,
                            ds.nb_nodes(),
                            ds.nb_nodes(),
                            ds.elements_of_domain().len()
                        );
                        for tid in ds.elements_of_domain() {
                            let e = ds.get_element(*tid);
                            let n = ds.element_nodes(&e);
                            eprintln!("PORTT PRE T {} {} {} {}", tid, n[0], n[1], n[2]);
                        }
                        for i in 1..=ds.nb_nodes() as i32 {
                            let v = ds.get_node(i);
                            eprintln!(
                                "PORTN PRE {} {:.17e} {:.17e} {}",
                                i,
                                v.location.x(),
                                v.location.y(),
                                v.state.index()
                            );
                        }
                        for i in 1..=ds.nb_elements() as i32 {
                            let e = ds.get_element(i);
                            let n = ds.element_nodes(&e);
                            eprintln!("PORTEL PRE E {} {} {} {} {}", i, n[0], n[1], n[2], 0 as i32);
                        }
                    }
                    self.optimize_mesh(&mut delaun, surf.as_ref(), params, face_deflection, _face_index);
                    if std::env::var("DCDBG").is_ok() {
                        eprintln!(
                            "DCDBG post_dc nodes={} tris={}",
                            delaun.result().nb_nodes(),
                            delaun.result().elements_of_domain().len()
                        );
                    }
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
        face_index: usize,
    ) {
        let (ru, rv, delta) = match self.splitter.as_ref() {
            Some(sp) => (sp.range_u(), sp.range_v(), sp.delta()),
            None => return,
        };
        let mut ctrl = DelaunayDeflectionControlMeshAlgo::new(params);
        let mut view = ParametricDelaun { delaun, ru, rv, delta };
        // TEMPORARY diagnostic (DELAUNDUMP_GRID): per-pass domain triangles, to
        // line up each `optimizeMesh` pass against the probe's `OCCTDCPASS`.
        let dump_dir = std::env::var("DELAUNDUMP_GRID").ok().filter(|_| {
            std::env::var("DELSTAGE").ok().and_then(|v| v.parse::<i32>().ok())
                == Some(face_index as i32)
        });
        let mut pass_no = 0usize;
        let ctrl_dump = std::env::var("CTRLDUMP").ok().filter(|_| {
            std::env::var("DELSTAGE").ok().and_then(|v| v.parse::<i32>().ok())
                == Some(face_index as i32)
        });
        let max_defl = ctrl.post_process(&mut view, surface, face_deflection, &mut |mesh, nodes| {
            pass_no += 1;
            let pre_nodes = mesh.delaun.result().nb_nodes();
            let pre_dom = mesh.delaun.result().elements_of_domain().len();
            let inserted = mesh.insert_parametric(self, surface, nodes);
            if std::env::var("DCDBG").is_ok() {
                eprintln!(
                    "DCDBG insert f={face_index} pass={pass_no} ctrl={} pre_nodes={} pre_dom={} post_nodes={} post_dom={} inserted={}",
                    nodes.len(),
                    pre_nodes,
                    pre_dom,
                    mesh.delaun.result().nb_nodes(),
                    mesh.delaun.result().elements_of_domain().len(),
                    inserted
                );
            }
            // TEMPORARY diagnostic (CTRLDUMP): the control node list this pass
            // hands to `insertNodes`, in face-basis coordinates.
            if let Some(dir) = ctrl_dump.as_deref() {
                use std::fmt::Write as _;
                let mut out = String::new();
                for n in nodes {
                    let _ = writeln!(out, "{:.17e} {:.17e}", n.x(), n.y());
                }
                let path = std::path::Path::new(dir)
                    .join(format!("ctrl_pass{pass_no}_f{face_index}.txt"));
                if std::fs::write(&path, out).is_err() {
                    eprintln!("CTRLDUMP: cannot write {}", path.display());
                }
            }
            if let Some(dir) = dump_dir.as_deref() {
                Self::dump_domain_state(
                    mesh.delaun.result(),
                    dir,
                    face_index,
                    &format!("pass{pass_no}"),
                );
            }
            inserted
        });
        // `hxx:117-120`: `if (!(myMaxSqDeflection < 0.)) SetDeflection(sqrt(...))`.
        // `post_process` returns 0 when the cache stayed at -1, so use the
        // returned value only when a sample actually updated it (`> 0` or a
        // zero-deflection flat pass that still ran geometry — both yield a
        // non-negative cache; flat geometry returns 0 from sqrt(0)).
        if ctrl.max_sq_deflection >= 0.0 {
            self.optimized_face_deflection = Some(max_defl);
        }
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
        if std::env::var("MBDIAG").is_ok() && uv.len() <= 1 {
            eprintln!(
                "MBDIAG thin face={face_index} uv={} wires={} edges={:?}",
                uv.len(),
                wires_uv.len(),
                wires_uv.iter().map(|w| w.len()).collect::<Vec<_>>()
            );
        }
        if uv.is_empty() {
            return Err(format!(
                "DelaunayNodeInsertionMeshAlgo::perform: face {face_index} has no boundary UV points"
            ));
        }
        self.boundary_uv = uv.clone();

        // TEMP diagnostic (TORSBD=<x1>,<y1>,<z1>,<x2>,<y2>,<z2>): per-wire-edge
        // pcurve point arrays of the matching face, so they can be diffed against
        // `occt_probe --healerdbg` (the same six bbox coordinates).
        if let Ok(spec) = std::env::var("TORSBD") {
            let t: Vec<f64> = spec
                .split(',')
                .filter_map(|s| s.trim().parse::<f64>().ok())
                .collect();
            if t.len() == 6 {
                let bb = crate::brep_bnd_lib::shape_bnd_box(&model.face(face_index)?.face().0);
                let (mn, mx) = (bb.corner_min(), bb.corner_max());
                {
                    let hit = (mn.x() - t[0]).abs() < 1e-6
                        && (mn.y() - t[1]).abs() < 1e-6
                        && (mn.z() - t[2]).abs() < 1e-6
                        && (mx.x() - t[3]).abs() < 1e-6
                        && (mx.y() - t[4]).abs() < 1e-6
                        && (mx.z() - t[5]).abs() < 1e-6;
                    if hit {                        let face = model.face(face_index)?;
                        for (wi, &wire_index) in face.wires().iter().enumerate() {
                            let wire = model.wire(wire_index)?;
                            for j in 0..wire.edges_nb() {
                                let edge_index = wire.edge(j)?;
                                let orientation = wire.edge_orientation(j)?;
                                let rev = wire.edge_reverse_walk(j);
                                let edge = model.edge(edge_index)?;
                                let Some(pc) = edge.pcurve_for(face_index, orientation) else {
                                    eprintln!("TORSBD face={face_index} w={wi} e={j} pc=none");
                                    continue;
                                };
                                let pts = pc.points();
                                let n = pts.len();
                                let f = pts.first().copied().unwrap_or(GpPnt2d::new(0.0, 0.0));
                                let l = pts.last().copied().unwrap_or(GpPnt2d::new(0.0, 0.0));
                                eprintln!(
                                    "TORSBD face={face_index} w={wi} e={j} nb={n} fwd={} rev={rev} \
same_param={} same_range={} free={} degen={} defl={:.15} ang={:.15} \
first=({:.15},{:.15}) last=({:.15},{:.15})",
                                    pc.is_forward() as u8,
                                    edge.same_param() as u8,
                                    edge.same_range() as u8,
                                    edge.is_free() as u8,
                                    edge.degenerated() as u8,
                                    edge.deflection(),
                                    edge.angular_deflection(),
                                    f.x(),
                                    f.y(),
                                    l.x(),
                                    l.y()
                                );
                                let mut s = format!("TORSBDV face={face_index} w={wi} e={j} v=");
                                for p in pts {
                                    s += &format!("{:.15},", p.y());
                                }
                                eprintln!("{s}");
                            }
                        }
                    }
                }
            }
        }

        // `BRepMesh_NodeInsertionMeshAlgo::initDataStructure` (`hxx:47-112`):
        // Reset the splitter, feed the discrete wire points (`collectWirePoints`
        // calls `myRangeSplitter.AddPoint` per pcurve sample, `hxx:179`),
        // AdjustRange, SetCellSize/SetTolerance, then RegisterWire with the
        // splitter range. Without AddPoint a plane's range stays the adaptor's
        // infinite bounds and the cell size collapses the whole face.
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
                    if std::env::var("MBDIAG").is_ok() {
                        let (su0, su1) = surface.as_ref().u_range();
                        let (sv0, sv1) = surface.as_ref().v_range();
                        let (ru0, ru1) = sp.range_u();
                        let (rv0, rv1) = sp.range_v();
                        let (du, dv) = sp.delta();
                        let (bu0, bu1, bv0, bv1) = uv_bounds(&uv);
                        eprintln!(
                            "MBDIAG rng face={face_index} surf={:?} sup={} svp={} su=[{su0:.6},{su1:.6}] sv=[{sv0:.6},{sv1:.6}] du=[{ru0:.6},{ru1:.6}] dv=[{rv0:.6},{rv1:.6}] delta=({du:.6},{dv:.6}) buv=[{bu0:.6},{bu1:.6}]x[{bv0:.6},{bv1:.6}] nuv={}",
                            classify_surface(surface.as_ref()),
                            surface.as_ref().is_u_periodic(),
                            surface.as_ref().is_v_periodic(),
                            uv.len()
                        );
                    }
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
        self.diag_dump_structure(face_index);
        if let Ok(spec) = std::env::var("DPROBE") {
            if spec.trim().parse::<i32>().ok() == Some(face_index as i32) {
                let nb = self.structure.nb_nodes() as i32;
                let (mut nfree, mut nfront, mut nfixed) = (0usize, 0usize, 0usize);
                for i in 1..=nb {
                    match self.structure.get_node(i).state {
                        VertexState::Free => nfree += 1,
                        VertexState::Frontier => nfront += 1,
                        VertexState::Fixed => nfixed += 1,
                        _ => {}
                    }
                }
                eprintln!(
                    "DPROBE face={face_index} nb_nodes={nb} nb_links={} uv={} free={nfree} frontier={nfront} fixed={nfixed}",
                    self.structure.nb_links(),
                    uv.len()
                );
                for i in 1..=nb {
                    let n = self.structure.get_node(i);
                    let p3d = self
                        .nodes_map
                        .get(i as usize - 1)
                        .copied()
                        .unwrap_or_else(GpPnt::zero);
                    eprintln!(
                        "DPROBE n i={i} st={} uv=({:.6},{:.6}) p=({:.6},{:.6},{:.6})",
                        n.state.to_str(),
                        n.location.x(),
                        n.location.y(),
                        p3d.x(),
                        p3d.y(),
                        p3d.z()
                    );
                }
            }
        }

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
        self.optimized_face_deflection = None;
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
            if let Some(d) = self.optimized_face_deflection.take() {
                model
                    .face_mut(face_index)
                    .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                    .set_deflection(d);
            }
            return result;
        }
        let insert = self.list_surface_nodes(model, face_index, params)?;
        let result = self.finish_mesh(&insert, params, face_deflection, face_index);
        if std::env::var("DPROBE").ok().and_then(|v| v.trim().parse::<i32>().ok())
            == Some(face_index as i32)
        {
            if let Ok(res) = &result {
                eprintln!(
                    "DPROBE result face={face_index} nodes={} tris={} insert={}",
                    res.nodes.len(),
                    res.triangles.len(),
                    insert.len()
                );
                for (i, n) in res.nodes.iter().enumerate() {
                    eprintln!(
                        "DPROBE r i={} st={} uv=({:.6},{:.6}) p=({:.6},{:.6},{:.6})",
                        i + 1,
                        n.state.to_str(),
                        n.location.x(),
                        n.location.y(),
                        n.p3d.x(),
                        n.p3d.y(),
                        n.p3d.z()
                    );
                }
            }
        }
        if let Some(d) = self.optimized_face_deflection.take() {
            model
                .face_mut(face_index)
                .map_err(|e| format!("DelaunayNodeInsertionMeshAlgo::perform: {e}"))?
                .set_deflection(d);
        }
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
                if std::env::var("BPUV").is_ok() {
                    let (mut minv, mut maxv) = (f64::INFINITY, f64::NEG_INFINITY);
                    let (mut minu, mut maxu) = (f64::INFINITY, f64::NEG_INFINITY);
                    for p in pts {
                        minv = minv.min(p.y());
                        maxv = maxv.max(p.y());
                        minu = minu.min(p.x());
                        maxu = maxu.max(p.x());
                    }
                    eprintln!(
                        "BPUV face={face_index} wire={wire_it} edge={edge_index} n={n} \
first=({:.17e},{:.17e}) last=({:.17e},{:.17e}) \
u=[{:.17e},{:.17e}] v=[{:.17e},{:.17e}]",
                        pts[0].x(),
                        pts[0].y(),
                        pts[n - 1].x(),
                        pts[n - 1].y(),
                        minu,
                        maxu,
                        minv,
                        maxv
                    );
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

                // `collectWirePoints` (`hxx:153-180`): walk the pcurve in the
                // effective edge orientation. `reverse_walk` models
                // `ShapeExtend_WireData::Edge(signed)` = `Edge.Reverse()` on a
                // seam; `BRepMesh_ShapeVisitor::addWire` stores that reversed
                // orientation into `IMeshData_Wire`, so `GetEdgeOrientation`
                // (and thus `IMeshData_PCurve::IsForward`) already includes it.
                // The port keeps the original orientation plus `reverse_walk`,
                // so the classifier polygon must XOR it exactly like the
                // frontier links above.
                if !reverse {
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
            let (skip, slot_edges) = {
                let wire = model.wire(wire_index)?;
                // `NodeInsertionMeshAlgo::initDataStructure` (`hxx:70-74`):
                // skip a self-intersecting wire, and skip an OpenWire unless
                // it is the outer wire (`aWireIt == 0`).
                if wire.is_status(MeshStatus::SELF_INTERSECTING_WIRE)
                    || (wire.is_status(MeshStatus::OPEN_WIRE) && wire_it != 0)
                {
                    (true, Vec::new())
                } else {
                    let mut edges = Vec::new();
                    for j in 0..wire.edges_nb() {
                        edges.push(wire.edge(j)?);
                    }
                    (false, edges)
                }
            };
            if skip {
                continue;
            }
            // `BRepMesh_BaseMeshAlgo::initDataStructure` (`cxx:87-125`): for every
            // wire slot, *every* pcurve of that edge on this face is registered,
            // in `aDEdge->GetPCurves(myDFace)` order (the faces' pcurve index
            // list). A seam edge the wire walks twice therefore contributes both
            // of its period copies, and the walk order is the pcurve's own point
            // order - OCCT never consults the CheckOrder slot orientation here.
            // The link orientation is `fixSeamEdgeOrientation` (`cxx:197-225`).
            for edge_index in slot_edges {
                let jobs: Vec<(usize, Orientation, Vec<GpPnt2d>, Vec<GpPnt>)> = {
                    let edge = model.edge(edge_index)?;
                    let pts3d = edge.discretization().points().to_vec();
                    let mut jobs = Vec::new();
                    for pc_index in edge.pcurves_for(face_index) {
                        let pc = edge.pcurve(pc_index)?;
                        let ori = Self::fix_seam_edge_orientation(&edge, pc_index, face_index)?;
                        jobs.push((pc_index, ori, pc.points().to_vec(), pts3d.clone()));
                    }
                    jobs
                };
                for (pc_index, ori, pts, pts3d) in jobs {
                    // TEMPORARY diagnostic (DELREGSKEL): the registration skeleton in the
                    // `occt_probe --structdump` `OCCTSTRUCT E` line format, so the pcurve
                    // visit order and the effective link orientation can be diffed against
                    // `BRepMesh_BaseMeshAlgo::initDataStructure`.
                    if std::env::var("DELREGSKEL").is_ok() {
                        let first = pts.first().copied().unwrap_or_else(GpPnt2d::zero);
                        let last = pts.last().copied().unwrap_or_else(GpPnt2d::zero);
                        eprintln!(
                            "REG f={} w={} edge={} pc={} ori={:?} n={} first=({:.9},{:.9}) last=({:.9},{:.9})",
                            face_index,
                            wire_it,
                            edge_index,
                            pc_index,
                            ori,
                            pts.len(),
                            first.x(),
                            first.y(),
                            last.x(),
                            last.y()
                        );
                    }
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
        }
        Ok(())
    }

    /// TEMPORARY diagnostic (DELAUNDUMP): structure dump in the text format
    /// consumed by `occt_probe --delaunstruct <dir>` (`delaun_in_f<face>.txt`):
    /// `TOL`, `CELLS`, `N i u v movability`, `L i first last movability`, `V ...`.
    fn dump_delaun_structure(
        structure: &DelaunDataStructure,
        dir: &str,
        face_index: usize,
        cells_u: i32,
        cells_v: i32,
        indices: &[i32],
        free_xy: &[(f64, f64)],
    ) {
        use std::fmt::Write as _;
        let (tol_u, tol_v) = structure.get_tolerance();
        let mut out = String::new();
        let _ = writeln!(out, "TOL {tol_u:.17e} {tol_v:.17e}");
        let _ = writeln!(out, "CELLS {cells_u} {cells_v}");
        // 3D bbox of the registered nodes: pairs this dump with the probe's
        // `--boundary` / `--facestats` face list (whose face order differs from
        // the port's, so the index alone is not a pairing key).
        {
            let (mut lo, mut hi) = (GpPnt::zero(), GpPnt::zero());
            let mut first = true;
            for i in 1..=structure.nb_nodes() as i32 {
                let p = structure.get_node(i).p3d;
                if first {
                    lo = p;
                    hi = p;
                    first = false;
                } else {
                    lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
                    hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
                }
            }
            if !first {
                let _ = writeln!(
                    out,
                    "BOX {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
                    lo.x(),
                    lo.y(),
                    lo.z(),
                    hi.x(),
                    hi.y(),
                    hi.z()
                );
            }
        }
        for i in 1..=structure.nb_nodes() as i32 {
            let n = structure.get_node(i);
            let _ = writeln!(
                out,
                "N {i} {:.17e} {:.17e} {}",
                n.location.x(),
                n.location.y(),
                n.state.index()
            );
        }
        for i in 1..=structure.nb_links() as i32 {
            let l = structure.get_link(i);
            let _ = writeln!(
                out,
                "L {i} {} {} {}",
                l.first_node(),
                l.last_node(),
                structure.link_movability(i).index()
            );
        }
        let mut vline = String::from("V");
        for &idx in indices {
            let _ = write!(vline, " {idx}");
        }
        let _ = writeln!(out, "{vline}");
        // TEMPORARY diagnostic (DELAUNDUMP): the surface nodes that
        // `AddVertices` will insert are not registered in the structure yet, so
        // they are appended as `N` lines with ids continuing after the
        // registered ones, together with the phase-1/phase-2 index lists.
        let base = structure.nb_nodes() as i32;
        for (k, (x, y)) in free_xy.iter().enumerate() {
            let _ = writeln!(out, "N {} {x:.17e} {y:.17e} {}", base + k as i32 + 1, VertexState::Free.index());
        }
        let mut v1 = String::from("V1");
        for &idx in indices {
            let _ = write!(v1, " {idx}");
        }
        let _ = writeln!(out, "{v1}");
        let mut v2 = String::from("V2");
        for k in 0..free_xy.len() {
            let _ = write!(v2, " {}", base + k as i32 + 1);
        }
        let _ = writeln!(out, "{v2}");
        let path = std::path::Path::new(dir).join(format!("delaun_in_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELAUNDUMP_GRID): dumps the pre-`AddVertices`
    /// structure (nodes and links of the base mesh plus the registered surface
    /// nodes) together with the two insertion phases as `V1` / `V2`, so the
    /// probe can replay `BRepMesh_Delaun(structure, V1, cellsU, cellsV)` followed
    /// by `AddVertices(V2)` -- the exact control flow of
    /// `BRepMesh_DelaunayNodeInsertionMeshAlgo::postProcessMesh`.
    fn dump_delaun_grid(
        structure: &DelaunDataStructure,
        dir: &str,
        face_index: usize,
        cells_u: i32,
        cells_v: i32,
        phase1: &[i32],
        phase2: &[i32],
    ) {
        use std::fmt::Write as _;
        let (tol_u, tol_v) = structure.get_tolerance();
        let mut out = String::new();
        let _ = writeln!(out, "TOL {tol_u:.17e} {tol_v:.17e}");
        let _ = writeln!(out, "CELLS {cells_u} {cells_v}");
        for i in 1..=structure.nb_nodes() as i32 {
            let n = structure.get_node(i);
            let _ = writeln!(
                out,
                "N {i} {:.17e} {:.17e} {}",
                n.location.x(),
                n.location.y(),
                n.state.index()
            );
        }
        for i in 1..=structure.nb_links() as i32 {
            let l = structure.get_link(i);
            let _ = writeln!(
                out,
                "L {i} {} {} {}",
                l.first_node(),
                l.last_node(),
                structure.link_movability(i).index()
            );
        }
        for (tag, ids) in [("V1", phase1), ("V2", phase2)] {
            let mut line = String::from(tag);
            for &idx in ids {
                let _ = write!(line, " {idx}");
            }
            let _ = writeln!(out, "{line}");
        }
        let path = std::path::Path::new(dir).join(format!("delaun_in_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP_GRID: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELAUNDUMP_GRID): per-pass state of the deflection
    /// control, node table (`ND <id> <x> <y>`) plus domain triangles
    /// (`T <tid> <a> <b> <c>`), for coordinate-level comparison with the probe's
    /// `OCCTDC P<n>` dump.
    fn dump_domain_state(structure: &DelaunDataStructure, dir: &str, face_index: usize, tag: &str) {
        use std::fmt::Write as _;
        let mut out = String::new();
        for i in 1..=structure.nb_nodes() as i32 {
            let n = structure.get_node(i);
            if n.state == VertexState::Deleted {
                continue;
            }
            let _ = writeln!(out, "ND {i} {:.17e} {:.17e}", n.location.x(), n.location.y());
        }
        let mut tris: Vec<(i32, [i32; 3])> = Vec::new();
        for &tid in structure.elements_of_domain() {
            let e = structure.get_element(tid);
            let mut n = structure.element_nodes(&e);
            n.sort_unstable();
            tris.push((tid, n));
        }
        tris.sort_unstable();
        let _ = writeln!(out, "TRIS {}", tris.len());
        for (tid, n) in &tris {
            let _ = writeln!(out, "T {} {} {} {}", tid, n[0], n[1], n[2]);
        }
        let path = std::path::Path::new(dir).join(format!("delaun_{tag}_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP_GRID: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELONESTEP): appends one insertion step's domain
    /// triangle set (`STEP <k>` header, then `T <a> <b> <c>` with sorted node
    /// ids) to the step file, for step-by-step diffing against the probe's
    /// per-vertex `BRepMesh_Delaun::AddVertices` replay.
    fn append_step_triangles(out: &mut String, step: usize, structure: &DelaunDataStructure) {
        use std::fmt::Write as _;
        let _ = writeln!(out, "STEP {step}");
        let mut tris: Vec<[i32; 3]> = Vec::new();
        for &tid in structure.elements_of_domain() {
            let e = structure.get_element(tid);
            let mut n = structure.element_nodes(&e);
            n.sort_unstable();
            tris.push(n);
        }
        tris.sort_unstable();
        for n in &tris {
            let _ = writeln!(out, "T {} {} {}", n[0], n[1], n[2]);
        }
    }

    /// TEMPORARY diagnostic (DELAUNDUMP_GRID): dumps the domain triangles of the
    /// structure as `T <a> <b> <c>` with sorted node ids, to be diffed against
    /// the probe's replay of the same structure.
    fn dump_domain_triangles(structure: &DelaunDataStructure, dir: &str, face_index: usize, tag: &str) {
        use std::fmt::Write as _;
        let mut out = String::new();
        let mut tris: Vec<(i32, [i32; 3])> = Vec::new();
        for &tid in structure.elements_of_domain() {
            let e = structure.get_element(tid);
            let mut n = structure.element_nodes(&e);
            n.sort_unstable();
            tris.push((tid, n));
        }
        tris.sort_unstable();
        let _ = writeln!(out, "TRIS {}", tris.len());
        for (tid, n) in &tris {
            let _ = writeln!(out, "T {} {} {} {}", tid, n[0], n[1], n[2]);
        }
        let path = std::path::Path::new(dir).join(format!("delaun_{tag}_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP_GRID: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELDIAG): dumps the pre-triangulation structure of
    /// the face whose boundary contains the node at (13.0176, -24.7899, 290).
    fn diag_dump_structure(&self, face_index: usize) {
        if std::env::var("DELDIAG").is_err() {
            return;
        }
        let nb = self.structure.nb_nodes() as i32;
        let is_target = (1..=nb).any(|i| {
            let p = self.structure.get_node(i).p3d;
            (p.x() - 13.017621944).abs() < 1e-3
                && (p.y() + 24.789947988).abs() < 1e-3
                && (p.z() - 290.0).abs() < 1e-3
        });
        if !is_target {
            return;
        }
        let nb_links = self.structure.nb_links() as i32;
        eprintln!("DELDIAG face={face_index} nodes={nb} links={nb_links}");
        for i in 1..=nb {
            let n = self.structure.get_node(i);
            eprintln!(
                "DELDIAG n i={i} st={} uv=({:.6},{:.6}) p=({:.6},{:.6},{:.6})",
                n.state.to_str(),
                n.location.x(),
                n.location.y(),
                n.p3d.x(),
                n.p3d.y(),
                n.p3d.z()
            );
        }
        for i in 1..=nb_links {
            let l = self.structure.get_link(i);
            let a = self.structure.get_node(l.first_node());
            let b = self.structure.get_node(l.last_node());
            eprintln!(
                "DELDIAG l i={i} st={} n=({},{}) p=[({:.6},{:.6},{:.6}),({:.6},{:.6},{:.6})]",
                self.structure.link_movability(i).to_str(),
                l.first_node(),
                l.last_node(),
                a.p3d.x(),
                a.p3d.y(),
                a.p3d.z(),
                b.p3d.x(),
                b.p3d.y(),
                b.p3d.z()
            );
        }
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
        let mut out_cls = 0usize;
        for &uv in nodes {
            if !algo.classifier.is_inside(&uv) {
                out_cls += 1;
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
        if std::env::var("DCDBG").is_ok() {
            eprintln!(
                "DCDBG insert in={} out_cls={} added={}",
                nodes.len(),
                out_cls,
                idxs.len()
            );
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
    if kind == MeshSurfaceType::Torus {
        cells_u = 2.0_f64.powf((du / delta.0).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / delta.1).log10()).ceil() as i32;
    } else if kind == MeshSurfaceType::Cylinder {
        cells_u = 2.0_f64.powf((du / delta.0 / dv).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / err_v).log10()).ceil() as i32;
    } else {
        cells_u = 2.0_f64.powf((du / delta.0 / err_u).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / delta.1 / err_v).log10()).ceil() as i32;
    }
    if kind == MeshSurfaceType::OtherSurface {
        return (-1, -1);
    }
    let logn = 2.0_f64.powf((vertices_nb as f64).log10()).ceil() as i32;
    match kind {
        MeshSurfaceType::Plane => {
            cells_u = logn;
            cells_v = logn;
        }
        MeshSurfaceType::Cylinder | MeshSurfaceType::Cone => {
            cells_v = logn;
        }
        MeshSurfaceType::SurfaceOfExtrusion | MeshSurfaceType::SurfaceOfRevolution => {
            if let Some(c) = basis_curve_for_err(surface, kind) {
                let planar_like = c.is_line()
                    || (c.bspline_poles().is_some()
                        && c.nurbs_degree().map(|d| d < 2).unwrap_or(false));
                if planar_like {
                    if kind == MeshSurfaceType::SurfaceOfExtrusion {
                        cells_u = logn;
                    } else {
                        cells_v = logn;
                    }
                }
            }
            if kind == MeshSurfaceType::SurfaceOfExtrusion {
                cells_v = logn;
            }
        }
        MeshSurfaceType::BezierSurface | MeshSurfaceType::BSplineSurface => {
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
fn compute_err_factors(surface: &dyn Surface, kind: MeshSurfaceType, deflection: f64) -> (f64, f64) {
    let mut err_u = deflection * 10.0;
    let mut err_v = deflection * 10.0;
    match kind {
        MeshSurfaceType::Cylinder | MeshSurfaceType::Cone | MeshSurfaceType::Sphere | MeshSurfaceType::Torus => {}
        MeshSurfaceType::SurfaceOfExtrusion | MeshSurfaceType::SurfaceOfRevolution => {
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
        MeshSurfaceType::BezierSurface => {
            let du = surface.u_degree();
            let dv = surface.v_degree();
            if du > 2 {
                err_u /= du as f64;
            }
            if dv > 2 {
                err_v /= dv as f64;
            }
        }
        MeshSurfaceType::BSplineSurface => {
            let du = surface.u_degree();
            let dv = surface.v_degree();
            if du > 2 {
                err_u /= du as f64 * (surface.nb_u_intervals(0) + 1) as f64;
            }
            if dv > 2 {
                err_v /= dv as f64 * (surface.nb_v_intervals(0) + 1) as f64;
            }
        }
        MeshSurfaceType::Plane | MeshSurfaceType::OffsetSurface | MeshSurfaceType::OtherSurface => {
            err_u = 1.0;
            err_v = 1.0;
        }
    }
    (err_u, err_v)
}

fn basis_curve_for_err(
    surface: &dyn Surface,
    kind: MeshSurfaceType,
) -> Option<std::sync::Arc<dyn occt_geom::Curve>> {
    match kind {
        MeshSurfaceType::SurfaceOfExtrusion => surface.extrusion_basis_curve(),
        MeshSurfaceType::SurfaceOfRevolution => surface.revolution_basis_curve(),
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
