//! `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` finish/optimize and triangle collection.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d, GpXY};
use occt_core::precision::SQUARE_CONFUSION;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::abs::Orientation;

use super::super::data_model::*;
use super::super::delaun::Delaun;
use super::super::delaun_data::DelaunDataStructure;
use super::super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::super::face_discret::Classifier;
#[allow(unused_imports)]
use super::super::mesh_algo::*;
use super::super::face_discret::PointState;
use super::super::deflection_control::{
    DelaunayDeflectionControlMeshAlgo, DeflectionMesh,
};
use super::super::parameters::MeshParameters;
use super::super::range_splitter::{
    mesh_surface_type, create_range_splitter, factory_uses_deflection_control,
    grab_params_of_internal_edges, RangeSplitter, MeshSurfaceType,
};
use super::*;
use super::parametric::{geom_tool_cells_count, ParametricDelaun};

use super::*;

impl DelaunayNodeInsertionMeshAlgo {
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
    pub(super) fn finish_mesh(
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
            .filter(|&i| structure.node(i).state != VertexState::Deleted)
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
        if delaun.is_failed() {
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
                        let p = ds.node(*i).location;
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
                    if crate::debug_flags::dc_dbg() {
                        let ds = delaun.result();
                        let mut edge_fp: u64 = 0;
                        for &tid in ds.elements_of_domain() {
                            let e = ds.element(tid);
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
                            let e = ds.element(*tid);
                            let n = ds.element_nodes(&e);
                            eprintln!("PORTT PRE T {} {} {} {}", tid, n[0], n[1], n[2]);
                        }
                        for i in 1..=ds.nb_nodes() as i32 {
                            let v = ds.node(i);
                            eprintln!(
                                "PORTN PRE {} {:.17e} {:.17e} {}",
                                i,
                                v.location.x(),
                                v.location.y(),
                                v.state.index()
                            );
                        }
                        for i in 1..=ds.nb_elements() as i32 {
                            let e = ds.element(i);
                            let n = ds.element_nodes(&e);
                            eprintln!("PORTEL PRE E {} {} {} {} {}", i, n[0], n[1], n[2], 0 as i32);
                        }
                    }
                    self.optimize_mesh(&mut delaun, surf.as_ref(), params, face_deflection, _face_index);
                    if crate::debug_flags::dc_dbg() {
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
    pub(super) fn optimize_mesh(
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
            if crate::debug_flags::dc_dbg() {
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
    pub(super) fn cells_count(&self, vertices_nb: i32) -> (i32, i32) {
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
    pub(super) fn collect_triangles(&mut self, delaun: Delaun) -> Result<TriangulationResult, String> {
        let ds = delaun.result();
        let mut used = self.used_nodes.clone();
        let mut triangles: Vec<DelaunTriangle> = Vec::new();
        for &id in ds.elements_of_domain() {
            let mut elem = ds.element(id);
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
            if ds.node(i).state == VertexState::Deleted {
                continue;
            }
            if let Some(&out) = used.get(&i) {
                let mut n = *ds.node(i);
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

}
