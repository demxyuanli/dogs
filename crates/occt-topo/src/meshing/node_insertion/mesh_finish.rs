//! `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` finish/optimize and triangle collection.


use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::Surface;


use super::super::delaun::Delaun;
use super::super::delaun_data::DelaunDataStructure;
use super::super::delaun_types::{DelaunTriangle, DelaunVertex, VertexState};
#[allow(unused_imports)]
use super::super::mesh_algo::*;
use super::super::deflection_control::DelaunayDeflectionControlMeshAlgo;
use super::super::parameters::MeshParameters;
use super::super::range_splitter::factory_uses_deflection_control;
use super::*;
use super::parametric::{geom_tool_cells_count, ParametricDelaun};


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
        if crate::debug_flags::post_dump() {
            crate::debug_flags::set_post_stage(1);
            crate::debug_flags::set_post_face(_face_index);
            let (tu, tv) = structure.tolerance();
            eprintln!(
                "POSTDUMP face={} struct_nodes={} indices={} tol=({:.17e},{:.17e}) sq={:.17e}",
                _face_index,
                structure.nb_nodes(),
                indices.len(),
                tu,
                tv,
                tu * tu + tv * tv
            );
        }
        let mut delaun = Delaun::new_with_data_cells(structure, &mut indices, cells_u, cells_v);
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

        // TEMP PREDCDUMP: the boundary-only base mesh, the state the probe dumps
        // as `OCCTDC POST_BASE` (`structdump_probe.cpp:235`). Same node ids as
        // OCCT only for the boundary nodes; the super-triangle nodes are still
        // registered at that point in OCCT and already gone here.
        if crate::debug_flags::pre_dc_dump() {
            let ds = delaun.result();
            eprintln!(
                "PREDCDUMP BASE f={} nodes={} domain={}",
                _face_index,
                ds.nb_nodes(),
                ds.elements_of_domain().len()
            );
            for &tid in ds.elements_of_domain() {
                let e = ds.element(tid);
                let n = ds.element_nodes(&e);
                eprintln!("PREDCDUMP BASE T {} {} {} {}", tid, n[0], n[1], n[2]);
            }
        }

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

        // TEMP PREDCDUMP: the state the deflection control starts from, i.e. the
        // probe's `OCCTDC POST_INSERT` stage (`structdump_probe.cpp` `dumpState`):
        // node 2D coordinates plus the domain triangles by structure index. The
        // port-side dump T-195 deleted only carried the input structure; this one
        // carries the triangulation, so the base mesh can be diffed triangle by
        // triangle before the control passes diverge.
        if crate::debug_flags::pre_dc_dump() {
            let ds = delaun.result();
            eprintln!(
                "PREDCDUMP f={} nodes={} domain={}",
                _face_index,
                ds.nb_nodes(),
                ds.elements_of_domain().len()
            );
            for i in 1..=ds.nb_nodes() as i32 {
                let v = ds.node(i);
                eprintln!(
                    "PREDCDUMP N {} {:.17e} {:.17e}",
                    i,
                    v.location.x(),
                    v.location.y()
                );
            }
            for &tid in ds.elements_of_domain() {
                let e = ds.element(tid);
                let n = ds.element_nodes(&e);
                eprintln!("PREDCDUMP T {} {} {} {}", tid, n[0], n[1], n[2]);
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

        self.collect_triangles(delaun, _face_index)
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
        if crate::debug_flags::post_dump() {
            crate::debug_flags::set_post_face(face_index);
        }
        let (ru, rv, delta) = match self.splitter.as_ref() {
            Some(sp) => (sp.range_u(), sp.range_v(), sp.delta()),
            None => return,
        };
        let mut ctrl = DelaunayDeflectionControlMeshAlgo::new(params);
        let mut view = ParametricDelaun { delaun, ru, rv, delta };
        let mut pass_no = 0usize;
        let max_defl = ctrl.post_process(&mut view, surface, face_deflection, &mut |mesh, nodes| {
            pass_no += 1;
            let pre_nodes = mesh.delaun.result().nb_nodes();
            let pre_dom = mesh.delaun.result().elements_of_domain().len();
            let inserted = mesh.insert_parametric(self, surface, nodes);
            // TEMP POSTDUMP: per-pass triangle set, the counterpart of the probe's
            // `OCCTDC P<n>` (`structdump_probe.cpp:369-373`, dumped right after
            // `insertNodes`). `POSTFACE` narrows it to one face index.
            if crate::debug_flags::post_dump()
                && crate::debug_flags::value("POSTFACE")
                    .and_then(|v| v.parse::<usize>().ok())
                    .map_or(true, |f| f == face_index)
            {
                let ds = mesh.delaun.result();
                eprintln!(
                    "POSTDUMP P{} face={} nodes={} domain={}",
                    pass_no,
                    face_index,
                    ds.nb_nodes(),
                    ds.elements_of_domain().len()
                );
                for &tid in ds.elements_of_domain() {
                    let e = ds.element(tid);
                    let n = ds.element_nodes(&e);
                    eprintln!("POSTDUMP P{} T {} {} {} {}", pass_no, tid, n[0], n[1], n[2]);
                }
                for id in 1..=ds.nb_nodes() as i32 {
                    let p = ds.node(id).location.coord;
                    eprintln!(
                        "POSTDUMP P{} N {} {:.17e} {:.17e}",
                        pass_no,
                        id,
                        p.x(),
                        p.y()
                    );
                }
            }
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
    pub(super) fn collect_triangles(
        &mut self,
        delaun: Delaun,
        face_index: usize,
    ) -> Result<TriangulationResult, String> {
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
        if crate::debug_flags::post_dump() {
            eprintln!(
                "POSTDUMP face={} nodes={} tris={}",
                face_index,
                result.nodes.len(),
                result.triangles.len()
            );
        }
        self.last_result = Some(result.clone());
        Ok(result)
    }

}
