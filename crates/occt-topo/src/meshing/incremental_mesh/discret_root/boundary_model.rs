//! `IncrementalMesh::build_boundary_model` and the per-face triangulation
//! driver (`triangulate_model_faces` / `map_triangulation`).

use super::*;

impl IncrementalMesh {
    /// Diagnostic entry point mirroring the public stage API of
    /// `IMeshTools_Context` (`IMeshTools_Context.hxx:44` `BuildModel`, `:53`
    /// `DiscretizeEdges`, `:69` `HealModel`, `:85` `PreProcessModel`): returns
    /// the discrete model exactly as it stands before `DiscretizeFaces`
    /// (`:99`), so a caller can dump boundary / grain state without meshing.
    pub fn build_boundary_model(
        shape: &TopoShape,
        parameters: &MeshParameters,
    ) -> Result<MeshModel, String> {
        let mut params = parameters.clone();
        init_parameters(&mut params)?;
        let mut model = ModelBuilder::build_model(shape, &params)?;
        let driver = IncrementalMesh {
            parameters: params,
            ..IncrementalMesh::new()
        };
        driver.prepare_model_boundaries(&mut model)?;
        Ok(model)
    }

    pub(in crate::meshing::incremental_mesh) fn triangulate_model_faces(&self, model: &mut MeshModel) -> Result<Vec<FaceTriangulation>, String> {
        self.prepare_model_boundaries(model)?;

        // `BRepMesh_FaceDiscret::process` / `BRepMesh_BaseMeshAlgo::Perform`
        // catch `Standard_Failure` per face. OCCT has **no** failure-ratio rule:
        // the "abort the shape when more than 10% of faces need the wireframe
        // fallback" threshold that used to sit here was invented and is gone
        // (audit A19). The per-face UV-grid rescue below is still UNPORTED
        // (see its comment).
        use std::panic::{catch_unwind, AssertUnwindSafe};

        struct PendingFace {
            index: usize,
            tri: Option<FaceTriangulation>,
        }

        let nfaces = model.faces_nb();
        let mut pending: Vec<PendingFace> = Vec::with_capacity(nfaces);
        for i in 0..nfaces {
            {
                let f = model
                    .face(i)
                    .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?;
                if f.is_status(MeshStatus::FAILURE) || f.is_status(MeshStatus::REUSED) {
                    if crate::debug_flags::mbdiag() {
                        eprintln!(
                            "MBDIAG skip face={i} bits={} selftint={} failure={} reused={}",
                            f.status().bits(),
                            f.is_status(MeshStatus::SELF_INTERSECTING_WIRE),
                            f.is_status(MeshStatus::FAILURE),
                            f.is_status(MeshStatus::REUSED)
                        );
                    }
                    continue;
                }
            }
            let surface = model
                .face(i)
                .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?
                .surface()
                .ok_or_else(|| {
                    format!("IncrementalMesh::triangulate_model_faces: face {i} has no surface")
                })?
                .clone();

            let params = &self.parameters;
            let hook = std::panic::take_hook();
            std::panic::set_hook(Box::new(|_| {}));
            let outcome = catch_unwind(AssertUnwindSafe(|| {
                let mut algo = DelaunayNodeInsertionMeshAlgo::new();
                algo.perform(model, i, params)
            }));
            std::panic::set_hook(hook);

            let tri = match outcome {
                Ok(Ok(result)) => {
                    if crate::debug_flags::mbdiag() {
                        eprintln!(
                            "MBDIAG ok face={i} nodes={} tris={}",
                            result.nodes.len(),
                            result.triangles.len()
                        );
                    }
                    crate::meshing::node_insertion::diag_dump_result(&result, i);
                    if std::env::var("SWEEP_FACE").is_ok() {
                        let o = surface.d0(0.0, 0.0);
                        eprintln!(
                            "[facedbg] face={i} surf_orgin=({:.3},{:.3},{:.3}) nodes={} tris={}",
                            o.x(), o.y(), o.z(), result.nodes.len(), result.triangles.len()
                        );
                        for (k, n) in result.nodes.iter().enumerate() {
                            eprintln!(
                                "[facedbg]   n{k} uv=({:.4},{:.4}) p=({:.4},{:.4},{:.4})",
                                n.location.x(), n.location.y(),
                                n.p3d.x(), n.p3d.y(), n.p3d.z()
                            );
                        }
                    }
                    Self::map_triangulation(&result, surface.as_ref(), i).ok()
                }
                Ok(Err(e)) => {
                    // `BRepMesh_BaseMeshAlgo::process` (`BRepMesh_BaseMeshAlgo.cxx:52-59`):
                    // when `initDataStructure` fails, the face is left with
                    // `IMeshData_Failure` and no triangulation, and the remaining
                    // faces of the shape are still meshed. A face the algorithm
                    // itself marked failed is therefore skipped, not fatal; any
                    // other error still propagates.
                    if model
                        .face(i)
                        .map(|f| f.is_status(MeshStatus::FAILURE))
                        .unwrap_or(false)
                    {
                        if crate::debug_flags::mbdiag() {
                            eprintln!("MBDIAG err face={i} FAILURE=1 msg={e}");
                        }
                        None
                    } else {
                        if crate::debug_flags::mbdiag() {
                            eprintln!("MBDIAG fatal face={i} msg={e}");
                        }
                        return Err(format!(
                            "IncrementalMesh::triangulate_model_faces: face {i}: {e}"
                        ));
                    }
                }
                Err(_) => {
                    if crate::debug_flags::mbdiag() {
                        eprintln!("MBDIAG panic face={i}");
                    }
                    None
                }
            };

            pending.push(PendingFace { index: i, tri });
        }

        let processed = pending.len();

        let mut out: Vec<FaceTriangulation> = Vec::with_capacity(processed);
        for p in pending {
            // OCCT has no per-face alternative tessellator — a face whose
            // Delaunay fails simply keeps `IMeshData_Failure`
            // (`BRepMesh_BaseMeshAlgo.cxx:40-62`) and stays unmeshed.
            let tri = p.tri;

            if let Some(tri) = tri {
                model
                    .face_mut(p.index)
                    .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?
                    .unset_status(MeshStatus::OUTDATED);
                out.push(tri);
            } else {
                model
                    .face_mut(p.index)
                    .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?
                    .set_status(MeshStatus::FAILURE);
            }
        }

        if out.is_empty() {
            return Err("IncrementalMesh::triangulate_model_faces: model has no faces".to_string());
        }
        Ok(out)
    }

    /// Build the 3D triangle mesh from Delaunay nodes. `TriangulationResult.nodes[i]`
    /// is the structure node with index `i + 1`, and every triangle vertex index
    /// is 1-based.
    ///
    /// Frontier (boundary) nodes keep the 3D polyline from `Tessellate3d`
    /// (`BRepMesh_BaseMeshAlgo::registerNode` passes `ICurve::GetPoint`, and
    /// `collectNodes` writes `myNodesMap`, not a re-evaluation of the surface).
    /// Interior nodes have no 3D polyline; they are `surface.d0` at their UV.
    pub(super) fn map_triangulation(
        result: &TriangulationResult,
        surface: &dyn Surface,
        face_index: usize,
    ) -> Result<FaceTriangulation, String> {
        let mut vertices = Vec::with_capacity(result.nodes.len());
        let mut uv = Vec::with_capacity(result.nodes.len());
        for n in &result.nodes {
            uv.push(n.location);
            let p = if n.state == VertexState::Frontier {
                n.p3d
            } else {
                surface.d0(n.location.x(), n.location.y())
            };
            vertices.push(p);
        }
        let triangles: Vec<Triangle> = result
            .triangles
            .iter()
            .map(|t| {
                Triangle::new(
                    (t.vertex_indices[0] - 1) as usize,
                    (t.vertex_indices[1] - 1) as usize,
                    (t.vertex_indices[2] - 1) as usize,
                )
            })
            .collect();
        if triangles.is_empty() {
            return Err(format!(
                "IncrementalMesh::map_triangulation: face {face_index}: Delaunay produced no triangles"
            ));
        }
        Ok(FaceTriangulation { face_index, vertices, triangles, uv })
    }

}
