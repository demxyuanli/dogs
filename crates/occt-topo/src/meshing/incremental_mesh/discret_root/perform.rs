//! `IncrementalMesh::perform` (`BRepMesh_IncrementalMesh.cxx`) and the
//! per-face assembly it drives (`build_shape_mesh`), with the triangle
//! orientation helpers and the model boundary preparation step.

use super::*;

impl IncrementalMesh {
    /// Performs meshing of the shape.
    ///
    /// Pipeline (OCCT `BRepMesh_IncrementalMesh::Perform`):
    /// `ModelBuilder::build_model` → `ModelPreProcessor::perform` → per-face
    /// discretization + triangulation → assembled `ShapeMesh`.
    ///
    /// Returns the triangle soup on success, or an error message (empty shape,
    /// invalid parameters, or no triangles produced). On success the result is
    /// cached in [`IncrementalMesh::mesh`] and the `IsDone` flag is set.
    pub fn perform(&mut self) -> Result<ShapeMesh, String> {
        self.root.set_not_done();
        self.status_flags = 0;

        let shape = match self.root.shape() {
            Some(s) => s.clone(),
            None => return Err("IncrementalMesh::perform: no shape set".to_string()),
        };

        // Validate/fix parameters (initParameters).
        init_parameters(&mut self.parameters)
            .map_err(|e| format!("IncrementalMesh::perform: {e}"))?;

        // 1. Build the discrete model from the topological shape.
        let mut model = ModelBuilder::build_model(&shape, &self.parameters)
            .map_err(|e| format!("IncrementalMesh::perform: {e}"))?;

        // 2. Discretize edges, heal, pre-process, and triangulate faces through
        //    the OCCT pipeline (EdgeDiscret → ModelHealer → ModelPreProcessor →
        //    FaceDiscret). A pipeline failure is reported, never rescued by an
        //    alternative tessellator: OCCT has none (`BRepMesh_BaseMeshAlgo.cxx`
        //    marks the face `IMeshData_Failure` and returns).
        let mesh = self.build_shape_mesh(&mut model)?;

        // 4. Accumulate status flags from faces and their wires
        //    (`BRepMesh_IncrementalMesh::Perform`).
        let mut flags = 0i32;
        for i in 0..model.faces_nb() {
            let f = model.face(i).expect("face index in range");
            flags |= f.status().bits() as i32;
            for &wi in f.wires() {
                flags |= model.wire(wi).expect("wire index in range").status().bits() as i32;
            }
        }
        self.status_flags = flags;

        self.modified = true;
        self.mesh = Some(mesh.clone());
        self.root.set_done();
        Ok(mesh)
    }

    /// Assemble a [`ShapeMesh`] through the OCCT-style pipeline.
    ///
    /// Every face is triangulated by [`Self::triangulate_model_faces`] and the
    /// per-face 3D nodes/triangles are concatenated into the flat triangle soup.
    pub(super) fn build_shape_mesh(&mut self, model: &mut MeshModel) -> Result<ShapeMesh, String> {
        let source_shape = model.shape().map(|s| s.shape_type()).unwrap_or(ShapeType::Shape);
        let reversed_faces = model.shape().map(Self::face_occurrence_reversed).unwrap_or_default();
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();
        self.face_stats.clear();
        self.face_triangulations.clear();

        for ft in self.triangulate_model_faces(model)? {
            self.face_triangulations.push(ft.clone());
            let ty = model
                .face(ft.face_index)
                .ok()
                .and_then(|f| f.surface())
                .map(|s| mesh_surface_type(s.as_ref()))
                .unwrap_or(MeshSurfaceType::OtherSurface);
            self.face_stats.push(FaceMeshStat {
                index: ft.face_index,
                shape_key: model
                    .face(ft.face_index)
                    .ok()
                    .map(|f| crate::tgeometry::GeometryRegistry::shape_key(&f.face().0))
                    .unwrap_or(0),
                surface: ty,
                vertices: ft.vertices.len(),
                triangles: ft.triangles.len(),
            });
            let reversed = reversed_faces.get(ft.face_index).copied().unwrap_or(false);
            let offset = vertices.len();
            vertices.extend(ft.vertices);
            for t in ft.triangles {
                let t = Self::oriented_triangle(&t, reversed);
                triangles.push(Triangle::new(offset + t.n0, offset + t.n1, offset + t.n2));
            }
        }

        if triangles.is_empty() {
            return Err("IncrementalMesh::perform: no triangles generated".to_string());
        }
        Ok(ShapeMesh { vertices, triangles, source_shape })
    }

    /// `RWMesh_FaceIterator::TriangleOriented` (`RWMesh_FaceIterator.hxx`):
    /// `Poly_Triangulation` holds the UV-Delaunay winding, i.e. the normal of the
    /// **surface**, not of the face occurrence, so a REVERSED face must be
    /// emitted with its last two nodes swapped (`(1,2,3) -> (1,3,2)`). The
    /// header's second term, `myIsMirrored` (an occurrence transform with a
    /// negative determinant), has no counterpart here: the port's shapes carry no
    /// location yet.
    pub(super) fn oriented_triangle(t: &Triangle, reversed: bool) -> Triangle {
        if reversed {
            Triangle::new(t.n0, t.n2, t.n1)
        } else {
            *t
        }
    }

    /// Per-face occurrence orientation of `shape`, in `faces_of` order.
    ///
    /// `ModelBuilder::build_model` adds faces by iterating `faces_of(shape)`
    /// (`wire_builder.rs:586`) and stores `Oriented(FORWARD)` copies
    /// (`wire_builder.rs:587` = `IMeshTools_ShapeExplorer.cxx:101`), so the
    /// discrete model keeps no trace of the occurrence orientation. Index `i`
    /// here is the same `i` as `MeshModel::face(i)`.
    pub(super) fn face_occurrence_reversed(shape: &TopoShape) -> Vec<bool> {
        crate::topo_tools::faces_of(shape)
            .iter()
            .map(|f| f.0.orientation().is_reversed())
            .collect()
    }

    /// Run the per-face pipeline over the model.
    ///
    /// For each face the boundary UV polygon is built from its discretized 3D
    /// edges ([`Self::build_face_uv_polygon`], the `EdgeDiscret` step) and the
    /// boundary UV point set is triangulated into 3D triangles
    /// ([`Triangulator::triangulate`], which runs the 2D Delaunay plus the
    /// OCCT deflection-controlled interior refinement). Faces that triangulate
    /// successfully are cleared of the `Outdated` marker.
    ///
    /// Returns the per-face triangulations so callers can inspect face coverage.
    ///
    /// The uniform interior UV grid (the range splitters' surface nodes) IS
    /// generated by the node-insertion step below: `perform` calls
    /// `list_surface_nodes`, which calls `RangeSplitter::generate_surface_nodes`
    /// (`node_insertion.rs:260` / `:269`). That matches OCCT, where
    /// `BRepMesh_DelaunayNodeInsertionMeshAlgo` always drives
    /// `getRangeSplitter().GenerateSurfaceNodes(...)`
    /// (`BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:65-69` when
    /// `IsPreProcessSurfaceNodes()` is true, `:94-100` otherwise). Only the base
    /// splitter yields no nodes (`range_splitter/stitching.rs:60-62`, OCCT's null list).
    /// Runs every stage `IMeshTools_Context` executes **before**
    /// `DiscretizeFaces`: `BuildModel` is the caller's `ModelBuilder` step,
    /// then `DiscretizeEdges` (`IMeshTools_Context.hxx:53`), `HealModel`
    /// (`:69`) and `PreProcessModel` (`:85`). Kept separate so the boundary /
    /// grain state can be inspected exactly as OCCT leaves it before faces are
    /// triangulated.
    pub(super) fn prepare_model_boundaries(&self, model: &mut MeshModel) -> Result<(), String> {
        // EdgeDiscret step (OCCT BRepMesh_EdgeDiscret): populate every edge
        // pcurve from the analytic MakePCurveOnFace (`make_pcurve_full`).
        self.discretize_pcurves(model, self.parameters.min_size)?;

        // ModelHealer step (`BRepMesh_ModelHealer::fixFaceBoundaries`): per wire
        // compute its deflection, snap the pcurve endpoints of adjacent wire edges
        // together (`connectClosestPoints`), then compute the face deflection.
        for i in 0..model.faces_nb() {
            let wire_indices: Vec<usize> = model.face(i)?.wires().to_vec();
            for wi in &wire_indices {
                ModelPreProcessor::compute_wire_deflection(model, *wi, &self.parameters)?;
            }
            ModelHealer::fix_face_boundaries(model, i)?;
            ModelPreProcessor::compute_face_deflection(model, i, &self.parameters)?;
        }

        // `BRepMesh_ModelHealer::process` FaceChecker + `amplifyEdges`
        // (`BRepMesh_ModelHealer.cxx:234-211`).
        self.heal_self_intersecting_wires(model)?;

        // ModelPreProcessor step (`BRepMesh_ModelPreProcessor::performInternal`):
        // seam-edge amplification + triangulation-consistency (no-ops on a fresh
        // model) before the faces are triangulated.
        ModelPreProcessor::perform(model, &self.parameters);
        Ok(())
    }

}
