//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]
use super::prelude::*;
use super::*;

/// `BRepMesh_IncrementalMesh::initParameters` — validates the meshing
/// parameters and fills the `Interior`/`MinSize` defaults from the boundary
/// values. Returns an error message for invalid parameter values.

pub(super) fn init_parameters(params: &mut MeshParameters) -> Result<(), String> {
    pub(super) const CONFUSION: f64 = 1e-7;
    if params.deflection < CONFUSION {
        return Err("invalid parameter value".to_string());
    }
    if params.deflection_interior < CONFUSION {
        params.deflection_interior = params.deflection;
    }
    if params.min_size < CONFUSION {
        params.min_size = (MeshParameters::rel_min_size()
            * params.deflection.min(params.deflection_interior))
            .max(CONFUSION);
    }
    pub(super) const ANGULAR: f64 = 1e-12;
    if params.angle < ANGULAR {
        return Err("invalid parameter value".to_string());
    }
    if params.angle_interior < ANGULAR {
        params.angle_interior = 2.0 * params.angle;
    }
    Ok(())
}

/// Common interface for meshing algorithms (port of `BRepMesh_DiscretRoot`).
///
/// Holds the shape to triangulate and the `IsDone` flag; `IncrementalMesh`
/// extends it with meshing parameters and the produced mesh.
#[derive(Debug, Default)]
pub struct DiscretRoot {
    pub(super) shape: Option<TopoShape>,
    pub(super) is_done: bool,
}

impl DiscretRoot {
    /// Default constructor.
    pub fn new() -> Self {
        Self { shape: None, is_done: false }
    }

    /// Set the shape to triangulate.
    pub fn set_shape(&mut self, shape: &TopoShape) {
        self.shape = Some(shape.clone());
    }

    /// The shape being triangulated, if set.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.shape.as_ref()
    }

    /// `true` if triangulation was performed and succeeded.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// Set the `IsDone` flag.
    pub fn set_done(&mut self) {
        self.is_done = true;
    }

    /// Clear the `IsDone` flag.
    pub fn set_not_done(&mut self) {
        self.is_done = false;
    }
}

/// Per-face vertex/triangle counts from the last `perform`.
#[derive(Clone, Copy, Debug)]
pub struct FaceMeshStat {
    pub index: usize,
    /// Pointer-identity key of the source face, so a caller can pair the stat
    /// with its face without relying on traversal order (the mesh model skips
    /// faces, so a positional lookup is off by the number of skipped faces).
    pub shape_key: usize,
    pub surface: SurfaceType,
    pub vertices: usize,
    pub triangles: usize,
}

/// Builds the mesh of a shape with respect to its correctly triangulated
/// parts (port of `BRepMesh_IncrementalMesh`).
///
/// Constructors follow OCCT and automatically call [`IncrementalMesh::perform`];
/// the resulting [`ShapeMesh`] is available via [`IncrementalMesh::mesh`].
#[derive(Debug)]
pub struct IncrementalMesh {
    pub(super) root: DiscretRoot,
    pub(super) parameters: MeshParameters,
    pub(super) modified: bool,
    pub(super) status_flags: i32,
    pub(super) mesh: Option<ShapeMesh>,
    pub(super) face_stats: Vec<FaceMeshStat>,
}

impl IncrementalMesh {
    /// Default constructor — performs nothing until a shape is set and
    /// [`IncrementalMesh::perform`] is called.
    pub fn new() -> Self {
        Self {
            root: DiscretRoot::new(),
            parameters: MeshParameters::default(),
            modified: false,
            status_flags: 0,
            mesh: None,
            face_stats: Vec::new(),
        }
    }

    /// Constructor from shape + linear/angular deflection. Automatically calls
    /// [`IncrementalMesh::perform`].
    ///
    /// `is_relative`: when `true`, the deflection used for each edge is
    /// `lin_deflection * size_of_edge`; the deflection used for faces is the
    /// maximum deflection of their edges.
    pub fn from_deflection(
        shape: &TopoShape,
        lin_deflection: f64,
        is_relative: bool,
        ang_deflection: f64,
    ) -> Self {
        let parameters = MeshParameters {
            deflection: lin_deflection,
            angle: ang_deflection,
            relative: is_relative,
            ..MeshParameters::default()
        };
        Self::from_parameters(shape, parameters)
    }

    /// Constructor from shape + full [`MeshParameters`]. Automatically calls
    /// [`IncrementalMesh::perform`].
    pub fn from_parameters(shape: &TopoShape, parameters: MeshParameters) -> Self {
        let mut m = Self {
            root: DiscretRoot::new(),
            parameters,
            modified: false,
            status_flags: 0,
            mesh: None,
            face_stats: Vec::new(),
        };
        m.root.set_shape(shape);
        let _ = m.perform();
        m
    }

    /// The meshing parameters.
    pub fn parameters(&self) -> &MeshParameters {
        &self.parameters
    }

    /// Mutable access to the meshing parameters.
    pub fn change_parameters(&mut self) -> &mut MeshParameters {
        &mut self.parameters
    }

    /// `true` after the shape has been (re)meshed.
    pub fn is_modified(&self) -> bool {
        self.modified
    }

    /// Accumulated status flags faced during meshing (`GetStatusFlags`).
    pub fn get_status_flags(&self) -> i32 {
        self.status_flags
    }

    /// The produced triangle mesh, if meshing has run successfully.
    pub fn mesh(&self) -> Option<&ShapeMesh> {
        self.mesh.as_ref()
    }

    /// Per-face vertex/triangle counts from the last successful `perform`.
    pub fn face_stats(&self) -> &[FaceMeshStat] {
        &self.face_stats
    }

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
        //    FaceDiscret). Any pipeline failure falls back to the wireframe
        //    UV-grid tessellator (the pre-pipeline behavior).
        let mesh = match self.build_shape_mesh(&mut model) {
            Ok(m) => m,
            Err(_) => self.build_shape_mesh_wireframe(&mut model)?,
        };

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
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();
        self.face_stats.clear();

        for ft in self.triangulate_model_faces(model)? {
            let ty = model
                .face(ft.face_index)
                .ok()
                .and_then(|f| f.surface())
                .map(|s| classify_surface(s.as_ref()))
                .unwrap_or(SurfaceType::OtherSurface);
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
            let offset = vertices.len();
            vertices.extend(ft.vertices);
            for t in ft.triangles {
                triangles.push(Triangle::new(offset + t.n0, offset + t.n1, offset + t.n2));
            }
        }

        if triangles.is_empty() {
            return Err("IncrementalMesh::perform: no triangles generated".to_string());
        }
        Ok(ShapeMesh { vertices, triangles, source_shape })
    }

    /// Wireframe UV-grid fallback — the pre-pipeline behavior. Faces are
    /// tessellated with `wireframe::face_to_triangles` (deflection-bounded UV
    /// grid) when the OCCT-style pipeline errors.
    pub(super) fn build_shape_mesh_wireframe(&mut self, model: &mut MeshModel) -> Result<ShapeMesh, String> {
        let source_shape = model.shape().map(|s| s.shape_type()).unwrap_or(ShapeType::Shape);
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();
        self.face_stats.clear();

        for i in 0..model.faces_nb() {
            let (face, deflection) = {
                let f = model.face(i).expect("face index in range");
                (f.face().clone(), f.deflection().max(self.parameters.deflection))
            };
            let (vs, ts) = wireframe::face_to_triangles(&face, deflection);
            if ts.is_empty() {
                model.face_mut(i).expect("face index in range").set_status(MeshStatus::FAILURE);
                continue;
            }
            // Successfully (re)meshed — drop the outdated marker.
            model.face_mut(i).expect("face index in range").unset_status(MeshStatus::OUTDATED);

            let ty = model
                .face(i)
                .ok()
                .and_then(|f| f.surface())
                .map(|s| classify_surface(s.as_ref()))
                .unwrap_or(SurfaceType::OtherSurface);
            self.face_stats.push(FaceMeshStat {
                index: i,
                shape_key: crate::tgeometry::GeometryRegistry::shape_key(&face.0),
                surface: ty,
                vertices: vs.len(),
                triangles: ts.len(),
            });
            let offset = vertices.len();
            vertices.extend(vs);
            for t in ts {
                triangles.push(Triangle::new(offset + t.n0, offset + t.n1, offset + t.n2));
            }
        }

        if triangles.is_empty() {
            self.status_flags |= MeshStatus::FAILURE.bits() as i32;
            return Err("IncrementalMesh::perform: no triangles generated".to_string());
        }
        Ok(ShapeMesh { vertices, triangles, source_shape })
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
    pub(super) fn triangulate_model_faces(&self, model: &mut MeshModel) -> Result<Vec<FaceTriangulation>, String> {
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

        // `BRepMesh_FaceDiscret::process` / `BRepMesh_BaseMeshAlgo::Perform`
        // catch `Standard_Failure` per face. OCCT has **no** failure-ratio rule:
        // the "abort the shape when more than 10% of faces need the wireframe
        // fallback" threshold that used to sit here was invented and is gone
        // (audit A19). The per-face UV-grid rescue below is still UNPORTED
        // (see its comment).
        use std::panic::{catch_unwind, AssertUnwindSafe};

        struct PendingFace {
            index: usize,
            topo_face: Face,
            deflection: f64,
            tri: Option<FaceTriangulation>,
            needs_fallback: bool,
        }

        let nfaces = model.faces_nb();
        let mut pending: Vec<PendingFace> = Vec::with_capacity(nfaces);
        for i in 0..nfaces {
            {
                let f = model
                    .face(i)
                    .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?;
                if f.is_status(MeshStatus::FAILURE) || f.is_status(MeshStatus::REUSED) {
                    continue;
                }
            }
            let deflection = model
                .face(i)
                .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?
                .deflection()
                .max(self.parameters.deflection);
            let topo_face = model
                .face(i)
                .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?
                .face()
                .clone();
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

            let (tri, needs_fallback) = match outcome {
                Ok(Ok(result)) => {
                    let mapped = Self::map_triangulation(&result, surface.as_ref(), i);
                    let tri = mapped.ok();
                    let needs_fallback = tri.is_none();
                    (tri, needs_fallback)
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
                        (None, false)
                    } else {
                        return Err(format!(
                            "IncrementalMesh::triangulate_model_faces: face {i}: {e}"
                        ));
                    }
                }
                Err(_) => {
                    (None, true)
                }
            };

            pending.push(PendingFace {
                index: i,
                topo_face,
                deflection,
                tri,
                needs_fallback,
            });
        }

        let processed = pending.len();

        let mut out: Vec<FaceTriangulation> = Vec::with_capacity(processed);
        for p in pending {
            let tri = p.tri.or_else(|| {
                if p.needs_fallback {
                    // UNPORTED (audit A19 / task T-68): OCCT has no per-face
                    // alternative tessellator — a failed face simply keeps
                    // `IMeshData_Failure` (`BRepMesh_BaseMeshAlgo.cxx:52-62`).
                    // This UV-grid rescue is kept only because the ported
                    // pipeline still fails 169 of 1772 faces of
                    // `data/occ/T0M.stp` (measured 2026-09-20); delete it once
                    // those faces mesh through the faithful path.
                    Self::wireframe_face_triangulation(&p.topo_face, p.deflection, p.index)
                } else {
                    None
                }
            });

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

    /// Per-face wireframe fallback when Delaunay fails for one face.
    fn wireframe_face_triangulation(
        face: &Face,
        deflection: f64,
        face_index: usize,
    ) -> Option<FaceTriangulation> {
        let (vertices, triangles) = wireframe::face_to_triangles(face, deflection);
        if triangles.is_empty() {
            return None;
        }
        Some(FaceTriangulation {
            face_index,
            vertices,
            triangles,
            uv: Vec::new(),
        })
    }

    /// `BRepMesh_ModelHealer::process` FaceChecker + `amplifyEdges`
    /// (`BRepMesh_ModelHealer.cxx:234-211`).
    fn heal_self_intersecting_wires(&self, model: &mut MeshModel) -> Result<(), String> {
        let mut intersecting: HashMap<usize, std::collections::HashSet<usize>> = HashMap::new();
        for i in 0..model.faces_nb() {
            if let Some(edges) = self.face_intersecting_edges(model, i)? {
                intersecting.insert(i, edges);
            }
        }

        const AMP_ITERS: usize = 5;
        for _ in 0..AMP_ITERS {
            let mut edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
            for set in intersecting.values() {
                edges.extend(set.iter().copied());
            }
            if edges.is_empty() {
                break;
            }
            // `EdgeAmplifier` (`cxx:50-83`): remesh only the intersecting
            // edges. Do not run `discretize_pcurves` on the whole model --
            // that rebuilds every pcurve with healer MinSize and wipes
            // `connectClosestPoints` ends (`Tessellate2d(false)` keeps them).
            self.amplify_remesh_edges(model, &edges)?;
            let faces: Vec<usize> = intersecting.keys().copied().collect();
            intersecting.clear();
            for i in faces {
                ModelHealer::fix_face_boundaries(model, i)?;
                if let Some(edges) = self.face_intersecting_edges(model, i)? {
                    intersecting.insert(i, edges);
                }
            }
        }

        for i in intersecting.keys().copied() {
            let f = model.face_mut(i)?;
            f.set_status(MeshStatus::SELF_INTERSECTING_WIRE);
            f.set_status(MeshStatus::FAILURE);
        }
        Ok(())
    }

    /// `BRepMesh_FaceChecker::Perform` plus the 2-edge / 2-point branch
    /// (`BRepMesh_ModelHealer.cxx:248-278`).
    fn face_intersecting_edges(
        &self,
        model: &MeshModel,
        face_index: usize,
    ) -> Result<Option<std::collections::HashSet<usize>>, String> {
        // `BRepMesh_FaceChecker` / `SegmentsFiller` (`cxx:61-81`) walks the
        // discretized pcurve points, not a re-projected 3D polyline.
        let uv = match Self::face_pcurve_uv_face(model, face_index) {
            Some(uv) => uv,
            None => return Ok(None),
        };
        let edge_pts = Self::face_pcurve_edge_polylines(model, face_index);
        let mut checker = FaceChecker::from_pcurve_edges(&uv, &edge_pts, 1e-7);
        let mut edges = std::collections::HashSet::new();
        // `BRepMesh_ModelHealer.cxx:248-278`: FaceChecker failure marks only
        // GetIntersectingEdges; the 2-edge / 2-point branch marks those two.
        if !checker.perform() {
            edges.extend(checker.intersecting_edges().iter().copied());
        } else if let Some((e0, e1)) = Self::two_edge_two_point_wire(model, face_index) {
            edges.insert(e0);
            edges.insert(e1);
        }
        if edges.is_empty() {
            Ok(None)
        } else {
            Ok(Some(edges))
        }
    }

    /// One UV chain per wire from stored pcurve samples.
    /// Source: `BRepMesh_FaceChecker.cxx:61-81`.
    fn face_pcurve_uv_face(model: &MeshModel, face_index: usize) -> Option<UvFace> {
        let face = model.face(face_index).ok()?;
        let mut wires: Vec<Vec<GpPnt2d>> = Vec::new();
        for &wi in face.wires() {
            let wire = model.wire(wi).ok()?;
            let mut chain = Vec::new();
            for j in 0..wire.edges_nb() {
                let ei = wire.edge(j).ok()?;
                let ori = wire.edge_orientation(j).ok()?;
                let reverse_walk = wire.edge_reverse_walk(j);
                let edge = model.edge(ei).ok()?;
                let Some(pc) = edge.pcurve_for(face_index, ori) else {
                    continue;
                };
                let pts = pc.points();
                let reverse = !pc.is_forward() ^ reverse_walk;
                let seq: Vec<GpPnt2d> = if reverse {
                    pts.iter().rev().copied().collect()
                } else {
                    pts.to_vec()
                };
                for p in seq {
                    if chain
                        .last()
                        .map(|q: &GpPnt2d| q.distance(&p) <= 1e-9)
                        .unwrap_or(false)
                    {
                        continue;
                    }
                    chain.push(p);
                }
            }
            if chain.len() >= 2 {
                wires.push(chain);
            }
        }
        let outer_wire = wires.first()?.clone();
        let inner_wires = wires.into_iter().skip(1).collect();
        Some(UvFace {
            outer_wire,
            inner_wires,
            deflection: face.deflection(),
            id: face_index,
        })
    }

    /// One polyline per edge pcurve in `GetPoint` order
    /// (`BRepMesh_FaceChecker.cxx:65-81`).
    fn face_pcurve_edge_polylines(
        model: &MeshModel,
        face_index: usize,
    ) -> Vec<(usize, usize, Vec<GpPnt2d>)> {
        let mut out = Vec::new();
        let Ok(face) = model.face(face_index) else {
            return out;
        };
        for (wi_local, &wi) in face.wires().iter().enumerate() {
            let Ok(wire) = model.wire(wi) else {
                continue;
            };
            for j in 0..wire.edges_nb() {
                let Ok(ei) = wire.edge(j) else {
                    continue;
                };
                let Ok(ori) = wire.edge_orientation(j) else {
                    continue;
                };
                let Ok(edge) = model.edge(ei) else {
                    continue;
                };
                let Some(pc) = edge.pcurve_for(face_index, ori) else {
                    continue;
                };
                let pts = pc.points().to_vec();
                if pts.len() >= 2 {
                    out.push((wi_local, ei, pts));
                }
            }
        }
        out
    }

    fn two_edge_two_point_wire(model: &MeshModel, face_index: usize) -> Option<(usize, usize)> {
        let face = model.face(face_index).ok()?;
        if face.wires_nb() != 1 {
            return None;
        }
        let wi = face.wire(0).ok()?;
        let w = model.wire(wi).ok()?;
        if w.edges_nb() != 2 {
            return None;
        }
        let (e0, e1) = (w.edge(0).ok()?, w.edge(1).ok()?);
        let (o0, o1) = (w.edge_orientation(0).ok()?, w.edge_orientation(1).ok()?);
        let pc0 = model
            .edge(e0)
            .ok()
            .and_then(|e| e.pcurve_for(face_index, o0).map(|p| p.parameters_nb()));
        let pc1 = model
            .edge(e1)
            .ok()
            .and_then(|e| e.pcurve_for(face_index, o1).map(|p| p.parameters_nb()));
        if pc0 == Some(2) && pc1 == Some(2) {
            Some((e0, e1))
        } else {
            None
        }
    }

    /// `BRepMesh_ModelHealer::EdgeAmplifier` (`cxx:40-84`).
    fn amplify_remesh_edges(
        &self,
        model: &mut MeshModel,
        edges: &std::collections::HashSet<usize>,
    ) -> Result<(), String> {
        let mut list: Vec<usize> = edges.iter().copied().collect();
        list.sort_unstable();
        for ei in list {
            self.amplify_one_edge(model, ei)?;
        }
        Ok(())
    }

    /// One `EdgeAmplifier::operator()` (`cxx:50-83`): keep endpoints, insert
    /// interiors, `MinSize = Confusion`, min points = previous `ParametersNb`
    /// (plus one when a face outer wire has at most two edges).
    fn amplify_one_edge(&self, model: &mut MeshModel, edge_index: usize) -> Result<(), String> {
        let mut a_points_nb = model.edge(edge_index)?.discretization().parameters_nb();
        let n_pc = model.edge(edge_index)?.pcurves_nb();
        for p in 0..n_pc {
            let face_index = model.edge(edge_index)?.pcurve(p)?.face();
            let face = model.face(face_index)?;
            if face.wires_nb() == 0 {
                continue;
            }
            let wi = face.wire(0)?;
            if model.wire(wi)?.edges_nb() <= 2 {
                a_points_nb += 1;
                break;
            }
        }
        let nd = (model.edge(edge_index)?.deflection() / 3.0).max(occt_core::precision::CONFUSION);
        {
            let e = model.edge_mut(edge_index)?;
            e.set_deflection(nd);
            e.set_status(MeshStatus::REMESH);
            e.clear(true);
        }
        let min_size = occt_core::precision::CONFUSION;
        let topo_edge = Edge(model.edge(edge_index)?.edge().0.oriented(Orientation::Forward));
        let (a, b) = BRepTool::edge_parameters(&topo_edge);
        if !(a.is_finite() && b.is_finite() && b - a >= 1e-15) {
            return Ok(());
        }
        let curve = model.edge(edge_index)?.curve().ok_or_else(|| {
            format!("EdgeAmplifier: edge {edge_index} has no 3D curve")
        })?;
        let (mut def, ang, degenerated, is_free, same_param, same_range, ori) = {
            let e = model.edge(edge_index)?;
            (
                (0.5 * e.deflection()).max(1e-9),
                0.5 * e.angular_deflection(),
                e.degenerated(),
                e.is_free(),
                e.same_param(),
                e.same_range(),
                e.edge().0.orientation(),
            )
        };
        if ori == Orientation::Internal {
            def = (0.5 * def).max(1e-9);
        }
        let mut cos_value: Option<(
            std::sync::Arc<dyn occt_geom2d::curve::Curve2d>,
            std::sync::Arc<dyn occt_geom::Surface>,
        )> = None;
        let (tess_curve, tess_first, tess_last) = if !is_free && !same_param {
            let (face_index, orientation) = {
                let pc = model.edge(edge_index)?.pcurve(0)?;
                (pc.face(), pc.orientation())
            };
            let topo_face = model.face(face_index)?.face().clone();
            let oriented = Edge(model.edge(edge_index)?.edge().0.oriented(orientation));
            match (
                crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                BRepTool::face_surface(&topo_face),
            ) {
                (Ok(c2), Some(s)) => {
                    let (pf, pl) = crate::boptools_2d::curve_on_surface_range(&oriented, &topo_face)
                        .map(|(_, f, l)| (f, l))
                        .unwrap_or((a, b));
                    cos_value = Some((c2.clone(), s.clone()));
                    (
                        std::sync::Arc::new(CurveOnSurface::new(c2, s, pf, pl))
                            as std::sync::Arc<dyn occt_geom::Curve>,
                        pf,
                        pl,
                    )
                }
                _ => (curve.clone(), a, b),
            }
        } else {
            (curve.clone(), a, b)
        };
        let min_pts = tessellator_min_points(tess_curve.as_ref()).max(a_points_nb.max(2));
        let tess_edge = model.edge(edge_index)?.edge().clone();
        let mut tess = CurveTessellator::from_range_angular_min(
            tess_curve,
            tess_first,
            tess_last,
            def,
            ang,
            min_pts,
            min_size,
        );
        tess.add_internal_vertices(&tess_edge);
        let mut tess_params = tess.params().to_vec();
        if !is_free && same_param && same_range && tess_params.len() > 1 {
            let mut pcs = Vec::new();
            for p in 0..model.edge(edge_index)?.pcurves_nb() {
                let (face_index, orientation) = {
                    let pc = model.edge(edge_index)?.pcurve(p)?;
                    (pc.face(), pc.orientation())
                };
                let topo_face = model.face(face_index)?.face().clone();
                let oriented = Edge(model.edge(edge_index)?.edge().0.oriented(orientation));
                if let (Ok(c2), Some(s)) = (
                    crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                    BRepTool::face_surface(&topo_face),
                ) {
                    pcs.push((c2, s));
                }
            }
            split_by_deflection2d(curve.as_ref(), &mut tess_params, &pcs, def, min_size);
        }
        if degenerated {
            return Ok(());
        }
        let edge_tol = BRepTool::edge_tolerance(&model.edge(edge_index)?.edge());
        let n = tess_params.len();
        for k in 1..n.saturating_sub(1) {
            let t = tess_params[k];
            if let Some((c2, s)) = cos_value.as_ref() {
                let p = curve.d0(t);
                if !curve_tessellator_value_ok(c2.as_ref(), s.as_ref(), t, &p, edge_tol) {
                    continue;
                }
            }
            let p = curve.d0(t);
            let pos = model.edge(edge_index)?.discretization().parameters_nb().saturating_sub(1);
            model
                .edge_mut(edge_index)?
                .discretization_mut()
                .insert_point(pos, p, t)?;
        }
        // `Tessellate2d(false)` (`cxx:307-325`): interiors only, keep snapped ends.
        let params3d = model.edge(edge_index)?.discretization().parameters().to_vec();
        let pts3d = model.edge(edge_index)?.discretization().points().to_vec();
        let n3 = params3d.len();
        if n3 < 3 {
            return Ok(());
        }
        let n_pc = model.edge(edge_index)?.pcurves_nb();
        for p in 0..n_pc {
            let face_index = model.edge(edge_index)?.pcurve(p)?.face();
            let orientation = model.edge(edge_index)?.pcurve(p)?.orientation();
            let topo_face = model.face(face_index)?.face().clone();
            let oriented = Edge(model.edge(edge_index)?.edge().0.oriented(orientation));
            let Ok(pc) = crate::pcurve_full::make_pcurve_full(&oriented, &topo_face) else {
                continue;
            };
            let Some(surface) = BRepTool::face_surface(&topo_face) else {
                continue;
            };
            let (pf, pl) = crate::boptools_2d::curve_on_surface_range(&oriented, &topo_face)
                .map(|(_, f, l)| (f, l))
                .filter(|(f, l)| f.is_finite() && l.is_finite())
                .unwrap_or((a, b));
            let mut provider = if same_param {
                EdgeParameterProvider::new(pf, pl)
            } else {
                let old_first = params3d.first().copied().unwrap_or(a);
                let old_last = params3d.last().copied().unwrap_or(b);
                EdgeParameterProvider::with_stored(pf, pl, old_first, old_last)
            };
            let cos = CurveOnSurface::new(pc.clone(), surface.clone(), pf, pl);
            for k in 1..n3 - 1 {
                let t = params3d[k];
                let u = provider.parameter_of(t, &pts3d[k], &cos);
                let pos = model.edge(edge_index)?.pcurve(p)?.parameters_nb().saturating_sub(1);
                model
                    .edge_mut(edge_index)?
                    .pcurve_mut(p)?
                    .insert_point(pos, pc.d0(u), t)?;
            }
        }
        Ok(())
    }

    /// Discretize every edge's 3D curve once (the shared 3D polyline, OCCT
    /// `BRepMesh_EdgeDiscret::Tessellate3d`) and mirror its parameters into each
    /// pcurve (`Tessellate2d`), so the 2D and 3D boundary points stay aligned and
    /// adjacent faces meet exactly on a shared edge.
    pub(crate) fn discretize_pcurves(
        &self,
        model: &mut MeshModel,
        min_size: f64,
    ) -> Result<(), String> {
        // Collect the pcurve jobs first (an immutable pass), then populate
        // (a mutable pass) so `model` is never borrowed twice at once.
        let mut jobs: Vec<(usize, usize, Edge, Face, f64, f64)> = Vec::new();
        for i in 0..model.edges_nb() {
            let edge = model.edge(i)?;
            // The pcurve is a property of the (edge, face) pair in the edge's
            // *natural* curve direction, independent of how a wire orients it;
            // the wire orientation is applied later when the boundary UV is
            // collected. The 3D parameter range is computed on the Forward edge
            // (`a -> b`), but each pcurve carries its own wire orientation — a
            // seam edge has a Forward and a Reversed pcurve on the same face,
            // and `make_pcurve_full` selects the matching side from it.
            let fwd_edge = Edge(edge.edge().0.oriented(Orientation::Forward));
            let (a, b) = BRepTool::edge_parameters(&fwd_edge);
            if !(a.is_finite() && b.is_finite() && b - a >= 1e-15) {
                continue;
            }
            for p in 0..edge.pcurves_nb() {
                let face_index = edge.pcurve(p)?.face();
                let orientation = edge.pcurve(p)?.orientation();
                let topo_face = model.face(face_index)?.face().clone();
                let oriented_edge = Edge(edge.edge().0.oriented(orientation));
                jobs.push((i, p, oriented_edge, topo_face, a, b));
            }
        }

        // Tessellate each edge's 3D curve once (cached), then reuse its parameter
        // sequence for every pcurve of that edge.
        let mut edge_params: HashMap<usize, Vec<f64>> = HashMap::new();
        for (edge_index, pcurve_index, topo_edge, topo_face, a, b) in jobs {
            let params = match edge_params.get(&edge_index) {
                Some(p) => p.clone(),
                None => {
                    // `BRepMesh_EdgeDiscret::process` computes the edge's
                    // deflection (`BRepMesh_Deflection::ComputeDeflection`) before
                    // tessellating its 3D curve. `EdgeAmplifier` already lowered
                    // the stored deflection (`BRepMesh_ModelHealer.cxx:56-57`);
                    // do not overwrite it.
                    if !model.edge(edge_index)?.is_status(MeshStatus::REMESH) {
                        ModelPreProcessor::compute_edge_deflection(
                            model,
                            edge_index,
                            &self.parameters,
                        )?;
                    }
                    // `BRepMesh_EdgeDiscret::process` (`cxx:113-126`):
                    // CheckAndUpdateFlags on every pcurve before tessellation.
                    let n_pc_flags = model.edge(edge_index)?.pcurves_nb();
                    for p in 0..n_pc_flags {
                        let face_index = model.edge(edge_index)?.pcurve(p)?.face();
                        let topo_face = model.face(face_index)?.face().clone();
                        ShapeTool::check_and_update_flags(
                            model.edge_mut(edge_index)?,
                            &topo_face,
                        );
                    }
                    let curve = model.edge(edge_index)?.curve().ok_or_else(|| {
                        format!("IncrementalMesh::discretize_pcurves: edge {edge_index} has no 3D curve")
                    })?;
                    let (mut def, ang, degenerated, is_free, same_param, same_range, ori) = {
                        let e = model.edge(edge_index)?;
                        // `BRepMesh_CurveTessellator::init`: both linear and angular
                        // deflections are halved (`BRepMesh_CurveTessellator.cxx:73-74`).
                        (
                            (0.5 * e.deflection()).max(1e-9),
                            0.5 * e.angular_deflection(),
                            e.degenerated(),
                            e.is_free(),
                            e.same_param(),
                            e.same_range(),
                            e.edge().0.orientation(),
                        )
                    };
                    // INTERNAL edges take an extra half on the linear term (`cxx:75-77`).
                    if ori == Orientation::Internal {
                        def = (0.5 * def).max(1e-9);
                    }
                    // `CreateEdgeTessellator` (`BRepMesh_EdgeDiscret.cxx:57-63`):
                    // !SameParam uses `BRepAdaptor_Curve(edge, face)`.
                    let mut cos_value: Option<(
                        std::sync::Arc<dyn occt_geom2d::curve::Curve2d>,
                        std::sync::Arc<dyn occt_geom::Surface>,
                    )> = None;
                    let (tess_curve, tess_first, tess_last) = if !is_free && !same_param {
                        let (face_index, orientation) = {
                            let pc = model.edge(edge_index)?.pcurve(0)?;
                            (pc.face(), pc.orientation())
                        };
                        let topo_face = model.face(face_index)?.face().clone();
                        let oriented = Edge(
                            model.edge(edge_index)?.edge().0.oriented(orientation),
                        );
                        match (
                            crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                            BRepTool::face_surface(&topo_face),
                        ) {
                            (Ok(c2), Some(s)) => {
                                // `BRepAdaptor_Curve(E,F)` loads the pcurve on
                                // the COS representation range (`BRep_GCurve`
                                // First/Last after CheckPCurves clamp).
                                let (pf, pl) = crate::boptools_2d::curve_on_surface_range(
                                    &oriented,
                                    &topo_face,
                                )
                                .map(|(_, f, l)| (f, l))
                                .unwrap_or((a, b));
                                cos_value = Some((c2.clone(), s.clone()));
                                (
                                    std::sync::Arc::new(CurveOnSurface::new(c2, s, pf, pl))
                                        as std::sync::Arc<dyn occt_geom::Curve>,
                                    pf,
                                    pl,
                                )
                            }
                            _ => (curve.clone(), a, b),
                        }
                    } else {
                        (curve.clone(), a, b)
                    };
                    let min_pts = tessellator_min_points(tess_curve.as_ref());
                    let tess_edge = model.edge(edge_index)?.edge().clone();
                    let mut tess = CurveTessellator::from_range_angular_min(
                        tess_curve,
                        tess_first,
                        tess_last,
                        def,
                        ang,
                        min_pts,
                        min_size,
                    );
                    tess.add_internal_vertices(&tess_edge);
                    let mut tess_params = tess.params().to_vec();
                    // `splitByDeflection2d` (`BRepMesh_CurveTessellator.cxx:157-188`).
                    if !is_free && same_param && same_range && tess_params.len() > 1 {
                        let mut pcs = Vec::new();
                        let n_pc = model.edge(edge_index)?.pcurves_nb();
                        for p in 0..n_pc {
                            let (face_index, orientation) = {
                                let pc = model.edge(edge_index)?.pcurve(p)?;
                                (pc.face(), pc.orientation())
                            };
                            let topo_face = model.face(face_index)?.face().clone();
                            let oriented = Edge(
                                model.edge(edge_index)?.edge().0.oriented(orientation),
                            );
                            if let (Ok(c2), Some(s)) = (
                                crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                                BRepTool::face_surface(&topo_face),
                            ) {
                                pcs.push((c2, s));
                            }
                        }
                        split_by_deflection2d(
                            curve.as_ref(),
                            &mut tess_params,
                            &pcs,
                            def,
                            min_size,
                        );
                    }
                    // Store the shared 3D polyline. OCCT `Tessellate3d(theUpdateEnds=true)`
                    // replaces the two end points with the exact vertex points
                    // (`BRep_Tool::Pnt`) of the edge in its natural (Forward)
                    // curve direction — not the pcurve's wire orientation.
                    // Degenerated edges skip interior tessellator points
                    // (`BRepMesh_EdgeDiscret.cxx:250-270`).
                    let natural = Edge(model.edge(edge_index)?.edge().0.oriented(Orientation::Forward));
                    let (first_vertex, last_vertex) =
                        crate::topo_tools_full::edge_vertices(&natural);
                    let first_pnt = first_vertex.as_ref().map(|v| BRepTool::vertex_point(v));
                    let last_pnt = last_vertex.as_ref().map(|v| BRepTool::vertex_point(v));
                    let edge_tol = BRepTool::edge_tolerance(&model.edge(edge_index)?.edge());
                    let em = model.edge_mut(edge_index)?;
                    em.discretization_mut().clear(false);
                    let n = tess_params.len();
                    let keep: Vec<usize> = if degenerated && n >= 2 {
                        vec![0, n - 1]
                    } else {
                        let mut idxs: Vec<usize> = (0..n).collect();
                        // `Tessellate3d` (`cxx:256-259`): drop interiors that
                        // `CurveTessellator::Value` rejects.
                        if let Some((c2, s)) = cos_value.as_ref() {
                            idxs.retain(|&k| {
                                if k == 0 || k + 1 == n {
                                    return true;
                                }
                                let t = tess_params[k];
                                let p = curve.d0(t);
                                curve_tessellator_value_ok(
                                    c2.as_ref(),
                                    s.as_ref(),
                                    t,
                                    &p,
                                    edge_tol,
                                )
                            });
                        }
                        idxs
                    };
                    let mut params = Vec::with_capacity(keep.len());
                    for &k in &keep {
                        let t = tess_params[k];
                        let p = if k == 0 {
                            first_pnt.unwrap_or_else(|| curve.d0(t))
                        } else if k == n - 1 {
                            last_pnt.unwrap_or_else(|| curve.d0(t))
                        } else {
                            curve.d0(t)
                        };
                        em.discretization_mut().add_point(p, t);
                        params.push(t);
                    }
                    edge_params.insert(edge_index, params.clone());
                    params
                }
            };

            let pc = crate::pcurve_full::make_pcurve_full(&topo_edge, &topo_face)
                .map_err(|e| format!("IncrementalMesh::discretize_pcurves: {e}"))?;
            // `BRepMesh_EdgeDiscret::Tessellate2d` + `EdgeParameterProvider`: the
            // pcurve parameter is found by projecting the edge's 3D polyline point
            // onto the surface and then onto the pcurve (a non-SameParameter STEP
            // edge's 3D curve may lie off the surface). The 3D polyline was stored
            // above, so read it back for the projection.
            let surface = BRepTool::face_surface(&topo_face).ok_or_else(|| {
                "IncrementalMesh::discretize_pcurves: face has no surface".to_string()
            })?;
            let pts3d: Vec<GpPnt> = model.edge(edge_index)?.discretization().points().to_vec();
            let same_param = model.edge(edge_index)?.same_param();
            let pcurve = model.edge_mut(edge_index)?.pcurve_mut(pcurve_index)?;
            pcurve.clear(false);
            let _pcurve_face = pcurve.face();
            // The pcurve's parameter range over this edge. A bounded curve is its
            // own `[first, last]`; an unbounded line's range is the arc length
            // between the edge's vertices projected onto the surface then the
            // pcurve (`BRep_Builder::UpdateEdge`).
            // `BRepAdaptor_Curve(E,F).First/Last` is COS range (`hxx:77-84`).
            // After RemovePCurves, `BRep_Tool.cxx:367-372` CurveOnPlane sets
            // First/Last to the 3D `f,l`, not `Geom2dCircle` `[0, 2pi]`.
            let (pf, pl) = if let Some((_, f, l)) =
                crate::boptools_2d::curve_on_surface_range(&topo_edge, &topo_face)
            {
                if f.is_finite() && l.is_finite() {
                    (f, l)
                } else if a.is_finite() && b.is_finite() && b > a {
                    (a, b)
                } else if pc.first_parameter().is_finite() && pc.last_parameter().is_finite() {
                    (pc.first_parameter(), pc.last_parameter())
                } else {
                    (a, b)
                }
            } else if a.is_finite() && b.is_finite() && b > a {
                // Bounded 3D curve (circle / B-spline / trimmed line) whose pcurve
                // is unbounded (a 2D line): the pcurve is parameterized over the
                // edge's own range (`SameParameter`), so use `[a, b]` directly.
                // Projecting the vertices collapses for a *closed* edge whose two
                // vertices coincide (a torus minor seam's inner point), which would
                // shrink the pcurve parameter range to a single value.
                (a, b)
            } else if pc.first_parameter().is_finite() && pc.last_parameter().is_finite() {
                (pc.first_parameter(), pc.last_parameter())
            } else {
                let proj = |v: Option<crate::shape::Vertex>| -> Option<f64> {
                    let p = BRepTool::vertex_point(&v?);
                    let pr =
                        occt_geom::geom_api::project_point_on_surface(surface.as_ref(), &p, 1e-7)?;
                    Some(crate::pcurve_full::project_uv_on_curve2d(
                        pc.as_ref(),
                        GpPnt2d::new(pr.u, pr.v),
                    ))
                };
                let (v1, v2) = crate::topo_tools_full::edge_vertices(&topo_edge);
                match (proj(v1), proj(v2)) {
                    (Some(u1), Some(u2)) => (u1.min(u2), u1.max(u2)),
                    _ => (a, b),
                }
            };
            let mut provider = if same_param {
                EdgeParameterProvider::new(pf, pl)
            } else {
                let old_first = params.first().copied().unwrap_or(a);
                let old_last = params.last().copied().unwrap_or(b);
                EdgeParameterProvider::with_stored(pf, pl, old_first, old_last)
            };
            let cos = CurveOnSurface::new(pc.clone(), surface.clone(), pf, pl);
            for (&t, p3d) in params.iter().zip(pts3d.iter()) {
                // `BRepMesh_EdgeParameterProvider::Parameter` (`hxx:109-139`)
                // + `Tessellate2d` (`BRepMesh_EdgeDiscret.cxx:312-317`).
                let u = provider.parameter_of(t, p3d, &cos);
                pcurve.add_point(pc.d0(u), t);
            }
        }
        Ok(())
    }

    /// Build the UV polygon of a model face from its boundary edges.
    ///
    /// Each wire edge is discretized in 3D ([`Self::discretize_edge`],
    /// `EdgeDiscret`) and its polyline is projected onto the face surface to
    /// give the boundary UV points; edges are stitched so the polygon is a
    /// continuous closed chain regardless of the wire's edge orientations.
    ///
    /// The first edge's direction is ambiguous (the wire stores every edge in
    /// its natural curve direction, not necessarily the traversal direction), so
    /// the chain is built twice — natural and flipped — and the closed one is
    /// kept.
    pub(super) fn build_face_uv_polygon(&self, model: &MeshModel, face_index: usize) -> Result<UvFace, String> {
        let face = model
            .face(face_index)
            .map_err(|e| format!("IncrementalMesh::build_face_uv_polygon: {e}"))?;
        let surface = face.surface().ok_or("IncrementalMesh::build_face_uv_polygon: face has no surface")?;

        let mut outer_wire: Vec<GpPnt2d> = Vec::new();
        let mut inner_wires: Vec<Vec<GpPnt2d>> = Vec::new();
        for (w_pos, &wire_index) in face.wires().iter().enumerate() {
            let mut chain = self.wire_uv_chain(model, wire_index, surface.as_ref(), false);
            if !chain_closed(&chain) {
                chain = self.wire_uv_chain(model, wire_index, surface.as_ref(), true);
            }
            if chain.len() >= 2 {
                if w_pos == 0 {
                    outer_wire = chain;
                } else {
                    inner_wires.push(chain);
                }
            }
        }

        if outer_wire.is_empty() {
            return Err(format!(
                "IncrementalMesh::build_face_uv_polygon: face {face_index} has no boundary UV points"
            ));
        }
        Ok(UvFace {
            outer_wire,
            inner_wires,
            deflection: face.deflection(),
            id: face_index,
        })
    }

    /// UV boundary chain of one wire on a face surface.
    ///
    /// Each edge is discretized in 3D and projected to UV, then stitched into a
    /// continuous chain. `reverse_first` flips the first edge's direction,
    /// resolving the wire-start ambiguity.
    pub(super) fn wire_uv_chain(
        &self,
        model: &MeshModel,
        wire_index: usize,
        surface: &dyn Surface,
        reverse_first: bool,
    ) -> Vec<GpPnt2d> {
        let wire = match model.wire(wire_index) {
            Ok(w) => w,
            Err(_) => return Vec::new(),
        };
        let mut chain: Vec<GpPnt2d> = Vec::new();
        for j in 0..wire.edges_nb() {
            let edge_index = match wire.edge(j) {
                Ok(e) => e,
                Err(_) => return Vec::new(),
            };
            let edge = match model.edge(edge_index) {
                Ok(e) => e,
                Err(_) => return Vec::new(),
            };
            let mut uv = self.edge_uv_points(edge, surface);
            if reverse_first && j == 0 {
                uv.reverse();
            }
            stitch_chain(&mut chain, uv);
        }
        chain
    }

    /// UV points of one edge on a face.
    ///
    /// The 3D edge is discretized into a deflection-bounded polyline
    /// ([`Self::discretize_edge`]) and every polyline point is projected onto
    /// the face surface.
    pub(super) fn edge_uv_points(&self, edge: &MeshEdge, surface: &dyn Surface) -> Vec<GpPnt2d> {
        let deflection = edge.deflection();
        let deflection = if deflection.is_finite() && deflection > 0.0 && deflection < f64::MAX {
            deflection
        } else {
            self.parameters.deflection
        };
        self.discretize_edge(edge.edge(), deflection)
            .into_iter()
            .map(|p| {
                // Analytic ProjLib inverse first (exact even on unbounded-v
                // surfaces); Newton fallback for non-analytic (B-spline) faces.
                crate::pcurve_full::project_point_on_surface(surface, &p).unwrap_or_else(|| {
                    let (u, v) = project_uv(surface, &p);
                    GpPnt2d::new(u, v)
                })
            })
            .collect()
    }

    /// Discretize a 3D edge into a deflection-bounded polyline.
    ///
    /// Port of `BRepMesh_EdgeDiscret`/`GCPnts_UniformDeflection`: the edge's
    /// curve is sampled adaptively so every chord deviates from the curve by at
    /// most `deflection`. Both endpoints are always included. Returns an empty
    /// polyline when the edge has no registered curve or an unbounded range.
    pub fn discretize_edge(&self, edge: &Edge, deflection: f64) -> Vec<GpPnt> {
        let Some(curve) = BRepTool::edge_curve(edge) else {
            return Vec::new();
        };
        let (first, last) = BRepTool::edge_parameters(edge);
        if !(first.is_finite() && last.is_finite() && last >= first) {
            return Vec::new();
        }
        EdgeDiscret::discretize_edge(&curve, first, last, deflection.max(1e-9))
    }

    /// Discretize a single face into a 3D triangle soup through the OCCT-style
    /// pipeline (a temporary single-face model is triangulated).
    ///
    /// Falls back to the wireframe UV-grid tessellator on any pipeline error.
    pub fn discretize_face(&mut self, face: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
        let mut params = self.parameters.clone();
        params.deflection = deflection;
        let mut model = match ModelBuilder::build_model(&face.0, &params) {
            Ok(m) => m,
            Err(_) => return wireframe::face_to_triangles(face, deflection),
        };
        ModelPreProcessor::perform(&mut model, &params);
        match self.build_shape_mesh(&mut model) {
            Ok(m) => (m.vertices, m.triangles),
            // ponytail: fall back to the wireframe UV-grid tessellator.
            Err(_) => wireframe::face_to_triangles(face, deflection),
        }
    }
}

impl Default for IncrementalMesh {
    fn default() -> Self {
        Self::new()
    }
}

// DiscretRoot accessors re-exposed on IncrementalMesh for callers that do not
// want to reach through `.root`.
impl IncrementalMesh {
    /// Set the shape to triangulate (before [`IncrementalMesh::perform`]).
    pub fn set_shape(&mut self, shape: &TopoShape) {
        self.root.set_shape(shape);
        self.modified = false;
    }

    /// The shape being triangulated, if set.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.root.shape()
    }

    /// `true` if triangulation was performed and succeeded.
    pub fn is_done(&self) -> bool {
        self.root.is_done()
    }
}

/// Mesh a shape into a [`ShapeMesh`] through the OCCT-style BRepMesh pipeline
/// with the given linear deflection.
///
/// Convenience entry wrapping [`IncrementalMesh::perform`]; the returned
/// [`ShapeMesh`] holds the assembled triangle soup (vertices + triangles).
pub fn incremental_mesh_to_shape_mesh(shape: &TopoShape, deflection: f64) -> Result<ShapeMesh, String> {
    let mut inc = IncrementalMesh::new();
    inc.set_shape(shape);
    inc.change_parameters().deflection = deflection;
    inc.perform()
}
