//! Port of OCCT incremental_mesh — BRepMesh integration.
//!
//! `BRepMesh_IncrementalMesh` is the meshing entry point: it builds the discrete
//! model (`ModelBuilder`), pre-processes it (`ModelPreProcessor`), then runs the
//! OCCT-style per-face pipeline — `EdgeDiscret` collects the boundary UV points
//! from the discretized edges, `Triangulator` triangulates them into 3D
//! triangles with deflection-controlled interior refinement — and assembles the
//! crate's `ShapeMesh` triangle soup.
//!
//! Reference: `BRepMesh_IncrementalMesh.hxx/.cxx`, `BRepMesh_DiscretRoot.hxx`.
//!
//! # Pipeline
//! `ModelBuilder::build_model` → `ModelPreProcessor::perform` →
//! per-face [`IncrementalMesh::triangulate_model_faces`]
//! (`EdgeDiscret` boundary UV → `Triangulator::triangulate`) → assembled
//! `ShapeMesh`.
//!
//! # Fallback
//! If any pipeline stage fails for a shape (unusual/degenerate geometry), the
//! pre-pipeline wireframe UV-grid tessellator
//! ([`IncrementalMesh::build_shape_mesh_wireframe`]) is used so `perform` never
//! panics.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::poly::triangulation::Triangle;
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::intpatch::refine_point_on_surface;
use crate::mesh::{mesh_surface_area, ShapeMesh};
use crate::shape::{Edge, Face, TopoShape};
use crate::wireframe;

use super::data_model::{MeshEdge, MeshModel, MeshStatus};
use super::edge_discret::{EdgeDiscret, MeshFace as UvFace};
use super::mesh_tool::MeshTool;
use super::model_builder::{ModelBuilder, ModelPreProcessor};
use super::parameters::MeshParameters;
use super::triangulator::{FaceTriangulation, Triangulator};

/// `BRepMesh_IncrementalMesh::initParameters` — validates the meshing
/// parameters and fills the `Interior`/`MinSize` defaults from the boundary
/// values. Returns an error message for invalid parameter values.
fn init_parameters(params: &mut MeshParameters) -> Result<(), String> {
    const CONFUSION: f64 = 1e-7;
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
    const ANGULAR: f64 = 1e-12;
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
    shape: Option<TopoShape>,
    is_done: bool,
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

/// Builds the mesh of a shape with respect to its correctly triangulated
/// parts (port of `BRepMesh_IncrementalMesh`).
///
/// Constructors follow OCCT and automatically call [`IncrementalMesh::perform`];
/// the resulting [`ShapeMesh`] is available via [`IncrementalMesh::mesh`].
#[derive(Debug)]
pub struct IncrementalMesh {
    root: DiscretRoot,
    parameters: MeshParameters,
    modified: bool,
    status_flags: i32,
    mesh: Option<ShapeMesh>,
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

        // 2. Pre-process: initialize per-entity deflection and status.
        ModelPreProcessor::perform(&mut model, &self.parameters);

        // 3. Discretize edges and triangulate faces through the OCCT-style
        //    pipeline (`EdgeDiscret` boundary UV → `Triangulator`).
        // ponytail: any pipeline failure falls back to the wireframe UV-grid
        // tessellator (the pre-pipeline behavior) so `perform` never panics for
        // unusual or degenerate shapes.
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
    fn build_shape_mesh(&self, model: &mut MeshModel) -> Result<ShapeMesh, String> {
        let source_shape = model.shape().map(|s| s.shape_type()).unwrap_or(ShapeType::Shape);
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();

        for ft in self.triangulate_model_faces(model)? {
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
    fn build_shape_mesh_wireframe(&mut self, model: &mut MeshModel) -> Result<ShapeMesh, String> {
        let source_shape = model.shape().map(|s| s.shape_type()).unwrap_or(ShapeType::Shape);
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();

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
    /// ponytail: the `FaceDiscret` uniform interior UV grid is not added here —
    /// the current `Delaun` port produces gapped triangulations when a dense set
    /// of Free interior points is combined with Frontier boundary vertices (a
    /// 3×3 grid works, anything denser leaves the boundary strips uncovered).
    /// `Triangulator::triangulate`'s own deflection-controlled refinement plays
    /// the interior-point role instead. Revisit once the Delaunay insertion
    /// handles dense Free point sets.
    fn triangulate_model_faces(&self, model: &mut MeshModel) -> Result<Vec<FaceTriangulation>, String> {
        let tool = MeshTool::new(self.parameters.clone());
        let triangulator = Triangulator::new(self.parameters.clone());

        let mut out: Vec<FaceTriangulation> = Vec::with_capacity(model.faces_nb());
        for i in 0..model.faces_nb() {
            let surface = {
                let f = model
                    .face(i)
                    .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?;
                f.surface()
                    .ok_or_else(|| format!("IncrementalMesh::triangulate_model_faces: face {i} has no surface"))?
            };

            // EdgeDiscret step: boundary UV polygon from the discretized edges.
            let uv_face = self.build_face_uv_polygon(model, i)?;
            let data = tool
                .extract_face(&uv_face, 1e-6)
                .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: face {i}: {e}"))?;

            // Triangulator step: 2D Delaunay + deflection-controlled refinement.
            let tri = triangulator
                .triangulate(surface.as_ref(), &data)
                .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: face {i}: {e}"))?;

            model
                .face_mut(i)
                .map_err(|e| format!("IncrementalMesh::triangulate_model_faces: {e}"))?
                .unset_status(MeshStatus::OUTDATED);
            out.push(tri);
        }

        if out.is_empty() {
            return Err("IncrementalMesh::triangulate_model_faces: model has no faces".to_string());
        }
        Ok(out)
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
    fn build_face_uv_polygon(&self, model: &MeshModel, face_index: usize) -> Result<UvFace, String> {
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
    fn wire_uv_chain(
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
    fn edge_uv_points(&self, edge: &MeshEdge, surface: &dyn Surface) -> Vec<GpPnt2d> {
        let deflection = edge.deflection();
        let deflection = if deflection.is_finite() && deflection > 0.0 && deflection < f64::MAX {
            deflection
        } else {
            self.parameters.deflection
        };
        self.discretize_edge(edge.edge(), deflection)
            .into_iter()
            .map(|p| {
                let (u, v) = project_uv(surface, &p);
                GpPnt2d::new(u, v)
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
    pub fn discretize_face(&self, face: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
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

/// Append an edge's UV points to a boundary chain, reversing the edge when
/// needed so the chain stays a continuous closed polygon.
fn stitch_chain(chain: &mut Vec<GpPnt2d>, mut edge_uv: Vec<GpPnt2d>) {
    if edge_uv.is_empty() {
        return;
    }
    if let Some(&last) = chain.last() {
        let first = edge_uv[0];
        let last_edge = edge_uv[edge_uv.len() - 1];
        // The edge's own last point joining the chain's last means the edge is
        // stored in the traversal direction already; otherwise flip it.
        if last.distance(&last_edge) < last.distance(&first) {
            edge_uv.reverse();
        }
    }
    for p in edge_uv {
        if chain.last().map_or(true, |q: &GpPnt2d| q.distance(&p) > 1e-9) {
            chain.push(p);
        }
    }
}

/// Whether a boundary chain is a closed loop (its last point equals its first).
fn chain_closed(chain: &[GpPnt2d]) -> bool {
    chain.len() >= 2 && chain[0].distance(&chain[chain.len() - 1]) < 1e-6
}

/// Invert a surface point to its `(u, v)` parameters via Newton iteration,
/// seeded from the surface range center (or the origin for unbounded ranges).
///
/// For analytic surfaces (planes, cylinders, …) the residual is near-linear, so
/// the iterate converges in a couple of steps.
fn project_uv(surface: &dyn Surface, p: &GpPnt) -> (f64, f64) {
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    let su = if u0.is_finite() && u1.is_finite() { 0.5 * (u0 + u1) } else { 0.0 };
    let sv = if v0.is_finite() && v1.is_finite() { 0.5 * (v0 + v1) } else { 0.0 };
    let (u, v, _) = refine_point_on_surface(surface, *p, su, sv, 12);
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::mesh_surface_area;
    use crate::primitives::BRepPrimBox;

    fn unit_box() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone()
    }

    #[test]
    fn perform_produces_triangles_for_box() {
        let shape = unit_box();
        let mut inc = IncrementalMesh::new();
        inc.set_shape(&shape);
        inc.change_parameters().deflection = 0.05;
        let mesh = inc.perform().expect("mesh produced");
        assert!(inc.is_done(), "IsDone set after perform");
        assert!(inc.is_modified(), "modified after perform");
        assert!(mesh.triangles.len() > 0, "triangles {}", mesh.triangles.len());
        assert!(mesh.vertices.len() >= 8, "vertices {}", mesh.vertices.len());
        // Unit box: surface area ≈ 6.
        let area = mesh_surface_area(&mesh);
        assert!((area - 6.0).abs() < 0.5, "box area {area}");
    }

    #[test]
    fn status_flags_are_no_error_for_clean_box() {
        let shape = unit_box();
        let params = MeshParameters { deflection: 0.05, ..MeshParameters::default() };
        let inc = IncrementalMesh::from_parameters(&shape, params);
        assert!(inc.mesh().is_some(), "mesh cached");
        assert_eq!(inc.get_status_flags(), 0, "flags {}", inc.get_status_flags());
    }

    #[test]
    fn constructor_performs_automatically() {
        let shape = unit_box();
        let inc = IncrementalMesh::from_deflection(&shape, 0.05, false, 0.5);
        assert!(inc.is_done(), "constructor performs");
        assert!(inc.mesh().is_some(), "constructor caches mesh");
        let m = inc.mesh().unwrap();
        assert!(m.triangles.len() > 0);
    }

    #[test]
    fn perform_without_shape_fails() {
        let mut inc = IncrementalMesh::new();
        assert!(inc.perform().is_err());
        assert!(!inc.is_done());
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let shape = unit_box();
        let mut inc = IncrementalMesh::new();
        inc.set_shape(&shape);
        inc.change_parameters().deflection = 0.0; // below Precision::Confusion()
        assert!(inc.perform().is_err());
        assert!(!inc.is_done());
    }

    #[test]
    fn incremental_mesh_to_shape_mesh_produces_triangles() {
        let shape = unit_box();
        let mesh = incremental_mesh_to_shape_mesh(&shape, 0.05).expect("mesh produced");
        assert!(mesh.triangles.len() > 0, "triangles {}", mesh.triangles.len());
        assert!(!mesh.vertices.is_empty(), "vertices non-empty");
        // The 6 box faces tile the unit surface: area ≈ 6.
        let area = mesh_surface_area(&mesh);
        assert!((area - 6.0).abs() < 0.5, "box area {area}");
    }

    #[test]
    fn pipeline_preserves_face_count() {
        let shape = unit_box();
        let mut inc = IncrementalMesh::new();
        inc.set_shape(&shape);
        inc.change_parameters().deflection = 0.05;
        // Build + pre-process the model and run the per-face pipeline directly.
        let mut model = ModelBuilder::build_model(&shape, inc.parameters()).expect("model built");
        ModelPreProcessor::perform(&mut model, inc.parameters());
        let tris = inc.triangulate_model_faces(&mut model).expect("pipeline ran");
        assert_eq!(tris.len(), 6, "box has 6 faces");
        for ft in &tris {
            assert!(ft.triangles.len() > 0, "face {} has triangles", ft.face_index);
        }
    }
}
