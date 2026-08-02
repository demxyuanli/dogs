//! Port of OCCT incremental_mesh — Wave 1 BRepMesh.
//!
//! `BRepMesh_IncrementalMesh` is the meshing entry point: it builds the discrete
//! model (`ModelBuilder`), pre-processes it (`ModelPreProcessor`), then
//! discretizes edges/faces and assembles the crate's `ShapeMesh` triangle soup.
//!
//! Reference: `BRepMesh_IncrementalMesh.hxx/.cxx`, `BRepMesh_DiscretRoot.hxx`.
//!
//! # Wave 1 dependency note
//! The edge/face discretization siblings (`edge_discret.rs`, `face_discret.rs`)
//! are filled by parallel Wave 1 agents. Until they land, `discretize_edge` /
//! `discretize_face` below call the existing `wireframe` tessellators as minimal
//! in-file placeholders (`// W1-D 占位，波2/3 替换`), keeping the pipeline shape
//! correct.

use occt_core::gp::GpPnt;
use occt_core::poly::triangulation::Triangle;

use crate::abs::ShapeType;
use crate::mesh::{mesh_surface_area, ShapeMesh};
use crate::shape::{Edge, Face, TopoShape};
use crate::wireframe;

use super::data_model::{MeshModel, MeshStatus};
use super::model_builder::{ModelBuilder, ModelPreProcessor};
use super::parameters::MeshParameters;

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
    /// `ModelBuilder::build_model` → `ModelPreProcessor::perform` →
    /// per-edge/face discretization → assembled `ShapeMesh`.
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

        // 3. Discretize every edge/face and assemble the ShapeMesh.
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

    /// Assemble a [`ShapeMesh`] by tessellating every face of the model.
    ///
    /// Faces that triangulate successfully are cleared of the `Outdated` marker
    /// (their mesh now matches the requested deflection); faces that produce no
    /// triangles get a `Failure` status.
    fn build_shape_mesh(&mut self, model: &mut MeshModel) -> Result<ShapeMesh, String> {
        let source_shape = model.shape().map(|s| s.shape_type()).unwrap_or(ShapeType::Shape);
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();

        for i in 0..model.faces_nb() {
            let (face, deflection) = {
                let f = model.face(i).expect("face index in range");
                (f.face().clone(), f.deflection().max(self.parameters.deflection))
            };
            let (vs, ts) = self.discretize_face(&face, deflection);
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

    /// W1-D 占位，波2/3 替换: 调用 `crate::meshing::face_discret` 的 FaceDiscret
    /// 管线。当前委托 `wireframe::face_to_triangles` 完成 deflection-bounded
    /// UV-grid 三角化。
    fn discretize_face(&self, face: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
        wireframe::face_to_triangles(face, deflection)
    }

    /// W1-D 占位，波2/3 替换: 调用 `crate::meshing::edge_discret` 的 EdgeDiscret
    /// 管线。当前委托 `wireframe::edge_to_polyline` 完成 deflection-bounded 折线。
    #[allow(dead_code)]
    fn discretize_edge(&self, edge: &Edge, deflection: f64) -> Vec<GpPnt> {
        wireframe::edge_to_polyline(edge, deflection)
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
}
