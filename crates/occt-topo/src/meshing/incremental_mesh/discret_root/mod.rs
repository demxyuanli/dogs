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
    pub surface: MeshSurfaceType,
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
    /// Per-face triangulation from the last successful `perform`, in
    /// `MeshModel::face` order.
    ///
    /// OCCT stores the `Poly_Triangulation` on the face itself
    /// (`BRepMesh_ShapeVisitor::Visit` -> `IMeshData_Face::SetTriangulation`,
    /// `BRep_Tool::Triangulation` reads it back), so keeping the per-face nodes,
    /// UV nodes and index triples here is the same store, not an extra one. The
    /// flat [`ShapeMesh`] soup is the concatenation of these.
    pub(super) face_triangulations: Vec<FaceTriangulation>,
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
            face_triangulations: Vec::new(),
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
            face_triangulations: Vec::new(),
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
    pub fn status_flags(&self) -> i32 {
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

    /// Per-face triangulation from the last successful `perform`, in
    /// `MeshModel::face` order.
    ///
    /// Same store as OCCT's `Poly_Triangulation` on the face: 3D nodes, the
    /// aligned UV nodes and 0-based index triples, per face.
    pub fn face_triangulations(&self) -> &[FaceTriangulation] {
        &self.face_triangulations
    }

}

mod amplify;
mod boundary_model;
mod discretize;
mod discretize_pcurves;
mod perform;
mod self_intersect;
mod traits;
mod uv_polygon;

pub use self::amplify::*;
pub use self::boundary_model::*;
pub use self::discretize::*;
pub use self::discretize_pcurves::*;
pub use self::perform::*;
pub use self::self_intersect::*;
pub use self::traits::*;
pub use self::uv_polygon::*;
