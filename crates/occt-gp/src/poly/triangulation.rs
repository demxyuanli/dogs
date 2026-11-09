//! Triangle mesh triangulation. Source: `Poly_Triangulation.hxx`
use crate::gp::GpPnt;

/// Triangle in a triangulation. Indices into the node array (0-based).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Triangle { pub n0: usize, pub n1: usize, pub n2: usize }

impl Triangle { pub fn new(n0: usize, n1: usize, n2: usize) -> Self { Self { n0, n1, n2 } } }

/// Initial set of parameters a triangulation was built for.
/// Source: `Poly_TriangulationParameters.hxx:23-67`.
///
/// A negative value means "undefined" (`:46-51`); the accessors keep that
/// convention instead of normalizing to `Option<f64>`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriangulationParameters {
    pub deflection: f64,
    pub angle: f64,
    pub min_size: f64,
}

impl TriangulationParameters {
    /// `Poly_TriangulationParameters(theDeflection, theAngle, theMinSize)`
    /// (`hxx:30-37`).
    pub const fn new(deflection: f64, angle: f64, min_size: f64) -> Self {
        Self { deflection, angle, min_size }
    }

    /// `HasDeflection` (`hxx:45`).
    pub fn has_deflection(&self) -> bool { !(self.deflection < 0.0) }
    /// `HasAngle` (`hxx:48`).
    pub fn has_angle(&self) -> bool { !(self.angle < 0.0) }
    /// `HasMinSize` (`hxx:51`).
    pub fn has_min_size(&self) -> bool { !(self.min_size < 0.0) }
}

/// The fields of `Poly_Triangulation` that OCCT rewrites on an already shared
/// object: `myDeflection` and `myParams` (`Poly_Triangulation.hxx:398-401`).
#[derive(Debug, Clone, Copy)]
struct TriangulationMeta {
    deflection: f64,
    parameters: Option<TriangulationParameters>,
}

/// Triangulation: nodes (3D points) + triangles (index triples).
#[derive(Debug)]
pub struct Triangulation {
    pub nodes: Vec<GpPnt>,
    pub triangles: Vec<Triangle>,
    pub normals: Option<Vec<GpPnt>>, // per-node normals
    pub uv_nodes: Option<Vec<(f64, f64)>>, // per-node UV coordinates
    /// `Poly_Triangulation::myDeflection` + `myParams` (`hxx:398-401`).
    ///
    /// A `Poly_Triangulation` is a shared object in OCCT: its handle is the key
    /// under which `BRep_PolygonOnTriangulation` representations are stored
    /// (`BRep_PolygonOnTriangulation.cxx:36-49`), so `Deflection(...)`
    /// (`hxx:113`) and `Parameters(...)` (`hxx:118-121`) mutate the object in
    /// place and every holder of the handle sees the new value. This crate
    /// shares node data through `Arc<Triangulation>`, hence the cell: the
    /// setters take `&self` and keep the handle - and therefore the
    /// polygon-on-triangulation key - valid.
    meta: std::sync::Mutex<TriangulationMeta>,
}

impl Clone for Triangulation {
    /// `Poly_Triangulation`'s copy constructor copies the mesh and leaves
    /// `myParams` null (the generated `Poly_Triangulation.cxx` never touches
    /// it); the port copies the value into a fresh cell so the clone never
    /// aliases the original's mutable fields.
    fn clone(&self) -> Self {
        Self {
            nodes: self.nodes.clone(),
            triangles: self.triangles.clone(),
            normals: self.normals.clone(),
            uv_nodes: self.uv_nodes.clone(),
            meta: std::sync::Mutex::new(*self.meta.lock().unwrap()),
        }
    }
}

impl Triangulation {
    pub fn new(nodes: Vec<GpPnt>, triangles: Vec<Triangle>) -> Self {
        Self {
            nodes,
            triangles,
            normals: None,
            uv_nodes: None,
            meta: std::sync::Mutex::new(TriangulationMeta {
                deflection: 0.0,
                parameters: None,
            }),
        }
    }

    pub fn with_normals(mut self, normals: Vec<GpPnt>) -> Self { self.normals = Some(normals); self }
    pub fn with_uv(mut self, uv: Vec<(f64, f64)>) -> Self { self.uv_nodes = Some(uv); self }

    /// `Poly_Triangulation::Deflection()` (`Poly_Triangulation.hxx:111`).
    pub fn deflection(&self) -> f64 {
        self.meta.lock().unwrap().deflection
    }

    /// `Poly_Triangulation::Deflection(theDeflection)` (`hxx:113`): rewrites the
    /// deviation of *this* triangulation object, so a handle already shared with
    /// an edge polygon keeps pointing at the updated value.
    pub fn set_deflection(&self, deflection: f64) {
        self.meta.lock().unwrap().deflection = deflection;
    }

    /// `Poly_Triangulation::Parameters()` (`Poly_Triangulation.hxx:116`).
    pub fn parameters(&self) -> Option<TriangulationParameters> {
        self.meta.lock().unwrap().parameters
    }

    /// `Poly_Triangulation::Parameters(theParams)` (`hxx:118-121`): rewrites the
    /// parameter set of *this* triangulation object, so a handle already shared
    /// with an edge polygon keeps pointing at the updated parameter set.
    pub fn set_parameters(&self, params: Option<TriangulationParameters>) {
        self.meta.lock().unwrap().parameters = params;
    }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn nb_triangles(&self) -> usize { self.triangles.len() }

    /// Compute bounding box
    pub fn bounding_box(&self) -> crate::bnd::BndBox {
        let mut b = crate::bnd::BndBox::new();
        for n in &self.nodes { b.add_point(n); }
        b
    }
}
