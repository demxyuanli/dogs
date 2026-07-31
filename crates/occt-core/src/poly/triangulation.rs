//! Triangle mesh triangulation. Source: `Poly_Triangulation.hxx`
use crate::gp::GpPnt;

/// Triangle in a triangulation. Indices into the node array (0-based).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Triangle { pub n0: usize, pub n1: usize, pub n2: usize }

impl Triangle { pub fn new(n0: usize, n1: usize, n2: usize) -> Self { Self { n0, n1, n2 } } }

/// Triangulation: nodes (3D points) + triangles (index triples).
#[derive(Debug, Clone)]
pub struct Triangulation {
    pub nodes: Vec<GpPnt>,
    pub triangles: Vec<Triangle>,
    pub normals: Option<Vec<GpPnt>>, // per-node normals
    pub uv_nodes: Option<Vec<(f64, f64)>>, // per-node UV coordinates
    pub deflection: f64,
}

impl Triangulation {
    pub fn new(nodes: Vec<GpPnt>, triangles: Vec<Triangle>) -> Self {
        Self { nodes, triangles, normals: None, uv_nodes: None, deflection: 0.0 }
    }

    pub fn with_normals(mut self, normals: Vec<GpPnt>) -> Self { self.normals = Some(normals); self }
    pub fn with_uv(mut self, uv: Vec<(f64, f64)>) -> Self { self.uv_nodes = Some(uv); self }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn nb_triangles(&self) -> usize { self.triangles.len() }

    /// Compute bounding box
    pub fn bounding_box(&self) -> crate::bnd::BndBox {
        let mut b = crate::bnd::BndBox::new();
        for n in &self.nodes { b.add_point(n); }
        b
    }
}
