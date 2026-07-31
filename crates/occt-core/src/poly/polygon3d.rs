//! 3D polygon — ordered list of points. Source: `Poly_Polygon3D.hxx`
use crate::gp::GpPnt;

#[derive(Debug, Clone)]
pub struct Polygon3D {
    pub nodes: Vec<GpPnt>,
    pub deflection: f64,
}

impl Polygon3D {
    pub fn new(nodes: Vec<GpPnt>) -> Self { Self { nodes, deflection: 0.0 } }
    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn bounding_box(&self) -> crate::bnd::BndBox {
        let mut b = crate::bnd::BndBox::new();
        for n in &self.nodes { b.add_point(n); }
        b
    }
}
