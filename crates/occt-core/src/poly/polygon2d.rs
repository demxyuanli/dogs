//! 2D polygon — ordered list of points. Source: `Poly_Polygon2D.hxx`
use crate::gp::GpPnt2d;

#[derive(Debug, Clone)]
pub struct Polygon2D {
    pub nodes: Vec<GpPnt2d>,
    pub deflection: f64,
}

impl Polygon2D {
    pub fn new(nodes: Vec<GpPnt2d>) -> Self { Self { nodes, deflection: 0.0 } }
    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
}
