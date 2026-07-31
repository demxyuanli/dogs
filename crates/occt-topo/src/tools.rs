//! Topology analysis tools. Source: `BRepTools`, `TopExp`
use crate::abs::ShapeType;
use crate::shape::{TopoShape, Vertex, Edge, Face, Solid};

/// Count occurrences of shapes by type in a compound.
/// In OCCT TopExp::MapShapes uses a map; here we use an explicit child list.
pub fn map_shapes_count(children: &[TopoShape]) -> std::collections::HashMap<ShapeType, usize> {
    let mut m = std::collections::HashMap::new();
    for c in children { *m.entry(c.shape_type()).or_insert(0) += 1; }
    m
}

/// Simple Euler characteristic for a connected mesh/solid region.
/// V - E + F (for 2-complex) or V - E + F - C (cells).
pub fn euler_characteristic(n_vertices: usize, n_edges: usize, n_faces: usize, n_cells: usize) -> i32 {
    (n_vertices as i32) - (n_edges as i32) + (n_faces as i32) - (n_cells as i32)
}

/// Classify a shape's bounding box (loose): vertices only → 0D; edges → 1D;
/// faces → 2D; solids → 3D.
pub fn shape_dimension(t: ShapeType) -> u8 {
    match t {
        ShapeType::Vertex => 0,
        ShapeType::Edge => 1,
        ShapeType::Wire => 1,
        ShapeType::Face => 2,
        ShapeType::Shell => 2,
        ShapeType::Solid => 3,
        ShapeType::CompSolid => 3,
        ShapeType::Compound => 3,
        ShapeType::Shape => 3,
    }
}

/// Check if a shape is closed (all boundary sub-shapes have even multiplicity).
/// Simplified: uses the `closed` flag stored on the TShape.
pub fn is_closed(s: &TopoShape) -> bool { s.closed() }

/// Vertex tools — access the 3D point (requires VertexShape data).
impl Vertex {
    /// Dummy point accessor — real impl stores point in VertexShape.
    pub fn point(&self) -> occt_core::gp::GpPnt { occt_core::gp::GpPnt::zero() }
}

/// Edge tools — parameter range.
impl Edge {
    pub fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
    pub fn last_parameter(&self) -> f64 { f64::INFINITY }
    pub fn is_degenerated(&self) -> bool { false }
    pub fn is_seam(&self) -> bool { false }
}

/// Face tools — natural restriction.
impl Face {
    pub fn is_natural_restriction(&self) -> bool { true }
}

/// Solid tools — volume via GProp later.
impl Solid { pub fn is_valid(&self) -> bool { true } }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euler_tetrahedron() {
        // V=4, E=6, F=4, C=1 → 4-6+4-1 = 1
        assert_eq!(euler_characteristic(4, 6, 4, 1), 1);
    }

    #[test]
    fn dimension() {
        assert_eq!(shape_dimension(ShapeType::Vertex), 0);
        assert_eq!(shape_dimension(ShapeType::Face), 2);
        assert_eq!(shape_dimension(ShapeType::Solid), 3);
    }
}
