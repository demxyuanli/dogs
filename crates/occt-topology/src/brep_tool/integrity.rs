use crate::abs::ShapeType;
use crate::shape::TopoShape;

use super::BRepTool;

impl BRepTool {
    /// Structural integrity of a shape list. Always true in this port:
    /// TShape types are fixed at construction.
    pub fn check_shape_integrity(shapes: &[TopoShape]) -> bool {
        let _ = shapes;
        true
    }

    /// Human-readable type name.
    pub fn shape_type_str(t: ShapeType) -> &'static str {
        t.to_str()
    }
}
