//! Type helper functions of the BOPDS package. Source: `BOPDS_Tools.hxx`.
use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{TopoShape, Vertex};

/// Returns the type of the shape.
pub fn shape_type(s: &TopoShape) -> ShapeType {
    s.shape_type()
}

/// True if the shape is a vertex.
pub fn is_vertex(s: &TopoShape) -> bool {
    s.shape_type() == ShapeType::Vertex
}

/// True if the shape is an edge.
pub fn is_edge(s: &TopoShape) -> bool {
    s.shape_type() == ShapeType::Edge
}

/// True if the shape is a face.
pub fn is_face(s: &TopoShape) -> bool {
    s.shape_type() == ShapeType::Face
}

/// True if the shape is a wire.
pub fn is_wire(s: &TopoShape) -> bool {
    s.shape_type() == ShapeType::Wire
}

/// True if the shape is a shell.
pub fn is_shell(s: &TopoShape) -> bool {
    s.shape_type() == ShapeType::Shell
}

/// True if the shape is a solid.
pub fn is_solid(s: &TopoShape) -> bool {
    s.shape_type() == ShapeType::Solid
}

/// The 3D point of a vertex (`BRep_Tool::Pnt`); `None` for non-vertices.
pub fn vertex_point(s: &TopoShape) -> Option<GpPnt> {
    if is_vertex(s) {
        Some(BRepTool::vertex_point(&Vertex(s.clone())))
    } else {
        None
    }
}

/// True if the type corresponds to a shape having a boundary
/// representation (vertex / edge / face).
pub fn has_brep(t: ShapeType) -> bool {
    matches!(t, ShapeType::Vertex | ShapeType::Edge | ShapeType::Face)
}

/// True if the type can be a participant of an interference.
pub fn is_interfering(t: ShapeType) -> bool {
    has_brep(t) || t == ShapeType::Solid
}

/// Converts a shape type to an integer (the OCCT `TopAbs_ShapeEnum` value).
pub fn type_to_integer(t: ShapeType) -> i32 {
    t as i32
}

/// Converts the combination of two shape types to the index of the
/// corresponding interference type (VV=0, VE=1, EE=2, VF=3, EF=4, FF=5,
/// VZ=6, EZ=7, FZ=8, ZZ=9), or `-1` for an incompatible combination.
pub fn type_to_integer2(t1: ShapeType, t2: ShapeType) -> i32 {
    let x = type_to_integer(t2) * 10 + type_to_integer(t1);
    match x {
        77 => 0,             // VV
        76 | 67 => 1,        // VE
        66 => 2,             // EE
        74 | 47 => 3,        // VF
        64 | 46 => 4,        // EF
        44 => 5,             // FF
        72 | 27 => 6,        // VZ
        62 | 26 => 7,        // EZ
        42 | 24 => 8,        // FZ
        22 => 9,             // ZZ
        _ => -1,
    }
}
