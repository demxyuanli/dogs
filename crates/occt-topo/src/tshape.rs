//! TShape — underlying shape data (geometric + topological content).
//! Source: `TopoDS_TShape`
use std::sync::{Arc, RwLock};
use std::fmt;
use crate::abs::{ShapeType, ShapeFlags};
use occt_core::toploc::TopLocLocation;

/// Thread-safe shared handle to a TShape (replaces Handle(TopoDS_TShape)).
pub type HandleTShape = Arc<RwLock<TShape>>;

/// Shared shape data. OCCT's TopoDS_TShape (reference-counted, shared).
/// Use `Arc<TShape>` + `RwLock` for interior mutability like OCCT's Handle(TShape).
pub struct TShape {
    pub shape_type: ShapeType,
    pub flags: ShapeFlags,
    pub location: TopLocLocation,
    /// Real children list (OCCT's `TopoDS_TShape::myShapes`). A wire holds its
    /// edges, a face its wires, a solid its shells, a compound arbitrary shapes.
    pub children: Vec<HandleTShape>,
}

impl fmt::Debug for TShape {
    // Do not recurse into children (that would loop on cyclic topology);
    // print type + child count instead.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "TShape({} children={})", self.shape_type.to_str(), self.children.len())
    }
}

impl TShape {
    pub fn new(shape_type: ShapeType) -> Self {
        Self { shape_type, flags: ShapeFlags::default(), location: TopLocLocation::identity(), children: Vec::new() }
    }
    pub fn shape_type(&self) -> ShapeType { self.shape_type }
    pub fn free(&self) -> bool { self.flags.free }
    pub fn set_free(&mut self, v: bool) { self.flags.free = v; }
    pub fn closed(&self) -> bool { self.flags.closed }
    pub fn set_closed(&mut self, v: bool) { self.flags.closed = v; }
    pub fn infinite(&self) -> bool { self.flags.infinite }
    pub fn set_infinite(&mut self, v: bool) { self.flags.infinite = v; }
    pub fn nb_children(&self) -> usize { self.children.len() }
    pub fn set_location(&mut self, l: &TopLocLocation) { self.location = l.clone(); }
    pub fn child(&self, i: usize) -> Option<HandleTShape> { self.children.get(i).cloned() }
    pub fn add_child(&mut self, c: HandleTShape) { self.children.push(c); }
}

/// Vertex shape data — 3D point (from BRep_TVertex).
#[derive(Debug)]
pub struct VertexShape {
    pub base: TShape,
    pub point: occt_core::gp::GpPnt,
    pub tolerance: f64,
}

/// Edge shape data — curve + parameter range (from BRep_TEdge).
#[derive(Debug)]
pub struct EdgeShape {
    pub base: TShape,
    pub first: f64,
    pub last: f64,
    pub tolerance: f64,
    pub same_parameter: bool,
    pub same_range: bool,
    pub degenerated: bool,
}

/// Wire shape data.
#[derive(Debug)]
pub struct WireShape { pub base: TShape }

/// Face shape data — surface (from BRep_TFace).
#[derive(Debug)]
pub struct FaceShape {
    pub base: TShape,
    pub tolerance: f64,
    pub natural_restriction: bool,
}

/// Shell shape data.
#[derive(Debug)]
pub struct ShellShape { pub base: TShape }

/// Solid shape data.
#[derive(Debug)]
pub struct SolidShape { pub base: TShape }

/// Compound shape data — groups arbitrary sub-shapes.
#[derive(Debug)]
pub struct CompoundShape { pub base: TShape }

impl VertexShape { pub fn new(p: occt_core::gp::GpPnt) -> Self { Self { base: TShape::new(ShapeType::Vertex), point: p, tolerance: 0.0 } } }
impl EdgeShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Edge), first: f64::NEG_INFINITY, last: f64::INFINITY, tolerance: 0.0, same_parameter: false, same_range: false, degenerated: false } } }
impl WireShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Wire) } } }
impl FaceShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Face), tolerance: 0.0, natural_restriction: false } } }
impl ShellShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Shell) } } }
impl SolidShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Solid) } } }
impl CompoundShape { pub fn new() -> Self { Self { base: TShape::new(ShapeType::Compound) } } }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_defaults() {
        let v = VertexShape::new(occt_core::gp::GpPnt::new(1.,2.,3.));
        assert_eq!(v.base.shape_type(), ShapeType::Vertex);
        assert_eq!(v.point.x(), 1.0);
    }

    #[test]
    fn edge_bounds() {
        let e = EdgeShape::new();
        assert!(e.first.is_infinite());
        assert!(!e.degenerated);
    }
}
