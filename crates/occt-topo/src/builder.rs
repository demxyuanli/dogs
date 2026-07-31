//! Topological shape builder — constructs shapes from sub-shapes.
//! Source: `BRep_Builder`
use crate::abs::ShapeType;
use crate::shape::{TopoShape, Vertex, Edge, Wire, Face, Shell, Solid, Compound};

/// Builds topological shapes (mirrors BRep_Builder).
#[derive(Debug, Default)]
pub struct TopoBuilder;

impl TopoBuilder {
    pub fn new() -> Self { Self }

    /// Make an empty shape of given type.
    pub fn make_shape(&self, shape_type: ShapeType) -> TopoShape {
        TopoShape::new(shape_type)
    }

    /// Add sub-shape to compound.
    pub fn add(&self, shape: &mut TopoShape, sub: &TopoShape) {
        if let Ok(mut t) = shape.tshape.write() {
            t.nb_children += 1;
        }
        // Real BRep stores children in a list; here we track count only.
        let _ = sub;
    }

    /// Add edge to wire.
    pub fn add_edge(&self, wire: &mut Wire, edge: &Edge) { self.add(&mut wire.0, &edge.0); }
    /// Add wire to face.
    pub fn add_wire(&self, face: &mut Face, wire: &Wire) { self.add(&mut face.0, &wire.0); }
    /// Add face to shell.
    pub fn add_face(&self, shell: &mut Shell, face: &Face) { self.add(&mut shell.0, &face.0); }
    /// Add shell to solid.
    pub fn add_shell(&self, solid: &mut Solid, shell: &Shell) { self.add(&mut solid.0, &shell.0); }
    /// Add any shape to compound.
    pub fn add_compound(&self, comp: &mut Compound, shape: &TopoShape) { self.add(&mut comp.0, shape); }

    /// Make a vertex with a point.
    pub fn make_vertex(&self, _p: occt_core::gp::GpPnt, tolerance: f64) -> Vertex {
        let v = Vertex::new();
        v.set_tolerance(tolerance);
        v.set_free(true);
        v
    }

    /// Set the free flag.
    pub fn make_free(&self, shape: &TopoShape) { shape.set_free(true); }

    /// Delete a sub-shape (OCCT BRep_Builder::Remove). Stub: no-op in flat model.
    pub fn remove(&self, _shape: &mut TopoShape, _sub: &TopoShape) {}
}

impl Vertex {
    pub fn set_tolerance(&self, tol: f64) { self.0.tshape.write().unwrap().flags.check = tol > 1e-15; }
    pub fn set_point(&self, _p: occt_core::gp::GpPnt) {}
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpPnt;

    #[test]
    fn build_vertex() {
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(1.,2.,3.), 0.001);
        assert!(v.is_vertex());
        assert!(v.free());
    }

    #[test]
    fn make_compound() {
        let b = TopoBuilder::new();
        let mut c = Compound::new();
        let v = b.make_vertex(GpPnt::zero(), 0.0);
        b.add_compound(&mut c, &v.0);
        assert_eq!(c.shape_type(), ShapeType::Compound);
    }
}
