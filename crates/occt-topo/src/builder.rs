//! Topological shape builder — constructs shapes from sub-shapes.
//! Source: `BRep_Builder`
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_geom::{Curve, GeomCircle, GeomLine, GeomPlane, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Compound, Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::tgeometry::{EdgeGeom, FaceGeom, GeometryRegistry, VertexGeom};

/// Builds topological shapes (mirrors BRep_Builder).
#[derive(Debug, Default)]
pub struct TopoBuilder;

impl TopoBuilder {
    pub fn new() -> Self { Self }

    /// Make an empty shape of given type.
    pub fn make_shape(&self, shape_type: ShapeType) -> TopoShape {
        TopoShape::new(shape_type)
    }

    /// Add sub-shape. Source: `TopoDS_Builder::Add` (TopoDS_Builder.cxx).
    ///
    /// Stores a full `TopoShape` (TShape + location + orientation). If the
    /// parent *view* is `REVERSED`, the stored child is reversed; if the
    /// parent location is not identity, the stored child is moved by
    /// `parent.location.inverted()`. The component TShape is frozen
    /// (`Free(false)`). Type compatibility is the OCCT `aTb` table; an
    /// incompatible pair is skipped (OCCT throws `TopoDS_UnCompatibleShapes`).
    pub fn add(&self, shape: &mut TopoShape, sub: &TopoShape) {
        // OCCT freezes the component first (also guards self-insertion).
        sub.set_free(false);
        if !shape.free() {
            return;
        }
        if !self.compatible(shape.shape_type(), sub.shape_type()) {
            return;
        }
        let mut child = sub.clone();
        if shape.orientation() == crate::abs::Orientation::Reversed {
            child.reverse();
        }
        if !shape.location().is_identity() {
            child.move_location(&shape.location().inverted());
        }
        if let Ok(mut t) = shape.tshape.write() {
            t.add_child(child);
            t.set_modified(true);
        }
    }

    /// OCCT `TopoDS_Builder::Add` compatibility bits: which parent types may
    /// contain a given component type.
    fn compatible(&self, parent: ShapeType, component: ShapeType) -> bool {
        use ShapeType::*;
        match component {
            Compound => parent == Compound,
            CompSolid => parent == Compound,
            Solid => parent == Compound || parent == CompSolid,
            Shell => parent == Compound || parent == Solid,
            Face => parent == Compound || parent == Shell,
            Wire => parent == Compound || parent == Face,
            Edge => parent == Compound || parent == Solid || parent == Wire,
            Vertex => parent == Compound || parent == Solid || parent == Face || parent == Edge,
            Shape => false,
        }
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

    /// Make a vertex with a point. The point and tolerance are registered in
    /// the geometry side-table (`BRep_TVertex`).
    pub fn make_vertex(&self, p: GpPnt, tolerance: f64) -> Vertex {
        let v = Vertex::new();
        v.set_point(p);
        v.set_tolerance(tolerance);
        v.set_free(true);
        v
    }

    /// Make an edge over `curve` on the parameter range `[first, last]`.
    /// Registers an `EdgeGeom` in the geometry side-table (`BRep_TEdge`).
    pub fn make_edge(&self, curve: Arc<dyn Curve>, first: f64, last: f64) -> Edge {
        let e = Edge::new();
        GeometryRegistry::global().set_edge(&e.0, EdgeGeom::new(curve, first, last));
        e
    }

    /// Make a straight segment from `p1` to `p2`, adding registered endpoint
    /// vertices as children (`BRepBuilderAPI_MakeEdge(P1, P2)`).
    pub fn make_edge_segment(&self, p1: &GpPnt, p2: &GpPnt) -> Edge {
        let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
            .expect("make_edge_segment: p1 and p2 must be distinct");
        let lin = GpLin::from_pnt_dir(*p1, dir);
        let mut e = self.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
        let v1 = self.make_vertex(*p1, 0.0);
        let v2 = self.make_vertex(*p2, 0.0);
        self.add_edge_vertices(&mut e, &v1, &v2);
        e
    }

    /// Make a straight segment edge from `p1` to `p2`, reusing the given
    /// endpoint vertices instead of creating fresh ones — so edges sharing an
    /// endpoint reference the same vertex `TShape` (needed for a meaningful
    /// Euler characteristic / vertex-sharing across faces).
    pub fn make_edge_segment_with_vertices(
        &self,
        p1: &GpPnt,
        p2: &GpPnt,
        v1: &Vertex,
        v2: &Vertex,
    ) -> Edge {
        let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
            .expect("make_edge_segment_with_vertices: p1 and p2 must be distinct");
        let lin = GpLin::from_pnt_dir(*p1, dir);
        let mut e = self.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
        self.add_edge_vertices(&mut e, v1, v2);
        e
    }

    /// `BRepLib_MakeEdge`: first vertex FORWARD, last vertex REVERSED.
    pub(crate) fn add_edge_vertices(&self, e: &mut Edge, v1: &Vertex, v2: &Vertex) {
        let mut a = v1.0.clone();
        a.set_orientation(Orientation::Forward);
        let mut b = v2.0.clone();
        b.set_orientation(Orientation::Reversed);
        self.add(&mut e.0, &a);
        self.add(&mut e.0, &b);
    }

    /// Make a circular arc edge in the plane `axis` with the given radius and
    /// parameter range (`BRepBuilderAPI_MakeEdge(gp_Circ, ...)`).
    pub fn make_edge_circle(&self, axis: &GpAx2, radius: f64, first: f64, last: f64) -> Edge {
        self.make_edge(Arc::new(GeomCircle::new(GpCirc::new(*axis, radius))), first, last)
    }

    /// Make a wire containing `edges` (`BRepBuilderAPI_MakeWire`).
    pub fn make_wire(&self, edges: &[Edge]) -> Wire {
        let mut wire = Wire::new();
        for e in edges {
            self.add_edge(&mut wire, e);
        }
        wire
    }

    /// Make a face over `surface`, adding `wires` as children. Registers a
    /// `FaceGeom` (`BRep_TFace`) with natural restriction (no trimming wires).
    pub fn make_face(&self, surface: Arc<dyn Surface>, wires: &[Wire]) -> Face {
        let mut face = Face::new();
        GeometryRegistry::global().set_face(&face.0, FaceGeom::new(surface));
        for w in wires {
            self.add_wire(&mut face, w);
        }
        face
    }

    /// Make an unbounded planar face (`BRepBuilderAPI_MakeFace(gp_Pln)`).
    pub fn make_face_plane(&self, pln: &GpPln) -> Face {
        self.make_face(Arc::new(GeomPlane::new(pln.clone())), &[])
    }

    /// Make a shell containing `faces` (`BRepBuilderAPI_MakeShell`).
    pub fn make_shell(&self, faces: &[Face]) -> Shell {
        let mut shell = Shell::new();
        for f in faces {
            self.add_face(&mut shell, f);
        }
        shell
    }

    /// Make a solid containing `shells` (`BRepBuilderAPI_MakeSolid`).
    pub fn make_solid(&self, shells: &[Shell]) -> Solid {
        let mut solid = Solid::new();
        for s in shells {
            self.add_shell(&mut solid, s);
        }
        solid
    }

    /// Make a compound of `shapes` (`BRepBuilderAPI_MakeCompound`).
    pub fn make_compound_of(&self, shapes: &[TopoShape]) -> Compound {
        let mut comp = Compound::new();
        for s in shapes {
            self.add_compound(&mut comp, s);
        }
        comp
    }

    /// Set the free flag.
    pub fn make_free(&self, shape: &TopoShape) { shape.set_free(true); }

    /// Delete a sub-shape (OCCT BRep_Builder::Remove). Stub: no-op in flat model.
    pub fn remove(&self, _shape: &mut TopoShape, _sub: &TopoShape) {}
}

impl Vertex {
    /// Set the vertex tolerance: sets the check flag and, when geometry is
    /// registered, keeps the side-table `VertexGeom` tolerance in sync.
    pub fn set_tolerance(&self, tol: f64) {
        self.0.tshape.write().unwrap().flags.check = tol > 1e-15;
        if let Some(mut g) = GeometryRegistry::global().vertex_geom(&self.0) {
            g.tolerance = tol;
            GeometryRegistry::global().set_vertex(&self.0, g);
        }
    }

    /// Register the vertex point in the geometry side-table, preserving any
    /// tolerance already registered (`BRep_TVertex::Pnt`).
    pub fn set_point(&self, p: GpPnt) {
        let tol = GeometryRegistry::global().vertex_geom(&self.0).map(|g| g.tolerance).unwrap_or(0.0);
        GeometryRegistry::global().set_vertex(&self.0, VertexGeom { point: p, tolerance: tol });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpAx3;

    fn nb_children(s: &TopoShape) -> usize {
        s.tshape.read().unwrap().children.len()
    }

    /// Release registry entries for a shape tree so tests don't leave stale
    /// geometry keyed by a freed Arc address in the process-wide side-table.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    #[test]
    fn build_vertex() {
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(1.,2.,3.), 0.001);
        assert!(v.is_vertex());
        assert!(v.free());
        clear_tree(&v.0);
    }

    #[test]
    fn make_compound() {
        let b = TopoBuilder::new();
        let mut c = Compound::new();
        let v = b.make_vertex(GpPnt::zero(), 0.0);
        b.add_compound(&mut c, &v.0);
        assert_eq!(c.shape_type(), ShapeType::Compound);
        clear_tree(&v.0);
        clear_tree(&c.0);
    }

    #[test]
    fn add_reversed_parent_stores_reversed_child() {
        let b = TopoBuilder::new();
        let mut wire = Wire::new();
        wire.0.set_orientation(crate::abs::Orientation::Reversed);
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        b.add_edge(&mut wire, &e);
        let stored = wire.0.tshape.read().unwrap().children[0].clone();
        assert!(stored.orientation().is_reversed());
        let viewed = crate::iterator::ShapeIterator::of_shape(&wire.0).next().unwrap();
        assert!(viewed.orientation().is_forward());
        assert!(!e.free());
        clear_tree(&wire.0);
    }

    #[test]
    fn wire_shell_solid_build_child_trees() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let wire = b.make_wire(&[e1, e2]);
        assert_eq!(nb_children(&wire.0), 2);

        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        assert_eq!(nb_children(&face.0), 0);

        let shell = b.make_shell(&[face.clone()]);
        assert_eq!(nb_children(&shell.0), 1);

        let solid = b.make_solid(&[shell.clone()]);
        assert_eq!(nb_children(&solid.0), 1);

        let comp = b.make_compound_of(&[wire.0.clone(), face.0.clone()]);
        assert_eq!(nb_children(&comp.0), 2);
        clear_tree(&comp.0);
    }
}
