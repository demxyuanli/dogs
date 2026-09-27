//! Geometric access to topology — a port of `BRep_Tool`.
//!
//! Reads real geometry (curves, surfaces, tolerances) back from the
//! `GeometryRegistry` side-table that `TopoBuilder` populates. Unregistered
//! shapes fall back to OCCT-like defaults (origin point, unbounded ranges,
//! no curve/surface).

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

/// Namespace for BRep_Tool-style static queries.
pub struct BRepTool;

impl BRepTool {
    /// `BRep_Tool::IsClosed(E, F)` (`BRep_Tool.cxx:795-805`): true when the
    /// edge has a CurveOnSurface representation with two pcurves on the face's
    /// surface (`BRep_Tool.cxx:814-841`).
    ///
    /// UNPORTED: the triangulation arm (`BRep_Tool.cxx:803-804`, `:849-...`);
    /// the port does not attach a `Poly_Triangulation` to faces for this query.
    pub fn is_closed_edge_face(edge: &Edge, face: &Face) -> bool {
        let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
            return false;
        };
        // `BRep_Tool.cxx:819-822`: a plane is never a closed surface.
        if surf.gp_pln().is_some() {
            return false;
        }
        let face_key = GeometryRegistry::shape_key(&face.0);
        GeometryRegistry::global().edge_pcurves(&edge.0, face_key).len() > 1
    }

    /// The location of a vertex (from `BRep_TVertex`), or the origin when
    /// the vertex has no registered geometry.
    pub fn vertex_point(v: &Vertex) -> GpPnt {
        GeometryRegistry::global().vertex_point(&v.0)
    }

    /// Vertex tolerance.
    pub fn vertex_tolerance(v: &Vertex) -> f64 {
        GeometryRegistry::global().vertex_tolerance(&v.0)
    }

    /// The curve underlying an edge (`Handle(Geom_Curve)`), if registered.
    pub fn edge_curve(e: &Edge) -> Option<Arc<dyn Curve>> {
        GeometryRegistry::global().edge_curve(&e.0)
    }

    /// Parameter range of an edge. Unregistered edges report the unbounded
    /// range `(-inf, +inf)`.
    pub fn edge_parameters(e: &Edge) -> (f64, f64) {
        GeometryRegistry::global().edge_parameters(&e.0)
    }

    /// Edge tolerance.
    pub fn edge_tolerance(e: &Edge) -> f64 {
        GeometryRegistry::global().edge_tolerance(&e.0)
    }

    /// The surface underlying a face (`Handle(Geom_Surface)`), if registered.
    pub fn face_surface(f: &Face) -> Option<Arc<dyn Surface>> {
        GeometryRegistry::global().face_surface(&f.0)
    }

    /// Face tolerance.
    pub fn face_tolerance(f: &Face) -> f64 {
        GeometryRegistry::global().face_tolerance(&f.0)
    }

    /// Endpoint points of the edge: the curve evaluated at the first and last
    /// parameters. `None` when the edge has no registered curve.
    pub fn edge_vertices(e: &Edge) -> Option<(GpPnt, GpPnt)> {
        let curve = Self::edge_curve(e)?;
        let (first, last) = Self::edge_parameters(e);
        Some((curve.d0(first), curve.d0(last)))
    }

    /// Whether the edge is same-parameter (`BRep_TEdge` flag).
    pub fn same_parameter(e: &Edge) -> bool {
        GeometryRegistry::global().same_parameter(&e.0)
    }

    /// Whether the edge is degenerated (`BRep_TEdge` flag).
    pub fn is_degenerated(e: &Edge) -> bool {
        GeometryRegistry::global().is_degenerated_edge(&e.0)
    }

    /// The first (`i == 0`) or last parameter of the edge.
    pub fn parameter_on_edge(e: &Edge, i: usize) -> f64 {
        let (first, last) = Self::edge_parameters(e);
        if i == 0 {
            first
        } else {
            last
        }
    }

    /// A closed edge has its first parameter equal to its last within tolerance.
    pub fn is_closed_edge(e: &Edge) -> bool {
        let (first, last) = Self::edge_parameters(e);
        (first - last).abs() <= Self::edge_tolerance(e)
    }

    /// Whether the face is bounded by its natural parametric curves
    /// (`BRep_TFace` flag).
    pub fn natural_restriction(f: &Face) -> bool {
        GeometryRegistry::global().natural_restriction(&f.0)
    }

    /// Parametric bounds of the face: `(u_min, u_max, v_min, v_max)` from the
    /// surface's `u_range`/`v_range`. Unbounded fallback for unregistered faces.
    pub fn uv_bounds(f: &Face) -> (f64, f64, f64, f64) {
        match Self::face_surface(f) {
            Some(s) => {
                let (u1, u2) = s.u_range();
                let (v1, v2) = s.v_range();
                (u1, u2, v1, v2)
            }
            None => (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    /// World-coordinate vertex point: the registered local point transformed
    /// by the shape's accumulated `location` (`TopoDS_Shape::Location`).
    pub fn vertex_point_world(v: &Vertex) -> GpPnt {
        let t = v.0.location().transformation();
        Self::vertex_point(v).transformed(&t)
    }

    /// World-coordinate edge curve: the registered local curve transformed by
    /// the edge's location.
    pub fn edge_curve_world(e: &Edge) -> Option<Arc<dyn Curve>> {
        let t = e.0.location().transformation();
        Self::edge_curve(e).map(|c| Arc::from(c.transformed(&t)))
    }

    /// World-coordinate face surface: the registered local surface transformed
    /// by the face's location.
    pub fn face_surface_world(f: &Face) -> Option<Arc<dyn Surface>> {
        let t = f.0.location().transformation();
        Self::face_surface(f).map(|s| Arc::from(s.transformed(&t)))
    }

    /// World-coordinate face UV bounds (from the transformed surface).
    pub fn uv_bounds_world(f: &Face) -> (f64, f64, f64, f64) {
        match Self::face_surface_world(f) {
            Some(s) => {
                let (u1, u2) = s.u_range();
                let (v1, v2) = s.v_range();
                (u1, u2, v1, v2)
            }
            None => (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY),
        }
    }

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use occt_core::gp::{GpAx3, GpPln};

    fn edge() -> Edge {
        Edge(TopoShape::new(ShapeType::Edge))
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
    fn edge_parameters_are_unbounded() {
        let (first, last) = BRepTool::edge_parameters(&edge());
        assert_eq!(first, f64::NEG_INFINITY);
        assert_eq!(last, f64::INFINITY);
    }

    #[test]
    fn unbounded_edge_is_not_closed() {
        assert!(!BRepTool::is_closed_edge(&edge()));
    }

    #[test]
    fn shape_type_str_names_vertex() {
        assert_eq!(BRepTool::shape_type_str(ShapeType::Vertex), "Vertex");
    }

    #[test]
    fn empty_shape_list_is_integral() {
        let shapes: Vec<TopoShape> = Vec::new();
        assert!(BRepTool::check_shape_integrity(&shapes));
    }

    #[test]
    fn make_vertex_stores_point() {
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(1.0, 2.0, 3.0), 1e-3);
        assert!(BRepTool::vertex_point(&v).is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
        assert_eq!(BRepTool::vertex_tolerance(&v), 1e-3);
        clear_tree(&v.0);
    }

    #[test]
    fn make_edge_stores_curve() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(10.0, 0.0, 0.0));
        assert_eq!(BRepTool::edge_parameters(&e), (0.0, 10.0));
        let c = BRepTool::edge_curve(&e).expect("curve present");
        let mid = c.d0(5.0);
        assert!(mid.is_equal(&GpPnt::new(5.0, 0.0, 0.0)));
        assert!(BRepTool::same_parameter(&e));
        assert!(!BRepTool::is_degenerated(&e));
        assert_eq!(BRepTool::parameter_on_edge(&e, 0), 0.0);
        assert_eq!(BRepTool::parameter_on_edge(&e, 1), 10.0);
        clear_tree(&e.0);
    }

    #[test]
    fn segment_has_two_vertices_and_endpoints() {
        let b = TopoBuilder::new();
        let p1 = GpPnt::new(0.0, 0.0, 0.0);
        let p2 = GpPnt::new(3.0, 4.0, 0.0);
        let e = b.make_edge_segment(&p1, &p2);
        assert_eq!(e.0.tshape.read().unwrap().children.len(), 2);
        let (a, z) = BRepTool::edge_vertices(&e).expect("endpoints");
        assert!(a.is_equal(&p1));
        assert!(z.is_equal(&p2));
        clear_tree(&e.0);
    }

    #[test]
    fn make_face_plane_registers_surface() {
        let b = TopoBuilder::new();
        let f = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        assert!(BRepTool::face_surface(&f).is_some());
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(&f);
        assert_eq!(u1, f64::NEG_INFINITY);
        assert_eq!(u2, f64::INFINITY);
        assert_eq!(v1, f64::NEG_INFINITY);
        assert_eq!(v2, f64::INFINITY);
        assert!(BRepTool::natural_restriction(&f));
        clear_tree(&f.0);
    }

    #[test]
    fn child_trees_build() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let wire = b.make_wire(&[e1, e2]);
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let shell = b.make_shell(&[face.clone()]);
        let solid = b.make_solid(&[shell.clone()]);
        assert_eq!(wire.0.tshape.read().unwrap().children.len(), 2);
        assert_eq!(face.0.tshape.read().unwrap().children.len(), 0);
        assert_eq!(shell.0.tshape.read().unwrap().children.len(), 1);
        assert_eq!(solid.0.tshape.read().unwrap().children.len(), 1);
        clear_tree(&solid.0);
        clear_tree(&wire.0);
    }

    #[test]
    fn unregistered_shapes_fall_back() {
        let v = Vertex::new();
        assert!(BRepTool::vertex_point(&v).is_equal(&GpPnt::zero()));
        assert_eq!(BRepTool::vertex_tolerance(&v), 0.0);
        assert!(BRepTool::edge_curve(&edge()).is_none());
        assert!(BRepTool::edge_vertices(&edge()).is_none());
        assert!(BRepTool::face_surface(&Face::new()).is_none());
        assert_eq!(
            BRepTool::uv_bounds(&Face::new()),
            (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY)
        );
        assert!(BRepTool::natural_restriction(&Face::new()));
    }
}
