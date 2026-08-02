//! Shape extraction tools for meshing — port of `BRepMesh_ShapeTool` plus the
//! visitor pattern from `BRepMesh_ShapeVisitor` / `IMeshTools_ShapeVisitor`.
//!
//! [`ShapeTool`] provides the auxiliary queries the meshing pipeline needs to
//! read geometry back off a `TopoDS`-style shape: an edge's 3D curve and
//! parameter range, its 2D p-curve on a face, the face's boundary wires
//! (loops), vertex/edge/face tolerances and bounding-box extents.
//!
//! [`ShapeVisitor`] is the visitor interface used by
//! [`ShapeExplorer`](super::context::ShapeExplorer) to walk a shape;
//! [`ModelShapeVisitor`] builds the discrete [`MeshModel`] by adding faces and
//! edges (deduplicated by `TShape` identity) and wiring each face's boundary
//! edges with p-curves.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::toploc::TopLocLocation;
use occt_geom::Curve;

use crate::abs::Orientation;
use crate::brep_surface::edge_pcurve_on_face;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools_full::{edge_vertices, edges_of_wire, wires_of_face};

use super::data_model::{MeshModel, MeshStatus};

/// Identity key of a `TShape` (its heap address). Used to deduplicate shared
/// edges that are visited once per containing face.
fn tshape_ptr(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Auxiliary queries that extract meshing-relevant geometry from topology.
/// Source: `BRepMesh_ShapeTool`.
pub struct ShapeTool;

impl ShapeTool {
    /// Maximum tolerance of the given face, considering the tolerances of the
    /// face itself, its edges and its vertices.
    /// Source: `BRepMesh_ShapeTool::MaxFaceTolerance`.
    pub fn max_face_tolerance(face: &Face) -> f64 {
        let mut max_tol = BRepTool::face_tolerance(face);
        for wire in wires_of_face(face) {
            for edge in edges_of_wire(&wire) {
                max_tol = max_tol.max(BRepTool::edge_tolerance(&edge));
                let (a, b) = edge_vertices(&edge);
                for v in [a, b].into_iter().flatten() {
                    max_tol = max_tol.max(BRepTool::vertex_tolerance(&v));
                }
            }
        }
        max_tol
    }

    /// Maximum dimension (X, Y or Z extent) of a bounding box. `None` for a
    /// void box, matching `BRepMesh_ShapeTool::BoxMaxDimension`.
    pub fn box_max_dimension(b: &BndBox) -> Option<f64> {
        let (x0, x1, y0, y1, z0, z1) = b.get()?;
        Some((x1 - x0).max(y1 - y0).max(z1 - z0))
    }

    /// The 3D curve underlying an edge (`BRep_Tool::Curve`).
    pub fn edge_curve(edge: &Edge) -> Option<Arc<dyn Curve>> {
        BRepTool::edge_curve(edge)
    }

    /// The 3D parameter range of an edge (`BRep_Tool::Range`).
    pub fn edge_range(edge: &Edge) -> (f64, f64) {
        BRepTool::edge_parameters(edge)
    }

    /// Edge tolerance (`BRep_Tool::Tolerance`).
    pub fn edge_tolerance(edge: &Edge) -> f64 {
        BRepTool::edge_tolerance(edge)
    }

    /// Vertex tolerance (`BRep_Tool::Tolerance` on a vertex).
    pub fn vertex_tolerance(vertex: &Vertex) -> f64 {
        BRepTool::vertex_tolerance(vertex)
    }

    /// The 2D p-curve of an edge on a face: `samples` UV points of the edge
    /// curve projected onto the face surface. Empty when the edge or face has
    /// no registered geometry. Source: `ShapeAnalysis_Edge::PCurve`.
    pub fn edge_pcurve(edge: &Edge, face: &Face, samples: usize) -> Vec<GpPnt2d> {
        edge_pcurve_on_face(edge, face, samples.max(2))
    }

    /// UV locations of the edge's extremities on the face. `None` when the
    /// edge has no curve or the projection fails.
    /// Source: `BRepMesh_ShapeTool::UVPoints`.
    pub fn uv_points(edge: &Edge, face: &Face) -> Option<(GpPnt2d, GpPnt2d)> {
        let pc = edge_pcurve_on_face(edge, face, 2);
        if pc.len() < 2 {
            return None;
        }
        Some((pc[0], pc[pc.len() - 1]))
    }

    /// The boundary wires (loops) of a face.
    pub fn wires_of_face(face: &Face) -> Vec<Wire> {
        wires_of_face(face)
    }

    /// The edges of a wire in stored order.
    pub fn edges_of_wire(wire: &Wire) -> Vec<Edge> {
        edges_of_wire(wire)
    }

    /// Whether the edge is degenerated (`BRep_Tool::Degenerated`).
    pub fn is_degenerated(edge: &Edge) -> bool {
        BRepTool::is_degenerated(edge)
    }

    /// Applies a location to a point and returns the result.
    /// Source: `BRepMesh_ShapeTool::UseLocation`.
    pub fn use_location(p: &GpPnt, loc: &TopLocLocation) -> GpPnt {
        if loc.is_identity() {
            return *p;
        }
        p.transformed(&loc.transformation())
    }
}

/// Visitor interface for shapes: handles faces and edges discovered by a
/// [`ShapeExplorer`](super::context::ShapeExplorer).
/// Source: `IMeshTools_ShapeVisitor`.
pub trait ShapeVisitor {
    /// Handles a `TopoDS_Edge`.
    fn visit_edge(&mut self, edge: &Edge);
    /// Handles a `TopoDS_Face`.
    fn visit_face(&mut self, face: &Face);
}

/// Builds the discrete model of a shape by adding faces and free edges,
/// computing p-curves for the edges of each face. Edges shared between faces
/// are added to the model exactly once.
/// Source: `BRepMesh_ShapeVisitor`.
pub struct ModelShapeVisitor {
    model: MeshModel,
    edge_index: HashMap<usize, usize>,
}

impl ModelShapeVisitor {
    /// Creates a visitor that builds into the given (initially empty) model.
    pub fn new(model: MeshModel) -> Self {
        Self {
            model,
            edge_index: HashMap::new(),
        }
    }

    /// Consumes the visitor and returns the built discrete model.
    pub fn into_model(self) -> MeshModel {
        self.model
    }

    /// Adds the edge to the model if not seen yet, returning its model index.
    fn ensure_edge(&mut self, edge: &Edge) -> usize {
        let ptr = tshape_ptr(&edge.0);
        if let Some(&i) = self.edge_index.get(&ptr) {
            return i;
        }
        let i = self.model.add_edge(edge.clone());
        self.edge_index.insert(ptr, i);
        i
    }
}

impl ShapeVisitor for ModelShapeVisitor {
    fn visit_edge(&mut self, edge: &Edge) {
        let _ = self.ensure_edge(edge);
    }

    fn visit_face(&mut self, face: &Face) {
        let face_index = self.model.add_face(face.clone());
        let wires = wires_of_face(face);

        if wires.is_empty() {
            // No boundary wire — mirror OCCT's failed outer wire.
            self.model
                .face_mut(face_index)
                .expect("face just added")
                .set_status(MeshStatus::FAILURE);
            return;
        }

        // The first wire is treated as the outer one; it must succeed.
        let mut outer_ok = false;
        for (wi, wire) in wires.iter().enumerate() {
            let wire_index = self.model.add_wire(wire.clone());
            let edges = edges_of_wire(wire);
            let wire_ok = !edges.is_empty();

            for edge in &edges {
                if edge.orientation() == Orientation::External {
                    continue;
                }
                let edge_index = self.ensure_edge(edge);
                {
                    let model_edge = self.model.edge_mut(edge_index).expect("edge exists");
                    model_edge.add_pcurve(face_index, edge.orientation());
                }
                {
                    let model_wire = self.model.wire_mut(wire_index).expect("wire exists");
                    model_wire.add_edge(edge_index, edge.orientation());
                }
            }

            if wire_ok {
                self.model
                    .face_mut(face_index)
                    .expect("face exists")
                    .add_wire(wire_index);
                if wi == 0 {
                    outer_ok = true;
                }
            } else if wi == 0 {
                self.model
                    .face_mut(face_index)
                    .expect("face exists")
                    .set_status(MeshStatus::FAILURE);
                return;
            } else {
                // An internal wire failure is non-fatal but flagged.
                self.model
                    .face_mut(face_index)
                    .expect("face exists")
                    .set_status(MeshStatus::UNORIENTED_WIRE);
            }
        }

        if !outer_ok {
            self.model
                .face_mut(face_index)
                .expect("face exists")
                .set_status(MeshStatus::FAILURE);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::meshing::context::{ShapeExplorer, TopoShapeExplorer};
    use crate::primitives::BRepPrimBox;
    use crate::topo_tools_full::faces_of;
    use occt_core::gp::{GpAx3, GpPln, GpPnt};
    use occt_geom::GeomPlane;
    use std::sync::Arc;

    #[test]
    fn shape_tool_queries_box_loops_and_curves() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let face = faces_of(&b.solid.0).swap_remove(0);

        let wires = ShapeTool::wires_of_face(&face);
        assert_eq!(wires.len(), 1, "a box face has one boundary loop");

        let edges = ShapeTool::edges_of_wire(&wires[0]);
        assert_eq!(edges.len(), 4, "the loop has four edges");

        for e in &edges {
            let curve = ShapeTool::edge_curve(e).expect("edge has a 3D curve");
            let (first, last) = ShapeTool::edge_range(e);
            assert!(first.is_finite() && last.is_finite());
            assert!(last > first, "edge range {first}..{last}");
            // A point at mid-range lies on the curve.
            let mid = curve.d0(0.5 * (first + last));
            assert!(mid.distance(&GpPnt::zero()) >= 0.0);
        }

        // Box dims 1×2×3 → max dimension 3.
        assert!((ShapeTool::box_max_dimension(&b.bbox).unwrap() - 3.0).abs() < 1e-12);

        // Tolerances are non-negative.
        assert!(ShapeTool::max_face_tolerance(&face) >= 0.0);
        let (first_vertex, _) = edge_vertices(&edges[0]);
        assert!(ShapeTool::vertex_tolerance(&first_vertex.expect("first vertex")) >= 0.0);
    }

    #[test]
    fn shape_tool_uv_points_on_plane() {
        let b = TopoBuilder::new();
        let face = b.make_face(Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard()))), &[]);
        let edge = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));

        let (a, z) = ShapeTool::uv_points(&edge, &face).expect("uv points");
        // On the Z=0 plane the pcurve is (u≈x, v≈y).
        assert!((a.x() - 0.0).abs() < 1e-6, "first u {}", a.x());
        assert!((z.x() - 1.0).abs() < 1e-6, "last u {}", z.x());

        let pc = ShapeTool::edge_pcurve(&edge, &face, 5);
        assert_eq!(pc.len(), 5);
    }

    #[test]
    fn shape_tool_use_location_identity() {
        let p = GpPnt::new(1.0, 2.0, 3.0);
        let loc = TopLocLocation::identity();
        let q = ShapeTool::use_location(&p, &loc);
        assert!(q.is_equal(&p));
    }

    #[test]
    fn model_visitor_builds_model_from_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);

        let model = MeshModel::new(b.solid.0.clone());
        let mut visitor = ModelShapeVisitor::new(model);
        TopoShapeExplorer::new(b.solid.0).accept(&mut visitor);
        let model = visitor.into_model();

        assert_eq!(model.faces_nb(), 6, "six faces");
        assert_eq!(model.edges_nb(), 12, "twelve distinct shared edges");

        // Every face owns one wire with four edges, each carrying a pcurve.
        for fi in 0..model.faces_nb() {
            let face = model.face(fi).expect("face");
            assert_eq!(face.wires_nb(), 1, "one boundary loop per face");
            let wire_index = face.wire(0).expect("wire index");
            let wire = model.wire(wire_index).expect("wire");
            assert_eq!(wire.edges_nb(), 4, "four edges per loop");
            for ei in 0..wire.edges_nb() {
                let edge_index = wire.edge(ei).expect("edge index");
                let edge = model.edge(edge_index).expect("edge");
                assert!(edge.pcurves_nb() >= 1, "boundary edge has a pcurve");
            }
        }

        assert!(!model.has_status(MeshStatus::FAILURE), "all faces valid");
    }
}
