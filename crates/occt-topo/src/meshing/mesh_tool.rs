//! Port of OCCT mesh_tool — BRepMesh_MeshTool.
//!
//! In the OCCT pipeline `BRepMesh_MeshTool` manipulates the
//! `BRepMesh_DataStructureOfDelaun` (triangle legalization, link cleanup,
//! frontier handling). That structure (`delaun_data.rs`) is still a stub, so
//! this Wave-2 `MeshTool` is its model-side counterpart: it extracts the 2D UV
//! vertex set and the constraint (boundary) edges of a face from the discrete
//! model — either the UV polygons produced by the face discretizer
//! (`edge_discret::MeshFace`) or the pcurves of a `data_model::MeshModel` —
//! ready to be fed to the Delaunay triangulator. The self-contained
//! [`NodeClassifier`] helper (left/right of a constraint) is ported as-is.

use occt_core::gp::GpPnt2d;

use super::data_model;
use super::delaun_index::{Box2d, DelaunVertex, VertexTool};
use super::edge_discret;

/// Extracted 2D Delaunay input for one face.
///
/// `vertices` are the deduplicated UV nodes (indices are positions in this
/// vector); `constraint_edges` reference them as `(first, last)` pairs.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceMeshData {
    /// Deduplicated UV vertices.
    pub vertices: Vec<DelaunVertex>,
    /// Constraint (boundary) edges as vertex-index pairs.
    pub constraint_edges: Vec<(usize, usize)>,
}

impl FaceMeshData {
    /// Number of extracted UV vertices.
    pub fn vertices_nb(&self) -> usize {
        self.vertices.len()
    }

    /// Number of extracted constraint edges.
    pub fn edges_nb(&self) -> usize {
        self.constraint_edges.len()
    }

    /// Vertex with the given index.
    pub fn vertex(&self, index: usize) -> Option<&DelaunVertex> {
        self.vertices.get(index)
    }
}

/// Tool extracting Delaunay input (UV vertices + constraint edges) from a
/// discrete mesh model. Port of `BRepMesh_MeshTool` (model side).
#[derive(Debug, Clone)]
pub struct MeshTool {
    params: edge_discret::MeshParameters,
}

impl MeshTool {
    /// Builds a tool using `params` for deflection computations.
    pub fn new(params: edge_discret::MeshParameters) -> Self {
        Self { params }
    }

    /// The parameters in effect.
    pub fn parameters(&self) -> &edge_discret::MeshParameters {
        &self.params
    }

    /// Extracts UV vertices and constraint edges of a face from its UV polygons.
    ///
    /// `outer_wire` and each `inner_wires` polygon are treated as closed chains;
    /// coincident vertices within `tolerance` collapse to a single node.
    pub fn extract_face(
        &self,
        face: &edge_discret::MeshFace,
        tolerance: f64,
    ) -> Result<FaceMeshData, String> {
        let mut chains: Vec<Vec<GpPnt2d>> = Vec::new();
        if face.outer_wire.len() >= 2 {
            chains.push(face.outer_wire.clone());
        }
        for wire in &face.inner_wires {
            if wire.len() >= 2 {
                chains.push(wire.clone());
            }
        }
        Self::build_uv_data(&chains, tolerance)
            .ok_or_else(|| "BRepMesh_MeshTool::extract_face: face has no UV points".to_string())
    }

    /// Extracts UV vertices and constraint edges of a model face from its pcurves.
    ///
    /// Walks the face's wires → edges → pcurves (in wire order, honoring each
    /// pcurve orientation) and deduplicates the collected UV points. Requires
    /// the edge discretizers to have populated the pcurves.
    pub fn extract_model_face(
        &self,
        model: &data_model::MeshModel,
        face_index: usize,
        tolerance: f64,
    ) -> Result<FaceMeshData, String> {
        let face = model.face(face_index)?;
        let mut chains: Vec<Vec<GpPnt2d>> = Vec::new();
        for &wire_index in face.wires() {
            let wire = model.wire(wire_index)?;
            let mut chain: Vec<GpPnt2d> = Vec::new();
            for j in 0..wire.edges_nb() {
                let edge_index = wire.edge(j)?;
                let orientation = wire.edge_orientation(j)?;
                let edge = model.edge(edge_index)?;
                let Some(pcurve) = edge.pcurve_for(face_index, orientation) else {
                    continue;
                };
                if pcurve.is_forward() {
                    chain.extend_from_slice(pcurve.points());
                } else {
                    chain.extend(pcurve.points().iter().rev());
                }
            }
            if chain.len() >= 2 {
                chains.push(chain);
            }
        }
        Self::build_uv_data(&chains, tolerance).ok_or_else(|| {
            format!("BRepMesh_MeshTool::extract_model_face: face {face_index} has no UV points")
        })
    }

    /// Deflection for an edge: the edge's own stored deflection when set, else
    /// the params-based absolute deflection for a shape of `max_shape_size`.
    pub fn compute_edge_deflection(
        &self,
        edge: &data_model::MeshEdge,
        max_shape_size: f64,
    ) -> f64 {
        // `MeshEdge::deflection` defaults to `RealLast()` (f64::MAX) meaning
        // "unset" — such an edge falls back to the params-based absolute value.
        let stored = edge.deflection();
        if stored.is_finite() && stored > 0.0 && stored < f64::MAX {
            stored
        } else {
            edge_discret::Deflection::compute_absolute_deflection(
                self.params.deflection,
                max_shape_size,
            )
        }
    }

    /// Deduplicates UV points across the given ordered chains into a vertex set
    /// and builds the constraint edges (consecutive pairs within a chain, plus
    /// the closing edge of an open chain). Degenerate (self-loop) edges are
    /// dropped. Returns `None` when no point was supplied.
    fn build_uv_data(chains: &[Vec<GpPnt2d>], tolerance: f64) -> Option<FaceMeshData> {
        let mut raw: Vec<GpPnt2d> = Vec::new();
        for chain in chains {
            raw.extend_from_slice(chain);
        }
        if raw.is_empty() {
            return None;
        }

        let (mut min_u, mut min_v, mut max_u, mut max_v) =
            (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        for p in &raw {
            min_u = min_u.min(p.x());
            min_v = min_v.min(p.y());
            max_u = max_u.max(p.x());
            max_v = max_v.max(p.y());
        }

        let box2d = Box2d::new(min_u, min_v, max_u, max_v);
        let (cells_u, cells_v) = default_cells(raw.len());
        let mut tool = VertexTool::new(box2d, cells_u, cells_v);
        let indices: Vec<usize> = raw
            .iter()
            .map(|p| tool.add_vertex(DelaunVertex::new(p.x(), p.y()), tolerance))
            .collect();

        let mut edges: Vec<(usize, usize)> = Vec::new();
        let mut k = 0;
        for chain in chains {
            let start = k;
            for t in 0..chain.len() - 1 {
                let (a, b) = (indices[start + t], indices[start + t + 1]);
                if a != b {
                    edges.push((a, b));
                }
            }
            // Close an open chain (a closed polygon's duplicate corner already
            // maps to the first vertex, making this a self-loop).
            let (first, last) = (indices[start], indices[k + chain.len() - 1]);
            if first != last {
                edges.push((last, first));
            }
            k += chain.len();
        }

        Some(FaceMeshData { vertices: tool.into_vertices(), constraint_edges: edges })
    }
}

/// Functor separating points to the left / right of a constraint edge.
///
/// Port of `BRepMesh_MeshTool::NodeClassifier`. `is_above` reports whether a
/// point lies on the "upper" side of the directed constraint, matching the
/// OCCT sign convention (`Direction().X() > 0` selects the `cross < 0` side).
#[derive(Debug, Clone)]
pub struct NodeClassifier {
    ax: f64,
    ay: f64,
    dx: f64,
    dy: f64,
    sign: bool,
}

impl NodeClassifier {
    /// Builds a classifier over the constraint from `(x1, y1)` to `(x2, y2)`.
    pub fn new(x1: f64, y1: f64, x2: f64, y2: f64) -> Self {
        Self { ax: x1, ay: y1, dx: x2 - x1, dy: y2 - y1, sign: (x2 - x1) > 0.0 }
    }

    /// Builds a classifier over the constraint between two UV vertices.
    pub fn from_vertices(a: &DelaunVertex, b: &DelaunVertex) -> Self {
        Self::new(a.u, a.v, b.u, b.v)
    }

    /// `true` when the point lies strictly on the classifier's upper side of
    /// the constraint; points on the constraint itself yield `false`.
    pub fn is_above(&self, u: f64, v: f64) -> bool {
        let nu = u - self.ax;
        let nv = v - self.ay;
        if nu * nu + nv * nv > 1e-9 {
            let cross = nu * self.dy - nv * self.dx;
            if cross.abs() > 1e-9 {
                return if self.sign { cross < 0.0 } else { cross > 0.0 };
            }
        }
        false
    }
}

/// Default grid split for `n` points: roughly square, at least `1 × 1`.
fn default_cells(n: usize) -> (usize, usize) {
    let s = (n as f64).sqrt().ceil().max(1.0) as usize;
    (s, s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::meshing::edge_discret::MeshFace;
    use crate::meshing::model_builder::ModelBuilder;
    use crate::meshing::parameters::MeshParameters;
    use crate::primitives::BRepPrimBox;
    use crate::tgeometry::GeometryRegistry;

    fn box_face_polygon() -> MeshFace {
        MeshFace {
            outer_wire: vec![
                GpPnt2d::new(0.0, 0.0),
                GpPnt2d::new(1.0, 0.0),
                GpPnt2d::new(1.0, 1.0),
                GpPnt2d::new(0.0, 1.0),
                GpPnt2d::new(0.0, 0.0),
            ],
            inner_wires: vec![],
            deflection: 0.1,
            id: 0,
        }
    }

    #[test]
    fn extract_box_face_uv_count() {
        let tool = MeshTool::new(MeshParameters::default());
        let data = tool.extract_face(&box_face_polygon(), 1e-6).expect("extract");
        // Four corners after dedup, four boundary edges.
        assert_eq!(data.vertices_nb(), 4, "box face has 4 unique corners");
        assert_eq!(data.edges_nb(), 4, "box face has 4 boundary edges");
        // Edges form the expected closed loop (orientation-insensitive compare).
        let mut sorted: Vec<(usize, usize)> = data
            .constraint_edges
            .iter()
            .map(|&(a, b)| (a.min(b), a.max(b)))
            .collect();
        sorted.sort_unstable();
        assert_eq!(
            sorted,
            vec![(0, 1), (0, 3), (1, 2), (2, 3)],
            "loop connects every corner to its neighbours"
        );
    }

    #[test]
    fn extract_face_dedups_with_tolerance_and_keeps_holes() {
        let tool = MeshTool::new(MeshParameters::default());
        let mut face = box_face_polygon();
        // A hole in the middle of the face.
        face.inner_wires.push(vec![
            GpPnt2d::new(0.25, 0.25),
            GpPnt2d::new(0.75, 0.25),
            GpPnt2d::new(0.75, 0.75),
            GpPnt2d::new(0.25, 0.75),
            GpPnt2d::new(0.25, 0.25),
        ]);
        let data = tool.extract_face(&face, 1e-6).expect("extract");
        assert_eq!(data.vertices_nb(), 8, "4 outer + 4 inner corners");
        assert_eq!(data.edges_nb(), 8, "4 outer + 4 inner boundary edges");
        // Distinct corners stay distinct under a tight tolerance.
        let tight = tool.extract_face(&face, 0.0).expect("extract");
        assert_eq!(tight.vertices_nb(), 8);
    }

    #[test]
    fn extract_model_face_unpopulated_errors() {
        let shape = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone();
        let params = MeshParameters::default();
        let model = ModelBuilder::build_model(&shape, &params).expect("model built");
        // Pcurves are empty until the edge discretizers run.
        let tool = MeshTool::new(params);
        let err = tool.extract_model_face(&model, 0, 1e-6);
        assert!(err.is_err(), "unpopulated pcurves yield no UV points");
        GeometryRegistry::global().clear_shape(&shape);
    }

    #[test]
    fn compute_edge_deflection_prefers_stored_value() {
        let b = TopoBuilder::new();
        let edge = b.make_edge_segment(
            &occt_core::gp::GpPnt::new(0.0, 0.0, 0.0),
            &occt_core::gp::GpPnt::new(1.0, 0.0, 0.0),
        );
        let mut mesh_edge = data_model::MeshEdge::new(edge);
        mesh_edge.set_deflection(0.25);

        let tool = MeshTool::new(MeshParameters::default());
        assert_eq!(tool.compute_edge_deflection(&mesh_edge, 2.0), 0.25);

        // Unset deflection falls back to the params-based absolute value.
        mesh_edge.set_deflection(f64::MAX);
        let abs = tool.compute_edge_deflection(&mesh_edge, 2.0);
        assert!(abs > 0.0 && abs <= 0.002, "params deflection 0.001 * coeff, got {abs}");
    }

    #[test]
    fn node_classifier_splits_sides_of_constraint() {
        // Horizontal constraint (0,0)→(1,0): "above" is the positive-V side.
        let nc = NodeClassifier::new(0.0, 0.0, 1.0, 0.0);
        assert!(nc.is_above(0.5, 0.5));
        assert!(!nc.is_above(0.5, -0.5));
        assert!(!nc.is_above(0.5, 0.0), "on the constraint is not above");

        // OCCT's `mySign = Direction().X() > 0` convention keeps the "above"
        // set identical when the horizontal constraint is reversed.
        let nc_rev = NodeClassifier::new(1.0, 0.0, 0.0, 0.0);
        assert!(nc_rev.is_above(0.5, 0.5));

        // A vertical constraint separates by the U side and reverses with the
        // constraint direction.
        let nc_v = NodeClassifier::new(0.0, 0.0, 0.0, 1.0);
        assert!(nc_v.is_above(0.5, 0.5));
        assert!(!nc_v.is_above(-0.5, 0.5));
        let nc_v_rev = NodeClassifier::new(0.0, 1.0, 0.0, 0.0);
        assert!(!nc_v_rev.is_above(0.5, 0.5));

        // from_vertices builds the same classifier from UV vertices.
        let a = DelaunVertex::new(0.0, 0.0);
        let b = DelaunVertex::new(1.0, 0.0);
        let nc2 = NodeClassifier::from_vertices(&a, &b);
        assert!(nc2.is_above(0.5, 0.5));
    }
}
