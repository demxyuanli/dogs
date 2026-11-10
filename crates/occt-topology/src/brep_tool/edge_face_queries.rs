use occt_core::gp::GpPnt;

use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

use super::BRepTool;

impl BRepTool {
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
}
