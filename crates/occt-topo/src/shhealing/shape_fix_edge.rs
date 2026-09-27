//! `ShapeFix_Edge` subset needed by `ShapeFix_ComposeShell::DispatchWires`:
//! `FixAddCurve3d` (`ShapeFix_Edge.cxx:618-638`).

use std::sync::Arc;

use occt_geom::Curve;

use crate::brep_tool::BRepTool;
use crate::boptools_2d;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, Face};
use crate::tgeometry::{EdgeGeom, GeometryRegistry};

/// `ShapeFix_Edge` (`ShapeFix_Edge.hxx`).
pub struct ShapeFixEdge;

impl ShapeFixEdge {
    /// `FixAddCurve3d(edge)` (`cxx:618-638`).
    ///
    /// UNPORTED: `TempSameRange` (`cxx:335-...`, the `!BRep_Tool::SameRange`
    /// branch at `cxx:626-629`) is not applied.
    pub fn fix_add_curve3d(&self, edge: &Edge, face: &Face) -> bool {
        let reg = GeometryRegistry::global();
        // cxx:622.
        if BRepTool::is_degenerated(edge) || reg.edge_curve(&edge.0).is_some() {
            return false;
        }
        // cxx:631-635: ShapeBuild_Edge::BuildCurve3d(edge). The port represents
        // the 3D image of a pcurve as `CurveOnSurface`, the adapter BRepAdaptor_Curve
        // resolves for a curve-less edge.
        let (c2d, f, l) = match boptools_2d::curve_on_surface_range(edge, face) {
            Some(v) => v,
            None => return false,
        };
        let surface = match BRepTool::face_surface(face) {
            Some(s) => s,
            None => return false,
        };
        let c3d: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(c2d, surface, f, l));
        let mut g = EdgeGeom::new(c3d, f, l);
        if let Some(o) = reg.edge_geom(&edge.0) {
            g.tolerance = o.tolerance;
            g.same_parameter = o.same_parameter;
            g.same_range = o.same_range;
            g.degenerated = o.degenerated;
        }
        reg.set_edge(&edge.0, g);
        true
    }
}
