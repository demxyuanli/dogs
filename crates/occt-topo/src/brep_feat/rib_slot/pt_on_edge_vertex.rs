//! `BRepFeat_RibSlot::PtOnEdgeVertex`
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:1327-1411`).
//!
//! Answers whether a point of the profile lies on an edge of the basis, and
//! whether it also lies on one of that edge's vertices.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::curve::Curve;
use occt_geom::extrema_pc::point_curve_extrema_all;
use occt_geom::trimmed::GeomTrimmedCurveBasis;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::locope::shape_maps::ShapeSet;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::shhealing::{first_vertex, last_vertex};
use crate::topo_tools::{explorer_shapes, vertex_position};

/// The four outputs of `PtOnEdgeVertex` (`cxx:1334-1337`).
#[derive(Clone, Debug, Default)]
pub struct PtOnEdgeVertex {
    /// `PtOnEdge`.
    pub pt_on_edge: bool,
    /// `OnEdge`.
    pub on_edge: Option<Edge>,
    /// `PtOnVertex`.
    pub pt_on_vertex: bool,
    /// `OnVertex`.
    pub on_vertex: Option<Vertex>,
}

/// `PtOnEdgeVertex(RevolRib, shape, point, FirstVertex, LastVertex, PtOnEdge,
/// OnEdge, PtOnVertex, OnVertex)` (`cxx:1327-1411`).
///
/// The two vertex arguments are unused in the reference as well: both are
/// commented out of the signature body (`cxx:1330-1331`).
///
/// The projection takes the `GeomAPI_ProjectPointOnCurve` route of
/// `cxx:1373-1391`: for a non-revolution rib the curve is trimmed to the edge
/// range and exactly one extremum is required, for a revolution rib the
/// untrimmed curve answers with at least one.
pub fn pt_on_edge_vertex(
    revol_rib: bool,
    shape: &TopoShape,
    point: &GpPnt,
) -> PtOnEdgeVertex {
    let mut out = PtOnEdgeVertex::default();
    // `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> Map`
    // (`cxx:1351`): `IsSame` keys, so a shared edge is visited once.
    let mut map = ShapeSet::new();
    for sh in explorer_shapes(shape, ShapeType::Edge) {
        let Some(e) = Edge::wrap(sh.clone()) else {
            continue;
        };
        if !map.add(&sh) {
            continue;
        }
        if !revol_rib && BRepTool::is_degenerated(&e) {
            continue;
        }
        let (fff, lll) = BRepTool::edge_parameters(&e);
        let Some(ccc) = BRepTool::edge_curve_world(&e) else {
            continue;
        };
        // `ccc = new Geom_TrimmedCurve(ccc, fff, lll)` for a non-revolution rib
        // (`cxx:1370`).
        let ccc: Arc<dyn Curve> = if revol_rib {
            ccc
        } else {
            Arc::new(GeomTrimmedCurveBasis::new(ccc, fff, lll))
        };
        let all = point_curve_extrema_all(ccc.as_ref(), point);
        // `proj.NbPoints() == 1` (`cxx:1378`) / `proj.NbPoints() >= 1`
        // (`cxx:1384`).
        let test_ok = if revol_rib {
            !all.is_empty()
        } else {
            all.len() == 1
        };
        if !test_ok {
            continue;
        }
        // `proj.Distance(1) <= BRep_Tool::Tolerance(e)` (`cxx:1389`).
        if all[0].distance > BRepTool::edge_tolerance(&e) {
            continue;
        }
        out.pt_on_edge = true;
        out.on_edge = Some(e.clone());
        let (Some(ev1), Some(ev2)) = (first_vertex(&e), last_vertex(&e)) else {
            break;
        };
        let ep1 = vertex_position(&ev1);
        let ep2 = vertex_position(&ev2);
        // `cxx:1397-1407`; the second test re-reads `Tolerance(ev1)` in the
        // reference.
        if point.distance(&ep1) <= BRepTool::vertex_tolerance(&ev1) {
            out.pt_on_vertex = true;
            out.on_vertex = Some(ev1);
        } else if point.distance(&ep2) <= BRepTool::vertex_tolerance(&ev1) {
            out.pt_on_vertex = true;
            out.on_vertex = Some(ev2);
        }
        break;
    }
    out
}
