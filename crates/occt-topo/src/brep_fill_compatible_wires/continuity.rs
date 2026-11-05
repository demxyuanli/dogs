//! The file-static `WireContinuity` of `BRepFill_CompatibleWires.cxx`.
//!
//! Source: `BRepFill_CompatibleWires.cxx:300-364`.

use occt_core::precision::ANGULAR;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::brep_lprop::{self, CN, C0};
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{Edge, Wire};

/// `static void WireContinuity(const TopoDS_Wire& W, GeomAbs_Shape& contW)`
/// (`cxx:300-364`): the minimum continuity over the wire's consecutive edge
/// pairs, closing the loop from the last edge back to the first.
///
/// A wire with a degenerated edge answers `GeomAbs_CN` (`cxx:307-320`): the
/// pair loop is skipped entirely.
///
/// The pairs are evaluated with `BRepLProp::Continuity(Curve1, Curve2, U1, U2,
/// Eps, Precision::Angular())`, whose `BRepAdaptor_Curve` arguments carry no
/// orientation, so both sides are evaluated FORWARD here.
pub fn wire_continuity(w: &Wire) -> u8 {
    let mut exp = WireExplorer::new();
    exp.init(w, None);
    let mut edges: Vec<Edge> = Vec::new();
    let mut is_degenerated = false;
    while exp.more() {
        if let Some(e) = exp.current_edge() {
            edges.push(e.clone());
            // `if (BRep_Tool::Degenerated(anExp.Current()))`
            if BRepTool::is_degenerated(&e) {
                is_degenerated = true;
            }
        }
        exp.next();
    }
    if is_degenerated {
        return CN;
    }

    let nb_edges = edges.len();
    let mut cont_w = CN;
    let mut testconti = true;
    for j in 1..=nb_edges {
        // `cxx:329-339`: the closing pair is (last, first).
        let (edge1, edge2) = if j == nb_edges {
            (edges[nb_edges - 1].clone(), edges[0].clone())
        } else {
            (edges[j - 1].clone(), edges[j].clone())
        };

        // `TopExp::Vertices(Edge1, Vbid, V1, true)` / `(Edge2, V2, Vbid, true)`
        // (`cxx:342-343`): with `CumOri = true`, V1 is the *last* vertex of
        // Edge1 and V2 the *first* of Edge2.
        let (e1_first, e1_last) = crate::brep_tools_wire_explorer::edge_vertices_cum_ori(&edge1);
        let (e2_first, e2_last) = crate::brep_tools_wire_explorer::edge_vertices_cum_ori(&edge2);
        let (Some(v1), Some(v2)) = (e1_last.or(e1_first), e2_first.or(e2_last)) else {
            // `BRep_Tool::Parameter` on a null vertex is undefined in OCCT.
            continue;
        };
        let u1 = BRepTool::parameter_of_vertex(&edge1.0, &v1.0);
        let u2 = BRepTool::parameter_of_vertex(&edge2.0, &v2.0);
        let Some(c1) = BRepTool::edge_curve_world(&edge1) else {
            continue;
        };
        let Some(c2) = BRepTool::edge_curve_world(&edge2) else {
            continue;
        };
        // `double Eps = BRep_Tool::Tolerance(V2) + BRep_Tool::Tolerance(V1)`
        // (`cxx:349`).
        let eps = crate::brep_tool::BRepTool::vertex_tolerance(&v2)
            + crate::brep_tool::BRepTool::vertex_tolerance(&v1);
        // `if (j == nbEdges) testconti = Curve1.Value(U1).IsEqual(
        //  Curve2.Value(U2), Eps)` (`cxx:351-354`).
        if j == nb_edges {
            testconti = c1.value(u1).is_equal_with_tolerance(&c2.value(u2), eps);
        }

        if testconti {
            // `cont = BRepLProp::Continuity(...)` (`cxx:358`).
            let cont = brep_lprop::continuity(
                c1.as_ref(),
                Orientation::Forward,
                u1,
                (c1.first_parameter(), c1.last_parameter()),
                c2.as_ref(),
                Orientation::Forward,
                u2,
                (c2.first_parameter(), c2.last_parameter()),
                false,
                eps,
                ANGULAR,
            )
            .unwrap_or(C0);
            if cont <= cont_w {
                cont_w = cont;
            }
        }
    }
    cont_w
}
