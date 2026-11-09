//! `UpdateVertex` (`BRepOffset_Inter2d.cxx:1680-1721`): move a vertex onto a
//! rebuilt edge when it sits at one of the original edge's ends.
//!
//! The reference tests the vertex point against both ends of `OE`; when it
//! matches an end, that end's parameter is accepted only if it lies inside
//! `NE`'s range (within `Precision::PConfusion()`) and `NE` evaluates back to
//! within `TolConf` of the point, and the parameter is then written on `NE`.
//!
//! No caller in the ported part yet: `ConnexIntByInt` (`cxx:1830-2095`) uses it
//! when it rebinds an edge rebuilt by `ExtentEdge`.

use occt_core::precision::PCONFUSION;

use crate::abs::Orientation;
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Vertex};
use crate::topo_tools::vertex_position;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

/// `UpdateVertex(V, OE, NE, TolConf)` (`cxx:1680-1721`).
pub(crate) fn update_vertex(v: &Vertex, oe: &Edge, ne: &Edge, tol_conf: f64) -> bool {
    let oc = CurveAdaptor::of_edge(oe);
    let nc = CurveAdaptor::of_edge(ne);
    let of = oc.first_parameter();
    let ol = oc.last_parameter();
    let nf = nc.first_parameter();
    let nl = nc.last_parameter();
    let par_tol = PCONFUSION;
    let p = vertex_position(v);

    let mut u = 0.0;
    let mut ok = false;

    // `cxx:1694-1701`: the first end.
    if let Some(pof) = oc.d0(of) {
        if p.distance(&pof) < tol_conf
            && of >= nf + par_tol
            && of <= nl + par_tol
            && nc.d0(of).is_some_and(|pn| p.distance(&pn) < tol_conf)
        {
            ok = true;
            u = of;
        }
    }
    // `cxx:1702-1709`: the last end (it overwrites `U` when both match).
    if let Some(pol) = oc.d0(ol) {
        if p.distance(&pol) < tol_conf
            && ol >= nf + par_tol
            && ol <= nl + par_tol
            && nc.d0(ol).is_some_and(|pn| p.distance(&pn) < tol_conf)
        {
            ok = true;
            u = ol;
        }
    }

    if ok {
        // `TopoDS_Edge EE = TopoDS::Edge(NE.Oriented(FORWARD))` (`cxx:1711-1713`)
        // is built and never read; then
        // `B.UpdateVertex(V.Oriented(INTERNAL), U, NE, Tolerance(NE))`.
        let _ = ne.0.oriented(Orientation::Forward);
        let v_int = v.0.oriented(Orientation::Internal);
        let _ = update_vertex_on_edge(&v_int, u, &ne.0, BRepTool::edge_tolerance(ne));
    }
    ok
}
