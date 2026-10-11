//! `ShapeAlgo_AlgoContainer::ConnectNextWire`
//! (`ShapeAlgo_AlgoContainer.cxx:61-113`) and the status bits
//! `ShapeAnalysis_Wire::CheckShapeConnect` writes
//! (`ShapeAnalysis_Wire.cxx:2111-2202`).

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::REAL_SMALL;
use occt_geom::curve::Curve;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Edge;
use crate::shape_fix_compose_shell::reverse_wire_data;
use crate::shhealing::{first_vertex, last_vertex};
use crate::topo_tools_full::is_same;

const DONE1: i32 = 0x0001;
const DONE2: i32 = 0x0002;
const DONE3: i32 = 0x0004;
const DONE4: i32 = 0x0008;
const DONE5: i32 = 0x0010;
const DONE6: i32 = 0x0020;
const FAIL1: i32 = 0x0100;
const FAIL2: i32 = 0x0200;
const FAIL: i32 = 0xff00;

fn has(bits: i32, mask: i32) -> bool {
    bits & mask != 0
}

/// First vertex of the first edge and last vertex of the last edge
/// (`ShapeAnalysis::FindBounds` on the wire `CheckShapeConnect` receives).
fn boundary_ends(edges: &[Edge]) -> Option<(GpPnt, GpPnt)> {
    let first = edges.first()?;
    let last = edges.last()?;
    let head = first_vertex(first)?;
    let tail = last_vertex(last)?;
    Some((BRepTool::vertex_point(&head), BRepTool::vertex_point(&tail)))
}

struct ConnectCheck {
    status: i32,
    tailhead: f64,
    tailtail: f64,
}

/// `CheckShapeConnect` (`cxx:2111-2202`). `prec` is the `maxtol` argument.
fn check_shape_connect(prev: &[Edge], next: &[Edge], precision: f64, prec: f64) -> ConnectCheck {
    let (Some((pf, pl)), Some((p1, p2))) = (boundary_ends(prev), boundary_ends(next)) else {
        return ConnectCheck { status: FAIL1, tailhead: 0.0, tailtail: 0.0 };
    };
    let tailhead = p1.distance(&pl);
    let tailtail = p2.distance(&pl);
    let headhead = p1.distance(&pf);
    let headtail = p2.distance(&pf);
    let mut dm1 = tailhead;
    let mut dm2 = headtail;
    let mut res1 = 0;
    let mut res2 = 0;
    if tailhead > tailtail {
        res1 = 1;
        dm1 = tailtail;
    }
    if headtail > headhead {
        res2 = 1;
        dm2 = headhead;
    }
    let mut result = res1;
    let min3d = dm1.min(dm2);
    if dm1 > dm2 {
        result = res2 + 2;
    }
    let mut status = match result {
        1 => DONE2,
        2 => DONE3,
        3 => DONE4,
        _ => DONE1,
    };
    if res1 == 0 {
        status |= DONE5;
    }
    if res2 == 0 {
        status |= DONE6;
    }
    if min3d > precision.max(prec) {
        status = FAIL2;
    }
    ConnectCheck { status, tailhead, tailtail }
}

/// Reverse `next` (and sometimes `sewd`) so the new segment meets the
/// accumulated wire, then append clones of it. An empty `next` succeeds
/// (`cxx:70-72`). An empty `sewd` takes `next` unchanged (`cxx:77-80`).
///
/// Returns false on `ShapeExtend_FAIL`. `IGESControl_IGESBoundary::Transfer`
/// then rebuilds both wires through `ShapeFix_Wire::Perform` and keeps the one
/// with the smaller max tolerance (`IGESControl_IGESBoundary.cxx:441-508`),
/// which the port does in `iges_to_brep_topo_curve.rs` (the
/// `if !connect_next_wire(...)` arm of `transfer_boundary`).
pub(in crate::iges_read) fn connect_next_wire(
    sewd: &mut Vec<Edge>,
    next: &mut Vec<Edge>,
    precision: f64,
    maxtol: f64,
) -> bool {
    if next.is_empty() {
        return true;
    }
    if sewd.is_empty() {
        sewd.extend(next.iter().cloned());
        return true;
    }
    let check = check_shape_connect(sewd, next, precision, maxtol);
    if check.tailhead > precision
        && check.tailtail > precision
        && (has(check.status, DONE4) || has(check.status, DONE3))
    {
        reverse_wire_data(sewd); // cxx:91
        if has(check.status, DONE3) {
            reverse_wire_data(next); // cxx:96
        }
    } else if !has(check.status, FAIL) && !has(check.status, DONE5) {
        reverse_wire_data(next); // cxx:103-105
    }
    let ok = !has(check.status, FAIL);
    if ok {
        sewd.extend(next.iter().cloned()); // cxx:110
    }
    ok
}

/// `IGESControl_IGESBoundary::Check` closure test (`cxx:69-89`).
///
/// `CheckConnected(1)` (`ShapeAnalysis_Wire.cxx:693-760`) with precision
/// `GetMaxTol` and the default `prec` 0. `LastCheckStatus(FAIL)` is true for
/// both `FAIL1` and `FAIL2`, so a gap above `maxtol` fails whether or not the
/// reversed end is close. The same vertex, or a gap within `maxtol`, passes.
pub(in crate::iges_read) fn closure_failed(edges: &[Edge], maxtol: f64) -> bool {
    if edges.is_empty() {
        return false; // cxx:695-697 status stays OK
    }
    let e1 = &edges[edges.len() - 1];
    let e2 = &edges[0];
    let (Some(v1), Some(v2)) = (last_vertex(e1), first_vertex(e2)) else {
        return true; // cxx:712-714 FAIL2
    };
    if is_same(&v1.0, &v2.0) {
        return false; // cxx:716-718
    }
    let dist = BRepTool::vertex_point(&v1).distance(&BRepTool::vertex_point(&v2));
    if dist <= REAL_SMALL || dist <= maxtol {
        return false; // cxx:724-731 DONE1 / DONE2. prec is 0, so DONE3 does not pass a positive gap.
    }
    true
}

/// `IGESToBRep_IGESBoundary::ReverseCurves3d` (`IGESToBRep_IGESBoundary.cxx:348-384`).
/// Edge order is reversed and each 3D curve is reversed. The new edge keeps
/// the orientation from before `ShapeExtend_WireData::Reverse`, which the
/// following `Reverse` cancels.
pub(in crate::iges_read) fn reverse_curves_3d(edges: &[Edge], precision: f64, maxtol: f64) -> Vec<Edge> {
    let b = TopoBuilder::new();
    let mut built = Vec::new();
    for edge in edges.iter().rev() {
        let Some(curve) = BRepTool::edge_curve(edge) else {
            built.push(edge.clone());
            continue;
        };
        let (p1, p2) = BRepTool::edge_parameters(edge);
        let rev: Arc<dyn Curve> = Arc::from(curve.reversed());
        let (first, last) = if curve.is_periodic() {
            (curve.reversed_parameter(p2), curve.reversed_parameter(p1))
        } else {
            (
                curve
                    .reversed_parameter(curve.last_parameter())
                    .max(curve.reversed_parameter(p2)),
                curve
                    .reversed_parameter(curve.first_parameter())
                    .min(curve.reversed_parameter(p1)),
            )
        };
        let mut neu = b.make_edge(rev.clone(), first, last);
        let v1 = b.make_vertex(rev.d0(first), 0.0);
        let v2 = b.make_vertex(rev.d0(last), 0.0);
        b.add_edge_vertices(&mut neu, &v1, &v2);
        neu.0.set_orientation(edge.0.orientation());
        if !edge.0.location().is_identity() {
            neu.0.move_location(&edge.0.location());
        }
        built.push(neu);
    }
    let mut wire = b.make_wire(&built);
    crate::shhealing::fix_connected_closed(&mut wire, precision, maxtol, true);
    crate::topo_tools_full::edges_of_wire(&wire)
}
