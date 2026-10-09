use super::*;

/// `ShapeFix_Wire::FixReorder` (`cxx:487-534`) then `FixReorder(sawo)`
/// (`cxx:1351-1399`). `Perform` calls `FixReorder()` with `theModeBoth=false`,
/// so CheckOrder is 3D (`cxx:504`) for every surface type. Status 0 leaves the
/// wire unchanged.
///
/// Returns the `ReorderOK` flag of `Perform` (`cxx:317-325`): false only when
/// the checked order could not be applied (`ShapeExtend_FAIL`, which then
/// disables FixShifted at `cxx:368-374`). A wire that is already ordered
/// returns true (`FixReorder(sawo)` returns false for status 0 without setting
/// FAIL).
pub fn fix_reorder_wire(wire: &Wire, _face: &Face) -> bool {
    let (ok, _status) = fix_reorder_wire_3d(wire);
    ok
}

/// The stored FORWARD/REVERSED edges of `wire` (`ShapeExtend_WireData::Init`).
pub(in crate::shhealing) fn stored_manifold_edges(wire: &Wire) -> Vec<Edge> {
    edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect()
}

/// `ShapeAnalysis_Wire::CheckOrder(..., isClosed, mode3d=true)` -> WireOrder
/// (`ShapeFix_Wire.cxx:496-505`). `None` is the FAIL2 short-circuit of
/// `ShapeAnalysis_Wire.cxx:617-621` (a missing edge vertex), for which the
/// caller sees ReorderOK = true and status 0.
pub(in crate::shhealing) fn build_wire_order_3d(wire: &Wire) -> Option<(WireOrder, Vec<Edge>)> {
    let stored = stored_manifold_edges(wire);
    if stored.len() < 2 {
        return None;
    }
    let mut order = WireOrder::new();
    for e in &stored {
        let (Some(v1), Some(v2)) = (first_vertex(e), last_vertex(e)) else {
            return None;
        };
        order.add_edge_xyz(BRepTool::vertex_point(&v1), BRepTool::vertex_point(&v2));
    }
    order.perform();
    Some((order, stored))
}

/// `ShapeFix_Wire::FixReorder()` (`ShapeFix_Wire.cxx:487-534`): reorder the
/// wire in 3D. Returns (ReorderOK, `ShapeAnalysis_WireOrder::Status`);
/// `Reversed` is the `myStatusReorder` DONE3 case of `cxx:524-527`.
pub fn fix_reorder_wire_3d(wire: &Wire) -> (bool, WireOrderStatus) {
    let Some((order, stored)) = build_wire_order_3d(wire) else {
        return (true, WireOrderStatus::Same);
    };
    let status = order.status();
    let ok = if status == WireOrderStatus::Same {
        true
    } else if order.nb_edges() != stored.len() {
        false
    } else {
        apply_wire_order(wire, &stored, &order)
    };
    (ok, status)
}

/// `ShapeFix_Wire::FixReorder(sawo)` (`ShapeFix_Wire.cxx:1351-1400`): apply a
/// prebuilt order. `true` is the DONE1 of `cxx:1398`; the status 0 guard
/// (`cxx:1360-1363`) and the FAIL1/2/3 guards (`cxx:1364-1385`) return false.
pub fn fix_reorder_wire_with_order(wire: &Wire, order: &WireOrder) -> bool {
    let stored = stored_manifold_edges(wire);
    if stored.len() != order.nb_edges() {
        return false; // cxx:1372-1376 FAIL2
    }
    for i in 1..=stored.len() {
        if order.ordered(i) == 0 {
            return false; // cxx:1380-1385 FAIL3
        }
    }
    apply_wire_order(wire, &stored, order)
}

/// The `ShapeFix_Wire::FixReorder` tail (`cxx:1387-1399` / `cxx:1431-1474`):
/// build the reordered edge list from the signed order and store it back,
/// keeping the non-edge children and the port's wire-orientation convention.
pub(in crate::shhealing) fn apply_wire_order(wire: &Wire, stored: &[Edge], order: &WireOrder) -> bool {
    let mut new_edges = Vec::with_capacity(stored.len());
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        if signed == 0 {
            return false;
        }
        let mut e = stored[signed.unsigned_abs() as usize - 1].clone();
        if signed < 0 {
            e.0.reverse();
        }
        new_edges.push(e);
    }
    let extras: Vec<TopoShape> = {
        let ts = wire.0.tshape.read().expect("poisoned TShape lock");
        ts.children
            .iter()
            .filter(|c| {
                !c.is_edge()
                    || (c.orientation() != Orientation::Forward
                        && c.orientation() != Orientation::Reversed)
            })
            .cloned()
            .collect()
    };
    // `ShapeExtend_WireData` stores the edges with the wire orientation already
    // composed (`ShapeExtend_WireData::Init`, `ShapeExtend_WireData.cxx:80-120`)
    // and `ShapeExtend_WireData::Wire()` rebuilds a fresh FORWARD wire from that
    // list (`ShapeExtend_WireData.cxx:651-685`). Our children stay under the
    // original wire, so undo the wire's own orientation before storing, keeping
    // the composed view identical to OCCT's WireData.
    let wire_reversed = wire.0.orientation() == Orientation::Reversed;
    {
        let mut ts = wire.0.tshape.write().expect("poisoned TShape lock");
        ts.children.clear();
        for e in &new_edges {
            let mut child = e.0.clone();
            if wire_reversed {
                child.reverse();
            }
            ts.children.push(child);
        }
        ts.children.extend(extras);
    }
    true
}

/// `ShapeAnalysis_Wire::CheckConnected(num, prec)` status
/// (`ShapeAnalysis_Wire.cxx:693-762`). The three `ShapeExtend_DONE` cases are
/// the only ones `FixConnected` acts on; `FAIL1` / `FAIL2` both return false
/// before any fix (`cxx:1487-1493`).
#[derive(Clone, Copy, PartialEq, Eq)]
pub(in crate::shhealing) enum ConnectedCheck {
    None,
    /// `cxx:727`: `myMin3d <= gp::Resolution()`, "absolutely confused".
    Done1,
    /// `cxx:730`: `myMin3d <= myPrecision`.
    Done2,
    /// `cxx:733`: `myMin3d <= prec`.
    Done3,
}

/// `ShapeAnalysis_Wire::CheckConnected(num, prec)`
/// (`ShapeAnalysis_Wire.cxx:693-762`). `my_precision` is the analyzer's
/// `myPrecision` (`ShapeAnalysis_Wire.cxx:203`), the value
/// `ShapeFix_Wire::Perform` was driven with.
pub(in crate::shhealing) fn check_connected(wire: &Wire, my_precision: f64, prec: f64, num: usize) -> ConnectedCheck {
    let nb = wire_edges_nb(wire);
    if nb < 1 {
        return ConnectedCheck::None;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let (Some(v1), Some(v2)) = (last_vertex(&edges[n1 - 1]), first_vertex(&edges[n2 - 1])) else {
        return ConnectedCheck::None;
    };
    if is_same(&v1.0, &v2.0) {
        return ConnectedCheck::None;
    }
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    let mut min3d = p1.distance(&p2);
    if min3d <= occt_core::precision::REAL_SMALL {
        return ConnectedCheck::Done1;
    }
    if min3d <= my_precision {
        return ConnectedCheck::Done2;
    }
    if min3d <= prec {
        return ConnectedCheck::Done3;
    }
    // `cxx:737-758`: the reverse test on the last vertex of E2.
    if n1 == n2 {
        return ConnectedCheck::None;
    }
    let Some(v2b) = last_vertex(&edges[n2 - 1]) else {
        return ConnectedCheck::None;
    };
    let p2b = BRepTool::vertex_point(&v2b);
    let dist = p1.distance(&p2b);
    if dist > my_precision {
        return ConnectedCheck::None;
    }
    min3d = dist;
    let _ = min3d;
    ConnectedCheck::None
}

/// `ShapeFix_Wire::FixConnected(num, prec, theUpdateWire)`
/// (`ShapeFix_Wire.cxx:1476-1608`), at read time (`theUpdateWire = false`,
/// `cxx:569`). `Context()` is not modelled: the two `CopyReplaceVertices`
/// writes to `WireData` (`cxx:1585-1601`) are applied directly to the wire's
/// child list, the same convention `fix_degenerated` uses.
pub(in crate::shhealing) fn fix_connected_one(wire: &mut Wire, my_precision: f64, prec: f64, num: usize) -> bool {
    let check = check_connected(wire, my_precision, prec, num);
    if check == ConnectedCheck::None {
        return false;
    }
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let v1 = last_vertex(&e1);
    let v2 = first_vertex(&e2);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return false;
    };
    // `cxx:1505-1542`: the merged vertex.
    let combined: Vertex = if check == ConnectedCheck::Done1 {
        let v2_is_last = last_vertex(&e2)
            .map(|l| is_same(&v2.0, &l.0))
            .unwrap_or(false);
        if v2_is_last {
            v2.clone()
        } else {
            v1.clone()
        }
    } else {
        combine_vertex(&v1, &v2, 1.0001)
    };

    // `cxx:1544-1601`: this port's edges are frozen by `TopoBuilder::add`, so
    // the `E1.Free() && E2.Free() && myTopoMode` arm never runs; the
    // `CopyReplaceVertices` arm does.
    if nb < 2 {
        // `cxx:1559-1565`: `sbe.CopyReplaceVertices(E2, V, V)` + `Set(tmpE, n2)`.
        let tmp_e = copy_replace_vertices_with(&e2, Some(&combined), Some(&combined));
        wire_set_edge_composed(wire, n2, &tmp_e);
    } else {
        // `cxx:1584-1599`: `CopyReplaceVertices(E2, V, null)` + `Set(tmpE2, n2)`,
        // then the same on E1 unless DONE1 with a non-closed E2.
        let tmp_e2 = copy_replace_vertices_with(&e2, Some(&combined), None);
        wire_set_edge_composed(wire, n2, &tmp_e2);
        let first_is_last = first_vertex(&e2)
            .zip(last_vertex(&e2))
            .map(|(f, l)| is_same(&f.0, &l.0))
            .unwrap_or(false);
        if check != ConnectedCheck::Done1 || first_is_last {
            let tmp_e1 = copy_replace_vertices_with(&e1, None, Some(&combined));
            wire_set_edge_composed(wire, n1, &tmp_e1);
        }
    }
    true
}

/// `ShapeFix_Wire::FixConnected(prec)` (`ShapeFix_Wire.cxx:557-596`), the API
/// driver `Perform` calls at `cxx:353-357` under `myFixConnectedMode`
/// (`FromSTEP.FixShape.FixConnectedMode`, `STEPControl_Controller.cxx:235`,
/// default `-1` -> `NeedFix(-1, ReorderOK)` = ReorderOK). `myClosedMode` is
/// true (`ShapeFix_Wire.cxx:842`), so `aStop = 0` and every junction is
/// visited, from `NbEdges()` down to 1. `Perform` passes the default
/// `prec = -1.`, which `cxx:1485` turns into `MaxTolerance()`.
/// `ShapeFix_Shape::SetMaxTolerance` propagates through
/// `ShapeFix_Solid`/`ShapeFix_Shell`/`ShapeFix_Face` to the wire tool
/// (`ShapeFix_Solid.cxx:745-749`, `ShapeFix_Shell.cxx:1709-1713`,
/// `ShapeFix_Face.cxx:173-177`), so this is `FromSTEP.FixShape.MaxTolerance3d`
/// (`STEPControl_Controller.cxx:204-206`, set from
/// `STEPControl_ActorRead::myMaxTol` at `cxx:1146`).
pub(in crate::shhealing) fn fix_connected_all(wire: &mut Wire, my_precision: f64, max_tolerance: f64) -> bool {
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    let mut done = false;
    for a_i in (1..=nb).rev() {
        done |= fix_connected_one(wire, my_precision, max_tolerance, a_i);
    }
    done
}
