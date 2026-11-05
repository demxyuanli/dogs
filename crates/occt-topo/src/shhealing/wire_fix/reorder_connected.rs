use super::*;
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};

/// `V.Oriented(anOrientation)` (`TopoDS_Shape::Oriented`): the same `TShape`
/// carried with another orientation, which is what `FixConnected` hands to
/// `Context()->Replace`.
fn oriented_like(v: &Vertex, orientation: Orientation) -> TopoShape {
    let mut s = v.0.clone();
    s.set_orientation(orientation);
    s
}

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

/// The FORWARD/REVERSED edge sequence built by `ShapeExtend_WireData::Init`
/// (`ShapeExtend_WireData.cxx:80-122`). `TopoDS_Iterator` composes each edge
/// orientation with the wire orientation; a REVERSED wire then prepends every
/// edge, reversing the resulting sequence.
pub(crate) fn stored_manifold_edges(wire: &Wire) -> Vec<Edge> {
    let mut edges: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if wire.0.orientation() == Orientation::Reversed {
        edges.reverse();
    }
    edges
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

/// `ShapeAnalysis_Wire::CheckOrder(sawo, closed)` over the 3D wire
/// (`ShapeFix_Wire.cxx:316-318`): `ShapeFix_Wire::Perform` reads it as
/// `ReorderOK = (CheckOrder(...) == 0)`, i.e. `ShapeAnalysis_Same`. The FAIL2
/// short-circuit (`ShapeAnalysis_Wire.cxx:617-621`, a missing edge vertex)
/// answers `true`, exactly like `build_wire_order_3d`'s `None`.
pub fn check_order_3d(wire: &Wire) -> bool {
    match build_wire_order_3d(wire) {
        Some((order, _stored)) => order.status() == WireOrderStatus::Same,
        None => true,
    }
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
    // original wire, so undo both the prepend and the composed orientation
    // before storing. The next `stored_manifold_edges` then reproduces this
    // WireData sequence.
    let wire_reversed = wire.0.orientation() == Orientation::Reversed;
    {
        let mut ts = wire.0.tshape.write().expect("poisoned TShape lock");
        ts.children.clear();
        let stored_edges: Vec<Edge> = if wire_reversed {
            new_edges.iter().rev().cloned().collect()
        } else {
            new_edges.clone()
        };
        for e in &stored_edges {
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
/// `ShapeExtend_WireData` index (1-based, the index `FixConnected` computes)
/// to the child slot `wire_set_edge` addresses. A REVERSED wire's WireData
/// list is the reverse of storage (`ShapeExtend_WireData.cxx:115-118`).
fn wiredata_storage_index(wire: &Wire, wiredata_index: usize, nb: usize) -> usize {
    if wire.0.orientation() == Orientation::Reversed {
        nb - wiredata_index + 1
    } else {
        wiredata_index
    }
}

pub(in crate::shhealing) fn check_connected(wire: &Wire, precision_par: f64, prec: f64, num: usize) -> ConnectedCheck {
    let nb = wire_edges_nb(wire);
    if nb < 1 {
        return ConnectedCheck::None;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    // `ShapeAnalysis_Wire::CheckConnected` reads `WireData()->Edge`
    // (`ShapeAnalysis_Wire.cxx:705-706`), the prepended sequence of a
    // REVERSED wire, not raw storage order.
    let edges = stored_manifold_edges(wire);
    if n1 == 0 || n2 == 0 || n1 > edges.len() || n2 > edges.len() {
        return ConnectedCheck::None;
    }
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
    if min3d <= precision_par {
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
    if dist > precision_par {
        return ConnectedCheck::None;
    }
    min3d = dist;
    let _ = min3d;
    ConnectedCheck::None
}

/// `ShapeFix_Wire::FixConnected(num, prec, theUpdateWire)`
/// (`ShapeFix_Wire.cxx:1476-1608`), at read time (`theUpdateWire = false`,
/// `cxx:569`).
///
/// `ctx` is `Context()`: when it is present, every substitution the `.cxx`
/// records goes through it - `Context()->Replace(V1, V.Oriented(...))` /
/// `Replace(V2, ...)` (`cxx:1513`, `:1525`, `:1537-1541`) and
/// `Context()->Replace(E1/E2, tmpE)` (`cxx:1561`, `:1587`, `:1594`). The
/// substitution is keyed by `TShape` identity (`TopTools_ShapeMapHasher`,
/// `IsSame`), so `UpdateWire` (`cxx:4293-4310`) later rewrites *every*
/// occurrence of the old vertex in the wire, not only the junction's two
/// endpoints; that is what keeps a closed edge closed when it shares its
/// vertex with a junction. `None` is the null-context tool (the naked
/// `ShapeFix_Wire` of `IGESControl_IGESBoundary::Transfer`), where the `.cxx`
/// skips both the `Replace` calls and `UpdateWire`.
pub(in crate::shhealing) fn fix_connected_one(
    wire: &mut Wire,
    ctx: Option<&mut SharedReShape>,
    precision: f64,
    prec: f64,
    num: usize,
) -> bool {
    let check = check_connected(wire, precision, prec, num);
    if check == ConnectedCheck::None {
        return false;
    }
    let mut ctx = ctx;
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    // Same WireData sequence as `CheckConnected` (`ShapeFix_Wire.cxx:1496-1500`).
    let edges = stored_manifold_edges(wire);
    if n1 == 0 || n2 == 0 || n1 > edges.len() || n2 > edges.len() {
        return false;
    }
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let v1 = last_vertex(&e1);
    let v2 = first_vertex(&e2);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return false;
    };
    // `cxx:1505-1542`: the merged vertex, plus the `Context()->Replace` pair
    // that makes the substitution global.
    let combined: Vertex = if check == ConnectedCheck::Done1 {
        let v2_is_last = last_vertex(&e2)
            .map(|l| is_same(&v2.0, &l.0))
            .unwrap_or(false);
        let (keep, drop) = if v2_is_last {
            (v2.clone(), v1.clone())
        } else {
            (v1.clone(), v2.clone())
        };
        if let Some(c) = ctx.as_deref_mut() {
            let repl = oriented_like(&keep, drop.0.orientation());
            c.replace(&drop.0, &repl);
        }
        keep
    } else {
        let merged = combine_vertex(&v1, &v2, 1.0001);
        if let Some(c) = ctx.as_deref_mut() {
            for src in [&v1, &v2] {
                let repl = oriented_like(&merged, src.0.orientation());
                c.replace(&src.0, &repl);
            }
        }
        merged
    };

    // `cxx:1544-1601`: this port's edges are frozen by `TopoBuilder::add`, so
    // the `E1.Free() && E2.Free() && myTopoMode` arm never runs; the
    // `CopyReplaceVertices` arm does.
    let s1 = wiredata_storage_index(wire, n1, nb);
    let s2 = wiredata_storage_index(wire, n2, nb);
    if nb < 2 {
        // `cxx:1559-1565`: `sbe.CopyReplaceVertices(E2, V, V)` + `Set(tmpE, n2)`.
        let tmp_e = copy_replace_vertices_with(&e2, Some(&combined), Some(&combined));
        if let Some(c) = ctx.as_deref_mut() {
            c.replace(&e2.0, &tmp_e.0);
        }
        wire_set_edge_composed(wire, s2, &tmp_e);
    } else {
        // `cxx:1584-1599`: `CopyReplaceVertices(E2, V, null)` + `Set(tmpE2, n2)`,
        // then the same on E1 unless DONE1 with a non-closed E2.
        let tmp_e2 = copy_replace_vertices_with(&e2, Some(&combined), None);
        if let Some(c) = ctx.as_deref_mut() {
            c.replace(&e2.0, &tmp_e2.0);
        }
        wire_set_edge_composed(wire, s2, &tmp_e2);
        let first_is_last = first_vertex(&e2)
            .zip(last_vertex(&e2))
            .map(|(f, l)| is_same(&f.0, &l.0))
            .unwrap_or(false);
        if check != ConnectedCheck::Done1 || first_is_last {
            let tmp_e1 = copy_replace_vertices_with(&e1, None, Some(&combined));
            if let Some(c) = ctx.as_deref_mut() {
                c.replace(&e1.0, &tmp_e1.0);
            }
            wire_set_edge_composed(wire, s1, &tmp_e1);
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
pub fn fix_connected_all(wire: &mut Wire, precision: f64, max_tolerance: f64) -> bool {
    fix_connected_all_in_context(wire, None, precision, max_tolerance, true)
}

/// `ShapeFix_Wire::FixConnected(prec)` with `Context()` possibly present
/// (`ShapeFix_Wire.cxx:557-596`).
///
/// With a context the `.cxx` does two things the null-context tool does not:
/// the per-junction refresh of the edge the next iteration reads as `n1`
/// (`cxx:570-588`, `Context()->Apply(aEPrev)` + `WireData()->Set`) and the
/// closing `UpdateWire()` (`cxx:590-593`, `cxx:4293-4310`), which rewrites
/// every WireData entry through the recorded substitutions.
pub fn fix_connected_all_in_context(
    wire: &mut Wire,
    mut ctx: Option<&mut SharedReShape>,
    precision: f64,
    max_tolerance: f64,
    closed: bool,
) -> bool {
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    // `cxx:565`: `aStop = myClosedMode ? 0 : 1`, `for (aI = NbEdges(); aI > aStop; aI--)`.
    let stop = if closed { 0 } else { 1 };
    let mut done = false;
    let mut a_i = nb;
    while a_i > stop {
        done |= fix_connected_one(wire, ctx.as_deref_mut(), precision, max_tolerance, a_i);
        // `cxx:570-588`.
        if ctx.is_some() && a_i - 1 > stop {
            let n1_next = if a_i - 1 > 1 { a_i - 2 } else { nb };
            if n1_next != a_i && n1_next != a_i - 1 {
                let nb_now = wire_edges_nb(wire);
                let edges = stored_manifold_edges(wire);
                if n1_next >= 1 && n1_next <= edges.len() && n1_next <= nb_now {
                    let eprev = edges[n1_next - 1].clone();
                    let c = ctx.as_deref_mut().expect("context checked above");
                    let refresh = c.apply(&eprev.0);
                    if refresh.shape_type() == ShapeType::Edge && !is_same(&refresh, &eprev.0) {
                        // `cxx:587`: `aWireSBWD->Set(TopoDS::Edge(aRefresh), aN1Next)`.
                        let slot = wiredata_storage_index(wire, n1_next, nb_now);
                        wire_set_edge_composed(wire, slot, &Edge(refresh));
                    }
                }
            }
        }
        a_i -= 1;
    }
    if let Some(c) = ctx.as_deref_mut() {
        update_wire_from_context(wire, c);
    }
    done
}

/// `ShapeFix_Wire::UpdateWire()` (`ShapeFix_Wire.cxx:4293-4310`): each
/// `WireData` entry is replaced by `Context()->Apply` of itself when the
/// substitution changed it. The `.cxx` unrolls an entry that `Apply` turned
/// into several edges (`cxx:4303-4306`); `CopyReplaceVertices` never changes
/// the edge count, so here the replacement is one edge.
pub(in crate::shhealing) fn update_wire_from_context(wire: &mut Wire, ctx: &mut SharedReShape) {
    let mut i = 1;
    loop {
        let nb = wire_edges_nb(wire);
        if i > nb {
            break;
        }
        let edges = stored_manifold_edges(wire);
        if i > edges.len() {
            break;
        }
        let e = edges[i - 1].clone();
        // `cxx:4299-4302`: `if (S == E) continue;`.
        let applied = ctx.apply(&e.0);
        if !is_same(&applied, &e.0) && applied.shape_type() == ShapeType::Edge {
            let slot = wiredata_storage_index(wire, i, nb);
            wire_set_edge_composed(wire, slot, &Edge(applied));
        }
        i += 1;
    }
}
