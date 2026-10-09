use super::*;

/// `ShapeFix_Root::MinTolerance()` at read time: `FromSTEP.FixShape.MinTolerance3d`
/// is `"1.e-7"` (`STEPControl_Controller.cxx:207`), `ShapeProcess_OperLibrary.cxx:807`
/// feeds it to `ShapeFix_Shape::SetMinTolerance`, which propagates through
/// `ShapeFix_Solid`/`_Shell`/`_Face` down to the wire tool
/// (`ShapeFix_Shape.cxx:339-342`, `ShapeFix_Solid.cxx:737-740`,
/// `ShapeFix_Shell.cxx:1701-1704`, `ShapeFix_Face.cxx:165-168`).
pub(in crate::shhealing) const SHAPE_FIX_MIN_TOLERANCE: f64 = CONFUSION;

/// Status of `ShapeAnalysis_Wire::CheckSmall(num, precsmall)`
/// (`ShapeAnalysis_Wire.cxx:765-835`): `ok` is that function's own return value,
/// `fail` and `done2` are the `ShapeExtend` bits it leaves in `myStatus`
/// (`DONE1` is `ok && !done2`, `DONE2` is `ok && done2`).
pub(in crate::shhealing) struct SmallCheck {
    pub(in crate::shhealing) ok: bool,
    /// `ShapeExtend_FAIL1` / `FAIL2`. OCCT only ORs these into
    /// `ShapeFix_Wire::myLastFixStatus` (`ShapeFix_Wire.cxx:1420-1423`,
    /// `1434-1443`), a status no caller of this port reads.
    #[allow(dead_code)]
    fail: bool,
    done2: bool,
}

/// `ShapeAnalysis_Wire::CheckSmall(num, precsmall)`
/// (`ShapeAnalysis_Wire.cxx:765-835`). `num == 0` means the last edge (`cxx:774`).
pub(in crate::shhealing) fn check_small(wire: &Wire, face: &Face, precsmall: f64, num: usize) -> SmallCheck {
    let nb = wire_edges_nb(wire);
    // `cxx:768-771`: `!IsLoaded() || NbEdges() <= 1`.
    if nb <= 1 {
        return SmallCheck { ok: false, fail: false, done2: false };
    }
    let edges = edges_of_wire(wire);
    let n = if num > 0 { num } else { nb };
    let Some(e) = edges.get(n - 1) else {
        return SmallCheck { ok: false, fail: false, done2: false };
    };
    // `cxx:777-786`: a degenerated edge that still carries a pcurve on the face
    // is left alone; one without a pcurve is a `FAIL1` and falls through.
    let mut fail = false;
    if BRepTool::is_degenerated(e) {
        if curve_on_surface_oriented(e, face, false).is_some() {
            return SmallCheck { ok: false, fail: false, done2: false };
        }
        fail = true;
    }
    // `cxx:787-794`: missing vertices are a `FAIL2`.
    let (Some(v1), Some(v2)) = (first_vertex(e), last_vertex(e)) else {
        return SmallCheck { ok: false, fail: true, done2: false };
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    // `cxx:796`: `Min(myPrecision, precsmall)` is commented out in OCCT, so the
    // caller's `precsmall` is used as-is.
    let prec = precsmall;
    // `cxx:797-800`: `dist > prec` is not small enough (any `FAIL1` survives).
    if p1.distance(&p2) > prec {
        return SmallCheck { ok: false, fail, done2: false };
    }
    // `cxx:806-826`: the midpoint of the 3d curve, else of the pcurve lifted by
    // the surface; with neither, `FAIL1` and the first vertex is used.
    let mid = match BRepTool::edge_curve(e) {
        Some(c3d) => {
            let (cf, cl) = BRepTool::edge_parameters(e);
            c3d.d0(0.5 * (cf + cl))
        }
        None => {
            let lifted = BRepTool::face_surface(face).and_then(|surf| {
                curve_on_surface_oriented(e, face, false).map(|(c2d, cf, cl)| {
                    let p2m = c2d.d0(0.5 * (cf + cl));
                    surf.d0(p2m.x(), p2m.y())
                })
            });
            match lifted {
                Some(m) => m,
                None => {
                    fail = true;
                    p1.clone()
                }
            }
        }
    };
    // `cxx:827-830`.
    if mid.distance(&p1) > prec || mid.distance(&p2) > prec {
        return SmallCheck { ok: false, fail, done2: false };
    }
    // `cxx:832`: `DONE1` when both vertices are the same shape, else `DONE2`.
    SmallCheck { ok: true, fail, done2: !is_same(&v1.0, &v2.0) }
}

/// `ShapeFix_Wire::FixSmall(num, lockvtx, precsmall)`
/// (`ShapeFix_Wire.cxx:1404-1472`).
pub(in crate::shhealing) fn fix_small_one(
    wire: &mut Wire,
    face: &Face,
    my_precision: f64,
    precsmall: f64,
    lockvtx: bool,
    num: usize,
) -> bool {
    let nb = wire_edges_nb(wire);
    // `cxx:1407-1417`: `!IsLoaded() || NbEdges() <= 1`, then the null analyzer
    // guard (this port always has the face, so the analyzer is never null).
    if nb <= 1 {
        return false;
    }
    let n = if num > 0 { num } else { nb };
    let check = check_small(wire, face, precsmall, n);
    // `cxx:1420-1428`: a `FAIL` is recorded but not fatal, and without a `DONE`
    // bit there is nothing to fix.
    if !check.ok {
        return false;
    }
    // `cxx:1433-1444`: a small edge whose vertices are not the same is only
    // removed when `myTopoMode` is on and the vertices are not locked.
    if check.done2 && lockvtx {
        return false;
    }
    // `cxx:1449-1456`: `Context()->Remove(WireData()->Edge(n))` then
    // `WireData()->Remove(n)`, plus the `FixAdvWire.FixSmall.MSG0` warning.
    // UNPORTED: this port has no `ShapeBuild_ReShape` context at read time, so
    // only the wire removal is applied (same convention as `fix_degenerated`);
    // `SendWarning` needs `Message_Msg`, which is not ported either.
    wire_remove_edge(wire, n);
    // `cxx:1459-1468`: `FixConnected(n <= NbEdges() ? n : 1, precsmall)` when the
    // removed edge had two distinct vertices. A `FAIL` there only sets `FAIL3`
    // in `myLastFixStatus`, which the driver ignores, so the result is dropped.
    if check.done2 {
        let nb2 = wire_edges_nb(wire);
        let cnum = if n <= nb2 { n } else { 1 };
        let _ = fix_connected_one(wire, my_precision, precsmall, cnum);
    }
    true
}

/// `ShapeFix_Wire::FixSmall(lockvtx, precsmall)` (`ShapeFix_Wire.cxx:538-553`),
/// the API driver `Perform` calls at `cxx:334` under
/// `NeedFix(myFixSmallMode, myTopoMode)`. `FromSTEP.FixShape.FixSmallMode` is
/// "-1" (`STEPControl_Controller.cxx:233`), so with `myTopoMode` true on the
/// face path (`ShapeFix_Shape.cxx:200`) the pass runs on every read-in wire,
/// with `lockvtx = !myTopoMode || !ReorderOK` and `precsmall = MinTolerance()`.
pub fn fix_small_all(
    wire: &mut Wire,
    face: &Face,
    my_precision: f64,
    precsmall: f64,
    lockvtx: bool,
) -> bool {
    // `cxx:546-550`: `for (int i = NbEdges(); i > 0; i--)`. `i` is initialised
    // only once, so after a removal the following indices address the shrunken
    // list exactly as OCCT's do.
    let mut done = false;
    let mut i = wire_edges_nb(wire);
    while i > 0 {
        done |= fix_small_one(wire, face, my_precision, precsmall, lockvtx, i);
        i -= 1;
    }
    done
}
