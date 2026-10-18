//! `FixAddCurve3d` arm of `ShapeFix_Wire::FixEdgeCurves`
//! (`ShapeFix_Wire.cxx:881-918`).
//!
//! `NeedFix(myFixAddCurve3dMode)` with the tool default `-1`
//! (`ShapeFix_Wire.cxx:181`) is true, so the arm runs on every wire of the
//! IGES and STEP reads. Each edge goes through
//! `ShapeFix_Edge::FixAddCurve3d` (`ShapeFix_Edge.cxx:618-638`), which leaves a
//! degenerate edge or one that already carries a 3D curve alone
//! (`ShapeAnalysis_Edge::HasCurve3d` is `BRep_Tool::Curve`'s null-ness) and
//! otherwise builds the 3D curve from the pcurve through
//! `ShapeBuild_Edge::BuildCurve3d` (`ShapeBuild_Edge.cxx:714-775` ->
//! `BRepLib::BuildCurve3d`, `BRepLib.cxx:301-455`).
//!
//! The arm also owns the failed-build clean-up (`cxx:888-917`): when the build
//! fails and the edge holds no usable pcurve, the edge is removed from the wire
//! and `FixConnected` runs at the removal index.

use crate::shhealing::shape_fix_edge::ShapeFixEdge;
use crate::shape::{Face, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::edges_of_wire;

use super::degenerated::fix_degenerated;
use super::lacking::fix_lacking_one;
use super::{fix_connected_at, remove_wire_edge};

/// `FixAddCurve3d` (`ShapeFix_Wire.cxx:881-918`). Reuses one `ShapeFix_Edge`
/// (OCCT's `myFixEdge`) and queries its `Status` after each call, mirroring
/// `ShapeFix_Shape::Perform`'s wire arm (`ShapeFix_Shape.cxx:228-232`), where a
/// `true` `Perform()` gates `Context()->Replace(S, sfw->Wire())`.
///
/// `Status(ShapeExtend_DONE)` (a 3D curve was built) sets `ShapeExtend_DONE5`;
/// `Status(ShapeExtend_FAIL)` (the build failed) removes the edge when it has no
/// usable pcurve, else keeps it, and sets `FAIL5` (`cxx:888-917`). An edge that
/// is degenerate or already carries a 3D curve leaves the status at `OK`
/// (`ShapeFix_Edge.cxx:622-625`), so it stays. The return value is the DONE5
/// accumulation (3D curve added or edge dropped).
///
/// `face` is the wire's own face: OCCT reaches `FixClosed` (`cxx:909`) whose
/// face-dependent sub-calls read `myAnalyzer->Surface()`. A free wire passes the
/// bare (unloaded) face, so they return false, as OCCT's null analyzer does.
pub fn fix_add_curve3d_or_drop(wire: &mut Wire, face: &Face, precision: f64) -> bool {
    use crate::shape_extend::{SHAPEEXTEND_DONE, SHAPEEXTEND_FAIL};
    let reg = GeometryRegistry::global();
    let mut sfe = ShapeFixEdge::new();
    let mut changed = false;
    let mut index = 1usize;
    let mut nb = edges_of_wire(wire).len();
    while index <= nb {
        let edges = edges_of_wire(wire);
        let edge = &edges[index - 1];
        // `cxx:883`: `myFixEdge->FixAddCurve3d(sbwd->Edge(i))`.
        sfe.fix_add_curve3d(edge);
        // `cxx:884-887`: `Status(ShapeExtend_DONE)` -> `DONE5`.
        if sfe.status(SHAPEEXTEND_DONE) {
            changed = true;
        }
        // `cxx:888-917`: `Status(ShapeExtend_FAIL)` -> maybe remove the edge.
        if sfe.status(SHAPEEXTEND_FAIL) {
            let reps = reg.edge_pcurve_reps(&edge.0);
            let missing = reps.is_empty()
                || reps
                    .iter()
                    .all(|(_, _, a, b)| (b - a).abs() < occt_core::precision::PCONFUSION);
            if missing {
                let removed_at = index;
                remove_wire_edge(wire, removed_at);
                changed = true; // `cxx:906` DONE5 on removal
                nb = edges_of_wire(wire).len();
                // `cxx:905`: `sbwd->Remove(i--)`, so the loop's own `i++` puts
                // `i` back on the same slot, which now holds the next edge. The
                // two arms below test that decremented index.
                let new_i = removed_at - 1;
                if new_i == nb {
                    // `cxx:909`: `FixClosed(Precision())`.
                    let _ = fix_closed(wire, face, occt_core::precision::CONFUSION);
                } else {
                    // `cxx:913`: `FixConnected(i + 1, Precision())`.
                    let _ = fix_connected_at(wire, precision, precision, new_i + 1);
                }
                continue; // same `i` as OCCT: `i--` then `i++`
            }
        }
        index += 1;
    }
    changed
}

/// `ShapeFix_Wire::FixClosed(prec)` (`ShapeFix_Wire.cxx:1304-1343`).
///
/// `FixConnected(1, prec)` (`cxx:1311`), `FixDegenerated(1)` (`cxx:1321`) and
/// `FixLacking(1)` (`cxx:1331`), in that order, each folded into
/// `myStatusClosed`. The caller at `cxx:909` discards the return value, so the
/// accumulated `ShapeExtend_DONE` is all that is reported.
fn fix_closed(wire: &mut Wire, face: &Face, prec: f64) -> bool {
    let mut done = false;
    // `cxx:1311-1320`.
    done |= fix_connected_at(wire, prec, prec, 1);
    // `cxx:1321-1330`: `FixDegenerated(1)`, precision `Precision()` inside the
    // method (`ShapeFix_Wire.cxx:2130-2205`).
    done |= fix_degenerated(wire, face, prec, 1).0;
    // `cxx:1331-1340`: `FixLacking(1)`, `force` false
    // (`ShapeFix_Wire.cxx:1275-1299`).
    done |= fix_lacking_one(wire, face, 1, false, prec, super::SHAPE_FIX_MAX_TOLERANCE);
    done
}
