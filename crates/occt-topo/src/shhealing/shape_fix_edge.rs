//! `ShapeFix_Edge` subset needed by `ShapeFix_ComposeShell::DispatchWires`:
//! `TempSameRange` (`ShapeFix_Edge.cxx:335-464`) and `FixAddCurve3d`
//! (`ShapeFix_Edge.cxx:618-638`).

use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::brep_lib_same_range::{rep_ranges, set_rep_ranges};
use crate::brep_tool::BRepTool;
use crate::shhealing::{geom_lib_same_range, ShapeBuildEdge};
use crate::shape::Edge;
use crate::tgeometry::GeometryRegistry;

/// `ShapeFix_Edge` (`ShapeFix_Edge.hxx`).
pub struct ShapeFixEdge;

/// File-static helper `TempSameRange(AnEdge, Tolerance)` in `ShapeFix_Edge.cxx:335` (not a `ShapeFix_Edge` member).
///
/// A copy of `BRepLib::SameRange` modified to be able to fix seam edges: the
/// reference range is the 3D curve range when the edge carries a 3D curve,
/// else the first representation's range, and the pcurves that disagree are
/// reparameterised onto it through `GeomLib::SameRange`.
///
/// UNPORTED: the periodic pcurve shift (`cxx:392-400`) and the
/// `Geom2d_BezierCurve::Segment` work-around before the closed-surface remap
/// (`cxx:405-415`); `geom_lib_same_range` performs the reparameterisation for
/// the types it covers.
pub(crate) fn temp_same_range(edge: &Edge) {
    let reg = GeometryRegistry::global();
    let reps = rep_ranges(&edge.0);
    // `cxx:349-353`. The port's `rep_ranges` puts the 3D curve representation
    // first when the edge carries one, matching `BRep_Tool::Curve` at
    // `cxx:349`.
    let Some(head) = reps.first() else {
        return;
    };
    let (current_first, current_last) = (head.first, head.last);
    if !current_first.is_finite() || !current_last.is_finite() {
        return;
    }
    for r in &reps[1..] {
        let Some(key) = r.key else {
            continue;
        };
        // `cxx:382-386`: `abs(first - current_first) > PConfusion || ...`.
        if (r.first - current_first).abs() <= PCONFUSION
            && (r.last - current_last).abs() <= PCONFUSION
        {
            continue;
        }
        // `cxx:387-455`: `GeomLib::SameRange` on the PCurve and, for a
        // closed-surface representation, on the PCurve2.
        let pcs = reg.edge_pcurves(&edge.0, key);
        if pcs.is_empty() {
            continue;
        }
        let new_pcs: Vec<Arc<dyn Curve2d>> = pcs
            .into_iter()
            .map(|pc| geom_lib_same_range(pc, r.first, r.last, current_first, current_last))
            .collect();
        if new_pcs.len() == 1 {
            reg.set_edge_pcurve(&edge.0, key, new_pcs.into_iter().next().unwrap());
        } else {
            reg.set_edge_pcurves(&edge.0, key, new_pcs);
        }
    }
    // `cxx:461-463`: `B.Range(edge, current_first, current_last)` +
    // `B.SameRange(edge, true)`.
    set_rep_ranges(&edge.0, current_first, current_last);
    reg.set_same_range(&edge.0, true);
}

impl ShapeFixEdge {
    /// `FixAddCurve3d(edge)` (`cxx:618-638`).
    ///
    /// Note that the OCCT method takes the edge alone: `ShapeBuild_Edge::
    /// BuildCurve3d(edge)` walks every `BRep_GCurve` representation of the edge
    /// (`BRepLib::BuildCurve3d`, `BRepLib.cxx:301-455`), so no face is needed.
    ///
    /// UNPORTED: the `ShapeExtend` status writes (`cxx:620` OK, `cxx:634` FAIL1,
    /// `cxx:637` DONE1); this port has no `myStatus` accumulator.
    pub fn fix_add_curve3d(&self, edge: &Edge) -> bool {
        let reg = GeometryRegistry::global();
        // `cxx:622`: `BRep_Tool::Degenerated(edge) || EA.HasCurve3d(edge)`.
        if BRepTool::is_degenerated(edge) || reg.edge_curve(&edge.0).is_some() {
            return false;
        }
        // `cxx:626-629`.
        if !reg.same_range(&edge.0) {
            temp_same_range(edge);
        }
        // `cxx:631-635`: `ShapeBuild_Edge().BuildCurve3d(edge)`.
        ShapeBuildEdge.build_curve3d(edge)
    }
}
