//! `ShapeFix_Edge` subset needed by `ShapeFix_ComposeShell::DispatchWires`:
//! `TempSameRange` (`ShapeFix_Edge.cxx:335-464`) and `FixAddCurve3d`
//! (`ShapeFix_Edge.cxx:618-638`).

use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::brep_lib_same_range::{rep_ranges, set_rep_ranges};
use crate::brep_tool::BRepTool;
use crate::shhealing::{geom_lib_same_range, ShapeBuildEdge};
use crate::shape::{Edge, TopoShape};
use crate::shape_extend::{
    decode_status, SHAPEEXTEND_DONE1, SHAPEEXTEND_FAIL1, SHAPEEXTEND_OK,
};
use crate::tgeometry::GeometryRegistry;

use super::{first_vertex, last_vertex};

/// `ShapeFix_Edge` (`ShapeFix_Edge.hxx`). `status` is `myStatus`
/// (`ShapeFix_Edge.hxx:52`), the `ShapeExtend` accumulator queried by callers
/// (e.g. `ShapeFix_Wire::FixAddCurve3d`, `ShapeFix_Wire.cxx:881-918`).
pub struct ShapeFixEdge {
    status: i32,
}

impl ShapeFixEdge {
    /// `ShapeFix_Edge()` (`ShapeFix_Edge.cxx:56-66`): `myStatus` starts unset.
    pub fn new() -> Self {
        Self { status: 0 }
    }

    /// `ShapeFix_Edge::Status(status)` (`ShapeFix_Edge.hxx:41-44`).
    pub fn status(&self, status: i32) -> bool {
        decode_status(self.status, status)
    }
}

impl Default for ShapeFixEdge {
    fn default() -> Self {
        Self::new()
    }
}

/// `ShapeAnalysis_Curve::IsPeriodic(Handle(Geom2d_Curve))`
/// (`ShapeAnalysis_Curve.cxx:1471-1490`): unwrap Offset/Trimmed bases, then
/// `IsPeriodic`.
fn is_periodic_curve2d(c: &dyn Curve2d) -> bool {
    let mut cur: Arc<dyn Curve2d> = Arc::from(c.clone_dyn());
    for _ in 0..32 {
        if let Some(basis) = cur.offset_basis() {
            cur = Arc::from(basis.clone_dyn());
            continue;
        }
        if let Some(basis) = cur.trimmed_basis() {
            cur = Arc::from(basis.clone_dyn());
            continue;
        }
        break;
    }
    cur.is_periodic()
}

/// Rebuild a concrete `Geom2d_BezierCurve` from trait accessors for
/// `IsKind(STANDARD_TYPE(Geom2d_BezierCurve))` + `Segment`.
fn clone_bezier2d(pc: &dyn Curve2d) -> Option<Geom2dBezierCurve> {
    if !pc.is_bezier2d() {
        return None;
    }
    let poles = pc.poles2d()?;
    match pc.bezier_weights2d() {
        Some(w) => Geom2dBezierCurve::rational(poles, w.to_vec()).ok(),
        None => Geom2dBezierCurve::new(poles).ok(),
    }
}

/// File-static helper `TempSameRange(AnEdge, Tolerance)` in `ShapeFix_Edge.cxx:335`
/// (not a `ShapeFix_Edge` member). Callers always pass `Precision::PConfusion()`
/// (`cxx:628`, `:809`, `:833`).
///
/// A copy of `BRepLib::SameRange` modified to be able to fix seam edges: the
/// reference range is the 3D curve range when the edge carries a 3D curve,
/// else the first representation's range, and the pcurves that disagree are
/// reparameterised onto it through `GeomLib::SameRange`.
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
        // `cxx:390-400`: range from the GCurve; periodic shift uses PCurve only
        // (`ShapeAnalysis_Curve::IsPeriodic` on Curve2dPtr), then both Cos share
        // the shifted `[oldFirst, oldLast]`.
        let mut old_first = r.first;
        let mut old_last = r.last;
        if let Some(pc0) = pcs.first() {
            if is_periodic_curve2d(pc0.as_ref()) {
                let tc = Geom2dTrimmedCurve::new(pc0.clone(), old_first, old_last);
                let shift = tc.first_parameter() - old_first;
                old_first += shift;
                old_last += shift;
            }
        }
        let new_pcs: Vec<Arc<dyn Curve2d>> = pcs
            .into_iter()
            .map(|pc| {
                let mut curve = pc;
                let mut first_on = old_first;
                let mut last_on = old_last;
                // `cxx:403-416` / `:432-445`: Bezier Segment onto [oldFirst, oldLast]
                // then SameRange from [0, 1].
                if let Some(mut bez) = clone_bezier2d(curve.as_ref()) {
                    let preci = PCONFUSION;
                    if old_first.abs() > preci || (old_last - 1.0).abs() > preci {
                        bez.segment(old_first, old_last);
                        curve = Arc::new(bez);
                    }
                    first_on = 0.0;
                    last_on = 1.0;
                }
                geom_lib_same_range(
                    PCONFUSION,
                    curve,
                    first_on,
                    last_on,
                    current_first,
                    current_last,
                )
            })
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
    /// Accumulates `myStatus` like OCCT: `OK` at entry (`cxx:620`), `FAIL1` when
    /// `BuildCurve3d` fails (`cxx:634`), else `DONE1` (`cxx:637`). A degenerate
    /// edge or one that already has a 3D curve returns `false` leaving the status
    /// at `OK` (`cxx:622-625`), so `Status(ShapeExtend_FAIL)` is false and
    /// `ShapeFix_Wire::FixAddCurve3d` keeps the edge.
    pub fn fix_add_curve3d(&mut self, edge: &Edge) -> bool {
        // `cxx:620`: `myStatus = ShapeExtend::EncodeStatus(ShapeExtend_OK)`.
        self.status = SHAPEEXTEND_OK;
        let reg = GeometryRegistry::global();
        // `cxx:622`: `BRep_Tool::Degenerated(edge) || EA.HasCurve3d(edge)`.
        if BRepTool::is_degenerated(edge) || reg.edge_curve(&edge.0).is_some() {
            return false;
        }
        // `cxx:626-629`.
        if !reg.same_range(&edge.0) {
            temp_same_range(edge);
        }
        // `cxx:631-637`: `ShapeBuild_Edge().BuildCurve3d(edge)`.
        if !ShapeBuildEdge.build_curve3d(edge) {
            self.status |= SHAPEEXTEND_FAIL1; // cxx:634
            return false;
        }
        self.status |= SHAPEEXTEND_DONE1; // cxx:637
        true
    }

    /// `FixReversed2d(edge, face)` (`ShapeFix_Edge.cxx:735-786`).
    ///
    /// `CheckCurve3dWithPCurve` (`ShapeAnalysis_Edge.cxx:379-425`) returns
    /// without `ShapeExtend_DONE` on a plane, a missing curve, or when the
    /// 3D ends already match the pcurve ends (`CheckPoints`, `cxx:429-445`).
    /// `DONE1` is the crossed pairing, and only that bit reverses the stored
    /// pcurve and rewrites its range.
    pub fn fix_reversed_2d(&self, edge: &Edge, face: &TopoShape) -> bool {
        let reg = GeometryRegistry::global();
        let Some(surf) = reg.face_surface(face) else {
            return false;
        };
        if surf.gp_pln().is_some() {
            return false; // cxx:385-387
        }
        let key = GeometryRegistry::shape_key(face);
        let Some(c2d) = reg.edge_pcurves(&edge.0, key).into_iter().next() else {
            return false; // cxx:392-395 FAIL1
        };
        let Some(c3d) = reg.edge_curve(&edge.0) else {
            return false; // cxx:400-403 FAIL2
        };
        let (Some(v1), Some(v2)) = (first_vertex(edge), last_vertex(edge)) else {
            return false; // cxx:409-411
        };
        let (f3d, l3d) = reg.edge_parameters(&edge.0);
        let (f2d, l2d) = reg.pcurve_range(&edge.0, key).unwrap_or((f3d, l3d));
        let preci1 = BRepTool::vertex_tolerance(&v1);
        let preci2 = BRepTool::vertex_tolerance(&v2);
        let p2d1 = c2d.d0(f2d);
        let p2d2 = c2d.d0(l2d);
        let p1a = c3d.d0(f3d);
        let p1b = c3d.d0(l3d);
        // Face location is identity here: entity 124 is baked into the surface.
        let p2a = surf.d0(p2d1.x(), p2d1.y());
        let p2b = surf.d0(p2d2.x(), p2d2.y());
        if p1a.square_distance(&p2a) <= preci1 * preci1
            && p1b.square_distance(&p2b) <= preci2 * preci2
        {
            return false; // cxx:437-439
        }
        let crossed = p1a.distance(&p2b) + p1b.distance(&p2a);
        let direct = p1a.distance(&p2a) + p1b.distance(&p2b);
        if !(crossed < direct) {
            return false; // status stays OK, so Status(DONE) is false (cxx:760)
        }
        let newf = c2d.reversed_parameter(l2d); // cxx:769
        let newl = c2d.reversed_parameter(f2d);
        let mut reversed = c2d.clone_dyn();
        reversed.reverse(); // cxx:770
        let stored: Arc<dyn Curve2d> = Arc::from(reversed);
        reg.set_edge_pcurve(&edge.0, key, stored);
        reg.set_pcurve_range(&edge.0, key, newf, newl); // cxx:773
        let (first, last) = reg.edge_parameters(&edge.0);
        if first != newf || last != newl {
            reg.set_same_range(&edge.0, false); // cxx:781
            reg.set_same_parameter(&edge.0, false); // cxx:782
        }
        true
    }
}
