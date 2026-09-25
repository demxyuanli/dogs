//! `XSAlgo_ShapeProcessor::CheckPCurve` — the "advanced check" the STEP reader
//! runs on every edge of a translated `EDGE_LOOP`.
//!
//! Call chain in OCCT:
//! `StepToTopoDS_TranslateEdgeLoop.cxx:875` → `CheckPCurves` (`:105-177`) →
//! `XSAlgo_ShapeProcessor::CheckPCurve` (`XSAlgo_ShapeProcessor.cxx:344`).
//!
//! Only the two *early* rejection tests are ported here (they are the ones that
//! drop a pcurve the file got wrong, after which `ShapeFix_Edge::FixAddPCurve`
//! re-projects it during the `ShapeFix` pass):
//!   * `:360-376` — a pcurve whose U or V span exceeds 6/8 of the surface span
//!     wraps around the surface (e.g. UV in degrees on a surface in radians), so
//!     it is discarded;
//!   * `:378-401` — the pcurve's end points must agree with the 3D curve's end
//!     points within the read precision.
//!
//! OCCT tests the single pcurve `ShapeAnalysis_Edge::PCurve(edge, face, ...,
//! false)` yields; a face can carry two (seam / two-sided edge), so every stored
//! pcurve is tested here and only the rejected ones are dropped — dropping the
//! good one as well would force a needless re-projection.
//!
//! UNPORTED: `:403-464` (deviation between pcurve and 3D curve over the whole
//! edge, then "best of the two" re-projection via `FixSameParameter` /
//! `FixAddPCurve`, `:466-490` write-back) and the `theIsSeam` two-pcurve write
//! (`:414-434`).

use super::prelude::*;
use crate::brep_tool::BRepTool;
use occt_geom2d::Curve2d;

/// Returns `true` when the edge keeps at least one pcurve on `face`.
///
/// Mirrors `XSAlgo_ShapeProcessor::CheckPCurve` (`XSAlgo_ShapeProcessor.cxx:344`)
/// with `thePrecision = preci`; rejected pcurves are removed
/// (`ShapeBuild_Edge().RemovePCurve`, `cxx:374` / `cxx:399`).
pub fn xsalgo_check_pcurve(edge: &Edge, face: &Face, preci: f64) -> bool {
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    let stored = reg.edge_pcurves(&edge.0, face_key);
    if stored.is_empty() {
        return false;
    }
    let Some(surface) = BRepTool::face_surface(face) else {
        return true;
    };
    let (u1, u2) = surface.u_range();
    let (v1, v2) = surface.v_range();
    let (a, b) = BRepTool::edge_parameters(edge);
    let curve3d = BRepTool::edge_curve_world(edge);
    let (vtx1, vtx2) = edge_vertices(edge);
    let pv1 = match curve3d.as_ref() {
        Some(c) if a.is_finite() => c.d0(a),
        _ => vtx1
            .as_ref()
            .map(vertex_position)
            .unwrap_or_else(GpPnt::zero),
    };
    let pv2 = match curve3d.as_ref() {
        Some(c) if b.is_finite() => c.d0(b),
        _ => vtx2
            .as_ref()
            .map(vertex_position)
            .unwrap_or_else(GpPnt::zero),
    };

    let params = |c2d: &Arc<dyn Curve2d>| -> Option<(f64, f64)> {
        let (mut t1, mut t2) = reg
            .pcurve_range(&edge.0, face_key)
            .unwrap_or_else(|| reg.edge_parameters(&edge.0));
        if !(t1.is_finite() && t2.is_finite()) {
            t1 = c2d.first_parameter();
            t2 = c2d.last_parameter();
        }
        (t1.is_finite() && t2.is_finite()).then_some((t1, t2))
    };

    let total = stored.len();
    let mut kept: Vec<Arc<dyn Curve2d>> = Vec::new();
    for c2d in stored {
        let Some((t1, t2)) = params(&c2d) else {
            // `cxx:355-358`: no parameters means OCCT returns without touching
            // the pcurve.
            kept.push(c2d);
            continue;
        };
        let q1 = c2d.d0(t1);
        let q2 = c2d.d0(t2);
        // `XSAlgo_ShapeProcessor.cxx:367-376`: `anEdgeSpanX / 8. > (U2 / 6. - U1 / 6.)`.
        let span_x = (q1.x() - q2.x()).abs();
        let span_y = (q1.y() - q2.y()).abs();
        if span_x / 8.0 > (u2 / 6.0 - u1 / 6.0) || span_y / 8.0 > (v2 / 6.0 - v1 / 6.0) {
            continue;
        }
        // `XSAlgo_ShapeProcessor.cxx:378-401`: `aCurve3DPoint = Surface->Value(uv)`
        // against `aCurve3D->Value(param)` (or the vertex point when the edge has
        // no 3D curve); a distance above `thePrecision` drops the pcurve.
        let s1 = surface.d0(q1.x(), q1.y());
        let s2 = surface.d0(q2.x(), q2.y());
        if pv1.distance(&s1) > preci || pv2.distance(&s2) > preci {
            continue;
        }
        kept.push(c2d);
    }

    if kept.len() < total {
        if kept.is_empty() {
            reg.remove_pcurves_on_surface(&edge.0, &face.0);
            return false;
        }
        reg.set_edge_pcurves(&edge.0, face_key, kept);
    }
    true
}
