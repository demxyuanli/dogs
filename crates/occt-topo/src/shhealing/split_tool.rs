//! `ShapeFix_SplitTool` (`ShapeFix_SplitTool.cxx`): the range-only edge cut.
//!
//! Ported: `CutEdge` (`ShapeFix_SplitTool.cxx:206-303`). It is the "cut" half
//! of `ShapeFix_Wire::FixIntersectingEdges` (`ShapeFix_Wire.cxx:3135`, `:3150`)
//! and is also reached from `ShapeFix_IntersectionTool::CutEdge`
//! (`ShapeFix_IntersectionTool.cxx:194-268`).
//!
//! `SplitEdge` / `SplitEdge1` / `SplitEdge2` (`ShapeFix_SplitTool.cxx:308-731`)
//! are not ported here: they build new topology through a `ShapeBuild_ReShape`
//! context and only `ShapeFix_IntersectionTool` calls them.

use super::prelude::*;

use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::shhealing::shape_analysis_curve::validate_range;
use crate::shhealing::wire_fix::fix_same_parameter;

/// `10. * Precision::PConfusion()`, the guard OCCT applies to both the
/// requested cut length and the edge's own span (`ShapeFix_SplitTool.cxx:216`,
/// `:220`, `:269`, `:273`).
fn ten_pconfusion() -> f64 {
    10.0 * PCONFUSION
}

/// `Geom2d_TrimmedCurve(Geom2d_Line)` test of `ShapeFix_SplitTool.cxx:232-235`.
fn is_trimmed_line(crv: &dyn Curve2d) -> bool {
    crv.trimmed_basis().map(|b| b.is_line()).unwrap_or(false)
}

/// `ShapeFix_SplitTool::CutEdge(edge, pend, cut, face, iscutline)`
/// (`ShapeFix_SplitTool.cxx:206-303`).
///
/// `pend` is the kept end of the current range and `cut` the parameter to cut
/// at; both are *pcurve* parameters of `face` when the pcurve is a trimmed
/// line, `[pend, cut]` is the sub-range otherwise. The edge is modified in
/// place (`BRep_Builder::Range` writes the `BRep_TEdge` of `edge`), and
/// `is_cut_line` reports the trimmed-line branch (`cxx:248`, `:258`).
pub fn cut_edge(edge: &Edge, pend: f64, cut: f64, face: &Face, is_cut_line: &mut bool) -> bool { if std::env::var("PROJDIAG").is_ok() && (pend - 2.078215).abs() < 1e-5 { eprintln!("BTCUT pend={:.6} cut={:.6}\n{}", pend, cut, std::backtrace::Backtrace::force_capture()); }
    // `cxx:214-217`.
    if (cut - pend).abs() < ten_pconfusion() {
        return false;
    }
    let a_range = (cut - pend).abs();
    let reg = GeometryRegistry::global();
    // `cxx:218`: `BRep_Tool::Range(edge, a, b)`, the edge's own 3D range.
    let (a, b) = reg.edge_parameters(&edge.0);
    *is_cut_line = false;
    // `cxx:220-223`.
    if a_range < ten_pconfusion() {
        return false;
    }

    // `cxx:226-265`: with `!SameParameter` the pcurve drives the cut and only
    // a trimmed *line* pcurve is handled; any other pcurve falls through to
    // `return true` at `cxx:265` with the range already written at `cxx:240`.
    if !reg.same_parameter(&edge.0) {
        if let Some((crv, fp, lp)) = curve_on_surface_oriented(edge, face, false) {
            if is_trimmed_line(crv.as_ref()) {
                // `cxx:240`: `B.Range(edge, min(pend,cut), max(pend,cut))` with
                // `Only3d = false`, i.e. *every* representation - the 3D curve
                // and the pcurve - gets the 2D cut interval.
                reg.set_ranges_all(&edge.0, pend.min(cut), pend.max(cut));
                if (pend - lp).abs() < PCONFUSION {
                    // `cxx:241-249`: cut from the beginning.
                    let cut3d = (cut - fp) * (b - a) / (lp - fp);
                    if cut3d <= PCONFUSION {
                        return false;
                    }
                    reg.set_edge_range(&edge.0, a + cut3d, b);
                    *is_cut_line = true;
                } else if (pend - fp).abs() < PCONFUSION {
                    // `cxx:250-259`: cut from the end.
                    let cut3d = (lp - cut) * (b - a) / (lp - fp);
                    if cut3d <= PCONFUSION {
                        return false;
                    }
                    reg.set_edge_range(&edge.0, a, b - cut3d);
                    *is_cut_line = true;
                }
            }
        }
        return true;
    }

    // `cxx:267-277`: the det-study guard and the re-read of the 3D range.
    if ((a - b).abs() - a_range).abs() < PCONFUSION {
        return false;
    }
    if a_range < ten_pconfusion() {
        return false;
    }
    let curve = reg.edge_curve(&edge.0);
    let a = pend.min(cut);
    let b = pend.max(cut);
    let (mut na, mut nb) = (a, b);
    // `cxx:283-298`.
    let fixed = curve.as_ref().is_some_and(|c| {
        !reg.is_degenerated_edge(&edge.0)
            && validate_range(c.as_ref(), &mut na, &mut nb, PCONFUSION)
            && (na != a || nb != b)
    });
    if fixed {
        reg.set_edge_range(&edge.0, na, nb);
        let face_key = GeometryRegistry::shape_key(&face.0);
        // `cxx:288-291`: `ShapeAnalysis_Edge::HasPCurve(edge, face)`.
        if !reg.edge_pcurves(&edge.0, face_key).is_empty() {
            reg.set_same_range(&edge.0, false);
        }
        // `cxx:294`: `ShapeFix_Edge::FixSameParameter(edge)`. OCCT walks every
        // face of `edge` (`ShapeFix_Edge.cxx:324-336`); the port holds only the
        // calling face here, so the fix runs on that one.
        fix_same_parameter(edge, face);
    } else {
        // `cxx:298`: `B.Range(edge, a, b, false)`, all representations.
        reg.set_ranges_all(&edge.0, a, b);
    }
    true
}
