//! `EdgeToBSpline` — the file static of `BRepFill_NSections.cxx:66-129`.
//!
//! Turns one section edge into a `Geom_BSplineCurve` parameterised on
//! `[0, 1]`, running against the edge occurrence (i.e. reversed when the edge
//! is stored `REVERSED`).
//!
//! Translation notes:
//! - `BRep_Tool::Curve(theEdge, aLoc, aFirst, aLast)` plus the later
//!   `aBSCurve->Transform(aLoc.Transformation())` (`cxx:90-116`) is
//!   [`BRepTool::edge_curve_world`]: the entry curve is already in world space,
//!   so the transform step has nothing left to do and no branch is skipped.
//! - `Geom_TrimmedCurve` of this port reparameterises the window onto `[0, 1]`
//!   and reports the basis-space range through `untrimmed_basis()`, which is
//!   what `GeomConvert::CurveToBSplineCurve` cuts against (`GeomConvert.cxx:171-189`).
//! - `Geom_BSplineCurve::Reverse()` (`cxx:124-128`) mutates in place; the port
//!   owns the curve, so [`occt_geom::curve::Curve::reverse`] is the same call.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::convert_approx_curve::GeomConvertApproxCurve;
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::curve::Curve;
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::edge_vertices_cum_ori;
use crate::shape::Edge;

/// `EdgeToBSpline(theEdge)` (`cxx:66-129`).
pub fn edge_to_bspline(edge: &Edge) -> Result<GeomBSplineCurve, String> {
    // `cxx:69-88`: a degenerated edge answers the two-pole point curve
    // `Geom_BSplineCurve(aPoles, {0, 1}, {2, 2}, 1)`, whose poles are the
    // `TopExp::Vertices(theEdge, vl, vf)` pair in reverse: pole 1 is `V2`
    // (`vf`), pole 2 is `V1` (`vl`).
    if BRepTool::is_degenerated(edge) {
        let (v1, v2) = edge_vertices_cum_ori(edge);
        let (Some(v1), Some(v2)) = (v1, v2) else {
            return Err(
                "Standard_NullObject: EdgeToBSpline on a degenerated edge without two vertices"
                    .to_string(),
            );
        };
        let poles = vec![
            BRepTool::vertex_point_world(&v2),
            BRepTool::vertex_point_world(&v1),
        ];
        return GeomBSplineCurve::from_poles_knots_mults(poles, vec![0.0, 1.0], vec![2, 2], 1)
            .map_err(|e| format!("Geom_BSplineCurve: {e}"));
    }

    // `cxx:90-96`: the used part of the edge curve as a trimmed curve.
    let curve = BRepTool::edge_curve_world(edge)
        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
    let (first, last) = BRepTool::edge_parameters(edge);
    let trim: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(Arc::clone(&curve), first, last));

    // `cxx:100-107`: `GeomConvert_ApproxCurve(aCurveTemp, Precision::Confusion(),
    // GeomAbs_C1, 16, 14)`, falling back to
    // `GeomConvert::CurveToBSplineCurve(aTrimCurve)` (default `Convert_QuasiAngular`).
    let approx = GeomConvertApproxCurve::new(&*trim, CONFUSION, Shape::C1, 16, 14);
    let mut bspline = match approx.curve() {
        Some(curve) if approx.has_result() => curve.clone(),
        _ => curve_to_bspline_curve(&*trim, ParameterisationType::QuasiAngular)
            .map_err(|e| format!("GeomConvert::CurveToBSplineCurve: {e:?}"))?,
    };

    // `cxx:118-121`: `BSplCLib::Reparametrize(0., 1., aKnots)` + `SetKnots`.
    let (mut knots, _) = bspline.distinct_knots_and_mults();
    occt_core::bspl::knots::reparametrize(0.0, 1.0, &mut knots);
    bspline
        .set_knots(&knots)
        .map_err(|e| format!("Geom_BSplineCurve::SetKnots: {e}"))?;

    // `cxx:124-128`: a `REVERSED` edge carries the reversed curve.
    if edge.0.orientation() == Orientation::Reversed {
        bspline.reverse();
    }

    Ok(bspline)
}

/// The pole pair of the point curve built for a punctual section edge
/// (`cxx:156-160` and `:236-240`): `Extremities(1) = Pnt(vf)`,
/// `Extremities(2) = Pnt(vl)`, i.e. `TopExp::Vertices(edge, vl, vf)` reversed.
pub fn point_curve_of_edge(edge: &Edge) -> Result<GeomBSplineCurve, String> {
    let (v1, v2) = edge_vertices_cum_ori(edge);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return Err(
            "Standard_NullObject: BRepFill_NSections punctual section without two vertices"
                .to_string(),
        );
    };
    let poles = vec![
        BRepTool::vertex_point_world(&v2),
        BRepTool::vertex_point_world(&v1),
    ];
    GeomBSplineCurve::from_poles_knots_mults(poles, vec![0.0, 1.0], vec![2, 2], 1)
        .map_err(|e| format!("Geom_BSplineCurve: {e}"))
}
