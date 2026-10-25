//! `isTangentToArc` — whether the vector of a common point is tangent to the
//! arc it sits on.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:192-209`
//! (OCCT 8.0.0).
//!
//! `StartSol` asks this before it accepts the support face it walked to as the
//! second face of the band (`_2.cxx:1419-1425`): a common point whose vector is
//! parallel to the arc's tangent leaves through the face it came from.
//!
//! Divergences that are not heuristics:
//! - `BRep_Tool::Curve(Arc, f, l)` is [`BRepTool::edge_curve`]; the two
//!   parameters OCCT discards are not read here either,
//! - `Geom_Curve::D1` is the [`Curve`] trait's `d1`,
//! - `HasVector()` / `Vector()` / `Arc()` / `ParameterOnArc()` raise in OCCT
//!   when the common point does not carry them; the port's `Option` accessors
//!   answer `false` there, which is the "no vector, no answer" of `cxx:194-197`.

use crate::brep_tool::BRepTool;
use crate::chfi_ds::CommonPoint;
use crate::shape::Edge;

/// `isTangentToArc(aCommonPoint, anAngularTolerance)` (`cxx:192-209`).
pub fn is_tangent_to_arc(common_point: &CommonPoint, angular_tolerance: f64) -> bool {
    // `if (!aCommonPoint.HasVector()) return false;` (`cxx:194-197`).
    let Some(vector) = common_point.vector() else {
        return false;
    };
    let Some(arc) = common_point.arc().cloned() else {
        return false;
    };
    let Some(arc) = Edge::wrap(arc) else {
        return false;
    };
    // `const handle(Geom_Curve) anArcCurve = BRep_Tool::Curve(aCommonPoint.Arc(),
    // aDummyParam1, aDummyParam2);` (`cxx:200-201`).
    let Some(curve) = BRepTool::edge_curve(&arc) else {
        return false;
    };
    let Some(param) = common_point.parameter_on_arc() else {
        return false;
    };
    // `anArcCurve->D1(aCommonPoint.ParameterOnArc(), aDummyPoint,
    // anArcTangent)` (`cxx:204-206`).
    let (_, tangent) = curve.d1(param);
    // `aCommonPointVector.IsParallel(anArcTangent, anAngularTolerance)`
    // (`cxx:208`).
    vector.is_parallel_ang(&tangent, angular_tolerance)
}
