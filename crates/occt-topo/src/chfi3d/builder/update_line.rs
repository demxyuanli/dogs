//! `UpdateLine` — refreshes the extremities of a filling line after a partial
//! invalidation.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:344-387`
//! (OCCT 8.0.0), the static helper of `ChFi3d_Builder` used by
//! `PerformMoreSurfdata` / `PerformIntersectionAtEnd` before `CompleteData`.
//!
//! The `BRepBlend_Extremity::SetValue(P, U, V, Tolerance)` call of
//! `cxx:356` / `:363` / `:375` / `:382` resolves to the four-argument overload
//! `SetValue(const gp_Pnt&, const double W, const double Param, const double
//! Tol)` (`BRepBlend_Extremity.hxx:34-43`): the *surface* V coordinate of the
//! blend point is stored as the parameter on guide. The port keeps that, so
//! [`Extremity::on_curve`] receives `(U, V)` in that order.

use occt_core::precision::CONFUSION;

use crate::brep_blend::extremity::Extremity;
use crate::brep_blend::line::Line;

/// `static void UpdateLine(BRepBlend_Line& Line, const bool isfirst)`
/// (`ChFi3d_Builder_6.cxx:344-387`).
///
/// `isfirst` selects the `Point(1)` / `StartPoint*` half; the other one uses
/// `Point(NbPoints())` / `EndPoint*` and inverts the comparisons. The
/// `Standard_OutOfRange` of `Point` and the `Standard_DomainError` of
/// `ParametersOnS1` / `ParametersOnS2` are answered as `Err`.
pub fn update_line(line: &mut Line, isfirst: bool) -> Result<(), String> {
    let out_of_range = || "BRepBlend_Line::Point: Standard_OutOfRange".to_string();
    if isfirst {
        // `const Blend_Point& BP = Line.Point(1)` (`cxx:350-351`); copied, so the
        // `SetStartPoints` below can borrow the line mutably.
        let bp = *line.point(1).ok_or_else(out_of_range)?;
        let tguide = bp.parameter();
        if line.start_point_on_first().parameter_on_guide() < tguide {
            let (u, v) = bp.parameters_on_s1()?;
            let be = Extremity::on_curve(*bp.point_on_s1(), u, v, CONFUSION);
            let second = line.start_point_on_second().clone();
            line.set_start_points(be, second);
        }
        if line.start_point_on_second().parameter_on_guide() < tguide {
            let (u, v) = bp.parameters_on_s2()?;
            let be = Extremity::on_curve(*bp.point_on_s2(), u, v, CONFUSION);
            let first = line.start_point_on_first().clone();
            line.set_start_points(first, be);
        }
    } else {
        // `const Blend_Point& BP = Line.Point(Line.NbPoints())` (`cxx:370-371`).
        let bp = *line.point(line.nb_points()).ok_or_else(out_of_range)?;
        let tguide = bp.parameter();
        if line.end_point_on_first().parameter_on_guide() > tguide {
            let (u, v) = bp.parameters_on_s1()?;
            let be = Extremity::on_curve(*bp.point_on_s1(), u, v, CONFUSION);
            let second = line.end_point_on_second().clone();
            line.set_end_points(be, second);
        }
        if line.end_point_on_second().parameter_on_guide() > tguide {
            let (u, v) = bp.parameters_on_s2()?;
            let be = Extremity::on_curve(*bp.point_on_s2(), u, v, CONFUSION);
            let first = line.end_point_on_first().clone();
            line.set_end_points(first, be);
        }
    }
    Ok(())
}
