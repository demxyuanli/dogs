//! `BuildTangents` / `BuildPeriodicTangent` / `ScaleTangents` — the free
//! functions of `GeomAPI_Interpolate.cxx:113-278` that fill the end tangents
//! of a degree-3 interpolation from a Lagrange derivative of the point data.
//!
//! `PLib::EvalLagrange` is called with `DerivativeRequest = 1`,
//! `Dimension = 3` and `Degree` 1/2/3 depending on the point count, and the
//! callers read `eval_result[1][*]` — the first derivative — from it.

use super::*;
use occt_core::bspl::plib::eval_lagrange;

/// `(degree + 1) * 3` coordinates of the point window starting at `start`,
/// i.e. the raw `gp_Pnt` storage OCCT hands to `PLib::EvalLagrange` as
/// `(double*)&PointsArray.Value(index)`.
fn point_window(points: &[GpPnt], start: usize, degree: usize) -> Vec<f64> {
    let mut values = Vec::with_capacity((degree + 1) * 3);
    for i in 0..=degree {
        let p = &points[start + i];
        values.push(p.x());
        values.push(p.y());
        values.push(p.z());
    }
    values
}

/// `PLib::EvalLagrange(Parameter, 1, Degree, 3, ...)` over the window, with
/// `eval_result[2][3]` flattened into a 6-element array.
fn lagrange_derivative(
    points: &[GpPnt],
    parameters: &[f64],
    start: usize,
    degree: usize,
    parameter: f64,
) -> Option<[f64; 6]> {
    let values = point_window(points, start, degree);
    let params = &parameters[start..start + degree + 1];
    let mut results = [0.0f64; 6];
    if eval_lagrange(parameter, 1, degree, 3, &values, params, &mut results) != 0 {
        return None;
    }
    Some(results)
}

/// `BuildTangents` (`GeomAPI_Interpolate.cxx:155-213`): the end tangents of a
/// non-periodic degree-3 interpolation, from a degree-3 (or degree-2 for three
/// points) Lagrange derivative.
pub(super) fn build_tangents(
    points: &[GpPnt],
    tangents: &mut [GpVec],
    flags: &mut [bool],
    parameters: &[f64],
) -> Result<(), String> {
    let n = points.len();
    if n < 3 {
        return Err("Standard_ConstructionError: GeomAPI_Interpolate".into());
    }
    let degree = if n == 3 { 2 } else { 3 };

    if !flags[0] {
        // `cxx:177-196`: the first point, `parameter_array = &Parameters(1)`.
        let res = lagrange_derivative(points, parameters, 0, degree, parameters[0])
            .ok_or("Standard_ConstructionError: GeomAPI_Interpolate")?;
        flags[0] = true;
        tangents[0] = GpVec::new(res[3], res[4], res[5]);
    }
    if !flags[n - 1] {
        // `cxx:197-212`: `point_array = &Points(Upper - degree)`.
        let start = n - 1 - degree;
        let res = lagrange_derivative(points, parameters, start, degree, parameters[n - 1])
            .ok_or("Standard_ConstructionError: GeomAPI_Interpolate")?;
        flags[n - 1] = true;
        tangents[n - 1] = GpVec::new(res[3], res[4], res[5]);
    }
    Ok(())
}

/// `BuildPeriodicTangent` (`GeomAPI_Interpolate.cxx:113-153`): only the tangent
/// at the first point is needed, since a periodic curve closes on itself.
pub(super) fn build_periodic_tangent(
    points: &[GpPnt],
    tangents: &mut [GpVec],
    flags: &mut [bool],
    parameters: &[f64],
) -> Result<(), String> {
    let n = points.len();
    if n < 3 {
        return Err("Standard_ConstructionError: GeomAPI_Interpolate".into());
    }
    if !flags[0] {
        let degree = if n == 3 { 2 } else { 3 };
        let res = lagrange_derivative(points, parameters, 0, degree, parameters[0])
            .ok_or("Standard_ConstructionError: GeomAPI_Interpolate")?;
        flags[0] = true;
        tangents[0] = GpVec::new(res[3], res[4], res[5]);
    }
    Ok(())
}

/// `ScaleTangents` (`GeomAPI_Interpolate.cxx:221-278`): rescale every flagged
/// tangent so its coordinate sum matches the Lagrange derivative's.
pub(super) fn scale_tangents(
    points: &[GpPnt],
    tangents: &mut [GpVec],
    flags: &[bool],
    parameters: &[f64],
) {
    let num_points = points.len();
    let degree = if num_points == 2 {
        1
    } else if num_points >= 3 {
        2
    } else {
        0
    };
    let mut index = 0usize;
    let max_index = points.len().saturating_sub(1 + degree);
    for ii in 0..flags.len() {
        if !flags[ii] {
            continue;
        }
        if index > max_index {
            index = max_index;
        }
        let Some(res) = lagrange_derivative(points, parameters, index, degree, parameters[ii])
        else {
            continue;
        };
        let value0 = tangents[ii].x().abs() + tangents[ii].y().abs() + tangents[ii].z().abs();
        let value1 = res[3].abs() + res[4].abs() + res[5].abs();
        let ratio = value1 / value0;
        tangents[ii] = tangents[ii].multiplied_scalar(ratio);
        if ii != 0 {
            index += 1;
        }
        if index > max_index {
            index = max_index;
        }
    }
}

impl GeomAPIInterpolate {
    /// `ScaleTangents(myPoints, myTangents, myTangentFlags, myParameters)`
    /// (`GeomAPI_Interpolate.cxx:377-383`, `:409-415`).
    pub(super) fn scale_tangents(&mut self) {
        scale_tangents(
            &self.points,
            &mut self.tangents,
            &self.tangent_flags,
            &self.parameters,
        );
    }
}
