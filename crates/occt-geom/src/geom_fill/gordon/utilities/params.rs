//! Parameter helpers of `GeomFill_GordonUtilities.pxx` (`pxx:38-118`, `:194-210`,
//! `:1296-1337`).

use crate::bspline_curve::MAX_DEGREE;
use occt_core::precision::{APPROXIMATION, INTERSECTION, PCONFUSION};

/// `interpolationDegree` (`pxx:38-42`).
pub fn interpolation_degree() -> i32 {
    3i32.min(MAX_DEGREE as i32)
}

/// One-arg `Precision::Parametric(P)` (`Precision.hxx:328`) = `P * 0.01`.
#[inline]
fn parametric_1(p: f64) -> f64 {
    p * PCONFUSION / occt_core::precision::CONFUSION
}

/// `Precision::PIntersection(T)` (`pxx:45-48` via `Precision.hxx`).
#[inline]
fn p_intersection(t: f64) -> f64 {
    parametric_1(INTERSECTION) / t
}

/// `parameterResolution` (`pxx:44-48`).
pub fn parameter_resolution(first: f64, last: f64) -> f64 {
    p_intersection((last - first).abs().max(1.0))
}

/// `boundaryParameterResolution` (`pxx:50-57`).
pub fn boundary_parameter_resolution(first: f64, last: f64, tolerance: f64) -> f64 {
    let range = (last - first).abs().max(1.0);
    parametric_1(tolerance.max(APPROXIMATION)) * range
}

/// `appendUnitParameter` (`pxx:59-63`).
pub fn append_unit_parameter(parameters: &mut Vec<f64>, parameter: f64) {
    parameters.push(parameter.clamp(0.0, 1.0));
}

/// `compactParameters` (`pxx:65-89`).
pub fn compact_parameters(parameters: &mut Vec<f64>) {
    if parameters.is_empty() {
        return;
    }
    parameters.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let tol = parameter_resolution(0.0, 1.0);
    let mut write_idx = 0usize;
    for read_idx in 1..parameters.len() {
        if (parameters[read_idx] - parameters[write_idx]).abs() > tol {
            write_idx += 1;
            parameters[write_idx] = parameters[read_idx];
        }
    }
    parameters.truncate(write_idx + 1);
}

/// `addIntervalMidpoints` (`pxx:91-118`).
pub fn add_interval_midpoints(parameters: &mut Vec<f64>, min_size: usize) {
    compact_parameters(parameters);
    while parameters.len() < min_size && parameters.len() > 1 {
        let mut split_idx = 0usize;
        let mut max_gap = -1.0;
        for param_idx in 1..parameters.len() {
            let gap = parameters[param_idx] - parameters[param_idx - 1];
            if gap > max_gap {
                split_idx = param_idx;
                max_gap = gap;
            }
        }
        if max_gap <= parameter_resolution(0.0, 1.0) {
            return;
        }
        append_unit_parameter(
            parameters,
            0.5 * (parameters[split_idx - 1] + parameters[split_idx]),
        );
        compact_parameters(parameters);
    }
}

/// `OrderedCurve` (`pxx:187-192`).
#[derive(Clone, Copy, Debug)]
pub struct OrderedCurve {
    pub index: usize,
    pub rank: f64,
}

/// `columnMeans` (`pxx:194-201`): average over profile rows for each guide column.
pub fn column_means(params: &[Vec<f64>]) -> Vec<f64> {
    if params.is_empty() {
        return Vec::new();
    }
    let nb_cols = params[0].len();
    let nb_rows = params.len() as f64;
    let mut sums = vec![0.0; nb_cols];
    for row in params {
        for (col, &v) in row.iter().enumerate() {
            sums[col] += v;
        }
    }
    for s in &mut sums {
        *s /= nb_rows;
    }
    sums
}

/// `rowMeans` (`pxx:203-210`): average over guide columns for each profile row.
pub fn row_means(params: &[Vec<f64>]) -> Vec<f64> {
    if params.is_empty() {
        return Vec::new();
    }
    let nb_cols = params[0].len() as f64;
    params
        .iter()
        .map(|row| row.iter().sum::<f64>() / nb_cols)
        .collect()
}

/// `pushConstraint` (`pxx:1296-1315`).
pub fn push_constraint(
    source: &mut Vec<f64>,
    target: &mut Vec<f64>,
    source_param: f64,
    target_param: f64,
    source_tol: f64,
    target_tol: f64,
) {
    if let (Some(&last_s), Some(&last_t)) = (source.last(), target.last()) {
        if (source_param - last_s).abs() <= source_tol || (target_param - last_t).abs() <= target_tol
        {
            return;
        }
    }
    source.push(source_param);
    target.push(target_param.clamp(0.0, 1.0));
}

/// `samplingParameters` (`pxx:1317-1337`).
pub fn sampling_parameters(targets: &[f64], nb_source_poles: usize) -> Vec<f64> {
    let mut samples = Vec::new();
    append_unit_parameter(&mut samples, 0.0);
    append_unit_parameter(&mut samples, 1.0);
    for &t in targets {
        append_unit_parameter(&mut samples, t);
    }
    compact_parameters(&mut samples);
    let min_size = (interpolation_degree() as usize + 1)
        .max(nb_source_poles * 2)
        .max(targets.len() + 2);
    add_interval_midpoints(&mut samples, min_size);
    samples
}
