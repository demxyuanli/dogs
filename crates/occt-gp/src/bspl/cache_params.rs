//! `BSplCLib_CacheParams`. Source: `BSplCLib_CacheParams.hxx:25-121`.

use crate::bspl::locate::locate_parameter;
use occt_foundation::precision::epsilon;

/// One-direction cache span state for Bezier / B-spline evaluation.
/// Source: `BSplCLib_CacheParams`.
#[derive(Debug, Clone)]
pub struct CacheParams {
    /// `Degree`.
    pub degree: i32,
    /// `IsPeriodic`.
    pub is_periodic: bool,
    /// `FirstParameter`.
    pub first_parameter: f64,
    /// `LastParameter`.
    pub last_parameter: f64,
    /// `SpanIndexMin`.
    pub span_index_min: i32,
    /// `SpanIndexMax`.
    pub span_index_max: i32,
    /// `SpanStart`.
    pub span_start: f64,
    /// `SpanLength`.
    pub span_length: f64,
    /// `SpanIndex` (1-based flat-knot index).
    pub span_index: i32,
}

impl CacheParams {
    /// `BSplCLib_CacheParams(Degree, Periodic, FlatKnots)` (`hxx:43-56`).
    pub fn new(degree: i32, is_periodic: bool, flat_knots: &[f64]) -> Self {
        let lower = 1i32;
        let upper = flat_knots.len() as i32;
        let first = flat_knots
            .get((lower + degree - 1) as usize)
            .copied()
            .unwrap_or(0.0);
        let last = flat_knots
            .get((upper - degree - 1) as usize)
            .copied()
            .unwrap_or(0.0);
        Self {
            degree,
            is_periodic,
            first_parameter: first,
            last_parameter: last,
            span_index_min: lower + degree,
            span_index_max: upper - degree - 1,
            span_start: 0.0,
            span_length: 0.0,
            span_index: 0,
        }
    }

    /// `PeriodicNormalization` (`hxx:60-77`).
    pub fn periodic_normalization(&self, parameter: f64) -> f64 {
        if !self.is_periodic {
            return parameter;
        }
        let period = self.last_parameter - self.first_parameter;
        if parameter < self.first_parameter {
            let scale = ((self.first_parameter - parameter) / period).trunc();
            return parameter + period * (scale + 1.0);
        }
        if parameter > self.last_parameter {
            let scale = ((parameter - self.last_parameter) / period).trunc();
            return parameter - period * (scale + 1.0);
        }
        parameter
    }

    /// `IsCacheValid` (`hxx:81-99`).
    pub fn is_cache_valid(&self, parameter: f64) -> bool {
        let a_new = self.periodic_normalization(parameter);
        let delta = a_new - self.span_start;
        if (delta < 0.0 && self.span_index != self.span_index_min)
            || (delta >= self.span_length && self.span_index != self.span_index_max)
        {
            return false;
        }
        if self.span_index == self.span_index_max {
            return true;
        }
        let an_eps = epsilon(self.last_parameter.abs().min(a_new.abs()));
        let delta_to_next = (delta - self.span_length).abs();
        delta_to_next > an_eps
    }

    /// `LocateParameter` (`hxx:104-116`). Updates span fields; returns the
    /// periodically-normalized parameter.
    pub fn locate_parameter(&mut self, parameter: f64, flat_knots: &[f64]) -> f64 {
        let mut u = parameter;
        self.span_index = 0;
        let (idx, new_u) = locate_parameter(
            self.degree,
            flat_knots,
            None,
            u,
            self.is_periodic,
            self.span_index,
        );
        self.span_index = idx;
        u = new_u;
        let i = self.span_index as usize;
        self.span_start = flat_knots.get(i.wrapping_sub(1)).copied().unwrap_or(0.0);
        let next = flat_knots.get(i).copied().unwrap_or(self.span_start);
        self.span_length = next - self.span_start;
        u
    }
}
