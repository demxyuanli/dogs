//! `BSplSLib_Cache`. Source: `BSplSLib_Cache.cxx` / `BSplSLib_Cache.hxx`.

use crate::bspl::build_cache::build_cache_array;
use crate::bspl::cache_params::CacheParams;
use crate::bspl::coefs::no_derivative_eval_polynomial;
use crate::gp::GpPnt;

/// Cached span evaluator for Bezier / B-spline surfaces.
/// Source: `BSplSLib_Cache`.
#[derive(Debug, Clone)]
pub struct SurfaceCache {
    is_rational: bool,
    params_u: CacheParams,
    params_v: CacheParams,
    /// Flat poles/weights: `rows = max(deg)+1`, each row
    /// `dim * (min(deg)+1)` floats (`hxx:154-159`).
    poles_weights: Vec<f64>,
    /// Row length in floats (`myPolesWeights->RowLength()`).
    cache_cols: usize,
    /// Number of rows (`max(deg)+1`).
    cache_rows: usize,
}

impl SurfaceCache {
    /// `BSplSLib_Cache::BSplSLib_Cache` (`cxx:232-248`).
    pub fn new(
        degree_u: i32,
        periodic_u: bool,
        flat_knots_u: &[f64],
        degree_v: i32,
        periodic_v: bool,
        flat_knots_v: &[f64],
        weights: Option<&[Vec<f64>]>,
    ) -> Self {
        let is_rational = weights.is_some();
        let min_degree = degree_u.min(degree_v);
        let max_degree = degree_u.max(degree_v);
        let pw_col = if is_rational { 4 } else { 3 };
        let cache_cols = (pw_col * (min_degree + 1)) as usize;
        let cache_rows = (max_degree + 1) as usize;
        Self {
            is_rational,
            params_u: CacheParams::new(degree_u, periodic_u, flat_knots_u),
            params_v: CacheParams::new(degree_v, periodic_v, flat_knots_v),
            poles_weights: vec![0.0; cache_rows * cache_cols],
            cache_cols,
            cache_rows,
        }
    }

    /// `IsCacheValid` (`cxx:252-255`).
    pub fn is_cache_valid(&self, u: f64, v: f64) -> bool {
        self.params_u.is_cache_valid(u) && self.params_v.is_cache_valid(v)
    }

    /// `BuildCache` (`cxx:259-297`).
    pub fn build_cache(
        &mut self,
        u: f64,
        v: f64,
        flat_knots_u: &[f64],
        flat_knots_v: &[f64],
        poles: &[Vec<GpPnt>],
        weights: Option<&[Vec<f64>]>,
    ) {
        // `cxx:266-271`: PeriodicNormalization before LocateParameter.
        let mut new_u = self.params_u.periodic_normalization(u);
        let mut new_v = self.params_v.periodic_normalization(v);
        new_u = self.params_u.locate_parameter(new_u, flat_knots_u);
        new_v = self.params_v.locate_parameter(new_v, flat_knots_v);

        // `cxx:276-279`: BSplSLib span mid / half-length convention.
        let span_length_u = 0.5 * self.params_u.span_length;
        let span_start_u = self.params_u.span_start + span_length_u;
        let span_length_v = 0.5 * self.params_v.span_length;
        let span_start_v = self.params_v.span_start + span_length_v;

        let flat = build_cache_array(
            span_start_u,
            span_start_v,
            span_length_u,
            span_length_v,
            self.params_u.is_periodic,
            self.params_v.is_periodic,
            self.params_u.degree,
            self.params_v.degree,
            self.params_u.span_index,
            self.params_v.span_index,
            flat_knots_u,
            flat_knots_v,
            poles,
            weights,
            self.is_rational,
            self.cache_rows,
            self.cache_cols,
        );
        self.poles_weights = flat;
        let _ = new_u;
        let _ = new_v;
    }

    /// `D0` (`cxx:301-305`).
    pub fn d0(&self, u: f64, v: f64) -> GpPnt {
        let (local_u, local_v) = self.to_local_params_d0(u, v);
        self.d0_local(local_u, local_v)
    }

    /// `toLocalParamsD0` (`cxx:90-102`).
    fn to_local_params_d0(&self, u: f64, v: f64) -> (f64, f64) {
        let new_u = self.params_u.periodic_normalization(u);
        let new_v = self.params_v.periodic_normalization(v);
        let span_length_u = 0.5 * self.params_u.span_length;
        let span_start_u = self.params_u.span_start + span_length_u;
        let span_length_v = 0.5 * self.params_v.span_length;
        let span_start_v = self.params_v.span_start + span_length_v;
        (
            (new_u - span_start_u) / span_length_u,
            (new_v - span_start_v) / span_length_v,
        )
    }

    /// `D0Local` (`cxx:309-345`).
    pub fn d0_local(&self, local_u: f64, local_v: f64) -> GpPnt {
        let dimension = if self.is_rational { 4 } else { 3 };
        let is_max_u = self.params_u.degree > self.params_v.degree;
        let (min_degree, max_degree) = if self.params_u.degree <= self.params_v.degree {
            (self.params_u.degree, self.params_v.degree)
        } else {
            (self.params_v.degree, self.params_u.degree)
        };
        let (min_param, max_param) = if is_max_u {
            (local_v, local_u)
        } else {
            (local_u, local_v)
        };

        let transient = no_derivative_eval_polynomial(
            max_param,
            max_degree,
            self.cache_cols as i32,
            &self.poles_weights,
        );
        let point = no_derivative_eval_polynomial(
            min_param,
            min_degree,
            dimension,
            &transient,
        );
        let mut p = GpPnt::new(
            point.first().copied().unwrap_or(0.0),
            point.get(1).copied().unwrap_or(0.0),
            point.get(2).copied().unwrap_or(0.0),
        );
        if self.is_rational {
            let w = point.get(3).copied().unwrap_or(1.0);
            if w.abs() > 0.0 {
                p = GpPnt::new(p.x() / w, p.y() / w, p.z() / w);
            }
        }
        let _ = self.cache_rows;
        p
    }
}
