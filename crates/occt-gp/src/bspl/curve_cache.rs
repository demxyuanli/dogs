//! `BSplCLib_Cache` for 2D curves.
//! Source: `BSplCLib_Cache.cxx` / `BSplCLib_Cache.hxx`,
//! `BSplCLib_BuildCache` Array2 (`BSplCLib_CurveComputation.pxx:1613-1672`).

use crate::bspl::cache_params::CacheParams;
use crate::bspl::coefs::no_derivative_eval_polynomial;
use crate::bspl::curve_eval::prepare_curve_eval_2d;
use crate::bspl::bohm::bohm;
use crate::gp::GpPnt2d;

/// Cached span evaluator for Bezier / B-spline curves (`BSplCLib_Cache`).
#[derive(Debug, Clone)]
pub struct CurveCache {
    is_rational: bool,
    params: CacheParams,
    /// Row length in floats (`myRowLength`): 2 or 3 (with weight).
    row_length: usize,
    /// Flat poles/weights: `(degree+1) * row_length` floats.
    poles_weights: Vec<f64>,
}

impl CurveCache {
    /// `BSplCLib_Cache` 2d constructor (`cxx:28-38`).
    pub fn new_2d(
        degree: i32,
        periodic: bool,
        flat_knots: &[f64],
        weights: Option<&[f64]>,
    ) -> Self {
        let is_rational = weights.is_some();
        let row_length = if is_rational { 3 } else { 2 };
        Self {
            is_rational,
            params: CacheParams::new(degree, periodic, flat_knots),
            row_length,
            poles_weights: vec![0.0; ((degree + 1) as usize) * row_length],
        }
    }

    /// `IsCacheValid` (`cxx:56-59`).
    pub fn is_cache_valid(&self, parameter: f64) -> bool {
        self.params.is_cache_valid(parameter)
    }

    /// `BuildCache` for 2d poles (`cxx:63-89`).
    pub fn build_cache_2d(
        &mut self,
        parameter: f64,
        flat_knots: &[f64],
        poles: &[GpPnt2d],
        weights: Option<&[f64]>,
    ) {
        // `BSplCLib_Cache::BuildCache` (`cxx:63-89`): PeriodicNormalization,
        // LocateParameter, then `BSplCLib::BuildCache(SpanStart, SpanLength, ...)`.
        // Curve cache uses the span start / full length — NOT the BSplSLib mid-span
        // / half-length frame used by `SurfaceCache`.
        let mut new_u = self.params.periodic_normalization(parameter);
        new_u = self.params.locate_parameter(new_u, flat_knots);
        self.poles_weights = build_curve_cache_array_2d(
            self.params.span_start,
            self.params.span_length,
            self.params.is_periodic,
            self.params.degree,
            self.params.span_index,
            flat_knots,
            poles,
            weights,
            self.is_rational,
        );
        let _ = new_u;
    }

    /// `D0` for `gp_Pnt2d` (`cxx:190-195`).
    pub fn d0_2d(&self, parameter: f64) -> GpPnt2d {
        // `BSplCLib_Cache::D0` (`cxx:190-195`):
        // `(param - SpanStart) / SpanLength` — full span, not mid-span.
        let local = self.params.periodic_normalization(parameter);
        let local = (local - self.params.span_start) / self.params.span_length;
        self.d0_local_2d(local)
    }

    /// `D0Local` for `gp_Pnt2d` (`cxx:334-350`).
    fn d0_local_2d(&self, local_param: f64) -> GpPnt2d {
        let point = no_derivative_eval_polynomial(
            local_param,
            self.params.degree,
            self.row_length as i32,
            &self.poles_weights,
        );
        let mut p = GpPnt2d::new(
            point.first().copied().unwrap_or(0.0),
            point.get(1).copied().unwrap_or(0.0),
        );
        if self.is_rational {
            let w = point.get(2).copied().unwrap_or(1.0);
            if w.abs() > 0.0 {
                p = GpPnt2d::new(p.x() / w, p.y() / w);
            }
        }
        p
    }
}

/// `BSplCLib_BuildCache` Array2 for `Dimension == 2` (`pxx:1613-1672`).
fn build_curve_cache_array_2d(
    parameter: f64,
    span_domain: f64,
    periodic: bool,
    degree: i32,
    span_index: i32,
    flat_knots: &[f64],
    poles: &[GpPnt2d],
    weights: Option<&[f64]>,
    is_rational_on_param: bool,
) -> Vec<f64> {
    let mut prep = prepare_curve_eval_2d(
        parameter,
        span_index,
        degree,
        periodic,
        poles,
        weights,
        flat_knots,
        None,
    );
    let cache_shift = if is_rational_on_param && !prep.rational {
        prep.dim + 1
    } else {
        prep.dim
    };
    bohm(
        prep.u,
        degree,
        degree,
        &prep.knots,
        prep.dim,
        &mut prep.poles,
    );
    let rows = (degree + 1) as usize;
    let cols = cache_shift as usize;
    let mut out = vec![0.0f64; rows * cols];
    let mut coeff = 1.0e0;
    let mut poly = 0usize;
    let mut cache = 0usize;
    for _i in 0..=degree {
        for j in 0..prep.dim as usize {
            if cache + j < out.len() && poly + j < prep.poles.len() {
                out[cache + j] = prep.poles[poly + j] * coeff;
            }
        }
        coeff *= span_domain / ((_i + 1) as f64);
        poly += prep.dim as usize;
        cache += prep.dim as usize;
        if cache_shift > prep.dim {
            if cache < out.len() {
                out[cache] = 0.0;
            }
            cache += 1;
        }
    }
    if cache_shift > prep.dim && !out.is_empty() {
        // `pxx:1666-1670`: weight of the constant term is 1.
        let w_col = cache_shift as usize - 1;
        if w_col < cols {
            out[w_col] = 1.0;
        }
    }
    out
}
