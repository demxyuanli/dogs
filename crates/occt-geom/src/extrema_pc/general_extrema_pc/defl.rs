//! `Extrema_CurveTool::DeflCurvIntervals` (`Extrema_CurveTool.cxx:41-91`).

use super::*;

/// Adapter for `GCPnts_TangentialDeflection::PerformCurve`.
pub(crate) struct CurveTd<'a>(&'a dyn Curve);

impl CurveSecondDeriv for CurveTd<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.d0(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.0.d2(u)
    }
}

/// `Extrema_CurveTool::DeflCurvIntervals` (`Extrema_CurveTool.cxx:41-91`).
///
/// Builds curvature-deflection sample breakpoints via
/// `GCPnts_TangentialDeflection` when a single C2 span would otherwise leave
/// `[First, Last]` as the only ExtPC interval.
pub(crate) fn defl_curv_intervals(curve: &dyn Curve) -> Vec<f64> {
    let tf = curve.first_parameter();
    let tl = curve.last_parameter();
    if !(tf.is_finite() && tl.is_finite()) || (tl - tf).abs() <= 0.0 {
        return vec![tf, tl];
    }

    // Chord-sum from First to each of 22 interior/end samples (`cxx:51-57`).
    // OCCT does not update the running point — distances are from `Value(tf)`.
    let p = curve.d0(tf);
    let mut length = 0.0;
    let nbpnts = DEFL_LENGTH_SAMPLES;
    for i in 2..=nbpnts {
        let t = (tf * (nbpnts - i) as f64 + (i - 1) as f64 * tl) / (nbpnts - 1) as f64;
        length += p.distance(&curve.d0(t));
    }

    let span = tl - tf;
    let dldt = length / span;
    if length <= CONFUSION || dldt < EPS_DL_DT || span > 10000.0 {
        return vec![tf, tl];
    }

    let defl = (0.01 * length / (2.0 * std::f64::consts::PI)).max(MIN_DEFL);
    if defl > MAX_DEFL {
        return vec![tf, tl];
    }
    let min_len = (0.00001 * length).max(CONFUSION);
    let tol = (0.00001 * span).max(PCONFUSION);

    // `GCPnts_TangentialDeflection(C, PI/6, aDefl, 2, aTol, aMinLen)` then
    // `PerformCurve` CN intervals (`TangentialDeflection.cxx:554-556`).
    let cn = {
        let mut iv = curve.parameter_intervals(6);
        if iv.len() < 2 {
            iv = vec![tf, tl];
        }
        iv
    };
    // Offset / OtherCurve stay at MinimumOfPoints=2; BSpline/Bezier bump is
    // inside PerformCurve when GetType matches (`cxx:568-578`).
    let degree_min_nb = curve
        .nurbs_degree()
        .map(|d| (d as usize + 1).max(2))
        .unwrap_or(2);

    let (params, _) = perform_tangential_curve(
        &CurveTd(curve),
        tf,
        tl,
        std::f64::consts::PI / 6.0,
        defl,
        2,
        tol,
        min_len,
        &cn,
        degree_min_nb,
    );
    if params.len() < 2 {
        vec![tf, tl]
    } else {
        params
    }
}

