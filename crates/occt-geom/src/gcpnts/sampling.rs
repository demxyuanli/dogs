use occt_core::gcpnts::{perform_tangential_curve, CurveSecondDeriv};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::curve::Curve;

use super::abscissa::CpntsAbscissaPoint;
use super::length::curve_length_range;

/// `n + 1` parameters at equal arc-length spacing across the curve
/// (`GCPnts_UniformAbscissa`). A zero-length curve returns `[a; n + 1]`.
pub fn uniform_abscissa(c: &dyn Curve, n: usize) -> Result<Vec<f64>, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("uniform_abscissa: unbounded curve".to_string());
    }
    let n = n.max(1);
    let total = curve_length_range(c, a, b, CONFUSION * 0.1);
    if total <= CONFUSION {
        // Degenerate: every parameter maps to the same point.
        return Ok(vec![a; n + 1]);
    }
    let step = total / n as f64;
    let resolution = c.resolution(CONFUSION);
    let mut params = Vec::with_capacity(n + 1);
    params.push(a);
    let mut prev = a;
    for _ in 1..n {
        // Find u with arc length `step` from `prev`, inside `[prev, b]`,
        // through the faithful `CPnts_AbscissaPoint::Init(C, U1, U2)` +
        // `Perform(Abscissa, U0, Resolution)` (`CPnts_AbscissaPoint.cxx:301-313`,
        // `:374-391`), i.e. `math_FunctionRoot` on the integral of the speed.
        let mut computer = CpntsAbscissaPoint::new(c);
        computer.init_range(c, prev, b);
        computer.perform(step, prev, resolution);
        if !computer.is_done() {
            return Err(
                "uniform_abscissa: math_FunctionRoot did not converge".to_string(),
            );
        }
        let u = computer.parameter().min(b);
        params.push(u);
        prev = u;
    }
    params.push(b);
    Ok(params)
}

/// `n` parameters distributed at equal *chord* intervals over a polyline
/// approximation of the curve (port of `GCPnts_QuasiUniformAbscissa`, which
/// uses a 2n-sample chord-length table).
pub fn quasi_uniform_abscissa(c: &dyn Curve, n: usize) -> Result<Vec<f64>, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("quasi_uniform_abscissa: unbounded curve".to_string());
    }
    let n = n.max(2);
    let total = curve_length_range(c, a, b, CONFUSION * 0.1);
    if total <= CONFUSION {
        return Ok(vec![a; n]);
    }
    // Chord-length vs parameter table over 2n samples (last lands on b).
    let samples = 2 * n;
    let du = (b - a) / (samples - 1) as f64;
    let mut cum = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    cum.push(0.0);
    params.push(a);
    let mut prev = c.d0(a);
    for k in 1..samples {
        let u = a + k as f64 * du;
        let p = c.d0(u);
        let l = cum[k - 1] + prev.distance(&p);
        cum.push(l);
        params.push(u);
        prev = p;
    }
    let total = *cum.last().unwrap();
    let dcorde = total / (n - 1) as f64;
    let mut out = Vec::with_capacity(n);
    out.push(a);
    let mut idx = 1usize;
    for i in 1..n - 1 {
        let target = dcorde * i as f64;
        while idx < samples - 1 && cum[idx] < target {
            idx += 1;
        }
        let denom = cum[idx] - cum[idx - 1];
        let alpha = if denom.abs() > 1e-30 {
            (target - cum[idx - 1]) / denom
        } else {
            0.0
        };
        out.push(params[idx - 1] + alpha * (params[idx] - params[idx - 1]));
    }
    out.push(b);
    Ok(out)
}

/// Adapter from `&dyn Curve` to the `GCPnts_TangentialDeflection` engine.
struct Adapter<'a>(&'a dyn Curve);

impl CurveSecondDeriv for Adapter<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.d0(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.0.d2(u)
    }
}

/// Parameters refined by both chord deviation `tol` and tangent-angle change
/// `angle_tol`, via the faithful `GCPnts_TangentialDeflection` engine
/// (`GCPnts_TangentialDeflection.cxx:522-916`, ported in
/// `occt-core/src/gcpnts_perform.rs`).
///
/// The CN breakpoint set and the BSpline/Bezier minimum-point bump follow the
/// same derivation as `meshing::edge_discret::CurveTessellator::initialize`
/// (`GCPnts_TangentialDeflection::initialize`, `cxx:415-453`); `u_tol`/`min_len`
/// are the OCCT `Initialize` defaults (`cxx:302-322`, `1.0e-9` / `CONFUSION`).
pub fn tangential_deflection(c: &dyn Curve, tol: f64, angle_tol: f64) -> Vec<f64> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return vec![a, b];
    }
    let ang = if angle_tol.is_finite() && angle_tol > 0.0 {
        angle_tol
    } else {
        std::f64::consts::PI
    };
    let mut intervals = c.parameter_intervals(6);
    if intervals.len() < 2 {
        intervals = vec![a, b];
    }
    let degree_min_nb = c.nurbs_degree().map(|d| (d + 1).max(2)).unwrap_or(2);
    let (params, _) = perform_tangential_curve(
        &Adapter(c),
        a,
        b,
        ang,
        tol,
        2,
        PCONFUSION,
        CONFUSION,
        &intervals,
        degree_min_nb,
    );
    if params.len() < 2 {
        vec![a, b]
    } else {
        params
    }
}
