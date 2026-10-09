//! B-spline curve utility functions.
use crate::gp::GpPnt;

/// Reparameterize curve from [a,b] to [c,d].
pub fn reparameterize(knots: &mut [f64], old_a: f64, old_b: f64, new_a: f64, new_b: f64) {
    let scale = (new_b - new_a) / (old_b - old_a).max(1e-30);
    let offset = new_a - old_a * scale;
    for k in knots.iter_mut() { *k = *k * scale + offset; }
}

/// Compute total arc length via chord length approximation.
pub fn chord_length(poles: &[GpPnt]) -> f64 {
    poles.windows(2).map(|w| {
        let dx = w[1].x()-w[0].x(); let dy = w[1].y()-w[0].y(); let dz = w[1].z()-w[0].z();
        (dx*dx + dy*dy + dz*dz).sqrt()
    }).sum()
}

/// Compute parameter values by chord length parameterization.
pub fn chord_length_params(poles: &[GpPnt]) -> Vec<f64> {
    let n = poles.len(); if n < 2 { return vec![0.0]; }
    let mut params = vec![0.0f64; n];
    for i in 1..n {
        let dx = poles[i].x()-poles[i-1].x(); let dy = poles[i].y()-poles[i-1].y(); let dz = poles[i].z()-poles[i-1].z();
        params[i] = params[i-1] + (dx*dx+dy*dy+dz*dz).sqrt();
    }
    let total = params[n-1];
    if total > 1e-30 { for p in &mut params { *p /= total; } }
    params
}

/// Compute centripetal parameterization (smoother than chord length).
pub fn centripetal_params(poles: &[GpPnt]) -> Vec<f64> {
    let n = poles.len(); if n < 2 { return vec![0.0]; }
    let mut params = vec![0.0f64; n];
    for i in 1..n {
        let dx = poles[i].x()-poles[i-1].x(); let dy = poles[i].y()-poles[i-1].y(); let dz = poles[i].z()-poles[i-1].z();
        params[i] = params[i-1] + (dx*dx+dy*dy+dz*dz).sqrt().sqrt(); // sqrt of distance
    }
    let total = params[n-1];
    if total > 1e-30 { for p in &mut params { *p /= total; } }
    params
}

/// `Geom_BSplineCurve::Reverse()` (`Geom_BSplineCurve.cxx:496-511`): poles in
/// reverse order and the knot values mirrored by `BSplCLib::Reverse`
/// (`BSplCLib.cxx:802-824`) through the curve's own ends,
/// `K -> K_first + K_last - K`. `FirstParameter()`/`LastParameter()` of the
/// result are therefore `-Last`/`-First` of the original, which is what
/// `Geom_BSplineCurve::ReversedParameter` (`:520-523`, `First + Last - U`)
/// assumes; mirroring through the last knot alone (`K -> K_last - K`) shifts
/// every parameter by `First` and desynchronizes the two.
pub fn reverse_curve(poles: &mut [GpPnt], knots: &mut [f64]) {
    poles.reverse();
    let n = knots.len();
    if n < 2 {
        return;
    }
    let mirror = knots[0] + knots[n - 1];
    for i in 0..n / 2 {
        knots.swap(i, n - 1 - i);
    }
    for k in knots.iter_mut() {
        *k = mirror - *k;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn chord_length_unit() {
        let pts = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.)];
        let len = chord_length(&pts);
        assert!((len - (1.0 + 1.0)).abs() < 1e-14);
    }
    #[test]
    fn reverse_preserves_span() {
        let mut poles = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.)];
        let mut knots = vec![0.,0.,1.,1.];
        reverse_curve(&mut poles, &mut knots);
        assert!((knots[0] - 0.0).abs() < 1e-14);
        assert!((knots[3] - 1.0).abs() < 1e-14);
    }
}
