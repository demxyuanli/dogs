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
    mirror_flat_knots(knots);
}

/// `BSplCLib::Reverse(Knots)` (`BSplCLib.cxx:802-824`).
///
/// `Knots(Lower)` and the original `Knots(Upper)` are the range ends. The
/// walk rewrites every later knot from the reversed successive gaps and never
/// stores into the first slot, so a first knot that does not survive
/// `K_first + K_last` in f64 stays put. Multiplicities are reversed by the
/// caller (`BSplCLib::Reverse(Mults)`, `cxx:828`) before the flat sequence is
/// rebuilt.
pub fn reverse_distinct_knots(knots: &mut [f64]) {
    if knots.len() < 2 {
        return;
    }
    let mut first: i32 = 0;
    let mut last: i32 = knots.len() as i32 - 1;
    let mut kfirst = knots[0];
    let mut klast = knots[last as usize];
    let mut tfirst = kfirst;
    let mut tlast = klast;
    first += 1;
    last -= 1;
    while first <= last {
        let fi = first as usize;
        let li = last as usize;
        tfirst += klast - knots[li];
        tlast -= knots[fi] - kfirst;
        kfirst = knots[fi];
        klast = knots[li];
        knots[fi] = tfirst;
        knots[li] = tlast;
        first += 1;
        last -= 1;
    }
}

/// `BSplCLib_Reverse` (`BSplCLib_CurveComputation.pxx:275-288`) for a periodic
/// curve. `Geom_BSplineCurve::Reverse` (`Geom_BSplineCurve.cxx:502-510`) passes
/// `theL = myFlatKnots.Upper() - myDeg - 1`, not the last pole. When the
/// periodic end multiplicity equals the degree that window is a single pole,
/// so the seam pole stays put and only the rest of the array is reversed.
pub fn reverse_periodic_span<T>(items: &mut [T], flat_knot_len: usize, degree: usize) {
    let n = items.len() as i32;
    if n <= 0 {
        return;
    }
    let last = flat_knot_len as i32 - degree as i32 - 1;
    let a_l = 1 + (last - 1).rem_euclid(n);
    let a = a_l as usize;
    items[..a].reverse();
    if a < items.len() {
        items[a..].reverse();
    }
}

/// Affine reflection of a flat knot vector, `K -> K_first + K_last - K`.
///
/// In exact arithmetic this matches [`reverse_distinct_knots`]. In f64,
/// `K_first + K_last` drops a first knot smaller than half an ulp of the
/// last knot, so the reflected value is `0` (or one ulp). `BSplCLib::Reverse`
/// never writes `Knots(Lower)` and keeps that knot. Callers that reverse a
/// curve must use [`reverse_distinct_knots`] on the distinct knots.
pub fn mirror_flat_knots(knots: &mut [f64]) {
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
