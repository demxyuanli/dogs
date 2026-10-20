//! 3D B-spline curve. Source: `Geom_BSplineCurve.hxx`

mod curve;
mod eval;
mod knot_ops;
mod periodic;
#[cfg(test)]
mod tests;

pub use self::curve::GeomBSplineCurve;

use occt_core::bspl::knots;

/// `BSplCLib::MaxDegree()` (`BSplCLib.lxx:24-27`).
pub const MAX_DEGREE: usize = 25;

/// Recover the stored distinct knots / multiplicities from a flat knot vector
/// (`Geom_BSplineCurve::updateKnots`, `Geom_BSplineCurve.cxx:1171-1188`, via
/// `BSplCLib::KnotAnalysis`).
///
/// For a periodic curve the period-extension knots lie outside
/// `[FirstParameter, LastParameter]` (`stored_knot ± period`), so the run
/// lengths of the flat array restricted to that interval are exactly
/// `myMults`; the end knots are included because they equal `FirstParameter`
/// and `LastParameter`.
fn distinct_knots_and_mults_of(flat: &[f64], degree: usize, periodic: bool) -> (Vec<f64>, Vec<i32>) {
    let (uknots, umults) = knots::unique_knots_mults(flat);
    if !periodic || uknots.is_empty() {
        return (uknots, umults);
    }
    if flat.len() <= degree {
        // Malformed input (fewer knots than one clamped end): there is no
        // parameter window to restrict to, so keep the run lengths as they are.
        return (uknots, umults);
    }
    let (first, last) = (flat[degree], flat[flat.len() - 1 - degree]);
    let mut out_knots = Vec::new();
    let mut out_mults = Vec::new();
    for (k, m) in uknots.iter().zip(umults.iter()) {
        if *k < first || *k > last {
            continue;
        }
        out_knots.push(*k);
        out_mults.push(*m);
    }
    if out_knots.is_empty() {
        (uknots, umults)
    } else {
        (out_knots, out_mults)
    }
}
