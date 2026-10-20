//! Knot-vector helpers of `Geom_BSplineSurface` (`myFlatKnots` ->
//! `myUKnots`/`myUMults`), split out of the representation head.

use super::*;

impl GeomBSplineSurface {
    /// Distinct knots (`Geom_BSplineSurface` knot array, multiplicity collapsed).
    pub(super) fn unique_knots(knots: &[f64]) -> Vec<f64> {
        let mut out = Vec::new();
        for &k in knots {
            if out.last().map_or(true, |p: &f64| (k - *p).abs() > 1e-14) {
                out.push(k);
            }
        }
        out
    }
}

/// Distinct knots / multiplicities of a direction, restricted to the parameter
/// window for a periodic flat sequence (the extension knots lie outside it).
pub(super) fn distinct_km(flat: &[f64], periodic: bool, range: (f64, f64)) -> (Vec<f64>, Vec<i32>) {
    let (uks, ums) = GeomBSplineSurface::unique_knots_mults(flat);
    if !periodic {
        return (uks, ums);
    }
    let mut ok = Vec::new();
    let mut om = Vec::new();
    for (k, m) in uks.iter().zip(ums.iter()) {
        if *k < range.0 || *k > range.1 {
            continue;
        }
        ok.push(*k);
        om.push(*m);
    }
    if ok.is_empty() {
        (uks, ums)
    } else {
        (ok, om)
    }
}
