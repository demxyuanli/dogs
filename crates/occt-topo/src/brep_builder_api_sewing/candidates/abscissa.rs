//! `GCPnts_UniformAbscissa` sampling used by `EvaluateAngulars` (GCPnts_UniformAbscissa.cxx:505-515).

use occt_core::precision::CONFUSION;
use occt_geom::gcpnts::{abscissa_point_with_tolerance, curve_length_range};
use occt_geom::Curve;




#[allow(dead_code)] // OCCT calls it only from EvaluateAngulars, whose only caller is commented out.
pub(super) fn uniform_abscissa_range(
    c: &dyn Curve,
    npt: usize,
    first: f64,
    last: f64,
) -> Option<Vec<f64>> {
    if npt < 2 {
        return Some(vec![first]);
    }
    let total = curve_length_range(c, first, last, CONFUSION);
    if total <= CONFUSION {
        return Some(vec![first; npt]);
    }
    let step = total / (npt - 1) as f64;
    let res = c.resolution(CONFUSION);
    let mut out = Vec::with_capacity(npt);
    out.push(first);
    let mut prev = first;
    for _ in 1..npt - 1 {
        let u = abscissa_point_with_tolerance(c, step, prev, res).ok()?;
        let u = u.min(last);
        out.push(u);
        prev = u;
    }
    out.push(last);
    Some(out)
}
