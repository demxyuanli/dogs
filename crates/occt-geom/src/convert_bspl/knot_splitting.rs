// ---------------------------------------------------------------------------
// Knot splitting (GeomConvert_BSplineCurveKnotSplitting /
// BSplineSurfaceKnotSplitting).
// ---------------------------------------------------------------------------

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;

/// Distinct knot values and their multiplicities from a full knot vector.
fn distinct_knots(knots: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut vals: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for &k in knots {
        if let Some(last) = vals.last() {
            if (k - last).abs() < 1e-15 {
                *mults.last_mut().unwrap() += 1;
                continue;
            }
        }
        vals.push(k);
        mults.push(1);
    }
    (vals, mults)
}

/// Split indices for `degree`/`mults` at the requested continuity range.
///
/// Returns 0-based indices into the distinct-knot array (the OCCT
/// `splitIndexes`), starting at the first and ending at the last distinct
/// knot. Port of `GeomConvert_BSplineCurveKnotSplitting` / the surface
/// counterpart.
fn knot_splits(degree: usize, mults: &[usize], continuity_range: usize) -> Vec<usize> {
    let first = 0usize;
    let last = mults.len() - 1;
    if continuity_range == 0 {
        return vec![first, last];
    }
    let mmax = *mults.iter().max().unwrap_or(&0);
    if degree >= mmax + continuity_range {
        return vec![first, last];
    }
    let mut out = vec![first];
    for i in 1..last {
        // Continuity at knot i is (degree − multiplicity); split where it
        // drops below `continuity_range`.
        if degree < mults[i] + continuity_range {
            out.push(i);
        }
    }
    out.push(last);
    out
}

/// Split indices of `c` where its continuity drops below C1 (the OCCT default
/// `ContinuityRange = 1`). Indices are into the curve's distinct-knot array.
pub fn knot_splitting_curve(c: &GeomBSplineCurve) -> Vec<usize> {
    knot_splitting_curve_with_continuity(c, 1)
}

/// `knot_splitting_curve` with an explicit continuity range.
pub fn knot_splitting_curve_with_continuity(c: &GeomBSplineCurve, continuity_range: usize) -> Vec<usize> {
    let (_, mults) = distinct_knots(&c.knots);
    knot_splits(c.degree, &mults, continuity_range)
}

/// Split indices of `s` in u and v (each a `Vec<usize>` into the respective
/// distinct-knot array), at the OCCT default continuity range 1.
pub fn knot_splitting_surface(s: &GeomBSplineSurface) -> (Vec<usize>, Vec<usize>) {
    knot_splitting_surface_with_continuity(s, 1, 1)
}

/// `knot_splitting_surface` with explicit u/v continuity ranges.
pub fn knot_splitting_surface_with_continuity(
    s: &GeomBSplineSurface,
    u_continuity: usize,
    v_continuity: usize,
) -> (Vec<usize>, Vec<usize>) {
    let (_, um) = distinct_knots(&s.knots_u);
    let (_, vm) = distinct_knots(&s.knots_v);
    (
        knot_splits(s.deg_u, &um, u_continuity),
        knot_splits(s.deg_v, &vm, v_continuity),
    )
}
