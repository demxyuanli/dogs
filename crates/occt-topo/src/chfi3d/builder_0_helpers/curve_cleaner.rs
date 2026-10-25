//! `CurveCleaner` — file static of `ChFi3d_Builder_0.cxx:3659-3682`
//! (OCCT 8.0.0).
//!
//! Makes a B-spline "as continued as possible at a given tolerance": from
//! `Mult = Degree()` down to `MultMin + 1`, it halves the tolerance and tries
//! to lower every knot of that exact multiplicity by one
//! (`RemoveKnot(ii, Mult - 1, tol)`).

use occt_geom::bspline_curve::GeomBSplineCurve;

/// `CurveCleaner(BS, Tol, MultMin)` (`cxx:3663-3682`).
///
/// `NbK` is captured once before the multiplicity loop (`cxx:3669`), so the
/// inner scan always starts at the knot count the curve had on entry. When a
/// `RemoveKnot(.., 0, ..)` drops a knot the array shrinks under that captured
/// bound; OCCT would then read `Multiplicity(ii)` out of range, the port skips
/// such indices instead (the removal itself is already done).
pub fn curve_cleaner(bs: &mut GeomBSplineCurve, tol: f64, mult_min: i32) {
    let mut tol = tol;
    // `const int NbK = BS->NbKnots();` (`cxx:3669`).
    let nb_k = bs.distinct_knots().len() as i32;

    let mut mult = bs.degree() as i32;
    while mult > mult_min {
        // `tol *= 0.5; // Progressive reduction` (`cxx:3672`).
        tol *= 0.5;
        let mut ii = nb_k;
        while ii > 1 {
            let current = (ii - 1) as usize;
            let multiplicity = match bs.multiplicities().get(current) {
                Some(m) => *m,
                None => {
                    ii -= 1;
                    continue;
                }
            };
            if multiplicity == mult {
                // `BS->RemoveKnot(ii, Mult - 1, tol);` (`cxx:3677`): the result
                // is not read by OCCT.
                let _ = bs.remove_knot(ii, mult - 1, tol);
            }
            ii -= 1;
        }
        mult -= 1;
    }
}
