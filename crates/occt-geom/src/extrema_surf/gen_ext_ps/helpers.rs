use super::super::prelude::*;

// ---------------------------------------------------------------------------
// Helpers (Extrema_GenExtPS.cxx:321-369, :755-854)
// ---------------------------------------------------------------------------

/// The distinct knot values of a flat knot vector - Geom_BSplineSurface::UKnots()
/// / Geom_BSplineCurve::Knots(), which OCCT stores as distinct knots plus
/// multiplicities (Geom_BSplineSurface.cxx:1148-1165). The port's surfaces keep
/// the flat sequence, so collapsing equal neighbors recovers the OCCT array.
pub(super) fn distinct_knots(flat: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in flat {
        if out.last() != Some(&k) {
            out.push(k);
        }
    }
    out
}

/// fillParams (Extrema_GenExtPS.cxx:321-369): parametric samples derived from
/// the knot spans of a B-spline / Bezier geometry. Returns None where OCCT
/// leaves theParams null (too few points, cxx:360-363).
pub(super) fn fill_params(
    knots: &[f64],
    degree: i32,
    par_min: f64,
    par_max: f64,
    sample: i32,
) -> Option<Vec<f64>> {
    let mut params: Vec<f64> = Vec::new();
    let mut prev_par = par_min;
    params.push(prev_par);
    let mut i = 1usize;
    while i < knots.len() && knots[i - 1] < (par_max - PCONFUSION) {
        if knots[i] < par_min + PCONFUSION {
            i += 1;
            continue;
        }
        let step = (knots[i] - knots[i - 1]) / degree.max(2) as f64;
        let mut k = 1;
        while k <= degree {
            let par = knots[i - 1] + k as f64 * step;
            if par > par_max - PCONFUSION {
                break;
            }
            if par > prev_par + PCONFUSION {
                params.push(par);
                prev_par = par;
            }
            k += 1;
        }
        i += 1;
    }
    params.push(par_max);
    let nb_par = params.len();
    if nb_par < sample as usize {
        return None;
    }
    Some(params)
}

/// LengthOfIso (Extrema_GenExtPS.cxx:755-791): the polyline length of the iso
/// line theIso at parameter thePar, sampled in theNbPnts points from thePar1 to
/// thePar2.
fn length_of_iso(
    s: &dyn Surface,
    iso_u: bool,
    par1: f64,
    par2: f64,
    nb_pnts: i32,
    par: f64,
) -> f64 {
    let mut len = 0.0;
    let dpar = (par2 - par1) / (nb_pnts - 1) as f64;
    let mut par_val = par1 + dpar;
    let mut p1 = if iso_u { s.d0(par, par1) } else { s.d0(par1, par) };
    for _i in 2..=nb_pnts {
        let p2 = if iso_u { s.d0(par, par_val) } else { s.d0(par_val, par) };
        len += p1.distance(&p2);
        p1 = p2;
        par_val += dpar;
    }
    len
}

/// CorrectNbSamples (Extrema_GenExtPS.cxx:793-854). Faithful, including the
/// cxx:846-853 branch which scales theNbV under the aRatio < 0.1 test as the
/// OCCT source does.
pub(super) fn correct_nb_samples(
    s: &dyn Surface,
    u1: f64,
    u2: f64,
    nb_u: &mut i32,
    v1: f64,
    v2: f64,
    nb_v: &mut i32,
) {
    let min_len = 1.0e-3;
    let nbp = 23.min(*nb_v);
    let mut len_u1 = length_of_iso(s, true, v1, v2, nbp, u1);
    if len_u1 <= min_len {
        let l = length_of_iso(s, true, v1, v2, nbp, 0.7 * u1 + 0.3 * u2);
        len_u1 = l.max(len_u1);
    }
    let mut len_u2 = length_of_iso(s, true, v1, v2, nbp, u2);
    if len_u2 <= min_len {
        let l = length_of_iso(s, true, v1, v2, nbp, 0.3 * u1 + 0.7 * u2);
        len_u2 = l.max(len_u2);
    }
    let nbp = 23.min(*nb_v);
    let mut len_v1 = length_of_iso(s, false, u1, u2, nbp, v1);
    if len_v1 <= min_len {
        let l = length_of_iso(s, false, u1, u2, nbp, 0.7 * v1 + 0.3 * v2);
        len_v1 = l.max(len_v1);
    }
    let mut len_v2 = length_of_iso(s, false, u1, u2, nbp, v2);
    if len_v2 <= min_len {
        let l = length_of_iso(s, false, u1, u2, nbp, 0.3 * v1 + 0.7 * v2);
        len_v2 = l.max(len_v2);
    }

    let step_v1 = len_u1 / *nb_v as f64;
    let step_v2 = len_u2 / *nb_v as f64;
    let step_u1 = len_v1 / *nb_u as f64;
    let step_u2 = len_v2 / *nb_u as f64;

    let max_step_v = step_v1.max(step_v2);
    let max_step_u = step_u1.max(step_u2);

    let ratio = max_step_v / max_step_u;
    if ratio > 10.0 {
        let mult = real_to_int(ratio.ln());
        if mult > 1 {
            *nb_v *= mult;
        }
    } else if ratio < 0.1 {
        let mult = real_to_int(-ratio.ln());
        if mult > 1 {
            *nb_v *= mult;
        }
    }
}

/// RealToInt (Standard_Real.hxx:319-328): truncation toward zero with the int
/// range clamp.
fn real_to_int(v: f64) -> i32 {
    if v < i32::MIN as f64 {
        i32::MIN
    } else if v > i32::MAX as f64 {
        i32::MAX
    } else {
        v as i32
    }
}
