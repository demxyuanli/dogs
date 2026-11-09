//! BSplCLib::IncreaseDegree / BSplCLib::IncreaseDegreeCountKnots and the
//! pole/weight wrapper BSplCLib_IncreaseDegree.
//!
//! Source: BSplCLib.cxx:2546-2588 (IncreaseDegreeCountKnots),
//! BSplCLib.cxx:2592-2963 (dimension form, Prautzsch degree elevation,
//! CAGD 1 (1984)) and BSplCLib_CurveComputation.pxx:537-589
//! (PLib::SetPoles / PLib::GetPoles wrapper).
//!
//! Translation notes:
//! - OCCT indexes every NCollection_Array1 1-based with Lower() == 1; the
//!   port keeps that arithmetic via the get1 / set1 helpers.
//! - InsertKnots is called with Mults and NewMults aliased to the same
//!   wmults array (BSplCLib.cxx:2858-2869); the port snapshots the input to
//!   mults_in (OCCT only reads positions it has not yet overwritten, so the
//!   snapshot is behaviour-identical).
//! - The local wpoles read at NewPoles.Lower() + i * Dimension + k
//!   (BSplCLib.cxx:2809, :2817) uses the *output* array lower bound; OCCT's
//!   output is always Lower() == 1, which is what the 0-based port assumes.
//! - In the final copy (BSplCLib.cxx:2948-2953) the index k can walk past
//!   wknots/wmults for a non-clamped curve whose FirstUKnotIndex is not 1
//!   (e.g. mults [2,3,3,2], degree 3 -> 5): OCCT 8.0.0 reads the slot out of
//!   bounds (verified with DRAWEXE: the 4th knot/mult are garbage). The port
//!   clamps the index, exactly like super::insert_knots does for its own
//!   out-of-bounds walks.

use crate::gp::GpPnt;

use super::insert_knots::insert_knots_dim;
use super::locate::{first_u_knot_index, last_u_knot_index};

/// 1-based read of arr clamped into range (OCCT would read out of bounds on a
/// malformed input; the clamp keeps the port panic-free).
#[inline]
fn get1(arr: &[f64], i: i32) -> f64 {
    if arr.is_empty() {
        return 0.0;
    }
    arr[(i.clamp(1, arr.len() as i32) - 1) as usize]
}

#[inline]
fn set1(arr: &mut [f64], i: i32, v: f64) {
    if arr.is_empty() {
        return;
    }
    let n = arr.len() as i32;
    arr[(i.clamp(1, n) - 1) as usize] = v;
}

#[inline]
fn get1_i32(arr: &[i32], i: i32) -> i32 {
    if arr.is_empty() {
        return 0;
    }
    arr[(i.clamp(1, arr.len() as i32) - 1) as usize]
}

#[inline]
fn set1_i32(arr: &mut [i32], i: i32, v: i32) {
    if arr.is_empty() {
        return;
    }
    let n = arr.len() as i32;
    arr[(i.clamp(1, n) - 1) as usize] = v;
}

/// BSplCLib::IncreaseDegreeCountKnots (BSplCLib.cxx:2546-2588): the number of
/// distinct knots of the curve after raising degree to new_degree.
pub fn increase_degree_count_knots(
    degree: i32,
    new_degree: i32,
    periodic: bool,
    mults: &[i32],
) -> i32 {
    if periodic {
        return mults.len() as i32;
    }
    let degree = degree.max(0);
    let f = first_u_knot_index(degree, mults);
    let l = last_u_knot_index(degree, mults);
    let step = new_degree - degree;
    let mut removed = 0i32;

    // Lower end (BSplCLib.cxx:2559-2571).
    let mut i = 1i32;
    let mut m = degree + (f - i + 1) * step + 1;
    while m > new_degree + 1 {
        removed += 1;
        m -= get1_i32(mults, i) + step;
        i += 1;
    }
    if m < new_degree + 1 {
        removed -= 1;
    }

    // Upper end (BSplCLib.cxx:2573-2585).
    let upper = mults.len() as i32;
    let mut i = upper;
    let mut m = degree + (i - l + 1) * step + 1;
    while m > new_degree + 1 && i >= 1 {
        removed += 1;
        m -= get1_i32(mults, i) + step;
        i -= 1;
    }
    if m < new_degree + 1 {
        removed -= 1;
    }

    mults.len() as i32 - removed
}

/// BSplCLib::IncreaseDegree(Degree, NewDegree, Periodic, Dimension, Poles,
/// Knots, Mults, NewPoles, NewKnots, NewMults) (BSplCLib.cxx:2592-2963).
///
/// poles is dimension-major (Poles.Length() == Dimension * n_poles);
/// new_poles must already carry the size Geom_BSplineCurve::IncreaseDegree
/// computes as myPoles.Length() + Step * (ToK2 - FromK1)
/// (Geom_BSplineCurve.cxx:261).
#[allow(clippy::too_many_arguments)]
pub fn increase_degree_dim(
    degree: i32,
    new_degree: i32,
    periodic: bool,
    dim: i32,
    poles: &[f64],
    knots: &[f64],
    mults: &[i32],
    new_poles: &mut [f64],
    new_knots: &mut [f64],
    new_mults: &mut [i32],
) {
    let dim_us = dim.max(1) as usize;
    let degree = degree.max(0);
    let n_orig = knots.len() as i32;
    let nm = mults.len().min(knots.len());
    let f = first_u_knot_index(degree, &mults[..nm]);
    let l = last_u_knot_index(degree, &mults[..nm]);

    let mut pf = 0i32; // number of null poles added at the beginning
    let mut pl = 0i32; // number of null poles added at the end
    let mut nbwknots = n_orig;

    if periodic {
        // Periodic curves are transformed in non-periodic curves
        // (BSplCLib.cxx:2640-2661).
        nbwknots += f - 1; // f - Mults.Lower()
        pf = -degree - 1;
        for i in 1..=f {
            pf += get1_i32(mults, i);
        }
        nbwknots += n_orig - l; // Mults.Upper() - l
        pl = -degree - 1;
        for i in l..=n_orig {
            pl += get1_i32(mults, i);
        }
    }

    // Copy the knots and multiplicities (BSplCLib.cxx:2663-2697).
    let nwk = nbwknots.max(0) as usize;
    let mut wknots = vec![0.0f64; nwk];
    let mut wmults = vec![0i32; nwk];
    if !periodic {
        for i in 1..=n_orig {
            set1(&mut wknots, i, get1(knots, i));
            set1_i32(&mut wmults, i, get1_i32(mults, i));
        }
    } else if n_orig > 0 {
        let period = get1(knots, n_orig) - get1(knots, 1);
        let mut i = 0i32;
        for k in l..n_orig {
            i += 1;
            set1(&mut wknots, i, get1(knots, k) - period);
            set1_i32(&mut wmults, i, get1_i32(mults, k));
        }
        for k in 1..=n_orig {
            i += 1;
            set1(&mut wknots, i, get1(knots, k));
            set1_i32(&mut wmults, i, get1_i32(mults, k));
        }
        for k in 2..=f {
            i += 1;
            set1(&mut wknots, i, get1(knots, k) + period);
            set1_i32(&mut wmults, i, get1_i32(mults, k));
        }
    }

    // Set the first and last mults to Degree+1 and add null poles
    // (BSplCLib.cxx:2699-2705).
    pf += degree + 1 - get1_i32(&wmults, 1);
    set1_i32(&mut wmults, 1, degree + 1);
    pl += degree + 1 - get1_i32(&wmults, nbwknots);
    set1_i32(&mut wmults, nbwknots, degree + 1);

    // Poles of the working curve (BSplCLib.cxx:2711-2743).
    let mut nbwpoles = 0i32;
    for i in 1..=nbwknots {
        nbwpoles += get1_i32(&wmults, i);
    }
    nbwpoles -= degree + 1;

    let nbwp_max = nbwpoles + (nbwknots - 1) * (new_degree - degree);
    let mut wpoles = vec![0.0f64; (nbwp_max.max(0) as usize) * dim_us];

    if pf > 0 {
        let to = (pf * dim).min(wpoles.len() as i32);
        for i in 1..=to {
            set1(&mut wpoles, i, 0.0);
        }
    }
    let mut k = 1i32;
    let from = pf * dim + 1;
    let to = (nbwpoles - pl) * dim;
    if from <= to {
        for i in from..=to {
            set1(&mut wpoles, i, get1(poles, k));
            k += 1;
            if k > poles.len() as i32 {
                k = 1;
            }
        }
    }
    let from = ((nbwpoles - pl) * dim + 1).max(1);
    let to = nbwpoles * dim;
    if from <= to {
        for i in from..=to {
            set1(&mut wpoles, i, 0.0);
        }
    }

    // Loop on degree incrementation (BSplCLib.cxx:2768-2895).
    let mut nbp = nbwpoles;
    let mut nbwp = nbp;
    for cur_deg in degree..new_degree {
        nbp = nbwp; // current number of poles
        nbwp = nbp + nbwknots - 1; // new number of poles

        let mut nwpoles = vec![0.0f64; (nbwp.max(0) as usize) * dim_us];

        for step in 0..=cur_deg {
            if step != 0 {
                for i in 1..=nbwknots {
                    let v = get1_i32(&wmults, i) - 1;
                    set1_i32(&mut wmults, i, v);
                }
            }

            // Poles are the current poles but the poles congruent to step
            // are duplicated (BSplCLib.cxx:2801-2820).
            let mut tempc1 = vec![0.0f64; (nbwp.max(0) as usize) * dim_us];
            let mut offset = 0i32;
            for i in 0..nbp {
                offset += 1;
                for kk in 0..dim {
                    tempc1[((offset - 1) * dim + kk) as usize] = wpoles[(i * dim + kk) as usize];
                }
                if i % (cur_deg + 1) == step {
                    offset += 1;
                    for kk in 0..dim {
                        tempc1[((offset - 1) * dim + kk) as usize] =
                            wpoles[(i * dim + kk) as usize];
                    }
                }
            }

            // Knot multiplicities are increased (BSplCLib.cxx:2825-2844).
            let mut stepmult = step + 1;
            let mut smult = 0i32;
            let mut iknots: Vec<f64> = Vec::new();
            for kk in 1..=nbwknots {
                smult += get1_i32(&wmults, kk);
                if smult >= stepmult {
                    stepmult += cur_deg + 1;
                    let v = get1_i32(&wmults, kk) + 1;
                    set1_i32(&mut wmults, kk, v);
                } else {
                    iknots.push(get1(&wknots, kk));
                }
            }

            if !iknots.is_empty() {
                // InsertKnots(curDeg+1, false, Dimension, curve, wknots, wmults,
                // aknots, NoMults(), ncurve, nknots, wmults, 0.0) -- Add
                // defaults to true (BSplCLib.hxx:568-580).
                let mut ncurve = vec![0.0f64; (nbwp.max(0) as usize) * dim_us];
                let mut nknots_buf = vec![0.0f64; nwk];
                let mults_in = wmults.clone();
                let lim = ((offset * dim).max(0) as usize).min(tempc1.len());
                insert_knots_dim(
                    cur_deg + 1,
                    false,
                    dim,
                    &tempc1[..lim],
                    &wknots,
                    &mults_in,
                    &iknots,
                    None,
                    &mut ncurve,
                    &mut nknots_buf,
                    &mut wmults,
                    0.0,
                    true,
                );
                for i in 0..(nbwp.max(0) as usize) * dim_us {
                    nwpoles[i] += ncurve[i];
                }
            } else {
                for i in 0..(nbwp.max(0) as usize) * dim_us {
                    nwpoles[i] += tempc1[i];
                }
            }
        }

        // The result is the average (BSplCLib.cxx:2889-2894).
        for i in 0..(nbwp.max(0) as usize) * dim_us {
            wpoles[i] = nwpoles[i] / (cur_deg + 1) as f64;
        }
    }

    // Copy the results (BSplCLib.cxx:2897-2962).
    let firstknot = if periodic { n_orig - l + 1 } else { f };
    let mut m = 0i32;
    for kk in 1..=firstknot {
        m += get1_i32(&wmults, kk);
    }
    let mut k = 1i32;
    let mut pf_out = 0i32;
    while m > new_degree + 1 {
        k += 1;
        m -= get1_i32(&wmults, k);
        pf_out += get1_i32(&wmults, k);
    }
    if m < new_degree + 1 {
        k -= 1;
        let v = get1_i32(&wmults, k) + (m - new_degree - 1);
        set1_i32(&mut wmults, k, v);
        pf_out += m - new_degree - 1;
    }
    if periodic {
        k = firstknot;
    }
    let nk_out = new_knots.len() as i32;
    for i in 1..=nk_out {
        set1(new_knots, i, get1(&wknots, k));
        set1_i32(new_mults, i, get1_i32(&wmults, k));
        k += 1;
    }
    let pf_scaled = (pf_out * dim).max(0) as usize;
    for (j, out) in new_poles.iter_mut().enumerate() {
        *out = if pf_scaled + j < wpoles.len() {
            wpoles[pf_scaled + j]
        } else {
            0.0
        };
    }
}

/// Result of increase_degree: the degree-elevated curve.
pub struct IncreaseDegreeResult {
    pub poles: Vec<GpPnt>,
    pub weights: Option<Vec<f64>>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
}

/// BSplCLib_IncreaseDegree (BSplCLib_CurveComputation.pxx:537-589) with
/// PLib::SetPoles / PLib::GetPoles (PLib.cxx:115-201): poles and (optional)
/// weights in, poles and weights out.
#[allow(clippy::too_many_arguments)]
pub fn increase_degree(
    degree: i32,
    new_degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    nb_new_poles: usize,
    nb_new_knots: usize,
) -> IncreaseDegreeResult {
    let rational = weights.is_some();
    let dim = if rational { 4 } else { 3 };

    // PLib::SetPoles (PLib.cxx:115-156).
    let mut fpoles = vec![0.0f64; dim * poles.len()];
    for (i, p) in poles.iter().enumerate() {
        let w = weights.map(|w| w[i]).unwrap_or(1.0);
        fpoles[i * dim] = p.x() * w;
        fpoles[i * dim + 1] = p.y() * w;
        fpoles[i * dim + 2] = p.z() * w;
        if rational {
            fpoles[i * dim + 3] = w;
        }
    }

    let mut new_fpoles = vec![0.0f64; dim * nb_new_poles];
    let mut new_knots = vec![0.0f64; nb_new_knots];
    let mut new_mults = vec![0i32; nb_new_knots];

    increase_degree_dim(
        degree,
        new_degree,
        periodic,
        dim as i32,
        &fpoles,
        knots,
        mults,
        &mut new_fpoles,
        &mut new_knots,
        &mut new_mults,
    );

    // PLib::GetPoles (PLib.cxx:160-201).
    let mut out_poles = Vec::with_capacity(nb_new_poles);
    let mut out_weights = if rational {
        Some(Vec::with_capacity(nb_new_poles))
    } else {
        None
    };
    for i in 0..nb_new_poles {
        let (x, y, z) = (
            new_fpoles[i * dim],
            new_fpoles[i * dim + 1],
            new_fpoles[i * dim + 2],
        );
        if let Some(w) = out_weights.as_mut() {
            let ww = new_fpoles[i * dim + 3];
            w.push(ww);
            out_poles.push(GpPnt::new(x / ww, y / ww, z / ww));
        } else {
            out_poles.push(GpPnt::new(x, y, z));
        }
    }

    IncreaseDegreeResult {
        poles: out_poles,
        weights: out_weights,
        knots: new_knots,
        mults: new_mults,
    }
}
