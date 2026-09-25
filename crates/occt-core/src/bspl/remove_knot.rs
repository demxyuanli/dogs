//! BSplCLib::AntiBoorScheme / BSplCLib::RemoveKnot and the pole/weight
//! wrapper BSplCLib_RemoveKnot.
//!
//! Source: BSplCLib.cxx:1075-1162 (AntiBoorScheme), BSplCLib.cxx:2355-2542
//! (RemoveKnot, dimension form) and BSplCLib_CurveComputation.pxx:303-363
//! (PLib::SetPoles / PLib::GetPoles wrapper).
//!
//! The helpers BoorIndex / Copy / GetPole live next to the InsertKnots port in
//! super::insert_knots and are reused here (BSplCLib.cxx:1810-1845, :2028-2059).
//!
//! Translation notes:
//! - AntiBoorScheme's pointer wound is signed: firstpole starts at
//!   (Depth - 1) * Dimension and is decremented before use
//!   (BSplCLib.cxx:1088, :1114), so for Depth == 1 the base offset is
//!   negative while every dereference lands inside the buffer. The port keeps
//!   the signed walk exactly.
//! - The wrapper's local pole buffer is (2 * Degree + 1) * Dimension
//!   (BSplCLib.cxx:2406) and the local knot window is 4 * Degree
//!   (BSplCLib.cxx:2405).

use crate::gp::GpPnt;

use super::build_knots::build_knots;
use super::insert_knots::{boor_index, copy, get_pole};
use super::knots::pole_index;
use super::locate::{first_u_knot_index, last_u_knot_index};

#[inline]
fn get1(arr: &[f64], i: i32) -> f64 {
    if arr.is_empty() {
        return 0.0;
    }
    arr[(i.clamp(1, arr.len() as i32) - 1) as usize]
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

/// BSplCLib::AntiBoorScheme (BSplCLib.cxx:1075-1162): the Boor scheme
/// reverted. Returns false when the removal cannot meet Tolerance.
///
/// knots is the local window of length 2 * Degree; poles is the
/// (2 * Degree + 1) * Dimension work buffer built by RemoveKnot.
pub fn anti_boor_scheme(
    u: f64,
    degree: i32,
    knots: &[f64],
    dim: i32,
    poles: &mut [f64],
    depth: i32,
    length: i32,
    tolerance: f64,
) -> bool {
    let d = dim.max(1) as isize;
    let du = dim.max(1) as usize;
    let mut firstpole: isize = (depth as isize - 1) * d;

    // Special case length = 1: only verify the central point
    // (BSplCLib.cxx:1090-1107).
    if length == 1 {
        let idx = degree.max(0) as usize;
        let x = (knots[idx] - u) / (knots[idx] - knots[0]);
        let y = 1.0 - x;
        for k in 0..du {
            let a = (firstpole + k as isize) as usize;
            let b = (firstpole + (2 * d) + k as isize) as usize;
            let c = (firstpole + d + k as isize) as usize;
            let z = x * poles[a] + y * poles[b];
            if (z - poles[c]).abs() > tolerance {
                return false;
            }
        }
        return true;
    }

    // General case (BSplCLib.cxx:1112-1160).
    for step in (0..depth).rev() {
        firstpole -= d;
        let mut pole = firstpole;

        // First step from left to right.
        let mut i = step;
        while i < length - 1 {
            pole += 2 * d;
            let idx = (i + degree - step) as usize;
            let x = (knots[idx] - u) / (knots[idx] - knots[i as usize]);
            let y = 1.0 - x;
            for k in 0..du {
                let p0 = (pole + k as isize) as usize;
                let pm = (pole - d + k as isize) as usize;
                let pp = (pole + d + k as isize) as usize;
                poles[pp] = (poles[p0] - x * poles[pm]) / y;
            }
            i += 1;
        }

        // Second step from right to left; only half of the way to avoid
        // overflows (BSplCLib.cxx:1132-1159).
        pole += 4 * d;
        let half_length = (length - 1 + step) / 2;
        let mut i = length - 1;
        while i > half_length {
            pole -= 2 * d;
            let idx = (i + degree - step) as usize;
            let x = (knots[idx] - u) / (knots[idx] - knots[i as usize]);
            let y = 1.0 - x;
            for k in 0..du {
                let p0 = (pole + k as isize) as usize;
                let pp = (pole + d + k as isize) as usize;
                let pm = (pole - d + k as isize) as usize;
                let z = (poles[p0] - y * poles[pp]) / x;
                if (z - poles[pm]).abs() > tolerance {
                    return false;
                }
                poles[pm] += z;
                poles[pm] /= 2.0;
            }
            i -= 1;
        }
    }
    true
}

/// BSplCLib::RemoveKnot(Index, Mult, Degree, Periodic, Dimension, Poles,
/// Knots, Mults, NewPoles, NewKnots, NewMults, Tolerance)
/// (BSplCLib.cxx:2355-2542), the dimension-major flat form.
///
/// Returns false where OCCT returns false (the protection of
/// BSplCLib.cxx:2372-2396 or a failed AntiBoorScheme).
#[allow(clippy::too_many_arguments)]
pub fn remove_knot_dim(
    index: i32,
    mult: i32,
    degree: i32,
    periodic: bool,
    dim: i32,
    poles: &[f64],
    knots: &[f64],
    mults: &[i32],
    new_poles: &mut [f64],
    new_knots: &mut [f64],
    new_mults: &mut [i32],
    tolerance: f64,
) -> bool {
    let degree = degree.max(0);
    let d = dim.max(1);
    let mut the_index = index;

    // Protection (BSplCLib.cxx:2372-2396).
    let (first, last) = if periodic {
        (1i32, knots.len() as i32)
    } else {
        (
            first_u_knot_index(degree, mults) + 1,
            last_u_knot_index(degree, mults) - 1,
        )
    };
    if index < first {
        return false;
    }
    if index > last {
        return false;
    }
    if periodic && index == first {
        the_index = last;
    }

    let depth = get1_i32(mults, the_index) - mult;
    let length = degree - mult;

    // Local arrays (BSplCLib.cxx:2404-2406).
    let mut lk = vec![0.0f64; (4 * degree) as usize];
    let mut lp = vec![0.0f64; ((2 * degree + 1) * d) as usize];

    // Build the knots for the anti Boor scheme (BSplCLib.cxx:2408-2429).
    let kb1 = build_knots(degree, the_index - 1, periodic, knots, Some(mults));
    for i in 0..(2 * degree) as usize {
        lk[i] = kb1.get(i).copied().unwrap_or(0.0);
    }
    let mut p_index = pole_index(degree, the_index - 1, periodic, mults);
    let kb2 = build_knots(degree, the_index, periodic, knots, Some(mults));
    for i in 0..(2 * degree) as usize {
        lk[(2 * degree) as usize + i] = kb2.get(i).copied().unwrap_or(0.0);
    }
    p_index += mult;

    for i in 0..(degree - mult) {
        lk[i as usize] = lk[(i + mult) as usize];
    }
    for i in (degree - mult)..(2 * degree) {
        lk[i as usize] = lk[(2 * degree + i) as usize];
    }

    // Build the poles for the anti Boor scheme (BSplCLib.cxx:2431-2450).
    let mut p = 1 + p_index * d;
    for i in 0..=(length + depth) {
        let j = d * boor_index(i, length, depth);
        for k in 0..d {
            lp[(j + k) as usize] = get1(poles, p + k);
        }
        p += d;
        if p > poles.len() as i32 {
            p = 1;
        }
    }

    // Anti Boor scheme (BSplCLib.cxx:2452-2457).
    if !anti_boor_scheme(
        get1(knots, the_index),
        degree,
        &lk,
        d,
        &mut lp,
        depth,
        length,
        tolerance,
    ) {
        return false;
    }

    // Copy the results (BSplCLib.cxx:2459-2539).
    let mut p = 1i32;
    let mut np = 1i32;
    copy((p_index + 1) * d, &mut p, poles, &mut np, new_poles);
    for i in 1..=length {
        get_pole(&lp, i, length, 0, d, new_poles, new_poles.len() as i32, &mut np);
    }
    p += (length + depth) * d;
    if p != 1 {
        let n = poles.len() as i32 - p + 1;
        copy(n, &mut p, poles, &mut np, new_poles);
    }

    // Knots and mults (BSplCLib.cxx:2489-2538).
    if mult > 0 {
        let n = knots.len().min(new_knots.len());
        new_knots[..n].copy_from_slice(&knots[..n]);
        let n = mults.len().min(new_mults.len());
        new_mults[..n].copy_from_slice(&mults[..n]);
        set1_i32(new_mults, the_index, mult);
        if periodic {
            if the_index == first {
                set1_i32(new_mults, last, mult);
            }
            if the_index == last {
                set1_i32(new_mults, first, mult);
            }
        }
    } else if !periodic || (the_index != first && the_index != last) {
        for i in 1..the_index {
            let (s, d2) = ((i - 1) as usize, (i - 1) as usize);
            if s < knots.len() && d2 < new_knots.len() {
                new_knots[d2] = knots[s];
            }
            if s < mults.len() && d2 < new_mults.len() {
                new_mults[d2] = mults[s];
            }
        }
        for i in (the_index + 1)..=knots.len() as i32 {
            let (s, d2) = ((i - 1) as usize, (i - 2) as usize);
            if s < knots.len() && d2 < new_knots.len() {
                new_knots[d2] = knots[s];
            }
            if s < mults.len() && d2 < new_mults.len() {
                new_mults[d2] = mults[s];
            }
        }
    } else {
        // The interesting case of a periodic curve where the first and last
        // knot is removed.
        for i in first..(last - 1) {
            let (s, d2) = (i as usize, (i - 1) as usize);
            if s < knots.len() && d2 < new_knots.len() {
                new_knots[d2] = knots[s];
            }
            if s < mults.len() && d2 < new_mults.len() {
                new_mults[d2] = mults[s];
            }
        }
        let li = (last - 2) as usize;
        let fi = (first - 1) as usize;
        if li < new_knots.len() && fi < new_knots.len() {
            new_knots[li] = new_knots[fi] + get1(knots, last) - get1(knots, first);
        }
        if li < new_mults.len() && fi < new_mults.len() {
            new_mults[li] = new_mults[fi];
        }
    }

    true
}

/// Result of remove_knot: the curve with the knot multiplicity lowered.
pub struct RemoveKnotResult {
    pub poles: Vec<GpPnt>,
    pub weights: Option<Vec<f64>>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
}

/// BSplCLib_RemoveKnot (BSplCLib_CurveComputation.pxx:303-363) with
/// PLib::SetPoles / PLib::GetPoles (PLib.cxx:115-201). None mirrors OCCT's
/// false return.
#[allow(clippy::too_many_arguments)]
pub fn remove_knot(
    index: i32,
    mult: i32,
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    nb_new_poles: usize,
    nb_new_knots: usize,
    tolerance: f64,
) -> Option<RemoveKnotResult> {
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

    if !remove_knot_dim(
        index,
        mult,
        degree,
        periodic,
        dim as i32,
        &fpoles,
        knots,
        mults,
        &mut new_fpoles,
        &mut new_knots,
        &mut new_mults,
        tolerance,
    ) {
        return None;
    }

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

    Some(RemoveKnotResult {
        poles: out_poles,
        weights: out_weights,
        knots: new_knots,
        mults: new_mults,
    })
}
