//! `BSplCLib::PrepareUnperiodize` / `BSplCLib::Unperiodize` on the port's
//! flat-knot representation, plus the surface wrapper `BSplSLib::Unperiodize`.
//!
//! Source: `BSplCLib.cxx:2967-3020` (`PrepareUnperiodize`), `:3024-3080`
//! (`Unperiodize`), `BSplSLib.cxx:2331-2380` (`Unperiodize`).
//!
//! OCCT keeps a periodic B-spline in the "periodic" representation: distinct
//! knots + multiplicities whose sum is `nb_poles + degree + 1`. Unperiodizing
//! raises the end multiplicities to `degree + 1`, moves the wrapped knots one
//! period outwards and extends the pole array cyclically
//! (`NewPoles(k) = Poles((k - 1) % Poles.Length() + 1)`, `BSplCLib.cxx:3076-3079`).
//!
//! `BSplSLib::Unperiodize` flattens the 2D pole array with `SetPoles`
//! (`BSplSLib.cxx:1916-1970`): row-major `(u, v)` for the U direction, and
//! column-major `(v, u)` for V. Because the flat wrap's modulus is the whole
//! array, the wrapped index is exactly `old_index mod n` **along that
//! direction** in both layouts, which is what [`unperiodize_direction`]
//! returns as a pole map (poles and weights follow the same map because OCCT
//! unperiodizes homogeneous poles).

/// Distinct knots of a periodic flat sequence, restricted to the fundamental
/// window. An open sequence (`len == nb_poles + degree + 1`) is returned whole.
fn periodic_window_knots(degree: i32, flat_knots: &[f64], nb_poles: usize) -> (Vec<f64>, Vec<i32>) {
    let (knots, mults) = distinct_knots_and_mults(flat_knots);
    if degree < 0 || nb_poles == 0 {
        return (knots, mults);
    }
    let deg = degree as usize;
    if flat_knots.len() <= deg.saturating_mul(2)
        || flat_knots.len() == nb_poles + deg + 1
    {
        return (knots, mults);
    }
    let lo = flat_knots[deg];
    let hi = flat_knots[flat_knots.len() - 1 - deg];
    let mut wk = Vec::new();
    let mut wm = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        if *k < lo || *k > hi {
            continue;
        }
        wk.push(*k);
        wm.push(*m);
    }
    if wk.len() < 2 {
        (knots, mults)
    } else {
        (wk, wm)
    }
}

/// Distinct knots and their multiplicities of a flat knot vector.
pub fn distinct_knots_and_mults(flat_knots: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for k in flat_knots {
        match knots.last() {
            Some(last) if *last == *k => {
                if let Some(m) = mults.last_mut() {
                    *m += 1;
                }
            }
            _ => {
                knots.push(*k);
                mults.push(1);
            }
        }
    }
    (knots, mults)
}

/// The flat knot vector of a distinct-knot/multiplicity pair.
pub fn flat_knots_from_mults(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut out = Vec::new();
    for (k, m) in knots.iter().zip(mults) {
        for _ in 0..(*m).max(0) {
            out.push(*k);
        }
    }
    out
}

/// `BSplCLib::PrepareUnperiodize` (`BSplCLib.cxx:2967-3020`): the sizes of the
/// unperiodized knot and pole arrays.
pub fn prepare_unperiodize(degree: i32, mults: &[i32]) -> (usize, usize) {
    let mut nb_knots = mults.len() as i32;
    let mut nb_poles = -degree - 1;
    for m in mults {
        nb_poles += *m;
    }

    // Knots added in front until the end multiplicity reaches `degree + 1`.
    let mut sigma = mults[0];
    let mut k = mults.len() as i32 - 2;
    while sigma < degree + 1 {
        sigma += mults[k as usize];
        nb_poles += mults[k as usize];
        k -= 1;
        nb_knots += 1;
    }
    if sigma > degree + 1 {
        nb_poles -= sigma - degree - 1;
    }

    // Same at the end.
    let mut sigma = mults[mults.len() - 1];
    let mut k = 1i32;
    while sigma < degree + 1 {
        sigma += mults[k as usize];
        nb_poles += mults[k as usize];
        k += 1;
        nb_knots += 1;
    }
    if sigma > degree + 1 {
        nb_poles -= sigma - degree - 1;
    }

    (nb_knots.max(0) as usize, nb_poles.max(0) as usize)
}

/// `BSplCLib::Unperiodize` (`BSplCLib.cxx:3024-3080`) on distinct
/// knots/multiplicities. Returns `(new_knots, new_mults, index)` where `index`
/// is the number of knots prepended (the knot-loop offset of the algorithm).
pub fn unperiodize_knots(degree: i32, knots: &[f64], mults: &[i32]) -> (Vec<f64>, Vec<i32>, usize) {
    let (nb_knots, _nb_poles) = prepare_unperiodize(degree, mults);
    if knots.is_empty() || mults.is_empty() || nb_knots == 0 {
        return (knots.to_vec(), mults.to_vec(), 0);
    }

    let mut index = 0usize;
    let mut sigma = mults[0];
    let mut k = mults.len() as i32 - 2;
    while sigma < degree + 1 {
        sigma += mults[k as usize];
        k -= 1;
        index += 1;
    }

    let n = knots.len();
    let period = knots[n - 1] - knots[0];
    let mut new_knots = vec![0.0; nb_knots];
    let mut new_mults = vec![0i32; nb_knots];

    // Interior: the original knots shifted by `index`.
    for i in 0..n {
        new_knots[i + index] = knots[i];
        new_mults[i + index] = mults[i];
    }
    // Starting knots: one period before the wrapped end.
    for i in 0..index {
        new_knots[i] = new_knots[i + n - 1] - period;
        new_mults[i] = new_mults[i + n - 1];
    }
    new_mults[0] -= sigma - degree - 1;

    // Ending knots: one period after the wrapped start.
    let mut sigma_end = new_mults[index + n - 1];
    for i in (n + index)..nb_knots {
        new_knots[i] = new_knots[i - n + 1] + period;
        new_mults[i] = new_mults[i - n + 1];
        sigma_end += new_mults[i - n + 1];
    }
    new_mults[nb_knots - 1] -= sigma_end - degree - 1;

    (new_knots, new_mults, index)
}

/// `BSplSLib::Unperiodize` (`BSplSLib.cxx:2331-2380`) for one direction of the
/// port's flat-knot surface representation.
///
/// Returns the unperiodized flat knot vector and the pole map:
/// `pole_map[k]` is the index of the old pole (along that direction) that
/// supplies the new pole `k`. The map is the cyclic wrap of
/// `BSplCLib::Unperiodize`, valid for both `SetPoles` layouts as explained in
/// the module header.
///
/// `nb_poles` is the pole count of the periodic representation along that
/// direction - `Poles.Length()` in OCCT, i.e.
/// `BSplCLib::NbPoles(degree, true, mults)` (`BSplCLib.cxx:392-453`), which is
/// `sum(mults) - Mults(1)`, **not** `flat_knots.len() - degree - 1`: a periodic
/// flat vector is `KnotSequenceLength` long
/// (`BSplCLib.cxx:455-474`, `sum(mults) + 2 * (degree + 1 - Mults(1))`), so the
/// open-layout identity does not hold there. `Geom_BSplineSurface::SetUPeriodic`
/// (`Geom_BSplineSurface_1.cxx:940-981`) trims the poles to `NbPoles` while
/// `updateUKnots` rebuilds the flat sequence, which is exactly this state.
pub fn unperiodize_direction(
    degree: i32,
    flat_knots: &[f64],
    nb_poles: usize,
) -> (Vec<f64>, Vec<usize>) {
    // `SetUNotPeriodic` reads `myUKnots` (`Geom_BSplineSurface_1.cxx:1238-1288`),
    // the distinct knots of the periodic window. The flat sequence also carries
    // the one-period extension (`BSplCLib::KnotSequence`, `BSplCLib.cxx:517-546`);
    // treating those copies as knots unperiodizes a second time
    // (`PrepareUnperiodize` then reports twice the poles). The window is the
    // values from `flat[degree]` through `flat[len - degree - 1]`, which is
    // where the first and last fundamental knots sit (`M1 = degree + 1 - Mults(1)`).
    let (knots, mults) = periodic_window_knots(degree, flat_knots, nb_poles);
    let (new_knots, new_mults, _index) = unperiodize_knots(degree, &knots, &mults);
    let new_flat = flat_knots_from_mults(&new_knots, &new_mults);

    let n_old = nb_poles as i64;
    let n_new = new_flat.len() as i64 - degree as i64 - 1;
    if n_old <= 0 || n_new <= 0 {
        return (flat_knots.to_vec(), Vec::new());
    }
    let pole_map = (0..n_new as usize).map(|k| k % n_old as usize).collect();
    (new_flat, pole_map)
}
