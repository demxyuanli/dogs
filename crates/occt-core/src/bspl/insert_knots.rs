//! `BSplCLib::PrepareInsertKnots` / `BSplCLib::InsertKnots` in the
//! **distinct knot + multiplicity** form used by `Geom_BSplineCurve`, together
//! with the Boor-scheme helpers they run on.
//!
//! Source: `BSplCLib.cxx:1010-1071` (`BoorScheme`), `:1783-1806` (`BuildBoor`),
//! `:1810-1821` (`BoorIndex`), `:1825-1845` (`GetPole`), `:1849-2024`
//! (`PrepareInsertKnots`), `:2028-2059` (file-static `Copy`), `:2063-2351`
//! (`InsertKnots`, dimension form), `BSplCLib_CurveComputation.pxx:380-438`
//! (`BSplCLib_InsertKnots`, the pole/weight wrapper) and `PLib::SetPoles` /
//! `PLib::GetPoles` (`PLib.cxx:115-201`).
//!
//! This is *not* `bspl::bezier::boehm_insert` (a single Boehm step used on a
//! flat knot window): the routine below is the general one — it takes distinct
//! knots and multiplicities, supports periodic curves (the first/last knot is
//! one knot and the poles wrap) and implements both OCCT `Add` semantics
//! ("raise the multiplicity to M" when `Add == false`, "add M" when `true`).
//!
//! Translation notes:
//! - OCCT indexes the flat pole arrays 1-based; the port keeps the same
//!   arithmetic with `i32` indices and `Lower() == 1` for every array.
//! - OCCT's `NCollection_LocalArray`/pointer walk in `BoorScheme` and
//!   `GetPole` reaches outside the "array" on purpose (strides of
//!   `2 * Dimension` and a wrap); the walk is reproduced verbatim.
//! - Where OCCT would read out of bounds on a malformed (but never produced by
//!   `Geom_BSplineCurve`) input, the port clamps the index and says so.

use crate::gp::GpPnt;
use crate::precision::{epsilon, REAL_SMALL};

use super::build_knots::build_knots;
use super::locate::{first_u_knot_index, last_u_knot_index};

/// 1-based read of `arr` with the index clamped into range.
#[inline]
fn get1(arr: &[f64], i: i32) -> f64 {
    if arr.is_empty() {
        return 0.0;
    }
    arr[(i.clamp(1, arr.len() as i32) - 1) as usize]
}

/// 1-based write of `arr` with the index clamped into range.
#[inline]
fn set1(arr: &mut [f64], i: i32, v: f64) {
    if arr.is_empty() {
        return;
    }
    let n = arr.len() as i32;
    arr[(i.clamp(1, n) - 1) as usize] = v;
}

/// `BSplCLib::PrepareInsertKnots` (`BSplCLib.cxx:1849-2024`).
///
/// Returns `(NbPoles, NbKnots)` for the curve that results from inserting
/// `add_knots` with `add_mults` (multiplicity `1` each when `add_mults` is
/// `None`, i.e. OCCT's `addflat`). `None` mirrors the `false` return that
/// `Geom_BSplineCurve::InsertKnots` turns into `Standard_ConstructionError`.
pub fn prepare_insert_knots(
    degree: i32,
    periodic: bool,
    knots: &[f64],
    mults: &[i32],
    add_knots: &[f64],
    add_mults: Option<&[i32]>,
    tolerance: f64,
    add: bool,
) -> Option<(i32, i32)> {
    if knots.is_empty() || mults.is_empty() || add_knots.is_empty() {
        return None;
    }
    let addflat = add_mults.is_none();

    // `first`/`last` are the first/last **unique** knots of the curve
    // (`cxx:1862-1872`).
    let (first, last) = if periodic {
        (1i32, knots.len() as i32)
    } else {
        (
            first_u_knot_index(degree, mults),
            last_u_knot_index(degree, mults),
        )
    };
    let a_delta_k1 = get1(knots, first) - add_knots[0];
    let a_delta_k2 = add_knots[add_knots.len() - 1] - get1(knots, last);
    if a_delta_k1 > tolerance {
        return None;
    }
    if a_delta_k2 > tolerance {
        return None;
    }

    let mut sigma = 0i32;
    let mut nb_knots = 0i32;
    let knots_upper = knots.len() as i32;
    let mut k = 0i32; // `Knots.Lower() - 1`
    let mut ak = 0usize; // `AddKnots.Lower()` (0-based)
    if periodic
        && add_knots.len() > 1
        && a_delta_k1.abs() <= REAL_SMALL
        && a_delta_k2.abs() <= REAL_SMALL
    {
        // gka: a full-period segment may add only one knot, at the end.
        ak += 1;
    }

    let mut a_last_knot_mult = mults[mults.len() - 1];
    let mut oldau = add_knots[ak];

    while ak < add_knots.len() {
        let au = add_knots[ak];
        if au < oldau {
            return None;
        }
        oldau = au;

        let eps = tolerance.max(epsilon(au));

        while k < knots_upper && get1(knots, k + 1) - au <= eps {
            k += 1;
            nb_knots += 1;
            sigma += get1_i32(mults, k);
        }

        let mut amult = if addflat {
            1
        } else {
            add_mults.unwrap()[ak].max(0)
        };

        while ak + 1 < add_knots.len() && (au - add_knots[ak + 1]).abs() <= eps {
            ak += 1;
            if add {
                if addflat {
                    amult += 1;
                } else {
                    amult += add_mults.unwrap()[ak].max(0);
                }
            }
        }

        // `k >= first >= 1` here: `PrepareInsertKnots` already rejected an
        // `AddKnots` below `Knots(first)`, so every earlier knot satisfies
        // `Knots(k + 1) - au <= tolerance <= Eps` and `k` advanced past them.
        let kk = k.max(1);
        if (au - get1(knots, kk)).abs() <= eps {
            // identical to an existing knot
            let mult = get1_i32(mults, kk);
            if add {
                if mult + amult > degree {
                    amult = (degree - mult).max(0);
                }
                sigma += amult;
            } else if amult > mult {
                if amult > degree {
                    amult = degree;
                }
                if k == knots_upper && periodic {
                    a_last_knot_mult = amult.max(mult);
                    sigma += 2 * (a_last_knot_mult - mult);
                } else {
                    sigma += amult - mult;
                }
            }
        } else if amult > 0 {
            // not identical to an existing knot
            if amult > degree {
                amult = degree;
            }
            nb_knots += 1;
            sigma += amult;
        }

        ak += 1;
    }

    while k < knots_upper {
        k += 1;
        nb_knots += 1;
        sigma += get1_i32(mults, k);
    }

    let nb_poles = if periodic {
        sigma - a_last_knot_mult
    } else {
        sigma - degree - 1
    };
    Some((nb_poles, nb_knots))
}

#[inline]
fn get1_i32(arr: &[i32], i: i32) -> i32 {
    if arr.is_empty() {
        return 0;
    }
    arr[(i.clamp(1, arr.len() as i32) - 1) as usize]
}

/// `BSplCLib::BoorIndex` (`BSplCLib.cxx:1810-1821`).
pub(super) fn boor_index(index: i32, length: i32, depth: i32) -> i32 {
    if index <= depth {
        return index;
    }
    if index <= length {
        return 2 * index - depth;
    }
    length + index - depth
}

/// `BSplCLib::BuildBoor` (`BSplCLib.cxx:1783-1806`): copies `length + 1` poles
/// of the view `[1, view_upper]` into the local Boor work array at a stride of
/// `2 * Dimension`, wrapping at `view_upper`.
fn build_boor(pole_view: &[f64], view_upper: i32, index: i32, length: i32, dim: i32, local: &mut [f64]) {
    let mut ip = 1 + index * dim;
    for i in 0..=length {
        for k in 0..dim {
            local[(i * 2 * dim + k) as usize] = get1(pole_view, ip);
            ip += 1;
            if ip > view_upper {
                ip = 1;
            }
        }
    }
}

/// `BSplCLib::BoorScheme` (`BSplCLib.cxx:1010-1071`).
///
/// `local` is the 0-based Boor work array of `(2 * Degree + 1) * Dimension`
/// doubles; `firstpole = &Poles - 2 * Dimension` (`cxx:1047`) is kept as an
/// `i64` index exactly as OCCT's pointer walk.
fn boor_scheme(u: f64, degree: i32, knots: &[f64], dim: i32, local: &mut [f64], depth: i32, length: i32) {
    let dim_i = dim as i64;
    let mut firstpole: i64 = -2 * dim_i;
    for step in 0..depth as i64 {
        firstpole += dim_i;
        let mut pole = firstpole;
        for i in step..length as i64 {
            pole += 2 * dim_i;
            let kd = (i + degree as i64 - step) as usize;
            let ki = i as usize;
            let x = (knots[kd] - u) / (knots[kd] - knots[ki]);
            let y = 1.0 - x;
            for k in 0..dim_i {
                let p = (pole + k) as usize;
                local[p] = x * local[p - dim_i as usize] + y * local[p + dim_i as usize];
            }
        }
    }
}

/// File-static `Copy` (`BSplCLib.cxx:2028-2059`) with both arrays' `Lower() == 1`.
/// `nb` counts scalars (a dimension-major flat array), not poles.
pub(super) fn copy(nb: i32, old_first: &mut i32, old: &[f64], new_first: &mut i32, new: &mut [f64]) {
    let old_len = old.len() as i32;
    let new_len = new.len() as i32;
    if old_len <= 0 || new_len <= 0 || nb <= 0 {
        return;
    }
    *old_first = 1 + (*old_first - 1).rem_euclid(old_len);
    *new_first = 1 + (*new_first - 1).rem_euclid(new_len);
    for _ in 0..nb {
        new[(*new_first - 1) as usize] = old[(*old_first - 1) as usize];
        *old_first += 1;
        if *old_first > old_len {
            *old_first = 1;
        }
        *new_first += 1;
        if *new_first > new_len {
            *new_first = 1;
        }
    }
}

/// `BSplCLib::GetPole` (`BSplCLib.cxx:1825-1845`): reads the pole
/// `BoorIndex(index, …)` out of the Boor work array into `pole` at `position`,
/// wrapping `position` at `pole_upper`.
#[allow(clippy::too_many_arguments)]
pub(super) fn get_pole(
    local: &[f64],
    index: i32,
    length: i32,
    depth: i32,
    dim: i32,
    pole: &mut [f64],
    pole_upper: i32,
    position: &mut i32,
) {
    let base = (boor_index(index, length, depth) * dim) as usize;
    for k in 0..dim as usize {
        let pos = (*position).clamp(1, pole.len() as i32);
        pole[(pos - 1) as usize] = local[base + k];
        *position += 1;
        if *position > pole_upper {
            *position = 1;
        }
    }
}

/// `BSplCLib::InsertKnots(Degree, Periodic, Dimension, Poles, Knots, Mults,
/// AddKnots, AddMults, NewPoles, NewKnots, NewMults, Tolerance, Add)`
/// (`BSplCLib.cxx:2063-2351`), the dimension-major flat form.
///
/// The output arrays must already have the sizes `PrepareInsertKnots` reports
/// (`NewPoles = Dimension * NbPoles`, `NewKnots`/`NewMults = NbKnots`), exactly
/// as OCCT requires.
#[allow(clippy::too_many_arguments)]
pub fn insert_knots_dim(
    degree: i32,
    periodic: bool,
    dim: i32,
    poles: &[f64],
    knots: &[f64],
    mults: &[i32],
    add_knots: &[f64],
    add_mults: Option<&[i32]>,
    new_poles: &mut [f64],
    new_knots: &mut [f64],
    new_mults: &mut [i32],
    tolerance: f64,
    add: bool,
) {
    let addflat = add_mults.is_none();
    let deg = degree.max(0);
    let mut local_knots: Vec<f64> = vec![0.0; (2 * deg) as usize];
    let mut local_poles: Vec<f64> = vec![0.0; ((2 * deg + 1) * dim.max(1)) as usize];

    let knots_upper = knots.len() as i32;
    let mut curk = 0i32; // `Knots.Lower() - 1`
    let mut curnk = 0i32; // `NewKnots.Lower() - 1`
    let mut curp = 1i32; // `Poles.Lower()`
    let mut curnp = 1i32; // `NewPoles.Lower()`

    let mut index = if periodic { -mults[0] } else { -degree - 1 };
    let mut firstmult = 0i32;

    for kn in 0..add_knots.len() {
        let u = add_knots[kn];
        let eps = tolerance.max(epsilon(u));

        // find the position in the old knots and copy to the new knots
        while curk < knots_upper && get1(knots, curk + 1) - u <= eps {
            curk += 1;
            curnk += 1;
            set1(new_knots, curnk, get1(knots, curk));
            let m = get1_i32(mults, curk);
            set1_i32(new_mults, curnk, m);
            index += m;
        }

        // slice the knots and mults to the current size of the new curve
        let view_upper = (curnk + knots_upper - curk).max(0);
        let view_len = (view_upper as usize).min(new_knots.len());

        // copy enough knots to compute the insertion schema
        let mut k = curk;
        let mut i = curnk;
        let mut mult = 0i32;
        while mult < degree && k < knots_upper {
            k += 1;
            i += 1;
            set1(new_knots, i, get1(knots, k));
            let m = get1_i32(mults, k);
            set1_i32(new_mults, i, m);
            mult += m;
        }

        // copy knots at the end for periodic curve
        if periodic {
            mult = 0;
            k = knots_upper;
            i = view_upper;
            while mult < degree && i > curnk {
                set1(new_knots, i, get1(knots, k));
                let m = get1_i32(mults, k);
                set1_i32(new_mults, i, m);
                mult += m;
                k -= 1;
                i -= 1;
            }
            let first = get1_i32(new_mults, 1);
            set1_i32(new_mults, view_upper, first);
        }

        // Boor scheme on the new curve to insert the new knot
        let sameknot = (u - get1(new_knots, curnk)).abs() <= eps;

        let length = if sameknot {
            (degree - get1_i32(new_mults, curnk)).max(0)
        } else {
            degree
        };

        let mut depth = if addflat { 1 } else { degree.min(add_mults.unwrap()[kn]) };

        if sameknot {
            if add {
                let m = get1_i32(new_mults, curnk);
                if m + depth > degree {
                    depth = degree - m;
                }
            } else {
                depth = (depth - get1_i32(new_mults, curnk)).max(0);
            }

            if periodic && (curk == 1 || curk == knots_upper) {
                // the first and last knot are the same knot and are delayed
                // to the end of the routine
                if firstmult == 0 {
                    firstmult += depth;
                }
                depth = 0;
            }
        }
        if depth <= 0 {
            continue;
        }

        local_knots = build_knots(
            degree,
            curnk,
            periodic,
            &new_knots[..view_len],
            Some(&new_mults[..view_len]),
        );

        // copy the poles
        let mut need = 1 + (index + length + 1) * dim - curnp;
        need = need.min(poles.len() as i32 - curp + 1);

        let mut p = curp;
        let mut np = curnp;
        copy(need, &mut p, poles, &mut np, new_poles);
        curp += need;
        curnp += need;

        // slice the poles to the current number of poles in case of periodic
        let boor_upper = curnp - 1;
        build_boor(new_poles, boor_upper, index, length, dim, &mut local_poles);
        boor_scheme(u, degree, &local_knots, dim, &mut local_poles, depth, length);

        // copy the new poles
        curnp += depth * dim;
        let the_upper = curnp - 1;
        let mut np = 1 + (index + 1) * dim;
        for i in 1..=(length + depth) {
            get_pole(
                &local_poles,
                i,
                length,
                depth,
                dim,
                new_poles,
                the_upper,
                &mut np,
            );
        }

        // insert the knot
        index += depth;
        if sameknot {
            let m = get1_i32(new_mults, curnk) + depth;
            set1_i32(new_mults, curnk, m);
        } else {
            curnk += 1;
            set1(new_knots, curnk, u);
            set1_i32(new_mults, curnk, depth);
        }
    }

    // copy the last poles and knots
    let mut p = curp;
    let mut np = curnp;
    copy(poles.len() as i32 - curp + 1, &mut p, poles, &mut np, new_poles);

    while curk < knots_upper {
        curk += 1;
        curnk += 1;
        set1(new_knots, curnk, get1(knots, curk));
        let m = get1_i32(mults, curk);
        set1_i32(new_mults, curnk, m);
    }

    // process the first-last knot on periodic curves
    if firstmult > 0 {
        let curnk = 1i32; // `NewKnots.Lower()`
        let mut firstmult = firstmult;
        if get1_i32(new_mults, curnk) + firstmult > degree {
            firstmult = degree - get1_i32(new_mults, curnk);
        }
        if firstmult > 0 {
            let length = degree - get1_i32(new_mults, curnk);
            let depth = firstmult;

            let lk = build_knots(degree, curnk, periodic, new_knots, Some(new_mults));
            let npoles_upper = (new_poles.len() as i32 - depth * dim).max(0);
            build_boor(new_poles, npoles_upper, 0, length, dim, &mut local_poles);
            boor_scheme(
                get1(new_knots, curnk),
                degree,
                &lk,
                dim,
                &mut local_poles,
                depth,
                length,
            );

            // copy the new poles but rotate them with depth
            let mut np = 1;
            for i in depth..(length + depth) {
                get_pole(
                    &local_poles,
                    i,
                    length,
                    depth,
                    dim,
                    new_poles,
                    new_poles.len() as i32,
                    &mut np,
                );
            }
            let mut np = new_poles.len() as i32 - depth * dim + 1;
            for i in 0..depth {
                get_pole(
                    &local_poles,
                    i,
                    length,
                    depth,
                    dim,
                    new_poles,
                    new_poles.len() as i32,
                    &mut np,
                );
            }

            let lower = get1_i32(new_mults, 1) + depth;
            set1_i32(new_mults, 1, lower);
            let n = new_mults.len() as i32;
            let upper = get1_i32(new_mults, n) + depth;
            set1_i32(new_mults, n, upper);
        }
    }
}

#[inline]
fn set1_i32(arr: &mut [i32], i: i32, v: i32) {
    if arr.is_empty() {
        return;
    }
    let n = arr.len() as i32;
    arr[(i.clamp(1, n) - 1) as usize] = v;
}

/// Result of [`insert_knots`]: the re-represented curve.
pub struct InsertKnotsResult {
    pub poles: Vec<GpPnt>,
    pub weights: Option<Vec<f64>>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
}

/// `BSplCLib_InsertKnots` (`BSplCLib_CurveComputation.pxx:380-438`) with
/// `PLib::SetPoles`/`PLib::GetPoles` (`PLib.cxx:115-201`): poles and (optional)
/// weights in, poles and weights out. `None` when `PrepareInsertKnots`
/// refuses the request or the result is degenerate.
#[allow(clippy::too_many_arguments)]
pub fn insert_knots(
    degree: i32,
    periodic: bool,
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    add_knots: &[f64],
    add_mults: Option<&[i32]>,
    tolerance: f64,
    add: bool,
) -> Option<InsertKnotsResult> {
    let (nb_poles, nb_knots) =
        prepare_insert_knots(degree, periodic, knots, mults, add_knots, add_mults, tolerance, add)?;
    if nb_poles <= 0 || nb_knots <= 0 {
        return None;
    }
    let nb_poles = nb_poles as usize;
    let nb_knots = nb_knots as usize;
    let rational = weights.is_some();
    let dim = if rational { 4 } else { 3 };

    // `PLib::SetPoles` (`PLib.cxx:115-156`).
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

    let mut new_fpoles = vec![0.0f64; dim * nb_poles];
    let mut new_knots = vec![0.0f64; nb_knots];
    let mut new_mults = vec![0i32; nb_knots];

    insert_knots_dim(
        degree,
        periodic,
        dim as i32,
        &fpoles,
        knots,
        mults,
        add_knots,
        add_mults,
        &mut new_fpoles,
        &mut new_knots,
        &mut new_mults,
        tolerance,
        add,
    );

    // `PLib::GetPoles` (`PLib.cxx:160-201`).
    let mut new_poles = Vec::with_capacity(nb_poles);
    let mut new_weights = if rational {
        Some(Vec::with_capacity(nb_poles))
    } else {
        None
    };
    for i in 0..nb_poles {
        let (x, y, z) = (
            new_fpoles[i * dim],
            new_fpoles[i * dim + 1],
            new_fpoles[i * dim + 2],
        );
        if let Some(w) = new_weights.as_mut() {
            let ww = new_fpoles[i * dim + 3];
            w.push(ww);
            new_poles.push(GpPnt::new(x / ww, y / ww, z / ww));
        } else {
            new_poles.push(GpPnt::new(x, y, z));
        }
    }

    Some(InsertKnotsResult {
        poles: new_poles,
        weights: new_weights,
        knots: new_knots,
        mults: new_mults,
    })
}
