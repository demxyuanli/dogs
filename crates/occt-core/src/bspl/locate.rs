//! `BSplCLib::Hunt` / `LocateParameter` / `ElCLib::InPeriod`.
//! Knot slices are 0-based; indices follow OCCT `Array1` Lower=1.

use crate::precision::{epsilon, Precision};

/// `BSplCLib::Hunt` (`BSplCLib.cxx:74-107`). Result is 1-based.
pub fn hunt_occt(knots: &[f64], x: f64) -> i32 {
    if knots.is_empty() {
        return 0;
    }
    if knots[0] > x {
        return 0;
    }
    if knots[knots.len() - 1] < x {
        return knots.len() as i32 + 1;
    }
    let mut pos = 1i32;
    if knots.len() <= 1 {
        return pos;
    }
    let mut hi = knots.len() as i32;
    while hi - pos != 1 {
        let mid = (hi + pos) / 2;
        if knots[(mid as usize) - 1] < x {
            pos = mid;
        } else {
            hi = mid;
        }
    }
    pos
}

/// `ElCLib::InPeriod` (`ElCLib.cxx:95-111`).
pub fn in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    if Precision::is_infinite(u) || Precision::is_infinite(u_first) || Precision::is_infinite(u_last)
    {
        return u;
    }
    let period = u_last - u_first;
    if period < epsilon(u_last) {
        return u;
    }
    u_first.max(u + period * ((u_first - u) / period).ceil())
}

fn knot_1(knots: &[f64], i: i32) -> f64 {
    knots[(i - 1) as usize]
}

/// `BSplCLib::LocateParameter(Knots, U, Periodic, FromK1, ToK2, ..., UFirst, ULast)`
/// (`BSplCLib.cxx:218-317`).
pub fn locate_parameter_range(
    knots: &[f64],
    u: f64,
    is_periodic: bool,
    from_k1: i32,
    to_k2: i32,
    u_first: f64,
    u_last: f64,
) -> (i32, f64) {
    let (first, last) = if from_k1 < to_k2 {
        (from_k1, to_k2)
    } else {
        (to_k2, from_k1)
    };
    let last1 = last - 1;
    let mut new_u = u;
    if is_periodic && (new_u < u_first || new_u > u_last) {
        new_u = in_period(new_u, u_first, u_last);
    }
    let mut knot_index = hunt_occt(knots, new_u);
    let k_upper = knots.len() as i32;
    let eps = epsilon(knots[k_upper as usize - 1].abs().min(u.abs()));
    if knot_index < k_upper {
        let mut val = new_u - knot_1(knots, knot_index + 1);
        if val < 0.0 {
            val = -val;
        }
        if val <= eps {
            knot_index += 1;
        }
    }
    if knot_index < first {
        knot_index = first;
    }
    if knot_index > last1 {
        knot_index = last1;
    }
    if knot_index != last1 {
        let mut k1 = knot_1(knots, knot_index);
        let mut k2 = knot_1(knots, knot_index + 1);
        let mut val = (k2 - k1).abs();
        while val <= eps {
            knot_index += 1;
            if knot_index >= k_upper {
                break;
            }
            k1 = k2;
            k2 = knot_1(knots, knot_index + 1);
            val = (k2 - k1).abs();
        }
    }
    (knot_index, new_u)
}

/// `BSplCLib::FirstUKnotIndex` (`BSplCLib.cxx:111-121`). `mults` is 0-based.
pub fn first_u_knot_index(degree: i32, mults: &[i32]) -> i32 {
    let mut index = 1i32;
    let mut sigma = *mults.first().unwrap_or(&0);
    while sigma <= degree {
        index += 1;
        let mi = (index - 1) as usize;
        if mi >= mults.len() {
            break;
        }
        sigma += mults[mi];
    }
    index
}

/// `BSplCLib::LastUKnotIndex` (`BSplCLib.cxx:126-137`).
pub fn last_u_knot_index(degree: i32, mults: &[i32]) -> i32 {
    let mut index = mults.len() as i32;
    let mut sigma = *mults.last().unwrap_or(&0);
    while sigma <= degree {
        index -= 1;
        if index < 1 {
            break;
        }
        sigma += mults[(index - 1) as usize];
    }
    index
}

/// `BSplCLib::LocateParameter(Degree, Knots, Mults, U, Periodic, KnotIndex, NewU)`
/// (`BSplCLib.cxx:321-364`). `knot_index` is the incoming guess (0 from EvalDN).
pub fn locate_parameter(
    degree: i32,
    knots: &[f64],
    mults: Option<&[i32]>,
    u: f64,
    periodic: bool,
    knot_index: i32,
) -> (i32, f64) {
    let lower = 1i32;
    let upper = knots.len() as i32;
    let (first, last) = if let Some(m) = mults {
        if periodic {
            (lower, upper)
        } else {
            (first_u_knot_index(degree, m), last_u_knot_index(degree, m))
        }
    } else {
        (lower + degree, upper - degree)
    };
    if knot_index < first || knot_index > last {
        locate_parameter_range(
            knots,
            u,
            periodic,
            first,
            last,
            knot_1(knots, first),
            knot_1(knots, last),
        )
    } else {
        (knot_index, u)
    }
}
