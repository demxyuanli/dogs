//! `BSplCLib::Intervals` and `GeomAdaptor_Curve::LocalContinuity`.
//!
//! Source: `BSplCLib.cxx:4824-4961`, `GeomAdaptor_Curve.cxx:143-225`.

use crate::precision::PCONFUSION;

/// Compress a flat knot vector. Same grouping as `knots::unique_knots_mults`.
fn unique_of(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    crate::bspl::knots::unique_knots_mults(flat)
}

/// 0-based left knot of `u` in a unique knot sequence.
fn locate_left(knots: &[f64], u: f64) -> usize {
    if knots.is_empty() {
        return 0;
    }
    if u <= knots[0] {
        return 0;
    }
    let n = knots.len();
    if u >= knots[n - 1] {
        return n - 1;
    }
    crate::bspl::knots::hunt(knots, u)
}

/// `GeomAdaptor_Curve::LocalContinuity` (`cxx:143-225`).
///
/// Continuity on `(u1, u2)` is `C(degree - max interior multiplicity)`.
/// No interior knot (one span) is `CN`. Encoded as `GeomAbs_Shape`
/// (`cxx:205-224`): C0=0, C1=2, C2=4, C3=5, CN=6.
pub fn local_continuity(flat_knots: &[f64], degree: usize, periodic: bool, u1: f64, u2: f64) -> u8 {
    let (tk, tm) = unique_of(flat_knots);
    let nb = tk.len();
    if nb < 2 {
        return 6;
    }
    let mut index1 = locate_left(&tk, u1);
    let mut index2 = locate_left(&tk, u2);
    if index1 + 1 < nb && (tk[index1 + 1] - u1).abs() < PCONFUSION && index1 + 1 < nb {
        index1 += 1;
    }
    if (tk[index2] - u2).abs() < PCONFUSION {
        index2 = index2.saturating_sub(1);
    }
    if periodic && index1 + 1 == nb {
        index1 = 0;
    }
    let mult_max = if (index2 as isize) - (index1 as isize) <= 0 && !periodic {
        100
    } else {
        let lo = (index1 + 1).min(nb.saturating_sub(1));
        let hi = index2.min(nb.saturating_sub(1));
        let mut m = tm.get(lo).copied().unwrap_or(0);
        for i in lo..=hi {
            m = m.max(tm[i]);
        }
        degree as i32 - m
    };
    if mult_max <= 0 {
        0
    } else if mult_max == 1 {
        2
    } else if mult_max == 2 {
        4
    } else if mult_max == 3 {
        5
    } else {
        6
    }
}

/// `BSplCLib::Intervals` (`cxx:4824-4961`) plus the GeomAdaptor CN mapping:
/// requested `CN` (`>= 6`) uses `aCont = degree`.
pub fn intervals(
    flat_knots: &[f64],
    degree: usize,
    periodic: bool,
    continuity: u8,
    first: f64,
    last: f64,
    tolerance: f64,
) -> Vec<f64> {
    let (first, last) = if first <= last {
        (first, last)
    } else {
        (last, first)
    };
    let (tk, tm) = unique_of(flat_knots);
    if tk.len() < 2 {
        return vec![first, last];
    }
    let a_cont = if continuity >= 6 {
        degree as i32
    } else {
        continuity as i32
    };
    let thresh = degree as i32 - a_cont;
    let mut kept = Vec::new();
    for i in 0..tk.len() {
        if tm[i] > thresh || i == 0 || i + 1 == tk.len() {
            kept.push(tk[i]);
        }
    }
    if kept.len() < 2 {
        return vec![first, last];
    }
    let n = kept.len();
    let (mut cur_first, mut cur_last) = (first, last);
    let mut first_period = 0i32;
    let mut last_period = 0i32;
    let period = if periodic {
        let lower = tk[0];
        let upper = tk[tk.len() - 1];
        let p = upper - lower;
        if p.abs() > 0.0 {
            while cur_first < lower {
                cur_first += p;
                first_period -= 1;
            }
            while cur_last < lower {
                cur_last += p;
                last_period -= 1;
            }
            while cur_first >= upper {
                cur_first -= p;
                first_period += 1;
            }
            while cur_last >= upper {
                cur_last -= p;
                last_period += 1;
            }
            p
        } else {
            0.0
        }
    } else {
        0.0
    };

    let mut index1 = locate_left(&kept, cur_first);
    let mut index2 = locate_left(&kept, cur_last);
    if index1 + 1 < n && (kept[index1 + 1] - cur_first).abs() < tolerance {
        index1 += 1;
    }
    if n > 0 && (kept[index2] - cur_last).abs() < tolerance {
        index2 = index2.saturating_sub(1);
    }
    let period_extra = if periodic {
        (last_period - first_period) * (n as i32 - 1)
    } else {
        0
    };
    let nb = ((index2 as i32) - (index1 as i32) + 1 + period_extra).max(1) as usize;

    let mut out = Vec::with_capacity(nb + 1);
    if periodic && last_period != first_period && period != 0.0 {
        for i in index1..n {
            out.push(kept[i] + first_period as f64 * period);
        }
        let mut pnum = first_period + 1;
        while pnum < last_period {
            for i in 0..n.saturating_sub(1) {
                out.push(kept[i] + pnum as f64 * period);
            }
            pnum += 1;
        }
        for i in 0..=index2 {
            out.push(kept[i] + last_period as f64 * period);
        }
    } else {
        let hi = index2.min(n.saturating_sub(1));
        if index1 <= hi {
            for i in index1..=hi {
                out.push(kept[i] + first_period as f64 * period);
            }
        }
    }
    if out.is_empty() {
        return vec![first, last];
    }
    out[0] = first;
    if out.len() == nb + 1 {
        out[nb] = last;
    } else {
        out.push(last);
    }
    out
}

/// `GeomAdaptor_Curve::NbIntervals` / `Intervals` on `[first, last]` (`cxx:371-413`).
///
/// `tolerance` is `min(Resolution(Confusion), PConfusion)` (`cxx:402`, `cxx:499`).
pub fn adaptor_intervals(
    flat_knots: &[f64],
    degree: usize,
    periodic: bool,
    continuity: u8,
    first: f64,
    last: f64,
    tolerance: f64,
) -> Vec<f64> {
    if continuity == 0
        || (!periodic && continuity <= local_continuity(flat_knots, degree, periodic, first, last))
    {
        return vec![first, last];
    }
    let eps = if tolerance > 0.0 { tolerance } else { PCONFUSION };
    intervals(flat_knots, degree, periodic, continuity, first, last, eps)
}
