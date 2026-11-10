//! `ShapeAnalysis::AdjustToPeriod` and `fixPeriodicityTroubles` (`cxx:286-412`).

use super::super::prelude::*;
use super::*;

/// `ShapeAnalysis::AdjustToPeriod` (`ShapeAnalysis.cxx:66-69`).
pub fn adjust_to_period(val: f64, val_min: f64, val_max: f64) -> f64 {
    adjust_by_period(val, 0.5 * (val_min + val_max), val_max - val_min)
}

/// Coordinate accessor for `fixPeriodicityTroubles`: OCCT's `theIdx` is 1 for
/// U (the `x` coordinate) and 2 for V (`cxx:285`).
fn iso_coord(p: &GpPnt2d, idx: usize) -> f64 {
    if idx == 0 {
        p.x()
    } else {
        p.y()
    }
}

fn set_iso_coord(p: &mut GpPnt2d, idx: usize, v: f64) {
    if idx == 0 {
        p.set_x(v);
    } else {
        p.set_y(v);
    }
}

/// `fixPeriodicityTroubles` (`cxx:286-412`): fold the four probe points onto one
/// period window around `theSavedParam` and undo a period jump between them.
/// Returns `true` when a jump was found and fixed (`theIsRecompute`).
pub(crate) fn fix_periodicity_troubles(
    pnt: &mut [GpPnt2d],
    idx: usize,
    period: f64,
    saved_point: i32,
    saved_param: f64,
) -> bool {
    let mut min_param = 0.0;
    let mut max_param = period;
    // `cxx:295-315`: with no cached point the reference is the middle of the
    // period and the walk starts at index 0.
    let (saved_param_val, saved_point_val) = if saved_point < 0 {
        (0.5 * period, 0usize)
    } else {
        while min_param > saved_param {
            min_param -= period;
            max_param -= period;
        }
        while max_param < saved_param {
            min_param += period;
            max_param += period;
        }
        (saved_param, saved_point as usize)
    };

    // `cxx:317-324`: an iso line sits on the window border and is pinned there.
    let mut fix_iso_param = min_param;
    let mut is_iso_line = false;
    if max_param - saved_param_val < PCONFUSION || saved_param_val - min_param < PCONFUSION {
        fix_iso_param = saved_param_val;
        is_iso_line = true;
    }

    for i in 0..4 {
        let mut param = iso_coord(&pnt[i], idx);
        param += adjust_to_period(param, min_param, max_param);
        if is_iso_line {
            if max_param - param < PCONFUSION || param - min_param < PCONFUSION {
                param = fix_iso_param;
            }
        } else {
            if max_param - param < PCONFUSION {
                param = max_param;
            }
            if param - min_param < PCONFUSION {
                param = min_param;
            }
        }
        set_iso_coord(&mut pnt[i], idx, param);
    }

    // `cxx:355-371`: a monotonicity break is a period jump.
    let mut is_jump = false;
    let mut prev_diff = 0.0;
    let mut sum_diff = 1.0;
    for i in 0..3 {
        let diff = iso_coord(&pnt[i + 1], idx) - iso_coord(&pnt[i], idx);
        if diff < -PCONFUSION {
            sum_diff *= -1.0;
        }
        if diff * prev_diff < -PCONFUSION {
            is_jump = true;
        }
        prev_diff = diff;
    }
    if !is_jump {
        return false;
    }

    let sp = saved_point_val.min(3);
    if sum_diff > 0.0 {
        for i in (1..=sp).rev() {
            if iso_coord(&pnt[i], idx) > iso_coord(&pnt[i - 1], idx) {
                let v = iso_coord(&pnt[i - 1], idx) + period;
                set_iso_coord(&mut pnt[i - 1], idx, v);
            }
        }
        for i in sp..3 {
            if iso_coord(&pnt[i], idx) < iso_coord(&pnt[i + 1], idx) {
                let v = iso_coord(&pnt[i + 1], idx) - period;
                set_iso_coord(&mut pnt[i + 1], idx, v);
            }
        }
    } else {
        for i in (1..=sp).rev() {
            if iso_coord(&pnt[i], idx) < iso_coord(&pnt[i - 1], idx) {
                let v = iso_coord(&pnt[i - 1], idx) - period;
                set_iso_coord(&mut pnt[i - 1], idx, v);
            }
        }
        for i in sp..3 {
            if iso_coord(&pnt[i], idx) > iso_coord(&pnt[i + 1], idx) {
                let v = iso_coord(&pnt[i + 1], idx) + period;
                set_iso_coord(&mut pnt[i + 1], idx, v);
            }
        }
    }
    true
}
