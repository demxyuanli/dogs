//! `BSplCLib::Resolution` for 3D curves (`BSplCLib.cxx:4316-4820`, dim=3).
//!
//! Upper-bounds the first derivative from the pole grid, then
//! `UTolerance = tolerance_3d / (max_derivative * degree)`.

use crate::gp::GpPnt;
use crate::precision::REAL_SMALL;

/// `BSplCLib::Resolution` (3D). `flat_knots` is the expanded knot sequence.
pub fn bspline_curve_resolution(
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    flat_knots: &[f64],
    degree: i32,
    tolerance_3d: f64,
) -> f64 {
    let num_poles = poles.len() as i32;
    // Rust-only bounds guards: OCCT (`BSplCLib.cxx:4337`) computes
    // `num_poles = FlatKnots.Length() - Deg1` and leaves `max_derivative = 0`
    // for these degenerate inputs, which falls into the `RealSmall()` branch
    // at `BSplCLib.cxx:4812-4817` -- the same result returned here.
    if num_poles < 1 || degree < 1 || flat_knots.is_empty() {
        return tolerance_3d / REAL_SMALL;
    }
    let deg1 = degree + 1;
    let deg2 = (degree << 1) + 1;
    let num_poles_flat = flat_knots.len() as i32 - deg1;
    if num_poles_flat < 2 {
        return tolerance_3d / REAL_SMALL;
    }

    let mut max_derivative = 0.0;
    if let Some(w) = weights {
        let min_weights = w.iter().copied().fold(f64::INFINITY, f64::min);
        for ii in 1..num_poles_flat {
            let ii_index = (ii.rem_euclid(num_poles)) as usize;
            let ii_minus = ((ii - 1).rem_euclid(num_poles)) as usize;
            let span = flat_knots[(ii + degree) as usize] - flat_knots[ii as usize];
            let inverse = 1.0 / span;
            let mut lower = ii - deg1;
            if lower < 0 {
                lower = 0;
            }
            let mut upper = deg2 + ii;
            if upper > num_poles_flat {
                upper = num_poles_flat;
            }
            let pi = poles[ii_index];
            let pm = poles[ii_minus];
            let wi = w[ii_index];
            let wm = w[ii_minus];
            for jj in lower..upper {
                let jj_index = (jj.rem_euclid(num_poles)) as usize;
                let pj = poles[jj_index];
                let mut value = 0.0;
                let mut factor = (pj.x() - pi.x()) * wi - (pj.x() - pm.x()) * wm;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = (pj.y() - pi.y()) * wi - (pj.y() - pm.y()) * wm;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = (pj.z() - pi.z()) * wi - (pj.z() - pm.z()) * wm;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                value *= inverse;
                if max_derivative < value {
                    max_derivative = value;
                }
            }
        }
        max_derivative /= min_weights;
    } else {
        for ii in 1..num_poles_flat {
            let ii_index = (ii.rem_euclid(num_poles)) as usize;
            let ii_minus = ((ii - 1).rem_euclid(num_poles)) as usize;
            let span = flat_knots[(ii + degree) as usize] - flat_knots[ii as usize];
            let inverse = 1.0 / span;
            let pi = poles[ii_index];
            let pm = poles[ii_minus];
            let mut value = 0.0;
            let mut factor = pi.x() - pm.x();
            if factor < 0.0 {
                factor = -factor;
            }
            value += factor;
            factor = pi.y() - pm.y();
            if factor < 0.0 {
                factor = -factor;
            }
            value += factor;
            factor = pi.z() - pm.z();
            if factor < 0.0 {
                factor = -factor;
            }
            value += factor;
            value *= inverse;
            if max_derivative < value {
                max_derivative = value;
            }
        }
    }
    max_derivative *= degree as f64;
    // `BSplCLib.cxx:4811-4818`: `max_derivative *= Degree;
    // if (max_derivative > RealSmall()) UTolerance = Tolerance3D / max_derivative;
    // else UTolerance = Tolerance3D / RealSmall();`
    if max_derivative > REAL_SMALL {
        tolerance_3d / max_derivative
    } else {
        tolerance_3d / REAL_SMALL
    }
}
