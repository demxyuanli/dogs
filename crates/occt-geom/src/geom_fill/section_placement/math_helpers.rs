//! Scalar helpers of `GeomFill_SectionPlacement` (`GeomFill_SectionPlacement.cxx`
//! `:55-109` and `:955-983`).

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, RESOLUTION};

use crate::curve::Curve;

/// `Tangente(Path, Param, P, Tang)` (`cxx:55-69`): the point and the first
/// non-vanishing derivative (up to order 11), normalized when it clears
/// `100 * gp::Resolution()`.
pub(crate) fn tangente(path: &dyn Curve, param: f64) -> (GpPnt, GpVec) {
    let (p, mut tang) = path.d1(param);
    let mut norm = tang.magnitude();

    let mut ii = 2;
    while ii < 12 && norm < CONFUSION {
        tang = path.eval_dn(param, ii);
        norm = tang.magnitude();
        ii += 1;
    }

    if norm > 100.0 * RESOLUTION {
        tang = tang.divided(norm);
    }
    (p, tang)
}

/// `Penalite(angle, dist)` (`cxx:71-97`).
pub(crate) fn penalite(angle: f64, dist: f64) -> f64 {
    let mut penal = if dist < 1.0 {
        dist.sqrt()
    } else if dist < 2.0 {
        dist * dist
    } else {
        dist + 2.0
    };

    if angle > 1.0e-3 {
        penal += 1.0 / angle - 2.0 / std::f64::consts::PI;
    } else {
        penal += 1.0e3;
    }

    penal
}

/// `EvalAngle(V1, V2)` (`cxx:100-109`).
pub(crate) fn eval_angle(v1: &GpVec, v2: &GpVec) -> f64 {
    let mut angle = v1.angle(v2);
    if angle > std::f64::consts::FRAC_PI_2 {
        angle = std::f64::consts::PI - angle;
    }
    angle
}

/// `GeomFill_SectionPlacement::Choix(dist, angle)` (`cxx:955-983`): is the
/// candidate `(dist, angle)` better than the current `(cur_dist, cur_angle)`?
pub(crate) fn choix(
    dist: f64,
    angle: f64,
    cur_dist: f64,
    cur_angle: f64,
    gabarit: f64,
) -> bool {
    let evoldist = dist - cur_dist;
    let evolangle = angle - cur_angle;

    // (1) Si la gain en distance est > que le gabarit, on prend
    if evoldist < -gabarit {
        return true;
    }

    // (2) si l'ecart en distance est de l'ordre du gabarit
    if evoldist.abs() < gabarit {
        // (2.1) si le gain en angle est important on garde
        if evolangle > 0.5 {
            return true;
        }
        // (2.2) si la variation d'angle est moderee on evalue une fonction de
        // penalite
        if penalite(angle, dist / gabarit) < penalite(cur_angle, cur_dist / gabarit) {
            return true;
        }
    }

    false
}
