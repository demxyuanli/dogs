//! Surface helpers from `ShapeUpgrade_UnifySameDomain.cxx`.

use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom::Surface;

/// `ClearRts` (`ShapeUpgrade_UnifySameDomain.cxx:1435-1440`).
pub(super) fn clear_rectangular_trim(surface: &Arc<dyn Surface>) -> Arc<dyn Surface> {
    surface
        .rectangular_trimmed_basis()
        .unwrap_or_else(|| surface.clone())
}

/// `SameSurf` (`ShapeUpgrade_UnifySameDomain.cxx:922-1053`).
pub(super) fn same_surface(first: &Arc<dyn Surface>, second: &Arc<dyn Surface>) -> bool {
    let (mut uf1, mut ul1) = first.u_range();
    let (mut vf1, mut vl1) = first.v_range();
    let (uf2, ul2) = second.u_range();
    let (vf2, vl2) = second.v_range();

    if !matching_lower_bound(&mut uf1, ul1, uf2)
        || !matching_lower_bound(&mut vf1, vl1, vf2)
        || !matching_upper_bound(&mut ul1, uf1, ul2)
        || !matching_upper_bound(&mut vl1, vf1, vl2)
    {
        return false;
    }

    const COEFFICIENTS: [f64; 2] = [0.3399811, 0.7745966];
    let du = ul1 - uf1;
    let dv = vl1 - vf1;
    for u_coefficient in COEFFICIENTS {
        let u = uf1 + u_coefficient * du;
        for v_coefficient in COEFFICIENTS {
            let v = vf1 + v_coefficient * dv;
            if first.d0(u, v).distance(&second.d0(u, v)) > PCONFUSION {
                return false;
            }
        }
    }
    true
}

fn matching_lower_bound(first: &mut f64, first_upper: f64, second: f64) -> bool {
    if first.is_infinite() && first.is_sign_negative() {
        if !(second.is_infinite() && second.is_sign_negative()) {
            return false;
        }
        *first = (-1.0_f64).min(first_upper - 1.0);
        true
    } else {
        second.is_finite() && (*first - second).abs() <= PCONFUSION
    }
}

fn matching_upper_bound(first: &mut f64, first_lower: f64, second: f64) -> bool {
    if first.is_infinite() && first.is_sign_positive() {
        if !(second.is_infinite() && second.is_sign_positive()) {
            return false;
        }
        *first = 1.0_f64.max(first_lower + 1.0);
        true
    } else {
        second.is_finite() && (*first - second).abs() <= PCONFUSION
    }
}
