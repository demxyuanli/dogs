//! Periodic-seam fix-ups of `IntPatch_PrmPrmIntersection`.
//!
//! Source: `IntPatch/IntPatch_PrmPrmIntersection.cxx`:
//! `AdjustOnPeriodic` (`:2329-2448`) and `MakeNewPoint` (`:2452-2483`).

use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::PntOn2S;

/// `AdjustOnPeriodic(Surf1, Surf2, aSLin)` (`cxx:2329-2448`).
///
/// Moves the first and the last point of every long walking line off the
/// periodic seam: a parameter sitting exactly on `0` or on the period is
/// rewritten to the other end when the neighbour point is more than a quarter
/// period away.
pub(crate) fn adjust_on_periodic(
    surf1: &dyn Surface,
    surf2: &dyn Surface,
    a_slin: &mut [GeomIntLine],
) {
    let a_eps = CONFUSION;

    let mut periodic = [false; 4];
    let mut period = [0.0f64; 4];
    periodic[0] = surf1.is_u_periodic();
    periodic[1] = surf1.is_v_periodic();
    periodic[2] = surf2.is_u_periodic();
    periodic[3] = surf2.is_v_periodic();
    if periodic[0] {
        period[0] = surf1.u_period();
    }
    if periodic[1] {
        period[1] = surf1.v_period();
    }
    if periodic[2] {
        period[2] = surf2.u_period();
    }
    if periodic[3] {
        period[3] = surf2.v_period();
    }

    let mut d_period = [0.0f64; 4];
    for k in 0..4 {
        d_period[k] = 0.25 * period[k];
    }

    for line in a_slin.iter_mut() {
        let Some(a_l) = line.as_wline_mut() else {
            continue;
        };
        let a_nb_px = a_l.nb_pnts();
        if a_nb_px < 10 {
            continue;
        }
        for j in 0..2 {
            let mut modified = false;
            let (a_indx, a_indq) = if j == 0 {
                (1, 2)
            } else {
                (a_nb_px, a_nb_px - 1)
            };
            let psx = *a_l.point(a_indx);
            let psq = *a_l.point(a_indq);

            let mut ux = [psx.u1, psx.v1, psx.u2, psx.v2];
            let uq = [psq.u1, psq.v1, psq.u2, psq.v2];
            for k in 0..4 {
                let mut is_null = false;
                let mut is_period = false;
                if !periodic[k] {
                    continue;
                }
                if ux[k].abs() < a_eps {
                    modified = true;
                    is_null = true;
                } else if (ux[k] - period[k]).abs() < a_eps {
                    modified = true;
                    is_period = true;
                }
                if modified {
                    let du = (ux[k] - uq[k]).abs();
                    if du > d_period[k] {
                        if is_null {
                            ux[k] = period[k];
                        }
                        if is_period {
                            ux[k] = 0.0;
                        }
                    }
                }
            }
            if modified {
                let mut a_pnt_on_2s = psx;
                a_pnt_on_2s.set_parameters(ux[0], ux[1], ux[2], ux[3]);
                *a_l.point_mut(a_indx) = a_pnt_on_2s;
            }
        }
    }
}

/// `MakeNewPoint(replacePnt, oldPnt, Periods)` (`cxx:2452-2483`).
///
/// Takes the 3D point of `replace_pnt` and the parameters of `replace_pnt`
/// shifted by one period when they are further than half a period from
/// `old_pnt` (the neighbour already on the line), so the parameter stays
/// continuous along the line.
pub(crate) fn make_new_point(
    replace_pnt: &PntOn2S,
    old_pnt: &PntOn2S,
    periods: &[f64; 4],
) -> PntOn2S {
    let mut new_point = PntOn2S {
        p: replace_pnt.p,
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
    };

    let old_params = [old_pnt.u1, old_pnt.v1, old_pnt.u2, old_pnt.v2];
    let mut new_params = [replace_pnt.u1, replace_pnt.v1, replace_pnt.u2, replace_pnt.v2];

    for i in 0..4 {
        if periods[i] != 0.0 && (new_params[i] - old_params[i]).abs() >= 0.5 * periods[i] {
            if new_params[i] < old_params[i] {
                new_params[i] += periods[i];
            } else {
                new_params[i] -= periods[i];
            }
        }
    }

    new_point.set_parameters(new_params[0], new_params[1], new_params[2], new_params[3]);
    new_point
}
