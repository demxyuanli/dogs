//! "Equations exprimant le meme rapport de tangence sur chaque courbe"
//! (`AppParCurves_ResolConstraint.gxx:288-406`).
//!
//! For every tangency point and every adjacent pair of curves, one row ties
//! the two tangency coordinates together through the most significant
//! component chosen in `Ibont`. The index arithmetic (`--Inc3` before the
//! loop, `++Inc3` at each curve pair, `IncCol += IP`) is copied from OCCT.
//!
//! Translation notes:
//! - OCCT's `switch` has no `default`, so an `Ibont` value outside `1..=3`
//!   would read an uninitialised `T1` / `IP` / `addIndex_1`; the port panics
//!   on that unreachable arm instead of reading garbage.
//! - The 2d tangency vector of column `k` is `tabV2d(k - nb3d)`
//!   (`gxx:313`), while the neighbouring column reads `tabV2d(k + 1 - nb3d)`
//!   (`gxx:365`): both spellings are kept.

use occt_core::gp::{GpVec, GpVec2d};
use occt_math::MathMatrix;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::resol_constraint::{ibont, ResolConstraint};

use super::perform::{tangency_2d, tangency_3d};

/// `gxx:288-406`. `inc3` is the running row index left by the tangency block.
pub(super) fn fill_same_ratio<L: LineTool>(
    this: &mut ResolConstraint,
    ssp: &L,
    tab_v: &mut [GpVec],
    tab_v2d: &mut [GpVec2d],
    ibonts: &[i32],
    derivative_bern: &MathMatrix,
    nb_cu: i32,
    nb3d: i32,
    npol: i32,
    mut inc3: i32,
) {
    let inc_tan = this.inc_tan;
    // `--Inc3;` (`gxx:290`).
    inc3 -= 1;

    for i in 1..=inc_tan as i32 {
        let mut inc_col = 0;
        let npt = this.itan[(i - 1) as usize];
        for k in 1..=nb_cu - 1 {
            inc3 += 1;
            // First relation variable (`gxx:298-348`).
            let a_val = ibont(ibonts, inc_tan, k, i);
            let (t1, ip, add_index_1) = match a_val {
                1 => {
                    if k <= nb3d {
                        tangency_3d(ssp, npt, tab_v);
                        (tab_v[(k - 1) as usize].x(), 3 * npol, 0)
                    } else {
                        tangency_2d(ssp, npt, tab_v2d);
                        (tab_v2d[(k - nb3d - 1) as usize].x(), 2 * npol, 0)
                    }
                }
                2 => {
                    if k <= nb3d {
                        tangency_3d(ssp, npt, tab_v);
                        (tab_v[(k - 1) as usize].y(), 3 * npol, npol)
                    } else {
                        tangency_2d(ssp, npt, tab_v2d);
                        (tab_v2d[(k - nb3d - 1) as usize].y(), 2 * npol, npol)
                    }
                }
                3 => {
                    tangency_3d(ssp, npt, tab_v);
                    (tab_v[(k - 1) as usize].z(), 3 * npol, 2 * npol)
                }
                _ => panic!("AppParCurves_ResolConstraint: Ibont out of range"),
            };
            // Second relation variable (`gxx:350-396`).
            let a_next = ibont(ibonts, inc_tan, k + 1, i);
            let (t2, add_index_2) = match a_next {
                1 => {
                    if k + 1 <= nb3d {
                        tangency_3d(ssp, npt, tab_v);
                        (tab_v[k as usize].x(), 0)
                    } else {
                        tangency_2d(ssp, npt, tab_v2d);
                        (tab_v2d[(k + 1 - nb3d - 1) as usize].x(), 0)
                    }
                }
                2 => {
                    if k + 1 <= nb3d {
                        tangency_3d(ssp, npt, tab_v);
                        (tab_v[k as usize].y(), npol)
                    } else {
                        tangency_2d(ssp, npt, tab_v2d);
                        (tab_v2d[(k + 1 - nb3d - 1) as usize].y(), npol)
                    }
                }
                3 => {
                    tangency_3d(ssp, npt, tab_v);
                    (tab_v[k as usize].z(), 2 * npol)
                }
                _ => panic!("AppParCurves_ResolConstraint: Ibont out of range"),
            };
            // Relations between T1 and T2 (`gxx:398-405`).
            for j in 1..=npol {
                let daij = derivative_bern.value(npt as usize, j as usize);
                this.cont.set_value(
                    inc3 as usize,
                    (j + inc_col + add_index_1) as usize,
                    daij * t2,
                );
                this.cont.set_value(
                    inc3 as usize,
                    (j + ip + inc_col + add_index_2) as usize,
                    -daij * t1,
                );
            }
            inc_col += ip;
        }
    }
}




