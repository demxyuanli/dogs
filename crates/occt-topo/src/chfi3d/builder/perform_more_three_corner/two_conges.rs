//! The `deuxconges` / `deuxcgnontg` / `isOnSameDiff` flags of
//! `PerformMoreThreeCorner`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:1583-1676`
//! (OCCT 8.0.0). The `PerformOneCorner(Jndex, Standard_True)` call of
//! `cxx:1671-1675` stays in the builder entry point, because it mutates the
//! builder and not the arrays.
//!
//! Notes:
//! - `trouve` is the out-parameter of `ChFi3d_cherche_vertex` in the topology
//!   stage and OCCT reuses it here as the "stop the loop" flag (`cxx:1587`,
//!   `cxx:1593`). The port keeps a fresh local: the two uses never overlap,
//!   because `trouve` is reset at `cxx:1587` before the scan.
//! - `ChFiDS_State stat;` (`cxx:1650`) is only assigned in the two branches of
//!   `sens(ic)`, so no default is needed.

use std::f64::consts::PI;

use crate::chfi3d::builder_helpers::angle_edge;
use crate::chfi_ds::State;

use super::surf_index::indices;
use super::MoreThreeCorner;

impl MoreThreeCorner {
    /// `cxx:1583-1676`: detects the "two fillets on two tangent living edges"
    /// configuration and, for three edges with two fillets, the
    /// `OnSame + OnDiff` one.
    ///
    /// `PerformOneCorner(Jndex, true)` is **not** issued here — the caller does
    /// it when [`Self::is_on_same_diff`] is set, exactly as `cxx:1671-1675`.
    pub(super) fn two_conges(&mut self, nconges: usize) -> Result<(), String> {
        if self.nedge < 1 {
            return Ok(());
        }
        let nedge = self.nedge as usize;
        let vertex = crate::shape::Vertex::wrap(self.vtx.clone())
            .ok_or_else(|| "PerformMoreThreeCorner: V1 is not a vertex".to_string())?;

        // `bool deuxconges, deuxcgnontg; deuxconges = false; trouve = false;`
        // (`cxx:1586-1588`).
        let mut deuxconges = false;

        // `cxx:1590-1605`: two fillets on two non-sharp edges that are tangent.
        if nconges == 2 {
            for ic in 0..nedge {
                let (icplus, _icmoins) = indices(nedge, ic);
                if !self.sharp[ic] && !self.sharp[icplus] {
                    let e1 = self.edge_at(ic)?;
                    let e2 = self.edge_at(icplus)?;
                    deuxconges = angle_edge(&vertex, &e1, &e2)?.abs() < 0.01;
                    if deuxconges {
                        break;
                    }
                }
            }
        }

        // `cxx:1608-1626`: same thing for two sharp edges, when none of them is
        // a free border (cts60072).
        if nconges == 2 && nedge == 4 {
            for ic in 0..nedge {
                if deuxconges {
                    break;
                }
                let (icplus, _icmoins) = indices(nedge, ic);
                if self.sharp[ic] && self.sharp[icplus] {
                    let e1 = self.edge_at(ic)?;
                    let e2 = self.edge_at(icplus)?;
                    if !self.is_free_border(&e1) && !self.is_free_border(&e2) {
                        let ang = angle_edge(&vertex, &e1, &e2)?.abs();
                        deuxconges = ang < 0.01 || (ang - PI).abs() < 0.01;
                    }
                }
            }
        }

        // `deuxcgnontg = nconges == 2 && nedge == 3 && !deuxconges;`
        // (`cxx:1628`).
        let deuxcgnontg = nconges == 2 && nedge == 3 && !deuxconges;

        // `if (deuxconges) for (ic = 0; ic < nedge; ic++) regul.SetValue(ic,
        // false);` (`cxx:1630-1636`).
        if deuxconges {
            for regul in self.regul.iter_mut().take(nedge) {
                *regul = false;
            }
        }

        // `cxx:1640-1670`: 3 edges & 2 fillets, one `OnSame` and one `OnDiff`.
        let mut is_on_same_diff = false;
        if deuxcgnontg {
            let mut is_on_same = false;
            let mut is_on_diff = false;
            for ic in 0..nedge {
                if self.sharp[ic] {
                    continue;
                }
                let stripe = self.cd[ic]
                    .clone()
                    .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
                let first = self.sens[ic] == 1;
                let stat = {
                    let guard = stripe
                        .read()
                        .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                    let spine = guard
                        .spine()
                        .cloned()
                        .ok_or_else(|| "PerformMoreThreeCorner: the stripe has no spine".to_string())?;
                    drop(guard);
                    let guard = spine
                        .read()
                        .map_err(|_| "PerformMoreThreeCorner: poisoned spine".to_string())?;
                    if first {
                        guard.first_status()
                    } else {
                        guard.last_status()
                    }
                };
                if stat == State::OnSame {
                    is_on_same = true;
                } else if stat == State::OnDiff {
                    is_on_diff = true;
                }
            }
            is_on_same_diff = is_on_same && is_on_diff;
        }

        self.deuxconges = deuxconges;
        self.deuxcgnontg = deuxcgnontg;
        self.is_on_same_diff = is_on_same_diff;
        Ok(())
    }
}
