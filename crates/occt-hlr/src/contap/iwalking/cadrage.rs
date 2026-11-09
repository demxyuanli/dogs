//! `IntWalk_IWalking::Cadrage` (`IntWalk_IWalking.gxx:405-582`): locate the
//! approached point in the natural domain of the surface, clip it on the
//! border and re-estimate the step.

use occt_core::precision::PCONFUSION;

use super::types::uv_on;
use super::ContapTheIWalking;

impl ContapTheIWalking {
    /// `Cadrage(BornInf, BornSup, UVap, Step, StepSign)` (`gxx:405-582`).
    ///
    /// Returns `true` when the point was put back on a border (the caller then
    /// restores the domain bounds and treats the step as a framed one).
    pub(super) fn cadrage(
        &self,
        born_inf: &mut [f64; 2],
        born_sup: &mut [f64; 2],
        uvap: &mut [f64; 2],
        step: &mut f64,
        step_sign: i32,
    ) -> bool {
        let duvx = self.previous_d2d.x();
        let duvy = self.previous_d2d.y();

        let (pu, pv) = uv_on(&self.previous_point, self.reversed);
        uvap[0] = pu;
        uvap[1] = pv;

        let step_sign = f64::from(step_sign);
        let u1 = uvap[0] + *step * duvx * step_sign;
        let v1 = uvap[1] + *step * duvy * step_sign;

        let infu = u1 <= born_inf[0] + PCONFUSION;
        let supu = u1 >= born_sup[0] - PCONFUSION;
        let infv = v1 <= born_inf[1] + PCONFUSION;
        let supv = v1 >= born_sup[1] - PCONFUSION;

        if !infu && !supu && !infv && !supv {
            uvap[0] = u1;
            uvap[1] = v1;
            return false;
        }

        if (infu || supu) && (infv || supv) {
            // `gxx:449-526`: both parameters leave the domain, clip on the
            // border reached first.
            let the_step_u = if infu {
                if duvx != 0.0 {
                    ((born_inf[0] - uvap[0]) / duvx).abs()
                } else {
                    *step
                }
            } else if duvx != 0.0 {
                ((born_sup[0] - uvap[0]) / duvx).abs()
            } else {
                *step
            };

            let the_step_v = if infv {
                if duvy != 0.0 {
                    ((born_inf[1] - uvap[1]) / duvy).abs()
                } else {
                    *step
                }
            } else if duvy != 0.0 {
                ((born_sup[1] - uvap[1]) / duvy).abs()
            } else {
                *step
            };

            if the_step_u <= the_step_v {
                *step = the_step_u;
                if infu {
                    uvap[0] = born_inf[0];
                    born_sup[0] = born_inf[0];
                } else {
                    uvap[0] = born_sup[0];
                    born_inf[0] = born_sup[0];
                }
                uvap[1] += *step * duvy * step_sign;
            } else {
                *step = the_step_v;
                if infv {
                    uvap[1] = born_inf[1];
                    born_sup[1] = born_inf[1];
                } else {
                    uvap[1] = born_sup[1];
                    born_inf[1] = born_sup[1];
                }
                uvap[0] += *step * duvx * step_sign;
            }
            return true;
        }

        if infu {
            // `gxx:529-541`.
            if duvx != 0.0 {
                let a_step = ((born_inf[0] - uvap[0]) / duvx).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_sup[0] = born_inf[0];
            uvap[0] = born_inf[0];
            uvap[1] += *step * duvy * step_sign;
            return true;
        }
        if supu {
            // `gxx:542-554`.
            if duvx != 0.0 {
                let a_step = ((born_sup[0] - uvap[0]) / duvx).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_inf[0] = born_sup[0];
            uvap[0] = born_sup[0];
            uvap[1] += *step * duvy * step_sign;
            return true;
        }
        if infv {
            // `gxx:556-567`.
            if duvy != 0.0 {
                let a_step = ((born_inf[1] - uvap[1]) / duvy).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_sup[1] = born_inf[1];
            uvap[0] += *step * duvx * step_sign;
            uvap[1] = born_inf[1];
            return true;
        }
        if supv {
            // `gxx:569-580`.
            if duvy != 0.0 {
                let a_step = ((born_sup[1] - uvap[1]) / duvy).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_inf[1] = born_sup[1];
            uvap[0] += *step * duvx * step_sign;
            uvap[1] = born_sup[1];
            return true;
        }
        true
    }
}
