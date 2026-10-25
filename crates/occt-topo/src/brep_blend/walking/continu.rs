//! `BRepBlend_Walking::Continu` (both overloads) and `Complete`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Walking.cxx:624-755`
//! (OCCT 8.0.0).
//!
//! Each of the three entry points re-seeds `previousP` / `param` / `sol` from a
//! point of the line already walked and then runs `InternalPerform` again.
//! `Continu` (`cxx:624-649`) extends towards `P`; the `OnS1` / `OnS2` overload
//! (`cxx:651-724`) switches the classification off on one side and drops the
//! points walked while it was off if the other side did not reach the
//! restriction; `Complete` (`cxx:726-755`) walks the other way, once.
//!
//! Translation notes:
//! - `throw StdFail_NotDone()` becomes `Err`; OCCT's `return false` stays
//!   `Ok(false)`.

use crate::blend::{BlendFunction, BlendFuncInv};
use crate::brep_blend::Extremity;

use super::Walking;

impl Walking {
    /// `Continu(Func, FuncInv, P)` (`cxx:624-649`): walk further from the end of
    /// the line the parameter `P` lies beyond, or from its start when `P` lies
    /// before it.
    pub fn continu(
        &mut self,
        func: &mut dyn BlendFunction,
        func_inv: &mut dyn BlendFuncInv,
        p: f64,
    ) -> Result<bool, String> {
        self.require_done()?;

        {
            let line = self
                .line
                .as_ref()
                .ok_or_else(|| "BRepBlend_Walking::Continu: StdFail_NotDone".to_string())?;
            let first_bp = *line
                .point(1)
                .ok_or_else(|| "BRepBlend_Walking::Continu: StdFail_NotDone".to_string())?;
            let last_bp = line
                .point(line.nb_points())
                .copied()
                .ok_or_else(|| "BRepBlend_Walking::Continu: StdFail_NotDone".to_string())?;

            if p < first_bp.parameter() {
                self.sens = -1.0;
                self.previous_p = first_bp;
            } else if p > last_bp.parameter() {
                self.sens = 1.0;
                self.previous_p = last_bp;
            }
        }

        self.seed_sol_from_previous();
        self.internal_perform(func, func_inv, p);
        Ok(true)
    }

    /// `Continu(Func, FuncInv, P, OnS1)` (`cxx:651-724`): extend towards `P`
    /// with the classification switched off on one surface; the points walked
    /// are kept only when the other surface reaches a restriction.
    pub fn continu_on(
        &mut self,
        func: &mut dyn BlendFunction,
        func_inv: &mut dyn BlendFuncInv,
        p: f64,
        on_s1: bool,
    ) -> Result<bool, String> {
        self.require_done()?;

        let (ext1, ext2) = {
            let line = self
                .line
                .as_ref()
                .ok_or_else(|| "BRepBlend_Walking::Continu: StdFail_NotDone".to_string())?;
            if self.sens < 0.0 {
                let ext1: Extremity = line.start_point_on_first().clone();
                let ext2: Extremity = line.start_point_on_second().clone();
                if (on_s1 && ext1.nb_point_on_rst() == 0) || (!on_s1 && ext2.nb_point_on_rst() == 0) {
                    return Ok(false);
                }
                let first_bp = line
                    .point(1)
                    .copied()
                    .ok_or_else(|| "BRepBlend_Walking::Continu: StdFail_NotDone".to_string())?;
                self.previous_p = first_bp;
                (ext1, ext2)
            } else {
                let ext1: Extremity = line.end_point_on_first().clone();
                let ext2: Extremity = line.end_point_on_second().clone();
                if (on_s1 && ext1.nb_point_on_rst() == 0) || (!on_s1 && ext2.nb_point_on_rst() == 0) {
                    return Ok(false);
                }
                let last_bp = line
                    .point(line.nb_points())
                    .copied()
                    .ok_or_else(|| "BRepBlend_Walking::Continu: StdFail_NotDone".to_string())?;
                self.previous_p = last_bp;
                (ext1, ext2)
            }
        };

        let length = self.line.as_ref().map(|l| l.nb_points()).unwrap_or(0);
        self.seed_sol_from_previous();

        if on_s1 {
            self.clasons1 = false;
        } else {
            self.clasons2 = false;
        }

        self.internal_perform(func, func_inv, p);

        self.clasons1 = true;
        self.clasons2 = true;

        let newlength = self.line.as_ref().map(|l| l.nb_points()).unwrap_or(0);
        let Some(line) = self.line.as_mut() else {
            return Ok(true);
        };
        if self.sens < 0.0 {
            if (on_s1 && line.start_point_on_second().nb_point_on_rst() == 0)
                || (!on_s1 && line.start_point_on_first().nb_point_on_rst() == 0)
            {
                line.remove(1, newlength - length);
                line.set_start_points(ext1, ext2);
                return Ok(false);
            }
        } else if (on_s1 && line.end_point_on_second().nb_point_on_rst() == 0)
            || (!on_s1 && line.end_point_on_first().nb_point_on_rst() == 0)
        {
            line.remove(length, newlength);
            line.set_end_points(ext1, ext2);
            return Ok(false);
        }
        Ok(true)
    }

    /// `Complete(Func, FuncInv, Pmin)` (`cxx:726-755`): walk the other way,
    /// once (`iscomplete`).
    pub fn complete(
        &mut self,
        func: &mut dyn BlendFunction,
        func_inv: &mut dyn BlendFuncInv,
        pmin: f64,
    ) -> Result<bool, String> {
        self.require_done()?;
        if self.iscomplete {
            return Ok(true);
        }

        let picked = {
            let line = self
                .line
                .as_ref()
                .ok_or_else(|| "BRepBlend_Walking::Complete: StdFail_NotDone".to_string())?;
            if self.sens > 0.0 {
                line.point(1).copied()
            } else {
                line.point(line.nb_points()).copied()
            }
        };
        if let Some(bp) = picked {
            self.previous_p = bp;
        }

        self.sens = -self.sens;
        self.seed_sol_from_previous();
        self.internal_perform(func, func_inv, pmin);

        self.iscomplete = true;
        Ok(true)
    }

    /// `if (!done) throw StdFail_NotDone()` (`cxx:626`, `:657`, `:728`).
    fn require_done(&self) -> Result<(), String> {
        if !self.done {
            return Err("BRepBlend_Walking: StdFail_NotDone".to_string());
        }
        Ok(())
    }

    /// `param = previousP.Parameter(); previousP.ParametersOnS1(sol(1), sol(2));
    /// previousP.ParametersOnS2(sol(3), sol(4));` (`cxx:645-648`, `:684-687`,
    /// `:746-749`).
    fn seed_sol_from_previous(&mut self) {
        self.param = self.previous_p.parameter();
        if let Ok((u, v)) = self.previous_p.parameters_on_s1() {
            self.sol.set_value(1, u);
            self.sol.set_value(2, v);
        }
        if let Ok((u, v)) = self.previous_p.parameters_on_s2() {
            self.sol.set_value(3, u);
            self.sol.set_value(4, v);
        }
    }
}
