//! `BRepBlend_RstRstLineBuilder::PerformFirstSection` and `Complete`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstLineBuilder.cxx:284-569`
//! (OCCT 8.0.0).
//!
//! `PerformFirstSection` (`cxx:284-533`) solves the initial inverse problem and
//! then reframes it on up to four supports at once -- `rst1`, `surf1` as a
//! point, `rst2`, `surf2` as a point -- selected by `RecRst1` / `RecP1` /
//! `RecRst2` / `RecP2`. The reframing data is first reduced per side
//! (`:368-435`), then the earliest exit decides the stop state (`:437-529`).
//!
//! Translation notes:
//! - The `recadp1 && recadrst1` branch writes `IsVtxrst2` / `Vtxrst2` from the
//!   `p1` data (`cxx:379-381`); that is what OCCT 8.0.0 does, and the sibling
//!   `recadp1` branch right below writes `IsVtxrst1` / `Vtxrst1` from the same
//!   source (`:394-396`). Both are kept verbatim.
//! - The `OnRst12` arms assign `sol(1) = trst11` and `sol(2) = trst22`, not
//!   `trst12` / `trst21` (`:441-444`, `:488-491`); that asymmetry is OCCT's.
//! - `pasmax` and `fleche` are not set here (`cxx:302-310`): OCCT's initial list
//!   leaves them at whatever `Perform` last wrote.

use occt_math::MathVector;

use crate::blend::{
    BlendCurvPointFuncInv, BlendRstRstFunction, BlendStatus, BlendSurfCurvFuncInv,
};
use crate::brep_blend::BlendTool;
use crate::brep_top_adaptor_hvertex::HVertex;

use super::RstRstLineBuilder;

impl RstRstLineBuilder {
    /// `PerformFirstSection(Func, Finv1, FinvP1, Finv2, FinvP2, Pdep, Pmax,
    /// Soldep, Tol3d, TolGuide, RecRst1, RecP1, RecRst2, RecP2, Psol, ParSol)`
    /// (`cxx:284-533`).
    ///
    /// The vertex flags OCCT accumulates here (`IsVtxrst1/2`, `Vtxrst1/2`,
    /// `:319`) are never read before returning: only `param` and `sol` leave
    /// through `Psol` / `ParSol`, so the writes are dead in OCCT as well.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn perform_first_section(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        pdep: f64,
        pmax: f64,
        sol_dep: &MathVector,
        tol3d: f64,
        tol_guide: f64,
        rec_rst1: bool,
        rec_p1: bool,
        rec_rst2: bool,
        rec_p2: bool,
        psol: &mut f64,
        par_sol: &mut MathVector,
    ) -> bool {
        self.done = false;
        self.iscomplete = false;
        self.comptra = false;
        self.line = Some(crate::brep_blend::Line::new());
        self.tolpoint3d = tol3d;
        self.tolgui = tol_guide.abs();
        self.rebrou = false;

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        let mut state = BlendStatus::OnRst12;
        let mut trst11 = 0.0;
        let mut trst12 = 0.0;
        let mut trst21 = 0.0;
        let mut trst22 = 0.0;
        let mut infbound = MathVector::with_init(1, 2, 0.0);
        let mut supbound = MathVector::with_init(1, 2, 0.0);
        let mut tolerance = MathVector::with_init(1, 2, 0.0);
        let mut solinvp1 = MathVector::new(1, 2);
        let mut solinvp2 = MathVector::new(1, 2);
        let mut solinvrst1 = MathVector::new(1, 3);
        let mut solinvrst2 = MathVector::new(1, 3);
        let mut vtxp1: Option<HVertex> = None;
        let mut vtxp2: Option<HVertex> = None;
        let mut vtxrst1: Option<HVertex> = None;
        let mut vtxrst2: Option<HVertex> = None;
        let mut is_vtxp1 = false;
        let mut is_vtxp2 = false;
        let mut is_vtxrst1 = false;
        let mut is_vtxrst2 = false;
        let mut wp1;
        let mut wp2;
        let mut wrst1;
        let mut wrst2;

        wp1 = pmax;
        wp2 = pmax;
        wrst1 = pmax;
        wrst2 = pmax;

        self.param = pdep;
        func.set(self.param);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        let Some(root) = super::super::walking::solve_root_set(
            func,
            sol_dep,
            &tolerance,
            &infbound,
            &supbound,
            30,
        ) else {
            return false;
        };
        self.sol = root;

        let mut recadrst1 = rec_rst1
            && self.recadre1_surf_curv(
                func,
                finv1,
                &mut solinvrst1,
                &mut is_vtxrst1,
                &mut vtxrst1,
            );
        if recadrst1 {
            wrst1 = solinvrst1.value(1);
        }

        let mut recadp1 = rec_p1
            && self.recadre1_curv_point(finv_p1, &mut solinvp1, &mut is_vtxp1, &mut vtxp1);
        if recadp1 {
            wp1 = solinvp1.value(1);
        }

        let mut recadrst2 = rec_rst2
            && self.recadre2_surf_curv(
                func,
                finv2,
                &mut solinvrst2,
                &mut is_vtxrst2,
                &mut vtxrst2,
            );
        if recadrst2 {
            wrst2 = solinvrst2.value(1);
        }

        let mut recadp2 = rec_p2
            && self.recadre2_curv_point(finv_p2, &mut solinvp2, &mut is_vtxp2, &mut vtxp2);
        if recadp2 {
            wp2 = solinvp2.value(1);
        }

        if !recadrst1 && !recadp1 && !recadrst2 && !recadp2 {
            return false;
        }

        // It is checked if the contact was lost or domain 1 was left
        // (`cxx:370-402`).
        //
        // `BRepBlend_BlendTool::Parameter(Vtxp1, rst1)` (`:378`, `:395`) is
        // called by OCCT whether or not a vertex was found; a null handle would
        // dereference there, so the port only writes the parameter when
        // `Recadre` actually returned one and leaves `trst11` untouched
        // otherwise.
        if recadp1 && recadrst1 {
            if self.sens * (wrst1 - wp1) > self.tolgui {
                // At first one leaves the domain.
                wrst1 = wp1;
                trst12 = solinvp1.value(2);
                if let Some(vtx) = vtxp1.as_ref() {
                    trst11 = BlendTool::parameter(vtx, &self.rst1);
                }
                is_vtxrst2 = is_vtxp1;
                vtxrst2 = vtxp1.clone();
                recadrst1 = false;
            } else {
                // The contact is lost.
                trst11 = solinvrst1.value(3);
                trst12 = solinvrst1.value(2);
                recadp1 = false;
            }
        } else if recadp1 {
            wrst1 = wp1;
            trst12 = solinvp1.value(2);
            if let Some(vtx) = vtxp1.as_ref() {
                trst11 = BlendTool::parameter(vtx, &self.rst1);
            }
            is_vtxrst1 = is_vtxp1;
            vtxrst1 = vtxp1.clone();
        } else if recadrst1 {
            trst11 = solinvrst1.value(3);
            trst12 = solinvrst1.value(2);
        }

        // It is checked if the contact was lost or domain 2 was left
        // (`cxx:405-435`).
        if recadp2 && recadrst2 {
            if self.sens * (wrst2 - wp2) > self.tolgui {
                wrst2 = wp2;
                trst21 = solinvp2.value(2);
                if let Some(vtx) = vtxp2.as_ref() {
                    trst22 = BlendTool::parameter(vtx, &self.rst2);
                }
                is_vtxrst2 = is_vtxp2;
                vtxrst2 = vtxp2.clone();
                recadrst2 = false;
            } else {
                trst22 = solinvrst2.value(3);
                trst21 = solinvrst2.value(2);
                recadp2 = false;
            }
        } else if recadp2 {
            wrst2 = wp2;
            trst21 = solinvp2.value(2);
            if let Some(vtx) = vtxp2.as_ref() {
                trst22 = BlendTool::parameter(vtx, &self.rst2);
            }
            is_vtxrst2 = is_vtxp2;
            vtxrst2 = vtxp2.clone();
        } else if recadrst2 {
            trst22 = solinvrst2.value(3);
            trst21 = solinvrst2.value(2);
        }

        // It is checked on which curve the contact is lost earlier
        // (`cxx:437-463`).
        if recadrst1 && recadrst2 {
            if (wrst1 - wrst2).abs() < self.tolgui {
                state = BlendStatus::OnRst12;
                self.param = 0.5 * (wrst1 + wrst2);
                self.sol.set_value(1, trst11);
                self.sol.set_value(2, trst22);
            } else if self.sens * (wrst1 - wrst2) < 0.0 {
                // Contact lost on Rst1.
                state = BlendStatus::OnRst1;
                self.param = wrst1;
                self.sol.set_value(1, trst11);
                self.sol.set_value(2, trst12);
            } else {
                // Contact lost on rst2.
                state = BlendStatus::OnRst2;
                self.param = wrst2;
                self.sol.set_value(1, trst21);
                self.sol.set_value(2, trst22);
            }
            func.set(self.param);
        } else if recadrst1 {
            state = BlendStatus::OnRst1;
            self.param = wrst1;
            self.sol.set_value(1, trst11);
            self.sol.set_value(2, trst12);
            func.set(self.param);
        } else if recadrst2 {
            state = BlendStatus::OnRst2;
            self.param = wrst2;
            self.sol.set_value(1, trst21);
            self.sol.set_value(2, trst22);
            func.set(self.param);
        }
        // It is checked on which curves one leaves first (`cxx:485-511`).
        else if recadp1 && recadp2 {
            if (wrst1 - wrst2).abs() < self.tolgui {
                state = BlendStatus::OnRst12;
                self.param = 0.5 * (wrst1 + wrst2);
                self.sol.set_value(1, trst11);
                self.sol.set_value(2, trst22);
            } else if self.sens * (wrst1 - wrst2) < 0.0 {
                state = BlendStatus::OnRst1;
                self.param = wrst1;
                self.sol.set_value(1, trst11);
                self.sol.set_value(2, trst12);
            } else {
                state = BlendStatus::OnRst2;
                self.param = wrst2;
                self.sol.set_value(1, trst21);
                self.sol.set_value(2, trst22);
            }
            func.set(self.param);
        } else if recadp1 {
            state = BlendStatus::OnRst1;
            self.param = wrst1;
            self.sol.set_value(1, trst11);
            self.sol.set_value(2, trst12);
            func.set(self.param);
        } else if recadp2 {
            state = BlendStatus::OnRst2;
            self.param = wrst2;
            self.sol.set_value(1, trst21);
            self.sol.set_value(2, trst22);
            func.set(self.param);
        }

        let _state = self.test_arret(func, false, state);
        *psol = self.param;
        *par_sol = self.sol.clone();
        true
    }

    /// `Complete(Func, Finv1, FinvP1, Finv2, FinvP2, Pmin)` (`cxx:537-569`).
    pub fn complete(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        pmin: f64,
    ) -> Result<bool, String> {
        if !self.done {
            return Err("BRepBlend_RstRstLineBuilder::Complete: StdFail_NotDone".to_string());
        }
        if self.iscomplete {
            return Ok(true);
        }

        let nb = self.line_nb_points();
        let point = if self.sens > 0.0 {
            self.line_point(1)
        } else {
            self.line_point(nb)
        };
        let Some(point) = point else {
            return Ok(false);
        };
        self.previous_p = point;

        self.sens = -self.sens;
        self.param = self.previous_p.parameter();
        self.sol
            .set_value(1, self.previous_p.parameter_on_c1().unwrap_or(0.0));
        self.sol
            .set_value(2, self.previous_p.parameter_on_c2().unwrap_or(0.0));

        self.internal_perform(func, finv1, finv_p1, finv2, finv_p2, pmin);
        self.iscomplete = true;
        Ok(true)
    }
}
