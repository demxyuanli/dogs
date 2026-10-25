//! `BRepBlend_SurfRstLineBuilder::PerformFirstSection` and `Complete`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstLineBuilder.cxx:326-529`
//! (OCCT 8.0.0).
//!
//! `PerformFirstSection` (`cxx:326-495`) solves the initial inverse problem and
//! then reframes it on up to three supports at once, the selection of which is
//! driven by `RecRst` / `RecP` / `RecS`. The earliest of the three reframings
//! wins (`cxx:443-465`), which is what tells the caller where the section
//! actually starts.
//!
//! `Complete` (`cxx:498-529`) walks the other way from the first point of the
//! line. Unlike [`super::cs_walking`] it does not rebuild `previousP` from
//! `firstsol`: it takes the extremity point of the existing line and reloads
//! `sol` from its parameters (`cxx:513-524`).

use crate::blend::{
    BlendFuncInv, BlendStatus, BlendSurfCurvFuncInv, BlendSurfPointFuncInv, BlendSurfRstFunction,
};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::BlendTool;
use crate::brep_top_adaptor_hvertex::HVertex;
use occt_math::MathVector;

use super::SurfRstLineBuilder;

impl SurfRstLineBuilder {
    /// `PerformFirstSection(F, Finv, FinvP, FinvC, Pdep, Pmax, ParDep, Tol3d,
    /// Tol2d, TolGuide, RecRst, RecP, RecS, Psol, ParSol)` (`cxx:326-495`).
    ///
    /// `pasmax` and `fleche` are deliberately left untouched: OCCT's initial
    /// list (`cxx:344-351`) does not set them either.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_first_section(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        pdep: f64,
        pmax: f64,
        par_dep: &MathVector,
        tol3d: f64,
        tol2d: f64,
        tol_guide: f64,
        rec_rst: bool,
        rec_p: bool,
        rec_s: bool,
        psol: &mut f64,
        par_sol: &mut MathVector,
    ) -> bool {
        self.done = false;
        self.iscomplete = false;
        self.comptra = false;
        self.line = Some(crate::brep_blend::Line::new());
        self.tolpoint3d = tol3d;
        self.tolpoint2d = tol2d;
        self.tolgui = tol_guide.abs();
        self.rebrou = false;

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        let mut state = BlendStatus::OnRst12;
        let mut trst = 0.0;
        let mut u = 0.0;
        let mut v = 0.0;
        let mut infbound = MathVector::with_init(1, 3, 0.0);
        let mut supbound = MathVector::with_init(1, 3, 0.0);
        let mut tolerance = MathVector::with_init(1, 3, 0.0);
        let mut solinvp = MathVector::new(1, 3);
        let mut solinvrst = MathVector::new(1, 4);
        let mut solinvs = MathVector::new(1, 3);
        let mut vtxp: Option<HVertex> = None;
        let mut vtxrst: Option<HVertex> = None;
        let mut vtxs: Option<HVertex> = None;
        let mut is_vtxp = false;
        let mut is_vtxrst = false;
        let mut is_vtxs = false;
        let mut arc: Option<Curve2dAdaptor> = None;

        let mut ws = pmax;
        let mut wrst = pmax;
        let mut wp = pmax;

        self.param = pdep;
        func.set(self.param);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        let Some(root) = super::super::walking::solve_root_set(
            func,
            par_dep,
            &tolerance,
            &infbound,
            &supbound,
            30,
        ) else {
            return false;
        };
        self.sol = root;

        let mut recads = rec_s
            && self.recadre_on_curv(finv_c, &mut solinvs, &mut arc, &mut is_vtxs, &mut vtxs);
        if recads {
            ws = solinvs.value(1);
        }
        let mut recadp = rec_p
            && self.recadre_on_point(finv_p, &mut solinvp, &mut is_vtxp, &mut vtxp);
        if recadp {
            wp = solinvp.value(1);
        }
        let mut recadrst = rec_rst
            && self.recadre_surf_rst(func, finv, &mut solinvrst, &mut is_vtxrst, &mut vtxrst);
        if recadrst {
            wrst = solinvrst.value(2);
        }

        if !recads && !recadp && !recadrst {
            return false;
        }

        if recadp && recadrst {
            if self.sens * (wrst - wp) > self.tolgui {
                // The first one leaves the domain.
                wrst = wp;
                u = solinvp.value(2);
                v = solinvp.value(3);
                if let Some(vtx) = vtxp.as_ref() {
                    trst = BlendTool::parameter(vtx, &self.rst);
                }
                is_vtxrst = is_vtxp;
                vtxrst = vtxp.clone();
            } else {
                u = solinvrst.value(3);
                v = solinvrst.value(4);
                trst = solinvrst.value(1);
            }
        } else if recadp {
            wrst = wp;
            u = solinvp.value(2);
            v = solinvp.value(3);
            if let Some(vtx) = vtxp.as_ref() {
                trst = BlendTool::parameter(vtx, &self.rst);
            }
            is_vtxrst = is_vtxp;
            vtxrst = vtxp.clone();
            recadrst = true;
        } else if recadrst {
            u = solinvrst.value(3);
            v = solinvrst.value(4);
            trst = solinvrst.value(1);
        }

        if recads && recadrst {
            if (ws - wrst).abs() < self.tolgui {
                state = BlendStatus::OnRst12;
                self.param = 0.5 * (ws + wrst);
                self.sol.set_value(1, u);
                self.sol.set_value(2, v);
                self.sol.set_value(3, solinvs.value(2));
            } else if self.sens * (ws - wrst) < 0.0 {
                // Ground on the surface.
                state = BlendStatus::OnRst1;
                self.param = ws;
                if let Some(p2d) = arc.as_ref().and_then(|a| a.d0(solinvs.value(3))) {
                    u = p2d.x();
                    v = p2d.y();
                }
                self.sol.set_value(1, u);
                self.sol.set_value(2, v);
                self.sol.set_value(3, solinvs.value(2));
            } else {
                // Ground on the restriction.
                state = BlendStatus::OnRst2;
                self.param = wrst;
                self.sol.set_value(1, u);
                self.sol.set_value(2, v);
                self.sol.set_value(3, trst);
            }
            func.set(self.param);
        } else if recads {
            // Ground on the surface.
            state = BlendStatus::OnRst1;
            self.param = ws;
            if let Some(p2d) = arc.as_ref().and_then(|a| a.d0(solinvs.value(3))) {
                u = p2d.x();
                v = p2d.y();
            }
            self.sol.set_value(1, u);
            self.sol.set_value(2, v);
            self.sol.set_value(3, solinvs.value(2));
            func.set(self.param);
        } else if recadrst {
            // Ground on the restriction.
            state = BlendStatus::OnRst2;
            self.param = wrst;
            self.sol.set_value(1, u);
            self.sol.set_value(2, v);
            self.sol.set_value(3, trst);
            func.set(self.param);
        }

        let _state = self.test_arret(func, false, state);
        *psol = self.param;
        *par_sol = self.sol.clone();
        true
    }

    /// `Complete(F, Finv, FinvP, FinvC, Pmin)` (`cxx:498-529`).
    pub fn complete(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        pmin: f64,
    ) -> Result<bool, String> {
        if !self.done {
            return Err("BRepBlend_SurfRstLineBuilder::Complete: StdFail_NotDone".to_string());
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
        let (u, v) = self.previous_p.parameters_on_s().unwrap_or((0.0, 0.0));
        self.sol.set_value(1, u);
        self.sol.set_value(2, v);
        self.sol.set_value(3, self.previous_p.parameter_on_c().unwrap_or(0.0));

        self.internal_perform(func, finv, finv_p, finv_c, pmin);
        self.iscomplete = true;
        Ok(true)
    }
}
