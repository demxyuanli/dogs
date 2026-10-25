//! `BRepBlend_SurfRstLineBuilder::InternalPerform`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstLineBuilder.cxx:532-995`
//! (OCCT 8.0.0).
//!
//! The marching loop. Each iteration solves the inverse problem at `param`,
//! checks whether the solution is still inside both supports (`CheckInside`),
//! reframes the section on the neighbouring arc when it left `domain1` or the
//! range of `rst` (`Recadre`), and maps the resulting `Blend_Status` onto one of
//! seven outcomes:
//!
//! - `Blend_OK` / `Blend_StepTooSmall`: append `previousP`, advance;
//! - `Blend_StepTooLarge`: halve the step, or stop when it falls under
//!   `tolgui`;
//! - `Blend_OnRst1` / `Blend_OnRst2` / `Blend_OnRst12`: stop on the surface,
//!   on the restriction, or on both;
//! - `Blend_SamePoints`: stop, taking both extremities from `previousP`.
//!
//! Translation notes:
//! - `decroch` is deliberately not reset between iterations: OCCT only clears
//!   it inside the `bonpoint` block (`cxx:685`), so a failed solve keeps the
//!   value `CheckInside` produced and that value ends up in `decrochfin` /
//!   `decrochdeb` (`cxx:986-991`).
//! - `State` is declared outside the loop (`cxx:557`) but assigned on every
//!   path that reaches the switch, so the port keeps it as a loop-carried
//!   variable for fidelity.
//! - `BRepBlend_BlendTool::Parameter(Vtxp, rst)` (`cxx:703`, `:717`) would
//!   raise on the null handle OCCT can carry there; the port keeps `trst` at
//!   its previous value in that case.

use occt_core::gp::GpVec;

use crate::blend::{
    BlendFuncInv, BlendStatus, BlendSurfCurvFuncInv, BlendSurfPointFuncInv, BlendSurfRstFunction,
};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::{BlendTool, Extremity};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::fclass2d::FaceState;
use occt_math::MathVector;

use super::SurfRstLineBuilder;

impl SurfRstLineBuilder {
    /// `InternalPerform(Func, Finv, FinvP, FinvC, Bound)` (`cxx:532-995`).
    pub(super) fn internal_perform(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        bound: f64,
    ) {
        let mut stepw = self.pasmax;
        let nbp = self.line_nb_points();
        if nbp >= 2 {
            // The last step is reproduced when it is not too small.
            if self.sens < 0.0 {
                if let (Some(p2), Some(p1)) = (self.line_point(2), self.line_point(1)) {
                    stepw = p2.parameter() - p1.parameter();
                }
            } else if let (Some(pn), Some(pn1)) =
                (self.line_point(nbp), self.line_point(nbp - 1))
            {
                stepw = pn.parameter() - pn1.parameter();
            }
            stepw = stepw.max(100.0 * self.tolgui);
        }

        let mut parprec = self.param;
        if self.sens * (parprec - bound) >= -self.tolgui {
            return;
        }

        let mut state = BlendStatus::OnRst12;
        let mut situonc = FaceState::Unknown;
        let mut situons = FaceState::Unknown;
        let mut decroch = false;
        let mut arrive = false;

        let mut infbound = MathVector::with_init(1, 3, 0.0);
        let mut supbound = MathVector::with_init(1, 3, 0.0);
        let mut parinit = MathVector::with_init(1, 3, 0.0);
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
        let mut extrst = Extremity::new();
        let mut exts = Extremity::new();

        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        parinit = self.sol.clone();

        self.param = parprec + self.sens * stepw;
        if self.sens * (self.param - bound) > 0.0 {
            stepw = self.sens * (bound - parprec) * 0.5;
            self.param = parprec + self.sens * stepw;
        }

        while !arrive {
            let mut bonpoint = true;

            func.set(self.param);
            match super::super::walking::solve_root_set(
                func,
                &parinit,
                &tolerance,
                &infbound,
                &supbound,
                30,
            ) {
                Some(root) => {
                    self.sol = root;
                    if !self.check_inside(func, &mut situonc, &mut situons, &mut decroch)
                        && self.line_nb_points() == 1
                    {
                        state = BlendStatus::StepTooLarge;
                        bonpoint = false;
                    }
                }
                None => {
                    state = BlendStatus::StepTooLarge;
                    bonpoint = false;
                }
            }

            if bonpoint {
                let mut wp = bound;
                let mut wrst = bound;
                let mut ws = bound;
                let mut recadp = false;
                let mut recadrst = false;
                let mut recads = false;
                let mut echecrecad = false;
                let mut u = 0.0;
                let mut v = 0.0;
                let mut trst = 0.0;

                if situons == FaceState::Out || situons == FaceState::On {
                    // Inverse problem on the restriction of the second
                    // surface.
                    recads =
                        self.recadre_on_curv(finv_c, &mut solinvs, &mut arc, &mut is_vtxs, &mut vtxs);
                    if recads {
                        ws = solinvs.value(1);
                        // The deviation has to be re-evaluated (BUC60360).
                        let mut t = GpVec::zero();
                        let mut n = GpVec::zero();
                        func.set(ws);
                        if let Some(p2d) = arc.as_ref().and_then(|a| a.d0(solinvs.value(3))) {
                            u = p2d.x();
                            v = p2d.y();
                        }
                        self.sol.set_value(1, u);
                        self.sol.set_value(2, v);
                        self.sol.set_value(3, solinvs.value(2));
                        decroch = func.decroch(&self.sol, &mut t, &mut n);
                    } else {
                        echecrecad = true;
                    }
                }
                if situonc == FaceState::Out || situonc == FaceState::On {
                    // Inverse problem point/surface.
                    recadp = self.recadre_on_point(finv_p, &mut solinvp, &mut is_vtxp, &mut vtxp);
                    if recadp {
                        wp = solinvp.value(1);
                    } else {
                        echecrecad = true;
                    }
                }
                if decroch {
                    // Inverse problem restriction/surface.
                    recadrst = self.recadre_surf_rst(
                        func,
                        finv,
                        &mut solinvrst,
                        &mut is_vtxrst,
                        &mut vtxrst,
                    );
                    if recadrst {
                        wrst = solinvrst.value(2);
                    } else {
                        echecrecad = true;
                    }
                }
                decroch = false;
                if recadp || recads || recadrst {
                    echecrecad = false;
                }

                if !echecrecad {
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
                            decroch = true;
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
                        decroch = true;
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
                            self.sol.set_value(3, solinvs.value(3));
                        } else if self.sens * (ws - wrst) < 0.0 {
                            // Ground on the surface.
                            decroch = false;
                            state = BlendStatus::OnRst1;
                            self.param = ws;
                            if let Some(p2d) =
                                arc.as_ref().and_then(|a| a.d0(solinvs.value(3)))
                            {
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
                    } else {
                        state = BlendStatus::Ok;
                    }
                    state = self.test_arret(func, true, state);
                } else {
                    // Failed reframing: leave with PointsConfondus
                    // (`cxx:789-794`).
                    state = BlendStatus::SamePoints;
                }
            }

            match state {
                BlendStatus::Ok => {
                    self.push_previous();
                    parinit = self.sol.clone();
                    parprec = self.param;

                    if self.param == bound {
                        arrive = true;
                        exts = Extremity::on_surface(
                            *self.previous_p.point_on_s(),
                            self.sol.value(1),
                            self.sol.value(2),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        self.make_extremity(
                            &mut extrst,
                            false,
                            &self.rst.clone(),
                            self.sol.value(3),
                            is_vtxrst,
                            vtxrst.as_ref(),
                        );
                        // Indicate end on Bound.
                    } else {
                        self.param += self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                }
                BlendStatus::StepTooLarge => {
                    stepw /= 2.0;
                    if stepw.abs() < self.tolgui {
                        let (u, v) = self.previous_p.parameters_on_s().unwrap_or((0.0, 0.0));
                        exts = Extremity::on_surface(
                            *self.previous_p.point_on_s(),
                            u,
                            v,
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        extrst = Extremity::on_curve(
                            *self.previous_p.point_on_c(),
                            self.previous_p.parameter_on_c().unwrap_or(0.0),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        arrive = true;
                        // `cxx:848-856`: when the line already has two points
                        // the stop happens during the processing.
                    } else {
                        // No risk to exceed `Bound`.
                        self.param = parprec + self.sens * stepw;
                    }
                }
                BlendStatus::StepTooSmall => {
                    self.push_previous();
                    parinit = self.sol.clone();
                    parprec = self.param;

                    stepw = (1.5 * stepw).min(self.pasmax);
                    if self.param == bound {
                        arrive = true;
                        exts = Extremity::on_surface(
                            *self.previous_p.point_on_s(),
                            self.sol.value(1),
                            self.sol.value(2),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        self.make_extremity(
                            &mut extrst,
                            false,
                            &self.rst.clone(),
                            self.sol.value(3),
                            is_vtxrst,
                            vtxrst.as_ref(),
                        );
                        // Indicate end on Bound.
                    } else {
                        self.param += self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                }
                BlendStatus::OnRst1 => {
                    self.push_previous();
                    if let Some(a) = arc.clone() {
                        self.make_extremity(
                            &mut exts,
                            true,
                            &a,
                            solinvs.value(3),
                            is_vtxs,
                            vtxs.as_ref(),
                        );
                    }
                    self.make_extremity(
                        &mut extrst,
                        false,
                        &self.rst.clone(),
                        self.sol.value(3),
                        is_vtxrst,
                        vtxrst.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::OnRst2 => {
                    self.push_previous();
                    exts = Extremity::on_surface(
                        *self.previous_p.point_on_s(),
                        self.sol.value(1),
                        self.sol.value(2),
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    self.make_extremity(
                        &mut extrst,
                        false,
                        &self.rst.clone(),
                        self.sol.value(3),
                        is_vtxrst,
                        vtxrst.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::OnRst12 => {
                    self.push_previous();
                    if let Some(a) = arc.clone() {
                        self.make_extremity(
                            &mut exts,
                            true,
                            &a,
                            solinvs.value(1),
                            is_vtxs,
                            vtxs.as_ref(),
                        );
                    }
                    self.make_extremity(
                        &mut extrst,
                        false,
                        &self.rst.clone(),
                        self.sol.value(3),
                        is_vtxrst,
                        vtxrst.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::SamePoints => {
                    let (u, v) = self.previous_p.parameters_on_s().unwrap_or((0.0, 0.0));
                    exts = Extremity::on_surface(
                        *self.previous_p.point_on_s(),
                        u,
                        v,
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    extrst = Extremity::on_curve(
                        *self.previous_p.point_on_c(),
                        self.previous_p.parameter_on_c().unwrap_or(0.0),
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    arrive = true;
                }
                // OCCT's `default:` arm (`cxx:977-978`) does nothing and
                // leaves `Arrive` untouched. `Blend_Backward` is the only
                // remaining state and `TestArret` never returns it, so this
                // arm is unreachable.
                BlendStatus::Backward => {}
            }

            if arrive {
                if self.sens > 0.0 {
                    if let Some(line) = self.line.as_mut() {
                        line.set_end_points(exts.clone(), extrst.clone());
                    }
                    self.decrochfin = decroch;
                } else {
                    if let Some(line) = self.line.as_mut() {
                        line.set_start_points(exts.clone(), extrst.clone());
                    }
                    self.decrochdeb = decroch;
                }
            }
        }
    }

    /// The `line->Append(previousP)` / `line->Prepend(previousP)` pair that
    /// every accepting arm performs (`cxx:808-815`).
    fn push_previous(&mut self) {
        let p = self.previous_p;
        if let Some(line) = self.line.as_mut() {
            if self.sens > 0.0 {
                line.append(p);
            } else {
                line.prepend(p);
            }
        }
    }
}
