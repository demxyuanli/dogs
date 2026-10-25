//! `BRepBlend_RstRstLineBuilder::InternalPerform`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstLineBuilder.cxx:571-1136`
//! (OCCT 8.0.0).
//!
//! The marching loop over the guide parameter. Each iteration solves the
//! inverse problem at `param`, checks both restrictions with `CheckInside`, and
//! reframes the section on the neighbouring data when either restriction falls
//! outside its range (`situonc1/2 != TopAbs_IN`) or when the contact with a
//! surface is lost (`decroch`):
//!
//! | status | reframes |
//! |---|---|
//! | `situonc1 != IN` | `Recadre1(FinvP1)` -- point/surface on `rst1` |
//! | `situonc2 != IN` | `Recadre2(FinvP2)` -- point/surface on `rst2` |
//! | `DecrochRst1` / `Both` | `Recadre1(Func, Finv1)` -- rst/surf on `rst1` |
//! | `DecrochRst2` / `Both` | `Recadre2(Func, Finv2)` -- rst/surf on `rst2` |
//!
//! The reduction of the four reframing results into `(trst11, trst12,
//! trst21, trst22)` is the same block as in `PerformFirstSection`
//! (`:743-810` here, `:368-435` there) with the extra `decroch` writes that
//! only this function performs (`:813-860`).
//!
//! Translation notes:
//! - `decroch` is *not* reset per iteration at the top of the loop; OCCT clears
//!   it after the reframing block (`:737`) and sets it again in the
//!   `recadrst1` / `recadrst2` arms, so it is the value that reaches
//!   `decrochfin` / `decrochdeb` (`:1126-1135`).
//! - `BRepBlend_BlendTool::Parameter(Vtxp1, rst1)` (`:754`, `:769`) and its
//!   `rst2` sibling (`:786`, `:803`) are called by OCCT without a null check;
//!   the port writes the parameter only when `Recadre` returned a vertex.
//! - `Blend_OK` and `Blend_StepTooSmall` at `param == Bound` set `Extrst1`
//!   directly but build `Extrst2` through `MakeExtremity` (`:948-956`,
//!   `:1017-1025`); `Blend_OnRst1` / `OnRst2` / `OnRst12` build both through
//!   `MakeExtremity` (`:1053-1056`, `:1073-1076`, `:1094-1097`) and are
//!   identical to each other, as in OCCT.

use occt_math::MathVector;

use crate::blend::{
    BlendCurvPointFuncInv, BlendDecrochStatus, BlendRstRstFunction, BlendStatus,
    BlendSurfCurvFuncInv,
};
use crate::brep_blend::{BlendTool, Extremity};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::fclass2d::FaceState;

use super::RstRstLineBuilder;

impl RstRstLineBuilder {
    /// `InternalPerform(Func, Finv1, FinvP1, Finv2, FinvP2, Bound)`
    /// (`cxx:571-1136`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn internal_perform(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        bound: f64,
    ) {
        let mut stepw = self.pasmax;
        let nbp = self.line_nb_points();
        if nbp >= 2 {
            // The last step is redone when it is not too small.
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
        let mut trst11 = 0.0;
        let mut trst12 = 0.0;
        let mut trst21 = 0.0;
        let mut trst22 = 0.0;
        let mut situonc1 = FaceState::Unknown;
        let mut situonc2 = FaceState::Unknown;
        let mut decroch = BlendDecrochStatus::NoDecroch;
        let mut arrive = false;

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
        let mut extrst1 = Extremity::new();
        let mut extrst2 = Extremity::new();

        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        let mut parinit = self.sol.clone();

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
                    if !self.check_inside(func, &mut situonc1, &mut situonc2, &mut decroch)
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
                let mut wp1 = bound;
                let mut wp2 = bound;
                let mut wrst1 = bound;
                let mut wrst2 = bound;
                let mut recadp1 = false;
                let mut recadp2 = false;
                let mut recadrst1 = false;
                let mut recadrst2 = false;
                let mut echecrecad = false;

                if situonc1 != FaceState::In {
                    // Inverse problem rst/rst on the first restriction.
                    recadp1 = self.recadre1_curv_point(finv_p1, &mut solinvp1, &mut is_vtxp1, &mut vtxp1);
                    if recadp1 {
                        wp1 = solinvp1.value(1);
                    } else {
                        echecrecad = true;
                    }
                }

                if situonc2 != FaceState::In {
                    // Inverse problem point/surface on the second restriction.
                    recadp2 = self.recadre2_curv_point(finv_p2, &mut solinvp2, &mut is_vtxp2, &mut vtxp2);
                    if recadp2 {
                        wp2 = solinvp2.value(1);
                    } else {
                        echecrecad = true;
                    }
                }

                if decroch == BlendDecrochStatus::DecrochRst1
                    || decroch == BlendDecrochStatus::DecrochBoth
                {
                    // Inverse problem rst1/surf1.
                    recadrst1 = self.recadre1_surf_curv(
                        func,
                        finv1,
                        &mut solinvrst1,
                        &mut is_vtxrst1,
                        &mut vtxrst1,
                    );
                    if recadrst1 {
                        wrst1 = solinvrst1.value(1);
                    } else {
                        echecrecad = true;
                    }
                }

                if decroch == BlendDecrochStatus::DecrochRst2
                    || decroch == BlendDecrochStatus::DecrochBoth
                {
                    // Inverse problem rst2/surf2.
                    recadrst2 = self.recadre2_surf_curv(
                        func,
                        finv2,
                        &mut solinvrst2,
                        &mut is_vtxrst2,
                        &mut vtxrst2,
                    );
                    if recadrst2 {
                        wrst2 = solinvrst2.value(1);
                    } else {
                        echecrecad = true;
                    }
                }

                decroch = BlendDecrochStatus::NoDecroch;
                if recadp1 || recadp2 || recadrst1 || recadrst2 {
                    echecrecad = false;
                }

                if !echecrecad {
                    // It is checked if the contact was lost or domain 1 was
                    // left (`cxx:743-777`).
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

                    // It is checked if the contact was lost or domain 2 was
                    // left (`cxx:779-810`).
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

                    // It is checked on which curve the contact is lost
                    // earlier (`cxx:812-862`).
                    if recadrst1 && recadrst2 {
                        if (wrst1 - wrst2).abs() < self.tolgui {
                            state = BlendStatus::OnRst12;
                            decroch = BlendDecrochStatus::DecrochBoth;
                            self.param = 0.5 * (wrst1 + wrst2);
                            self.sol.set_value(1, trst11);
                            self.sol.set_value(2, trst22);
                        } else if self.sens * (wrst1 - wrst2) < 0.0 {
                            // The contact is lost on Rst1.
                            state = BlendStatus::OnRst1;
                            decroch = BlendDecrochStatus::DecrochRst1;
                            self.param = wrst1;
                            self.sol.set_value(1, trst11);
                            self.sol.set_value(2, trst12);
                        } else {
                            // The contact is lost on rst2.
                            state = BlendStatus::OnRst2;
                            decroch = BlendDecrochStatus::DecrochRst2;
                            self.param = wrst2;
                            self.sol.set_value(1, trst21);
                            self.sol.set_value(2, trst22);
                        }
                        func.set(self.param);
                    } else if recadrst1 {
                        // Ground on rst1.
                        state = BlendStatus::OnRst1;
                        decroch = BlendDecrochStatus::DecrochRst1;
                        self.param = wrst1;
                        self.sol.set_value(1, trst11);
                        self.sol.set_value(2, trst12);
                        func.set(self.param);
                    } else if recadrst2 {
                        // Ground on rst2.
                        state = BlendStatus::OnRst2;
                        decroch = BlendDecrochStatus::DecrochRst2;
                        self.param = wrst2;
                        self.sol.set_value(1, trst21);
                        self.sol.set_value(2, trst22);
                        func.set(self.param);
                    }
                    // It is checked on which curve the contact is lost
                    // earlier (`cxx:863-908`).
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
                        // Ground on rst1.
                        state = BlendStatus::OnRst1;
                        self.param = wrst1;
                        self.sol.set_value(1, trst11);
                        self.sol.set_value(2, trst12);
                        func.set(self.param);
                    } else if recadp2 {
                        // Ground on rst2.
                        state = BlendStatus::OnRst2;
                        self.param = wrst2;
                        self.sol.set_value(1, trst21);
                        self.sol.set_value(2, trst22);
                        func.set(self.param);
                    } else {
                        state = BlendStatus::Ok;
                    }

                    state = self.test_arret(func, true, state);
                } else {
                    // Reframing failed: leave with points confoundus
                    // (`cxx:914-923`).
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
                        extrst1 = Extremity::on_curve(
                            *self.previous_p.point_on_c1(),
                            self.previous_p.parameter_on_c1().unwrap_or(0.0),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        self.make_extremity(
                            &mut extrst2,
                            false,
                            &self.rst2.clone(),
                            self.sol.value(2),
                            is_vtxrst2,
                            vtxrst2.as_ref(),
                        );
                        // Show that the end is on Bound.
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
                        extrst1 = Extremity::on_curve(
                            *self.previous_p.point_on_c1(),
                            self.previous_p.parameter_on_c1().unwrap_or(0.0),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        extrst2 = Extremity::on_curve(
                            *self.previous_p.point_on_c2(),
                            self.previous_p.parameter_on_c2().unwrap_or(0.0),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        arrive = true;
                    } else {
                        // There is no risk to exceed Bound.
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
                        extrst1 = Extremity::on_curve(
                            *self.previous_p.point_on_c1(),
                            self.previous_p.parameter_on_c1().unwrap_or(0.0),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        self.make_extremity(
                            &mut extrst2,
                            false,
                            &self.rst2.clone(),
                            self.sol.value(2),
                            is_vtxrst2,
                            vtxrst2.as_ref(),
                        );
                        // Indicate that the end is on Bound.
                    } else {
                        self.param += self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                }
                BlendStatus::OnRst1 => {
                    self.push_previous();
                    self.make_extremity(
                        &mut extrst1,
                        true,
                        &self.rst1.clone(),
                        self.sol.value(1),
                        is_vtxrst1,
                        vtxrst1.as_ref(),
                    );
                    self.make_extremity(
                        &mut extrst2,
                        false,
                        &self.rst2.clone(),
                        self.sol.value(2),
                        is_vtxrst2,
                        vtxrst2.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::OnRst2 => {
                    self.push_previous();
                    self.make_extremity(
                        &mut extrst1,
                        true,
                        &self.rst1.clone(),
                        self.sol.value(1),
                        is_vtxrst1,
                        vtxrst1.as_ref(),
                    );
                    self.make_extremity(
                        &mut extrst2,
                        false,
                        &self.rst2.clone(),
                        self.sol.value(2),
                        is_vtxrst2,
                        vtxrst2.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::OnRst12 => {
                    self.push_previous();
                    self.make_extremity(
                        &mut extrst1,
                        true,
                        &self.rst1.clone(),
                        self.sol.value(1),
                        is_vtxrst1,
                        vtxrst1.as_ref(),
                    );
                    self.make_extremity(
                        &mut extrst2,
                        false,
                        &self.rst2.clone(),
                        self.sol.value(2),
                        is_vtxrst2,
                        vtxrst2.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::SamePoints => {
                    // Stop.
                    extrst1 = Extremity::on_curve(
                        *self.previous_p.point_on_c1(),
                        self.previous_p.parameter_on_c1().unwrap_or(0.0),
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    extrst2 = Extremity::on_curve(
                        *self.previous_p.point_on_c2(),
                        self.previous_p.parameter_on_c2().unwrap_or(0.0),
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    arrive = true;
                }
                // OCCT's `default:` arm (`cxx:1118-1119`) does nothing and
                // leaves `Arrive` untouched. `Blend_Backward` is the only
                // remaining state and `TestArret` never returns it.
                BlendStatus::Backward => {}
            }

            if arrive {
                if self.sens > 0.0 {
                    if let Some(line) = self.line.as_mut() {
                        line.set_end_points(extrst1.clone(), extrst2.clone());
                    }
                    self.decrochfin = decroch;
                } else {
                    if let Some(line) = self.line.as_mut() {
                        line.set_start_points(extrst1.clone(), extrst2.clone());
                    }
                    self.decrochdeb = decroch;
                }
            }
        }
    }

    /// The `line->Append(previousP)` / `line->Prepend(previousP)` pair that
    /// every accepting arm performs (`cxx:934-943`, `:1002-1011`).
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
