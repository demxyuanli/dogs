//! `BRepBlend_RstRstLineBuilder::TestArret`, `CheckDeflectionOnRst1`,
//! `CheckDeflectionOnRst2` and `CheckInside`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstLineBuilder.cxx:1588-1964`
//! (OCCT 8.0.0).
//!
//! Both deflection checks (`cxx:1588-1674` for `rst1`, `cxx:1678-1764` for
//! `rst2`) are the 3d rule of `Blend_CSWalking` applied to the restriction
//! side: `CosRef3D = 0.98`, the chord guard `0.01 * tolpoint3d`, and the
//! sagitta estimate `|Tg - Tgsurf|^2 * Norme / 64`. They differ only in whether
//! they read `PointOnC1` / `TangentOnC1` or `PointOnC2` / `TangentOnC2`.
//!
//! `TestArret` (`cxx:1775-1917`) builds the candidate section, runs both
//! deflection checks when `TestDeflection` asks for them, and computes the two
//! line transitions the first time it can (`comptra == false`, non-tangency
//! point): each transition compares the 2d tangent of the restriction against
//! the 2d tangent of the section on that restriction, taking the orientation of
//! the restriction in its domain into account (`cxx:1852-1888`).
//!
//! `CheckInside` (`cxx:1921-1964`) reports the position of `sol` against both
//! restriction ranges and the `Decroch` status of the function.
//!
//! Translation notes:
//! - The transition gate uses `tolpoint3d` for *both* restrictions
//!   (`cxx:1859`, `:1874`), where the surface/restriction builder uses a
//!   separate `1e-8` for the second one; that is OCCT's difference, kept here.
//! - The `rst2` half of the transition block lives inside the `rst1` half's
//!   successful test (`cxx:1870-1886`): a null 2d tangent on `rst1` leaves both
//!   transitions undecided. Kept.
//! - `Decroch` writes four vectors `(NRst1, TgRst1, NRst2, TgRst2)`
//!   (`Blend_RstRstFunction.hxx:154-157`); the caller at `cxx:1961` names them
//!   `tgrst1, norst1, tgrst2, norst2` and passes them positionally, which is
//!   what the port does.

use occt_core::gp::GpVec;

use crate::abs::Orientation;
use crate::blend::{BlendDecrochStatus, BlendPoint, BlendRstRstFunction, BlendStatus};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::fclass2d::FaceState;
use crate::intsurf_transition::IntSurfTypeTrans;

use super::RstRstLineBuilder;

/// `CosRef3D` of both deflection checks (`cxx:1593`, `cxx:1683`).
const COS_REF3D: f64 = 0.98;

/// `IntSurf_TypeTrans ConvOrToTra(TopAbs_Orientation)` (`cxx:1766-1772`).
fn conv_or_to_tra(o: Orientation) -> IntSurfTypeTrans {
    if o == Orientation::Forward {
        IntSurfTypeTrans::In
    } else {
        IntSurfTypeTrans::Out
    }
}

impl RstRstLineBuilder {
    /// `CheckDeflectionOnRst1(CurPoint)` (`cxx:1588-1674`).
    fn check_deflection_on_rst1(&self, cur_point: &BlendPoint) -> BlendStatus {
        let cur_tangent = cur_point.is_tangency_point();
        let prev_tangent = self.previous_p.is_tangency_point();

        let psurf = *cur_point.point_on_c1();
        let tgsurf = cur_point
            .tangent_on_c1()
            .ok()
            .copied()
            .unwrap_or_else(GpVec::zero);
        let prev_p = *self.previous_p.point_on_c1();
        let prev_tg = self
            .previous_p
            .tangent_on_c1()
            .ok()
            .copied()
            .unwrap_or_else(GpVec::zero);

        let corde = GpVec::from_pnts(&prev_p, &psurf);
        let norme = corde.square_magnitude();
        let prev_norme = if prev_tangent {
            0.0
        } else {
            prev_tg.square_magnitude()
        };

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // It can be necessary to force the same point.
            return BlendStatus::SamePoints;
        }
        if !prev_tangent {
            if prev_norme <= toler3d * toler3d {
                return BlendStatus::SamePoints;
            }
            let cosi = self.sens * corde.dot(&prev_tg);
            if cosi < 0.0 {
                // 3d angle greater than pi/2: going backwards.
                return BlendStatus::Backward;
            }
            let cosi2 = cosi * cosi / prev_norme / norme;
            if cosi2 < COS_REF3D {
                return BlendStatus::StepTooLarge;
            }
        }

        if !cur_tangent {
            let cosi = self.sens * corde.dot(&tgsurf);
            let cosi2 = cosi * cosi / tgsurf.square_magnitude() / norme;
            if cosi2 < COS_REF3D || cosi < 0.0 {
                return BlendStatus::StepTooLarge;
            }
        }

        if !cur_tangent && !prev_tangent {
            let diff = prev_tg.normalized().subtracted(&tgsurf.normalized());
            let fleche_courante = diff.square_magnitude() * norme / 64.0;
            if fleche_courante <= 0.25 * self.fleche * self.fleche {
                return BlendStatus::StepTooSmall;
            }
            if fleche_courante > self.fleche * self.fleche {
                return BlendStatus::StepTooLarge;
            }
        }
        BlendStatus::Ok
    }

    /// `CheckDeflectionOnRst2(CurPoint)` (`cxx:1678-1764`).
    fn check_deflection_on_rst2(&self, cur_point: &BlendPoint) -> BlendStatus {
        let cur_tangent = cur_point.is_tangency_point();
        let prev_tangent = self.previous_p.is_tangency_point();

        let psurf = *cur_point.point_on_c2();
        let tgsurf = cur_point
            .tangent_on_c2()
            .ok()
            .copied()
            .unwrap_or_else(GpVec::zero);
        let prev_p = *self.previous_p.point_on_c2();
        let prev_tg = self
            .previous_p
            .tangent_on_c2()
            .ok()
            .copied()
            .unwrap_or_else(GpVec::zero);

        let corde = GpVec::from_pnts(&prev_p, &psurf);
        let norme = corde.square_magnitude();
        let prev_norme = if prev_tangent {
            0.0
        } else {
            prev_tg.square_magnitude()
        };

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // It can be necessary to force the same point.
            return BlendStatus::SamePoints;
        }
        if !prev_tangent {
            if prev_norme <= toler3d * toler3d {
                return BlendStatus::SamePoints;
            }
            let cosi = self.sens * corde.dot(&prev_tg);
            if cosi < 0.0 {
                // 3d angle greater than pi/2: going backwards.
                return BlendStatus::Backward;
            }
            let cosi2 = cosi * cosi / prev_norme / norme;
            if cosi2 < COS_REF3D {
                return BlendStatus::StepTooLarge;
            }
        }

        if !cur_tangent {
            let cosi = self.sens * corde.dot(&tgsurf);
            let cosi2 = cosi * cosi / tgsurf.square_magnitude() / norme;
            if cosi2 < COS_REF3D || cosi < 0.0 {
                return BlendStatus::StepTooLarge;
            }
        }

        if !cur_tangent && !prev_tangent {
            let diff = prev_tg.normalized().subtracted(&tgsurf.normalized());
            let fleche_courante = diff.square_magnitude() * norme / 64.0;
            if fleche_courante <= 0.25 * self.fleche * self.fleche {
                return BlendStatus::StepTooSmall;
            }
            if fleche_courante > self.fleche * self.fleche {
                return BlendStatus::StepTooLarge;
            }
        }
        BlendStatus::Ok
    }

    /// `TestArret(Func, TestDeflection, State)` (`cxx:1775-1917`).
    pub(super) fn test_arret(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        test_deflection: bool,
        state: BlendStatus,
    ) -> BlendStatus {
        if !func.is_solution(&self.sol, self.tolpoint3d) {
            return BlendStatus::StepTooLarge;
        }

        let cur_tangent = func.is_tangency_point();
        let ptrst1 = *func.point_on_rst1();
        let ptrst2 = *func.point_on_rst2();
        let pt2drst1 = *func.pnt2d_on_rst1();
        let pt2drst2 = *func.pnt2d_on_rst2();

        let mut curpoint = BlendPoint::new();
        if cur_tangent {
            curpoint.set_on_two_pcurves(
                ptrst1,
                ptrst2,
                self.param,
                pt2drst1.x(),
                pt2drst1.y(),
                pt2drst2.x(),
                pt2drst2.y(),
                self.sol.value(1),
                self.sol.value(2),
            );
        } else {
            curpoint.set_on_two_pcurves_with_tangents(
                ptrst1,
                ptrst2,
                self.param,
                pt2drst1.x(),
                pt2drst1.y(),
                pt2drst2.x(),
                pt2drst2.y(),
                self.sol.value(1),
                self.sol.value(2),
                *func.tangent_on_rst1(),
                *func.tangent_on_rst2(),
                *func.tangent2d_on_rst1(),
                *func.tangent2d_on_rst2(),
            );
        }

        let (mut state_rst1, mut state_rst2) = if test_deflection {
            (
                self.check_deflection_on_rst1(&curpoint),
                self.check_deflection_on_rst2(&curpoint),
            )
        } else {
            (BlendStatus::Ok, BlendStatus::Ok)
        };

        if state_rst1 == BlendStatus::Backward {
            state_rst1 = BlendStatus::StepTooLarge;
            self.rebrou = true;
        }
        if state_rst2 == BlendStatus::Backward {
            state_rst2 = BlendStatus::StepTooLarge;
            self.rebrou = true;
        }
        if state_rst1 == BlendStatus::StepTooLarge || state_rst2 == BlendStatus::StepTooLarge {
            return BlendStatus::StepTooLarge;
        }

        if !self.comptra && !cur_tangent {
            let mut trarst1 = IntSurfTypeTrans::Undecided;
            let mut trarst2 = IntSurfTypeTrans::Undecided;

            if let Some(tg2drstref) = self.rst1.d1(self.sol.value(1)).map(|(_, d)| d) {
                let mut testra = func.tangent2d_on_rst1().dot(&tg2drstref);
                let mut or = TopolTool::orientation_of_curve(&self.rst1);

                if testra.abs() > self.tolpoint3d {
                    if testra < 0.0 {
                        trarst1 = conv_or_to_tra(or.reversed());
                    } else if testra > 0.0 {
                        trarst1 = conv_or_to_tra(or);
                    }

                    if let Some(tg2drstref2) =
                        self.rst2.d1(self.sol.value(2)).map(|(_, d)| d)
                    {
                        testra = func.tangent2d_on_rst2().dot(&tg2drstref2);

                        or = TopolTool::orientation_of_curve(&self.rst2);
                        if testra.abs() > self.tolpoint3d {
                            if testra < 0.0 {
                                trarst2 = conv_or_to_tra(or.reversed());
                            } else if testra > 0.0 {
                                trarst2 = conv_or_to_tra(or);
                            }
                            self.comptra = true;
                            if let Some(line) = self.line.as_mut() {
                                line.set_transitions(trarst1, trarst2);
                            }
                        }
                    }
                }
            }
        }

        if state_rst1 == BlendStatus::Ok || state_rst2 == BlendStatus::Ok {
            self.previous_p = curpoint;
            return state;
        }
        if state_rst1 == BlendStatus::StepTooSmall && state_rst2 == BlendStatus::StepTooSmall {
            self.previous_p = curpoint;
            if state == BlendStatus::Ok {
                return BlendStatus::StepTooSmall;
            }
            return state;
        }
        if state == BlendStatus::Ok {
            BlendStatus::SamePoints
        } else {
            state
        }
    }

    /// `CheckInside(Func, SituOnC1, SituOnC2, Decroch)` (`cxx:1921-1964`).
    pub(super) fn check_inside(
        &self,
        func: &dyn BlendRstRstFunction,
        situ_on_c1: &mut FaceState,
        situ_on_c2: &mut FaceState,
        decroch: &mut BlendDecrochStatus,
    ) -> bool {
        let mut tolerance = occt_math::MathVector::with_init(1, 2, 0.0);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);

        // Face pcurve 1.
        let v = self.sol.value(1);
        let (first1, last1) = (self.rst1.first_parameter(), self.rst1.last_parameter());
        if v < first1 - tolerance.value(2) || v > last1 + tolerance.value(2) {
            *situ_on_c1 = FaceState::Out;
        } else if v > first1 && v < last1 {
            *situ_on_c1 = FaceState::In;
        } else {
            *situ_on_c1 = FaceState::On;
        }

        // Face pcurve 2.
        let v = self.sol.value(2);
        let (first2, last2) = (self.rst2.first_parameter(), self.rst2.last_parameter());
        if v < first2 - tolerance.value(2) || v > last2 + tolerance.value(2) {
            *situ_on_c2 = FaceState::Out;
        } else if v > first2 && v < last2 {
            *situ_on_c2 = FaceState::In;
        } else {
            *situ_on_c2 = FaceState::On;
        }

        // Lost contact.
        let mut tg_rst1 = GpVec::zero();
        let mut n_rst1 = GpVec::zero();
        let mut tg_rst2 = GpVec::zero();
        let mut n_rst2 = GpVec::zero();
        *decroch = func.decroch(
            &self.sol,
            &mut n_rst1,
            &mut tg_rst1,
            &mut n_rst2,
            &mut tg_rst2,
        );

        *situ_on_c1 == FaceState::In
            && *situ_on_c2 == FaceState::In
            && *decroch == BlendDecrochStatus::NoDecroch
    }
}
