//! `BRepBlend_SurfRstLineBuilder::TestArret`, `CheckDeflectionOnSurf`,
//! `CheckDeflectionOnRst` and `CheckInside`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstLineBuilder.cxx:1427-1777`
//! (OCCT 8.0.0).
//!
//! `TestArret` (`cxx:1612-1741`) decides whether the current section is
//! admissible. Unlike [`super::cs_walking`] both supports are checked: the
//! surface side by `CheckDeflectionOnSurf` (`cxx:1427-1511`) and the
//! restriction side by `CheckDeflectionOnRst` (`cxx:1515-1599`), each with
//! `CosRef3D = 0.98` (`cxx:1432`, `:1521`) and the sagitta estimate
//! `|Tg - Tgsurf|^2 * Norme / 64`. The chord guard uses `0.01 * tolpoint3d`
//! (`cxx:1458`, `:1546`), a hundredth of the 3d tolerance, where
//! `BRepBlend_CSWalking` uses the tolerance itself.
//!
//! The line transition is computed once, when `comptra` is still false and the
//! point is not a tangency point (`cxx:1679-1713`); `tras` comes from the mixed
//! product of the section tangent, the surface normal and the surface tangent,
//! `trarst` from the sign of the restriction tangents against the orientation
//! of `rst` in `domain2`.
//!
//! `CheckInside` (`cxx:1745-1777`) is the initial-section gate used by
//! `Perform` in the `Appro` case: the restriction parameter must be inside its
//! range, the `(U, V)` must classify `In` on `domain1`, and `Decroch` must be
//! `false`.

use occt_core::gp::{GpPnt2d, GpVec};
use occt_math::MathVector;

use crate::abs::Orientation;
use crate::blend::{BlendPoint, BlendStatus, BlendSurfRstFunction};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::fclass2d::FaceState;
use crate::intsurf_transition::IntSurfTypeTrans;

use super::SurfRstLineBuilder;

/// `CosRef3D` of both deflection checks (`cxx:1432`, `cxx:1521`).
const COS_REF3D: f64 = 0.98;

/// `IntSurf_TypeTrans ConvOrToTra(TopAbs_Orientation)` (`cxx:1602-1608`).
fn conv_or_to_tra(o: Orientation) -> IntSurfTypeTrans {
    if o == Orientation::Forward {
        IntSurfTypeTrans::In
    } else {
        IntSurfTypeTrans::Out
    }
}

impl SurfRstLineBuilder {
    /// `CheckDeflectionOnSurf(CurPoint)` (`cxx:1427-1511`).
    fn check_deflection_on_surf(&self, cur_point: &BlendPoint) -> BlendStatus {
        let cur_tangent = cur_point.is_tangency_point();
        let prev_tangent = self.previous_p.is_tangency_point();

        let psurf = *cur_point.point_on_s();
        let tgsurf = cur_point
            .tangent_on_s()
            .ok()
            .copied()
            .unwrap_or_else(GpVec::zero);
        let prev_p = *self.previous_p.point_on_s();
        let prev_tg = self
            .previous_p
            .tangent_on_s()
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

    /// `CheckDeflectionOnRst(CurPoint)` (`cxx:1515-1599`).
    fn check_deflection_on_rst(&self, cur_point: &BlendPoint) -> BlendStatus {
        let cur_tangent = cur_point.is_tangency_point();
        let prev_tangent = self.previous_p.is_tangency_point();

        let psurf = *cur_point.point_on_c();
        let tgsurf = cur_point
            .tangent_on_c()
            .ok()
            .copied()
            .unwrap_or_else(GpVec::zero);
        let prev_p = *self.previous_p.point_on_c();
        let prev_tg = self
            .previous_p
            .tangent_on_c()
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

    /// `TestArret(Func, TestDeflection, State)` (`cxx:1612-1741`).
    pub(super) fn test_arret(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        test_deflection: bool,
        state: BlendStatus,
    ) -> BlendStatus {
        if !func.is_solution(&self.sol, self.tolpoint3d) {
            return BlendStatus::StepTooLarge;
        }

        let cur_tangent = func.is_tangency_point();
        let pts = *func.point_on_s();
        let ptrst = *func.point_on_rst();
        let pt2drst = *func.pnt2d_on_rst();

        let mut curpoint = BlendPoint::new();
        if cur_tangent {
            curpoint.set_on_surface_and_pcurve(
                pts,
                ptrst,
                self.param,
                self.sol.value(1),
                self.sol.value(2),
                pt2drst.x(),
                pt2drst.y(),
                self.sol.value(3),
            );
        } else {
            curpoint.set_on_surface_and_pcurve_with_tangents(
                pts,
                ptrst,
                self.param,
                self.sol.value(1),
                self.sol.value(2),
                pt2drst.x(),
                pt2drst.y(),
                self.sol.value(3),
                *func.tangent_on_s(),
                *func.tangent_on_rst(),
                *func.tangent2d_on_s(),
                *func.tangent2d_on_rst(),
            );
        }

        let (mut state_s, mut state_rst) = if test_deflection {
            (
                self.check_deflection_on_surf(&curpoint),
                self.check_deflection_on_rst(&curpoint),
            )
        } else {
            (BlendStatus::Ok, BlendStatus::Ok)
        };

        if state_s == BlendStatus::Backward {
            state_s = BlendStatus::StepTooLarge;
            self.rebrou = true;
        }
        if state_rst == BlendStatus::Backward {
            state_rst = BlendStatus::StepTooLarge;
            self.rebrou = true;
        }
        if state_s == BlendStatus::StepTooLarge || state_rst == BlendStatus::StepTooLarge {
            return BlendStatus::StepTooLarge;
        }

        if !self.comptra && !cur_tangent {
            let mut tgsecs = GpVec::zero();
            let mut nors = GpVec::zero();
            func.decroch(&self.sol, &mut nors, &mut tgsecs);
            let nors = nors.normalized();
            let tgs = *func.tangent_on_s();
            let mut testra = tgsecs.dot(&nors.crossed(&tgs));
            if testra.abs() > self.tolpoint3d {
                let tras = if testra < 0.0 {
                    IntSurfTypeTrans::In
                } else if testra > 0.0 {
                    IntSurfTypeTrans::Out
                } else {
                    IntSurfTypeTrans::Undecided
                };
                let mut trarst = IntSurfTypeTrans::Undecided;
                let mut done = false;
                if let Some(tg2drstref) = self.rst.d1(self.sol.value(3)).map(|(_, d)| d) {
                    testra = func.tangent2d_on_rst().dot(&tg2drstref);
                    let or = TopolTool::orientation_of_curve(&self.rst);
                    if testra.abs() > 1.0e-8 {
                        if testra < 0.0 {
                            trarst = conv_or_to_tra(or.reversed());
                        } else if testra > 0.0 {
                            trarst = conv_or_to_tra(or);
                        }
                        done = true;
                    }
                }
                if done {
                    self.comptra = true;
                    if let Some(line) = self.line.as_mut() {
                        line.set_transitions(tras, trarst);
                    }
                }
            }
        }

        if state_s == BlendStatus::Ok || state_rst == BlendStatus::Ok {
            self.previous_p = curpoint;
            return state;
        }
        if state_s == BlendStatus::StepTooSmall && state_rst == BlendStatus::StepTooSmall {
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

    /// `CheckInside(Func, SituOnC, SituOnS, Decroch)` (`cxx:1745-1777`).
    pub(super) fn check_inside(
        &self,
        func: &dyn BlendSurfRstFunction,
        situ_on_c: &mut FaceState,
        situ_on_s: &mut FaceState,
        decroch: &mut bool,
    ) -> bool {
        let mut tolerance = MathVector::with_init(1, 3, 0.0);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);

        // Face pcurve.
        let w = self.sol.value(3);
        let first = self.rst.first_parameter();
        let last = self.rst.last_parameter();
        if w < first - tolerance.value(3) || w > last + tolerance.value(3) {
            *situ_on_c = FaceState::Out;
        } else if w > first && w < last {
            *situ_on_c = FaceState::In;
        } else {
            *situ_on_c = FaceState::On;
        }

        // Face surface.
        let p2d = GpPnt2d::new(self.sol.value(1), self.sol.value(2));
        *situ_on_s = self
            .domain1
            .borrow_mut()
            .classify(p2d, self.tolpoint2d, false);

        // Lost contact.
        let mut tg = GpVec::zero();
        let mut no = GpVec::zero();
        *decroch = func.decroch(&self.sol, &mut tg, &mut no);

        *situ_on_c == FaceState::In && *situ_on_s == FaceState::In && !*decroch
    }
}
