//! `BRepBlend_Walking::TestArret` and `CheckDeflection`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Walking.cxx:777-1195`
//! (OCCT 8.0.0).
//!
//! `TestArret` decides whether the current section is admissible: it solves
//! `Function.IsSolution`, checks the deflection criterion on both supports,
//! sets the line transitions once (`comptra`, `cxx:949-972`) and maps the two
//! per-support states onto `Blend_Status`.
//!
//! `CheckDeflection` (`cxx:1046-1195`) is the criterion itself: chord length
//! against `0.01 * Tol3d`, the 3d angle against `CosRef3D = 0.98`
//! (`cxx:1048-1050`), the 2d angle against `CosRef2D = 0.88` (`cxx:1052`) and
//! the sagitta estimate `|Tg - Tgsurf|^2 * Norme / 64` (`cxx:1169-1171`).
//! The constants are the ones a test campaign set in U4; they are not
//! parameters of the algorithm and the port does not invent new ones.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_math::MathVector;

use crate::blend::{BlendFunction, BlendPoint, BlendStatus};
use crate::intsurf_transition::IntSurfTypeTrans;

use super::Walking;

impl Walking {
    /// `TestArret(F, State, TestDefl, TestSolu, TestLengthStep)`
    /// (`cxx:777-1044`).
    pub(super) fn test_arret(
        &mut self,
        func: &mut dyn BlendFunction,
        state: BlendStatus,
        test_defl: bool,
        test_solu: bool,
        test_length_step: bool,
    ) -> BlendStatus {
        let mut tolsolu = self.tolpoint3d;
        if !test_solu {
            tolsolu *= 1000.0;
        }
        let sol = self.sol.clone();
        if !func.is_solution(&sol, tolsolu) {
            return BlendStatus::StepTooLarge;
        }

        let mut curpoint = BlendPoint::new();
        let curpointistangent = func.is_tangency_point();
        let pt1: GpPnt = *func.point_on_s1();
        let pt2: GpPnt = *func.point_on_s2();
        let mut v1 = GpVec::zero();
        let mut v2 = GpVec::zero();
        let mut loctwist1 = false;
        let mut loctwist2 = false;
        if curpointistangent {
            curpoint.set_on_two_surfaces(
                pt1,
                pt2,
                self.param,
                sol.value(1),
                sol.value(2),
                sol.value(3),
                sol.value(4),
            );
        } else {
            v1 = *func.tangent_on_s1();
            v2 = *func.tangent_on_s2();
            let v12d = *func.tangent2d_on_s1();
            let v22d = *func.tangent2d_on_s2();
            curpoint.set_on_two_surfaces_with_tangents(
                pt1,
                pt2,
                self.param,
                sol.value(1),
                sol.value(2),
                sol.value(3),
                sol.value(4),
                v1,
                v2,
                v12d,
                v22d,
            );
            if func.twist_on_s1() {
                loctwist1 = true;
            }
            if func.twist_on_s2() {
                loctwist2 = true;
            }
        }

        let mut state1;
        let mut state2;
        if test_defl && self.check {
            state1 = self.check_deflection(true, &curpoint);
            state2 = self.check_deflection(false, &curpoint);
        } else {
            state1 = BlendStatus::Ok;
            state2 = BlendStatus::Ok;
            if test_length_step {
                // Only the step length is checked: `cxx:877-902`.
                let mut inf = MathVector::with_init(1, 4, 0.0);
                let mut sup = MathVector::with_init(1, 4, 0.0);
                func.get_bounds(&mut inf, &mut sup);
                sup.subtract(&inf);
                sup.multiply_scalar(0.05); // Max step: 5% of the domain

                let (curu, curv) = curpoint.parameters_on_s1().unwrap_or((0.0, 0.0));
                let (prevu, prevv) = self.previous_p.parameters_on_s1().unwrap_or((0.0, 0.0));
                if (curu - prevu).abs() > sup.value(1) {
                    state1 = BlendStatus::StepTooLarge;
                }
                if (curv - prevv).abs() > sup.value(2) {
                    state1 = BlendStatus::StepTooLarge;
                }
                let (curu, curv) = curpoint.parameters_on_s2().unwrap_or((0.0, 0.0));
                let (prevu, prevv) = self.previous_p.parameters_on_s2().unwrap_or((0.0, 0.0));
                if (curu - prevu).abs() > sup.value(3) {
                    state2 = BlendStatus::StepTooLarge;
                }
                if (curv - prevv).abs() > sup.value(4) {
                    state2 = BlendStatus::StepTooLarge;
                }
            }
        }

        if state1 == BlendStatus::Backward {
            state1 = BlendStatus::StepTooLarge;
            self.rebrou = true;
        }
        if state2 == BlendStatus::Backward {
            state2 = BlendStatus::StepTooLarge;
            self.rebrou = true;
        }
        if state1 == BlendStatus::StepTooLarge || state2 == BlendStatus::StepTooLarge {
            return BlendStatus::StepTooLarge;
        }

        // The twist can only be stated here, once the rejections are done
        // (BUC60322, `cxx:926-935`).
        if loctwist1 {
            self.twistflag1 = true;
        }
        if loctwist2 {
            self.twistflag2 = true;
        }

        if !self.comptra && !curpointistangent {
            let mut tgp1 = GpVec::zero();
            let mut tgp2 = GpVec::zero();
            let mut nor1 = GpVec::zero();
            let mut nor2 = GpVec::zero();
            func.tangent(
                sol.value(1),
                sol.value(2),
                sol.value(3),
                sol.value(4),
                &mut tgp1,
                &mut tgp2,
                &mut nor1,
                &mut nor2,
            );
            nor1.normalize();
            nor2.normalize();

            let mut testra = tgp1.dot(&nor1.crossed(&v1));
            if testra.abs() > CONFUSION {
                let mut tras1 = IntSurfTypeTrans::In;
                if (testra > 0.0 && !loctwist1) || (testra < 0.0 && loctwist1) {
                    tras1 = IntSurfTypeTrans::Out;
                }

                testra = tgp2.dot(&nor2.crossed(&v2));
                if testra.abs() > CONFUSION {
                    let mut tras2 = IntSurfTypeTrans::Out;
                    if (testra > 0.0 && !loctwist2) || (testra < 0.0 && loctwist2) {
                        tras2 = IntSurfTypeTrans::In;
                    }
                    self.comptra = true;
                    if let Some(line) = self.line.as_mut() {
                        line.set_transitions(tras1, tras2);
                    }
                }
            }
        }

        if state1 == BlendStatus::Ok || state2 == BlendStatus::Ok {
            self.previous_p = curpoint;
            return state;
        }
        if state1 == BlendStatus::StepTooSmall && state2 == BlendStatus::StepTooSmall {
            self.previous_p = curpoint;
            return if state == BlendStatus::Ok {
                BlendStatus::StepTooSmall
            } else {
                state
            };
        }
        if state == BlendStatus::Ok {
            BlendStatus::SamePoints
        } else {
            state
        }
    }

    /// `CheckDeflection(OnFirst, CurPoint)` (`cxx:1046-1195`).
    pub(super) fn check_deflection(
        &self,
        on_first: bool,
        cur_point: &BlendPoint,
    ) -> BlendStatus {
        // Rule of the U4 tests: `CosRef3D` is 11.478 degrees, `CosRef2D` 25
        // degrees (`cxx:1048-1052`).
        const COS_REF3D: f64 = 0.98;
        const COS_REF2D: f64 = 0.88;

        let curpointistangent = cur_point.is_tangency_point();
        let prevpointistangent = self.previous_p.is_tangency_point();

        let (psurf, tgsurf, prev_p, prev_tg, tolu, tolv) = if on_first {
            let tgsurf = if !curpointistangent {
                cur_point.tangent_on_s1().ok().copied()
            } else {
                None
            };
            let prev_tg = if !prevpointistangent {
                self.previous_p.tangent_on_s1().ok().copied()
            } else {
                None
            };
            (
                *cur_point.point_on_s1(),
                tgsurf,
                *self.previous_p.point_on_s1(),
                prev_tg,
                u_resolution(self.surf1.as_ref(), self.tolpoint3d),
                v_resolution(self.surf1.as_ref(), self.tolpoint3d),
            )
        } else {
            let tgsurf = if !curpointistangent {
                cur_point.tangent_on_s2().ok().copied()
            } else {
                None
            };
            let prev_tg = if !prevpointistangent {
                self.previous_p.tangent_on_s2().ok().copied()
            } else {
                None
            };
            (
                *cur_point.point_on_s2(),
                tgsurf,
                *self.previous_p.point_on_s2(),
                prev_tg,
                u_resolution(self.surf2.as_ref(), self.tolpoint3d),
                v_resolution(self.surf2.as_ref(), self.tolpoint3d),
            )
        };

        let corde = GpVec::from_pnts(&prev_p, &psurf);
        let norme = corde.square_magnitude();
        let mut prev_norme = 0.0;
        if let Some(prev_tg) = prev_tg.as_ref() {
            prev_norme = prev_tg.square_magnitude();
        }

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            return BlendStatus::SamePoints;
        }
        if let Some(prev_tg) = prev_tg.as_ref() {
            if prev_norme <= toler3d * toler3d {
                return BlendStatus::SamePoints;
            }
            let cosi = self.sens * corde.dot(prev_tg);
            if cosi < 0.0 {
                // 3d angle greater than pi/2: going backwards.
                return BlendStatus::Backward;
            }
            let cosi2 = cosi * cosi / prev_norme / norme;
            if cosi2 < COS_REF3D {
                return BlendStatus::StepTooLarge;
            }
        }

        if let Some(tgsurf) = tgsurf.as_ref() {
            let cosi = self.sens * corde.dot(tgsurf);
            let cosi2 = cosi * cosi / tgsurf.square_magnitude() / norme;
            if cosi2 < COS_REF3D || cosi < 0.0 {
                return BlendStatus::StepTooLarge;
            }
        }

        if self.check2d {
            let (curu, curv, tgon, prevu, prevv, prevd2d) = if on_first {
                let tgon = if !curpointistangent {
                    cur_point.tangent2d_on_s1().ok()
                } else {
                    None
                };
                let prevd2d = if !prevpointistangent {
                    self.previous_p.tangent2d_on_s1().ok()
                } else {
                    None
                };
                let (curu, curv) = cur_point.parameters_on_s1().unwrap_or((0.0, 0.0));
                let (prevu, prevv) = self.previous_p.parameters_on_s1().unwrap_or((0.0, 0.0));
                (curu, curv, tgon, prevu, prevv, prevd2d)
            } else {
                let tgon = if !curpointistangent {
                    cur_point.tangent2d_on_s2().ok()
                } else {
                    None
                };
                let prevd2d = if !prevpointistangent {
                    self.previous_p.tangent2d_on_s2().ok()
                } else {
                    None
                };
                let (curu, curv) = cur_point.parameters_on_s2().unwrap_or((0.0, 0.0));
                let (prevu, prevv) = self.previous_p.parameters_on_s2().unwrap_or((0.0, 0.0));
                (curu, curv, tgon, prevu, prevv, prevd2d)
            };

            let du = curu - prevu;
            let dv = curv - prevv;
            let duv = du * du + dv * dv;
            if du.abs() < tolu && dv.abs() < tolv {
                // The same point might have to be forced.
                return BlendStatus::SamePoints;
            }
            if let Some(prevd2d) = prevd2d.as_ref() {
                if prevd2d.coord.x.abs() < tolu && prevd2d.coord.y.abs() < tolv {
                    return BlendStatus::SamePoints;
                }
                let cosi = self.sens * (du * prevd2d.coord.x + dv * prevd2d.coord.y);
                if cosi < 0.0 {
                    return BlendStatus::Backward;
                }
            }
            if let Some(tgon) = tgon.as_ref() {
                let cosi = self.sens * (du * tgon.coord.x + dv * tgon.coord.y) / tgon.magnitude();
                let cosi2 = cosi * cosi / duv;
                if cosi2 < COS_REF2D || cosi < 0.0 {
                    return BlendStatus::StepTooLarge;
                }
            }
        }

        if let (Some(prev_tg), Some(tgsurf)) = (prev_tg.as_ref(), tgsurf.as_ref()) {
            // Current sagitta estimate (`cxx:1169-1171`).
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
}
