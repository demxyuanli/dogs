//! `BRepBlend_CSWalking::TestArret`, `CheckDeflectionOnSurf` and
//! `CheckDeflectionOnCurv`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_CSWalking.cxx:258-587`
//! (OCCT 8.0.0).
//!
//! `TestArret` (`cxx:258-405`) decides whether the current section is
//! admissible: it asks `IsSolution`, checks the deflection criterion on the
//! surface, sets the line transition once (`comptra`, `cxx:332-351`) and maps
//! the per-support states onto `Blend_Status`. The curve-side deflection is
//! hard-coded to `Blend_StepTooSmall` there (`cxx:308`) because the point/face
//! paths had problems; `CheckDeflectionOnCurv` (`cxx:504-587`) is kept but no
//! longer called.
//!
//! `CheckDeflectionOnSurf` (`cxx:407-502`) uses `CosRef3D = 0.98`
//! (`cxx:412`) and `CosRef2D = 0.88` (`cxx:414`), the values a test campaign
//! set in U4, plus the sagitta estimate `|Tg - Tgsurf|^2 * Norme / 64`
//! (`cxx:488-489`). Unlike `BRepBlend_Walking::CheckDeflection`, the chord
//! guard uses `tolpoint3d` itself, not a hundredth of it (`cxx:427`).

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};

use crate::blend::{BlendCsFunction, BlendStatus};
use crate::intsurf_transition::IntSurfTypeTrans;

use super::CsWalking;

/// `CosRef3D` of `CheckDeflectionOnSurf` (`cxx:412`).
const COS_REF3D: f64 = 0.98;
/// `CosRef2D` of `CheckDeflectionOnSurf` (`cxx:414`).
const COS_REF2D: f64 = 0.88;
/// `CosRef3D` of `CheckDeflectionOnCurv` (`cxx:509`).
const COS_REF3D_ON_CURV: f64 = 0.98;

impl CsWalking {
    /// `TestArret(Function, Sol, TestDefl, State)` (`cxx:258-405`).
    pub(super) fn test_arret(
        &mut self,
        function: &mut dyn BlendCsFunction,
        sol: &occt_math::MathVector,
        test_defl: bool,
        state: BlendStatus,
    ) -> BlendStatus {
        if !function.is_solution(sol, self.tolpoint3d) {
            return BlendStatus::StepTooLarge;
        }

        let pt1 = *function.point_on_s();
        // `pt2` is only consumed by the curve-side check OCCT comments out
        // (`cxx:285`, `:304-308`).
        let _pt2 = *function.point_on_c();
        let pt2d = *function.pnt2d();
        let p_on_c = function.parameter_on_c();
        let v1 = *function.tangent_on_s();
        let v12d = *function.tangent2d();

        let (mut state1, mut state2) = if test_defl {
            // The curve-side check is put aside on purpose (`cxx:304-308`).
            (
                self.check_deflection_on_surf(&pt1, &pt2d, &v1, &v12d),
                BlendStatus::StepTooSmall,
            )
        } else {
            (BlendStatus::Ok, BlendStatus::Ok)
        };

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

        if !self.comptra {
            let mut tgp1 = GpVec::zero();
            let mut nor1 = GpVec::zero();
            function.tangent(pt2d.x(), pt2d.y(), &mut tgp1, &mut nor1);
            let testra = tgp1.dot(&nor1.crossed(&v1));
            if testra.abs() > CONFUSION {
                let tras = if testra < 0.0 {
                    IntSurfTypeTrans::In
                } else if testra > 0.0 {
                    IntSurfTypeTrans::Out
                } else {
                    IntSurfTypeTrans::Undecided
                };
                self.comptra = true;
                if let Some(line) = self.line.as_mut() {
                    line.set_transition(tras);
                }
            }
        }

        if state1 == BlendStatus::Ok || state2 == BlendStatus::Ok {
            self.previous_p.set_on_surface_and_curve_with_tangents(
                *function.point_on_s(),
                *function.point_on_c(),
                self.param,
                pt2d.x(),
                pt2d.y(),
                p_on_c,
                v1,
                *function.tangent_on_c(),
                v12d,
            );
            return state;
        }
        if state1 == BlendStatus::StepTooSmall && state2 == BlendStatus::StepTooSmall {
            self.previous_p.set_on_surface_and_curve_with_tangents(
                *function.point_on_s(),
                *function.point_on_c(),
                self.param,
                pt2d.x(),
                pt2d.y(),
                p_on_c,
                v1,
                *function.tangent_on_c(),
                v12d,
            );
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

    /// `CheckDeflectionOnSurf(Psurf, Ponsurf, Tgsurf, Tgonsurf)`
    /// (`cxx:407-502`).
    fn check_deflection_on_surf(
        &self,
        psurf: &GpPnt,
        ponsurf: &GpPnt2d,
        tgsurf: &GpVec,
        tgon_surf: &GpVec2d,
    ) -> BlendStatus {
        let prev_p = *self.previous_p.point_on_s();
        let prev_tg = self.previous_tangent();
        let tolu = u_resolution(self.surf.as_ref(), self.tolpoint3d);
        let tolv = v_resolution(self.surf.as_ref(), self.tolpoint3d);

        let corde = GpVec::from_pnts(&prev_p, psurf);
        let norme = corde.square_magnitude();
        let prev_norme = prev_tg.square_magnitude();

        let toler3d = self.tolpoint3d;
        if norme <= toler3d * toler3d || prev_norme <= toler3d * toler3d {
            return BlendStatus::SamePoints;
        }
        let mut cosi = self.sens * corde.dot(&prev_tg);
        if cosi < 0.0 {
            // 3d angle greater than pi/2: going backwards.
            return BlendStatus::Backward;
        }
        let cosi2 = cosi * cosi / prev_norme / norme;
        if cosi2 < COS_REF3D {
            return BlendStatus::StepTooLarge;
        }

        let (paramu, paramv) = self.previous_parameters_on_s();
        let previousd2d = self.previous_tangent2d();

        let du = ponsurf.x() - paramu;
        let dv = ponsurf.y() - paramv;
        let duv = du * du + dv * dv;
        if (du.abs() < tolu && dv.abs() < tolv)
            || (previousd2d.coord.x.abs() < tolu && previousd2d.coord.y.abs() < tolv)
        {
            // The same point might have to be forced.
            return BlendStatus::SamePoints;
        }
        cosi = self.sens * (du * previousd2d.coord.x + dv * previousd2d.coord.y);
        if cosi < 0.0 {
            return BlendStatus::Backward;
        }

        cosi = self.sens * corde.dot(tgsurf);
        let cosi2 = cosi * cosi / tgsurf.square_magnitude() / norme;
        if cosi2 < COS_REF3D || cosi < 0.0 {
            return BlendStatus::StepTooLarge;
        }

        let cosi = self.sens * (du * tgon_surf.coord.x + dv * tgon_surf.coord.y) / tgon_surf.magnitude();
        let cosi2 = cosi * cosi / duv;
        if cosi2 < COS_REF2D || cosi < 0.0 {
            return BlendStatus::StepTooLarge;
        }

        let diff = prev_tg.normalized().subtracted(&tgsurf.normalized());
        let fleche_courante = diff.square_magnitude() * norme / 64.0;

        if fleche_courante <= 0.25 * self.fleche * self.fleche {
            return BlendStatus::StepTooSmall;
        }
        if fleche_courante > self.fleche * self.fleche {
            return BlendStatus::StepTooLarge;
        }

        BlendStatus::Ok
    }

    /// `CheckDeflectionOnCurv(Pcurv, Param, Tgcurv)` (`cxx:504-587`).
    ///
    /// Defined in OCCT but no longer called: `TestArret` puts the curve-side
    /// check aside (`cxx:304-308`). Kept so the missing call stays visible.
    #[allow(dead_code)]
    fn check_deflection_on_curv(
        &self,
        pcurv: &GpPnt,
        param: f64,
        tgcurv: &GpVec,
    ) -> BlendStatus {
        let prev_p = *self.previous_p.point_on_c();
        let prev_tg = self
            .previous_p
            .tangent_on_c()
            .ok()
            .copied()
            .unwrap_or_else(GpVec::zero);
        let tolu = self.curv.resolution(self.tolpoint3d);

        let corde = GpVec::from_pnts(&prev_p, pcurv);
        let norme = corde.square_magnitude();
        let prev_norme = prev_tg.square_magnitude();

        let toler3d = self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // Changed on 95.01.10 (`cxx:529-531`).
            return BlendStatus::SamePoints;
        } else if prev_norme > toler3d * toler3d {
            let cosi = self.sens * corde.dot(&prev_tg);
            if cosi < 0.0 {
                // 3d angle greater than pi/2: going backwards.
                return BlendStatus::Backward;
            }
            let cosi2 = cosi * cosi / prev_norme / norme;
            if cosi2 < COS_REF3D_ON_CURV {
                return BlendStatus::StepTooLarge;
            }
        }

        let paramu = self.previous_p.parameter_on_c().unwrap_or(0.0);
        let du = param - paramu;
        if du.abs() < tolu {
            return BlendStatus::SamePoints;
        }

        if tgcurv.magnitude() <= self.tolpoint3d {
            // "Big bluff in the meantime" (`cxx:561`).
            return BlendStatus::SamePoints;
        }

        let cosi = self.sens * corde.dot(tgcurv);
        let cosi2 = cosi * cosi / tgcurv.square_magnitude() / norme;
        if cosi2 < COS_REF3D_ON_CURV || cosi < 0.0 {
            return BlendStatus::StepTooLarge;
        }

        if prev_norme > toler3d * toler3d {
            let diff = prev_tg.normalized().subtracted(&tgcurv.normalized());
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
