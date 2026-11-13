//! `IntWalk_IWalking::TestDeflection` (`IntWalk_IWalking.gxx:2600-2887`):
//! accept the step, divide it, or stop the line, from the 3D and 2D angles
//! between the last two points and from the estimated arrow.

use occt_core::gp::GpVec;
use occt_core::precision::{PCONFUSION, SQUARE_CONFUSION};

use super::types::{uv_on, StatusDeflection};
use super::ContapTheIWalking;

/// `CosRef3D` (`gxx:2607`).
const COS_REF_3D: f64 = 0.98;
/// `CosRef2D` (`gxx:2609`).
const COS_REF_2D: f64 = 0.88;
/// `MaxDivision` (`gxx:39` / formerly `:2611`, OCCT 8.0.1): max consecutive
/// bad-root subdivisions before `ComputeOpenLine` / `ComputeCloseLine` stop.
pub(super) const MAX_DIVISION: i32 = 60;

impl ContapTheIWalking {
    /// `TestDeflection(sp, Finished, UV, StatusPrecedent, NbDivision, Step,
    /// StepSign)` (`gxx:2614-2887`).
    pub(super) fn test_deflection(
        &mut self,
        sp: &mut super::super::surf_function::ContapSurfFunction,
        finished: bool,
        uv: &[f64; 2],
        status_precedent: StatusDeflection,
        nb_division: &mut i32,
        step: &mut f64,
        step_sign: i32,
    ) -> StatusDeflection {
        let mut a_status = StatusDeflection::Ok;

        let (paramu, paramv) = uv_on(&self.previous_point, self.reversed);

        let du = uv[0] - paramu;
        let dv = uv[1] - paramv;
        let duv = du * du + dv * dv;

        let corde = GpVec::from_pnts(&self.previous_point.p, sp.point());
        let norme = corde.square_magnitude();

        if norme <= 4.0 * SQUARE_CONFUSION
            && (duv <= PCONFUSION * PCONFUSION || status_precedent != StatusDeflection::Ok)
        {
            a_status = StatusDeflection::PointConfondu;
            if status_precedent == StatusDeflection::PasTropGrand {
                return StatusDeflection::ArretSurPointPrecedent;
            }
        } else {
            let cosi = corde.dot(&self.previous_d3d);
            let mut cosi2 = 0.0;

            if cosi * f64::from(step_sign) >= 0.0 {
                let a_div = self.previous_d3d.square_magnitude() * norme;
                if a_div == 0.0 {
                    return a_status;
                }
                cosi2 = cosi * cosi / a_div;
            }
            if cosi2 < COS_REF_3D {
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x()).abs();
                let step_v = (*step * self.previous_d2d.y()).abs();
                if step_u < self.vec_tolerance[0] && step_v < self.vec_tolerance[1] {
                    a_status = StatusDeflection::ArretSurPointPrecedent;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                return a_status;
            }
        }

        // `gxx:2703-2722`.
        let a_min_tol_u = 0.1 * (*step * self.previous_d2d.x()).abs();
        let a_min_tol_v = 0.1 * (*step * self.previous_d2d.y()).abs();
        let a_tol_u = if a_min_tol_u > 0.0 {
            self.vec_tolerance[0].min(a_min_tol_u)
        } else {
            self.vec_tolerance[0]
        };
        let a_tol_v = if a_min_tol_v > 0.0 {
            self.vec_tolerance[1].min(a_min_tol_v)
        } else {
            self.vec_tolerance[1]
        };

        if du.abs() < a_tol_u && dv.abs() < a_tol_v {
            return StatusDeflection::ArretSurPointPrecedent;
        }

        let mut cosi = f64::from(step_sign) * (du * self.previous_d2d.x() + dv * self.previous_d2d.y());

        if cosi < 0.0 && a_status == StatusDeflection::PointConfondu {
            return StatusDeflection::ArretSurPointPrecedent;
        }

        if sp.is_tangent() {
            return StatusDeflection::ArretSurPoint;
        }

        if *nb_division < MAX_DIVISION
            && a_status != StatusDeflection::PointConfondu
            && status_precedent != StatusDeflection::PointConfondu
        {
            let mut cosi2 = cosi * cosi / duv;
            if cosi2 < COS_REF_2D || cosi < 0.0 {
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x()).abs();
                let step_v = (*step * self.previous_d2d.y()).abs();
                if step_u < self.vec_tolerance[0] && step_v < self.vec_tolerance[1] {
                    a_status = StatusDeflection::ArretSurPointPrecedent;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                *nb_division += 1;
                return a_status;
            }

            let dir3d = match sp.direction3d() {
                Ok(v) => *v,
                Err(_) => return StatusDeflection::ArretSurPoint,
            };
            cosi = corde.dot(&dir3d);
            cosi2 = cosi * cosi / dir3d.square_magnitude() / norme;
            if cosi2 < COS_REF_3D {
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x()).abs();
                let step_v = (*step * self.previous_d2d.y()).abs();
                if step_u < self.vec_tolerance[0] && step_v < self.vec_tolerance[1] {
                    a_status = StatusDeflection::ArretSurPoint;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                return a_status;
            }

            let dir2d = match sp.direction2d() {
                Ok(v) => v,
                Err(_) => return StatusDeflection::ArretSurPoint,
            };
            cosi = du * dir2d.x() + dv * dir2d.y();
            cosi2 = cosi * cosi / duv;
            let dot_prev = dir2d.x() * self.previous_d2d.x() + dir2d.y() * self.previous_d2d.y();
            if cosi2 < COS_REF_2D || dot_prev < 0.0 {
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x()).abs();
                let step_v = (*step * self.previous_d2d.y()).abs();
                if step_u < self.vec_tolerance[0] && step_v < self.vec_tolerance[1] {
                    a_status = StatusDeflection::ArretSurPointPrecedent;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                return a_status;
            }
        }

        if !finished {
            if a_status == StatusDeflection::PointConfondu {
                // `gxx:2783-2803`.
                let step_u = (1.5 * du)
                    .abs()
                    .min(self.pas * (self.umax - self.um));
                let step_v = (1.5 * dv)
                    .abs()
                    .min(self.pas * (self.vmax - self.vm));

                let d2dx = self.previous_d2d.x().abs();
                let d2dy = self.previous_d2d.y().abs();

                if d2dx < self.vec_tolerance[0] {
                    *step = step_v / d2dy;
                } else if d2dy < self.vec_tolerance[1] {
                    *step = step_u / d2dx;
                } else {
                    *step = (step_u / d2dx).min(step_v / d2dy);
                }
            } else {
                // `gxx:2806-2885`.
                let dir3d = match sp.direction3d() {
                    Ok(v) => *v,
                    Err(_) => return StatusDeflection::ArretSurPoint,
                };
                let dir2d = match sp.direction2d() {
                    Ok(v) => v,
                    Err(_) => return StatusDeflection::ArretSurPoint,
                };

                let a = self.previous_d3d.normalized();
                let b = dir3d.normalized();
                let fleche_courante =
                    (a.x() - b.x()).powi(2) + (a.y() - b.y()).powi(2) + (a.z() - b.z()).powi(2);
                let fleche_courante = fleche_courante * norme / 64.0;

                if fleche_courante <= 0.25 * self.fleche * self.fleche {
                    let d2dx = dir2d.x().abs();
                    let d2dy = dir2d.y().abs();

                    let step_u = (1.5 * du)
                        .abs()
                        .min(self.pas * (self.umax - self.um));
                    let step_v = (1.5 * dv)
                        .abs()
                        .min(self.pas * (self.vmax - self.vm));

                    if d2dx < self.vec_tolerance[0] {
                        *step = step_v / d2dy;
                    } else if d2dy < self.vec_tolerance[1] {
                        *step = step_u / d2dx;
                    } else {
                        *step = (step_u / d2dx).min(step_v / d2dy);
                    }
                } else if fleche_courante > self.fleche * self.fleche {
                    *step /= 2.0;
                    let step_u = (*step * self.previous_d2d.x()).abs();
                    let step_v = (*step * self.previous_d2d.y()).abs();
                    if step_u < self.vec_tolerance[0] && step_v < self.vec_tolerance[1] {
                        a_status = StatusDeflection::ArretSurPointPrecedent;
                    } else {
                        a_status = StatusDeflection::PasTropGrand;
                    }
                } else {
                    let d2dx = dir2d.x().abs();
                    let d2dy = dir2d.y().abs();

                    let step_u = (1.5 * du)
                        .abs()
                        .min(self.pas * (self.umax - self.um));
                    let step_v = (1.5 * dv)
                        .abs()
                        .min(self.pas * (self.vmax - self.vm));

                    if d2dx < self.vec_tolerance[0] {
                        *step = (*step).min(step_v / d2dy);
                    } else if d2dy < self.vec_tolerance[1] {
                        *step = (*step).min(step_u / d2dx);
                    } else {
                        *step = (*step).min((step_u / d2dx).min(step_v / d2dy));
                    }
                }
            }
        }

        a_status
    }
}
