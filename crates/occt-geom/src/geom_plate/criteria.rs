//! `GeomPlate_PlateG0Criterion` and `GeomPlate_PlateG1Criterion`.
//!
//! Source: `TKGeomAlgo/GeomPlate/GeomPlate_PlateG0Criterion.cxx:33-122` and
//! `GeomPlate_PlateG1Criterion.cxx:31-133` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_Sequence<gp_XY>` is `Vec<(f64, f64)>`, `gp_XYZ` is
//!   [`GpXyz`].
//! - `P.Coefficients(1, C)->ChangeArray1()(Lower())` is
//!   [`CriterionPatch::equation`]; the `Patch` working array the criteria build
//!   from it is [`patch_coefficients`].
//! - `PLib::EvalPoly2Var` is `occt_core::bspl::plib_eval::eval_poly2var_deriv`.
//! - `P.SetCritValue(x)` is `patch.crit_value = x`; `P.CritValue()` is read
//!   back by `IsSatisfied` (`G0.cxx:118-122`, `G1.cxx:129-133`).

use occt_core::bspl::plib_eval::eval_poly2var_deriv;
use occt_core::gp::GpVec;
use occt_core::gp::GpXyz;

use crate::adv_app2var::App2VarContext;
use crate::adv_app2var::{
    patch_coefficients, Criterion, CriterionBase, CriterionPatch, CriterionRepartition,
    CriterionType,
};

/// `GeomPlate_PlateG0Criterion` (`GeomPlate_PlateG0Criterion.hxx:35-52`).
pub struct PlateG0Criterion {
    base: CriterionBase,
    /// `myData` (`hxx:50`).
    data: Vec<(f64, f64)>,
    /// `myXYZ` (`hxx:51`).
    xyz: Vec<GpXyz>,
}

impl PlateG0Criterion {
    /// `GeomPlate_PlateG0Criterion(Data, G0Data, Maximum, Type, Repart)`
    /// (`G0.cxx:33-48`).
    pub fn new(
        data: Vec<(f64, f64)>,
        g0_data: Vec<GpXyz>,
        maximum: f64,
        criterion_type: CriterionType,
        repartition: CriterionRepartition,
    ) -> Self {
        Self {
            base: CriterionBase {
                max_value: maximum,
                criterion_type,
                repartition,
            },
            data,
            xyz: g0_data,
        }
    }

    /// The header's default arguments (`hxx:42-43`): `AdvApp2Var_Absolute`,
    /// `AdvApp2Var_Regular`.
    pub fn with_defaults(data: Vec<(f64, f64)>, g0_data: Vec<GpXyz>, maximum: f64) -> Self {
        Self::new(
            data,
            g0_data,
            maximum,
            CriterionType::Absolute,
            CriterionRepartition::Regular,
        )
    }
}

impl Criterion for PlateG0Criterion {
    fn base(&self) -> &CriterionBase {
        &self.base
    }

    /// `Value(P, C)` (`G0.cxx:52-116`).
    fn value(&self, patch: &mut CriterionPatch<'_>, ctx: &App2VarContext) {
        let nb_coeff = [
            patch.nb_coeff_u.max(0) as usize,
            patch.nb_coeff_v.max(0) as usize,
        ];
        let dimension = 3 * nb_coeff[1];
        let mut coeffs = vec![0.0; nb_coeff[0] * dimension];
        patch_coefficients(patch, ctx, &mut coeffs);

        let mut dist = 0.0f64;
        for i in 0..self.data.len() {
            let (x2d, y2d) = self.data[i];
            if patch.u0 < x2d && x2d < patch.u1 && patch.v0 < y2d && y2d < patch.v1 {
                // `u, v` rescaled onto `(-1, 1)` (`G0.cxx:98-100`).
                let up = (2.0 * x2d - patch.u0 - patch.u1) / (patch.u1 - patch.u0);
                let vp = (2.0 * y2d - patch.v0 - patch.v1) / (patch.v1 - patch.v0);
                let mut digit = [0.0f64; 3];
                eval_poly2var_deriv(
                    up,
                    vp,
                    0,
                    0,
                    nb_coeff[0] as i32 - 1,
                    nb_coeff[1] as i32 - 1,
                    3,
                    &coeffs,
                    &mut digit,
                );
                let Some(point) = self.xyz.get(i) else {
                    continue;
                };
                let dx = digit[0] - point.x();
                let dy = digit[1] - point.y();
                let dz = digit[2] - point.z();
                let dist_tmp = dx * dx + dy * dy + dz * dz;
                if dist_tmp > dist {
                    dist = dist_tmp;
                }
            }
        }
        patch.crit_value = dist.sqrt();
    }

    /// `IsSatisfied(P)` (`G0.cxx:118-122`).
    fn is_satisfied(&self, patch: &CriterionPatch<'_>) -> bool {
        patch.crit_value < self.max_value()
    }
}

/// `GeomPlate_PlateG1Criterion` (`GeomPlate_PlateG1Criterion.hxx:33-50`).
pub struct PlateG1Criterion {
    base: CriterionBase,
    /// `myData` (`hxx:50`).
    data: Vec<(f64, f64)>,
    /// `myXYZ` (`hxx:51`).
    xyz: Vec<GpXyz>,
}

impl PlateG1Criterion {
    /// `GeomPlate_PlateG1Criterion(Data, G1Data, Maximum, Type, Repart)`
    /// (`G1.cxx:31-43`).
    pub fn new(
        data: Vec<(f64, f64)>,
        g1_data: Vec<GpXyz>,
        maximum: f64,
        criterion_type: CriterionType,
        repartition: CriterionRepartition,
    ) -> Self {
        Self {
            base: CriterionBase {
                max_value: maximum,
                criterion_type,
                repartition,
            },
            data,
            xyz: g1_data,
        }
    }

    /// The header's default arguments (`hxx:42-43`).
    pub fn with_defaults(data: Vec<(f64, f64)>, g1_data: Vec<GpXyz>, maximum: f64) -> Self {
        Self::new(
            data,
            g1_data,
            maximum,
            CriterionType::Absolute,
            CriterionRepartition::Regular,
        )
    }
}

impl Criterion for PlateG1Criterion {
    fn base(&self) -> &CriterionBase {
        &self.base
    }

    /// `Value(P, C)` (`G1.cxx:45-127`).
    fn value(&self, patch: &mut CriterionPatch<'_>, ctx: &App2VarContext) {
        let nb_coeff = [
            patch.nb_coeff_u.max(0) as usize,
            patch.nb_coeff_v.max(0) as usize,
        ];
        let dimension = 3 * nb_coeff[1];
        let mut coeffs = vec![0.0; nb_coeff[0] * dimension];
        patch_coefficients(patch, ctx, &mut coeffs);

        let mut ang = 0.0f64;
        for i in 0..self.data.len() {
            let Some(normal) = self.xyz.get(i) else {
                continue;
            };
            let v3h = GpVec::new(normal.x(), normal.y(), normal.z());
            let (x2d, y2d) = self.data[i];
            if patch.u0 < x2d && x2d < patch.u1 && patch.v0 < y2d && y2d < patch.v1 {
                let up = (2.0 * x2d - patch.u0 - patch.u1) / (patch.u1 - patch.u0);
                let vp = (2.0 * y2d - patch.v0 - patch.v1) / (patch.v1 - patch.v0);

                let mut digit = [0.0f64; 3];
                eval_poly2var_deriv(
                    up,
                    vp,
                    1,
                    0,
                    nb_coeff[0] as i32 - 1,
                    nb_coeff[1] as i32 - 1,
                    3,
                    &coeffs,
                    &mut digit,
                );
                let v1s = GpVec::new(digit[0], digit[1], digit[2]);

                eval_poly2var_deriv(
                    up,
                    vp,
                    0,
                    1,
                    nb_coeff[0] as i32 - 1,
                    nb_coeff[1] as i32 - 1,
                    3,
                    &coeffs,
                    &mut digit,
                );
                let v2s = GpVec::new(digit[0], digit[1], digit[2]);

                let v3s = v1s.crossed(&v2s);
                let angle = v3s.angle(&v3h);
                if angle > std::f64::consts::PI / 2.0 {
                    if (std::f64::consts::PI - angle) > ang {
                        ang = std::f64::consts::PI - angle;
                    }
                } else if angle > ang {
                    ang = angle;
                }
            }
        }
        patch.crit_value = ang;
    }

    /// `IsSatisfied(P)` (`G1.cxx:129-133`).
    fn is_satisfied(&self, patch: &CriterionPatch<'_>) -> bool {
        patch.crit_value < self.max_value()
    }
}
