//! `AdvApp2Var_Criterion` (`AdvApp2Var_Criterion.hxx:34-54`,
//! `AdvApp2Var_Criterion.cxx:21-41`).
//!
//! Translation notes:
//! - The abstract base becomes the [`Criterion`] trait; its three protected
//!   members (`hxx:50-52`) live in [`CriterionBase`], which each criterion
//!   embeds and returns from [`Criterion::base`].
//! - `AdvApp2Var_CriterionType` (`AdvApp2Var_CriterionType.hxx:26-30`) and
//!   `AdvApp2Var_CriterionRepartition`
//!   (`AdvApp2Var_CriterionRepartition.hxx:26-30`) are the two enums below.
//! - `AdvApp2Var_Patch&` / `const AdvApp2Var_Patch&` become
//!   [`CriterionPatch`], the subset of patch state the criteria touch.

/// `AdvApp2Var_CriterionType` (`AdvApp2Var_CriterionType.hxx:26-30`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CriterionType {
    /// `AdvApp2Var_Absolute`: cut when the criterion is not satisfied and
    /// deactivate the maximum error computation.
    Absolute,
    /// `AdvApp2Var_Relative`: cut when the maximum error is not good, or when
    /// it is good and the criterion is not satisfied.
    Relative,
}

/// `AdvApp2Var_CriterionRepartition`
/// (`AdvApp2Var_CriterionRepartition.hxx:26-30`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CriterionRepartition {
    /// `AdvApp2Var_Regular`: all new cutting points at each step.
    Regular,
    /// `AdvApp2Var_Incremental`: one new cutting point at each step.
    Incremental,
}

/// The members of `AdvApp2Var_Criterion` (`hxx:50-52`).
#[derive(Clone, Copy, Debug)]
pub struct CriterionBase {
    /// `myMaxValue`.
    pub max_value: f64,
    /// `myType`.
    pub criterion_type: CriterionType,
    /// `myRepartition`.
    pub repartition: CriterionRepartition,
}

/// The `AdvApp2Var_Patch` state a criterion reads or writes.
///
/// Readers: `U0()`/`U1()`/`V0()`/`V1()` (`AdvApp2Var_Patch.hxx:75-81`),
/// `NbCoeffInU()`/`NbCoeffInV()` (`:91-93`), `Coefficients(1, C)`
/// (`AdvApp2Var_Patch.cxx:1150-1164`, i.e. `myEquation`).
/// Writers: `SetCritValue()` (`:113`), read back by `IsSatisfied` through
/// `CritValue()` (`:111`).
pub struct CriterionPatch<'a> {
    /// `U0()`.
    pub u0: f64,
    /// `U1()`.
    pub u1: f64,
    /// `V0()`.
    pub v0: f64,
    /// `V1()`.
    pub v1: f64,
    /// `NbCoeffInU()` (`myNbCoeffInU`).
    pub nb_coeff_u: i32,
    /// `NbCoeffInV()` (`myNbCoeffInV`).
    pub nb_coeff_v: i32,
    /// `myEquation` (`AdvApp2Var_Patch.cxx:942-957`): `u`-major, then `v`, then
    /// the `ctx.total_dimension()` coordinates, padded to `ctx.lim_u` x
    /// `ctx.lim_v` coefficients.
    pub equation: &'a [f64],
    /// `myCritValue`.
    pub crit_value: f64,
    /// `IsApproximated()` (`myApprIsDone`), read by
    /// `CutSense(Crit, NumDec)` (`AdvApp2Var_Patch.cxx:1047`).
    pub is_approximated: bool,
}

/// `class AdvApp2Var_Criterion` (`hxx:34-54`).
pub trait Criterion {
    /// `base()`: the protected members (`hxx:50-52`).
    fn base(&self) -> &CriterionBase;

    /// `Value(P, C)` (`hxx:39`).
    fn value(&self, patch: &mut CriterionPatch<'_>, ctx: &super::context::App2VarContext);

    /// `IsSatisfied(P)` (`hxx:41`).
    fn is_satisfied(&self, patch: &CriterionPatch<'_>) -> bool;

    /// `MaxValue()` (`AdvApp2Var_Criterion.cxx:23-26`).
    fn max_value(&self) -> f64 {
        self.base().max_value
    }

    /// `Type()` (`AdvApp2Var_Criterion.cxx:30-33`).
    fn criterion_type(&self) -> CriterionType {
        self.base().criterion_type
    }

    /// `Repartition()` (`AdvApp2Var_Criterion.cxx:38-41`).
    fn repartition(&self) -> CriterionRepartition {
        self.base().repartition
    }
}

/// `P.Coefficients(1, C)` copied into the `Patch` working array of the two
/// criteria (`GeomPlate_PlateG0Criterion.cxx:74-88`,
/// `GeomPlate_PlateG1Criterion.cxx:71-83`).
///
/// The criteria hard-code `Dimension == 3` (`G0.cxx:69`, `G1.cxx:66`), so the
/// `u` stride is `3 * ctx.lim_v` (`G0.cxx:80`). `out` is filled
/// `u`-major / `v` / coordinate, ready for [`eval_poly2var_deriv`].
/// `AdvApp2Var_Patch` always allocates `lim_u * lim_v * dimension` values
/// (`AdvApp2Var_Patch.cxx:836`), so the reads stay in range; anything beyond a
/// shorter array reads as zero, as `eval_poly0` does.
///
/// [`eval_poly2var_deriv`]: occt_core::bspl::plib_eval::eval_poly2var_deriv
pub fn patch_coefficients(
    patch: &CriterionPatch<'_>,
    ctx: &super::context::App2VarContext,
    out: &mut [f64],
) {
    let stride = 3 * ctx.lim_v;
    let mut ll = 0usize;
    for k1 in 0..patch.nb_coeff_u.max(0) {
        let mut pos = stride * k1;
        for _k2 in 0..patch.nb_coeff_v.max(0) {
            for d in 0..3 {
                if ll >= out.len() {
                    return;
                }
                out[ll] = patch
                    .equation
                    .get((pos + d) as usize)
                    .copied()
                    .unwrap_or(0.0);
                ll += 1;
            }
            pos += 3;
        }
    }
}
