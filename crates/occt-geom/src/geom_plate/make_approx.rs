//! `GeomPlate_MakeApprox` (`GeomPlate_MakeApprox.hxx:31-87`,
//! `GeomPlate_MakeApprox.cxx:268-558`): converts a `GeomPlate_Surface` into a
//! `Geom_BSplineSurface` through `AdvApp2Var_ApproxAFunc2Var`.
//!
//! | OCCT | here |
//! |---|---|
//! | `GeomPlate_MakeApprox_Eval::Evaluate` (`cxx:37-265`) | `adv_app2var::eval_surface::evaluate` |
//! | ctor `(SurfPlate, PlateCrit, Tol3d, Nbmax, dgmax, Continuity, EnlargeCoeff)` (`cxx:268-336`) | [`MakeApprox::from_criterion`] |
//! | ctor `(SurfPlate, Tol3d, Nbmax, dgmax, dmax, CritOrder, Continuity, EnlargeCoeff)` (`cxx:338-538`) | [`MakeApprox::from_crit_order`] |
//! | `Surface` / `ApproxError` / `CriterionError` (`cxx:540-558`) | [`MakeApprox::surface`] / [`MakeApprox::approx_error`] / [`MakeApprox::criterion_error`] |
//!
//! Translation notes:
//! - `AdvApp2Var_ApproxAFunc2Var` is `adv_app2var::approximate_until_cdi`
//!   (`approximate_until_cdi` takes the surface the evaluator would sample, the
//!   way `AdvApp2Var_ApproxAFunc2Var` takes `AdvApp2Var_EvaluatorFunc2Var`).
//!   `AdvApp2Var_Context` is built by `context_until_roots` with
//!   `ifav = 2` / `iu = iv = Continuity` / `ndu = ndv = max(dgmax + 1, 2*iu + 2)`
//!   / `nb3 = 1`, which is `Init` (`AdvApp2Var_ApproxAFunc2Var.cxx:164-234`,
//!   `:207-215`).
//! - The two criterion-taking constructors pass `AdvApp2Var_Criterion&` down to
//!   `ComputePatches` (`AdvApp2Var_ApproxAFunc2Var.cxx:480-644`), which is
//!   `adv_app2var::approximate_until_cdi_with_criterion`.
//! - `nul1` / `nul2` (`cxx:288-290`) are the all-zero 1D / 2D tolerance arrays of
//!   the `nb1 = nb2 = 0` subspaces, so only `eps3D = Tol3d` and
//!   `epsfr = Tol3d` reach `AdvApp2Var_Context`.
//! - `mySurface` is a handle in OCCT and may stay null when the approximation
//!   fails; here it is an `Option`, and `surface()` mirrors the handle read.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_core::gp::GpXyz;

use crate::adv_app2var::{
    approximate_until_cdi, approximate_until_cdi_with_criterion, context_until_roots, Criterion,
    CriterionRepartition, CriterionType,
};
use crate::bspline_surface::GeomBSplineSurface;
use crate::geom_plate::{PlateG0Criterion, PlateG1Criterion, PlateSurface};
use crate::Surface;

/// `GeomPlate_MakeApprox` (`hxx:31-87`).
#[derive(Clone)]
pub struct MakeApprox {
    /// `myPlate` (`hxx:81`).
    plate: Arc<PlateSurface>,
    /// `mySurface` (`hxx:82`).
    surface: Option<GeomBSplineSurface>,
    /// `myAppError` (`hxx:83`).
    app_error: f64,
    /// `myCritError` (`hxx:84`).
    crit_error: f64,
}

/// `AdvApp2Var_ApproxAFunc2Var::Init`'s `myContInU` / `myContInV` switch
/// (`AdvApp2Var_ApproxAFunc2Var.cxx:178-206`).
fn continuity_order(continuity: Shape) -> Result<i32, String> {
    match continuity {
        Shape::C0 => Ok(0),
        Shape::C1 => Ok(1),
        Shape::C2 => Ok(2),
        _ => Err("AdvApp2Var_ApproxAFunc2Var : UContinuity Error".to_string()),
    }
}

impl MakeApprox {
    /// `GeomPlate_MakeApprox(SurfPlate, PlateCrit, Tol3d, Nbmax, dgmax,
    /// Continuity, EnlargeCoeff)` (`cxx:268-336`).
    ///
    /// `myPrec` is 0 (`cxx:305`) and the criterion reaches the approximation
    /// loop through `ComputePatches`.
    pub fn from_criterion(
        plate: Arc<PlateSurface>,
        plate_crit: &dyn Criterion,
        tol3d: f64,
        nb_max: i32,
        dg_max: i32,
        continuity: Shape,
        enlarge_coeff: f64,
    ) -> Result<Self, String> {
        let mut this = Self {
            plate,
            surface: None,
            app_error: 0.0,
            crit_error: 0.0,
        };
        this.run(
            tol3d,
            nb_max,
            dg_max,
            continuity,
            enlarge_coeff,
            0,
            Some(plate_crit),
        )?;
        Ok(this)
    }

    /// `GeomPlate_MakeApprox(SurfPlate, Tol3d, Nbmax, dgmax, dmax, CritOrder,
    /// Continuity, EnlargeCoeff)` (`cxx:338-538`).
    ///
    /// `crit_order == -1` (`cxx:424-458`) runs without a criterion; `0` and `1`
    /// build a `GeomPlate_PlateG0Criterion` / `GeomPlate_PlateG1Criterion`
    /// (`cxx:460-538`).
    pub fn from_crit_order(
        plate: Arc<PlateSurface>,
        tol3d: f64,
        nb_max: i32,
        dg_max: i32,
        dmax: f64,
        crit_order: i32,
        continuity: Shape,
        enlarge_coeff: f64,
    ) -> Result<Self, String> {
        // `cxx:350-381`: the 2D constraints of the plate and their 3D values
        // (G0 point or G1 normal). Only built when the criterion is asked for.
        let mut seq2d: Vec<(f64, f64)> = Vec::new();
        let mut seq3d: Vec<GpXyz> = Vec::new();
        if crit_order >= 0 {
            let constraints = plate.constraints();
            for p2d in &constraints {
                let (x, y) = (p2d.x(), p2d.y());
                let value = if crit_order == 0 {
                    // `cxx:359-365`: `myPlate->D0(P2d.X(), P2d.Y(), PP)`.
                    let p = plate.d0(x, y);
                    GpXyz::new(p.x(), p.y(), p.z())
                } else {
                    // `cxx:367-376`: `v3h = v1h ^ v2h`.
                    let (_, v1h, v2h) = plate.d1(x, y);
                    let v3h = v1h.crossed(&v2h);
                    GpXyz::new(v3h.x(), v3h.y(), v3h.z())
                };
                seq2d.push((x, y));
                seq3d.push(value);
            }
        }
        // `cxx:388-393`, `:396-410`: `seuil` is `10 * dmax` whenever `Tol3d`
        // is below it, for both the G0 and the G1 criterion.
        let seuil = if crit_order >= 0 && tol3d < 10.0 * dmax {
            10.0 * dmax
        } else {
            tol3d
        };
        let mut this = Self {
            plate,
            surface: None,
            app_error: 0.0,
            crit_error: 0.0,
        };
        match crit_order {
            // `cxx:424-458`: `myPrec = 1`.
            -1 => this.run(tol3d, nb_max, dg_max, continuity, enlarge_coeff, 1, None)?,
            // `cxx:460-497`.
            0 => {
                let crit = PlateG0Criterion::new(
                    seq2d,
                    seq3d,
                    seuil,
                    CriterionType::Absolute,
                    CriterionRepartition::Regular,
                );
                this.run(
                    tol3d,
                    nb_max,
                    dg_max,
                    continuity,
                    enlarge_coeff,
                    0,
                    Some(&crit),
                )?;
            }
            // `cxx:499-537`.
            1 => {
                let crit = PlateG1Criterion::new(
                    seq2d,
                    seq3d,
                    seuil,
                    CriterionType::Absolute,
                    CriterionRepartition::Regular,
                );
                this.run(
                    tol3d,
                    nb_max,
                    dg_max,
                    continuity,
                    enlarge_coeff,
                    0,
                    Some(&crit),
                )?;
            }
            // `cxx:424-537` covers only `-1`, `0` and `1`; anything else leaves
            // `mySurface` null.
            _ => {}
        }
        Ok(this)
    }

    /// The shared body of the `AdvApp2Var_ApproxAFunc2Var` construction
    /// (`cxx:281-331`, `cxx:382-458`).
    fn run(
        &mut self,
        tol3d: f64,
        nb_max: i32,
        dg_max: i32,
        continuity: Shape,
        enlarge_coeff: f64,
        precis_code: i32,
        crit: Option<&dyn Criterion>,
    ) -> Result<(), String> {
        // `cxx:279-285`, `:386-388`: `RealBounds` then the `EnlargeCoeff` scale
        // on all four parameters.
        let (u0, u1, v0, v1) = self.plate.real_bounds();
        let (u0, u1, v0, v1) = (
            enlarge_coeff * u0,
            enlarge_coeff * u1,
            enlarge_coeff * v0,
            enlarge_coeff * v1,
        );

        // `cxx:289-296`: one 3D subspace, `nul1` / `nul2` are all zeros.
        let three_d_tol = [tol3d];
        let three_d_tol_front = [[tol3d; 4]];

        // `AdvApp2Var_ApproxAFunc2Var.cxx:168-206`.
        let ord = continuity_order(continuity)?;
        let nd = (dg_max + 1).max(2 * ord + 2);
        if nd < 2 * ord + 2 {
            return Err("AdvApp2Var_ApproxAFunc2Var : UMaxDegree Error".to_string());
        }
        // `AdvApp2Var_ApproxAFunc2Var.cxx:217`.
        let precis_code = precis_code.clamp(0, 3);

        // `AdvApp2Var_ApproxAFunc2Var.cxx:218-233`.
        let ctx = context_until_roots(
            2,
            ord,
            ord,
            nd,
            nd,
            precis_code,
            0,
            0,
            1,
            &[],
            &[],
            &three_d_tol,
            &[],
            &[],
            &three_d_tol_front,
        )?;

        let result = match crit {
            Some(c) => approximate_until_cdi_with_criterion(
                self.plate.as_ref(),
                u0,
                u1,
                v0,
                v1,
                &ctx,
                nb_max,
                c,
            )?,
            None => approximate_until_cdi(self.plate.as_ref(), u0, u1, v0, v1, &ctx, nb_max)?,
        };
        // `cxx:327-329`, `:451-453`: `Surface(1)`, `MaxError(3, 1)`,
        // `CritError(3, 1)`; the `CritOrder == -1` arm stores `0` for the last.
        self.surface = Some(result.surface);
        self.app_error = result.max_error;
        self.crit_error = if crit.is_some() {
            result.criterion_error
        } else {
            0.0
        };
        Ok(())
    }

    /// `Surface()` (`cxx:540-544`).
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        self.surface.as_ref()
    }

    /// `Surface()` by value.
    pub fn into_surface(self) -> Option<GeomBSplineSurface> {
        self.surface
    }

    /// `ApproxError()` (`cxx:546-551`).
    pub fn approx_error(&self) -> f64 {
        self.app_error
    }

    /// `CriterionError()` (`cxx:553-558`).
    pub fn criterion_error(&self) -> f64 {
        self.crit_error
    }

    /// The `myPlate` handle (`hxx:81`), kept for the caller that wants the
    /// source surface back.
    pub fn plate(&self) -> &Arc<PlateSurface> {
        &self.plate
    }
}
