//! `BlendFunc` — the package of generic section helpers the blending
//! functions share, plus `BlendFunc_SectionShape`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc.hxx:31-71`,
//! `BlendFunc.cxx:40-152`, `BlendFunc_SectionShape.hxx:20-26` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `GetShape` (`cxx:40-87`) and `GetMinimalWeights` (`cxx:95-136`) are the
//!   `BlendFunc_SectionShape` dispatch in front of `GeomFill::GetShape` /
//!   `GeomFill::GetMinimalWeights`, which the port already has in
//!   [`occt_geom::geom_fill::section_convert`]; the two functions here keep
//!   OCCT's `switch` and delegate per arm.
//! - `NextShape` (`cxx:138-152`) maps the requested continuity to the next
//!   level, saturating at `GeomAbs_CN`.
//! - `Knots` / `Mults` are declared in `hxx:46-50` but defined nowhere in
//!   `BlendFunc.cxx` (dead declarations in OCCT 8.0.0), so nothing is ported.
//! - `ComputeNormal` (`cxx:155-215`) and `ComputeDNormal` (`cxx:219-280`) are
//!   ported on top of [`Surface::try_eval_dn`] (`GeomAdaptor_Surface::EvalDN`,
//!   `GeomAdaptor_Surface.cxx:1697-1814`) and the `CSLib` port in
//!   [`occt_geom::cs_lib`]. OCCT lets a failing `Surf->DN` throw
//!   (`Geom_UndefinedDerivative`), so both answer `Result<_, String>` with the
//!   message `try_eval_dn` reports instead of swallowing the failure.

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt2d, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_geom::cs_lib::{self, Array2Vec, NormalStatus};
use occt_geom::geom_fill::section_convert;
use occt_geom::Surface;

/// `enum BlendFunc_SectionShape` (`BlendFunc_SectionShape.hxx:20-26`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendFuncSectionShape {
    /// `BlendFunc_Rational`.
    Rational,
    /// `BlendFunc_QuasiAngular`.
    QuasiAngular,
    /// `BlendFunc_Polynomial`.
    Polynomial,
    /// `BlendFunc_Linear`.
    Linear,
}

/// `BlendFunc::GetShape(SShape, MaxAng, NbPoles, NbKnots, Degree, TConv)`
/// (`BlendFunc.cxx:40-87`). `TConv` is in-out: only the `Rational` and
/// `QuasiAngular` / `Polynomial` arms overwrite it, and `BlendFunc_Linear`
/// leaves it untouched.
pub fn get_shape(
    shape: BlendFuncSectionShape,
    max_ang: f64,
    nb_poles: &mut i32,
    nb_knots: &mut i32,
    degree: &mut i32,
    tconv: &mut ParameterisationType,
) {
    match shape {
        BlendFuncSectionShape::Rational => {
            let nb_span = (3.0 * max_ang.abs() / 2.0 / std::f64::consts::PI).ceil() as i32;
            *nb_poles = 2 * nb_span + 1;
            *nb_knots = nb_span + 1;
            *degree = 2;
            if nb_span == 1 {
                *tconv = ParameterisationType::TgtThetaOver2_1;
            } else {
                // `// QuasiAngular affin d'etre C1 (et meme beaucoup plus)`
                // (`cxx:58`).
                *nb_poles = 7;
                *nb_knots = 2;
                *degree = 6;
                *tconv = ParameterisationType::QuasiAngular;
            }
        }
        BlendFuncSectionShape::QuasiAngular => {
            *nb_poles = 7;
            *nb_knots = 2;
            *degree = 6;
            *tconv = ParameterisationType::QuasiAngular;
        }
        BlendFuncSectionShape::Polynomial => {
            *nb_poles = 8;
            *nb_knots = 2;
            *degree = 7;
            *tconv = ParameterisationType::Polynomial;
        }
        BlendFuncSectionShape::Linear => {
            *nb_poles = 2;
            *nb_knots = 2;
            *degree = 1;
        }
    }
}

/// `BlendFunc::GetMinimalWeights(SShape, TConv, MinAng, MaxAng, Weights)`
/// (`BlendFunc.cxx:95-136`): the polynomial and linear sections are not
/// rational, so every weight is one; the rational and quasi-angular ones take
/// the minimum weight of the `[0, MaxAng]` and `[0, max(PConfusion, MinAng)]`
/// arcs. OCCT raises out of `GeomConvert` when the angles are out of range and
/// the port reports that as [`occt_core::convert::ConvertError`].
pub fn get_minimal_weights(
    shape: BlendFuncSectionShape,
    tconv: ParameterisationType,
    min_ang: f64,
    max_ang: f64,
    weights: &mut [f64],
) -> Result<(), occt_core::convert::ConvertError> {
    match shape {
        BlendFuncSectionShape::Polynomial | BlendFuncSectionShape::Linear => {
            weights.fill(1.0);
            Ok(())
        }
        BlendFuncSectionShape::Rational | BlendFuncSectionShape::QuasiAngular => {
            section_convert::get_minimal_weights(tconv, min_ang, max_ang, weights)
        }
    }
}

/// `BlendFunc::NextShape(S)` (`BlendFunc.cxx:138-152`): `C0 -> C1 -> C2 -> C3`,
/// everything else (`C3`, `CN`) -> `CN`.
pub fn next_shape(s: Shape) -> Shape {
    match s {
        Shape::C0 => Shape::C1,
        Shape::C1 => Shape::C2,
        Shape::C2 => Shape::C3,
        _ => Shape::CN,
    }
}

/// `BlendFunc::ComputeNormal(Surf, p2d, Normal)` (`BlendFunc.cxx:155-215`).
///
/// `MaxOrder = 3` (`cxx:160`): `DerSurf` holds `DN(U, V, i, j)` for
/// `i = 0..=4`, `j = 0..=4` (`cxx:167-177`), `DerNUV` holds
/// `CSLib::DNNUV(i, j, DerSurf)` for `i, j = 0..=3` (`cxx:179-187`), and the
/// status comes from `CSLib::Normal(MaxOrder, DerNUV, 1.0e-9, ...)` with the
/// surface's `First/LastUParameter` / `First/LastVParameter` box
/// (`cxx:192-208`). Answers the normal only for `CSLib_Defined`
/// (`cxx:209-214`); `None` stands for OCCT's `return false`.
pub fn compute_normal(surf: &dyn Surface, p2d: &GpPnt2d) -> Result<Option<GpVec>, String> {
    const MAX_ORDER: i32 = 3;
    let u = p2d.x();
    let v = p2d.y();

    let mut der_surf = Array2Vec::new(0, MAX_ORDER + 1, 0, MAX_ORDER + 1);
    for i in 1..=MAX_ORDER + 1 {
        let d = surf.try_eval_dn(u, v, i, 0)?;
        der_surf.set_value(i, 0, d);
    }
    for i in 0..=MAX_ORDER + 1 {
        for j in 1..=MAX_ORDER + 1 {
            let d = surf.try_eval_dn(u, v, i, j)?;
            der_surf.set_value(i, j, d);
        }
    }

    let mut der_nuv = Array2Vec::new(0, MAX_ORDER, 0, MAX_ORDER);
    for i in 0..=MAX_ORDER {
        for j in 0..=MAX_ORDER {
            der_nuv.set_value(i, j, cs_lib::dnnuv(i, j, &der_surf));
        }
    }

    let (umin, umax) = surf.u_range();
    let (vmin, vmax) = surf.v_range();
    let (stat, thenormal, _order_u, _order_v) = cs_lib::normal_higher(
        MAX_ORDER, &der_nuv, 1.0e-9, u, v, umin, umax, vmin, vmax,
    );

    if stat == NormalStatus::Defined {
        Ok(thenormal.map(|d| GpVec::from_xyz(d.xyz())))
    } else {
        Ok(None)
    }
}

/// `BlendFunc::ComputeDNormal(Surf, p2d, Normal, DNu, DNv)`
/// (`BlendFunc.cxx:219-280`).
///
/// Same `DerSurf` / `DerNUV` build-up as [`compute_normal`], then
/// `DNu = CSLib::DNNormal(1, 0, DerNUV, OrderU, OrderV)` and
/// `DNv = CSLib::DNNormal(0, 1, DerNUV, OrderU, OrderV)` (`cxx:274-278`).
/// Answers `(normal, dnormal_du, dnormal_dv)` for `CSLib_Defined`.
pub fn compute_dnormal(
    surf: &dyn Surface,
    p2d: &GpPnt2d,
) -> Result<Option<(GpVec, GpVec, GpVec)>, String> {
    const MAX_ORDER: i32 = 3;
    let u = p2d.x();
    let v = p2d.y();

    let mut der_surf = Array2Vec::new(0, MAX_ORDER + 1, 0, MAX_ORDER + 1);
    for i in 1..=MAX_ORDER + 1 {
        let d = surf.try_eval_dn(u, v, i, 0)?;
        der_surf.set_value(i, 0, d);
    }
    for i in 0..=MAX_ORDER + 1 {
        for j in 1..=MAX_ORDER + 1 {
            let d = surf.try_eval_dn(u, v, i, j)?;
            der_surf.set_value(i, j, d);
        }
    }

    let mut der_nuv = Array2Vec::new(0, MAX_ORDER, 0, MAX_ORDER);
    for i in 0..=MAX_ORDER {
        for j in 0..=MAX_ORDER {
            der_nuv.set_value(i, j, cs_lib::dnnuv(i, j, &der_surf));
        }
    }

    let (umin, umax) = surf.u_range();
    let (vmin, vmax) = surf.v_range();
    let (stat, thenormal, order_u, order_v) = cs_lib::normal_higher(
        MAX_ORDER, &der_nuv, 1.0e-9, u, v, umin, umax, vmin, vmax,
    );

    if stat == NormalStatus::Defined {
        if let Some(d) = thenormal {
            let normal = GpVec::from_xyz(d.xyz());
            let dnu = cs_lib::dnnorm(1, 0, &der_nuv, order_u, order_v);
            let dnv = cs_lib::dnnorm(0, 1, &der_nuv, order_u, order_v);
            return Ok(Some((normal, dnu, dnv)));
        }
    }
    Ok(None)
}
