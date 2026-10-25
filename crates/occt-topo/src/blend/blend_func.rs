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
//! - `ComputeNormal` (`cxx:155-216`) and `ComputeDNormal`
//!   (`cxx:218-279`) need `CSLib::DNNUV` / `CSLib::DN` / `CSLib::Normal`, which
//!   have no port yet. UNPORTED.

use occt_core::convert::ParameterisationType;
use occt_core::kernel::geomabs::Shape;
use occt_geom::geom_fill::section_convert;

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
