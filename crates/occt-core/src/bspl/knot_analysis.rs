//! `BSplCLib::KnotForm` / `BSplCLib::MultForm` / `BSplCLib::KnotAnalysis`
//! (`BSplCLib.cxx:602-753`).
//!
//! The result drives `Geom_BSplineCurve::myKnotSet` /
//! `Geom_BSplineSurface::myUKnotSet` / `myVKnotSet` (`Geom_BSplineCurve.cxx:1176`,
//! `Geom_BSplineSurface.cxx:1153,1201`), which is what
//! `GeomToStep_MakeBSplineCurveWithKnots_gen.pxx:71-87` and
//! `GeomToStep_MakeBSplineSurfaceWithKnots.cxx:127-150` turn into the STEP
//! `knot_spec` attribute.

use crate::kernel::geomabs::BSplKnotDistribution;
use crate::precision::epsilon;

/// `BSplCLib_KnotDistribution` (`BSplCLib.hxx`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum KnotDistribution {
    Uniform,
    NonUniform,
}

/// `BSplCLib_MultDistribution` (`BSplCLib.hxx`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MultDistribution {
    Constant,
    QuasiConstant,
    NonConstant,
}

/// `BSplCLib::KnotForm(Knots, FromK1, ToK2)` (`BSplCLib.cxx:602-632`) with the
/// call both `KnotAnalysis` users make, `FromK1 = 1`, `ToK2 = Knots.Length()`
/// (`cxx:701`).
pub fn knot_form(knots: &[f64]) -> KnotDistribution {
    let from_k1 = 1usize;
    if from_k1 + 1 > knots.len() {
        return KnotDistribution::Uniform;
    }
    let mut a_ui = knots[from_k1 - 1].abs();
    let mut a_uj = knots[from_k1].abs();
    let mut a_du0 = (a_uj - a_ui).abs();
    let mut an_eps = epsilon(a_ui) + epsilon(a_uj) + epsilon(a_du0);
    let mut i = from_k1 + 1;
    while i < knots.len() {
        a_ui = knots[i - 1].abs();
        a_uj = knots[i].abs();
        let a_du1 = (a_uj - a_ui).abs();
        if (a_du1 - a_du0).abs() > an_eps {
            return KnotDistribution::NonUniform;
        }
        a_du0 = a_du1;
        an_eps = epsilon(a_ui) + epsilon(a_uj) + epsilon(a_du0);
        i += 1;
    }
    KnotDistribution::Uniform
}

/// `BSplCLib::MultForm(Mults, FromK1, ToK2)` (`BSplCLib.cxx:636-688`) with
/// `FromK1 = 1`, `ToK2 = Mults.Length()` (`cxx:705`).
pub fn mult_form(mults: &[i32]) -> MultDistribution {
    let a_first = 1usize;
    let a_last = mults.len();
    if a_first + 1 > mults.len() {
        return MultDistribution::Constant;
    }
    let a_first_mult = mults[a_first - 1];
    let mut a_form = MultDistribution::Constant;
    let mut a_mult = mults[a_first];
    let mut i = a_first + 1;
    while i <= a_last && a_form != MultDistribution::NonConstant {
        if i == a_first + 1 {
            if a_mult != a_first_mult {
                a_form = MultDistribution::QuasiConstant;
            }
        } else if i == a_last {
            if a_form == MultDistribution::QuasiConstant {
                if a_first_mult != mults[i - 1] {
                    a_form = MultDistribution::NonConstant;
                }
            } else if a_mult != mults[i - 1] {
                a_form = MultDistribution::NonConstant;
            }
        } else {
            if a_mult != mults[i - 1] {
                a_form = MultDistribution::NonConstant;
            }
            a_mult = mults[i - 1];
        }
        i += 1;
    }
    a_form
}

/// `BSplCLib::KnotAnalysis(Degree, Periodic, CKnots, CMults, KnotForm,
/// MaxKnotMult)` (`BSplCLib.cxx:692-753`).
///
/// `knots` are the distinct knots and `mults` their multiplicities, exactly the
/// arrays `Geom_BSplineCurve::myKnots` / `myMults` and
/// `Geom_BSplineSurface::myUKnots` / `myUMults` hold. The distribution starts as
/// `GeomAbs_NonUniform` (`cxx:699`) and is only promoted inside the
/// `KSet == Uniform` arm, so an unmatched multiplicity pattern stays
/// `NonUniform` (`cxx:708-709`).
///
/// The second value is `MaxKnotMult`, the largest multiplicity strictly between
/// the first and last *u*-knot index (`cxx:741-752`), which OCCT feeds to
/// `mySmooth`.
pub fn knot_analysis(
    degree: i32,
    periodic: bool,
    knots: &[f64],
    mults: &[i32],
) -> (BSplKnotDistribution, i32) {
    let mut knot_form_result = BSplKnotDistribution::NonUniform;

    if knot_form(knots) == KnotDistribution::Uniform {
        match mult_form(mults) {
            MultDistribution::NonConstant => {}
            MultDistribution::Constant => {
                if knots.len() == 2 {
                    knot_form_result = BSplKnotDistribution::PiecewiseBezier;
                } else if mults.first().copied() == Some(1) {
                    knot_form_result = BSplKnotDistribution::Uniform;
                }
            }
            MultDistribution::QuasiConstant => {
                if mults.first().copied() == Some(degree + 1) {
                    let m = mults.get(1).copied().unwrap_or(0);
                    if m == degree {
                        knot_form_result = BSplKnotDistribution::PiecewiseBezier;
                    } else if m == 1 {
                        knot_form_result = BSplKnotDistribution::QuasiUniform;
                    }
                }
            }
        }
    }

    let (first_km, last_km) = if periodic {
        (1i32, knots.len() as i32)
    } else {
        (
            crate::bspl::locate::first_u_knot_index(degree, mults),
            crate::bspl::locate::last_u_knot_index(degree, mults),
        )
    };
    let mut max_knot_mult = 0i32;
    if last_km - first_km != 1 {
        let mut i = first_km + 1;
        while i < last_km {
            let multi = mults.get((i - 1) as usize).copied().unwrap_or(0);
            max_knot_mult = max_knot_mult.max(multi);
            i += 1;
        }
    }

    (knot_form_result, max_knot_mult)
}
