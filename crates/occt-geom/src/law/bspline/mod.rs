//! `Law_BSpline` (the 1D B-spline evolution) and `Law_BSpFunc` (its trimmed
//! wrapper).
//!
//! Source: `TKGeomAlgo/Law/Law_BSpline.hxx`/`Law_BSpline.cxx:199-1766` and
//! `Law_BSpFunc.hxx`/`Law_BSpFunc.cxx:33-406`.
//!
//! Split by responsibility (repo rule: a module file stays <=800 lines):
//! - this file — the `LawBSpline` fields, the constructors, `UpdateKnots` and
//!   the `PLib::SetPoles` / `PLib::GetPoles` helpers,
//! - [`modify`] — degree / multiplicity / knot / pole / weight editing,
//!   `Reverse`, `Segment` and `MovePointAndTangent`,
//! - [`query`] — the property readers (knots, poles, weights, continuities,
//!   bounds) and `LocateU`,
//! - [`eval`] — `D0`/`D1`/`D2`/`D3`/`DN`, the `Local*` window evaluators and
//!   `Resolution`,
//! - [`bsp_func`] — `Law_BSpFunc`.
//!
//! Translation notes:
//! - Every public index follows OCCT and is 1-based (`SetPole`, `Knot`,
//!   `Multiplicity`, `RemoveKnot`, ...). The stored vectors are 0-based, so the
//!   OCCT index `I` is the vector index `I - 1`.
//! - OCCT's `BSplCLib::NoMults()` / `NoWeights()` (the null handles passed by
//!   `D0`/`D1`/`D2`/`D3`/`DN`, `Law_BSpline.cxx:1236-1315`) become `None`.
//! - `UpdateKnots` shares the knot array with `FlatKnots` when the curve is
//!   uniform and non-periodic (`Law_BSpline.cxx:1092-1094`); the port stores a
//!   copy, which is the same sequence of values.
//! - `occ::handle<Law_BSpline>` becomes an owned [`LawBSpline`]; `Copy()`
//!   returns a value instead of a new handle.
//! - OCCT raises `Standard_ConstructionError` / `Standard_OutOfRange` /
//!   `Standard_DomainError` / `Standard_NoSuchObject` / `Standard_RangeError`;
//!   the port returns `Err` with that class name as the message prefix.

mod bsp_func;
mod eval;
mod modify;
mod query;

use occt_core::bspl::knot_analysis::knot_analysis;
use occt_core::bspl::knots::{knot_sequence_length, knot_sequence_periodic, nb_poles};
use occt_core::bspl::unperiodize::flat_knots_from_mults;
use occt_core::kernel::geomabs::{BSplKnotDistribution, Shape};
use occt_core::precision::{epsilon, REAL_SMALL};

pub use bsp_func::LawBSpFunc;


/// `PLib::SetPoles` (`PLib.cxx:115-137`, as called by
/// `Law_BSpline.cxx:47-59`): `(P(i) * W(i), W(i))` pairs.
fn set_poles(poles: &[f64], weights: &[f64]) -> Vec<f64> {
    let mut fp = Vec::with_capacity(2 * poles.len());
    for (i, p) in poles.iter().enumerate() {
        let w = weights[i];
        fp.push(p * w);
        fp.push(w);
    }
    fp
}

/// `PLib::GetPoles` (`PLib.cxx:63-78`, as called by `Law_BSpline.cxx:62-77`).
fn get_poles(fp: &[f64], nb_poles: usize) -> (Vec<f64>, Vec<f64>) {
    let mut poles = Vec::with_capacity(nb_poles);
    let mut weights = Vec::with_capacity(nb_poles);
    for i in 0..nb_poles {
        let w = fp[2 * i + 1];
        weights.push(w);
        poles.push(fp[2 * i] / w);
    }
    (poles, weights)
}

/// `CheckCurveData` (`Law_BSpline.cxx:82-112`).
fn check_curve_data(
    poles: &[f64],
    knots: &[f64],
    mults: &[i32],
    degree: i32,
    periodic: bool,
) -> Result<(), String> {
    if !(1..=LawBSpline::MAX_DEGREE).contains(&degree) {
        return Err("Standard_ConstructionError: Law_BSpline".into());
    }
    if poles.len() < 2 {
        return Err("Standard_ConstructionError: Law_BSpline".into());
    }
    if knots.len() != mults.len() {
        return Err("Standard_ConstructionError: Law_BSpline".into());
    }
    for i in 0..knots.len().saturating_sub(1) {
        if knots[i + 1] - knots[i] <= epsilon(knots[i].abs()) {
            return Err("Standard_ConstructionError: Law_BSpline".into());
        }
    }
    if poles.len() as i32 != nb_poles(degree, periodic, mults) {
        return Err("Standard_ConstructionError: Law_BSpline".into());
    }
    Ok(())
}

/// `Rational(const Array1OfDouble&)` (`Law_BSpline.cxx:184-196`).
fn weights_are_rational(w: &[f64]) -> bool {
    for i in 0..w.len().saturating_sub(1) {
        if (w[i] - w[i + 1]).abs() > REAL_SMALL {
            return true;
        }
    }
    false
}

/// `Law_BSpline` (`Law_BSpline.hxx:38-535`): a 1D B-spline evolution.

/// `Law_BSpline` (`Law_BSpline.hxx:38-535`): a 1D B-spline evolution.
#[derive(Debug, Clone)]
pub struct LawBSpline {
    rational: bool,
    periodic: bool,
    knot_set: BSplKnotDistribution,
    smooth: Shape,
    deg: i32,
    poles: Vec<f64>,
    weights: Option<Vec<f64>>,
    flat_knots: Vec<f64>,
    knots: Vec<f64>,
    mults: Vec<i32>,
}

impl LawBSpline {
    /// `Law_BSpline::MaxDegree()` (`Law_BSpline.cxx:304-307`) =
    /// `BSplCLib::MaxDegree()` = 25.
    pub const MAX_DEGREE: i32 = 25;

    /// `Law_BSpline::Law_BSpline` (`Law_BSpline.cxx:222-247`), non-rational.
    pub fn new(
        poles: Vec<f64>,
        knots: Vec<f64>,
        mults: Vec<i32>,
        degree: i32,
        periodic: bool,
    ) -> Result<Self, String> {
        check_curve_data(&poles, &knots, &mults, degree, periodic)?;
        let mut curve = Self {
            rational: false,
            periodic,
            knot_set: BSplKnotDistribution::NonUniform,
            smooth: Shape::CN,
            deg: degree,
            poles,
            weights: None,
            flat_knots: Vec::new(),
            knots,
            mults,
        };
        curve.update_knots();
        Ok(curve)
    }

    /// `Law_BSpline::Law_BSpline(Poles, Weights, Knots, Mults, Degree,
    /// Periodic)` (`Law_BSpline.cxx:250-300`).
    pub fn new_rational(
        poles: Vec<f64>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        mults: Vec<i32>,
        degree: i32,
        periodic: bool,
    ) -> Result<Self, String> {
        check_curve_data(&poles, &knots, &mults, degree, periodic)?;
        if weights.len() != poles.len() {
            return Err("Standard_ConstructionError: Law_BSpline".into());
        }
        for &w in &weights {
            if w <= REAL_SMALL {
                return Err("Standard_ConstructionError: Law_BSpline".into());
            }
        }
        let rational = weights_are_rational(&weights);
        let mut curve = Self {
            rational,
            periodic,
            knot_set: BSplKnotDistribution::NonUniform,
            smooth: Shape::CN,
            deg: degree,
            poles,
            weights: if rational { Some(weights) } else { None },
            flat_knots: Vec::new(),
            knots,
            mults,
        };
        curve.update_knots();
        Ok(curve)
    }

    /// `Law_BSpline::Copy` (`Law_BSpline.cxx:200-217`).
    pub fn copy(&self) -> Self {
        self.clone()
    }

    /// `Law_BSpline::UpdateKnots` (`Law_BSpline.cxx:1083-1134`): rebuilds
    /// `rational`, `knotSet`, `flatknots` and `smooth`.
    fn update_knots(&mut self) {
        self.rational = self.weights.is_some();
        let (knot_set, max_knot_mult) =
            knot_analysis(self.deg, self.periodic, &self.knots, &self.mults);
        self.knot_set = knot_set;
        if knot_set == BSplKnotDistribution::Uniform && !self.periodic {
            self.flat_knots = self.knots.clone();
        } else if self.periodic {
            let mut flat = knot_sequence_periodic(&self.knots, &self.mults, self.deg);
            flat.resize(knot_sequence_length(&self.mults, self.deg, true), 0.0);
            self.flat_knots = flat;
        } else {
            self.flat_knots = flat_knots_from_mults(&self.knots, &self.mults);
        }
        self.smooth = if max_knot_mult == 0 {
            Shape::CN
        } else {
            match self.deg - max_knot_mult {
                0 => Shape::C0,
                1 => Shape::C1,
                2 => Shape::C2,
                _ => Shape::C3,
            }
        };
    }

}