use crate::curve::Curve;
use occt_core::gp::GpPnt;

use super::distinct_knots_and_mults_of;
use occt_core::bspl::knots;

/// Non-rational or rational B-spline curve in 3D.
///
/// The distinct knots and their multiplicities are stored, as OCCT stores them
/// in `myKnots` / `myMults` (`Geom_BSplineCurve.hxx:151-152`) and rebuilds them
/// in `updateKnots()` (`Geom_BSplineCurve.cxx:1171-1188`). Every evaluation
/// (`EvalD0..D3`) feeds those arrays to `BSplCLib`; the port used to recover
/// them by rescanning `knots` on each call.
#[derive(Clone)]
pub struct GeomBSplineCurve {
    pub poles: Vec<GpPnt>,
    pub weights: Option<Vec<f64>>,
    /// Flat knot vector (`myFlatKnots`).
    pub knots: Vec<f64>,
    pub degree: usize,
    /// `myPeriodic`. Read through [`Self::periodic`]; write through
    /// [`Self::set_periodic_flag`], which refreshes [`Self::distinct_knots`].
    pub(super) periodic: bool,
    /// `myKnots` — distinct knots of the stored representation.
    pub(super) distinct: Vec<f64>,
    /// `myMults` — multiplicities, one per entry of `distinct`.
    pub(super) mults: Vec<i32>,
}

impl GeomBSplineCurve {
    /// Build a non-rational B-spline. Knot count must be poles + degree + 1.
    pub fn new(poles: Vec<GpPnt>, knots: Vec<f64>, degree: usize) -> Result<Self, &'static str> {
        knots::check_degree(poles.len(), degree, knots.len())?;
        Ok(Self::from_parts(poles, None, knots, degree, false))
    }

    /// `Geom_BSplineCurve(Poles, Knots, Mults, Degree)` — unique knots + multiplicities.
    pub fn from_poles_knots_mults(
        poles: Vec<GpPnt>,
        knots: Vec<f64>,
        mults: Vec<i32>,
        degree: usize,
    ) -> Result<Self, &'static str> {
        let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, degree as i32);
        Self::new(poles, flat, degree)
    }

    /// Build a rational B-spline (weights length must equal pole count).
    pub fn rational(poles: Vec<GpPnt>, weights: Vec<f64>, knots: Vec<f64>, degree: usize) -> Result<Self, &'static str> {
        knots::check_degree(poles.len(), degree, knots.len())?;
        if weights.len() != poles.len() {
            return Err("GeomBSplineCurve: weight count mismatch");
        }
        Ok(Self::from_parts(poles, Some(weights), knots, degree, false))
    }

    /// Assemble a curve from its flat representation — the single place that
    /// derives the stored distinct knots / multiplicities (`updateKnots()`).
    /// Used by callers that build a periodic curve from an already
    /// period-extended flat knot vector.
    pub fn from_parts(
        poles: Vec<GpPnt>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        degree: usize,
        periodic: bool,
    ) -> Self {
        let (distinct, mults) = distinct_knots_and_mults_of(&knots, degree, periodic);
        Self { poles, weights, knots, degree, periodic, distinct, mults }
    }

    /// `myPeriodic`.
    pub fn periodic(&self) -> bool {
        self.periodic
    }

    /// `myPeriodic` setter. OCCT never exposes the flag as a plain field, so
    /// the callers that flip it on a built curve go through
    /// `Geom_BSplineCurve::SetPeriodic` / `SetNotPeriodic` or restore it on a
    /// copy; either way the knot representation must follow.
    pub fn set_periodic_flag(&mut self, periodic: bool) {
        if self.periodic == periodic {
            return;
        }
        self.periodic = periodic;
        self.refresh_knot_rep();
    }

    /// `myKnots` — the distinct knots (`Geom_BSplineCurve::Knots()`).
    pub fn distinct_knots(&self) -> &[f64] {
        &self.distinct
    }

    /// `myMults` — the multiplicities (`Geom_BSplineCurve::Multiplicities()`).
    pub fn multiplicities(&self) -> &[i32] {
        &self.mults
    }

    /// Rebuild the stored distinct knots / multiplicities after `knots`,
    /// `degree` or `periodic` changed (`updateKnots()`).
    pub(super) fn refresh_knot_rep(&mut self) {
        let (distinct, mults) = distinct_knots_and_mults_of(&self.knots, self.degree, self.periodic);
        self.distinct = distinct;
        self.mults = mults;
    }

    pub fn set_pole(&mut self, i: usize, p: GpPnt) { self.poles[i] = p; }
    pub fn set_weight(&mut self, i: usize, w: f64) {
        if let Some(weights) = self.weights.as_mut() { weights[i] = w; }
    }
    pub fn pole(&self, i: usize) -> &GpPnt { &self.poles[i] }
    pub fn nb_poles(&self) -> usize { self.poles.len() }
    pub fn nb_knots(&self) -> usize { self.knots.len() }
    pub fn degree(&self) -> usize { self.degree }
    pub fn is_rational(&self) -> bool { self.weights.is_some() }

    /// `Geom_BSplineCurve::IsClosed()` (`Geom_BSplineCurve_1.cxx:146-149`):
    /// `StartPoint().SquareDistance(EndPoint()) <= Precision::Computational()`.
    pub fn is_closed(&self) -> bool {
        let d = self.d0(self.first_parameter()).square_distance(&self.d0(self.last_parameter()));
        d <= occt_core::precision::COMPUTATIONAL
    }

    /// `Geom_BSplineCurve::KnotDistribution()` (`Geom_BSplineCurve_1.cxx:363-366`),
    /// i.e. `myKnotSet` as `updateKnots()` sets it
    /// (`Geom_BSplineCurve.cxx:1176`: `BSplCLib::KnotAnalysis(myDeg, myPeriodic,
    /// myKnots, myMults, myKnotSet, MaxKnotMult)`).
    pub fn knot_distribution(&self) -> occt_core::kernel::geomabs::BSplKnotDistribution {
        occt_core::bspl::knot_analysis::knot_analysis(
            self.degree as i32,
            self.periodic,
            &self.distinct,
            &self.mults,
        )
        .0
    }

    /// `Geom_BSplineCurve::Knots()` / `Multiplicities()`
    /// (`Geom_BSplineCurve_1.cxx:380-383`, `:148-151`): the stored distinct
    /// knots and their multiplicities.
    ///
    /// The port stores the flat sequence plus the derived distinct knots /
    /// multiplicities. For a periodic curve the period-extension knots live
    /// **outside** `[FirstParameter, LastParameter]` (they are `stored_knot ±
    /// period`), so the run lengths of the flat array restricted to that
    /// interval are exactly `myMults` — the equivalent of OCCT's separately
    /// stored arrays.
    ///
    /// Returns owned copies; hot paths read [`Self::distinct_knots`] and
    /// [`Self::multiplicities`] instead.
    pub fn distinct_knots_and_mults(&self) -> (Vec<f64>, Vec<i32>) {
        (self.distinct.clone(), self.mults.clone())
    }

    /// `Geom_BSplineCurve::SetKnots(K)` (`Geom_BSplineCurve.cxx:758-765`):
    /// `CheckCurveData(myPoles, K, myMults, myDeg, myPeriodic)` followed by
    /// `updateKnots()`. `k` holds the **distinct** knots.
    pub fn set_knots(&mut self, distinct_knots: &[f64]) -> Result<(), &'static str> {
        let (_, umults) = self.distinct_knots_and_mults();
        if umults.len() != distinct_knots.len() {
            return Err("GeomBSplineCurve::set_knots: knot count mismatch");
        }
        // `CheckCurveData` (`Geom_BSplineCurve.cxx:91-94`).
        if knots::nb_poles(self.degree as i32, self.periodic, &umults) as usize != self.poles.len()
        {
            return Err("GeomBSplineCurve::set_knots: pole/degree mismatch");
        }
        self.knots = if self.periodic {
            knots::knot_sequence_periodic(distinct_knots, &umults, self.degree as i32)
        } else {
            occt_core::bspl::banded_interp::knot_sequence(
                distinct_knots,
                &umults,
                self.degree as i32,
            )
        };
        self.refresh_knot_rep();
        Ok(())
    }

    /// `Geom_BSplineCurve::FirstUKnotIndex()` / `LastUKnotIndex()`
    /// (`Geom_BSplineCurve_1.cxx:334-344`, `:404-414`): `1` / `NbKnots()` on a
    /// periodic curve, otherwise the `BSplCLib` indices of the multiplicities.
    pub(super) fn u_knot_index_range(&self, umults: &[i32]) -> (i32, i32) {
        if self.periodic {
            (1, umults.len() as i32)
        } else {
            (
                occt_core::bspl::locate::first_u_knot_index(self.degree as i32, umults).max(1),
                occt_core::bspl::locate::last_u_knot_index(self.degree as i32, umults).max(1),
            )
        }
    }

    /// `Geom_BSplineCurve::updateKnots()` (`Geom_BSplineCurve.cxx:1171-1188`):
    /// rebuild the flat knot vector from the distinct knots and their
    /// multiplicities (`BSplCLib::KnotSequence`; the periodic representation
    /// uses the period-extended sequence).
    ///
    /// OCCT's `KnotSet == GeomAbs_Uniform && !Periodic` shortcut assigns the
    /// distinct array directly; that array is what the non-periodic expansion
    /// reproduces for an all-multiplicity-one curve, so the port expands in
    /// both cases (see `set_periodic`, which relies on the same equivalence).
    pub(super) fn update_flat_knots(&mut self, uknots: &[f64], umults: &[i32]) {
        self.knots = if self.periodic {
            knots::knot_sequence_periodic(uknots, umults, self.degree as i32)
        } else {
            occt_core::bspl::banded_interp::knot_sequence(uknots, umults, self.degree as i32)
        };
        // The stored representation is always derived from the flat vector, so
        // it is exactly what a fresh rescan of `knots` yields.
        self.refresh_knot_rep();
    }

    /// `Geom_BSplineCurve::Weight(Index)` (`Geom_BSplineCurve_1.cxx:148-151`):
    /// the stored weight, unit when the curve is not rational
    /// (`Geom_BSplineCurve` keeps `BSplCLib::UnitWeights` in that case).
    pub fn weight(&self, index: i32) -> f64 {
        match &self.weights {
            Some(w) => w.get((index - 1).max(0) as usize).copied().unwrap_or(1.0),
            None => 1.0,
        }
    }
}
