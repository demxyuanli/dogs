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

    /// `Geom_BSplineCurve::Geom_BSplineCurve(Poles, Knots, Mults, Degree,
    /// Periodic)` (`Geom_BSplineCurve.cxx:139-169`): the distinct knots and
    /// multiplicities are expanded by `updateKnots()` with the periodic
    /// extension when the flag is set
    /// ([`occt_core::bspl::knots::knot_sequence_periodic`]).
    pub fn from_poles_knots_mults_periodic(
        poles: Vec<GpPnt>,
        knots: Vec<f64>,
        mults: Vec<i32>,
        degree: usize,
        periodic: bool,
    ) -> Result<Self, &'static str> {
        let flat = if periodic {
            knots::knot_sequence_periodic(&knots, &mults, degree as i32)
        } else {
            occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, degree as i32)
        };
        knots::check_degree(poles.len(), degree, flat.len())?;
        Ok(Self::from_parts(poles, None, flat, degree, periodic))
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

    /// `Geom_BSplineCurve::IsEqual(theOther, thePreci)`
    /// (`Geom_BSplineCurve_1.cxx:802-852`): same degree, rationality, knot and
    /// pole counts, poles within `thePreci`, distinct knots within
    /// `Precision::Parametric(thePreci)` and identical multiplicities; a
    /// rational curve also compares its weights.
    pub fn is_equal(&self, other: &GeomBSplineCurve, preci: f64) -> bool {
        if self.distinct.is_empty()
            || self.poles.is_empty()
            || self.mults.is_empty()
            || self.degree != other.degree()
            || self.is_rational() != other.is_rational()
            || self.distinct.len() != other.nb_knots()
            || self.poles.len() != other.nb_poles()
        {
            return false;
        }
        for i in 0..self.poles.len() {
            let a = &self.poles[i];
            let b = other.pole(i);
            if (a.x() - b.x()).abs() > preci
                || (a.y() - b.y()).abs() > preci
                || (a.z() - b.z()).abs() > preci
            {
                return false;
            }
        }
        for i in 0..self.distinct.len() {
            if (self.distinct[i] - other.distinct_knots()[i]).abs()
                > occt_core::precision::Precision::parametric(preci, 100.0)
            {
                return false;
            }
        }
        for i in 0..self.mults.len() {
            if self.mults[i] != other.multiplicities()[i] {
                return false;
            }
        }
        if !self.is_rational() {
            return true;
        }
        for i in 0..self.poles.len() {
            let w = self.weight((i + 1) as i32);
            if (w - other.weight((i + 1) as i32)).abs() > occt_core::precision::epsilon(w) {
                return false;
            }
        }
        true
    }

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

    /// `Geom_BSplineCurve::Multiplicity(Index)` (`Geom_BSplineCurve_1.cxx:567-571`):
    /// 1-based index into `Multiplicities()`. An out-of-range index answers `0`
    /// where OCCT raises `Standard_OutOfRange`.
    pub fn multiplicity(&self, index: i32) -> i32 {
        if index < 1 || index as usize > self.mults.len() {
            return 0;
        }
        self.mults[(index - 1) as usize]
    }

    /// `Geom_BSplineCurve::FirstUKnotIndex()` (`Geom_BSplineCurve_1.cxx:334-344`).
    pub fn first_u_knot_index(&self) -> i32 {
        let mults = self.mults.clone();
        self.u_knot_index_range(&mults).0
    }

    /// `Geom_BSplineCurve::LastUKnotIndex()` (`Geom_BSplineCurve_1.cxx:404-414`).
    pub fn last_u_knot_index(&self) -> i32 {
        let mults = self.mults.clone();
        self.u_knot_index_range(&mults).1
    }

    /// `Geom_BSplineCurve::LocateU(U, ParametricTolerance, I1, I2, false)`
    /// (`Geom_BSplineCurve_1.cxx:702-750`) without knot repetition, i.e. over
    /// the stored distinct knots. Returns `(I1, I2)` in 1-based indices.
    pub fn locate_u(&self, u: f64, parametric_tolerance: f64) -> (i32, i32) {
        let uknots = &self.distinct;
        if uknots.is_empty() {
            return (0, 0);
        }
        let mut new_u = u;
        if self.periodic {
            new_u = occt_core::bspl::locate::in_period(
                new_u,
                self.first_parameter(),
                self.last_parameter(),
            );
        }
        let u_first = uknots[0];
        let u_last = uknots[uknots.len() - 1];
        let tol = parametric_tolerance.abs();
        if (new_u - u_first).abs() <= tol {
            return (1, 1);
        }
        if (new_u - u_last).abs() <= tol {
            let n = uknots.len() as i32;
            return (n, n);
        }
        if new_u < u_first {
            return (0, 1);
        }
        if new_u > u_last {
            let n = uknots.len() as i32;
            return (n, n + 1);
        }
        let mut i1 = occt_core::bspl::locate::hunt_occt(uknots, new_u);
        i1 = i1.max(1).min(uknots.len() as i32);
        while (i1 as usize) < uknots.len() && (uknots[i1 as usize] - new_u).abs() <= tol {
            i1 += 1;
        }
        if (uknots[(i1 - 1) as usize] - new_u).abs() <= tol {
            (i1, i1)
        } else {
            (i1, i1 + 1)
        }
    }

    /// `Geom_BSplineCurve::Resolution(Tolerance3D, UTolerance)`
    /// (`Geom_BSplineCurve_1.cxx:756-796`): the parametric tolerance matching a
    /// 3D tolerance. The periodic arm wraps the poles/weights through
    /// `BSplCLib::PrepareUnperiodize` first.
    pub fn resolution(&self, tolerance_3d: f64) -> f64 {
        if self.periodic {
            let (_, np) =
                occt_core::bspl::unperiodize::prepare_unperiodize(self.degree as i32, &self.mults);
            let n_old = self.poles.len();
            if n_old == 0 || np == 0 {
                return tolerance_3d / occt_core::precision::REAL_SMALL;
            }
            let new_poles: Vec<GpPnt> = (0..np).map(|i| self.poles[i % n_old].clone()).collect();
            let new_weights = self
                .weights
                .as_ref()
                .map(|w| (0..np).map(|i| w[i % n_old]).collect::<Vec<_>>());
            tolerance_3d
                * occt_core::bspl::bspline_curve_resolution(
                    &new_poles,
                    new_weights.as_deref(),
                    &self.knots,
                    self.degree as i32,
                    1.0,
                )
        } else {
            tolerance_3d
                * occt_core::bspl::bspline_curve_resolution(
                    &self.poles,
                    self.weights.as_deref(),
                    &self.knots,
                    self.degree as i32,
                    1.0,
                )
        }
    }
}
