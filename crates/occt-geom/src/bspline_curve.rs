//! 3D B-spline curve. Source: `Geom_BSplineCurve.hxx`

use crate::curve::Curve;
use occt_core::gp::{GpPnt, GpVec, GpTrsf};
use occt_core::bspl::{knots, eval, poles, curve_tools};

/// `BSplCLib::MaxDegree()` (`BSplCLib.lxx:24-27`).
pub const MAX_DEGREE: usize = 25;

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
    periodic: bool,
    /// `myKnots` — distinct knots of the stored representation.
    distinct: Vec<f64>,
    /// `myMults` — multiplicities, one per entry of `distinct`.
    mults: Vec<i32>,
}

/// Recover the stored distinct knots / multiplicities from a flat knot vector
/// (`Geom_BSplineCurve::updateKnots`, `Geom_BSplineCurve.cxx:1171-1188`, via
/// `BSplCLib::KnotAnalysis`).
///
/// For a periodic curve the period-extension knots lie outside
/// `[FirstParameter, LastParameter]` (`stored_knot ± period`), so the run
/// lengths of the flat array restricted to that interval are exactly
/// `myMults`; the end knots are included because they equal `FirstParameter`
/// and `LastParameter`.
fn distinct_knots_and_mults_of(flat: &[f64], degree: usize, periodic: bool) -> (Vec<f64>, Vec<i32>) {
    let (uknots, umults) = knots::unique_knots_mults(flat);
    if !periodic || uknots.is_empty() {
        return (uknots, umults);
    }
    if flat.len() <= degree {
        // Malformed input (fewer knots than one clamped end): there is no
        // parameter window to restrict to, so keep the run lengths as they are.
        return (uknots, umults);
    }
    let (first, last) = (flat[degree], flat[flat.len() - 1 - degree]);
    let mut out_knots = Vec::new();
    let mut out_mults = Vec::new();
    for (k, m) in uknots.iter().zip(umults.iter()) {
        if *k < first || *k > last {
            continue;
        }
        out_knots.push(*k);
        out_mults.push(*m);
    }
    if out_knots.is_empty() {
        (uknots, umults)
    } else {
        (out_knots, out_mults)
    }
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
    fn refresh_knot_rep(&mut self) {
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

    /// `Geom_BSplineCurve::SetNotPeriodic()` (`Geom_BSplineCurve.cxx:974-1019`):
    /// `BSplCLib::PrepareUnperiodize` + `BSplCLib::Unperiodize` (`BSplCLib.cxx:2967-3080`)
    /// followed by `updateKnots()`.
    ///
    /// `Unperiodize` raises the end multiplicities to `degree + 1` by prepending
    /// and appending one period of knots, and re-indexes the poles cyclically:
    /// `NewPoles(k) = Poles((k - 1) % n_old + 1)` (`cxx:3076-3079`).
    pub fn set_not_periodic(&mut self) {
        if !self.periodic {
            return;
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let degree = self.degree as i32;
        let (new_knots, new_mults, _index) =
            occt_core::bspl::unperiodize::unperiodize_knots(degree, &uknots, &umults);
        let n_new = (new_mults.iter().sum::<i32>() - degree - 1).max(0) as usize;
        let n_old = self.poles.len();
        if n_old == 0 || n_new == 0 {
            return;
        }
        self.poles = (0..n_new).map(|k| self.poles[k % n_old]).collect();
        if let Some(w) = self.weights.as_ref() {
            let old = w.clone();
            self.weights = Some((0..n_new).map(|k| old[k % n_old]).collect());
        }
        self.knots = occt_core::bspl::unperiodize::flat_knots_from_mults(&new_knots, &new_mults);
        self.periodic = false;
        self.refresh_knot_rep();
    }

    /// `Geom_BSplineCurve::SetPeriodic()` (`Geom_BSplineCurve.cxx:777-815`):
    /// convert a non-periodic representation into the periodic one.
    ///
    /// The kept knots are `FirstUKnotIndex()..LastUKnotIndex()` of the distinct
    /// array (`Geom_BSplineCurve_1.cxx:334-344`, `:404-414`), the end
    /// multiplicities are clamped to `degree`, the pole count becomes
    /// `BSplCLib::NbPoles(degree, true, mults)` and the flat knot vector is
    /// rebuilt with the periodic `BSplCLib::KnotSequence` (`updateKnots()`).
    ///
    /// OCCT's `myPoles.Resize(1, nbp, true)` keeps the leading poles when the
    /// count shrinks; when it grows, OCCT leaves the new poles
    /// default-constructed — reproduced here as the origin (poles) and `0.0`
    /// (weights).
    ///
    /// `ClearEvalRepresentation()` has no counterpart: this port stores no
    /// evaluation cache.
    pub fn set_periodic(&mut self) {
        let (uknots, umults) = self.distinct_knots_and_mults();
        if uknots.is_empty() || umults.is_empty() {
            return;
        }
        let degree = self.degree as i32;
        // `FirstUKnotIndex()` / `LastUKnotIndex()` (`Geom_BSplineCurve_1.cxx:334-344`,
        // `:404-414`): `1` / `NbKnots()` for a periodic curve, otherwise the
        // `BSplCLib` indices of the distinct multiplicities.
        let (first, last) = if self.periodic {
            (1usize, uknots.len())
        } else {
            (
                occt_core::bspl::locate::first_u_knot_index(degree, &umults).max(1) as usize,
                occt_core::bspl::locate::last_u_knot_index(degree, &umults)
                    .max(1) as usize,
            )
        };
        let first = first.min(uknots.len());
        let last = last.min(uknots.len()).max(first);
        let uknots = uknots[first - 1..last].to_vec();
        let mut umults = umults[first - 1..last].to_vec();
        let last_idx = umults.len() - 1;
        let m = degree.min(umults[0].max(umults[last_idx]));
        umults[0] = m;
        umults[last_idx] = m;
        let nbp = knots::nb_poles(degree, true, &umults).max(0) as usize;
        if nbp < self.poles.len() {
            self.poles.truncate(nbp);
            if let Some(w) = self.weights.as_mut() {
                w.truncate(nbp);
            }
        } else if nbp > self.poles.len() {
            self.poles.resize(nbp, GpPnt::zero());
            if let Some(w) = self.weights.as_mut() {
                w.resize(nbp, 0.0);
            }
        }
        self.knots = knots::knot_sequence_periodic(&uknots, &umults, degree);
        self.periodic = true;
        self.refresh_knot_rep();
    }

    /// `Geom_BSplineCurve::FirstUKnotIndex()` / `LastUKnotIndex()`
    /// (`Geom_BSplineCurve_1.cxx:334-344`, `:404-414`): `1` / `NbKnots()` on a
    /// periodic curve, otherwise the `BSplCLib` indices of the multiplicities.
    fn u_knot_index_range(&self, umults: &[i32]) -> (i32, i32) {
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
    fn update_flat_knots(&mut self, uknots: &[f64], umults: &[i32]) {
        self.knots = if self.periodic {
            knots::knot_sequence_periodic(uknots, umults, self.degree as i32)
        } else {
            occt_core::bspl::banded_interp::knot_sequence(uknots, umults, self.degree as i32)
        };
        // The stored representation is always derived from the flat vector, so
        // it is exactly what a fresh rescan of `knots` yields.
        self.refresh_knot_rep();
    }

    /// `Geom_BSplineCurve::InsertKnots(Knots, Mults, ParametricTolerance, Add)`
    /// (`Geom_BSplineCurve.cxx:351-416`).
    ///
    /// `add_knots` are **distinct** parameters, `add_mults[i]` their
    /// multiplicities (`None` = flat insertion, one per knot). With
    /// `add = true` an existing knot's multiplicity grows by `M`; with
    /// `add = false` it is raised to `M` (OCCT's default for `InsertKnots`).
    /// The tolerance for knot equality is `max(Epsilon(U), parametric_tolerance)`.
    pub fn insert_knots(
        &mut self,
        add_knots: &[f64],
        add_mults: Option<&[i32]>,
        parametric_tolerance: f64,
        add: bool,
    ) -> Result<(), &'static str> {
        let (uknots, umults) = self.distinct_knots_and_mults();
        // `BSplCLib::PrepareInsertKnots` returning `false` is OCCT's
        // `Standard_ConstructionError`.
        let (nbpoles, _nbknots) = match occt_core::bspl::insert_knots::prepare_insert_knots(
            self.degree as i32,
            self.periodic,
            &uknots,
            &umults,
            add_knots,
            add_mults,
            parametric_tolerance,
            add,
        ) {
            Some(sizes) => sizes,
            None => return Err("Geom_BSplineCurve::InsertKnots"),
        };
        if nbpoles as usize == self.poles.len() {
            return Ok(());
        }
        let out = match occt_core::bspl::insert_knots::insert_knots(
            self.degree as i32,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &uknots,
            &umults,
            add_knots,
            add_mults,
            parametric_tolerance,
            add,
        ) {
            Some(out) => out,
            None => return Err("Geom_BSplineCurve::InsertKnots"),
        };
        self.poles = out.poles;
        self.weights = out.weights;
        self.update_flat_knots(&out.knots, &out.mults);
        Ok(())
    }

    /// `Geom_BSplineCurve::InsertKnot(U, M, ParametricTolerance, Add)`
    /// (`Geom_BSplineCurve.cxx:337-347`), defaults `M = 1`,
    /// `ParametricTolerance = 0.0`, `Add = true` (`Geom_BSplineCurve.hxx:239-242`).
    pub fn insert_knot(
        &mut self,
        u: f64,
        m: i32,
        parametric_tolerance: f64,
        add: bool,
    ) -> Result<(), &'static str> {
        self.insert_knots(&[u], Some(&[m]), parametric_tolerance, add)
    }

    /// `Geom_BSplineCurve::IncreaseMultiplicity(Index, M)`
    /// (`Geom_BSplineCurve.cxx:302-309`): `InsertKnots({Knots(Index)},
    /// {M - Mults(Index)}, Epsilon(1.), true)`.
    pub fn increase_multiplicity(&mut self, index: i32, m: i32) -> Result<(), &'static str> {
        let (uknots, umults) = self.distinct_knots_and_mults();
        let i = index - 1;
        if i < 0 || i as usize >= uknots.len() {
            return Err("GeomBSplineCurve::increase_multiplicity: index out of range");
        }
        let k = uknots[i as usize];
        let mm = m - umults[i as usize];
        self.insert_knots(&[k], Some(&[mm]), occt_core::precision::epsilon(1.0), true)
    }

    /// `Geom_BSplineCurve::IncreaseMultiplicity(I1, I2, M)`
    /// (`Geom_BSplineCurve.cxx:313-323`).
    pub fn increase_multiplicity_range(
        &mut self,
        i1: i32,
        i2: i32,
        m: i32,
    ) -> Result<(), &'static str> {
        let (uknots, umults) = self.distinct_knots_and_mults();
        if i1 < 1 || i2 < i1 || i2 as usize > uknots.len() {
            return Err("GeomBSplineCurve::increase_multiplicity_range: index out of range");
        }
        let ks = uknots[(i1 - 1) as usize..i2 as usize].to_vec();
        let ms: Vec<i32> = umults[(i1 - 1) as usize..i2 as usize]
            .iter()
            .map(|u| m - u)
            .collect();
        self.insert_knots(&ks, Some(&ms), occt_core::precision::epsilon(1.0), true)
    }

    /// `Geom_BSplineCurve::IncrementMultiplicity(I1, I2, Step)`
    /// (`Geom_BSplineCurve.cxx:327-333`).
    pub fn increment_multiplicity(
        &mut self,
        i1: i32,
        i2: i32,
        step: i32,
    ) -> Result<(), &'static str> {
        let (uknots, _) = self.distinct_knots_and_mults();
        if i1 < 1 || i2 < i1 || i2 as usize > uknots.len() {
            return Err("GeomBSplineCurve::increment_multiplicity: index out of range");
        }
        let ks = uknots[(i1 - 1) as usize..i2 as usize].to_vec();
        let ms = vec![step; ks.len()];
        self.insert_knots(&ks, Some(&ms), occt_core::precision::epsilon(1.0), true)
    }

    /// `Geom_BSplineCurve::SetOrigin(Index)` (`Geom_BSplineCurve.cxx:819-909`):
    /// rotate a periodic curve so that the knot `Index` becomes the origin
    /// (the knots after it stay, the ones before it move one period up and the
    /// poles are rotated to match).
    pub fn set_origin(&mut self, index: i32) -> Result<(), &'static str> {
        if !self.periodic {
            return Err("Geom_BSplineCurve::SetOrigin");
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let (first, last) = self.u_knot_index_range(&umults);
        if index < first || index > last {
            return Err("Geom_BSplineCurve::SetOrigin");
        }
        let nbknots = uknots.len() as i32;
        let nbpoles = self.poles.len() as i32;
        if nbknots == 0 || nbpoles == 0 {
            return Err("Geom_BSplineCurve::SetOrigin");
        }

        // set the knots and mults
        let period = uknots[(last - 1) as usize] - uknots[(first - 1) as usize];
        let mut newknots = Vec::with_capacity(nbknots as usize);
        let mut newmults = Vec::with_capacity(nbknots as usize);
        for i in index..=last {
            newknots.push(uknots[(i - 1) as usize]);
            newmults.push(umults[(i - 1) as usize]);
        }
        for i in (first + 1)..=index {
            newknots.push(uknots[(i - 1) as usize] + period);
            newmults.push(umults[(i - 1) as usize]);
        }

        let mut pole_first = 1i32;
        for i in (first + 1)..=index {
            pole_first += umults[(i - 1) as usize];
        }

        // set the poles and weights
        let mut newpoles = Vec::with_capacity(nbpoles as usize);
        for i in pole_first..=nbpoles {
            newpoles.push(self.poles[(i - 1) as usize]);
        }
        for i in 1..pole_first {
            newpoles.push(self.poles[(i - 1) as usize]);
        }
        if let Some(w) = self.weights.as_ref() {
            let mut newweights = Vec::with_capacity(nbpoles as usize);
            for i in pole_first..=nbpoles {
                newweights.push(w[(i - 1) as usize]);
            }
            for i in 1..pole_first {
                newweights.push(w[(i - 1) as usize]);
            }
            self.weights = Some(newweights);
        }

        self.poles = newpoles;
        self.update_flat_knots(&newknots, &newmults);
        Ok(())
    }

    /// `Geom_BSplineCurve::SetOrigin(U, Tol)` (`Geom_BSplineCurve.cxx:913-970`):
    /// move the origin of a periodic curve to the parameter `U` (translating
    /// the whole knot vector if `U` differs from the folded value by more than
    /// `Tol`, inserting a knot at `U` when needed).
    pub fn set_origin_u(&mut self, u: f64, tol: f64) -> Result<(), &'static str> {
        if !self.periodic {
            return Err("Geom_BSplineCurve::SetOrigin");
        }
        // Is U within the period?
        let mut uf = self.first_parameter();
        let mut ul = self.last_parameter();
        let period = ul - uf;
        let mut uu = u;
        while tol < (uf - uu) {
            uu += period;
        }
        while tol > (ul - uu) {
            uu -= period;
        }

        if (u - uu).abs() > tol {
            // Reparametrize the curve
            let delta = u - uu;
            uf += delta;
            ul += delta;
            for k in self.knots.iter_mut() {
                *k += delta;
            }
            self.refresh_knot_rep();
        }
        // For a periodic curve, uf and ul represent the same point
        if (u - uf).abs() < tol || (u - ul).abs() < tol {
            return Ok(());
        }

        let (uknots, _) = self.distinct_knots_and_mults();
        let mut ik = 0i32;
        let mut delta = f64::MAX;
        for (i, k) in uknots.iter().enumerate() {
            let dki = k - u;
            if dki.abs() < delta.abs() {
                ik = i as i32 + 1;
                delta = dki;
            }
        }
        if delta.abs() > tol {
            // `InsertKnot(U)`: `M = 1`, `ParametricTolerance = 0.0`, `Add = true`.
            self.insert_knot(u, 1, 0.0, true)?;
            if delta < 0.0 {
                ik += 1;
            }
        }
        self.set_origin(ik)
    }

    /// `Geom_BSplineCurve::Segment(U1, U2, theTolerance)`
    /// (`Geom_BSplineCurve.cxx:527-715`): restrict the curve to `[U1, U2]`
    /// (`U2 < U1` and, for a periodic curve, `(U2 - U1) - Period >
    /// Precision::PConfusion()` are `Standard_DomainError`).
    ///
    /// Rust has no default arguments: callers reproduce OCCT's
    /// `theTolerance = Precision::PConfusion()` default explicitly.
    pub fn segment(&mut self, u1: f64, u2: f64, tolerance: f64) -> Result<(), &'static str> {
        if u2 < u1 {
            return Err("Geom_BSplineCurve::Segment");
        }

        let was_periodic = self.periodic;
        let new_u1;
        let mut new_u2;
        let mut du = 0.0;
        let mut ddu = 0.0;

        // define param distance to keep (eap, Apr 18 2002, occ311)
        if self.periodic {
            let period = self.last_parameter() - self.first_parameter();
            du = u2 - u1;
            if du - period > occt_core::precision::PCONFUSION {
                return Err("Geom_BSplineCurve::Segment");
            }
            if du > period {
                du = period;
            }
            ddu = du;
        }

        // `BSplCLib::LocateParameter(Degree, Knots, Mults, U, Periodic,
        // Knots.Lower(), Knots.Upper(), index, NewU)` (`BSplCLib.cxx:168-185`,
        // flat-knot form) — the port's `locate_parameter_range` with the
        // **distinct** knots and their full range, **not** the convenience
        // `locate_parameter` (which derives the range from the multiplicities).
        let (uknots, _umults) = self.distinct_knots_and_mults();
        if uknots.is_empty() {
            return Err("Geom_BSplineCurve::Segment");
        }
        let (lo, hi) = (1i32, uknots.len() as i32);
        let (_, nu1) = occt_core::bspl::locate::locate_parameter_range(
            &uknots,
            u1,
            self.periodic,
            lo,
            hi,
            uknots[0],
            uknots[uknots.len() - 1],
        );
        let (_, nu2) = occt_core::bspl::locate::locate_parameter_range(
            &uknots,
            u2,
            self.periodic,
            lo,
            hi,
            uknots[0],
            uknots[uknots.len() - 1],
        );
        new_u1 = nu1;
        new_u2 = nu2;

        let nu2_val = new_u2;

        let knots_pair = [new_u1.min(new_u2), new_u1.max(new_u2)];
        let mults_pair = [self.degree as i32, self.degree as i32];

        let abs_u_max = new_u1
            .abs()
            .max(new_u2.abs())
            .max(self.first_parameter().abs())
            .max(self.last_parameter().abs());
        let eps = occt_core::precision::epsilon(abs_u_max).max(tolerance);

        // `InsertKnots(Knots, Mults, Eps)` — `Add` defaults to **false**
        // (`Geom_BSplineCurve.hxx:262-265`).
        self.insert_knots(&knots_pair, Some(&mults_pair), eps, false)?;

        if self.periodic {
            // set the origin at NewU1
            let (uk, _) = self.distinct_knots_and_mults();
            let (mut index, u) = occt_core::bspl::locate::locate_parameter_range(
                &uk,
                u1,
                true,
                lo,
                uk.len() as i32,
                uk[0],
                uk[uk.len() - 1],
            );
            // Test if the insertion is OK, shift otherwise.
            if (index as usize) < uk.len() && (uk[index as usize] - u).abs() <= eps {
                index += 1;
            }
            self.set_origin(index)?;
            self.set_not_periodic();
            new_u2 = new_u1 + du;
        }

        // compute index1 and index2 to set the new knots and mults
        let (uk2, um2) = self.distinct_knots_and_mults();
        if uk2.is_empty() {
            return Err("Geom_BSplineCurve::Segment");
        }
        let from_u1 = 1i32;
        let to_u2 = uk2.len() as i32;
        let (mut index1, ua) = occt_core::bspl::locate::locate_parameter_range(
            &uk2,
            new_u1,
            self.periodic,
            from_u1,
            to_u2,
            uk2[0],
            uk2[uk2.len() - 1],
        );
        if (index1 as usize) < uk2.len() && (uk2[index1 as usize] - ua).abs() <= eps {
            index1 += 1;
        }

        let (mut index2, ub) = occt_core::bspl::locate::locate_parameter_range(
            &uk2,
            new_u2,
            self.periodic,
            from_u1,
            to_u2,
            uk2[0],
            uk2[uk2.len() - 1],
        );
        if (index2 as usize) < uk2.len() && (uk2[index2 as usize] - ub).abs() <= eps
            || index2 == index1
        {
            index2 += 1;
        }

        let nbknots = index2 - index1 + 1;
        if nbknots <= 0 {
            return Err("Geom_BSplineCurve::Segment");
        }
        let mut nknots = Vec::with_capacity(nbknots as usize);
        let mut nmults = Vec::with_capacity(nbknots as usize);

        // to restore changed U1
        if du > 0.0 {
            // if was periodic
            du = new_u1 - u1;
        }

        for i in index1..=index2 {
            nknots.push(uk2[(i - 1) as usize] - du);
            nmults.push(um2[(i - 1) as usize]);
        }
        let last_k = nbknots as usize - 1;
        nmults[0] = self.degree as i32 + 1;
        nmults[last_k] = self.degree as i32 + 1;

        // compute index1 and index2 to set the new poles and weights
        let mut pindex1 = knots::pole_index(self.degree as i32, index1, self.periodic, &um2);
        let mut pindex2 = knots::pole_index(self.degree as i32, index2, self.periodic, &um2);
        pindex1 += 1;
        pindex2 = (pindex2 + 1).min(self.poles.len() as i32);
        if pindex1 < 1 || pindex2 < pindex1 {
            return Err("Geom_BSplineCurve::Segment");
        }

        let mut newpoles = Vec::with_capacity((pindex2 - pindex1 + 1) as usize);
        for i in pindex1..=pindex2 {
            newpoles.push(self.poles[(i - 1) as usize]);
        }
        let mut newweights = self.weights.as_ref().map(|w| {
            let mut v = Vec::with_capacity((pindex2 - pindex1 + 1) as usize);
            for i in pindex1..=pindex2 {
                v.push(w[(i - 1) as usize]);
            }
            v
        });

        if was_periodic {
            nknots[0] = u1;
            if nu2_val < u2 {
                nknots[last_k] = u1 + ddu;
            }
        }

        self.poles = newpoles;
        if self.weights.is_some() {
            self.weights = newweights.take();
        }        self.update_flat_knots(&nknots, &nmults);
        Ok(())
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

    /// `Geom_BSplineCurve::IncreaseDegree(Degree)`
    /// (`Geom_BSplineCurve.cxx:243-298`): raise the degree to `degree` using
    /// `BSplCLib::IncreaseDegree` (the Prautzsch degree elevation,
    /// `BSplCLib.cxx:2592-2963`), wrapped by `BSplCLib_IncreaseDegree`
    /// (`BSplCLib_CurveComputation.pxx:537-589`).
    ///
    /// `degree == self.degree` is a no-op, `degree < self.degree` or
    /// `degree > MaxDegree()` is OCCT's `Standard_ConstructionError`
    /// (`Geom_BSplineCurve.cxx:252-255`). `MaxDegree()` is
    /// `BSplCLib::MaxDegree()` = 25 (`BSplCLib.lxx:24-27`).
    /// `ClearEvalRepresentation()` has no counterpart: no eval cache is kept.
    pub fn increase_degree(&mut self, degree: usize) -> Result<(), &'static str> {
        if degree == self.degree {
            return Ok(());
        }
        if degree < self.degree || degree > MAX_DEGREE {
            return Err("GeomBSplineCurve::increase_degree: bad degree value");
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let (from_k1, to_k2) = self.u_knot_index_range(&umults);
        let step = degree as i32 - self.degree as i32;
        // `Geom_BSplineCurve.cxx:256-261`: `FromK1`/`ToK2` and the new pole
        // count `myPoles.Length() + Step * (ToK2 - FromK1)`.
        let nb_new_poles = self.poles.len() as i32 + step * (to_k2 - from_k1);
        // `Geom_BSplineCurve.cxx:263`: `BSplCLib::IncreaseDegreeCountKnots`.
        let nb_new_knots = occt_core::bspl::increase_degree::increase_degree_count_knots(
            self.degree as i32,
            degree as i32,
            self.periodic,
            &umults,
        );
        if nb_new_poles <= 0 || nb_new_knots <= 0 {
            return Err("GeomBSplineCurve::increase_degree: bad degree value");
        }
        let out = occt_core::bspl::increase_degree::increase_degree(
            self.degree as i32,
            degree as i32,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &uknots,
            &umults,
            nb_new_poles as usize,
            nb_new_knots as usize,
        );
        self.poles = out.poles;
        self.weights = out.weights;
        self.degree = degree;
        self.update_flat_knots(&out.knots, &out.mults);
        Ok(())
    }

    /// `Geom_BSplineCurve::RemoveKnot(Index, M, Tolerance)`
    /// (`Geom_BSplineCurve.cxx:420-495`): lower the multiplicity of the
    /// distinct knot `index` to `m` (`m == 0` removes the knot entirely).
    ///
    /// Returns `Ok(false)` exactly where `BSplCLib::RemoveKnot`
    /// (`BSplCLib.cxx:2355-2542`) reports failure; the out-of-range index is
    /// OCCT's `Standard_OutOfRange` (`Geom_BSplineCurve.cxx:430-437`).
    pub fn remove_knot(&mut self, index: i32, m: i32, tolerance: f64) -> Result<bool, &'static str> {
        if m < 0 {
            return Ok(true);
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let (i1, i2) = self.u_knot_index_range(&umults);
        if !self.periodic && (index <= i1 || index >= i2) {
            return Err("GeomBSplineCurve::remove_knot: index out of range");
        } else if self.periodic && (index < i1 || index > i2) {
            return Err("GeomBSplineCurve::remove_knot: index out of range");
        }
        let idx = (index - 1) as usize;
        if idx >= umults.len() {
            return Err("GeomBSplineCurve::remove_knot: index out of range");
        }
        let step = umults[idx] - m;
        if step <= 0 {
            return Ok(true);
        }
        let nb_new_poles = self.poles.len() as i32 - step;
        let nb_new_knots = uknots.len() as i32 - if m == 0 { 1 } else { 0 };
        if nb_new_poles < 0 || nb_new_knots < 0 {
            return Err("GeomBSplineCurve::remove_knot: index out of range");
        }
        let out = occt_core::bspl::remove_knot::remove_knot(
            index,
            m,
            self.degree as i32,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &uknots,
            &umults,
            nb_new_poles as usize,
            nb_new_knots as usize,
            tolerance,
        );
        match out {
            Some(out) => {
                self.poles = out.poles;
                self.weights = out.weights;
                self.update_flat_knots(&out.knots, &out.mults);
                Ok(true)
            }
            None => Ok(false),
        }
    }

    /// Simple degree reduction: drop to degree-1 by removing end knots and
    /// re-interpolating at the new Greville abscissae.
    pub fn decrease_degree(&mut self, _tolerance: f64) {
        if self.degree <= 1 { return; }
        let new_degree = self.degree - 1;
        let new_knots = self.knots[1..self.knots.len() - 1].to_vec();
        let n_new = self.nb_poles() - 1;
        let params = poles::greville_abscissae(&new_knots, new_degree, n_new);
        let mut new_poles = Vec::with_capacity(n_new);
        for &u in &params {
            new_poles.push(self.d0(u));
        }
        self.poles = new_poles;
        self.knots = new_knots;
        self.degree = new_degree;
        if self.weights.is_some() {
            self.weights = Some(vec![1.0; n_new]);
        }
        self.refresh_knot_rep();
    }

    /// `Geom_BSplineCurve::EvalD0/D1/D2/D3` preamble
    /// (`Geom_BSplineCurve_1.cxx:182-190`).
    ///
    /// OCCT hands `BSplCLib::D0/D1/D2/D3` the *distinct* knots plus `&myMults`
    /// (not `myFlatKnots` / `NoMults`), together with a span index located in
    /// that distinct-knot array and post-corrected by `aSpanIndex--` when
    /// `aNewU` sits below `myKnots(aSpanIndex)`. The flat-knot arm is only used
    /// by `LocalD0..LocalD3` and `EvalDN`.
    fn eval_span(&self, u: f64) -> (f64, i32, &[f64], &[i32]) {
        let uknots = self.distinct.as_slice();
        let umults = self.mults.as_slice();
        let (mut span, new_u) = occt_core::bspl::locate::locate_parameter(
            self.degree as i32,
            uknots,
            Some(umults),
            u,
            self.periodic,
            0,
        );
        if span >= 1 {
            if let Some(&k) = uknots.get(span as usize - 1) {
                if new_u < k {
                    span -= 1;
                }
            }
        }
        // Debug builds re-derive the pair and compare, so a mutation path that
        // forgets to refresh the representation shows up in the existing parity
        // tests instead of silently shifting geometry.
        debug_assert!(
            {
                let (d, m) = distinct_knots_and_mults_of(&self.knots, self.degree, self.periodic);
                d.as_slice() == uknots && m.as_slice() == umults
            },
            "GeomBSplineCurve: stored distinct knots out of date"
        );
        (new_u, span, uknots, umults)
    }

    /// Finite-difference first derivative: (f(u+h) - f(u-h)) / (2h).
    fn fd_d1(&self, u: f64) -> GpVec {
        let h = 1e-6;
        let p1 = self.d0(u + h);
        let p2 = self.d0(u - h);
        GpVec::new(
            (p1.x() - p2.x()) / (2.0 * h),
            (p1.y() - p2.y()) / (2.0 * h),
            (p1.z() - p2.z()) / (2.0 * h),
        )
    }

    /// Finite-difference second derivative: (f(u+h) - 2f(u) + f(u-h)) / h^2.
    fn fd_d2(&self, u: f64) -> GpVec {
        let h = 1e-6;
        let p1 = self.d0(u + h);
        let p0 = self.d0(u);
        let p2 = self.d0(u - h);
        GpVec::new(
            (p1.x() - 2.0 * p0.x() + p2.x()) / (h * h),
            (p1.y() - 2.0 * p0.y() + p2.y()) / (h * h),
            (p1.z() - 2.0 * p0.z() + p2.z()) / (h * h),
        )
    }

    /// Private de Boor evaluation. The public d0 delegates to `eval`; this is
    /// kept as a self-contained reference.
    #[allow(dead_code)]
    fn de_boor(&self, u: f64) -> GpPnt {
        let n = self.poles.len();
        if n == 0 { return GpPnt::zero(); }
        let idx = knots::hunt(&self.knots, u).max(self.degree).min(n - 1);
        let mut pts = vec![GpPnt::zero(); self.degree + 1];
        for k in 0..=self.degree {
            let pi = idx - self.degree + k;
            pts[k] = if pi < n { self.poles[pi] } else { self.poles[n - 1] };
        }
        for r in 1..=self.degree {
            for i in (r..=self.degree).rev() {
                let k0 = idx + i - self.degree;
                let k1 = k0 + self.degree + 1 - r;
                let alpha = (u - self.knots[k0]) / (self.knots[k1] - self.knots[k0]);
                if alpha.is_finite() {
                    pts[i] = GpPnt::new(
                        (1.0 - alpha) * pts[i - 1].x() + alpha * pts[i].x(),
                        (1.0 - alpha) * pts[i - 1].y() + alpha * pts[i].y(),
                        (1.0 - alpha) * pts[i - 1].z() + alpha * pts[i].z(),
                    );
                }
            }
        }
        pts[self.degree]
    }
}

impl Curve for GeomBSplineCurve {
    fn d0(&self, u: f64) -> GpPnt {
        // `Geom_BSplineCurve::EvalD0` (`Geom_BSplineCurve_1.cxx:174-194`) ->
        // `BSplCLib::D0` over `myKnots` + `&myMults` (distinct knots).
        let (new_u, span, uknots, umults) = self.eval_span(u);
        occt_core::bspl::curve_eval::curve_d0(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        )
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        // `Geom_BSplineCurve::EvalD1` (`:198-227`) -> `BSplCLib::D1`.
        let (new_u, span, uknots, umults) = self.eval_span(u);
        occt_core::bspl::curve_eval::curve_d1(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        )
    }

    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        // `Geom_BSplineCurve::EvalDN` (`Geom_BSplineCurve_1.cxx:300-316`).
        // Illegal N<1 returns zero instead of throw. Eval-rep is empty.
        if n < 1 {
            return GpVec::zero();
        }
        occt_core::bspl::curve_dn::curve_dn(
            u,
            n,
            0,
            self.degree as i32,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.knots,
            None,
        )
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        // `Geom_BSplineCurve::EvalD2` (`:231-261`) -> `BSplCLib::D2`.
        let (new_u, span, uknots, umults) = self.eval_span(u);
        let r = occt_core::bspl::curve_eval::curve_d2(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        );
        r
    }

    /// `Geom_BSplineCurve::EvalD3` (`:265-296`) -> `BSplCLib::D3`, which
    /// `Geom_OffsetCurve` calls through `basisCurve->D3(U, …)` in
    /// `CalculateD2`. Without this override the trait default returned a
    /// **zero** third derivative for every B-spline basis.
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (new_u, span, uknots, umults) = self.eval_span(u);
        occt_core::bspl::curve_eval::curve_d3(
            new_u, span, self.degree as i32, self.periodic, &self.poles,
            self.weights.as_deref(), uknots, Some(umults),
        )
    }

    fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }
    fn is_periodic(&self) -> bool { self.periodic }
    fn continuity(&self) -> u8 {
        occt_core::bspl::local_continuity(
            &self.knots,
            self.degree,
            self.periodic,
            self.first_parameter(),
            self.last_parameter(),
        )
    }

    fn transform(&mut self, t: &GpTrsf) {
        for p in self.poles.iter_mut() {
            *p = p.transformed(t);
        }
    }

    fn reverse(&mut self) {
        // `Geom_BSplineCurve::Reverse` (`cxx:496-511`): `BSplCLib::Reverse` on
        // the distinct knots and on the multiplicities, then the poles
        // (periodic window `myFlatKnots.Upper() - Deg - 1`, `cxx:502-510`),
        // then `updateKnots()`.
        let (mut uknots, mut umults) = self.distinct_knots_and_mults();
        curve_tools::reverse_distinct_knots(&mut uknots);
        umults.reverse();
        if self.periodic {
            curve_tools::reverse_periodic_span(&mut self.poles, self.knots.len(), self.degree);
            if let Some(w) = self.weights.as_mut() {
                curve_tools::reverse_periodic_span(w, self.knots.len(), self.degree);
            }
            self.knots = knots::knot_sequence_periodic(&uknots, &umults, self.degree as i32);
        } else {
            self.poles.reverse();
            if let Some(w) = self.weights.as_mut() {
                w.reverse();
            }
            self.knots =
                occt_core::bspl::banded_interp::knot_sequence(&uknots, &umults, self.degree as i32);
        }
        self.refresh_knot_rep();
    }

    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
    fn bspline_poles(&self) -> Option<&[GpPnt]> { Some(&self.poles) }
    fn bspline_weights(&self) -> Option<&[f64]> { self.weights.as_deref() }
    fn bspline_knots(&self) -> Option<&[f64]> { Some(&self.knots) }
    fn nurbs_degree(&self) -> Option<usize> { Some(self.degree) }
    fn resolution(&self, r3d: f64) -> f64 {
        occt_core::bspl::bspline_curve_resolution(
            &self.poles,
            self.weights.as_deref(),
            &self.knots,
            self.degree as i32,
            r3d,
        )
    }
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        let eps = self
            .resolution(occt_core::precision::CONFUSION)
            .min(occt_core::precision::PCONFUSION);
        occt_core::bspl::adaptor_intervals(
            &self.knots,
            self.degree,
            self.periodic,
            continuity,
            self.first_parameter(),
            self.last_parameter(),
            eps,
        )
    }
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity).len().saturating_sub(1).max(1) as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_midpoint() {
        let c = GeomBSplineCurve::new(
            vec![GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.)],
            vec![0., 0., 1., 1.],
            1,
        ).unwrap();
        let p = c.d0(0.5);
        assert!((p.x() - 0.5).abs() < 1e-12);
        assert!((p.y() - 0.0).abs() < 1e-12);
        assert!((p.z() - 0.0).abs() < 1e-12);
    }
}
