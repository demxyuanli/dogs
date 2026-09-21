//! 3D B-spline curve. Source: `Geom_BSplineCurve.hxx`

use crate::curve::Curve;
use occt_core::gp::{GpPnt, GpVec, GpTrsf};
use occt_core::bspl::{knots, eval, poles, curve_tools};

/// Non-rational or rational B-spline curve in 3D.
#[derive(Clone)]
pub struct GeomBSplineCurve {
    pub poles: Vec<GpPnt>,
    pub weights: Option<Vec<f64>>,
    pub knots: Vec<f64>,
    pub degree: usize,
    pub periodic: bool,
}

impl GeomBSplineCurve {
    /// Build a non-rational B-spline. Knot count must be poles + degree + 1.
    pub fn new(poles: Vec<GpPnt>, knots: Vec<f64>, degree: usize) -> Result<Self, &'static str> {
        knots::check_degree(poles.len(), degree, knots.len())?;
        Ok(Self { poles, weights: None, knots, degree, periodic: false })
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
        Ok(Self { poles, weights: Some(weights), knots, degree, periodic: false })
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

    /// `Geom_BSplineCurve::Knots()` / `Multiplicities()`
    /// (`Geom_BSplineCurve_1.cxx:380-383`, `:148-151`): the stored distinct
    /// knots and their multiplicities.
    ///
    /// The port stores the flat sequence only. For a periodic curve the
    /// period-extension knots live **outside** `[FirstParameter,
    /// LastParameter]` (they are `stored_knot ± period`), so the run lengths of
    /// the flat array restricted to that interval are exactly `myMults` — the
    /// equivalent of OCCT's separately stored arrays.
    pub fn distinct_knots_and_mults(&self) -> (Vec<f64>, Vec<i32>) {
        let (uknots, umults) = knots::unique_knots_mults(&self.knots);
        if !self.periodic || uknots.is_empty() {
            return (uknots, umults);
        }
        let (first, last) = (self.first_parameter(), self.last_parameter());
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
    pub fn segment(&mut self, u1: f64, u2: f64, the_tolerance: f64) -> Result<(), &'static str> {
        if u2 < u1 {
            return Err("Geom_BSplineCurve::Segment");
        }

        let was_periodic = self.periodic;
        let new_u1;
        let mut new_u2;
        let mut du = 0.0;
        let mut a_ddu = 0.0;

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
            a_ddu = du;
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

        let a_nu2 = new_u2;

        let knots_pair = [new_u1.min(new_u2), new_u1.max(new_u2)];
        let mults_pair = [self.degree as i32, self.degree as i32];

        let abs_u_max = new_u1
            .abs()
            .max(new_u2.abs())
            .max(self.first_parameter().abs())
            .max(self.last_parameter().abs());
        let eps = occt_core::precision::epsilon(abs_u_max).max(the_tolerance);

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
            if a_nu2 < u2 {
                nknots[last_k] = u1 + a_ddu;
            }
        }

        self.poles = newpoles;
        if self.weights.is_some() {
            self.weights = newweights.take();
        }        self.update_flat_knots(&nknots, &nmults);
        Ok(())
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
        if self.periodic {
            // Periodic flat knot sequence: `BSplCLib::D0` goes through
            // `PrepareEval` (`LocateParameter` maps the parameter into the
            // period, `BuildEval` wraps the pole window) and `Bohm(…, 0, …)`.
            // `curve_dn::dn` is exactly that machinery with an explicit
            // derivative order; order 0 returns the point.
            let v = occt_core::bspl::curve_dn::dn(
                u, 0, 0, self.degree as i32, true, &self.poles,
                self.weights.as_deref(), &self.knots, None,
            );
            return GpPnt::new(v.x(), v.y(), v.z());
        }
        match &self.weights {
            Some(w) => eval::eval_curve_rational(&self.poles, w, &self.knots, self.degree, u),
            None => eval::eval_curve(&self.poles, &self.knots, self.degree, u),
        }
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        if self.periodic {
            return (self.d0(u), self.eval_dn(u, 1));
        }
        // `Geom_BSplineCurve::D1` / `BSplCLib::D1`. Rational uses the
        // homogeneous quotient already computed by `eval_curve_rational_d2`.
        match &self.weights {
            Some(w) => {
                let (p, d1, _) =
                    eval::eval_curve_rational_d2(&self.poles, w, &self.knots, self.degree, u);
                (p, d1)
            }
            None => eval::eval_curve_d1(&self.poles, &self.knots, self.degree, u),
        }
    }

    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        // `Geom_BSplineCurve::EvalDN` (`Geom_BSplineCurve_1.cxx:300-316`).
        // Illegal N<1 returns zero instead of throw. Eval-rep is empty.
        if n < 1 {
            return GpVec::zero();
        }
        occt_core::bspl::curve_dn::dn(
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
        if self.periodic {
            return (self.d0(u), self.eval_dn(u, 1), self.eval_dn(u, 2));
        }
        match &self.weights {
            Some(w) => eval::eval_curve_rational_d2(&self.poles, w, &self.knots, self.degree, u),
            None => eval::eval_curve_d2(&self.poles, &self.knots, self.degree, u),
        }
    }

    /// `Geom_BSplineCurve::D3` (`BSplCLib::DN(..., 3)`), which `Geom_OffsetCurve`
    /// calls through `basisCurve->D3(U, …)` in `CalculateD2`. Without this
    /// override the trait default returned a **zero** third derivative for every
    /// B-spline basis.
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (p, d1, d2) = self.d2(u);
        (p, d1, d2, self.eval_dn(u, 3))
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
        curve_tools::reverse_curve(&mut self.poles, &mut self.knots);
        if let Some(w) = self.weights.as_mut() {
            w.reverse();
        }
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
