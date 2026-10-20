use crate::curve::Curve;
use occt_core::bspl::{knots, poles};

use super::curve::GeomBSplineCurve;
use super::MAX_DEGREE;

impl GeomBSplineCurve {
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
}
