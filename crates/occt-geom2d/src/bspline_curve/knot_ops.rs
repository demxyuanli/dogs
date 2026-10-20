use crate::curve::Curve2d;
use occt_core::bspl::knots;

use super::curve::Geom2dBSplineCurve;
use super::MAX_DEGREE;

impl Geom2dBSplineCurve {
    /// `Geom2d_BSplineCurve::InsertKnots(Knots, Mults, ParametricTolerance, Add)`
    /// (`Geom2d_BSplineCurve.cxx:343-393`).
    ///
    /// `add_knots` are **distinct** parameters, `add_mults[i]` their
    /// multiplicities (`None` = flat insertion, one per knot). With
    /// `add = true` an existing knot's multiplicity grows by `M`; with
    /// `add = false` it is raised to `M` (the header default,
    /// `Geom2d_BSplineCurve.hxx:302-304`). The tolerance for knot equality is
    /// `max(Epsilon(U), parametric_tolerance)`.
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
            None => return Err("Geom2d_BSplineCurve::InsertKnots"),
        };
        if nbpoles as usize == self.xs.len() {
            return Ok(());
        }
        // The port stores (x, y) separately; `BSplCLib::InsertKnots` is
        // dimension-agnostic, so lift the poles to `z = 0` exactly as the
        // evaluators do (weights are carried alongside for a rational curve).
        let poles3d = self.poles_3d();
        let out = match occt_core::bspl::insert_knots::insert_knots(
            self.degree as i32,
            self.periodic,
            &poles3d,
            self.weights.as_deref(),
            &uknots,
            &umults,
            add_knots,
            add_mults,
            parametric_tolerance,
            add,
        ) {
            Some(out) => out,
            None => return Err("Geom2d_BSplineCurve::InsertKnots"),
        };
        // `BSplCLib::InsertKnots` output: refit poles and, for a rational curve,
        // the matching weights.
        self.xs = out.poles.iter().map(|p| p.x()).collect();
        self.ys = out.poles.iter().map(|p| p.y()).collect();
        if self.weights.is_some() {
            self.weights = out.weights;
        }
        self.update_flat_knots(&out.knots, &out.mults);
        Ok(())
    }

    /// `Geom2d_BSplineCurve::InsertKnot(U, M, ParametricTolerance)`
    /// (`Geom2d_BSplineCurve.cxx:332-341`), defaults `M = 1`,
    /// `ParametricTolerance = 0.0` (`Geom2d_BSplineCurve.hxx:280-282`).
    pub fn insert_knot(
        &mut self,
        u: f64,
        m: i32,
        parametric_tolerance: f64,
    ) -> Result<(), &'static str> {
        self.insert_knots(&[u], Some(&[m]), parametric_tolerance, true)
    }

    /// `Geom2d_BSplineCurve::IncreaseDegree(Degree)`
    /// (`Geom2d_BSplineCurve.cxx:236-292`): raise the degree to `degree`
    /// (`degree < self.degree` or `degree > MaxDegree()` is OCCT's
    /// `Standard_ConstructionError`, `:243-246`).
    pub fn increase_degree(&mut self, degree: usize) -> Result<(), &'static str> {
        if degree == self.degree {
            return Ok(());
        }
        if degree < self.degree || degree > MAX_DEGREE {
            return Err("Geom2dBSplineCurve::increase_degree: bad degree value");
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let (from_k1, to_k2) = self.u_knot_index_range(&umults);
        let step = degree as i32 - self.degree as i32;
        // `cxx:251`: the new pole count `myPoles.Length() + Step * (ToK2 - FromK1)`.
        let nb_new_poles = self.xs.len() as i32 + step * (to_k2 - from_k1);
        // `cxx:255`: `BSplCLib::IncreaseDegreeCountKnots`.
        let nb_new_knots = occt_core::bspl::increase_degree::increase_degree_count_knots(
            self.degree as i32,
            degree as i32,
            self.periodic,
            &umults,
        );
        if nb_new_poles <= 0 || nb_new_knots <= 0 {
            return Err("Geom2dBSplineCurve::increase_degree: bad degree value");
        }
        // `BSplCLib::IncreaseDegree` is dimension-agnostic; lift the (x, y)
        // poles to `z = 0` exactly as the evaluators do.
        let poles3d = self.poles_3d();
        let out = occt_core::bspl::increase_degree::increase_degree(
            self.degree as i32,
            degree as i32,
            self.periodic,
            &poles3d,
            self.weights.as_deref(),
            &uknots,
            &umults,
            nb_new_poles as usize,
            nb_new_knots as usize,
        );
        self.xs = out.poles.iter().map(|p| p.x()).collect();
        self.ys = out.poles.iter().map(|p| p.y()).collect();
        self.weights = out.weights;
        self.degree = degree;
        self.update_flat_knots(&out.knots, &out.mults);
        Ok(())
    }

    /// `Geom2d_BSplineCurve::RemoveKnot(Index, M, Tolerance)`
    /// (`Geom2d_BSplineCurve.cxx:410-478`): lower the multiplicity of the
    /// distinct knot `index` to `m` (`m == 0` removes the knot entirely).
    ///
    /// Returns `Ok(false)` exactly where `BSplCLib::RemoveKnot` reports failure;
    /// the out-of-range index is OCCT's `Standard_OutOfRange` (`:420-423`).
    pub fn remove_knot(&mut self, index: i32, m: i32, tolerance: f64) -> Result<bool, &'static str> {
        if m < 0 {
            return Ok(true);
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let (i1, i2) = self.u_knot_index_range(&umults);
        if index < i1 || index > i2 {
            return Err("Geom2dBSplineCurve::remove_knot: index out of range");
        }
        let idx = (index - 1) as usize;
        if idx >= umults.len() {
            return Err("Geom2dBSplineCurve::remove_knot: index out of range");
        }
        let step = umults[idx] - m;
        if step <= 0 {
            return Ok(true);
        }
        let nb_new_poles = self.xs.len() as i32 - step;
        let nb_new_knots = uknots.len() as i32 - if m == 0 { 1 } else { 0 };
        if nb_new_poles < 0 || nb_new_knots < 0 {
            return Err("Geom2dBSplineCurve::remove_knot: index out of range");
        }
        let poles3d = self.poles_3d();
        let out = occt_core::bspl::remove_knot::remove_knot(
            index,
            m,
            self.degree as i32,
            self.periodic,
            &poles3d,
            self.weights.as_deref(),
            &uknots,
            &umults,
            nb_new_poles as usize,
            nb_new_knots as usize,
            tolerance,
        );
        match out {
            Some(out) => {
                self.xs = out.poles.iter().map(|p| p.x()).collect();
                self.ys = out.poles.iter().map(|p| p.y()).collect();
                self.weights = out.weights;
                self.update_flat_knots(&out.knots, &out.mults);
                Ok(true)
            }
            None => Ok(false),
        }
    }

    /// `Geom2d_BSplineCurve::Segment(U1, U2, theTolerance)`
    /// (`Geom2d_BSplineCurve.cxx:707-889`): restrict the curve to `[U1, U2]`
    /// (`U2 < U1` and, for a periodic curve, `(U2 - U1) - Period >
    /// Precision::PConfusion()` are `Standard_DomainError`).
    ///
    /// Rust has no default arguments: callers reproduce OCCT's
    /// `theTolerance = Precision::PConfusion()` default explicitly.
    pub fn segment(&mut self, u1: f64, u2: f64, tolerance: f64) -> Result<(), &'static str> {
        if u2 < u1 {
            return Err("Geom2d_BSplineCurve::Segment");
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
                return Err("Geom2d_BSplineCurve::Segment");
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
            return Err("Geom2d_BSplineCurve::Segment");
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
        // (`Geom2d_BSplineCurve.hxx:302-304`).
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
            return Err("Geom2d_BSplineCurve::Segment");
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
            return Err("Geom2d_BSplineCurve::Segment");
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

        // compute index1 and index2 to set the new poles
        let mut pindex1 = knots::pole_index(self.degree as i32, index1, self.periodic, &um2);
        let mut pindex2 = knots::pole_index(self.degree as i32, index2, self.periodic, &um2);
        pindex1 += 1;
        pindex2 = (pindex2 + 1).min(self.xs.len() as i32);
        if pindex1 < 1 || pindex2 < pindex1 {
            return Err("Geom2d_BSplineCurve::Segment");
        }

        let mut newx = Vec::with_capacity((pindex2 - pindex1 + 1) as usize);
        let mut newy = Vec::with_capacity((pindex2 - pindex1 + 1) as usize);
        for i in pindex1..=pindex2 {
            newx.push(self.xs[(i - 1) as usize]);
            newy.push(self.ys[(i - 1) as usize]);
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

        self.xs = newx;
        self.ys = newy;
        if self.weights.is_some() {
            self.weights = newweights.take();
        }
        self.update_flat_knots(&nknots, &nmults);
        Ok(())
    }
}
