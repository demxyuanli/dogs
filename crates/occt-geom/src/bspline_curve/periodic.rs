use crate::curve::Curve;
use occt_core::bspl::knots;
use occt_core::gp::GpPnt;

use super::curve::GeomBSplineCurve;

impl GeomBSplineCurve {
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
}
