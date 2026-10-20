use occt_core::bspl::knots;

use super::curve::Geom2dBSplineCurve;

impl Geom2dBSplineCurve {
    /// `Geom2d_BSplineCurve::SetPeriodic()` (`Geom2d_BSplineCurve.cxx:948-…`):
    /// same construction as `Geom_BSplineCurve::SetPeriodic`
    /// (`Geom_BSplineCurve.cxx:777-815`).
    pub fn set_periodic(&mut self) {
        let (uknots, umults) = self.distinct_knots_and_mults();
        if uknots.is_empty() || umults.is_empty() {
            return;
        }
        let degree = self.degree as i32;
        let (first, last) = if self.periodic {
            (1usize, uknots.len())
        } else {
            (
                occt_core::bspl::locate::first_u_knot_index(degree, &umults).max(1) as usize,
                occt_core::bspl::locate::last_u_knot_index(degree, &umults).max(1) as usize,
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
        if nbp < self.xs.len() {
            self.xs.truncate(nbp);
            self.ys.truncate(nbp);
            if let Some(w) = self.weights.as_mut() {
                w.truncate(nbp);
            }
        } else if nbp > self.xs.len() {
            // OCCT's `Resize` leaves the new poles default-constructed (`gp_Pnt2d()`).
            self.xs.resize(nbp, 0.0);
            self.ys.resize(nbp, 0.0);
            if let Some(w) = self.weights.as_mut() {
                w.resize(nbp, 0.0);
            }
        }
        self.knots = knots::knot_sequence_periodic(&uknots, &umults, degree);
        self.periodic = true;
    }

    /// `Geom2d_BSplineCurve::SetNotPeriodic()` (`Geom2d_BSplineCurve.cxx:1087-…`):
    /// `BSplCLib::PrepareUnperiodize` + `BSplCLib::Unperiodize`
    /// (`BSplCLib.cxx:2967-3080`), poles re-indexed as
    /// `NewPoles(k) = Poles((k - 1) % n_old + 1)`.
    pub fn set_not_periodic(&mut self) {
        if !self.periodic {
            return;
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let degree = self.degree as i32;
        let (new_knots, new_mults, _index) =
            occt_core::bspl::unperiodize::unperiodize_knots(degree, &uknots, &umults);
        let n_new = (new_mults.iter().sum::<i32>() - degree - 1).max(0) as usize;
        let n_old = self.xs.len();
        if n_old == 0 || n_new == 0 {
            return;
        }
        self.xs = (0..n_new).map(|k| self.xs[k % n_old]).collect();
        self.ys = (0..n_new).map(|k| self.ys[k % n_old]).collect();
        if let Some(w) = self.weights.as_ref() {
            let old = w.clone();
            self.weights = Some((0..n_new).map(|k| old[k % n_old]).collect());
        }
        self.knots = occt_core::bspl::unperiodize::flat_knots_from_mults(&new_knots, &new_mults);
        self.periodic = false;
    }

    /// `Geom2d_BSplineCurve::SetOrigin(Index)` (`Geom2d_BSplineCurve.cxx:989-1084`):
    /// rotate a periodic curve so that the knot `Index` becomes the origin
    /// (the knots after it stay, the ones before it move one period up and the
    /// poles are rotated to match).
    pub fn set_origin(&mut self, index: i32) -> Result<(), &'static str> {
        if !self.periodic {
            return Err("Geom2d_BSplineCurve::SetOrigin");
        }
        let (uknots, umults) = self.distinct_knots_and_mults();
        let (first, last) = self.u_knot_index_range(&umults);
        if index < first || index > last {
            return Err("Geom2d_BSplineCurve::SetOrigin");
        }
        let nbknots = uknots.len() as i32;
        let nbpoles = self.xs.len() as i32;
        if nbknots == 0 || nbpoles == 0 {
            return Err("Geom2d_BSplineCurve::SetOrigin");
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

        // set the poles
        let mut newx = Vec::with_capacity(nbpoles as usize);
        let mut newy = Vec::with_capacity(nbpoles as usize);
        for i in pole_first..=nbpoles {
            newx.push(self.xs[(i - 1) as usize]);
            newy.push(self.ys[(i - 1) as usize]);
        }
        for i in 1..pole_first {
            newx.push(self.xs[(i - 1) as usize]);
            newy.push(self.ys[(i - 1) as usize]);
        }

        self.xs = newx;
        self.ys = newy;
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
        self.update_flat_knots(&newknots, &newmults);
        Ok(())
    }
}
