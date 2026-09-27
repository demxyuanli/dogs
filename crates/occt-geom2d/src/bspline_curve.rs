//! 2D B-spline curve. Source: `Geom2d_BSplineCurve.hxx`

use crate::curve::Curve2d;
use occt_core::bspl::{eval, knots};
use occt_core::gp::{GpPnt, GpPnt2d, GpTrsf2d, GpVec2d};

/// 2D B-spline curve (non-rational), stored as separate x/y pole arrays.
#[derive(Clone)]
pub struct Geom2dBSplineCurve {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub knots: Vec<f64>,
    pub degree: usize,
    /// `Geom2d_BSplineCurve::IsPeriodic()`: with `true` the flat knot vector is
    /// the periodic `BSplCLib::KnotSequence` (extended by one period each side).
    pub periodic: bool,
}

impl Geom2dBSplineCurve {
    /// Build a 2D B-spline. xs/ys lengths must match and knot count must be
    /// poles + degree + 1.
    pub fn new(xs: Vec<f64>, ys: Vec<f64>, knots: Vec<f64>, degree: usize) -> Result<Self, &'static str> {
        if xs.len() != ys.len() {
            return Err("Geom2dBSplineCurve: xs/ys length mismatch");
        }
        knots::check_degree(xs.len(), degree, knots.len())?;
        Ok(Self { xs, ys, knots, degree, periodic: false })
    }

    /// `Geom2d_BSplineCurve::Knots()` / `Multiplicities()`: stored distinct knots
    /// and multiplicities (see `GeomBSplineCurve::distinct_knots_and_mults` in
    /// `occt-geom` for why the periodic period-extension knots are excluded).
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
        } else if nbp > self.xs.len() {
            // OCCT's `Resize` leaves the new poles default-constructed (`gp_Pnt2d()`).
            self.xs.resize(nbp, 0.0);
            self.ys.resize(nbp, 0.0);
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
        self.knots = occt_core::bspl::unperiodize::flat_knots_from_mults(&new_knots, &new_mults);
        self.periodic = false;
    }

    pub fn nb_poles(&self) -> usize { self.xs.len() }
    pub fn degree(&self) -> usize { self.degree }
    pub fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    pub fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }

    /// Lift `(x, y)` poles to `z = 0` for `BSplCLib` evaluators.
    fn poles_3d(&self) -> Vec<GpPnt> {
        self.xs
            .iter()
            .zip(self.ys.iter())
            .map(|(&x, &y)| GpPnt::new(x, y, 0.0))
            .collect()
    }

    /// De Boor triangular evaluation on (x, y) pole pairs.
    fn de_boor(&self, u: f64) -> GpPnt2d {
        let n = self.xs.len();
        if n == 0 { return GpPnt2d::zero(); }
        let idx = knots::hunt(&self.knots, u).max(self.degree).min(n - 1);
        let mut x = vec![0.0f64; self.degree + 1];
        let mut y = vec![0.0f64; self.degree + 1];
        for k in 0..=self.degree {
            let pi = (idx - self.degree + k).min(n - 1);
            x[k] = self.xs[pi];
            y[k] = self.ys[pi];
        }
        for r in 1..=self.degree {
            for i in (r..=self.degree).rev() {
                let k0 = idx + i - self.degree;
                let k1 = k0 + self.degree + 1 - r;
                let alpha = (u - self.knots[k0]) / (self.knots[k1] - self.knots[k0]);
                if alpha.is_finite() {
                    x[i] = (1.0 - alpha) * x[i - 1] + alpha * x[i];
                    y[i] = (1.0 - alpha) * y[i - 1] + alpha * y[i];
                }
            }
        }
        GpPnt2d::new(x[self.degree], y[self.degree])
    }
}

impl Curve2d for Geom2dBSplineCurve {
    fn d0(&self, u: f64) -> GpPnt2d {
        if self.periodic {
            // `Geom2d_BSplineCurve::D0` on a periodic flat knot sequence:
            // `BSplCLib::PrepareEval` wraps the pole window and `LocateParameter`
            // maps the parameter into the period; `curve_dn::dn` with order 0 is
            // `BSplCLib::D0`.
            let v = occt_core::bspl::curve_dn::dn(
                u, 0, 0, self.degree as i32, true, &self.poles_3d(), None, &self.knots, None,
            );
            return GpPnt2d::new(v.x(), v.y());
        }
        self.de_boor(u)
    }

    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        // `Geom2d_BSplineCurve::D1` / `BSplCLib::D1`.
        let poles = self.poles_3d();
        if self.periodic {
            let d = occt_core::bspl::curve_dn::dn(
                u, 1, 0, self.degree as i32, true, &poles, None, &self.knots, None,
            );
            return (self.d0(u), GpVec2d::new(d.x(), d.y()));
        }
        let (p, d) = eval::eval_curve_d1(&poles, &self.knots, self.degree, u);
        (GpPnt2d::new(p.x(), p.y()), GpVec2d::new(d.x(), d.y()))
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        // `Geom2d_BSplineCurve::D2` / `BSplCLib::D2`.
        let poles = self.poles_3d();
        if self.periodic {
            let d1 = occt_core::bspl::curve_dn::dn(
                u, 1, 0, self.degree as i32, true, &poles, None, &self.knots, None,
            );
            let d2 = occt_core::bspl::curve_dn::dn(
                u, 2, 0, self.degree as i32, true, &poles, None, &self.knots, None,
            );
            return (
                self.d0(u),
                GpVec2d::new(d1.x(), d1.y()),
                GpVec2d::new(d2.x(), d2.y()),
            );
        }
        let (p, d1, d2) = eval::eval_curve_d2(&poles, &self.knots, self.degree, u);
        (
            GpPnt2d::new(p.x(), p.y()),
            GpVec2d::new(d1.x(), d1.y()),
            GpVec2d::new(d2.x(), d2.y()),
        )
    }

    fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }
    fn is_periodic(&self) -> bool { self.periodic }
    fn period(&self) -> f64 {
        if self.periodic {
            self.last_parameter() - self.first_parameter()
        } else {
            0.0
        }
    }
    fn continuity(&self) -> u8 {
        // `Geom2d_BSplineCurve::Continuity` / `GeomAdaptor` LocalContinuity.
        occt_core::bspl::local_continuity(
            &self.knots,
            self.degree,
            self.periodic,
            self.first_parameter(),
            self.last_parameter(),
        )
    }
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        occt_core::bspl::adaptor_intervals(
            &self.knots,
            self.degree,
            self.periodic,
            continuity,
            self.first_parameter(),
            self.last_parameter(),
            occt_core::precision::PCONFUSION,
        )
    }
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }

    fn transform(&mut self, t: &GpTrsf2d) {
        for i in 0..self.xs.len() {
            let mut p = GpPnt2d::new(self.xs[i], self.ys[i]);
            p.transform(t);
            self.xs[i] = p.x();
            self.ys[i] = p.y();
        }
    }

    fn reverse(&mut self) {
        // `Geom2d_BSplineCurve::Reverse` (`Geom2d_BSplineCurve.cxx:677-696`) runs
        // `BSplCLib::Reverse(myKnots)` + `BSplCLib::Reverse(myMults)` + reverse
        // the poles (+ weights) + `updateKnots()`. `BSplCLib::Reverse(
        // NCollection_Array1<double>& Knots)` (`BSplCLib.cxx:802-828`) maps every
        // knot to `kfirst + klast - k`, so the reversed curve keeps the *same*
        // parameter range; reversing the flat knot array already reverses the
        // multiplicities, so applying that affine map after the swap is
        // equivalent. `klast - k` alone (the former code) shifts the range by
        // `-kfirst` for any curve whose first knot is not 0 — e.g. the spherical
        // pcurves of `data/Offset.step`, whose knots start at `pi/2`:
        // `build_arc` then evaluated the reversed curve over a parameter window
        // the curve does not cover.
        self.xs.reverse();
        self.ys.reverse();
        let n = self.knots.len();
        if n == 0 {
            return;
        }
        let (kfirst, klast) = (self.knots[0], self.knots[n - 1]);
        for i in 0..n / 2 {
            self.knots.swap(i, n - 1 - i);
        }
        for k in self.knots.iter_mut() {
            *k = kfirst + klast - *k;
        }
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> { Box::new(self.clone()) }
    fn is_bspline2d(&self) -> bool {
        true
    }
    fn poles2d(&self) -> Option<Vec<GpPnt2d>> {
        Some(self.xs.iter().zip(self.ys.iter()).map(|(x, y)| GpPnt2d::new(*x, *y)).collect())
    }
    fn set_poles2d(&mut self, poles: &[GpPnt2d]) {
        for (i, p) in poles.iter().enumerate() {
            self.xs[i] = p.x();
            self.ys[i] = p.y();
        }
    }

    /// `Geom2d_BSplineCurve::NbKnots()` (`Geom2d_BSplineCurve_1.cxx:598-601`):
    /// the number of distinct knots. The curve stores the expanded knot
    /// sequence, so compress it as `Geom_BSplineCurve::Knots` does.
    fn bspline_nb_knots(&self) -> Option<usize> {
        Some(knots::unique_knots_mults(&self.knots).0.len())
    }

    /// `Geom2d_BSplineCurve::Knots()` flat knot sequence (`Knot(j)` in
    /// `ShapeAnalysis_TransferParametersProj::CorrectParameter`,
    /// `Proj.cxx:268-279`).
    fn bspline_knots2d(&self) -> Option<&[f64]> {
        Some(&self.knots)
    }

    /// `Geom2d_BSplineCurve::Degree()` (`Geom2d_BSplineCurve_1.cxx:168-171`).
    fn bspline_degree(&self) -> Option<usize> {
        Some(self.degree)
    }

    fn bspline_poles2d(&self) -> Option<(&[f64], &[f64])> {
        Some((&self.xs, &self.ys))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_midpoint() {
        let c = Geom2dBSplineCurve::new(
            vec![0.0, 1.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0, 1.0, 1.0],
            1,
        ).unwrap();
        let p = c.d0(0.5);
        assert!((p.x() - 0.5).abs() < 1e-12);
        assert!((p.y() - 0.0).abs() < 1e-12);
    }
}
