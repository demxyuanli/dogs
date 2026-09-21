//! 3D B-spline curve. Source: `Geom_BSplineCurve.hxx`

use crate::curve::Curve;
use occt_core::gp::{GpPnt, GpVec, GpTrsf};
use occt_core::bspl::{knots, eval, poles, bezier, curve_tools};

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
        if self.periodic {
            // OCCT's `SetPeriodic()` on an already periodic curve is a no-op:
            // `FirstUKnotIndex()`/`LastUKnotIndex()` return `1`/`myKnots.Length()`
            // (`Geom_BSplineCurve_1.cxx:334-344`, `:404-414`), the end
            // multiplicities are already `<= degree`, `NbPoles(degree, true, …)`
            // reproduces the pole count and `updateKnots()` rebuilds the same
            // flat sequence. The port stores only the flat sequence, whose
            // run-length decomposition would additionally expose the
            // period-extension knots (OCCT keeps those in `myFlatKnots` only),
            // so the no-op is taken explicitly rather than recomputed.
            return;
        }
        let (uknots, umults) = knots::unique_knots_mults(&self.knots);
        if uknots.is_empty() || umults.is_empty() {
            return;
        }
        let degree = self.degree as i32;
        let first = occt_core::bspl::locate::first_u_knot_index(degree, &umults).max(1) as usize;
        let last = occt_core::bspl::locate::last_u_knot_index(degree, &umults).max(1) as usize;
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

    /// Insert knot `u` with multiplicity `mult` (Boehm knot insertion).
    pub fn insert_knot(&mut self, u: f64, mult: usize) {
        for _ in 0..mult {
            let idx = knots::hunt(&self.knots, u);
            bezier::boehm_insert(&mut self.poles, &self.knots, idx, u, self.degree, self.weights.as_mut());
            self.knots = knots::insert_knot(&self.knots, u, 1);
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
