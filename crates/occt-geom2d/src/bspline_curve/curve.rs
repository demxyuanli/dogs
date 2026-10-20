use crate::curve::Curve2d;
use occt_core::bspl::knots;
use occt_core::gp::GpPnt;

/// 2D B-spline curve (non-rational), stored as separate x/y pole arrays.
#[derive(Clone)]
pub struct Geom2dBSplineCurve {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub weights: Option<Vec<f64>>,
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
        Ok(Self { xs, ys, weights: None, knots, degree, periodic: false })
    }

    /// Build a rational 2D B-spline (`Geom2d_BSplineCurve(Poles, Weights, Knots,
    /// Mults, Degree)`); `weights` length must equal the pole count.
    pub fn rational(
        xs: Vec<f64>,
        ys: Vec<f64>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        degree: usize,
    ) -> Result<Self, &'static str> {
        if xs.len() != ys.len() {
            return Err("Geom2dBSplineCurve: xs/ys length mismatch");
        }
        if weights.len() != xs.len() {
            return Err("Geom2dBSplineCurve: weights/poles length mismatch");
        }
        knots::check_degree(xs.len(), degree, knots.len())?;
        Ok(Self { xs, ys, weights: Some(weights), knots, degree, periodic: false })
    }

    /// Internal constructor mirroring OCCT's
    /// `Geom2d_BSplineCurve(Poles, Weights, Knots, Mults, Degree, IsPeriodic)`
    /// constructor after `updateKnots()`: the caller passes the already-expanded
    /// flat knot sequence (the period extension included when `periodic`), so
    /// the `n_poles + degree + 1` count check only applies to the non-periodic
    /// representation.
    pub fn from_flat(
        xs: Vec<f64>,
        ys: Vec<f64>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        degree: usize,
        periodic: bool,
    ) -> Result<Self, &'static str> {
        if xs.len() != ys.len() {
            return Err("Geom2dBSplineCurve: xs/ys length mismatch");
        }
        if let Some(w) = &weights {
            if w.len() != xs.len() {
                return Err("Geom2dBSplineCurve: weights/poles length mismatch");
            }
        }
        if !periodic {
            knots::check_degree(xs.len(), degree, knots.len())?;
        }
        Ok(Self { xs, ys, weights, knots, degree, periodic })
    }

    /// `Geom2d_BSplineCurve::IsRational()`.
    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
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

    /// `Geom2d_BSplineCurve::SetKnots(K)` (`Geom2d_BSplineCurve.cxx:...`):
    /// `CheckCurveData(myPoles, K, myMults, myDeg, myPeriodic)` followed by
    /// `updateKnots()`. `distinct_knots` holds the distinct knots, exactly like
    /// the 3D `GeomBSplineCurve::set_knots`.
    pub fn set_knots(&mut self, distinct_knots: &[f64]) -> Result<(), &'static str> {
        let (_, umults) = self.distinct_knots_and_mults();
        if umults.len() != distinct_knots.len() {
            return Err("Geom2dBSplineCurve::set_knots: knot count mismatch");
        }
        if knots::nb_poles(self.degree as i32, self.periodic, &umults) as usize != self.xs.len() {
            return Err("Geom2dBSplineCurve::set_knots: pole/degree mismatch");
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

    /// `Geom2d_BSplineCurve::FirstUKnotIndex()` / `LastUKnotIndex()`
    /// (`Geom2d_BSplineCurve_1.cxx:335-347`, `:405-417`): `1` / `NbKnots()` on
    /// a periodic curve, otherwise the `BSplCLib` indices of the multiplicities.
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

    /// `Geom2d_BSplineCurve::updateKnots()` (`Geom2d_BSplineCurve.cxx:1280-1324`):
    /// rebuild the flat knot vector from the distinct knots and their
    /// multiplicities (`BSplCLib::KnotSequence`; the periodic representation
    /// uses the period-extended sequence).
    ///
    /// OCCT's `KnotSet == GeomAbs_Uniform && !Periodic` shortcut assigns the
    /// distinct array directly; that array is what the non-periodic expansion
    /// reproduces for an all-multiplicity-one curve, so the port expands in
    /// both cases (same equivalence as `GeomBSplineCurve::update_flat_knots`).
    pub(super) fn update_flat_knots(&mut self, uknots: &[f64], umults: &[i32]) {
        self.knots = if self.periodic {
            knots::knot_sequence_periodic(uknots, umults, self.degree as i32)
        } else {
            occt_core::bspl::banded_interp::knot_sequence(uknots, umults, self.degree as i32)
        };
    }

    pub fn nb_poles(&self) -> usize { self.xs.len() }
    pub fn degree(&self) -> usize { self.degree }
    pub fn first_parameter(&self) -> f64 { self.knots[self.degree] }
    pub fn last_parameter(&self) -> f64 { self.knots[self.knots.len() - 1 - self.degree] }

    /// `Geom2d_BSplineCurve::FirstUKnotIndex()` (`Geom2d_BSplineCurve_1.cxx:335-347`).
    pub fn first_u_knot_index(&self) -> i32 {
        let (_, umults) = self.distinct_knots_and_mults();
        self.u_knot_index_range(&umults).0
    }

    /// `Geom2d_BSplineCurve::LastUKnotIndex()` (`Geom2d_BSplineCurve_1.cxx:405-417`).
    pub fn last_u_knot_index(&self) -> i32 {
        let (_, umults) = self.distinct_knots_and_mults();
        self.u_knot_index_range(&umults).1
    }

    /// `Geom2d_BSplineCurve::NbKnots()` (`Geom2d_BSplineCurve_1.cxx:598-601`).
    pub fn nb_knots(&self) -> usize {
        self.distinct_knots_and_mults().0.len()
    }

    /// `Geom2d_BSplineCurve::Weight(Index)` 1-based (`Geom2d_BSplineCurve.cxx:194-206`).
    pub fn weight(&self, index: usize) -> f64 {
        self.weights
            .as_ref()
            .map(|w| w[index - 1])
            .unwrap_or(1.0)
    }

    /// `Geom2d_BSplineCurve::SetWeight(Index, Weight)` 1-based
    /// (`Geom2d_BSplineCurve.cxx:534-...`).
    pub fn set_weight(&mut self, index: usize, weight: f64) {
        if index == 0 || index > self.xs.len() {
            return;
        }
        if self.weights.is_none() {
            self.weights = Some(vec![1.0; self.xs.len()]);
        }
        if let Some(w) = self.weights.as_mut() {
            w[index - 1] = weight;
        }
    }

    /// `Geom2d_BSplineCurve::Resolution(ToleranceUV, UTolerance)`
    /// (`Geom2d_BSplineCurve_1.cxx:764-806`): parametric tolerance for a
    /// given UV tolerance. Poles are lifted to `z = 0` so the shared
    /// `BSplCLib::Resolution` (dim-3) matches the 2D L1 bound.
    pub fn resolution(&self, tolerance_uv: f64) -> f64 {
        let poles = self.poles_3d();
        if self.periodic {
            // `PrepareUnperiodize` pole wrap (`cxx:768-791`).
            let (_, umults) = self.distinct_knots_and_mults();
            let (_nk, np) = occt_core::bspl::unperiodize::prepare_unperiodize(
                self.degree as i32,
                &umults,
            );
            let n_old = poles.len();
            if n_old == 0 || np == 0 {
                return tolerance_uv / occt_core::precision::REAL_SMALL;
            }
            let new_poles: Vec<GpPnt> = (0..np).map(|i| poles[i % n_old]).collect();
            let new_weights = self
                .weights
                .as_ref()
                .map(|w| (0..np).map(|i| w[i % n_old]).collect::<Vec<_>>());
            occt_core::bspl::bspline_curve_resolution(
                &new_poles,
                new_weights.as_deref(),
                &self.knots,
                self.degree as i32,
                tolerance_uv,
            )
        } else {
            occt_core::bspl::bspline_curve_resolution(
                &poles,
                self.weights.as_deref(),
                &self.knots,
                self.degree as i32,
                tolerance_uv,
            )
        }
    }

    /// `Geom2d_BSplineCurve::LocateU` without knot repetition
    /// (`Geom2d_BSplineCurve_1.cxx:710-760`, `WithKnotRepetition = false`).
    /// Returns `(I1, I2)` in 1-based distinct-knot indices.
    pub fn locate_u(&self, u: f64, parametric_tolerance: f64) -> (i32, i32) {
        let (uknots, _) = self.distinct_knots_and_mults();
        if uknots.is_empty() {
            return (0, 0);
        }
        let mut new_u = u;
        if self.periodic {
            let first = self.first_parameter();
            let last = self.last_parameter();
            new_u = occt_core::bspl::locate::in_period(new_u, first, last);
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
        let mut i1 = occt_core::bspl::locate::hunt_occt(&uknots, new_u);
        i1 = i1.max(1).min(uknots.len() as i32);
        while (i1 as usize) < uknots.len()
            && (uknots[i1 as usize] - new_u).abs() <= tol
        {
            i1 += 1;
        }
        if (uknots[(i1 - 1) as usize] - new_u).abs() <= tol {
            (i1, i1)
        } else {
            (i1, i1 + 1)
        }
    }

    /// Lift `(x, y)` poles to `z = 0` for `BSplCLib` evaluators.
    pub(super) fn poles_3d(&self) -> Vec<GpPnt> {
        self.xs
            .iter()
            .zip(self.ys.iter())
            .map(|(&x, &y)| GpPnt::new(x, y, 0.0))
            .collect()
    }

    /// `Geom2d_BSplineCurve::EvalD0` span resolution, mirroring
    /// `GeomBSplineCurve::eval_span` in `occt-geom`
    /// (`Geom_BSplineCurve_1.cxx:174-194` -> `BSplCLib::D0`). Returns the
    /// parameter mapped into the period, the span index, and the stored
    /// distinct knots/multiplicities that `BSplCLib` consumes.
    pub(super) fn eval_span(&self, u: f64) -> (f64, i32, Vec<f64>, Vec<i32>) {
        let (uknots, umults) = self.distinct_knots_and_mults();
        let (mut span, new_u) = occt_core::bspl::locate::locate_parameter(
            self.degree as i32,
            &uknots,
            Some(&umults),
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
        (new_u, span, uknots, umults)
    }
}
