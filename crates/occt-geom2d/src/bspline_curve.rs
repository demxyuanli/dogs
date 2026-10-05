//! 2D B-spline curve (rational and non-rational). Source: `Geom2d_BSplineCurve.hxx`

use crate::curve::Curve2d;
use occt_core::bspl::{eval, knots};
use occt_core::gp::{GpPnt, GpPnt2d, GpTrsf2d, GpVec2d};

/// `Geom2d_BSplineCurve::MaxDegree()` = `BSplCLib::MaxDegree()`
/// (`BSplCLib.lxx:24-27`), same value as the 3D port's `MAX_DEGREE`.
const MAX_DEGREE: usize = 25;

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

    /// `Geom2d_BSplineCurve::FirstUKnotIndex()` / `LastUKnotIndex()`
    /// (`Geom2d_BSplineCurve_1.cxx:335-347`, `:405-417`): `1` / `NbKnots()` on
    /// a periodic curve, otherwise the `BSplCLib` indices of the multiplicities.
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

    /// `Geom2d_BSplineCurve::updateKnots()` (`Geom2d_BSplineCurve.cxx:1280-1324`):
    /// rebuild the flat knot vector from the distinct knots and their
    /// multiplicities (`BSplCLib::KnotSequence`; the periodic representation
    /// uses the period-extended sequence).
    ///
    /// OCCT's `KnotSet == GeomAbs_Uniform && !Periodic` shortcut assigns the
    /// distinct array directly; that array is what the non-periodic expansion
    /// reproduces for an all-multiplicity-one curve, so the port expands in
    /// both cases (same equivalence as `GeomBSplineCurve::update_flat_knots`).
    fn update_flat_knots(&mut self, uknots: &[f64], umults: &[i32]) {
        self.knots = if self.periodic {
            knots::knot_sequence_periodic(uknots, umults, self.degree as i32)
        } else {
            occt_core::bspl::banded_interp::knot_sequence(uknots, umults, self.degree as i32)
        };
    }

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
    pub fn segment(&mut self, u1: f64, u2: f64, the_tolerance: f64) -> Result<(), &'static str> {
        if u2 < u1 {
            return Err("Geom2d_BSplineCurve::Segment");
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
                return Err("Geom2d_BSplineCurve::Segment");
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
            if a_nu2 < u2 {
                nknots[last_k] = u1 + a_ddu;
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
                u, 0, 0, self.degree as i32, true, &self.poles_3d(), self.weights.as_deref(),
                &self.knots, None,
            );
            return GpPnt2d::new(v.x(), v.y());
        }
        match &self.weights {
            Some(w) => {
                let p = eval::eval_curve_rational(&self.poles_3d(), w, &self.knots, self.degree, u);
                GpPnt2d::new(p.x(), p.y())
            }
            None => self.de_boor(u),
        }
    }

    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        // `Geom2d_BSplineCurve::D1` / `BSplCLib::D1`.
        let poles = self.poles_3d();
        if self.periodic {
            let d = occt_core::bspl::curve_dn::dn(
                u, 1, 0, self.degree as i32, true, &poles, self.weights.as_deref(), &self.knots,
                None,
            );
            return (self.d0(u), GpVec2d::new(d.x(), d.y()));
        }
        if let Some(w) = self.weights.as_deref() {
            // Rational uses the homogeneous quotient already computed by
            // `eval_curve_rational_d2`.
            let (p, d1, _) = eval::eval_curve_rational_d2(&poles, w, &self.knots, self.degree, u);
            return (GpPnt2d::new(p.x(), p.y()), GpVec2d::new(d1.x(), d1.y()));
        }
        let (p, d) = eval::eval_curve_d1(&poles, &self.knots, self.degree, u);
        (GpPnt2d::new(p.x(), p.y()), GpVec2d::new(d.x(), d.y()))
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        // `Geom2d_BSplineCurve::D2` / `BSplCLib::D2`.
        let poles = self.poles_3d();
        if self.periodic {
            let d1 = occt_core::bspl::curve_dn::dn(
                u, 1, 0, self.degree as i32, true, &poles, self.weights.as_deref(), &self.knots,
                None,
            );
            let d2 = occt_core::bspl::curve_dn::dn(
                u, 2, 0, self.degree as i32, true, &poles, self.weights.as_deref(), &self.knots,
                None,
            );
            return (
                self.d0(u),
                GpVec2d::new(d1.x(), d1.y()),
                GpVec2d::new(d2.x(), d2.y()),
            );
        }
        if let Some(w) = self.weights.as_deref() {
            let (p, d1, d2) = eval::eval_curve_rational_d2(&poles, w, &self.knots, self.degree, u);
            return (
                GpPnt2d::new(p.x(), p.y()),
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
        if let Some(w) = self.weights.as_mut() {
            w.reverse();
        }
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
    /// `Geom2d_BSplineCurve::IsRational()` (`Geom2dAdaptor_Curve.cxx:1296-1298`).
    fn is_rational(&self) -> bool {
        self.weights.is_some()
    }
    /// `Geom2d_BSplineCurve::Weights()` (`GeomLib.cxx:598-620`).
    fn bspline_weights2d(&self) -> Option<&[f64]> {
        self.weights.as_deref()
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

    /// `Geom2d_BSplineCurve::Copy()` (`Geom2d_BSplineCurve.cxx:109-136`): a deep
    /// copy of the stored poles/knots. The port's clone is that copy.
    fn bspline_copy2d(&self) -> Option<Geom2dBSplineCurve> {
        Some(self.clone())
    }

    /// `Geom2d_BSplineCurve::Knots()` / `Multiplicities()` as the distinct
    /// knots + multiplicities pair that the `Geom_BSplineCurve` constructor
    /// takes (`GeomLib.cxx:610-618`). The port stores no weights, so the
    /// rational branch of `GeomLib::To3d` cannot be reached from here.
    fn bspline_distinct_knots_mults(&self) -> Option<(Vec<f64>, Vec<i32>)> {
        Some(self.distinct_knots_and_mults())
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
