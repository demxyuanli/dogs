//! 3D Bezier curve.

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

/// Rational-free polynomial Bezier curve in 3D.
#[derive(Clone)]
pub struct GeomBezierCurve {
    pub poles: Vec<GpPnt>,
    /// `Geom_BezierCurve` weights. `None` is the polynomial curve.
    pub weights: Option<Vec<f64>>,
    /// `myClosed` (`Geom_BezierCurve.cxx:774-775`, `:390`, `:440`).
    closed: bool,
}

fn lerp(a: &GpPnt, b: &GpPnt, t: f64) -> GpPnt {
    GpPnt::new(
        a.x() + (b.x() - a.x()) * t,
        a.y() + (b.y() - a.y()) * t,
        a.z() + (b.z() - a.z()) * t,
    )
}

fn lerp_vec(a: &GpVec, b: &GpVec, t: f64) -> GpVec {
    GpVec::new(
        a.x() + (b.x() - a.x()) * t,
        a.y() + (b.y() - a.y()) * t,
        a.z() + (b.z() - a.z()) * t,
    )
}

fn de_casteljau(poles: &[GpPnt], u: f64) -> GpPnt {
    let mut pts: Vec<GpPnt> = poles.to_vec();
    while pts.len() > 1 {
        for i in 0..(pts.len() - 1) {
            pts[i] = lerp(&pts[i], &pts[i + 1], u);
        }
        pts.pop();
    }
    pts[0]
}

fn de_casteljau_vec(vs: &[GpVec], u: f64) -> GpVec {
    let mut pts: Vec<GpVec> = vs.to_vec();
    while pts.len() > 1 {
        for i in 0..(pts.len() - 1) {
            pts[i] = lerp_vec(&pts[i], &pts[i + 1], u);
        }
        pts.pop();
    }
    pts[0]
}

/// `Geom_BezierCurve::Rational` (`Geom_BezierCurve.cxx:60-72`).
fn weights_differ(weights: &[f64]) -> bool {
    for i in 0..weights.len().saturating_sub(1) {
        if (weights[i] - weights[i + 1]).abs() > occt_core::precision::REAL_SMALL {
            return true;
        }
    }
    false
}

/// Binomial coefficient `C(n, k)` as `u64`.
fn binomial(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut c = 1u64;
    for i in 0..k {
        c = c * (n - i) as u64 / (i + 1) as u64;
    }
    c
}

/// Bezier → monomial matrix entry `(i, j)`: the coefficient of `t^i` in
/// `B_j^d(t) = C(d,j)·t^j·(1−t)^(d−j)`, i.e.
/// `(−1)^(i−j)·C(d,j)·C(d−j, i−j)` (zero for `i < j`). This is the expansion
/// `BSplCLib::BuildCache` produces for a Bezier on `[0, 1]`.
fn bezier_power_entry(degree: usize, i: usize, j: usize) -> f64 {
    if i < j {
        return 0.0;
    }
    let sign = if (i - j) % 2 == 0 { 1.0 } else { -1.0 };
    sign * binomial(degree, j) as f64 * binomial(degree - j, i - j) as f64
}

/// Inverse of the Bezier → monomial matrix — the monomial → Bezier write-back of
/// `PLib::CoefficientsPoles` (`PLib.cxx:1482-1493` + the `dim`-arm conversion).
/// The forward matrix is triangular with nonzero diagonal, so forward
/// substitution gives the inverse exactly in `f64`.
fn power_to_bezier_matrix(degree: usize) -> Vec<f64> {
    let n = degree + 1;
    let a: Vec<f64> = (0..n * n)
        .map(|k| bezier_power_entry(degree, k / n, k % n))
        .collect();
    let mut inv = vec![0.0f64; n * n];
    for row in 0..n {
        for col in 0..=row {
            let mut v = if row == col { 1.0 } else { 0.0 };
            for k in col..row {
                v -= a[row * n + k] * inv[k * n + col];
            }
            inv[row * n + col] = v / a[row * n + row];
        }
    }
    inv
}

impl GeomBezierCurve {
    /// `Geom_BezierCurve::MaxDegree` (`Geom_BezierCurve.cxx:196-199`).
    pub fn max_degree() -> i32 {
        25
    }

    pub fn new(poles: Vec<GpPnt>) -> Result<Self, &'static str> {
        if poles.len() < 2 {
            return Err("GeomBezierCurve requires at least 2 poles");
        }
        let closed = poles[0].distance(poles.last().unwrap()) <= occt_core::precision::CONFUSION;
        Ok(Self {
            poles,
            weights: None,
            closed,
        })
    }

    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }

    pub fn degree(&self) -> usize {
        self.poles.len() - 1
    }

    pub fn pole(&self, i: usize) -> &GpPnt {
        &self.poles[i]
    }

    /// `Geom_BezierCurve::Weight` (`Geom_BezierCurve.cxx:688-698`).
    /// `Index` is 1-based. A polynomial curve reports weight 1.
    pub fn weight(&self, index: i32) -> Result<f64, &'static str> {
        let nb = self.poles.len() as i32;
        if index < 1 || index > nb {
            return Err("Geom_BezierCurve::Weight");
        }
        Ok(self
            .weights
            .as_ref()
            .map(|w| w[(index as usize) - 1])
            .unwrap_or(1.0))
    }

    /// `Geom_BezierCurve::IsClosed` (`Geom_BezierCurve.cxx:499-502`).
    pub fn is_closed(&self) -> bool {
        self.closed
    }

    fn closed_from_poles(poles: &[GpPnt]) -> bool {
        match (poles.first(), poles.last()) {
            (Some(a), Some(b)) => a.distance(b) <= occt_core::precision::CONFUSION,
            _ => false,
        }
    }

    /// Equivalent B-spline on `[0, 1]` with full end multiplicities
    /// (`Geom_BezierCurve.cxx:595`, `:618`).
    fn equivalent_bspline(&self) -> Option<GeomBSplineCurve> {
        let n = self.poles.len();
        if n < 2 {
            return None;
        }
        let degree = n - 1;
        let mult = (degree as i32) + 1;
        if let Some(weights) = self.weights.as_ref() {
            if weights.len() != n {
                return None;
            }
            let knots = occt_core::bspl::banded_interp::knot_sequence(
                &[0.0, 1.0],
                &[mult, mult],
                degree as i32,
            );
            GeomBSplineCurve::rational(self.poles.clone(), weights.clone(), knots, degree).ok()
        } else {
            GeomBSplineCurve::from_poles_knots_mults(
                self.poles.clone(),
                vec![0.0, 1.0],
                vec![mult, mult],
                degree,
            )
            .ok()
        }
    }

    /// `Geom_BezierCurve::Segment(U1, U2)` (`Geom_BezierCurve.cxx:388-425`),
    /// non-rational branch:
    ///
    /// ```text
    /// BSplCLib::BuildCache(0., 1., false, aDeg, KnotSequence(), myPoles,
    ///                      BSplCLib::NoWeights(), coeffs, BSplCLib::NoWeights());
    /// PLib::Trimming(U1, U2, coeffs, PLib::NoWeights());
    /// PLib::CoefficientsPoles(coeffs, PLib::NoWeights(), myPoles, PLib::NoWeights());
    /// ```
    ///
    /// A Bezier curve lives on `[0, 1]`, so its `BuildCache` coefficients are the
    /// monomial (power-basis) expansion — `B_j^d(t) = C(d,j)·t^j·(1−t)^(d−j)`
    /// gives `M[i][j] = (−1)^(i−j)·C(d,j)·C(d−j, i−j)` — and
    /// `PLib::CoefficientsPoles` is the inverse expansion. `PLib::Trimming`
    /// (`PLib.cxx:1642-1716`) re-expresses the sub-range in the same basis.
    /// The rational arm (`Geom_BezierCurve.cxx:395-408`) applies that matrix
    /// to `myWeights` and trims poles and weights together.
    pub fn segment(&mut self, u1: f64, u2: f64) {
        // `myClosed` from the untrimmed ends (`Geom_BezierCurve.cxx:390`).
        self.closed = self.d0(u1).distance(&self.d0(u2)) <= occt_core::precision::CONFUSION;
        let degree = self.degree();
        let n = degree + 1;
        let rational = self.weights.as_ref().map(|w| w.len() == n).unwrap_or(false);
        let mut coefs = vec![0.0f64; 3 * n];
        let mut wcoefs = rational.then(|| vec![0.0; n]);
        for i in 0..n {
            let mut p = GpPnt::zero();
            let mut weight = 0.0;
            for j in 0..n {
                let b = bezier_power_entry(degree, i, j);
                p = GpPnt::new(
                    p.x() + b * self.poles[j].x(),
                    p.y() + b * self.poles[j].y(),
                    p.z() + b * self.poles[j].z(),
                );
                if let Some(src) = self.weights.as_ref() {
                    weight += b * src[j];
                }
            }
            coefs[3 * i] = p.x();
            coefs[3 * i + 1] = p.y();
            coefs[3 * i + 2] = p.z();
            if let Some(dst) = wcoefs.as_mut() {
                dst[i] = weight;
            }
        }
        if let Some(wcoefs) = wcoefs.as_mut() {
            occt_core::bspl::plib::trimming_w(u1, u2, 3, &mut coefs, Some(wcoefs));
        } else {
            occt_core::bspl::plib::trimming(u1, u2, 3, &mut coefs);
        }
        let to_poles = power_to_bezier_matrix(degree);
        let mut next_w = rational.then(|| vec![0.0; n]);
        for i in 0..n {
            let (mut x, mut y, mut z) = (0.0f64, 0.0f64, 0.0f64);
            let mut weight = 0.0;
            for j in 0..n {
                let b = to_poles[i * n + j];
                x += b * coefs[3 * j];
                y += b * coefs[3 * j + 1];
                z += b * coefs[3 * j + 2];
                if let Some(src) = wcoefs.as_ref() {
                    weight += b * src[j];
                }
            }
            self.poles[i] = GpPnt::new(x, y, z);
            if let Some(dst) = next_w.as_mut() {
                dst[i] = weight;
            }
        }
        if let Some(next_w) = next_w {
            self.weights = Some(next_w);
        }
    }

    /// `Geom_BezierCurve::SetPole` (`Geom_BezierCurve.cxx:429-443`).
    /// `i` is 0-based. End poles refresh `myClosed` (`:438-441`).
    pub fn set_pole(&mut self, i: usize, p: GpPnt) {
        self.poles[i] = p;
        if i == 0 || i + 1 == self.poles.len() {
            self.closed = Self::closed_from_poles(&self.poles);
        }
    }

    /// `Geom_BezierCurve::Increase` (`Geom_BezierCurve.cxx:141-192`).
    /// `init` after elevation may drop unit weights.
    pub fn increase(&mut self, deg: i32) -> Result<(), &'static str> {
        let old = self.degree() as i32;
        if deg == old {
            return Ok(());
        }
        if deg < old || deg > Self::max_degree() {
            return Err("Geom_BezierCurve::Increase");
        }
        let knots = [0.0, 1.0];
        let mults = [old + 1, old + 1];
        let nb_new_knots = occt_core::bspl::increase_degree::increase_degree_count_knots(
            old, deg, false, &mults,
        )
        .max(0) as usize;
        let nb_new_poles = (deg as usize) + 1;
        let raised = occt_core::bspl::increase_degree::increase_degree(
            old,
            deg,
            false,
            &self.poles,
            self.weights.as_deref(),
            &knots,
            &mults,
            nb_new_poles,
            nb_new_knots,
        );
        self.poles = raised.poles;
        if let Some(weights) = raised.weights {
            if weights_differ(&weights) {
                self.weights = Some(weights);
            } else {
                self.weights = None;
            }
        } else {
            self.weights = None;
        }
        Ok(())
    }

    /// `Geom_BezierCurve::InsertPoleAfter` (`Geom_BezierCurve.cxx:210-280`).
    /// `index` is 1-based; 0 inserts before the first pole (`InsertPoleBefore(1)`).
    pub fn insert_pole_after(
        &mut self,
        index: i32,
        point: GpPnt,
        weight: f64,
    ) -> Result<(), &'static str> {
        let nbpoles = self.poles.len() as i32;
        if nbpoles >= Self::max_degree() || weight <= occt_core::precision::REAL_SMALL {
            return Err("Geom_BezierCurve::InsertPoleAfter");
        }
        if index < 0 || index > nbpoles {
            return Err("Geom_BezierCurve::InsertPoleAfter");
        }
        let insert_at = index as usize;
        let mut poles = Vec::with_capacity(self.poles.len() + 1);
        poles.extend_from_slice(&self.poles[..insert_at]);
        poles.push(point);
        poles.extend_from_slice(&self.poles[insert_at..]);
        let rat = self.is_rational()
            || (weight - 1.0).abs() > occt_core::precision::REAL_SMALL;
        if rat {
            let mut weights = Vec::with_capacity(poles.len());
            for i in 0..insert_at {
                weights.push(self.weights.as_ref().map(|w| w[i]).unwrap_or(1.0));
            }
            weights.push(weight);
            for i in insert_at..self.poles.len() {
                weights.push(self.weights.as_ref().map(|w| w[i]).unwrap_or(1.0));
            }
            if weights_differ(&weights) {
                self.weights = Some(weights);
            } else {
                self.weights = None;
            }
        } else {
            self.weights = None;
        }
        self.poles = poles;
        Ok(())
    }

    /// `Geom_BezierCurve::InsertPoleBefore` (`Geom_BezierCurve.cxx:284-294`).
    pub fn insert_pole_before(
        &mut self,
        index: i32,
        point: GpPnt,
        weight: f64,
    ) -> Result<(), &'static str> {
        self.insert_pole_after(index - 1, point, weight)
    }

    /// `Geom_BezierCurve::RemovePole` (`Geom_BezierCurve.cxx:298-340`).
    pub fn remove_pole(&mut self, index: i32) -> Result<(), &'static str> {
        let nbpoles = self.poles.len() as i32;
        if nbpoles <= 2 {
            return Err("Geom_BezierCurve::RemovePole");
        }
        if index < 1 || index > nbpoles {
            return Err("Geom_BezierCurve::RemovePole");
        }
        let at = (index as usize) - 1;
        self.poles.remove(at);
        if let Some(weights) = self.weights.as_mut() {
            weights.remove(at);
            if !weights_differ(weights) {
                self.weights = None;
            }
        }
        Ok(())
    }

    /// `Geom_BezierCurve::Reverse` (`Geom_BezierCurve.cxx:352-377`).
    /// `myMaxDerivInvOk` (`:376`) is not stored.
    pub fn reverse(&mut self) {
        self.poles.reverse();
        if let Some(weights) = self.weights.as_mut() {
            weights.reverse();
        }
    }

    /// `Geom_BezierCurve::IsRational` (`Geom_BezierCurve.cxx:520-523`).
    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    /// Weighted constructor (`Geom_BezierCurve.cxx:104-136`).
    /// Weights that do not vary are stored as a polynomial curve.
    pub fn with_weights(poles: Vec<GpPnt>, weights: Vec<f64>) -> Result<Self, &'static str> {
        let nb = poles.len();
        if nb < 2 || nb > (Self::max_degree() as usize) + 1 {
            return Err("Geom_BezierCurve");
        }
        if weights.len() != nb {
            return Err("Geom_BezierCurve");
        }
        for &weight in &weights {
            if weight <= occt_core::precision::REAL_SMALL {
                return Err("Geom_BezierCurve");
            }
        }
        let mut curve = Self::new(poles)?;
        if weights_differ(&weights) {
            curve.weights = Some(weights);
        }
        Ok(curve)
    }

    /// `Geom_BezierCurve::SetWeight` (`Geom_BezierCurve.cxx:455-495`).
    /// `Index` is 1-based. `myMaxDerivInvOk` (`:494`) is not stored.
    pub fn set_weight(&mut self, index: i32, weight: f64) -> Result<(), &'static str> {
        let nbpoles = self.poles.len() as i32;
        if index < 1 || index > nbpoles {
            return Err("Geom_BezierCurve::SetWeight");
        }
        if weight <= occt_core::precision::REAL_SMALL {
            return Err("Geom_BezierCurve::SetWeight");
        }
        let was_rat = self.is_rational();
        if !was_rat && (weight - 1.0).abs() <= occt_core::precision::REAL_SMALL {
            return Ok(());
        }
        if !was_rat {
            self.weights = Some(vec![1.0; self.poles.len()]);
        }
        if let Some(weights) = self.weights.as_mut() {
            weights[(index as usize) - 1] = weight;
        }
        let still = self.weights.as_ref().map(|w| weights_differ(w)).unwrap_or(false);
        if was_rat && !still {
            self.weights = None;
        }
        Ok(())
    }

    fn tangent(&self, u: f64) -> GpVec {
        let deg = self.degree();
        if deg == 0 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let d = deg as f64;
        let mut dp: Vec<GpVec> = Vec::with_capacity(deg);
        for i in 0..deg {
            let p0 = &self.poles[i];
            let p1 = &self.poles[i + 1];
            dp.push(GpVec::new(
                (p1.x() - p0.x()) * d,
                (p1.y() - p0.y()) * d,
                (p1.z() - p0.z()) * d,
            ));
        }
        de_casteljau_vec(&dp, u)
    }

    fn second_deriv(&self, u: f64) -> GpVec {
        let deg = self.degree();
        if deg < 2 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let d = deg as f64;
        let mut dp: Vec<GpVec> = Vec::with_capacity(deg);
        for i in 0..deg {
            let p0 = &self.poles[i];
            let p1 = &self.poles[i + 1];
            dp.push(GpVec::new(
                (p1.x() - p0.x()) * d,
                (p1.y() - p0.y()) * d,
                (p1.z() - p0.z()) * d,
            ));
        }
        let mut ddp: Vec<GpVec> = Vec::with_capacity(deg - 1);
        for i in 0..(deg - 1) {
            ddp.push(GpVec::new(
                (dp[i + 1].x() - dp[i].x()) * (d - 1.0),
                (dp[i + 1].y() - dp[i].y()) * (d - 1.0),
                (dp[i + 1].z() - dp[i].z()) * (d - 1.0),
            ));
        }
        de_casteljau_vec(&ddp, u)
    }
}

impl GeomBezierCurve {
    fn rational_deriv(&self, u: f64) -> Option<(GpPnt, GpVec, GpVec)> {
        let w = self.weights.as_ref()?;
        if w.len() != self.poles.len() || self.poles.len() < 2 {
            return None;
        }
        let degree = self.degree();
        let mult = (degree as i32) + 1;
        let knots = occt_core::bspl::banded_interp::knot_sequence(
            &[0.0, 1.0],
            &[mult, mult],
            degree as i32,
        );
        Some(occt_core::bspl::eval::eval_curve_rational_d2(
            &self.poles,
            w,
            &knots,
            degree,
            u,
        ))
    }
}

impl Curve for GeomBezierCurve {
    fn d0(&self, u: f64) -> GpPnt {
        if let Some((point, _, _)) = self.rational_deriv(u) {
            return point;
        }
        de_casteljau(&self.poles, u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        if let Some((point, d1, _)) = self.rational_deriv(u) {
            return (point, d1);
        }
        (self.d0(u), self.tangent(u))
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        if let Some(deriv) = self.rational_deriv(u) {
            return deriv;
        }
        (self.d0(u), self.tangent(u), self.second_deriv(u))
    }

    /// `Geom_BezierCurve::EvalD3` (`Geom_BezierCurve.cxx:586-597`).
    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        if let Some(bs) = self.equivalent_bspline() {
            return bs.d3(u);
        }
        let (p, d1, d2) = self.d2(u);
        (p, d1, d2, GpVec::zero())
    }

    /// `Geom_BezierCurve::EvalDN` (`Geom_BezierCurve.cxx:601-620`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        if n < 1 {
            return GpVec::zero();
        }
        if let Some(bs) = self.equivalent_bspline() {
            return bs.eval_dn(u, n);
        }
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            3 => self.d3(u).3,
            _ => GpVec::zero(),
        }
    }

    fn first_parameter(&self) -> f64 {
        0.0
    }

    fn last_parameter(&self) -> f64 {
        1.0
    }

    /// `Geom_BezierCurve::IsPeriodic` (`Geom_BezierCurve.cxx:513-516`).
    fn is_periodic(&self) -> bool {
        false
    }

    /// `Geom_BezierCurve::Continuity` (`Geom_BezierCurve.cxx:527-530`).
    fn continuity(&self) -> u8 {
        occt_core::kernel::geomabs::Shape::CN as u8
    }

    /// `Geom_BezierCurve::IsCN` (`Geom_BezierCurve.cxx:506-509`).
    fn is_cn(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_Curve::IsCN: negative order".into());
        }
        Ok(true)
    }

    /// `Geom_BezierCurve::Transform` (`Geom_BezierCurve.cxx:729-739`).
    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.poles {
            *p = p.transformed(t);
        }
    }

    /// `Geom_BezierCurve::Reverse` (`Geom_BezierCurve.cxx:352-377`).
    fn reverse(&mut self) {
        GeomBezierCurve::reverse(self);
    }

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }

    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        Some(&self.poles)
    }
    fn bezier_weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }
    fn nurbs_degree(&self) -> Option<usize> {
        Some(self.poles.len().saturating_sub(1))
    }

    /// `Geom_BezierCurve::Resolution` (`Geom_BezierCurve.cxx`): the first call
    /// computes `myMaxDerivInv` with `BSplCLib::Resolution(poles, weights,
    /// nbpoles, KnotSequence(), degree, 1., inv)` and then
    /// `UTolerance = Tolerance3D * myMaxDerivInv`. `BSplCLib::Resolution` is
    /// linear in its tolerance argument, so calling the faithful port
    /// (`occt_core::bspl::bspline_curve_resolution`) with `Tolerance3D`
    /// directly is the same value. The Bezier knot sequence is `deg+1` zeros
    /// followed by `deg+1` ones (`Geom_BezierCurve::KnotSequence`).
    fn resolution(&self, r3d: f64) -> f64 {
        let n = self.poles.len();
        if n < 2 {
            // `Geom_BezierCurve` requires at least two poles; OCCT's degenerate
            // path (through `RealSmall()`) is what the helper returns for it.
            return occt_core::bspl::bspline_curve_resolution(&self.poles, None, &[], 0, r3d);
        }
        let degree = (n - 1) as i32;
        let mut flat_knots = vec![0.0; n];
        flat_knots.extend(std::iter::repeat(1.0).take(n));
        occt_core::bspl::bspline_curve_resolution(
            &self.poles,
            self.weights.as_deref(),
            &flat_knots,
            degree,
            r3d,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: &GpPnt, b: &GpPnt) -> bool {
        (a.x() - b.x()).abs() < 1e-12
            && (a.y() - b.y()).abs() < 1e-12
            && (a.z() - b.z()).abs() < 1e-12
    }

    #[test]
    fn quadratic_bezier() {
        let c = GeomBezierCurve::new(vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 2.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
        ])
        .unwrap();

        assert!(approx(&c.d0(0.5), &GpPnt::new(1.0, 1.0, 0.0)));
        assert!(approx(&c.d0(0.0), &GpPnt::new(0.0, 0.0, 0.0)));
        assert!(approx(&c.d0(1.0), &GpPnt::new(2.0, 0.0, 0.0)));
    }

    #[test]
    fn requires_two_poles() {
        assert!(GeomBezierCurve::new(vec![GpPnt::new(0.0, 0.0, 0.0)]).is_err());
    }
}
