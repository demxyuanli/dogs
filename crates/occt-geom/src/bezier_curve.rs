//! 3D Bezier curve.

use crate::curve::Curve;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

/// Rational-free polynomial Bezier curve in 3D.
#[derive(Clone)]
pub struct GeomBezierCurve {
    pub poles: Vec<GpPnt>,
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
    pub fn new(poles: Vec<GpPnt>) -> Result<Self, &'static str> {
        if poles.len() < 2 {
            return Err("GeomBezierCurve requires at least 2 poles");
        }
        Ok(Self { poles })
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
    /// UNPORTED: the rational branch (`PLib::Trimming` with weights) — this curve
    /// type carries no weights.
    pub fn segment(&mut self, u1: f64, u2: f64) {
        let degree = self.degree();
        let n = degree + 1;
        let mut coefs = vec![0.0f64; 3 * n];
        for (i, slot) in coefs.chunks_mut(3).enumerate() {
            let mut p = GpPnt::zero();
            for j in 0..n {
                let b = bezier_power_entry(degree, i, j);
                p = GpPnt::new(
                    p.x() + b * self.poles[j].x(),
                    p.y() + b * self.poles[j].y(),
                    p.z() + b * self.poles[j].z(),
                );
            }
            slot[0] = p.x();
            slot[1] = p.y();
            slot[2] = p.z();
        }
        occt_core::bspl::plib::trimming(u1, u2, 3, &mut coefs);
        let to_poles = power_to_bezier_matrix(degree);
        for i in 0..n {
            let (mut x, mut y, mut z) = (0.0f64, 0.0f64, 0.0f64);
            for j in 0..n {
                let b = to_poles[i * n + j];
                x += b * coefs[3 * j];
                y += b * coefs[3 * j + 1];
                z += b * coefs[3 * j + 2];
            }
            self.poles[i] = GpPnt::new(x, y, z);
        }
    }

    pub fn set_pole(&mut self, i: usize, p: GpPnt) {
        self.poles[i] = p;
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

impl Curve for GeomBezierCurve {
    fn d0(&self, u: f64) -> GpPnt {
        de_casteljau(&self.poles, u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        (self.d0(u), self.tangent(u))
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        (self.d0(u), self.tangent(u), self.second_deriv(u))
    }

    fn first_parameter(&self) -> f64 {
        0.0
    }

    fn last_parameter(&self) -> f64 {
        1.0
    }

    fn continuity(&self) -> u8 {
        3
    }

    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.poles {
            *p = p.transformed(t);
        }
    }

    fn reverse(&mut self) {
        self.poles.reverse();
    }

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }

    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        Some(&self.poles)
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
        occt_core::bspl::bspline_curve_resolution(&self.poles, None, &flat_knots, degree, r3d)
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
