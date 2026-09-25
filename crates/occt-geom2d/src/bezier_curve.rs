//! 2D Bezier curve (rational and non-rational). Source: `Geom2d_BezierCurve.hxx`

use crate::curve::Curve2d;
use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};

/// 2D Bezier curve. Non-rational when `weights` is `None`, rational otherwise.
/// Parameter range is `[0, 1]`; degree = `poles.len() - 1`.
#[derive(Clone)]
pub struct Geom2dBezierCurve {
    pub poles: Vec<GpPnt2d>,
    pub weights: Option<Vec<f64>>,
}

impl Geom2dBezierCurve {
    /// Non-rational Bezier from control points (OCCT constructor).
    pub fn new(poles: Vec<GpPnt2d>) -> Result<Self, &'static str> {
        if poles.len() < 2 {
            return Err("Geom2dBezierCurve: need at least 2 poles");
        }
        Ok(Self { poles, weights: None })
    }

    /// Rational Bezier from control points and positive weights.
    pub fn rational(poles: Vec<GpPnt2d>, weights: Vec<f64>) -> Result<Self, &'static str> {
        if poles.len() < 2 {
            return Err("Geom2dBezierCurve: need at least 2 poles");
        }
        if weights.len() != poles.len() {
            return Err("Geom2dBezierCurve: weights/poles length mismatch");
        }
        if weights.iter().any(|&w| w <= 1e-12) {
            return Err("Geom2dBezierCurve: weights must be positive");
        }
        Ok(Self { poles, weights: Some(weights) })
    }

    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }

    /// Polynomial degree = number of poles minus one.
    pub fn degree(&self) -> usize {
        self.poles.len() - 1
    }

    /// Pole by 0-based index.
    pub fn pole(&self, i: usize) -> GpPnt2d {
        self.poles[i]
    }
}

impl Curve2d for Geom2dBezierCurve {
    fn d0(&self, u: f64) -> GpPnt2d {
        match &self.weights {
            Some(w) => rational_de_casteljau2d(&self.poles, w, u),
            None => de_casteljau2d(&self.poles, u),
        }
    }

    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let p = self.d0(u);
        let n = self.degree();
        let v = match &self.weights {
            Some(w) => {
                // C' = (N' w - N w') / w^2 with N = weighted numerator curve.
                let h: Vec<GpPnt2d> = self
                    .poles
                    .iter()
                    .zip(w.iter())
                    .map(|(pp, wi)| pnt_scale(*pp, *wi))
                    .collect();
                let h1: Vec<GpPnt2d> = (0..n)
                    .map(|i| pnt_scale(pnt_sub(h[i + 1], h[i]), n as f64))
                    .collect();
                let w1: Vec<f64> = (0..n).map(|i| (w[i + 1] - w[i]) * n as f64).collect();
                let nv = de_casteljau2d(&h1, u);
                let wv = de_casteljau_scalar(&w1, u);
                let wt = de_casteljau_scalar(w, u);
                let denom = wt * wt;
                GpVec2d::new((nv.x() * wt - p.x() * wv) / denom, (nv.y() * wt - p.y() * wv) / denom)
            }
            None => {
                // Derivative control points: Q_i = n (P_{i+1} - P_i).
                let q1: Vec<GpPnt2d> = (0..n)
                    .map(|i| pnt_scale(pnt_sub(self.poles[i + 1], self.poles[i]), n as f64))
                    .collect();
                let d = de_casteljau2d(&q1, u);
                GpVec2d::new(d.x(), d.y())
            }
        };
        (p, v)
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, v1) = self.d1(u);
        let n = self.degree();
        let v2 = match &self.weights {
            Some(w) => {
                // C'' = (N'' w - 2 N' w' + C (2 w'^2 - w w'')) / w^2
                let h: Vec<GpPnt2d> = self
                    .poles
                    .iter()
                    .zip(w.iter())
                    .map(|(pp, wi)| pnt_scale(*pp, *wi))
                    .collect();
                let h1: Vec<GpPnt2d> = (0..n)
                    .map(|i| pnt_scale(pnt_sub(h[i + 1], h[i]), n as f64))
                    .collect();
                let h2: Vec<GpPnt2d> = if n >= 2 {
                    (0..n - 1)
                        .map(|i| {
                            let d = pnt_sub(pnt_add(h[i + 2], h[i]), pnt_scale(h[i + 1], 2.0));
                            pnt_scale(d, (n * (n - 1)) as f64)
                        })
                        .collect()
                } else {
                    Vec::new()
                };
                let w1: Vec<f64> = (0..n).map(|i| (w[i + 1] - w[i]) * n as f64).collect();
                let w2: Vec<f64> = if n >= 2 {
                    (0..n - 1)
                        .map(|i| (w[i + 2] - 2.0 * w[i + 1] + w[i]) * (n * (n - 1)) as f64)
                        .collect()
                } else {
                    Vec::new()
                };
                let n2v = if h2.is_empty() { GpPnt2d::zero() } else { de_casteljau2d(&h2, u) };
                let w2v = if w2.is_empty() { 0.0 } else { de_casteljau_scalar(&w2, u) };
                let nv = de_casteljau2d(&h1, u);
                let wv = de_casteljau_scalar(&w1, u);
                let wt = de_casteljau_scalar(w, u);
                let denom = wt * wt;
                let s = 2.0 * wv * wv - wt * w2v;
                let numer = pnt_add(pnt_add(pnt_scale(n2v, wt), pnt_scale(nv, -2.0 * wv)), pnt_scale(p, s));
                GpVec2d::new(numer.x() / denom, numer.y() / denom)
            }
            None => {
                if n < 2 {
                    return (p, v1, GpVec2d::zero());
                }
                // Second-difference control points: Q_i = n(n-1)(P_{i+2} - 2 P_{i+1} + P_i).
                let q2: Vec<GpPnt2d> = (0..n - 1)
                    .map(|i| {
                        let d = pnt_sub(
                            pnt_add(self.poles[i + 2], self.poles[i]),
                            pnt_scale(self.poles[i + 1], 2.0),
                        );
                        pnt_scale(d, (n * (n - 1)) as f64)
                    })
                    .collect();
                let d = de_casteljau2d(&q2, u);
                GpVec2d::new(d.x(), d.y())
            }
        };
        (p, v1, v2)
    }

    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        1.0
    }
    fn is_periodic(&self) -> bool {
        false
    }
    fn continuity(&self) -> u8 {
        3
    }
    fn transform(&mut self, t: &GpTrsf2d) {
        for p in self.poles.iter_mut() {
            p.transform(t);
        }
    }
    fn reverse(&mut self) {
        self.poles.reverse();
        if let Some(w) = self.weights.as_mut() {
            w.reverse();
        }
    }
    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }

    /// `Geom2d_BezierCurve::NbPoles()` (`Geom2d_BezierCurve.cxx:600-603`).
    fn bezier_nb_poles(&self) -> Option<usize> {
        Some(self.poles.len())
    }

    /// `Geom2dAdaptor_Curve::IsRational()` (`Geom2dAdaptor_Curve.cxx:1297-1298`):
    /// `Geom2d_BezierCurve::IsRational()` is true once weights were supplied.
    fn is_rational(&self) -> bool {
        self.weights.is_some()
    }
}

/// de Casteljau evaluation of a non-rational Bezier.
pub fn de_casteljau2d(poles: &[GpPnt2d], t: f64) -> GpPnt2d {
    if poles.is_empty() {
        return GpPnt2d::zero();
    }
    let mut pts = poles.to_vec();
    let mut m = pts.len();
    while m > 1 {
        for i in 0..m - 1 {
            pts[i] = GpPnt2d::new(
                (1.0 - t) * pts[i].x() + t * pts[i + 1].x(),
                (1.0 - t) * pts[i].y() + t * pts[i + 1].y(),
            );
        }
        m -= 1;
    }
    pts[0]
}

/// Rational de Casteljau evaluation (homogeneous form).
pub fn rational_de_casteljau2d(poles: &[GpPnt2d], weights: &[f64], t: f64) -> GpPnt2d {
    if poles.is_empty() {
        return GpPnt2d::zero();
    }
    let mut pts = poles.to_vec();
    let mut w = weights.to_vec();
    let mut m = pts.len();
    while m > 1 {
        for i in 0..m - 1 {
            let a = (1.0 - t) * w[i];
            let b = t * w[i + 1];
            let wi = a + b;
            if wi.abs() < 1e-300 {
                pts[i] = GpPnt2d::new(
                    (1.0 - t) * pts[i].x() + t * pts[i + 1].x(),
                    (1.0 - t) * pts[i].y() + t * pts[i + 1].y(),
                );
                w[i] = wi;
            } else {
                pts[i] = GpPnt2d::new(
                    (a * pts[i].x() + b * pts[i + 1].x()) / wi,
                    (a * pts[i].y() + b * pts[i + 1].y()) / wi,
                );
                w[i] = wi;
            }
        }
        m -= 1;
    }
    pts[0]
}

/// Sample a Bezier curve as a polyline by recursive subdivision until the
/// control polygon is within `tol` of its chord.
pub fn bezier2d_to_polyline(c: &Geom2dBezierCurve, tol: f64) -> Vec<GpPnt2d> {
    let mut out = Vec::new();
    if c.poles.is_empty() {
        return out;
    }
    rec(&c.poles, c.weights.as_deref(), tol, &mut out);
    // Consecutive leaves share the de Casteljau midpoint — dedupe.
    let mut result = Vec::with_capacity(out.len());
    for p in out {
        if result.last().map_or(true, |q: &GpPnt2d| q.distance(&p) > 1e-12) {
            result.push(p);
        }
    }
    result
}

fn rec(poles: &[GpPnt2d], weights: Option<&[f64]>, tol: f64, out: &mut Vec<GpPnt2d>) {
    let n = poles.len();
    if n == 1 {
        out.push(poles[0]);
        return;
    }
    if is_flat(poles, tol) {
        out.push(poles[0]);
        out.push(poles[n - 1]);
        return;
    }
    if let Some(w) = weights {
        let ((lp, lw), (rp, rw)) = split_rational(poles, w, 0.5);
        rec(&lp, Some(&lw), tol, out);
        rec(&rp, Some(&rw), tol, out);
    } else {
        let (l, r) = split_bezier(poles, 0.5);
        rec(&l, None, tol, out);
        rec(&r, None, tol, out);
    }
}

fn is_flat(poles: &[GpPnt2d], tol: f64) -> bool {
    let a = poles[0];
    let b = poles[poles.len() - 1];
    poles[1..poles.len() - 1]
        .iter()
        .all(|p| dist_to_line(*p, a, b) <= tol)
}

fn dist_to_line(p: GpPnt2d, a: GpPnt2d, b: GpPnt2d) -> f64 {
    let dx = b.x() - a.x();
    let dy = b.y() - a.y();
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-30 {
        return p.distance(&a);
    }
    let t = ((p.x() - a.x()) * dx + (p.y() - a.y()) * dy) / len2;
    let px = a.x() + t * dx;
    let py = a.y() + t * dy;
    ((p.x() - px).powi(2) + (p.y() - py).powi(2)).sqrt()
}

/// Split a non-rational Bezier at parameter t into (left, right) control nets.
fn split_bezier(poles: &[GpPnt2d], t: f64) -> (Vec<GpPnt2d>, Vec<GpPnt2d>) {
    let n = poles.len();
    let mut q: Vec<Vec<GpPnt2d>> = vec![poles.to_vec()];
    for level in 1..n {
        let prev = &q[level - 1];
        let cur: Vec<GpPnt2d> = (0..n - level)
            .map(|i| GpPnt2d::new(
                (1.0 - t) * prev[i].x() + t * prev[i + 1].x(),
                (1.0 - t) * prev[i].y() + t * prev[i + 1].y(),
            ))
            .collect();
        q.push(cur);
    }
    let left: Vec<GpPnt2d> = q.iter().map(|row| row[0]).collect();
    let mut right: Vec<GpPnt2d> = (0..n).map(|level| q[level][n - 1 - level]).collect();
    right.reverse();
    (left, right)
}

/// Split a rational Bezier at parameter t into (left, right) pole+weight nets.
fn split_rational(poles: &[GpPnt2d], weights: &[f64], t: f64) -> ((Vec<GpPnt2d>, Vec<f64>), (Vec<GpPnt2d>, Vec<f64>)) {
    let n = poles.len();
    let mut hp: Vec<Vec<GpPnt2d>> = vec![poles.to_vec()];
    let mut hw: Vec<Vec<f64>> = vec![weights.to_vec()];
    for level in 1..n {
        let prevp = &hp[level - 1];
        let prevw = &hw[level - 1];
        let mut cp = Vec::with_capacity(n - level);
        let mut cw = Vec::with_capacity(n - level);
        for i in 0..n - level {
            let a = (1.0 - t) * prevw[i];
            let b = t * prevw[i + 1];
            let w = a + b;
            cp.push(GpPnt2d::new(
                (a * prevp[i].x() + b * prevp[i + 1].x()) / w,
                (a * prevp[i].y() + b * prevp[i + 1].y()) / w,
            ));
            cw.push(w);
        }
        hp.push(cp);
        hw.push(cw);
    }
    let left_p: Vec<GpPnt2d> = hp.iter().map(|r| r[0]).collect();
    let left_w: Vec<f64> = hw.iter().map(|r| r[0]).collect();
    let mut right_p: Vec<GpPnt2d> = (0..n).map(|l| hp[l][n - 1 - l]).collect();
    let mut right_w: Vec<f64> = (0..n).map(|l| hw[l][n - 1 - l]).collect();
    right_p.reverse();
    right_w.reverse();
    ((left_p, left_w), (right_p, right_w))
}

fn de_casteljau_scalar(c: &[f64], t: f64) -> f64 {
    if c.is_empty() {
        return 0.0;
    }
    let mut v = c.to_vec();
    let mut m = v.len();
    while m > 1 {
        for i in 0..m - 1 {
            v[i] = (1.0 - t) * v[i] + t * v[i + 1];
        }
        m -= 1;
    }
    v[0]
}

fn pnt_add(a: GpPnt2d, b: GpPnt2d) -> GpPnt2d {
    GpPnt2d::new(a.x() + b.x(), a.y() + b.y())
}
fn pnt_sub(a: GpPnt2d, b: GpPnt2d) -> GpPnt2d {
    GpPnt2d::new(a.x() - b.x(), a.y() - b.y())
}
fn pnt_scale(a: GpPnt2d, s: f64) -> GpPnt2d {
    GpPnt2d::new(a.x() * s, a.y() * s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_midpoint() {
        let c = Geom2dBezierCurve::new(vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(0.5, 1.0),
            GpPnt2d::new(1.0, 0.0),
        ])
        .unwrap();
        let p = c.d0(0.5);
        assert!((p.x() - 0.5).abs() < 1e-12, "x={}", p.x());
        assert!((p.y() - 0.5).abs() < 1e-12, "y={}", p.y());
    }

    #[test]
    fn endpoints_match_poles() {
        let c = Geom2dBezierCurve::new(vec![
            GpPnt2d::new(1.0, 2.0),
            GpPnt2d::new(2.0, 3.0),
            GpPnt2d::new(4.0, 5.0),
        ])
        .unwrap();
        assert!(c.d0(0.0).distance(&c.pole(0)) < 1e-12);
        assert!(c.d0(1.0).distance(&c.pole(2)) < 1e-12);
    }

    #[test]
    fn reverse_flips_endpoints() {
        let mut c = Geom2dBezierCurve::new(vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(2.0, 0.0),
        ])
        .unwrap();
        let orig_start = c.d0(0.0);
        let orig_end = c.d0(1.0);
        c.reverse();
        assert!(c.d0(0.0).distance(&orig_end) < 1e-12);
        assert!(c.d0(1.0).distance(&orig_start) < 1e-12);
    }

    #[test]
    fn transform_translates_poles() {
        let mut c = Geom2dBezierCurve::new(vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
        ])
        .unwrap();
        let mut t = GpTrsf2d::identity();
        t.set_translation_vec(&GpVec2d::new(1.0, 2.0));
        c.transform(&t);
        assert!((c.pole(0).x() - 1.0).abs() < 1e-12 && (c.pole(0).y() - 2.0).abs() < 1e-12);
        assert!((c.pole(1).x() - 2.0).abs() < 1e-12 && (c.pole(1).y() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn rational_quarter_circle() {
        // Quadratic rational Bezier for a quarter unit circle.
        let c = Geom2dBezierCurve::rational(
            vec![
                GpPnt2d::new(1.0, 0.0),
                GpPnt2d::new(1.0, 1.0),
                GpPnt2d::new(0.0, 1.0),
            ],
            vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0],
        )
        .unwrap();
        let p = c.d0(0.5);
        // Midpoint of the quarter circle = (sqrt(2)/2, sqrt(2)/2).
        let s = std::f64::consts::FRAC_1_SQRT_2;
        assert!((p.x() - s).abs() < 1e-9, "x={}", p.x());
        assert!((p.y() - s).abs() < 1e-9, "y={}", p.y());
    }

    #[test]
    fn polyline_endpoints() {
        let c = Geom2dBezierCurve::new(vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(0.5, 1.0),
            GpPnt2d::new(1.0, 0.0),
        ])
        .unwrap();
        let poly = bezier2d_to_polyline(&c, 1e-3);
        assert!(poly.len() >= 3);
        assert!(poly.first().unwrap().distance(&GpPnt2d::new(0.0, 0.0)) < 1e-9);
        assert!(poly.last().unwrap().distance(&GpPnt2d::new(1.0, 0.0)) < 1e-9);
    }
}
