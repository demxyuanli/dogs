//! B-spline evaluation via de Boor algorithm. Source: `BSplCLib.cxx` — Eval, EvalDeriv
//!
//! Given control points (poles), weights, knots, and degree, evaluate
//! the rational B-spline curve/surface at parameter u.

use crate::gp::GpPnt;

/// Evaluate non-rational B-spline curve at u.
/// Returns point on curve. Uses de Boor recursion (Cox-de Boor).
/// poles: control points (n). knots: knot vector (n+deg+1). degree: d.
pub fn eval_curve(poles: &[GpPnt], knots: &[f64], degree: usize, u: f64) -> GpPnt {
    let n = poles.len();
    if n == 0 { return GpPnt::zero(); }
    let idx = super::knots::hunt(knots, u).max(degree).min(n - 1);
    // Extract local knot span
    let mut pts = vec![GpPnt::zero(); degree + 1];
    for k in 0..=degree {
        let pi = idx - degree + k;
        pts[k] = if pi < poles.len() { poles[pi] } else { poles[poles.len()-1] };
    }
    // de Boor triangle
    for r in 1..=degree {
        for i in (r..=degree).rev() {
            let k0 = idx + i - degree;
            let k1 = k0 + degree + 1 - r;
            let alpha = (u - knots[k0]) / (knots[k1] - knots[k0]);
            if alpha.is_finite() {
                let x = (1.0 - alpha) * pts[i-1].coord.x + alpha * pts[i].coord.x;
                let y = (1.0 - alpha) * pts[i-1].coord.y + alpha * pts[i].coord.y;
                let z = (1.0 - alpha) * pts[i-1].coord.z + alpha * pts[i].coord.z;
                pts[i] = GpPnt::new(x, y, z);
            }
        }
    }
    pts[degree]
}

/// Evaluate rational B-spline curve at u (with weights).
pub fn eval_curve_rational(poles: &[GpPnt], weights: &[f64], knots: &[f64], degree: usize, u: f64) -> GpPnt {
    let n = poles.len();
    if n == 0 { return GpPnt::zero(); }
    let idx = super::knots::hunt(knots, u).max(degree).min(n - 1);
    let mut pts = vec![(GpXyz::zero(), 0.0f64); degree + 1];
    for k in 0..=degree {
        let pi = (idx - degree + k).min(poles.len()-1);
        pts[k] = (poles[pi].coord, weights[pi]);
    }
    for r in 1..=degree {
        for i in (r..=degree).rev() {
            let k0 = idx + i - degree;
            let k1 = k0 + degree + 1 - r;
            let alpha = (u - knots[k0]) / (knots[k1] - knots[k0]);
            if alpha.is_finite() {
                pts[i].0 = crate::gp::GpXyz::new(
                    (1.0 - alpha) * pts[i-1].0.x + alpha * pts[i].0.x,
                    (1.0 - alpha) * pts[i-1].0.y + alpha * pts[i].0.y,
                    (1.0 - alpha) * pts[i-1].0.z + alpha * pts[i].0.z,
                );
                pts[i].1 = (1.0 - alpha) * pts[i-1].1 + alpha * pts[i].1;
            }
        }
    }
    if pts[degree].1.abs() < 1e-30 {
        GpPnt::new(pts[degree].0.x, pts[degree].0.y, pts[degree].0.z)
    } else {
        let w = pts[degree].1;
        GpPnt::new(pts[degree].0.x / w, pts[degree].0.y / w, pts[degree].0.z / w)
    }
}

/// Evaluate non-rational B-spline curve D0 and D1 at u.
/// Returns (point, first derivative).
pub fn eval_curve_d1(poles: &[GpPnt], knots: &[f64], degree: usize, u: f64) -> (GpPnt, crate::gp::GpVec) {
    let n = poles.len();
    if n < 2 { return (GpPnt::zero(), crate::gp::GpVec::zero()); }
    // Derivative by differencing degree-1 spline
    let deg1 = if degree > 0 { degree - 1 } else { 0 };
    let mut dpoles = Vec::with_capacity(n - 1);
    for i in 0..n-1 {
        let k0 = i + 1;
        let k1 = k0 + degree;
        let alpha = if knots[k1] > knots[k0] { degree as f64 / (knots[k1] - knots[k0]) } else { 1.0 };
        let dx = alpha * (poles[i+1].coord.x - poles[i].coord.x);
        let dy = alpha * (poles[i+1].coord.y - poles[i].coord.y);
        let dz = alpha * (poles[i+1].coord.z - poles[i].coord.z);
        dpoles.push(GpPnt::new(dx, dy, dz));
    }
    let pt = eval_curve(poles, knots, degree, u);
    let d1 = if dpoles.is_empty() { crate::gp::GpVec::zero() } else {
        let p = eval_curve(&dpoles, &knots[1..knots.len()-1], deg1, u);
        crate::gp::GpVec::new(p.x(), p.y(), p.z())
    };
    (pt, d1)
}

use crate::gp::GpXyz;
use crate::gp::GpVec;

/// Evaluate non-rational B-spline surface at (u,v).
/// poles: flat array of control points, row-major (n_u rows, n_v cols).
pub fn eval_surface(poles: &[GpPnt], n_u: usize, n_v: usize,
                    knots_u: &[f64], knots_v: &[f64],
                    degree_u: usize, degree_v: usize,
                    u: f64, v: f64) -> GpPnt {
    // Compute v-direction intermediate curve
    let n = n_v;
    let mut curve_u = Vec::with_capacity(n);
    for j in 0..n {
        let col: Vec<GpPnt> = (0..n_u).map(|i| poles[i * n + j]).collect();
        curve_u.push(eval_curve(&col, knots_u, degree_u, u));
    }
    eval_curve(&curve_u, knots_v, degree_v, v)
}

/// Evaluate non-rational B-spline surface D0 and D1 at (u,v).
pub fn eval_surface_d1(poles: &[GpPnt], n_u: usize, n_v: usize,
                       knots_u: &[f64], knots_v: &[f64],
                       degree_u: usize, degree_v: usize,
                       u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    let n = n_v;
    let mut curve_u = Vec::with_capacity(n);
    let mut dcurve_u = Vec::with_capacity(n);
    for j in 0..n {
        let col: Vec<GpPnt> = (0..n_u).map(|i| poles[i * n + j]).collect();
        let (pt, du) = eval_curve_d1(&col, knots_u, degree_u, u);
        curve_u.push(pt);
        dcurve_u.push(du);
    }
    let pt = eval_curve(&curve_u, knots_v, degree_v, v);
    let du = GpVec::new(
        eval_curve(&dcurve_u.iter().map(|p| GpPnt::new(p.x(), 0., 0.)).collect::<Vec<_>>(), knots_v, degree_v, v).x(),
        eval_curve(&dcurve_u.iter().map(|p| GpPnt::new(p.y(), 0., 0.)).collect::<Vec<_>>(), knots_v, degree_v, v).x(),
        eval_curve(&dcurve_u.iter().map(|p| GpPnt::new(p.z(), 0., 0.)).collect::<Vec<_>>(), knots_v, degree_v, v).x(),
    );
    let (_, dv) = eval_curve_d1(&curve_u, knots_v, degree_v, v);
    (pt, du, dv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::GpPnt;

    #[test]
    fn eval_linear() {
        let poles = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.)];
        let knots = vec![0.,0.,1.,1.];
        let p = eval_curve(&poles, &knots, 1, 0.5);
        assert!((p.x() - 0.5).abs() < 1e-14);
    }

    #[test]
    fn eval_cubic_endpoints() {
        let poles = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(3.,1.,0.)];
        let knots = crate::bspl::knots::build_uniform_knots(4, 3);
        let p0 = eval_curve(&poles, &knots, 3, 0.0);
        assert!((p0.x() - 0.0).abs() < 1e-14);
        let p1 = eval_curve(&poles, &knots, 3, 1.0);
        assert!((p1.x() - 3.0).abs() < 1e-14);
    }
}
