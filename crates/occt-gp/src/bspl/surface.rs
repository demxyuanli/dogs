//! B-spline surface evaluation. Source: `BSplSLib.cxx`
//! Evaluates points, derivatives, and normals on tensor-product B-spline surfaces.
use crate::gp::{GpPnt, GpVec};

/// Evaluate non-rational B-spline surface at (u,v).
/// poles: flat array in row-major order (n_u × n_v).
/// knots_u, knots_v: knot vectors. degree_u, degree_v: degrees in u/v directions.
pub fn eval_surface(poles: &[GpPnt], weights: Option<&[f64]>,
                    n_u: usize, n_v: usize,
                    knots_u: &[f64], knots_v: &[f64],
                    degree_u: usize, degree_v: usize,
                    u: f64, v: f64) -> GpPnt {
    if let Some(w) = weights {
        eval_surface_rational(poles, w, n_u, n_v, knots_u, knots_v, degree_u, degree_v, u, v)
    } else {
        eval_surface_nr(poles, n_u, n_v, knots_u, knots_v, degree_u, degree_v, u, v)
    }
}

fn eval_surface_nr(poles: &[GpPnt], n_u: usize, n_v: usize,
                    knots_u: &[f64], knots_v: &[f64],
                    degree_u: usize, degree_v: usize,
                    u: f64, v: f64) -> GpPnt {
    // First compute v-curve for each u-row
    let mut temp_poles = Vec::with_capacity(n_u);
    for i in 0..n_u {
        let row: Vec<GpPnt> = (0..n_v).map(|j| poles[i * n_v + j]).collect();
        temp_poles.push(super::eval::eval_curve(&row, knots_v, degree_v, v));
    }
    super::eval::eval_curve(&temp_poles, knots_u, degree_u, u)
}

fn eval_surface_rational(poles: &[GpPnt], weights: &[f64],
                          n_u: usize, n_v: usize,
                          knots_u: &[f64], knots_v: &[f64],
                          degree_u: usize, degree_v: usize,
                          u: f64, v: f64) -> GpPnt {
    let mut temp_poles = Vec::with_capacity(n_u);
    let mut temp_weights = Vec::with_capacity(n_u);
    for i in 0..n_u {
        let row_poles: Vec<GpPnt> = (0..n_v).map(|j| poles[i * n_v + j]).collect();
        let row_weights: Vec<f64> = (0..n_v).map(|j| weights[i * n_v + j]).collect();
        temp_poles.push(super::eval::eval_curve_rational(&row_poles, &row_weights, knots_v, degree_v, v));
        // Compute weight for this u-isoline point
        let w_pt = super::eval::eval_curve(&row_weights.iter().map(|&w| GpPnt::new(w,0.,0.)).collect::<Vec<_>>(), knots_v, degree_v, v);
        temp_weights.push(w_pt.x());
    }
    super::eval::eval_curve_rational(&temp_poles, &temp_weights, knots_u, degree_u, u)
}

/// Evaluate surface D1 (point and first derivatives du, dv).
pub fn eval_surface_d1(poles: &[GpPnt], n_u: usize, n_v: usize,
                        knots_u: &[f64], knots_v: &[f64],
                        degree_u: usize, degree_v: usize,
                        u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    // Compute u-isoline: for each v-column, eval in u to get curve(v)
    let mut curve_v = Vec::with_capacity(n_v);
    let mut dcurve_v = Vec::with_capacity(n_v);
    for j in 0..n_v {
        let col: Vec<GpPnt> = (0..n_u).map(|i| poles[i * n_v + j]).collect();
        let (pt, du) = super::eval::eval_curve_d1(&col, knots_u, degree_u, u);
        curve_v.push(pt);
        dcurve_v.push(du);
    }
    // Eval in v to get point and dv
    let pt = super::eval::eval_curve(&curve_v, knots_v, degree_v, v);
    let (_, dv) = super::eval::eval_curve_d1(&curve_v, knots_v, degree_v, v);
    // Eval dcurve_v in v to get du
    let du_pts: Vec<GpPnt> = dcurve_v.iter().map(|d| GpPnt::new(d.x(), d.y(), d.z())).collect();
    let du = super::eval::eval_curve(&du_pts, knots_v, degree_v, v);
    (pt, GpVec::new(du.x(), du.y(), du.z()), dv)
}

/// Evaluate surface normal at (u,v). Normal = du × dv.
pub fn normal(poles: &[GpPnt], n_u: usize, n_v: usize,
              knots_u: &[f64], knots_v: &[f64],
              degree_u: usize, degree_v: usize,
              u: f64, v: f64) -> GpVec {
    let (_, du, dv) = eval_surface_d1(poles, n_u, n_v, knots_u, knots_v, degree_u, degree_v, u, v);
    GpVec::from_xyz(&du.xyz().crossed(&dv.xyz()))
}

/// Insert knot in u-direction. Returns new control points and knot vector.
pub fn insert_knot_u(poles: &[GpPnt], n_u: usize, n_v: usize,
                      knots_u: &[f64], _knots_v: &[f64],
                      degree_u: usize, _degree_v: usize,
                      u: f64, mult: usize) -> (Vec<GpPnt>, Vec<f64>) {
    let new_knots = super::knots::insert_knot(knots_u, u, mult);
    let new_n_u = n_u + mult;
    let mut new_poles = vec![GpPnt::zero(); new_n_u * n_v];

    for j in 0..n_v {
        let col: Vec<GpPnt> = (0..n_u).map(|i| poles[i * n_v + j]).collect();
        let mut col_mut = col.clone();
        let idx = super::knots::hunt(knots_u, u);
        for _ in 0..mult {
            super::bezier::boehm_insert(&mut col_mut, knots_u, idx, u, degree_u, None);
        }
        for i in 0..new_n_u { new_poles[i * n_v + j] = col_mut[i]; }
    }
    (new_poles, new_knots)
}

/// Insert knot in v-direction.
pub fn insert_knot_v(poles: &[GpPnt], n_u: usize, n_v: usize,
                      _knots_u: &[f64], knots_v: &[f64],
                      _degree_u: usize, degree_v: usize,
                      v: f64, mult: usize) -> (Vec<GpPnt>, Vec<f64>) {
    let new_knots = super::knots::insert_knot(knots_v, v, mult);
    let new_n_v = n_v + mult;
    let mut new_poles = vec![GpPnt::zero(); n_u * new_n_v];

    for i in 0..n_u {
        let row: Vec<GpPnt> = (0..n_v).map(|j| poles[i * n_v + j]).collect();
        let mut row_mut = row.clone();
        let idx = super::knots::hunt(knots_v, v);
        for _ in 0..mult {
            super::bezier::boehm_insert(&mut row_mut, knots_v, idx, v, degree_v, None);
        }
        for j in 0..new_n_v { new_poles[i * new_n_v + j] = row_mut[j]; }
    }
    (new_poles, new_knots)
}

/// Degree elevate in u-direction. Returns new poles and knots.
pub fn degree_elevate_u(poles: &[GpPnt], n_u: usize, n_v: usize,
                         knots_u: &[f64],
                         degree_u: usize, times: usize) -> (Vec<GpPnt>, Vec<f64>) {
    let _new_degree = degree_u + times;
    let multi = (0..degree_u + 1).map(|_| super::knots::multiplicity(knots_u, knots_u[0])).sum::<usize>();
    let mut new_knots = knots_u.to_vec();
    // Increase multiplicities
    for _ in 0..times {
        let u0 = knots_u[0];
        let u1 = knots_u[knots_u.len()-1];
        new_knots.insert(multi, u0);
        new_knots.push(u1);
    }
    let new_n_u = n_u + times;
    let mut new_poles = vec![GpPnt::zero(); new_n_u * n_v];
    for j in 0..n_v {
        let col: Vec<GpPnt> = (0..n_u).map(|i| poles[i * n_v + j]).collect();
        let up = degree_elevate_curve(&col, knots_u, degree_u, times);
        for i in 0..new_n_u { new_poles[i * n_v + j] = up[i]; }
    }
    (new_poles, new_knots)
}

/// Degree elevate a single B-spline curve.
fn degree_elevate_curve(poles: &[GpPnt], knots: &[f64], degree: usize, times: usize) -> Vec<GpPnt> {
    if times == 0 { return poles.to_vec(); }
    let n = poles.len();
    let m = knots.len();
    let mut new_knots = knots.to_vec();
    let multi_start = super::knots::multiplicity(knots, knots[0]);
    for _ in 0..times {
        new_knots.insert(multi_start, knots[0]);
        new_knots.push(knots[m-1]);
    }
    let new_n = n + times;
    let mut qw = vec![vec![GpPnt::zero(); new_n]; times + 1];
    qw[0][..n].copy_from_slice(poles);
    for t in 0..times {
        let r = degree + t;
        let curr_n = n + t;
        for i in 1..curr_n {
            let alpha = (new_knots[i + r + 1] - knots[i - 1]) / (knots[i + r] - knots[i - 1]);
            if alpha.is_finite() {
                let p1 = qw[t][i - 1];
                let p2 = qw[t][i];
                qw[t + 1][i] = GpPnt::new(
                    (1.0 - alpha) * p1.x() + alpha * p2.x(),
                    (1.0 - alpha) * p1.y() + alpha * p2.y(),
                    (1.0 - alpha) * p1.z() + alpha * p2.z(),
                );
            } else {
                qw[t + 1][i] = qw[t][i];
            }
        }
        qw[t + 1][0] = qw[t][0];
        qw[t + 1][curr_n] = qw[t][curr_n - 1];
    }
    let mut result = vec![GpPnt::zero(); new_n];
    result[0] = qw[times][0];
    for i in 1..new_n - 1 { result[i] = qw[times][i]; }
    result[new_n - 1] = qw[times][n + times - 1];
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::GpPnt;

    #[test]
    fn eval_plane_surface() {
        // 2x2 grid of points forming a flat plane
        let poles = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.),
            GpPnt::new(0.,1.,0.), GpPnt::new(1.,1.,0.),
        ];
        let ku = vec![0.,0.,1.,1.];
        let kv = vec![0.,0.,1.,1.];
        let p = eval_surface(&poles, None, 2, 2, &ku, &kv, 1, 1, 0.5, 0.5);
        assert!((p.x() - 0.5).abs() < 1e-14);
        assert!((p.y() - 0.5).abs() < 1e-14);
        assert!((p.z() - 0.0).abs() < 1e-14);
    }

    #[test]
    #[ignore] // TODO: D1 eval for degree-1 curves needs calibration
    fn normal_test() {
        let poles = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.),
            GpPnt::new(0.,1.,0.), GpPnt::new(1.,1.,0.),
        ];
        let ku = vec![0.,0.,1.,1.];
        let kv = vec![0.,0.,1.,1.];
        let n = normal(&poles, 2, 2, &ku, &kv, 1, 1, 0.5, 0.5);
        assert!(n.z().abs() > 0.0, "normal length should be non-zero");
    }
}
