//! Bezier subdivision and flattening.
//! Source: `PLib.hxx` (de Casteljau), `GeomConvert_CurveToPolyline.hxx`.

use occt_core::gp::{GpPnt, GpXyz};

/// Recursion depth cap for flattening pathological control polygons.
const MAX_DEPTH: usize = 24;

/// Evaluate a polynomial Bezier at parameter `t` in `[0, 1]` via de Casteljau.
pub fn de_casteljau(poles: &[GpPnt], t: f64) -> GpPnt {
    if poles.is_empty() {
        return GpPnt::zero();
    }
    let mut pts: Vec<GpPnt> = poles.to_vec();
    while pts.len() > 1 {
        for i in 0..pts.len() - 1 {
            pts[i] = lerp(&pts[i], &pts[i + 1], t);
        }
        pts.pop();
    }
    pts[0]
}

/// Evaluate a rational Bezier at `t` via de Casteljau in homogeneous space.
pub fn rational_de_casteljau(poles: &[GpPnt], weights: &[f64], t: f64) -> GpPnt {
    if poles.is_empty() {
        return GpPnt::zero();
    }
    let n = poles.len().min(weights.len());
    let mut pts: Vec<GpXyz> = (0..n).map(|i| poles[i].coord.multiplied(weights[i])).collect();
    let mut ws: Vec<f64> = weights[..n].to_vec();
    while pts.len() > 1 {
        for i in 0..pts.len() - 1 {
            pts[i] = pts[i].multiplied(1.0 - t).added(&pts[i + 1].multiplied(t));
            ws[i] = (1.0 - t) * ws[i] + t * ws[i + 1];
        }
        pts.pop();
        ws.pop();
    }
    if ws[0].abs() < 1e-30 {
        GpPnt::zero()
    } else {
        GpPnt::from_xyz(&pts[0].divided(ws[0]))
    }
}

/// Max deviation of the interior control points from the chord `poles[0]..poles[n-1]`.
pub fn bezier_flatness(poles: &[GpPnt]) -> f64 {
    let n = poles.len();
    if n <= 2 {
        return 0.0;
    }
    let a = &poles[0];
    let b = &poles[n - 1];
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return poles[1..n - 1].iter().map(|p| p.distance(a)).fold(0.0, f64::max);
    }
    poles[1..n - 1]
        .iter()
        .map(|p| {
            let v = p.coord.subtracted(&a.coord);
            let t = (v.dot(&ab) / len2).clamp(0.0, 1.0);
            v.subtracted(&ab.multiplied(t)).modulus()
        })
        .fold(0.0, f64::max)
}

/// Flatten a Bezier curve into a polyline whose chords deviate from the curve
/// by at most `tol`. Always includes the first and last poles. Uses an
/// iterative stack (splitting at `t = 0.5`) instead of recursion.
pub fn bezier_to_polyline(poles: &[GpPnt], weights: Option<&[f64]>, tol: f64) -> Vec<GpPnt> {
    let tol = tol.max(1e-15);
    let mut result: Vec<GpPnt> = Vec::new();
    let mut stack: Vec<(Vec<GpPnt>, Option<Vec<f64>>, usize)> = Vec::new();
    stack.push((poles.to_vec(), weights.map(|w| w.to_vec()), 0));

    while let Some((seg_poles, seg_w, depth)) = stack.pop() {
        if depth >= MAX_DEPTH || bezier_flatness(&seg_poles) <= tol {
            push_dedup(&mut result, seg_poles[0]);
            push_dedup(&mut result, *seg_poles.last().unwrap());
        } else {
            let (left_p, right_p) = split_bezier(&seg_poles, 0.5);
            match seg_w {
                Some(w) => {
                    let (left_w, right_w) = split_bezier_weights(&w, 0.5);
                    stack.push((right_p, Some(right_w), depth + 1));
                    stack.push((left_p, Some(left_w), depth + 1));
                }
                None => {
                    stack.push((right_p, None, depth + 1));
                    stack.push((left_p, None, depth + 1));
                }
            }
        }
    }
    result
}

/// Approximate arc length of a Bezier by summing `n` chord samples.
pub fn bezier_arc_length(poles: &[GpPnt], weights: Option<&[f64]>, n: usize) -> f64 {
    let n = n.max(2);
    let mut total = 0.0;
    let mut prev = match weights {
        Some(w) => rational_de_casteljau(poles, w, 0.0),
        None => de_casteljau(poles, 0.0),
    };
    for i in 1..=n {
        let t = i as f64 / n as f64;
        let cur = match weights {
            Some(w) => rational_de_casteljau(poles, w, t),
            None => de_casteljau(poles, t),
        };
        total += prev.distance(&cur);
        prev = cur;
    }
    total
}

fn push_dedup(result: &mut Vec<GpPnt>, p: GpPnt) {
    if result.last().map_or(true, |q| q.distance(&p) > 1e-12) {
        result.push(p);
    }
}

/// De Casteljau split at `t`: returns `(left_poles, right_poles)`.
fn split_bezier(poles: &[GpPnt], t: f64) -> (Vec<GpPnt>, Vec<GpPnt>) {
    let mut pts = poles.to_vec();
    let mut left = vec![pts[0]];
    let mut right = vec![*pts.last().unwrap()];
    while pts.len() > 1 {
        let mut next = Vec::with_capacity(pts.len() - 1);
        for i in 0..pts.len() - 1 {
            next.push(lerp(&pts[i], &pts[i + 1], t));
        }
        left.push(next[0]);
        right.push(*next.last().unwrap());
        pts = next;
    }
    right.reverse();
    (left, right)
}

/// Weight split that mirrors [`split_bezier`] (homogeneous coordinate part).
fn split_bezier_weights(weights: &[f64], t: f64) -> (Vec<f64>, Vec<f64>) {
    let mut ws = weights.to_vec();
    let mut left = vec![ws[0]];
    let mut right = vec![*ws.last().unwrap()];
    while ws.len() > 1 {
        let mut next = Vec::with_capacity(ws.len() - 1);
        for i in 0..ws.len() - 1 {
            next.push((1.0 - t) * ws[i] + t * ws[i + 1]);
        }
        left.push(next[0]);
        right.push(*next.last().unwrap());
        ws = next;
    }
    right.reverse();
    (left, right)
}

fn lerp(a: &GpPnt, b: &GpPnt, t: f64) -> GpPnt {
    GpPnt::new(
        a.x() + t * (b.x() - a.x()),
        a.y() + t * (b.y() - a.y()),
        a.z() + t * (b.z() - a.z()),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quad() -> Vec<GpPnt> {
        vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.5, 1.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ]
    }

    #[test]
    fn de_casteljau_midpoint() {
        let p = de_casteljau(&quad(), 0.5);
        assert!((p.x() - 0.5).abs() < 1e-12, "x={}", p.x());
        assert!((p.y() - 0.5).abs() < 1e-12, "y={}", p.y());
    }

    #[test]
    fn polyline_endpoints_and_midpoint() {
        let pl = bezier_to_polyline(&quad(), None, 1e-3);
        assert!(pl.first().unwrap().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
        assert!(pl.last().unwrap().distance(&GpPnt::new(1., 0., 0.)) < 1e-9);
        // The curve midpoint (0.5, 0.5) is a subdivision boundary and is
        // therefore present in the polyline.
        assert!(
            pl.iter().any(|p| p.distance(&GpPnt::new(0.5, 0.5, 0.)) < 1e-9),
            "polyline missing midpoint: {pl:?}"
        );
    }

    #[test]
    fn flatness_and_refinement() {
        assert!(bezier_flatness(&quad()) > 0.0);
        let coarse = bezier_to_polyline(&quad(), None, 0.1);
        let fine = bezier_to_polyline(&quad(), None, 1e-4);
        assert!(fine.len() >= coarse.len());
        // Converges to the true arc length as tol -> 0.
        let arc = bezier_arc_length(&quad(), None, 1000);
        let pl = bezier_to_polyline(&quad(), None, 1e-6);
        let pl_len: f64 = pl.windows(2).map(|w| w[0].distance(&w[1])).sum();
        assert!((pl_len - arc).abs() < 1e-4, "polyline {pl_len} vs arc {arc}");
    }

    #[test]
    fn rational_matches_polynomial() {
        let p = rational_de_casteljau(&quad(), &[1.0, 1.0, 1.0], 0.5);
        assert!((p.x() - 0.5).abs() < 1e-12);
        assert!((p.y() - 0.5).abs() < 1e-12);
        // Rational flattening with unit weights equals the polynomial case.
        let pl = bezier_to_polyline(&quad(), Some(&[1.0, 1.0, 1.0]), 1e-4);
        assert!(pl.first().unwrap().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
        assert!(pl.last().unwrap().distance(&GpPnt::new(1., 0., 0.)) < 1e-9);
    }
}
