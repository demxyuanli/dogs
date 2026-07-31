//! B-spline -> Bezier decomposition (2D).
//! Source: `Geom2dConvert_BSplineCurveToBezierCurve.cxx` (and the 3D
//! `GeomConvert_BSplineCurveToBezierCurve.cxx`).

use crate::bspline_curve::Geom2dBSplineCurve;
use occt_core::bspl::knots;
use occt_core::gp::GpPnt2d;

/// One Bezier arc extracted from a B-spline.
#[derive(Clone, Debug)]
pub struct BezierSegment2d {
    pub poles: Vec<GpPnt2d>,
    pub weights: Option<Vec<f64>>,
    pub u0: f64,
    pub u1: f64,
}

/// Split a (non-rational) B-spline curve into adjacent Bezier arcs.
///
/// Algorithm (GeomConvert_BSplineCurveToBezierCurve): raise every interior
/// knot to full multiplicity (degree+1) by Boehm knot insertion, then each
/// non-empty knot span carries exactly `degree+1` control points forming a
/// Bezier arc. The result's arc `k` evaluates on `[u_k, u_{k+1}]`.
pub fn split_bspline2d_to_beziers(curve: &Geom2dBSplineCurve) -> Vec<BezierSegment2d> {
    let deg = curve.degree;
    let np = curve.nb_poles();
    if np == 0 || curve.knots.len() < 2 {
        return Vec::new();
    }

    let mut poles: Vec<GpPnt2d> = (0..np)
        .map(|i| GpPnt2d::new(curve.xs[i], curve.ys[i]))
        .collect();
    let mut knots = curve.knots.clone();

    // Distinct interior knots (strictly between first and last parameter).
    let first = knots[deg];
    let last = knots[knots.len() - 1 - deg];
    let mut interior: Vec<f64> = Vec::new();
    for &k in &knots {
        if k > first + 1e-15 && k < last - 1e-15 {
            if interior.last().map_or(true, |&p| (p - k).abs() > 1e-15) {
                interior.push(k);
            }
        }
    }

    // Insert each interior knot until it reaches multiplicity degree+1.
    for &u in &interior {
        let mult = knots::multiplicity(&knots, u);
        for _ in 0..(deg + 1).saturating_sub(mult) {
            boehm_insert_2d(&mut poles, &knots, u, deg);
            knots = knots::insert_knot(&knots, u, 1);
        }
    }

    // Distinct knots: each consecutive pair bounds one Bezier arc.
    let mut dk: Vec<f64> = Vec::new();
    for &k in &knots {
        if dk.last().map_or(true, |&p| (p - k).abs() > 1e-15) {
            dk.push(k);
        }
    }

    let mut segs = Vec::new();
    for i in 0..dk.len().saturating_sub(1) {
        let u0 = dk[i];
        let u1 = dk[i + 1];
        if u1 - u0 <= 1e-15 {
            continue;
        }
        let start = i * (deg + 1);
        let end = start + deg;
        if end >= poles.len() {
            break;
        }
        segs.push(BezierSegment2d {
            poles: poles[start..=end].to_vec(),
            weights: None,
            u0,
            u1,
        });
    }
    segs
}

/// Number of Bezier arcs `split_bspline2d_to_beziers` would produce
/// (= number of distinct knots minus one).
pub fn bspline2d_bezier_count(curve: &Geom2dBSplineCurve) -> usize {
    let mut count = 0usize;
    let mut last: Option<f64> = None;
    for &k in &curve.knots {
        if last.map_or(true, |p| (p - k).abs() > 1e-15) {
            count += 1;
            last = Some(k);
        }
    }
    count.saturating_sub(1)
}

/// 2D Boehm knot insertion (mirrors `occt_core::bspl::bezier::boehm_insert`
/// for `GpPnt2d`).
fn boehm_insert_2d(poles: &mut Vec<GpPnt2d>, knots: &[f64], u: f64, degree: usize) {
    let n = poles.len();
    if n == 0 {
        return;
    }
    let k = knots::hunt(knots, u).min(n - 1);
    let m = knots::multiplicity(knots, u);
    let p = degree;
    let mut new_poles = Vec::with_capacity(n + 1);
    for i in 0..=n {
        if i <= k.saturating_sub(p) {
            new_poles.push(poles[i]);
        } else if i <= k.saturating_sub(m) {
            let alpha = knot_alpha(knots, u, i, p);
            new_poles.push(lerp2(poles[i - 1], poles[i], alpha));
        } else {
            new_poles.push(poles[i - 1]);
        }
    }
    *poles = new_poles;
}

/// Interpolation factor `a_i = (u - U_i) / (U_{i+p} - U_i)`, 0 when degenerate.
fn knot_alpha(knots: &[f64], u: f64, i: usize, p: usize) -> f64 {
    if i + p < knots.len() && (knots[i + p] - knots[i]).abs() > 1e-30 {
        ((u - knots[i]) / (knots[i + p] - knots[i])).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

/// Linear interpolation `a + t * (b - a)`.
fn lerp2(a: GpPnt2d, b: GpPnt2d, t: f64) -> GpPnt2d {
    GpPnt2d::new(a.x() + t * (b.x() - a.x()), a.y() + t * (b.y() - a.y()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_splits_into_one_segment() {
        let c = Geom2dBSplineCurve::new(
            vec![0.0, 1.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0, 1.0, 1.0],
            1,
        )
        .unwrap();
        let segs = split_bspline2d_to_beziers(&c);
        assert_eq!(segs.len(), 1);
        assert_eq!(segs[0].poles.len(), 2);
        assert!((segs[0].poles[0].x() - 0.0).abs() < 1e-12);
        assert!((segs[0].poles[1].x() - 1.0).abs() < 1e-12);
        assert!((segs[0].u0 - 0.0).abs() < 1e-12);
        assert!((segs[0].u1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn degree2_three_spans() {
        let c = Geom2dBSplineCurve::new(
            vec![0.0, 0.5, 1.0, 1.5, 2.0],
            vec![0.0, 1.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0],
            2,
        )
        .unwrap();
        let segs = split_bspline2d_to_beziers(&c);
        assert_eq!(segs.len(), 3, "segments={segs:?}");
        assert!((segs[0].u0 - 0.0).abs() < 1e-12 && (segs[0].u1 - 1.0).abs() < 1e-12);
        assert!((segs[1].u0 - 1.0).abs() < 1e-12 && (segs[1].u1 - 2.0).abs() < 1e-12);
        assert!((segs[2].u0 - 2.0).abs() < 1e-12 && (segs[2].u1 - 3.0).abs() < 1e-12);
        assert_eq!(bspline2d_bezier_count(&c), 3);
    }

    #[test]
    fn segments_join_continuously() {
        let c = Geom2dBSplineCurve::new(
            vec![0.0, 0.5, 1.0, 1.5, 2.0],
            vec![0.0, 1.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0],
            2,
        )
        .unwrap();
        let segs = split_bspline2d_to_beziers(&c);
        for k in 0..segs.len() - 1 {
            let a = *segs[k].poles.last().unwrap();
            let b = segs[k + 1].poles[0];
            assert!(a.distance(&b) < 1e-9, "junction {k}: {a:?} vs {b:?}");
        }
    }
}
