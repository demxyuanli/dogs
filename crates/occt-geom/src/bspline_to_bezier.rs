//! Split a B-spline curve into piecewise Bezier segments.
//! Source: `GeomConvert_SplitBSplineCurve.hxx`

use crate::bspline_curve::GeomBSplineCurve;
use occt_core::bspl::bezier::boehm_insert;
use occt_core::bspl::knots::{hunt, insert_knot, multiplicity};
use occt_core::gp::GpPnt;

/// A single Bezier segment extracted from a B-spline curve, defined over the
/// parameter interval `[u0, u1]`.
#[derive(Debug, Clone)]
pub struct BezierSegment {
    pub poles: Vec<GpPnt>,
    pub weights: Option<Vec<f64>>,
    pub u0: f64,
    pub u1: f64,
}

/// Number of Bezier segments the curve splits into — the number of distinct
/// knot values inside the parameter range, minus one (i.e. the number of
/// non-degenerate knot spans).
pub fn bezier_segment_count(curve: &GeomBSplineCurve) -> usize {
    let nk = curve.knots.len();
    let degree = curve.degree;
    if nk < 2 * degree + 2 || curve.poles.len() < degree + 1 {
        return 0;
    }
    let first = curve.knots[degree];
    let last = curve.knots[nk - 1 - degree];
    let mut distinct = 0usize;
    let mut prev: Option<f64> = None;
    for &k in &curve.knots {
        if k < first - 1e-15 || k > last + 1e-15 {
            continue;
        }
        if prev.map_or(true, |p| (k - p).abs() > 1e-15) {
            distinct += 1;
            prev = Some(k);
        }
    }
    distinct.saturating_sub(1)
}

/// Split `curve` into Bezier segments by Boehm knot insertion.
///
/// Mirrors OCCT's `GeomConvert_BSplineCurveToBezierCurve`: every interior knot
/// is inserted until it has multiplicity `degree`. At that point each distinct
/// knot span is a Bezier segment and consecutive segments share one control
/// point — segment `s` (0-based) is `poles[s*degree ..= s*degree + degree]`
/// over `[breakpoints[s], breakpoints[s+1]]`.
///
/// The curve is expected to be clamped (end knots at multiplicity `degree+1`),
/// which is the representation `GeomBSplineCurve::new`/`rational` produce.
pub fn split_bspline_to_beziers(curve: &GeomBSplineCurve) -> Vec<BezierSegment> {
    let degree = curve.degree;
    let nk = curve.knots.len();
    let target = degree + 1;
    if nk < 2 * target || curve.poles.len() < target {
        return Vec::new();
    }
    let first = curve.knots[degree];
    let last = curve.knots[nk - 1 - degree];

    let mut poles = curve.poles.clone();
    let mut weights = curve.weights.clone();
    let mut knots = curve.knots.clone();

    // Raise interior knots to multiplicity `degree`.
    let mut i = degree;
    while i < knots.len() - degree - 1 {
        let u = knots[i];
        if u <= first || u >= last {
            i += 1;
            continue;
        }
        let mut m = multiplicity(&knots, u);
        while m < degree {
            let idx = hunt(&knots, u);
            boehm_insert(&mut poles, &knots, idx, u, degree, weights.as_mut());
            knots = insert_knot(&knots, u, 1);
            m += 1;
        }
        i += m;
    }

    // Distinct knot values inside the parameter range become the segment
    // breakpoints.
    let mut breakpoints: Vec<f64> = Vec::new();
    for &k in &knots {
        if k < first - 1e-15 || k > last + 1e-15 {
            continue;
        }
        if breakpoints.last().map_or(true, |&b| (b - k).abs() > 1e-15) {
            breakpoints.push(k);
        }
    }

    let mut segments = Vec::with_capacity(breakpoints.len().saturating_sub(1));
    for s in 0..breakpoints.len().saturating_sub(1) {
        let u0 = breakpoints[s];
        let u1 = breakpoints[s + 1];
        if u1 <= u0 {
            continue;
        }
        let start = (s * degree).min(poles.len().saturating_sub(target));
        let seg_poles = poles[start..start + target].to_vec();
        let seg_weights = weights.as_ref().map(|w| w[start..start + target].to_vec());
        segments.push(BezierSegment { poles: seg_poles, weights: seg_weights, u0, u1 });
    }
    segments
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bezier_to_polyline::rational_de_casteljau;
    use crate::curve::Curve;

    fn approx(a: &GpPnt, b: &GpPnt, tol: f64) -> bool {
        a.distance(b) < tol
    }

    #[test]
    fn quadratic_splits_into_four_segments() {
        let poles = vec![
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(2., 0., 0.),
            GpPnt::new(3., -1., 0.),
            GpPnt::new(4., 0., 0.),
            GpPnt::new(5., 1., 0.),
        ];
        let knots = vec![0., 0., 0., 1., 2., 3., 4., 4., 4.];
        let c = GeomBSplineCurve::new(poles, knots, 2).unwrap();

        assert_eq!(bezier_segment_count(&c), 4);
        let segs = split_bspline_to_beziers(&c);
        assert_eq!(segs.len(), 4);

        // The first/last pole of each segment matches the curve at the
        // segment boundaries (continuity of the split).
        for s in &segs {
            assert_eq!(s.poles.len(), 3);
            assert!(s.weights.is_none());
            assert!(s.u1 > s.u0);
            assert!(approx(&s.poles[0], &c.d0(s.u0), 1e-9));
            assert!(approx(&s.poles[2], &c.d0(s.u1), 1e-9));
        }
        // Clamped B-spline endpoints equal the extreme poles.
        assert!(approx(&segs[0].poles[0], &c.d0(0.0), 1e-9));
        assert!(approx(&segs[3].poles[2], &c.d0(4.0), 1e-9));
    }

    #[test]
    fn linear_spline_single_segment() {
        let poles = vec![GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.)];
        let knots = vec![0., 0., 1., 1.];
        let c = GeomBSplineCurve::new(poles.clone(), knots, 1).unwrap();
        let segs = split_bspline_to_beziers(&c);
        assert_eq!(segs.len(), 1);
        assert_eq!(segs[0].poles, poles);
        assert_eq!(segs[0].u0, 0.0);
        assert_eq!(segs[0].u1, 1.0);
    }

    #[test]
    fn rational_circle_splits_and_matches() {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let poles = vec![
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
            GpPnt::new(-1., 1., 0.),
            GpPnt::new(-1., 0., 0.),
            GpPnt::new(-1., -1., 0.),
            GpPnt::new(0., -1., 0.),
            GpPnt::new(1., -1., 0.),
            GpPnt::new(1., 0., 0.),
        ];
        let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
        let pi = std::f64::consts::PI;
        let knots = vec![
            0., 0., 0.,
            pi / 2., pi / 2.,
            pi, pi,
            3. * pi / 2., 3. * pi / 2.,
            2. * pi, 2. * pi, 2. * pi,
        ];
        let c = GeomBSplineCurve::rational(poles, weights, knots, 2).unwrap();

        let segs = split_bspline_to_beziers(&c);
        assert_eq!(segs.len(), 4);

        let bps = [0.0, pi / 2.0, pi, 3.0 * pi / 2.0, 2.0 * pi];
        for (s, seg) in segs.iter().enumerate() {
            let (u0, u1) = (bps[s], bps[s + 1]);
            assert_eq!(seg.poles.len(), 3);
            assert_eq!(seg.weights.as_ref().unwrap().len(), 3);
            // Segment endpoints match the original B-spline at breakpoints.
            assert!(approx(&seg.poles[0], &c.d0(u0), 1e-9), "s={s} start");
            assert!(approx(&seg.poles[2], &c.d0(u1), 1e-9), "s={s} end");
            // Midpoint lies on the unit circle. (occt-core's rational B-spline
            // evaluator interpolates unweighted coordinates and is incorrect at
            // interior parameters, so compare against the analytic circle.)
            let tm = 0.5 * (u0 + u1);
            let expected = GpPnt::new(tm.cos(), tm.sin(), 0.0);
            let p = rational_de_casteljau(&seg.poles, seg.weights.as_ref().unwrap(), 0.5);
            assert!(approx(&p, &expected, 1e-9), "s={s} mid {:?} vs {expected:?}", p);
        }
    }
}
