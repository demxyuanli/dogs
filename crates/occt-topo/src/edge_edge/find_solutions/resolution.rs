use occt_core::precision::RESOLUTION;
use occt_geom::Curve;

use crate::inttools_data::IntRange;

use super::curve_kind::{geom_curve_kind, GeomCurveKind};

// ---------------------------------------------------------------------------
// Resolution helpers
// ---------------------------------------------------------------------------

/// `ResolutionCoeff` (`IntTools_EdgeEdge.cxx:1486-1559`).
pub(in crate::edge_edge) fn resolution_coeff(c: &dyn Curve, kind: GeomCurveKind, range: IntRange) -> f64 {
    match kind {
        GeomCurveKind::Circle => match c.circle_radius() {
            Some(r) => 1.0 / (2.0 * r),
            None => 0.0,
        },
        GeomCurveKind::Ellipse => match c.gp_ellipse() {
            Some(e) => 1.0 / e.major_radius(),
            None => 0.0,
        },
        GeomCurveKind::OffsetCurve => {
            // `cxx:1501-1520`: offset-of-line leaves the switch with 0 (its
            // `Resolution` arm returns `theR3D`), offset-of-circle/ellipse
            // return an analytic coefficient, every other basis falls through
            // into the sampling arm.
            let Some((basis, offset)) = c.offset_curve() else {
                return 0.0;
            };
            match geom_curve_kind(&*basis) {
                GeomCurveKind::Line => 0.0,
                GeomCurveKind::Circle => match basis.circle_radius() {
                    Some(r) => 1.0 / (2.0 * (offset + r)),
                    None => 0.0,
                },
                GeomCurveKind::Ellipse => match basis.gp_ellipse() {
                    Some(e) => 1.0 / (offset + e.major_radius()),
                    None => 0.0,
                },
                _ => sample_resolution_coeff(c, range),
            }
        }
        GeomCurveKind::Hyperbola | GeomCurveKind::Parabola | GeomCurveKind::OtherCurve => {
            sample_resolution_coeff(c, range)
        }
        _ => 0.0,
    }
}

/// The shared sampling arm of `ResolutionCoeff` (`cxx:1524-1551`): 30 equal
/// steps over the range, the smallest `dt / chord` ratio (or the initial `10.`
/// when every step is degenerate).
fn sample_resolution_coeff(c: &dyn Curve, range: IntRange) -> f64 {
    let nb_p = 30;
    let (t1, t2) = (range.first, range.last);
    let dt = (t2 - t1) / nb_p as f64;
    let mut t = t1;
    let mut k_min = 10.0;
    let mut p1 = c.d0(t1);
    for _ in 1..=nb_p {
        t += dt;
        let p2 = c.d0(t);
        let k = dt / p1.distance(&p2);
        if k < k_min {
            k_min = k;
        }
        p1 = p2;
    }
    k_min
}

/// `Resolution` (`IntTools_EdgeEdge.cxx:1561-1607`).
pub(in crate::edge_edge) fn resolution_of(c: &dyn Curve, kind: GeomCurveKind, res_coeff: f64, r3d: f64) -> f64 {
    match kind {
        GeomCurveKind::Line => r3d,
        GeomCurveKind::Circle => {
            let dt = res_coeff * r3d;
            if dt <= 1.0 {
                2.0 * dt.asin()
            } else {
                2.0 * std::f64::consts::PI
            }
        }
        // `Geom_BezierCurve::Resolution` / `Geom_BSplineCurve::Resolution`
        // (both `BSplCLib::Resolution`), reached through the `Curve` trait.
        GeomCurveKind::BezierCurve | GeomCurveKind::BSplineCurve => c.resolution(r3d),
        GeomCurveKind::OffsetCurve => {
            let Some((basis, _)) = c.offset_curve() else {
                return res_coeff * r3d;
            };
            match geom_curve_kind(&*basis) {
                GeomCurveKind::Line => r3d,
                GeomCurveKind::Circle => {
                    let dt = res_coeff * r3d;
                    if dt <= 1.0 {
                        2.0 * dt.asin()
                    } else {
                        2.0 * std::f64::consts::PI
                    }
                }
                _ => res_coeff * r3d,
            }
        }
        _ => res_coeff * r3d,
    }
}

/// `CurveDeflection` (`IntTools_EdgeEdge.cxx:1611-1638`): the total turn of the
/// tangent over 10 equal steps, skipping steps with a degenerate derivative.
pub(in crate::edge_edge) fn curve_deflection(c: &dyn Curve, range: IntRange) -> f64 {
    let nb_p = 10;
    let (t1, t2) = (range.first, range.last);
    let dt = (t2 - t1) / nb_p as f64;
    let mut t = t1;
    let mut defl = 0.0;
    let mut v1 = c.d1(t1).1;
    for _ in 1..=nb_p {
        t += dt;
        let v2 = c.d1(t).1;
        if v1.magnitude() > RESOLUTION && v2.magnitude() > RESOLUTION {
            defl += v1.angle(&v2);
        }
        v1 = v2;
    }
    defl
}

/// `IsClosed` (`IntTools_EdgeEdge.cxx:1642-1659`).
pub(in crate::edge_edge) fn is_closed(c: &dyn Curve, t1: f64, t2: f64, tol: f64, res: f64) -> bool {
    if (t1 - t2).abs() < res {
        return false;
    }
    let p1 = c.d0(t1);
    let p2 = c.d0(t2);
    p1.distance(&p2) < tol
}
