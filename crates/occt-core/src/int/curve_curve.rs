//! Curve-curve intersection primitives (port-internal).
//!
//! **Provenance (audit A9)**: **not** a translation of
//! `IntCurveCurve_IntImpCurveCurve` / `IntTools`. The functions below are
//! closest-point / overlap primitives in the style of Ericson, *Real-Time
//! Collision Detection*; OCCT intersects curves through `IntImp_CurveCurve` +
//! `IntCurve_*` (with `Extrema_ExtCC` for proximity) and edge–edge proximity
//! through `IntTools_EdgeEdge`.
use crate::gp::{GpLin, GpPnt, GpVec};

const EPS: f64 = 1e-14;

/// 3D segment-segment closest-point + coincidence test within `Precision::Confusion`.
/// Returns the intersection point (or the midpoint of the closest pair of the
/// shortest transversal) when the segments are within tolerance; for collinear
/// overlapping segments returns the midpoint of the overlap region.
pub fn segment_segment_intersection_3d(a: &GpPnt, b: &GpPnt, c: &GpPnt, d: &GpPnt) -> Option<GpPnt> {
    let tol = crate::precision::CONFUSION;
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(c, d);
    let l1 = d1.magnitude();
    let l2 = d2.magnitude();
    if l1 <= tol && l2 <= tol {
        return (a.distance(c) <= tol).then_some(*a);
    }
    if l1 <= tol {
        let (pc, _) = point_segment_closest(a, c, d);
        return (a.distance(&pc) <= tol).then_some(*a);
    }
    if l2 <= tol {
        let (pc, _) = point_segment_closest(c, a, b);
        return (c.distance(&pc) <= tol).then_some(*c);
    }
    let cross = d1.crossed(&d2);
    if cross.square_magnitude() <= tol * tol * l1 * l1 * l2 * l2 {
        // Parallel: collinear overlap test.
        let r = GpVec::from_pnts(a, c);
        let dist = r.crossed(&d2).magnitude() / l2;
        if dist > tol { return None; }
        let u = d1.divided(l1);
        let tc = GpVec::from_pnts(a, c).dot(&u);
        let td = GpVec::from_pnts(a, d).dot(&u);
        let lo = tc.min(td).max(0.0);
        let hi = tc.max(td).min(l1);
        if lo > hi + tol { return None; }
        return Some(a.translated_vec(&u.multiplied_scalar((lo + hi) * 0.5)));
    }
    // Non-parallel: closest points on the two segments.
    let (p1, p2) = segment_segment_closest(a, b, c, d);
    if p1.distance(&p2) <= tol {
        Some(midpoint(&p1, &p2))
    } else {
        None
    }
}

/// Exact minimum distance between two 3D segments.
pub fn segment_segment_distance_3d(a: &GpPnt, b: &GpPnt, c: &GpPnt, d: &GpPnt) -> f64 {
    let (p1, p2) = segment_segment_closest(a, b, c, d);
    p1.distance(&p2)
}

/// Intersection of segment [a,b] with an infinite line. None if parallel or outside.
pub fn segment_line_intersection_3d(a: &GpPnt, b: &GpPnt, l: &GpLin) -> Option<GpPnt> {
    let d = GpVec::from_pnts(a, b);
    let u = GpVec::from_xyz(l.direction().xyz());
    let o = l.location();
    let r = GpVec::from_pnts(&o, a);
    let du = d.dot(&u);
    let uu = u.dot(&u);
    let cross = d.crossed(&u);
    let denom = cross.square_magnitude();
    if denom <= crate::precision::RESOLUTION { return None; } // parallel
    let rd = r.dot(&d);
    let ru = r.dot(&u);
    let s = (du * ru - rd * uu) / denom;
    if s < 0.0 || s > 1.0 { return None; }
    let pt = a.translated_vec(&d.multiplied_scalar(s));
    (l.distance(&pt) <= crate::precision::CONFUSION).then_some(pt)
}

/// Closest point on segment [a,b] to point `p`, plus its parameter t ∈ [0,1].
pub fn point_segment_closest(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> (GpPnt, f64) {
    let d = GpVec::from_pnts(a, b);
    let dd = d.square_magnitude();
    if dd <= crate::precision::RESOLUTION {
        return (*a, 0.0);
    }
    let t = (GpVec::from_pnts(a, p).dot(&d) / dd).clamp(0.0, 1.0);
    (a.translated_vec(&d.multiplied_scalar(t)), t)
}

/// Shortest distance from point `p` to infinite line `l`.
pub fn point_line_distance(p: &GpPnt, l: &GpLin) -> f64 {
    l.distance(p)
}

// --- internals ---

/// Closest points on two 3D segments (Ericson, Real-Time Collision Detection).
fn segment_segment_closest(a: &GpPnt, b: &GpPnt, c: &GpPnt, d: &GpPnt) -> (GpPnt, GpPnt) {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(c, d);
    let r = GpVec::from_pnts(c, a); // a - c
    let aa = d1.square_magnitude();
    let e = d2.square_magnitude();
    let f = d2.dot(&r);
    if aa <= EPS && e <= EPS {
        return (*a, *c);
    }
    let mut s = 0.0;
    let mut t = 0.0;
    if aa <= EPS {
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c1 = d1.dot(&r);
        if e <= EPS {
            s = (-c1 / aa).clamp(0.0, 1.0);
        } else {
            let bb = d1.dot(&d2);
            let denom = aa * e - bb * bb;
            if denom > EPS {
                s = ((bb * f - c1 * e) / denom).clamp(0.0, 1.0);
            }
            t = (bb * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c1 / aa).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((bb - c1) / aa).clamp(0.0, 1.0);
            }
        }
    }
    (
        a.translated_vec(&d1.multiplied_scalar(s)),
        c.translated_vec(&d2.multiplied_scalar(t)),
    )
}

fn midpoint(p: &GpPnt, q: &GpPnt) -> GpPnt {
    GpPnt::new((p.x() + q.x()) * 0.5, (p.y() + q.y()) * 0.5, (p.z() + q.z()) * 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::GpDir;

    #[test]
    fn crossing_segments_3d() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(2.0, 0.0, 0.0);
        let c = GpPnt::new(1.0, -1.0, 0.0);
        let d = GpPnt::new(1.0, 1.0, 0.0);
        let p = segment_segment_intersection_3d(&a, &b, &c, &d).unwrap();
        assert!(p.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
        assert!(segment_segment_distance_3d(&a, &b, &c, &d) < 1e-9);
    }

    #[test]
    fn skew_segments() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(2.0, 0.0, 0.0);
        let c = GpPnt::new(1.0, 0.0, 1.0);
        let d = GpPnt::new(1.0, 2.0, 1.0);
        let dist = segment_segment_distance_3d(&a, &b, &c, &d);
        assert!((dist - 1.0).abs() < 1e-9);
        assert!(dist > 0.0);
        assert!(segment_segment_intersection_3d(&a, &b, &c, &d).is_none());
    }

    #[test]
    fn collinear_overlap() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(2.0, 0.0, 0.0);
        let c = GpPnt::new(1.0, 0.0, 0.0);
        let d = GpPnt::new(3.0, 0.0, 0.0);
        let p = segment_segment_intersection_3d(&a, &b, &c, &d).unwrap();
        assert!(p.distance(&GpPnt::new(1.5, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn segment_line_cross() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(2.0, 0.0, 0.0);
        let l = GpLin::from_pnt_dir(GpPnt::new(1.0, -1.0, 0.0), GpDir::new(0.0, 1.0, 0.0).unwrap());
        let p = segment_line_intersection_3d(&a, &b, &l).unwrap();
        assert!(p.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
        assert!(point_line_distance(&GpPnt::new(1.0, 5.0, 0.0), &l) < 1e-12);
        assert!((point_line_distance(&GpPnt::new(3.0, 0.0, 0.0), &l) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn point_segment_closest_param() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(4.0, 0.0, 0.0);
        let p = GpPnt::new(2.0, 3.0, 0.0);
        let (cl, t) = point_segment_closest(&p, &a, &b);
        assert!(t > 0.49 && t < 0.51);
        assert!(cl.distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-12);
    }
}
