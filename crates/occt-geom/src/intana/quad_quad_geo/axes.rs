//! Axis and plane helpers of `IntAna_QuadQuadGeo` (`AxeOperator`,
//! `AxeOperator::Distance`, `gce_MakePln` and the circle/plane builders).

use super::super::prelude::*;
use super::*;

/// Relation between two axes — port of `AxeOperator`
/// (IntAna_QuadQuadGeo.cxx). Returns `(parallel, distance, coplanar,
/// intersection point)`; `coplanar` needs the axes to (nearly) meet
/// (`distance < eps_dist` and the triple product within tolerance), and the
/// point is only computed for concurrent non-parallel axes.
pub(super) fn axe_operator(a1: &GpAx1, a2: &GpAx1, eps_dist: f64, eps_para: f64) -> (bool, f64, bool, Option<GpPnt>) {
    let mut v1 = *a1.direction();
    let mut v2 = *a2.direction();
    refine_dir(&mut v1);
    refine_dir(&mut v2);
    let (p1, p2) = (*a1.location(), *a2.location());
    let (w1, w2) = (GpVec::from_xyz(v1.xyz()), GpVec::from_xyz(v2.xyz()));
    let parallel = w1.cross_magnitude(&w2) <= eps_para;
    let distance = if parallel {
        let rel = GpVec::from_pnts(&p1, &p2);
        rel.subtracted(&w1.multiplied_scalar(rel.dot(&w1))).magnitude()
    } else {
        w1.crossed(&w2).normalized().dot(&GpVec::from_pnts(&p1, &p2)).abs()
    };
    let mut coplanar = false;
    let mut pt_intersect = None;
    if distance < eps_dist {
        let det = det33(&w1, &w2, &GpVec::from_pnts(&p2, &p1)); // rows V1, V2, P1−P2
        if det.abs() <= eps_dist {
            coplanar = true;
            if !parallel {
                // Concurrent axes: intersection point P1 + A·V1.
                let sm = GpVec::from_pnts(&p1, &p2);
                let d1 = w1.y() * w2.x() - w1.x() * w2.y();
                let d2 = w1.z() * w2.y() - w1.y() * w2.z();
                let d3 = w1.z() * w2.x() - w1.x() * w2.z();
                let a = if d1 != 0.0 && d1.abs() >= d2.abs() && d1.abs() >= d3.abs() {
                    (sm.y() * w2.x() - sm.x() * w2.y()) / d1
                } else if d2 != 0.0 && d2.abs() >= d1.abs() && d2.abs() >= d3.abs() {
                    (sm.z() * w2.y() - sm.y() * w2.z()) / d2
                } else {
                    (sm.z() * w2.x() - sm.x() * w2.z()) / d3
                };
                pt_intersect = Some(p1.translated_vec(&w1.multiplied_scalar(a)));
            }
        }
    }
    (parallel, distance, coplanar, pt_intersect)
}

/// The perpendicular segment between two non-parallel axes: signed distance and
/// the parameters of the closest points on each axis. Port of
/// `AxeOperator::Distance`.
pub(super) fn axe_distance(a1: &GpAx1, a2: &GpAx1) -> (f64, f64, f64) {
    let w1 = GpVec::from_xyz(a1.direction().xyz());
    let w2 = GpVec::from_xyz(a2.direction().xyz());
    let o1o2 = GpVec::from_pnts(a1.location(), a2.location());
    let n = w1.crossed(&w2);
    if n.magnitude() < 1e-12 {
        return (0.0, 0.0, 0.0);
    }
    let n = n.normalized();
    let d = det33(&w1, &w2, &n);
    if d != 0.0 {
        let dist = det33(&w1, &w2, &o1o2) / d;
        let p1 = det33(&o1o2, &w2, &n) / (-d);
        let p2 = det33(&w1, &o1o2, &n) / d;
        (dist, p1, p2)
    } else {
        (0.0, 0.0, 0.0)
    }
}

/// `(A1A2.Intersect() && Pt.Distance(A1A2.PtIntersect()) == 0.0) || A1A2.Same()`
/// (`IntAna_QuadQuadGeo.cxx:1378-1379` for cylinder/sphere, `cxx:1924-1925` for
/// sphere/cone): the sphere center must lie exactly on the point where the two
/// axes meet, or the axes must be parallel and coincident within
/// `myEPSILON_DISTANCE`. `AxeOperator`'s defaults are `theEpsDistance = 1.e-14`
/// and `theEpsAxesPara = Precision::Angular()` (`cxx:72-73`). OCCT adds no
/// tolerance of its own, so neither does this predicate.
pub fn axes_meet_at_point(a1: &GpAx1, a2: &GpAx1, pt: &GpPnt) -> bool {
    let (parallel, distance, coplanar, pt_intersect) = axe_operator(a1, a2, 1.0e-14, ANGULAR);
    // `A1A2.Intersect()` is `thecoplanar && !theparallel` (`cxx:86`).
    if coplanar && !parallel {
        if let Some(pi) = pt_intersect {
            if pt.distance(&pi) == 0.0 {
                return true;
            }
        }
    }
    // `A1A2.Same()` is `theparallel && (thedistance < myEPSILON_DISTANCE)` (`cxx:83`).
    parallel && distance < 1.0e-14
}

pub(super) fn mid_pnt(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// A circle with the given center / plane normal / radius (the frame X
/// direction is arbitrary — it does not affect the curve).
pub(super) fn circle_with_normal(center: GpPnt, normal: GpDir, radius: f64) -> GpCirc {
    GpCirc::new(ax2_from_dirs(center, normal, perp_x_dir(&normal)), radius)
}

/// A plane through `pt` with the given normal.
pub(super) fn plane_normal_at(pt: GpPnt, normal: GpDir) -> Option<GpPln> {
    let ax3 = GpAx3::new(pt, normal, &perp_x_dir(&normal)).ok()?;
    Some(GpPln::new(ax3))
}

/// `gce_MakePln(P1, P2, P3)` (`gce_MakePln.cxx:76-92`): the plane through three
/// points, X direction along `P1P2`. `None` for collinear points
/// (`gce_ColinearPoints`).
pub fn plane_through_three(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> Option<GpPln> {
    let v1 = GpVec::from_pnts(p1, p2);
    let v2 = GpVec::from_pnts(p1, p3);
    let normal = v1.crossed(&v2);
    if normal.magnitude() < REAL_SMALL {
        return None;
    }
    let n = GpDir::from_vec(&normal).ok()?;
    let x = GpDir::from_vec(&v1).ok()?;
    let ax3 = GpAx3::new(*p1, n, &x).ok()?;
    Some(GpPln::new(ax3))
}
