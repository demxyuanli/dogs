//! Extrema between geometric objects — point/curve/surface minimum and
//! maximum distances.
//!
//! Facade over the Phase 13 analytic/Newton ports (source: `Extrema`,
//! TKGeomBase): exact + all-extrema implementations live in
//! `extrema_pc` (point–curve), `extrema_cc` (curve–curve),
//! `extrema_surf` (point/curve–surface), `extrema_ss` (surface–surface).
//! Public signatures are kept compatible with the original sampling API.

use occt_core::gp::{GpPnt, GpVec};

use crate::{Curve, Surface};

/// A solved extremum between two objects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtremaPair {
    pub p1: GpPnt,
    pub p2: GpPnt,
    pub distance: f64,
    /// Parameters on each object (curve u, or surface u/v).
    pub u1: f64,
    pub v1: Option<f64>,
    pub u2: f64,
    pub v2: Option<f64>,
}

/// Golden-section refine of a 1-D closest-point parameter around the initial
/// guess. Handles unbounded curves via window expansion. Used as the fallback
/// for degenerate curves by `extrema_pc`.
pub fn refine_curve_point(c: &dyn Curve, p: &GpPnt, a: f64, b: f64) -> (f64, GpPnt) {
    let f = |u: f64| {
        let q = c.d0(u);
        (q.x() - p.x()).powi(2) + (q.y() - p.y()).powi(2) + (q.z() - p.z()).powi(2)
    };
    let (mut lo, mut hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        // Unbounded: probe expanding windows until a strict minimum is
        // interior (both endpoints larger than the center).
        let mut w = 1.0;
        let (mut lo, mut hi) = (-w, w);
        let mut fm = f(0.0);
        for _ in 0..8 {
            let fl = f(lo);
            let fh = f(hi);
            if fl >= fm && fh >= fm {
                break; // minimum interior to [lo, hi]
            }
            w *= 8.0;
            lo = -w;
            hi = w;
            fm = f(0.0);
        }
        (lo, hi)
    };
    // Golden-section refinement (golden ratio conjugate 0.618...).
    let phi = (5.0f64.sqrt() - 1.0) * 0.5;
    for _ in 0..64 {
        let x1 = hi - phi * (hi - lo);
        let x2 = lo + phi * (hi - lo);
        if f(x1) < f(x2) {
            hi = x2;
        } else {
            lo = x1;
        }
    }
    let u = 0.5 * (lo + hi);
    (u, c.d0(u))
}

/// Minimum distance from point `p` to curve `c` (with the closest point).
/// Exact analytic dispatch for lines/circles; Newton refinement otherwise.
pub fn point_curve_extrema(c: &dyn Curve, p: &GpPnt) -> ExtremaPair {
    crate::extrema_pc::point_curve_extrema(c, p)
}

/// Maximum distance from point `p` to curve `c` over its parameter range.
/// `samples` is kept for signature compatibility but no longer drives the
/// computation.
pub fn point_curve_max_extrema(c: &dyn Curve, p: &GpPnt, _samples: usize) -> ExtremaPair {
    crate::extrema_pc::point_curve_max_extrema(c, p)
}

/// All local extrema between two curves, sorted by distance (min first).
/// `samples` is kept for signature compatibility.
pub fn curve_curve_extrema(c1: &dyn Curve, c2: &dyn Curve, _samples: usize) -> Vec<ExtremaPair> {
    crate::extrema_cc::curve_curve_extrema_all(c1, c2)
}

/// Minimum distance between a curve and a surface.
/// `samples` is kept for signature compatibility.
pub fn curve_surface_extrema(c: &dyn Curve, s: &dyn Surface, _samples: usize) -> ExtremaPair {
    crate::extrema_surf::curve_surface_extrema(c, s)
}

/// Minimum distance between two surfaces.
/// `samples` is kept for signature compatibility.
pub fn surface_surface_extrema(s1: &dyn Surface, s2: &dyn Surface, _samples: usize) -> ExtremaPair {
    crate::extrema_ss::surface_surface_extrema(s1, s2)
}

/// Minimum distance from point `p` to surface `s` (with the closest point).
pub fn point_surface_extrema(s: &dyn Surface, p: &GpPnt) -> ExtremaPair {
    crate::extrema_surf::point_surface_extrema(s, p)
}

/// Build an extrema pair from a point and a tangent vector on a curve —
/// convenience for callers that already have the closest parameter.
pub fn tangent_at(c: &dyn Curve, u: f64) -> GpVec {
    let (_, d) = c.d1(u);
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeomCircle, GeomLine, GeomPlane, GeomSphere};
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir, GpLin, GpPln, GpPnt, GpSphere as _GpSphere, GpTrsf, GpVec};

    fn line(p: GpPnt, d: GpDir) -> GeomLine {
        GeomLine::new(GpLin::new(occt_core::gp::GpAx1::new(p, d)))
    }

    fn sphere(r: f64) -> GeomSphere {
        GeomSphere::new(_GpSphere::new(GpAx3::standard(), r).unwrap())
    }

    #[test]
    fn point_line_min_distance() {
        let line = line(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e = point_curve_extrema(&line, &GpPnt::new(3.0, 4.0, 0.0));
        assert!((e.distance - 4.0).abs() < 1e-7, "dist {}", e.distance);
        assert!((e.p2.x() - 3.0).abs() < 1e-6, "closest x {}", e.p2.x());
    }

    #[test]
    fn point_circle_min_and_max() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let p = GpPnt::new(3.0, 0.0, 0.0);
        let e = point_curve_extrema(&circle, &p);
        assert!((e.distance - 2.0).abs() < 1e-6, "min dist {}", e.distance);
        let m = point_curve_max_extrema(&circle, &p, 64);
        assert!((m.distance - 4.0).abs() < 1e-6, "max dist {}", m.distance);
    }

    #[test]
    fn curve_curve_circles_min() {
        let c1 = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let mut c2 = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(3.0, 0.0, 0.0));
        c2.transform(&t);
        let es = curve_curve_extrema(&c1, &c2, 24);
        assert!(!es.is_empty());
        assert!((es[0].distance - 1.0).abs() < 1e-5, "min {}", es[0].distance);
    }

    #[test]
    fn curve_surface_sphere_line() {
        // Line at y=3, sphere radius 1 at origin → min distance 2.
        let line = line(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let sphere = sphere(1.0);
        let e = curve_surface_extrema(&line, &sphere, 16);
        assert!((e.distance - 2.0).abs() < 1e-4, "min {}", e.distance);
    }

    #[test]
    fn point_surface_sphere() {
        let sphere = sphere(1.0);
        let e = point_surface_extrema(&sphere, &GpPnt::new(3.0, 0.0, 0.0));
        assert!((e.distance - 2.0).abs() < 1e-5, "min {}", e.distance);
    }

    #[test]
    fn surface_surface_planes_parallel() {
        let ax = GpAx3::standard();
        let p1 = GeomPlane::new(GpPln::new(ax.clone()));
        let p2 = GeomPlane::new(GpPln::new(ax));
        let e = surface_surface_extrema(&p1, &p2, 3);
        assert!(e.distance.abs() < 1e-6, "coincident planes {}", e.distance);
    }

    #[test]
    fn surface_surface_sphere_plane() {
        let sphere = sphere(1.0);
        let ax = GpAx3::standard();
        let mut pl = GpPln::new(ax);
        pl.set_location(&GpPnt::new(0.0, 0.0, 5.0));
        let plane = GeomPlane::new(pl);
        let e = surface_surface_extrema(&sphere, &plane, 4);
        assert!((e.distance - 4.0).abs() < 1e-4, "min {}", e.distance);
    }

    #[test]
    fn curve_curve_skew_lines() {
        let l1 = line(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let l2 = line(GpPnt::new(0.0, 0.0, 3.0), GpDir::new(0.0, 1.0, 0.0).unwrap());
        let es = curve_curve_extrema(&l1, &l2, 8);
        assert!(!es.is_empty());
        assert!((es[0].distance - 3.0).abs() < 1e-5, "min {}", es[0].distance);
    }
}
