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
/// `None` mirrors OCCT's `IsDone() == false` (`Extrema_GGExtPC.hxx:531`).
///
/// `Extrema_GenExtPC` lives in `extrema_pc`; this re-export keeps the
/// `extrema::point_curve_extrema` path working.
pub use crate::extrema_pc::point_curve_extrema;

/// Maximum distance from point `p` to curve `c` over its parameter range.
/// `samples` is kept for signature compatibility but no longer drives the
/// computation. `None` as above.
pub fn point_curve_max_extrema(c: &dyn Curve, p: &GpPnt, _samples: usize) -> Option<ExtremaPair> {
    crate::extrema_pc::point_curve_max_extrema(c, p)
}

/// All local extrema between two curves, sorted by distance (min first).
/// `samples` is kept for signature compatibility.
pub fn curve_curve_extrema(c1: &dyn Curve, c2: &dyn Curve, _samples: usize) -> Vec<ExtremaPair> {
    crate::extrema_cc::curve_curve_extrema_all(c1, c2)
}

/// Minimum distance between a curve and a surface.
/// `samples` is kept for signature compatibility.
/// `None` = `Extrema_ExtCS::IsDone() == false` or `NbExt() == 0`.
pub fn curve_surface_extrema(c: &dyn Curve, s: &dyn Surface, _samples: usize) -> Option<ExtremaPair> {
    crate::extrema_surf::curve_surface_extrema(c, s)
}

/// Minimum distance between two surfaces — the nearest-pair answer of
/// `GeomAPI_ExtremaSurfaceSurface`.
///
/// Port of `GeomAPI_ExtremaSurfaceSurface::Init(S1, S2)`
/// (`GeomAPI_ExtremaSurfaceSurface.cxx:57-91`) plus `NearestPoints`
/// (`cxx:202-208`), `LowerDistanceParameters` (`cxx:210-221`) and
/// `LowerDistance` (`cxx:223-231`), which all read the `myIndex` minimum that
/// `Init` keeps (`cxx:75-89`) — not the first reported extremum. The engine is
/// `Extrema_ExtSS` over each surface's own range with
/// `Tol = Precision::PConfusion()` (`cxx:63-68`).
/// `samples` is kept for signature compatibility.
///
/// `None` mirrors the two OCCT exits that leave no pair to report: `!myIsDone`
/// (`cxx:68`), and the `Extrema_ExtElSS` parallel-plane arm, which records the
/// square distance but allocates no `myPOnS1` / `myPOnS2`
/// (`Extrema_ExtSS.cxx:228-235`), so `Points(myIndex, P1, P2)`
/// (`Extrema_ExtSS.cxx:292-300`) walks an empty sequence.
pub fn surface_surface_extrema(
    s1: &dyn Surface,
    s2: &dyn Surface,
    _samples: usize,
) -> Option<ExtremaPair> {
    let ex = crate::geom_api::ExtremaSurfaceSurface::new(s1, s2);
    if ex.nb_extrema() == 0 || ex.is_parallel() {
        return None;
    }
    let (p1, p2) = ex.nearest_points();
    let (u1, v1, u2, v2) = ex.lower_distance_parameters();
    Some(ExtremaPair {
        p1,
        p2,
        distance: ex.lower_distance(),
        u1,
        v1: Some(v1),
        u2,
        v2: Some(v2),
    })
}

/// Minimum distance from point `p` to surface `s` (with the closest point).
/// `None` = `Extrema_ExtPS::IsDone() == false` or `NbExt() == 0`.
pub use crate::extrema_surf::point_surface_extrema;

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
        let e = point_curve_extrema(&line, &GpPnt::new(3.0, 4.0, 0.0)).expect("line projection");
        assert!((e.distance - 4.0).abs() < 1e-7, "dist {}", e.distance);
        assert!((e.p2.x() - 3.0).abs() < 1e-6, "closest x {}", e.p2.x());
    }

    #[test]
    fn point_circle_min_and_max() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let p = GpPnt::new(3.0, 0.0, 0.0);
        let e = point_curve_extrema(&circle, &p).expect("circle projection");
        assert!((e.distance - 2.0).abs() < 1e-6, "min dist {}", e.distance);
        let m = point_curve_max_extrema(&circle, &p, 64).expect("circle max");
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
        let e = curve_surface_extrema(&line, &sphere, 16).expect("Extrema_ExtCS done");
        assert!((e.distance - 2.0).abs() < 1e-4, "min {}", e.distance);
    }

    #[test]
    fn point_surface_sphere() {
        let sphere = sphere(1.0);
        let e = point_surface_extrema(&sphere, &GpPnt::new(3.0, 0.0, 0.0)).expect("Extrema_ExtPS done");
        assert!((e.distance - 2.0).abs() < 1e-5, "min {}", e.distance);
    }

    #[test]
    fn surface_surface_planes_parallel() {
        let ax = GpAx3::standard();
        let p1 = GeomPlane::new(GpPln::new(ax.clone()));
        let p2 = GeomPlane::new(GpPln::new(ax));
        // Coincident planes are the parallel `Extrema_ExtElSS` arm: it stores
        // the square distance and no points (`Extrema_ExtElSS.cxx:46-53`), so
        // `surface_surface_extrema` has no pair to report
        // (`Extrema_ExtSS.cxx:229-233`) and the faithful class is read instead.
        assert!(surface_surface_extrema(&p1, &p2, 3).is_none());
        let ext = crate::extrema_ss::ExtremaExtSS::new(&p1, &p2, 0.0, 0.0);
        assert!(ext.is_done(), "Extrema_ExtSS done");
        assert!(ext.is_parallel(), "coincident planes are parallel");
        assert!(ext.square_distance(1).abs() < 1e-6, "coincident planes");
    }

    /// Sphere x an *unbounded* plane is not an analytic arm in OCCT: the pair
    /// falls into `Extrema_GenExtSS` with the plane window taken from
    /// `GeomAdaptor_Surface` (`GeomAdaptor_Surface.hxx:150-160` ->
    /// `Geom_Plane::Bounds` = `±Precision::Infinite()`), so OCCT 8.0.0 answers
    /// `sqDist = 7.992002e+200` for `sphere(r=1)` vs `plane z=5`
    /// (`specs\occt_probe\run_here.bat dummy --ssx`, `SSX sph-pln` and
    /// `SSX api-sph-pln dist=2.827012911183817e+100`). The 4.0 the deleted
    /// `plane_sphere_extrema` arm used to invent is not what
    /// `Extrema_ExtElSS::Perform(gp_Pln, gp_Sphere)` does -- it throws
    /// `Standard_NotImplemented` (`Extrema_ExtElSS.cxx:62-69`). The GenExtSS
    /// detail is pinned in `crate::extrema_ss::tests`.
    #[test]
    fn surface_surface_sphere_plane() {
        let sphere = sphere(1.0);
        let ax = GpAx3::standard();
        let mut pl = GpPln::new(ax);
        pl.set_location(&GpPnt::new(0.0, 0.0, 5.0));
        let plane = GeomPlane::new(pl);
        let e = surface_surface_extrema(&sphere, &plane, 4).expect("Extrema_ExtSS done");
        let expected = 7.992002e200f64.sqrt();
        assert!(
            (e.distance - expected).abs() <= 1e-12 * expected,
            "sphere vs unbounded plane {}",
            e.distance
        );
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
