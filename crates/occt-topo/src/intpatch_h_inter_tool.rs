//! `IntPatch_HInterTool` - the sampling / restriction-arc services the
//! `IntPatch` polyhedron and walking codes call.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_HInterTool.cxx`
//! and `.hxx` / `.lxx`.
//!
//! Ported here (OCCT 8.0.0):
//! * `SingularOnUMin` / `SingularOnUMax` / `SingularOnVMin` / `SingularOnVMax`
//!   (`lxx:18-37`): all four answer `false` in OCCT regardless of the surface.
//! * `NbSamplesU(S, u1, u2)` (`cxx:76-112`) and `NbSamplesV(S, v1, v2)`
//!   (`cxx:37-73`). `u1` / `u2` are unused by OCCT's bodies and are kept as
//!   parameters for signature parity only.
//! * `NbSamplesOnArc(A)` (`cxx:232-261`).
//! * `Bounds(C, Ufirst, Ulast)` (`cxx:263-269`).
//! * `Project(C, P, Paramproj, Ptproj)` (`cxx:271-304`) over the port's
//!   `Extrema_EPCOfExtPC2d`.
//! * `Tolerance(V, C)` (`cxx:306-310`) and `Parameter(V, C)` (`cxx:312-316`):
//!   the two forwards to `Adaptor3d_HVertex` (`BRepTopAdaptor_HVertex`).
//!
//! UNPORTED in this module:
//! * `NbSamplePoints` (`cxx:114-162`) and `SamplePoint` (`cxx:164-230`). They
//!   keep the instance state `uinf` / `vinf` / `usup` / `vsup` (`hxx:150-154`)
//!   and are only called by `IntPatch_TheSearchInside` / `IntPatch_PrmPrmIntersection`,
//!   neither of which consumes them yet; their port lands with those callers
//!   instead of as dead code here.
//! * the `HasBeenSeen` / `NbPoints` / `Value` / `IsVertex` / `Vertex` /
//!   `NbSegments` / `HasFirstPoint` / `HasLastPoint` / `IsAllSolution`
//!   (`cxx:318-368`) stubs: each answers a constant (`false` / `0` /
//!   `Standard_OutOfRange`) in OCCT and has no caller in the port.

use occt_core::gp::GpPnt2d;
use occt_core::precision::Precision;
use occt_geom2d::extrema2d::EPcOfExtPc2d;
use occt_geom2d::Curve2d;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::geom_abs_curve_type::{curve2d_kind, GeomAbsCurveType};
use occt_geom::Surface;

/// `NBMAXUV` cap applied by `IntPatch_Polyhedron::NbPOnU/NbPOnV`
/// (`IntPatch_Polyhedron.cxx:34`).
pub const NB_MAX_UV: i32 = 30;

/// `SingularOnUMin(S)` (`lxx:18-21`): `false`.
pub fn singular_on_u_min(_s: &dyn Surface) -> bool {
    false
}

/// `SingularOnUMax(S)` (`lxx:23-26`): `false`.
pub fn singular_on_u_max(_s: &dyn Surface) -> bool {
    false
}

/// `SingularOnVMin(S)` (`lxx:28-31`): `false`.
pub fn singular_on_v_min(_s: &dyn Surface) -> bool {
    false
}

/// `SingularOnVMax(S)` (`lxx:33-36`): `false`.
pub fn singular_on_v_max(_s: &dyn Surface) -> bool {
    false
}

/// `NbSamplesV(S, v1, v2)` (`cxx:37-73`).
pub fn nb_samples_v(s: &dyn Surface) -> i32 {
    if s.is_bezier_surface() {
        // `case GeomAbs_BezierSurface: return (3 + S->NbVPoles());`
        return 3 + s.nb_v_poles();
    }
    if s.is_bspline_surface() {
        // `case GeomAbs_BSplineSurface:` (`cxx:46-57`)
        let mut nbs = s.bspline_surface_vknots().map(|k| k.len() as i32).unwrap_or(0)
            * s.v_degree();
        if !bspline_v_rational(s) {
            nbs *= 2;
        }
        if nbs < 4 {
            nbs = 4;
        }
        return nbs;
    }
    // `GeomAbs_SurfaceOfRevolution` / `GeomAbs_SurfaceOfExtrusion` are not
    // named by `classify_surface`, but OCCT's V switch groups them with the
    // quadrics (15, `cxx:59-65`), while U groups them with offset / other
    // (10, `cxx:100-107`).
    if is_swept_surface(s) {
        return 15;
    }
    match classify_surface(s) {
        // `case GeomAbs_Plane: return 2;` (`cxx:41-42`)
        SurfaceKind::Plane => 2,
        // `cxx:59-65`: cylinder / cone / sphere / torus answer 15 for V.
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Torus => 15,
        // `cxx:67-69` offset / other.
        SurfaceKind::Other => 10,
    }
}

/// `NbSamplesU(S, u1, u2)` (`cxx:76-112`).
pub fn nb_samples_u(s: &dyn Surface) -> i32 {
    if s.is_bezier_surface() {
        // `case GeomAbs_BezierSurface: return (3 + S->NbUPoles());`
        return 3 + s.nb_u_poles();
    }
    if s.is_bspline_surface() {
        // `case GeomAbs_BSplineSurface:` (`cxx:85-96`)
        let mut nbs = s.bspline_surface_uknots().map(|k| k.len() as i32).unwrap_or(0)
            * s.u_degree();
        if !bspline_u_rational(s) {
            nbs *= 2;
        }
        if nbs < 4 {
            nbs = 4;
        }
        return nbs;
    }
    match classify_surface(s) {
        // `case GeomAbs_Plane: return 2;` (`cxx:79-80`)
        SurfaceKind::Plane => 2,
        // `case GeomAbs_Torus: return 20;` (`cxx:97-98`)
        SurfaceKind::Torus => 20,
        // `cxx:100-107`: cylinder / cone / sphere / revolution / extrusion /
        // offset / other answer 10 for U.
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Other => 10,
    }
}

/// `GeomAbs_SurfaceOfRevolution` / `GeomAbs_SurfaceOfExtrusion`. Both are
/// `Geom_SweptSurface` subclasses in OCCT, and `GeomAdaptor_Surface::Load`
/// reports them by class (`GeomAdaptor_Surface.cxx`).
fn is_swept_surface(s: &dyn Surface) -> bool {
    s.revolution_basis_curve().is_some() || s.extrusion_basis_curve().is_some()
}

/// `Geom_BSplineSurface::IsVRational()`. OCCT asks the adaptor's
/// `IsVRational()`; the port reads it off the surface.
fn bspline_v_rational(s: &dyn Surface) -> bool {
    s.bspline_rational_flags().map(|(_, v)| v).unwrap_or(false)
}

/// `Geom_BSplineSurface::IsURational()`.
fn bspline_u_rational(s: &dyn Surface) -> bool {
    s.bspline_rational_flags().map(|(u, _)| u).unwrap_or(false)
}

/// `NbSamplesOnArc(A)` (`cxx:232-261`).
pub fn nb_samples_on_arc(a: &dyn Curve2d) -> i32 {
    // `GeomAbs_CurveType CurveType = A->GetType();`
    if let Some(nb_poles) = a.bezier_nb_poles() {
        // `case GeomAbs_BezierCurve: nbsOnC = A->NbPoles();`
        return nb_poles as i32;
    }
    match curve2d_kind(a) {
        // `case GeomAbs_Line: nbsOnC = 2;`
        GeomAbsCurveType::Line => 2,
        // `case GeomAbs_Circle/Ellipse/Hyperbola/Parabola: nbsOnC = 10;`
        GeomAbsCurveType::Circle
        | GeomAbsCurveType::Ellipse
        | GeomAbsCurveType::Hyperbola
        | GeomAbsCurveType::Parabola => 10,
        // `case GeomAbs_BSplineCurve: nbsOnC = 2 + A->NbKnots() * A->Degree();`
        GeomAbsCurveType::BSplineCurve => {
            let knots = a.bspline_nb_knots().unwrap_or(0) as i32;
            let degree = a.bspline_degree().unwrap_or(0) as i32;
            2 + knots * degree
        }
        // `default: nbsOnC = 10;` (the initial `nbsOnC = 5` is always
        // overwritten by one of the two switches).
        _ => 10,
    }
}

/// `Bounds(C, Ufirst, Ulast)` (`cxx:263-269`).
pub fn bounds(c: &dyn Curve2d) -> (f64, f64) {
    (c.first_parameter(), c.last_parameter())
}

/// `Project(C, P, Paramproj, Ptproj)` (`cxx:271-304`).
///
/// `epsX = 1.e-8` / `Nbu = 20` / `Tol = 1.e-5` are the local constants of the
/// source (`cxx:279-281`). The nearest of the sampled extrema wins, ties going
/// to the smallest index (`cxx:292-303`).
pub fn project(c: &dyn Curve2d, p: &GpPnt2d) -> Option<(f64, GpPnt2d)> {
    let eps_x = 1.0e-8;
    let nb_u = 20;
    let tol = 1.0e-5;

    let extrema = EPcOfExtPc2d::new(p, c, nb_u, eps_x, tol);
    if !extrema.is_done() {
        return None;
    }
    let nb_ext = extrema.nb_ext();
    if nb_ext == 0 {
        return None;
    }
    let mut index_min = 1;
    let mut dist2 = extrema.square_distance(1);
    for i in 2..=nb_ext {
        if extrema.square_distance(i) < dist2 {
            index_min = i;
            dist2 = extrema.square_distance(i);
        }
    }
    let _ = dist2;
    Some(extrema.point(index_min))
}

/// `Tolerance(V, C)` (`cxx:306-310`): `V->Resolution(C)`.
pub fn tolerance(v: &HVertex, c: &Curve2dAdaptor) -> f64 {
    v.resolution(c)
}

/// `Parameter(V, C)` (`cxx:312-316`): `V->Parameter(C)`.
pub fn parameter(v: &HVertex, c: &Curve2dAdaptor) -> f64 {
    v.parameter(c)
}

/// The `Epsilon(100.)` initial deflection of `IntPatch_Polyhedron`
/// (`IntPatch_Polyhedron.cxx:72`), kept next to the tool because it is the
/// same `Precision` service.
pub fn initial_deflection() -> f64 {
    occt_core::precision::epsilon(100.0)
}

/// `Precision::PConfusion()` re-export for the `IntPatch_PolyLine` initial
/// deflection (`IntPatch_PolyLine.cxx:26`).
pub fn p_confusion() -> f64 {
    Precision::PCONFUSION
}
