//! Curve/surface type dispatch used by EdgeFace.
//!
//! Every predicate here is the OCCT `GetType()` test that `IntTools_EdgeFace`
//! branches on (`IntTools_EdgeFace.cxx:74,247,532,620,650,798,827`), answered
//! from the curve/surface's own `gp_*` accessors - the same source as
//! `GeomAdaptor_Curve::load` / `GeomAdaptor_Surface::load`. No geometric
//! sampling is involved.
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::PCONFUSION;
use occt_geom::{Curve, Surface};

use crate::brep_surface::{classify_surface, SurfaceKind};

// ---------------------------------------------------------------------------
// Curve classification helpers
// ---------------------------------------------------------------------------

/// Coarse analytic kind of a 3D curve, mirroring the `GeomAbs_CurveType` subset
/// the algorithm branches on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CurveKind {
    Line,
    Circle,
    BSpline,
    Other,
}

/// `(center, radius, unit axis direction)` of a `GeomAbs_Circle` curve.
///
/// `CurveTool::Circle` (`IntTools_EdgeFace.cxx:800,829`), i.e. the
/// `Geom_Circle::Circ()` downcast; `None` for every other curve type. The
/// previous body fitted a circle through 8 sampled points, so a B-spline whose
/// image was a circle also answered `Some`.
pub(crate) fn circle_geometry(c: &dyn Curve) -> Option<(GpPnt, f64, GpVec)> {
    let circ = c.gp_circ()?;
    let nrm = GpVec::from_xyz(circ.pos.direction().xyz()).normalized();
    if nrm.xyz().square_modulus() < 1e-30 {
        return None;
    }
    Some((circ.location(), circ.radius, nrm))
}

/// `Adaptor3d_Curve::GetType()` reduced to the arms `IntTools_EdgeFace`
/// branches on.
///
/// `GeomAdaptor_Curve::load` (`GeomAdaptor_Curve.cxx:252-311`) unwraps a
/// `Geom_TrimmedCurve` to its basis and then tests the concrete class -
/// `gp_line()` for `GeomAbs_Line`, `gp_circ()` for `GeomAbs_Circle`, and the
/// Bezier/BSpline classes for `GeomAbs_BezierCurve`/`GeomAbs_BSplineCurve`
/// (both surfaced by [`Curve::nurbs_degree`]). Everything else - ellipse,
/// parabola, hyperbola, offset curve, other - takes the `default` arm, which
/// `IntTools_EdgeFace` treats as "neither Line nor Circle nor B-spline":
/// `cxx:74,309,532,701` (Line), `cxx:650,798,827` (Circle) and
/// `cxx:532` (BSpline/Bezier).
///
/// The previous body sampled 8 points of the curve and tested collinearity /
/// equidistance, so a `Geom_BSplineCurve` whose image was a straight line or a
/// circle was reported as `Line`/`Circle` and sent into the analytic arms.
pub(crate) fn curve_kind(c: &dyn Curve) -> CurveKind {
    if c.gp_line().is_some() {
        CurveKind::Line
    } else if c.gp_circ().is_some() {
        CurveKind::Circle
    } else if c.nurbs_degree().is_some() {
        CurveKind::BSpline
    } else {
        CurveKind::Other
    }
}

/// `BRepAdaptor_Curve::Resolution` = `GeomAdaptor_Curve::Resolution`
/// (`GeomAdaptor_Curve.cxx:1116-1148`): the [`Curve::resolution`] of the curve
/// itself, which `IntTools_EdgeFace` calls as `myC.Resolution(myCriteria)`
/// (`IntTools_EdgeFace.cxx:326-327` in `MakeType`, `:373` in `CheckTouch`) for
/// its "whole range" tests. No sampling is involved.
pub(crate) fn curve_resolution(c: &dyn Curve, tol: f64) -> f64 {
    c.resolution(tol)
}

// ---------------------------------------------------------------------------
// Analytic surface geometry recovery (dyn Surface cannot be downcast)
// ---------------------------------------------------------------------------

pub(crate) fn surface_sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// `(origin, unit-u, unit-v)` frame of a planar surface (natural parameterization).
pub(crate) fn plane_frame(s: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let (o, x, y) = s.d1(0.0, 0.0);
    let mx = x.magnitude();
    let my = y.magnitude();
    if mx < 1e-30 || my < 1e-30 {
        return None;
    }
    Some((o, x.divided(mx), y.divided(my)))
}

/// `ElSLib::Parameters(gp_Pln, P)` plus the signed plane distance, read from the
/// surface's own `gp_Pln` placement. This is OCCT's analytic arm for a planar
/// surface (`IntCurveSurface_InterUtils.pxx:906-909`); it needs no reconstructed
/// frame, so it is also the fallback when `plane_frame`'s `d1(0, 0)` sampling is
/// degenerate.
fn analytic_plane_projection(s: &dyn Surface, p: &GpPnt) -> Option<(f64, f64, f64)> {
    let pln = s.gp_pln()?;
    let (u, v) = occt_core::elib::slib::plane_parameters(&pln.pos, p);
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    if n.square_magnitude() < 1e-30 {
        return None;
    }
    let d = GpVec::from_pnts(&pln.location(), p);
    Some((u, v, d.dot(&n).abs()))
}

/// Exact projection of `p` onto a planar surface: `(u, v, distance)`.
pub(crate) fn plane_projection(s: &dyn Surface, p: &GpPnt) -> (f64, f64, f64) {
    let (o, x, y) = match plane_frame(s) {
        Some(f) => f,
        None => {
            // Degenerate `d1(0, 0)` frame: take the analytic `gp_Pln` placement
            // (`ElSLib::Parameters`), which does not need the frame; else fall
            // back to OCCT's point-on-surface projection (`Extrema_ExtPS`).
            if let Some(proj) = analytic_plane_projection(s, p) {
                return proj;
            }
            return match occt_geom::geom_api::project_point_on_surface(
                s,
                p,
                occt_core::precision::CONFUSION,
            ) {
                Some(ps) => (ps.u, ps.v, ps.distance),
                None => (0.0, 0.0, f64::INFINITY),
            };
        }
    };
    let n = x.crossed(&y);
    let m = n.magnitude();
    if m < 1e-30 {
        if let Some(proj) = analytic_plane_projection(s, p) {
            return proj;
        }
        return match occt_geom::geom_api::project_point_on_surface(
            s,
            p,
            occt_core::precision::CONFUSION,
        ) {
            Some(ps) => (ps.u, ps.v, ps.distance),
            None => (0.0, 0.0, f64::INFINITY),
        };
    }
    let n = n.divided(m);
    let d = GpVec::from_pnts(&o, p);
    let dist = d.dot(&n).abs();
    let proj = p.translated_vec(&n.multiplied_scalar(-d.dot(&n)));
    let w = GpVec::from_pnts(&o, &proj);
    (w.dot(&x), w.dot(&y), dist)
}

/// `gp_Pln::Axis().Direction()` of a `GeomAbs_Plane` surface - the plane normal
/// `IsCoplanar` compares the circle axis against (`IntTools_EdgeFace.cxx:804-808`,
/// `SurfaceTool::Plane`). `GpVec::zero()` for every other surface type.
pub(crate) fn surface_plane_normal(s: &dyn Surface) -> GpVec {
    match s.gp_pln() {
        Some(pln) => GpVec::from_xyz(pln.pos.direction().xyz()).normalized(),
        None => GpVec::zero(),
    }
}

/// `Adaptor3d_Surface::GetType()` (`GeomAdaptor_Surface::load`,
/// `GeomAdaptor_Surface.cxx:422-513`), i.e. [`classify_surface`]. This is the
/// predicate `IntTools_EdgeFace` branches on - `cxx:74,620` (Cylinder),
/// `cxx:650,798,827` (Plane), and `IsEqDistance`'s `cxx:249,263,282`
/// (Cylinder / Cone / Torus).
///
/// The previous body was a two-stage sampler: recover a quadric's parameters
/// from `d0` evaluations and then re-verify them on an 8x8 grid against
/// home-made `1e-4` tolerances. That reported a B-spline patch that merely
/// resembled a cylinder / cone / torus as the analytic type, and could equally
/// reject a genuine quadric whose window shifted the sampled frame.
pub(crate) fn edge_face_surface_kind(s: &dyn Surface) -> SurfaceKind {
    classify_surface(s)
}

/// Shortcut distance for a point lying on the axis of an analytic surface.
///
/// Port of `IntTools_EdgeFace::IsEqDistance` (`cxx:240-299`): the branch is
/// `aBAS.GetType()`, the frame comes from `aBAS.Cylinder()` / `Cone()` /
/// `Torus()` - `SurfaceTool::Cylinder` etc. - and the tolerance is the literal
/// `1.e-7` the caller passes.
pub(crate) fn is_eq_distance(p: &GpPnt, s: &dyn Surface) -> Option<f64> {
    const TOL: f64 = 1e-7;
    match edge_face_surface_kind(s) {
        SurfaceKind::Cylinder => {
            let cyl = s.gp_cylinder()?;
            let z = GpVec::from_xyz(cyl.pos.direction().xyz());
            let dc = GpVec::from_pnts(&cyl.location(), p).crossed(&z).magnitude();
            if dc < TOL {
                Some(cyl.radius)
            } else {
                None
            }
        }
        SurfaceKind::Cone => {
            let cone = s.gp_cone()?;
            // `gp_Cone::Apex()` (`Local.para` of `gp_Cone.hxx`): the apex is
            // `Location - V * (RefRadius / tan(SemiAngle))`, since the radius is
            // `v * tan(SemiAngle) + RefRadius` (`ElSLib::ConeValue`).
            let apex = cone
                .location()
                .translated_vec(&GpVec::from_xyz(cone.pos.direction().xyz()).multiplied_scalar(
                    -(cone.radius / cone.semi_angle.tan()),
                ));
            let z = GpVec::from_xyz(cone.pos.direction().xyz());
            let dc = GpVec::from_pnts(&apex, p).crossed(&z).magnitude();
            if dc < TOL {
                Some(p.distance(&apex) * cone.semi_angle.tan())
            } else {
                None
            }
        }
        SurfaceKind::Torus => {
            let tor = s.gp_torus()?;
            let dc = (tor.location().distance(p) - tor.major_radius).abs();
            if dc < TOL {
                Some(tor.minor_radius)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// `IntTools_Tools::IsDirsCoinside`: unit directions equal or opposite within
/// `0.0002` (measured as the distance between the unit-sphere points).
pub(crate) fn dirs_coinside(d1: &GpVec, d2: &GpVec) -> bool {
    let d = d1.subtracted(d2).magnitude();
    d < 0.0002 || (2.0 - d).abs() < 0.0002
}

/// `IsCoplanar` (`IntTools_EdgeFace.cxx:787-811`): the curve axis direction
/// coincides with the plane normal, `aCType == GeomAbs_Circle &&
/// aSType == GeomAbs_Plane` being the gate.
pub(crate) fn is_coplanar(curve: &dyn Curve, surface: &dyn Surface) -> bool {
    if curve_kind(curve) != CurveKind::Circle || edge_face_surface_kind(surface) != SurfaceKind::Plane {
        return false;
    }
    let (_, _, n) = match circle_geometry(curve) {
        Some(g) => g,
        None => return false,
    };
    dirs_coinside(&n, &surface_plane_normal(surface))
}

/// `IsRadius` (`IntTools_EdgeFace.cxx:816-840`): the circle center lies at
/// distance ~ radius from the plane, `aCType == GeomAbs_Circle &&
/// aSType == GeomAbs_Plane` being the gate. `aPln.Distance(aCenter)` is measured
/// from the plane's own `gp_Pln` placement (`SurfaceTool::Plane`), not from a
/// sampled `d0`.
pub(crate) fn is_radius(curve: &dyn Curve, surface: &dyn Surface, criteria: f64) -> bool {
    if curve_kind(curve) != CurveKind::Circle || edge_face_surface_kind(surface) != SurfaceKind::Plane {
        return false;
    }
    let (center, r, _) = match circle_geometry(curve) {
        Some(g) => g,
        None => return false,
    };
    let Some(pln) = surface.gp_pln() else {
        return false;
    };
    let d = GpVec::from_pnts(&pln.location(), &center).dot(&surface_plane_normal(surface)).abs();
    (d - r).abs() < criteria
}

// ---------------------------------------------------------------------------
// EdgeFace
// ---------------------------------------------------------------------------
