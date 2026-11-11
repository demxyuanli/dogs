//! Surface/curve helpers used by FaceFace analytic pairs and result assembly.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{
    GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpDir2d, GpLin2d, GpPln, GpPnt,
    GpPnt2d, GpSphere, GpVec, GpVec2d,
};
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::brep_surface::SurfaceKind;
use crate::intpatch;
use crate::inttools_data::IntRange;
use crate::pcurve_full::cone_params;
use crate::shape::Face;

use super::FaceFaceCurve;

/// Default half-extent used when an intersection curve reports an unbounded
/// parameter range, so it can still be represented in a `FaceFaceCurve`.
pub(crate) const DEFAULT_WINDOW: f64 = 4.0;

/// OCCT `IntTools_FaceFace::IndexType` — a total ordering used to sort the two
/// faces so the "higher" analytic type becomes `Face1`.
pub(crate) fn surface_sort_index(kind: SurfaceKind) -> usize {
    match kind {
        SurfaceKind::Plane => 0,
        SurfaceKind::Cylinder => 1,
        SurfaceKind::Cone => 2,
        SurfaceKind::Sphere => 3,
        SurfaceKind::Torus => 4,
        SurfaceKind::Other => 10,
    }
}

/// Cylinder parameters `(center, unit axis, radius)` recovered from a surface
/// by sampling its invariants (mirrors `pcurve_full::cylinder_params`).
pub(crate) fn cylinder_params(s: &dyn Surface) -> Option<(GpPnt, GpVec, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let axis = GpVec::from_pnts(&p0, &p1);
    let m = axis.magnitude();
    if m < 1e-9 {
        return None;
    }
    let ax = axis.divided(m);
    let q0 = s.d0(0.0, 0.0);
    let q1 = s.d0(PI, 0.0);
    let r = q0.distance(&q1) / 2.0;
    if r <= 1e-9 {
        return None;
    }
    let center = mid(&q0, &q1);
    let (nu, nv) = (8, 6);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = 2.0 * PI * i as f64 / nu as f64;
            let v = (j as f64 - nv as f64 / 2.0) / 2.0;
            let p = s.d0(u, v);
            if (dist_to_axis(&p, &center, &ax) - r).abs() > 1e-3 * r.max(1.0) {
                return None;
            }
        }
    }
    Some((center, ax, r))
}

pub(crate) fn mid(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

fn dist_to_axis(p: &GpPnt, c: &GpPnt, ax: &GpVec) -> f64 {
    let rel = GpVec::from_pnts(c, p);
    rel.subtracted(&ax.multiplied_scalar(rel.dot(ax))).magnitude()
}

/// The analytic family of a plane×cylinder intersection, recomputed from the
/// plane normal and the cylinder axis (matches `intersect_plane_cylinder`).

/// A unit direction perpendicular to `z` (axis for a cone/torus frame).
fn perp_x_dir(z: &GpDir) -> GpDir {
    let base = if z.x().abs() < 0.9 {
        GpDir::new(1.0, 0.0, 0.0).expect("x")
    } else {
        GpDir::new(0.0, 1.0, 0.0).expect("y")
    };
    base.cross(z).unwrap_or(base)
}

/// Recover a `GpCone` from a cone surface by sampling its invariants
/// (apex + axis + semi-angle from `pcurve_full`). The cone is anchored at its
/// apex with radius 0 — `IntAna_QuadQuadGeo` only reads apex / axis / semi-angle
/// geometry (plus `radius` in one axis-flip heuristic, where 0 is safe).
pub(crate) fn gp_cone_from_surface(s: &dyn Surface) -> Option<GpCone> {
    let (apex, ax, alpha) = cone_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(apex, z, &x).ok()?;
    GpCone::new(ax3, 0.0, alpha).ok()
}

/// Recover a `GpCylinder` from a cylinder surface. Port of
/// `IntSurf_Quadric::SetValue(const gp_Cylinder&)` (`IntSurf_Quadric.cxx:123-131`):
/// the quadric keeps the surface's own `Position()`, including its X direction.
/// `IntAna_QuadQuadGeo` re-bases the section circle onto that same frame
/// (`IntPatch_ImpImpIntersection.cxx:3300` -> `AdjustToSeam`, `:3895-3903`), so
/// dropping the U origin here mis-parameterizes the section curve and shifts the
/// pave vertices that split the face. The sampled reconstruction is only a
/// fallback for surfaces that do not expose their cylinder.
pub(crate) fn gp_cylinder_from_surface(s: &dyn Surface) -> Option<GpCylinder> {
    if let Some(c) = s.gp_cylinder() {
        return Some(c);
    }
    let (center, ax, r) = cylinder_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(center, z, &x).ok()?;
    GpCylinder::new(ax3, r).ok()
}

/// Recover a `GpSphere` from a sphere surface (sphere frame is arbitrary).
pub(crate) fn gp_sphere_from_surface(s: &dyn Surface) -> Option<GpSphere> {
    let (c, r) = intpatch::sphere_params(s)?;
    let z = GpDir::new(0.0, 0.0, 1.0).ok()?;
    let x = GpDir::new(1.0, 0.0, 0.0).ok()?;
    let ax3 = GpAx3::new(c, z, &x).ok()?;
    GpSphere::new(ax3, r).ok()
}

/// Plane coefficients `A·x + B·y + C·z + D = 0` with unit normal (OCCT form).
fn plane_coeffs(p: &GpPln) -> (f64, f64, f64, f64) {
    let n = GpVec::from_xyz(p.axis().direction().xyz()).normalized();
    let loc = p.location();
    let d = -(n.x() * loc.x() + n.y() * loc.y() + n.z() * loc.z());
    (n.x(), n.y(), n.z(), d)
}

/// A `GpCirc` in the plane through `center` with normal `normal` and `radius`.
fn circle_gp(center: GpPnt, normal: GpDir, radius: f64) -> Option<GpCirc> {
    let x_dir = perp_x_dir(&normal);
    let ax2 = GpAx2::new(center, normal, x_dir).ok()?;
    Some(GpCirc::new(ax2, radius))
}

/// Plane ∩ torus circles. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln,
/// gp_Torus)` (IntAna_QuadQuadGeo.cxx): up to two circles when the torus axis
/// is parallel to the plane normal (perpendicular cut → radii `major ± dt`,
/// `dt = √(minor² − dist²)` from the center) or perpendicular to it (axis in
/// the plane through the center → two `minor` circles at `±major`). Returns
/// `None` for any other orientation — the general (non-planar) torus section,
/// which OCCT routes to the numeric walker.
pub(crate) fn plane_torus_circles(
    pln: &GpPln,
    center: GpPnt,
    ax: &GpVec,
    major: f64,
    minor: f64,
    tol: f64,
) -> Option<Vec<GpCirc>> {
    if minor >= major {
        return None; // degenerate torus → IntAna_NoGeometricSolution
    }
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let az = ax.normalized();
    let n_par = n.dot(&az);
    let (a, b, c, d) = plane_coeffs(pln);
    let dist = a * center.x() + b * center.y() + c * center.z() + d;

    if (n_par.abs() - 1.0).abs() <= 1e-12 {
        // Axis ∥ plane normal → perpendicular cut.
        let dr = dist.abs() - minor;
        if dr > 1e-13 {
            return None; // plane misses the tube → IntAna_Empty
        }
        let dist = if dr.abs() < 1e-13 {
            if dist < 0.0 { -minor } else { minor }
        } else {
            dist
        };
        let dt = (minor * minor - dist * dist).max(0.0).sqrt();
        let center_on_plane = center.translated_vec(&n.multiplied_scalar(-dist));
        let normal = GpDir::from_vec(&n).ok()?;
        let mut out = vec![circle_gp(center_on_plane, normal, major + dt)?];
        if dr < -1e-13 && dt > tol {
            out.push(circle_gp(center_on_plane, normal, (major - dt).max(0.0))?);
        }
        return Some(out);
    }

    if n_par.abs() > 1e-12 {
        return None; // oblique → IntAna_NoGeometricSolution (numeric)
    }
    // Axis ⊥ normal → plane must contain the torus axis through the center.
    if dist.abs() > 1e-14 {
        return None;
    }
    let dir = GpDir::from_vec(&az).ok()?;
    let normal = GpDir::from_vec(&n).ok()?;
    let e = dir.cross(&normal).ok()?;
    let c1 = center.translated_vec(&GpVec::from_xyz(e.xyz()).multiplied_scalar(major));
    let c2 = center.translated_vec(&GpVec::from_xyz(e.xyz()).multiplied_scalar(-major));
    Some(vec![
        circle_gp(c1, normal, minor)?,
        circle_gp(c2, normal, minor)?,
    ])
}

/// Parameter interval of the (unit-speed) 2D line that lies inside the UV
/// rectangle — a Liang–Barsky clip. Returns `None` when the line misses the
/// rectangle entirely. Unbounded rectangle dimensions impose no constraint.
pub(crate) fn line_in_uv_rect(lin: &GpLin2d, bounds: (f64, f64, f64, f64)) -> Option<(f64, f64)> {
    let (umin, umax, vmin, vmax) = bounds;
    let loc = lin.pos.loc;
    let dir = lin.pos.vdir;
    let (x0, y0) = (loc.x(), loc.y());
    let (dx, dy) = (dir.x, dir.y);
    let mut t0 = f64::NEG_INFINITY;
    let mut t1 = f64::INFINITY;
    for (lo, hi, p0, dp) in [(umin, umax, x0, dx), (vmin, vmax, y0, dy)] {
        if !lo.is_finite() || !hi.is_finite() {
            continue;
        }
        if dp.abs() < 1e-15 {
            if p0 < lo - 1e-9 || p0 > hi + 1e-9 {
                return None;
            }
            continue;
        }
        let ta = (lo - p0) / dp;
        let tb = (hi - p0) / dp;
        let (e, l) = if ta < tb { (ta, tb) } else { (tb, ta) };
        t0 = t0.max(e);
        t1 = t1.min(l);
        if t0 > t1 + 1e-9 {
            return None;
        }
    }
    if !t0.is_finite() || !t1.is_finite() || t1 - t0 <= 1e-9 {
        return None;
    }
    Some((t0, t1))
}

/// Unit-speed 2D line through two distinct points (`None` when coincident).
pub(crate) fn lin2d_through(q0: GpPnt2d, q1: GpPnt2d) -> Option<GpLin2d> {
    let v = GpVec2d::new(q1.x() - q0.x(), q1.y() - q0.y());
    let m = v.magnitude();
    if m < 1e-12 {
        return None;
    }
    let d = GpDir2d::from_vec2d(&v).ok()?;
    Some(GpLin2d::from_pnt_dir(q0, d))
}

/// The finite range of a curve, or a default window for unbounded curves.
pub(crate) fn curve_range(curve: &dyn Curve) -> IntRange {
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if a.is_finite() && b.is_finite() && b > a {
        IntRange::new_unchecked(a, b)
    } else {
        IntRange::new_unchecked(-DEFAULT_WINDOW, DEFAULT_WINDOW)
    }
}

/// Pcurve of `curve` on `face` over `range` —
/// `GeomInt_IntSS::BuildPCurves` (`GeomInt_IntSS_1.cxx:1172-1303`).
///
/// `cxx:1187-1202`: the 2D curve is `GeomProjLib::Curve2d(theCurve, theFirst,
/// theLast, theSurface, Tol)` — `ProjLib_ProjectedCurve` over exactly that
/// range with the tolerance it reports. `ShapeConstruct_ProjectCurveOnSurface`
/// (`make_pcurve_full`) is a different OCCT class: it interpolates a
/// *fixed-count* uniform sample set, so its chord error grows with the range
/// (`h^4 |f''''|`) and it is not what this call site must use.
///
/// The three `cxx` branches after it are ported below: the
/// `Tol += Confusion` retry when the UV-box projection comes back null
/// (`cxx:1200-1204`), the `Geom2d_BSplineCurve` knot re-parametrization
/// (`cxx:1206-1222`), and the degree-1 segment for ranges too short for
/// `ProjLib_Function` (`cxx:1223-1282`). The U-periodic window shift
/// (`cxx:1284-1303`) is [`adjust_u_periodic_on_face`].
pub(crate) fn pcurve_of_curve(
    curve: &Arc<dyn Curve>,
    range: IntRange,
    face: &Face,
) -> Option<Arc<dyn Curve2d>> {
    if !range.is_valid() {
        return None;
    }
    let surf = crate::brep_tool::BRepTool::face_surface(face)?;
    let mut tol = occt_core::precision::CONFUSION;
    let mut pc = if range.length() > 2.0e-9 {
        // `cxx:1189-1199`: `ProjLib_Function` shrinks the range by 1.e-09, so
        // only a wider span reaches `GeomProjLib::Curve2d`.
        let mut projected =
            occt_geom::geom_proj_lib::curve_2d_range(curve, range.first, range.last, &surf, &mut tol);
        if projected.is_none() {
            // `cxx:1200-1204`: retry against the surface's own UV box (the
            // five-argument `GeomProjLib::Curve2d`) with the widened tolerance.
            tol += occt_core::precision::CONFUSION;
            projected = occt_geom::geom_proj_lib::curve_2d_with_tolerance(curve, &surf, &mut tol);
        }
        projected?
    } else {
        short_range_pcurve(curve, range, &surf, &mut tol)?
    };
    // `cxx:1206-1222`: the approximation can land on knots that miss the
    // requested ends by more than `Precision::PConfusion()`; re-parameterize
    // them onto `[theFirst, theLast]` the way `Geom2d_TrimmedCurve` expects.
    if let Some(mut bs) = pc.bspline_copy2d() {
        if bs.first_parameter() - range.first > occt_core::precision::PCONFUSION
            || range.last - bs.last_parameter() > occt_core::precision::PCONFUSION
        {
            if let Some((mut knots, _)) = pc.bspline_distinct_knots_mults() {
                occt_core::bspl::knots::reparametrize(range.first, range.last, &mut knots);
                let _ = bs.set_knots(&knots);
                pc = Arc::new(bs);
            }
        }
    }
    adjust_u_periodic_on_face(&mut pc, face, range.first, range.last);
    Some(pc)
}

/// `GeomInt_IntSS::BuildPCurves` short-range branch
/// (`GeomInt_IntSS_1.cxx:1223-1282`).
///
/// `theLast - theFirst` is inside `ProjLib_Function`'s own epsilon, so the
/// reference treats the 2D image as a segment: project both 3D ends with
/// `Extrema_ExtPS` (`ParametersOfNearestPointOnSurface`, `cxx:151-175`), emit a
/// degree-1 B-spline with knots `[theFirst, theLast]` and multiplicity 2, and
/// fold the mid-point's 3D-to-2D distance into `theTol`.
fn short_range_pcurve(
    curve: &Arc<dyn Curve>,
    range: IntRange,
    surf: &Arc<dyn occt_geom::surface::Surface>,
    tol: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    // `cxx:1225`: `(theLast - theFirst) > Epsilon(std::abs(theFirst))`.
    if range.length() <= occt_core::precision::epsilon(range.first) {
        return None;
    }
    let conf = occt_core::precision::CONFUSION;
    let p1 = occt_geom::geom_api::project_point_on_surface(surf.as_ref(), &curve.d0(range.first), conf)?;
    let p2 = occt_geom::geom_api::project_point_on_surface(surf.as_ref(), &curve.d0(range.last), conf)?;
    let a_p2d1 = occt_core::gp::GpPnt2d::new(p1.u, p1.v);
    let a_p2d2 = occt_core::gp::GpPnt2d::new(p2.u, p2.v);
    // `cxx:1258`: a degenerate 2D span leaves the curve null.
    if a_p2d1.distance(&a_p2d2) <= occt_core::precision::REAL_SMALL {
        return None;
    }
    // `cxx:1266-1270`: knots `[theFirst, theLast]`, both multiplicity 2.
    let knots = vec![range.first, range.first, range.last, range.last];
    let bs = occt_geom2d::bspline_curve::Geom2dBSplineCurve::from_flat(
        vec![a_p2d1.x(), a_p2d2.x()],
        vec![a_p2d1.y(), a_p2d2.y()],
        None,
        knots,
        1,
        false,
    )
    .ok()?;
    // `cxx:1272-1278`: widen `theTol` by the mid-point's deviation.
    let p_mid = curve.d0(0.5 * (range.first + range.last));
    let pmid2d = occt_core::gp::GpPnt2d::new(
        0.5 * (a_p2d1.x() + a_p2d2.x()),
        0.5 * (a_p2d1.y() + a_p2d2.y()),
    );
    let a_pc = surf.d0(pmid2d.x(), pmid2d.y());
    *tol = p_mid.distance(&a_pc).max(*tol);
    Some(Arc::new(bs))
}

/// `GeomInt_IntSS::BuildPCurves` tail (`GeomInt_IntSS_1.cxx:1285-1305`):
/// "Recadre dans le domaine UV de la face" — on a U-periodic surface the
/// freshly projected 2D curve is translated by whole periods so that the point
/// at the middle of the 3D range `[first, last]` (`aTm = .5 * (aFirst +
/// aLast)`) lands inside `[umin, umax]`. The offset comes from
/// `GeomInt::AdjustPeriodic`, and the translation is applied to the whole
/// curve, so the pcurve keeps its analytic type.
///
/// `theUmin`/`theUmax` are the surface adaptor's first/last U parameter, i.e.
/// the face's UV box (`IntTools_FaceFace.cxx:1035-1060`).
fn adjust_u_periodic_on_face(
    pc: &mut Arc<dyn Curve2d>,
    face: &Face,
    first: f64,
    last: f64,
) {
    let Some(surf) = crate::brep_tool::BRepTool::face_surface(face) else {
        return;
    };
    if !surf.is_u_periodic() {
        return;
    }
    let (umin, umax, _, _) = crate::brep_tool::BRepTool::uv_bounds(face);
    let period = surf.u_period();
    if !period.is_finite() || period <= 0.0 {
        return;
    }
    let midpoint = pc.d0(0.5 * (first + last)).x();
    let (_, du) = crate::geom_int::geom_int_adjust_periodic(
        midpoint,
        umin,
        umax,
        period,
        occt_core::precision::PCONFUSION,
    );
    if du == 0.0 {
        return;
    }
    let mut t = occt_core::gp::GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(du, 0.0));
    *pc = Arc::from(pc.transformed(&t));
}

/// Whether two curves represent the same intersection (same kind, same range
/// and a coincident mid-point) — used to deduplicate the result.
pub(crate) fn same_curve(a: &FaceFaceCurve, b: &FaceFaceCurve) -> bool {
    if a.kind != b.kind {
        return false;
    }
    if (a.range.first - b.range.first).abs() > 1e-7 || (a.range.last - b.range.last).abs() > 1e-7 {
        return false;
    }
    let ta = 0.5 * (a.range.first + a.range.last);
    let tb = 0.5 * (b.range.first + b.range.last);
    a.curve.d0(ta).distance(&b.curve.d0(tb)) < 1e-6
}
