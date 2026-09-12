//! Surface/curve helpers used by FaceFace analytic pairs and result assembly.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpDir2d, GpLin2d, GpPln, GpPnt,
    GpPnt2d, GpSphere, GpVec, GpVec2d,
};
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::brep_surface::SurfaceKind;
use crate::builder::TopoBuilder;
use crate::intpatch;
use crate::inttools_data::{CurveKind, IntRange};
use crate::pcurve_full::{cone_params, make_pcurve_full};
use crate::shape::Face;
use crate::tgeometry::GeometryRegistry;

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
pub(crate) fn plane_cylinder_kind(pln: &GpPln, ax: &GpAx1) -> CurveKind {
    let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
    let az = GpVec::from_xyz(ax.direction().xyz()).normalized();
    let n_par = n.dot(&az);
    let n_perp = n.subtracted(&az.multiplied_scalar(n_par));
    if n_perp.magnitude() <= 1e-9 {
        CurveKind::Circle
    } else if n_par.abs() <= 1e-9 {
        CurveKind::Line
    } else {
        CurveKind::Ellipse
    }
}

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
pub(crate) fn cone_from_surface(s: &dyn Surface) -> Option<GpCone> {
    let (apex, ax, alpha) = cone_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(apex, z, &x).ok()?;
    GpCone::new(ax3, 0.0, alpha).ok()
}

/// Recover a `GpCylinder` from a cylinder surface via `cylinder_params`
/// (axis frame X direction is arbitrary — it does not affect the section
/// curve geometry, only its parametrization).
pub(crate) fn cylinder_from_surface(s: &dyn Surface) -> Option<GpCylinder> {
    let (center, ax, r) = cylinder_params(s)?;
    let z = GpDir::from_vec(&ax).ok()?;
    let x = perp_x_dir(&z);
    let ax3 = GpAx3::new(center, z, &x).ok()?;
    GpCylinder::new(ax3, r).ok()
}

/// Recover a `GpSphere` from a sphere surface (sphere frame is arbitrary).
pub(crate) fn sphere_from_surface(s: &dyn Surface) -> Option<GpSphere> {
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
        let a_dr = dist.abs() - minor;
        if a_dr > 1e-13 {
            return None; // plane misses the tube → IntAna_Empty
        }
        let dist = if a_dr.abs() < 1e-13 {
            if dist < 0.0 { -minor } else { minor }
        } else {
            dist
        };
        let a_dt = (minor * minor - dist * dist).max(0.0).sqrt();
        let center_on_plane = center.translated_vec(&n.multiplied_scalar(-dist));
        let normal = GpDir::from_vec(&n).ok()?;
        let mut out = vec![circle_gp(center_on_plane, normal, major + a_dt)?];
        if a_dr < -1e-13 && a_dt > tol {
            out.push(circle_gp(center_on_plane, normal, (major - a_dt).max(0.0))?);
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
    let a_dir = GpDir::from_vec(&az).ok()?;
    let normal = GpDir::from_vec(&n).ok()?;
    let e = a_dir.cross(&normal).ok()?;
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

/// Pcurve of `curve` on `face` over `range`, via a temporary edge registered in
/// the geometry side-table (cleaned up immediately after). Best-effort: returns
/// `None` when the pcurve cannot be computed.
pub(crate) fn pcurve_of_curve(
    curve: &Arc<dyn Curve>,
    range: IntRange,
    face: &Face,
) -> Option<Arc<dyn Curve2d>> {
    if !range.is_valid() || range.length() < 1e-12 {
        return None;
    }
    let builder = TopoBuilder::new();
    let edge = builder.make_edge(curve.clone(), range.first, range.last);
    let pc = make_pcurve_full(&edge, face).ok();
    GeometryRegistry::global().clear_shape(&edge.0);
    pc
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
