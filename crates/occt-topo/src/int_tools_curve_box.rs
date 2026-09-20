//! Remaining `IntTools_Tools` helpers used by `PerformFF` / `MakeBlocks`.
//!
//! Source:
//! - `IntTools_Tools::CheckCurve` at `IntTools_Tools.cxx:549`
//! - `IntTools_Tools::IsOnPave` at 579
//! - `IntTools_Tools::VertexParameters` at 593
//! - `IntTools_Tools::VertexParameter` at 616
//! - `IntTools_Tools::IsOnPave1` at 627
//! - `IntTools_Tools::IsInRange` at 643
//! - `IntTools_Tools::SegPln` at 661 (kept as a documented translation
//!   boundary: the plane-segment classifier is not required by PerformFF)
//! - `IntTools_Tools::ComputeIntRange` at 783
//! - `Bnd_Box::IsThin` (the CheckCurve size gate)
//!
//! `CheckCurve` builds a bounding box of the 3D curve (enlarged by
//! `max(Tolerance, TangentialTolerance)`) and rejects the curve when that
//! box is thin in every coordinate direction (`Bnd_Box::IsThin`,
//! `3 * Precision::Confusion()`). A line-sized box (thin in two axes) is
//! kept — that is the plane/plane section case.

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::Curve;

use crate::inttools_data::IntRange;

/// `IntTools_Tools::CheckCurve` comparison tolerance: two vertices at
/// Confusion plus Confusion as the minimal distance between them.
pub const CHECK_CURVE_THIN: f64 = 3.0 * CONFUSION;

/// `Bnd_Box::IsThin(aTol)` — true only when every raw extent is smaller than `tol`.
pub fn box_is_thin(box_: &BndBox, tol: f64) -> bool {
    box_.is_thin(tol)
}

/// Build a bounding box of `curve` over `[first, last]`, enlarged by `tol`.
///
/// `BndLib_Add3dCurve::Add(theAdaptor, U1, U2, Tol, B)` — in 8.0.0 that forwards
/// to `GeomBndLib_Curve`, whose faithful port (per-type dispatch, analytic arms,
/// tolerance as the box gap) is [`crate::geom_bnd_lib_curve3d::box_curve`]. The
/// previous body sampled 48 points, a stand-in with no OCCT counterpart
/// (audit A19/T-55 family).
pub fn add_curve_to_box(curve: &dyn Curve, first: f64, last: f64, tol: f64, box_: &mut BndBox) {
    let (a, b) = finite_or_unit(curve, first, last);
    let b_curve = crate::geom_bnd_lib_curve3d::box_curve(curve, a, b, tol);
    box_.add_box(&b_curve);
}

fn finite_or_unit(curve: &dyn Curve, first: f64, last: f64) -> (f64, f64) {
    let a = if first.is_finite() {
        first
    } else {
        curve.first_parameter()
    };
    let b = if last.is_finite() {
        last
    } else {
        curve.last_parameter()
    };
    if a.is_finite() && b.is_finite() && (b - a).abs() > PCONFUSION {
        return (a.min(b), a.max(b));
    }
    (-1.0, 1.0)
}

/// `IntTools_Tools::CheckCurve(theCurve, theBox)`.
///
/// Returns `false` when the 3D curve is missing or the box is thin in every
/// direction (`3 * Confusion`). On success `the_box` is the enlarged bounding
/// box of the curve.
pub fn check_curve(
    curve: Option<&dyn Curve>,
    first: f64,
    last: f64,
    tolerance: f64,
    tangential: f64,
    the_box: &mut BndBox,
) -> bool {
    let Some(c) = curve else {
        return false;
    };
    *the_box = BndBox::new();
    let enlarge = tolerance.max(tangential);
    add_curve_to_box(c, first, last, enlarge, the_box);
    !box_is_thin(the_box, CHECK_CURVE_THIN)
}

/// Convenience: check an `IntTools_Curve`-shaped payload (curve + range + tols).
pub fn check_curve_payload(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    tolerance: f64,
    tangential: f64,
) -> Option<BndBox> {
    let mut box_ = BndBox::new();
    if check_curve(Some(curve), first, last, tolerance, tangential, &mut box_) {
        Some(box_)
    } else {
        None
    }
}

/// `IntTools_Tools::IsOnPave(aT1, aRange, aTolerance)`.
pub fn is_on_pave(t1: f64, range: IntRange, tolerance: f64) -> bool {
    let first_is = (range.first - t1).abs() < tolerance;
    let last_is = (range.last - t1).abs() < tolerance;
    first_is || last_is
}

/// `IntTools_Tools::IsOnPave1(aTR, aCPRange, aTol)`.
///
/// True when `aTR` is inside the range expanded by `aTol`, or on either
/// bound within `aTol`.
pub fn is_on_pave1(a_tr: f64, range: IntRange, a_tol: f64) -> bool {
    if is_on_pave(a_tr, range, a_tol) {
        return true;
    }
    a_tr > range.first - a_tol && a_tr < range.last + a_tol
}

/// `IntTools_Tools::IsInRange(aR, aRA, aTol)` — overlap of two 1D ranges
/// expanded by `aTol`.
pub fn is_in_range(a_r: IntRange, a_ra: IntRange, a_tol: f64) -> bool {
    let f = a_r.first - a_tol;
    let l = a_r.last + a_tol;
    let fa = a_ra.first;
    let la = a_ra.last;
    (fa >= f && fa <= l) || (la >= f && la <= l) || (fa <= f && la >= l)
}

/// `IntTools_Tools::VertexParameters` — pick the vertex parameter of a
/// common part on each edge, falling back to the mid of Range1 / Ranges2(1).
pub fn vertex_parameters(
    range1: IntRange,
    vertex_parameter1: f64,
    range2: IntRange,
    vertex_parameter2: f64,
) -> (f64, f64) {
    let mut t1 = 0.5 * (range1.first + range1.last);
    if vertex_parameter1 >= range1.first && vertex_parameter1 <= range1.last {
        t1 = vertex_parameter1;
    }
    let mut t2 = 0.5 * (range2.first + range2.last);
    if vertex_parameter2 >= range2.first && vertex_parameter2 <= range2.last {
        t2 = vertex_parameter2;
    }
    (t1, t2)
}

/// `IntTools_Tools::VertexParameter` — Range1 mid, or VertexParameter1 when
/// it sits inside Range1.
pub fn vertex_parameter(range1: IntRange, vertex_parameter1: f64) -> f64 {
    let mut t = 0.5 * (range1.first + range1.last);
    if vertex_parameter1 >= range1.first && vertex_parameter1 <= range1.last {
        t = vertex_parameter1;
    }
    t
}

/// `IntTools_Tools::ComputeIntRange(Tol1, Tol2, Angle)`.
///
/// The extra 3D range needed around a plane/plane intersection line so that
/// the face tolerances still cover the wedge. Straight copy of
/// `IntTools_Tools.cxx:783`.
pub fn compute_int_range(tol1: f64, tol2: f64, angle: f64) -> f64 {
    use std::f64::consts::PI;
    if (PI * 0.5 - angle).abs() < ANGULAR {
        return tol2;
    }
    let an_angle = if angle > PI * 0.5 { PI - angle } else { angle };
    let a1 = tol1 * (PI * 0.5 - an_angle).tan();
    let a2 = tol2 / an_angle.sin();
    a1 + a2
}

/// Face/face tangential tolerance of an intersection curve.
///
/// Non-planar: `max(TolF1, TolF2)`. Plane/plane: `sqrt(Dt^2 + TolF1^2)`
/// floored at `max(TolF1, TolF2)` (`IntTools_FaceFace::ComputeTolReached3d`).
pub fn curve_tangential_tolerance(
    f1_planar: bool,
    f2_planar: bool,
    t1: f64,
    t2: f64,
    plane_angle: Option<f64>,
) -> f64 {
    let t_max = t1.max(t2);
    if f1_planar && f2_planar {
        if let Some(angle) = plane_angle {
            let dt = compute_int_range(t1, t2, angle);
            return (dt * dt + t1 * t1).sqrt().max(t_max);
        }
    }
    t_max
}

/// Sampled bounding box of a 3D curve over its stored range, with the
/// PerformFF expand (`aTolFF + MaxVertexTol`) applied afterwards by the
/// caller via [`BndBox::enlarge`].
pub fn curve_box_unexpanded(curve: &dyn Curve, first: f64, last: f64, curve_tol: f64) -> BndBox {
    let mut box_ = BndBox::new();
    add_curve_to_box(curve, first, last, curve_tol, &mut box_);
    box_
}


/// True when `p` is inside `box` expanded by `tol` (used by IsExistingVertex).
pub fn point_box_out(box_: &BndBox, p: &GpPnt, tol: f64) -> bool {
    let mut pb = BndBox::new();
    pb.add_point(p);
    pb.enlarge(tol);
    box_.is_out_box(&pb)
}

/// Enlarge `box` by `value` (OCCT `Bnd_Box::Enlarge`).
pub fn enlarge_box(box_: &mut BndBox, value: f64) {
    box_.enlarge(value);
}

/// Combined enlarge used by PerformFF after CheckCurve succeeds:
/// `aBox.Enlarge(aBoxExpandValue)` where `aBoxExpandValue = aTolFF + MaxVertexTol`.
pub fn enlarge_section_curve_box(box_: &mut BndBox, tol_ff: f64, max_vertex_tol: f64) {
    box_.enlarge(tol_ff + max_vertex_tol);
}

/// Diagnostic dump of a box's raw extents (gap subtracted).
pub fn box_raw_extents(box_: &BndBox) -> Option<(f64, f64, f64)> {
    if box_.is_void() {
        return None;
    }
    let (xmin, xmax, ymin, ymax, zmin, zmax) = box_.get()?;
    let g = 2.0 * box_.gap();
    Some((
        (xmax - xmin - g).max(0.0),
        (ymax - ymin - g).max(0.0),
        (zmax - zmin - g).max(0.0),
    ))
}
