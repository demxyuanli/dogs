//! Port of `GeomBndLib_Curve.cxx` for 3D curves (`Box` only).
//!
//! Dispatch follows `GeomBndLib_Curve.cxx:101-172` (the handle constructor):
//! Line / Circle / Ellipse / Hyperbola / Parabola / Bezier / BSpline / Offset,
//! everything else `GeomBndLib_OtherCurve`.
//!
//! PARKED: Ellipse / Hyperbola / Parabola analytic boxes. The `occt_geom::Curve`
//! trait does not expose the underlying `gp_Elips` / `gp_Hypr` / `gp_Parab`
//! (`GeomBndLib_Ellipse.cxx`, `GeomBndLib_Hyperbola.cxx`,
//! `GeomBndLib_Parabola.cxx`), so they take the `OtherCurve` sampling path
//! below. Also parked: the periodic arm of `GeomBndLib_BSplineCurve::Box`
//! (`BSplineCurve.cxx:44-49` + `:299-330`), which needs `AdjustPeriodic` plus
//! `Geom_BSplineCurve::Segment` (`Geom_BSplineCurve.cxx:527-660`, i.e.
//! `SetOrigin` / `SetNotPeriodic` reparameterization); `occt_geom::Curve` has
//! no `Segment`, so periodic B-splines sample instead. The NON-periodic
//! sub-range arm (`BSplineCurve.cxx:39-85`) is ported: for a non-periodic
//! curve `Segment` keeps the parameterization (the knot vector is truncated,
//! its `DU` shift stays 0), so sampling the original curve over the same knot
//! spans yields the identical point set without rebuilding a segment object.

use occt_core::bnd::BndBox;
use occt_core::bspl::knots as bspl_knots;
use occt_core::gp::{GpLin, GpPnt};
use occt_core::precision::{Precision, PCONFUSION};

use occt_geom::curve::Curve;

use crate::geom_bnd_lib_circle3d::box_circ_range;
use crate::geom_bnd_lib_elclib2d::adjust_periodic;
use crate::geom_bnd_lib_inf3d::{open_max, open_min, open_min_max};

const WEAKNESS: f64 = 1.5;

/// `GeomBndLib_OtherCurve.cxx` local `FillBox` (33 samples) sampled envelope.
fn fill_box3d(box_: &mut BndBox, curve: &dyn Curve, first: f64, last: f64, n: i32) -> f64 {
    let mut a_p1 = curve.d0(first);
    box_.add_point(&a_p1);
    let mut tol = 0.0_f64;
    let mut p = first;
    let dp = last - first;
    if dp.abs() > PCONFUSION {
        let step = dp / (2 * n) as f64;
        for _ in 1..=n {
            p += step;
            let a_p2 = curve.d0(p);
            box_.add_point(&a_p2);
            p += step;
            let a_p3 = curve.d0(p);
            box_.add_point(&a_p3);
            let a_pc = mid(&a_p1, &a_p3);
            tol = tol.max(a_pc.distance(&a_p2));
            a_p1 = a_p3;
        }
    } else {
        box_.add_point(&curve.d0(last));
    }
    tol
}

fn mid(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

/// `GeomBndLib_SplineHelpers::ReduceSplineBox` 3D.
fn reduce_spline_box3d(poles: &[GpPnt], orig: &BndBox) -> BndBox {
    let mut poles_box = BndBox::new();
    for p in poles {
        poles_box.add_point(p);
    }
    let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = orig.get() else {
        return BndBox::new();
    };
    let mut reduced = BndBox::new();
    if !poles_box.is_void() {
        if let Some((pxmin, pxmax, pymin, pymax, pzmin, pzmax)) = poles_box.get() {
            reduced.update(
                xmin.max(pxmin),
                ymin.max(pymin),
                zmin.max(pzmin),
                xmax.min(pxmax),
                ymax.min(pymax),
                zmax.min(pzmax),
            );
        }
    } else {
        reduced.update(xmin, ymin, zmin, xmax, ymax, zmax);
    }
    reduced
}

/// `ReduceSplineBox` with a pole index range `[min_idx, max_idx]` (1-based).
/// Retained for the parked periodic B-spline arm.
#[allow(dead_code)]
fn reduce_spline_box3d_range(poles: &[GpPnt], min_idx: i32, max_idx: i32, orig: &BndBox) -> BndBox {
    let mut poles_box = BndBox::new();
    let n = poles.len() as i32;
    if n <= 0 {
        return reduce_spline_box3d(&[], orig);
    }
    let mut idx = min_idx;
    while idx <= max_idx {
        let mut i = idx;
        if i > n {
            i -= n;
        }
        if i >= 1 && (i as usize) <= poles.len() {
            poles_box.add_point(&poles[(i - 1) as usize]);
        }
        idx += 1;
    }
    let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = orig.get() else {
        return BndBox::new();
    };
    let mut reduced = BndBox::new();
    if !poles_box.is_void() {
        if let Some((pxmin, pxmax, pymin, pymax, pzmin, pzmax)) = poles_box.get() {
            reduced.update(
                xmin.max(pxmin),
                ymin.max(pymin),
                zmin.max(pzmin),
                xmax.min(pxmax),
                ymax.min(pymax),
                zmax.min(pzmax),
            );
        }
    } else {
        reduced.update(xmin, ymin, zmin, xmax, ymax, zmax);
    }
    reduced
}

/// `GeomBndLib_OtherCurve::Box` (`OtherCurve.cxx:72-87`).
pub fn box_other(curve: &dyn Curve, u1: f64, u2: f64, tol: f64) -> BndBox {
    let mut a_b1 = BndBox::new();
    let t = fill_box3d(&mut a_b1, curve, u1, u2, 33);
    a_b1.enlarge(WEAKNESS * t);
    let mut a_box = BndBox::new();
    if let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = a_b1.get() {
        a_box.update(xmin, ymin, zmin, xmax, ymax, zmax);
    }
    a_box.enlarge(tol);
    a_box
}

/// `GeomBndLib_Line::Box(gp_Lin, ...)` (`GeomBndLib_Line.hxx:66-121`).
pub fn box_lin(lin: &GpLin, u1: f64, u2: f64, tol: f64) -> BndBox {
    let mut a_box = BndBox::new();
    let dir = lin.direction();
    if Precision::is_negative_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            return a_box;
        } else if Precision::is_positive_infinite(u2) {
            open_min_max(&dir, &mut a_box);
            a_box.add_point(&occt_core::elib::clib::line_value(lin, 0.0));
        } else {
            open_min(&dir, &mut a_box);
            a_box.add_point(&occt_core::elib::clib::line_value(lin, u2));
        }
    } else if Precision::is_positive_infinite(u1) {
        if Precision::is_negative_infinite(u2) {
            open_min_max(&dir, &mut a_box);
            a_box.add_point(&occt_core::elib::clib::line_value(lin, 0.0));
        } else if Precision::is_positive_infinite(u2) {
            return a_box;
        } else {
            open_max(&dir, &mut a_box);
            a_box.add_point(&occt_core::elib::clib::line_value(lin, u2));
        }
    } else {
        a_box.add_point(&occt_core::elib::clib::line_value(lin, u1));
        if Precision::is_negative_infinite(u2) {
            open_min(&dir, &mut a_box);
        } else if Precision::is_positive_infinite(u2) {
            open_max(&dir, &mut a_box);
        } else {
            a_box.add_point(&occt_core::elib::clib::line_value(lin, u2));
        }
    }
    a_box.enlarge(tol);
    a_box
}

/// Reconstruct a `gp_Lin` from an unbounded line curve.
fn lin_from_curve(curve: &dyn Curve) -> Option<GpLin> {
    let p = curve.d0(0.0);
    let (_, t) = curve.d1(0.0);
    let d = occt_core::gp::GpDir::from_vec(&t).ok()?;
    Some(GpLin::from_pnt_dir(p, d))
}

/// `GeomBndLib_BezierCurve::Box` (`BezierCurve.cxx:29-46`).
fn box_bezier(curve: &dyn Curve, poles: &[GpPnt], u1: f64, u2: f64, tol: f64) -> BndBox {
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    if u1 - first > PCONFUSION || last - u2 > PCONFUSION {
        // PARK: `Copy`+`Segment` arm (`GeomBndLib_BezierCurve` / `SplineHelpers.pxx:236-253`).
        return box_other(curve, u1, u2, tol);
    }
    let degree = curve.nurbs_degree().unwrap_or(1) as i32;
    let mut a_sampled = BndBox::new();
    let defl = fill_box3d(&mut a_sampled, curve, u1, u2, degree.max(1));
    a_sampled.enlarge(WEAKNESS * defl);
    let mut a_box = reduce_spline_box3d(poles, &a_sampled);
    a_box.enlarge(tol);
    a_box
}

/// `GeomBndLib_BSplineCurve::Box` (`BSplineCurve.cxx:31-113`).
///
/// Non-periodic: full range and sub-range share the knot-span `FillBox` loop
/// (`cxx:85-105`) and `ReduceSplineBox(myGeom->Poles(), ...)` (`cxx:110`,
/// the ORIGINAL full pole array, not a window). A strict sub-range is
/// segmented by OCCT (`cxx:39-85`) merely to truncate the knot vector; for a
/// non-periodic curve that leaves the parameterization unchanged, so clipping
/// the spans to `[a_u1, a_u2]` reproduces the segmented curve's samples.
fn box_bspline(curve: &dyn Curve, poles: &[GpPnt], knots: &[f64], u1: f64, u2: f64, tol: f64) -> BndBox {
    let degree = curve.nurbs_degree().unwrap_or(1) as i32;
    if curve.is_periodic() {
        // PARK: periodic arm (`BSplineCurve.cxx:44-49`) needs `AdjustPeriodic`
        // plus `Geom_BSplineCurve::Segment` (`Geom_BSplineCurve.cxx:527-660`,
        // `SetOrigin` / `SetNotPeriodic`); `occt_geom::Curve` has no `Segment`.
        let mut a_u1 = u1;
        let mut a_u2 = u2;
        adjust_periodic(curve.first_parameter(), curve.last_parameter(), PCONFUSION, &mut a_u1, &mut a_u2);
        return box_other(curve, a_u1, a_u2, tol);
    }
    let a_u1 = u1.max(curve.first_parameter());
    let a_u2 = u2.min(curve.last_parameter());
    let (uknots, _umults) = bspl_knots::unique_knots_mults(knots);
    if uknots.len() < 2 {
        let mut a_box = BndBox::new();
        a_box.add_point(&curve.d0(a_u1));
        a_box.add_point(&curve.d0(a_u2));
        a_box.enlarge(tol);
        return a_box;
    }
    let lower = 0usize;
    let upper = uknots.len() - 1;
    let mut a_k_min = bspl_knots::hunt(&uknots, a_u1);
    a_k_min = a_k_min.clamp(lower, upper.saturating_sub(1));
    let mut a_k_max = bspl_knots::hunt(&uknots, a_u2);
    a_k_max = (a_k_max + 1).clamp(lower, upper);

    let mut a_b1 = BndBox::new();
    let mut a_tol = 0.0_f64;
    let mut a_first = a_u1;
    let n = degree.max(1);
    for a_k in (a_k_min + 1)..=a_k_max {
        let a_last = if a_k < a_k_max { a_u2.min(uknots[a_k]) } else { a_u2 };
        if a_last > a_first + PCONFUSION {
            a_tol = a_tol.max(fill_box3d(&mut a_b1, curve, a_first, a_last, n));
        }
        a_first = a_last;
        if a_first >= a_u2 - PCONFUSION {
            break;
        }
    }
    if a_b1.is_void() {
        let mut a_box = BndBox::new();
        a_box.add_point(&curve.d0(a_u1));
        a_box.add_point(&curve.d0(a_u2));
        a_box.enlarge(tol);
        return a_box;
    }
    a_b1.enlarge(WEAKNESS * a_tol);
    let mut a_box = reduce_spline_box3d(poles, &a_b1);
    a_box.enlarge(tol);
    a_box
}

/// `GeomBndLib_Curve::Box(theU1, theU2, theTol)` for a `Curve` handle.
pub fn box_curve(curve: &dyn Curve, u1: f64, u2: f64, tol: f64) -> BndBox {
    // `BRepAdaptor_Curve::Initialize` -> `GeomAdaptor_Curve::load`
    // (`GeomAdaptor_Curve.cxx:252-254`): a `Geom_TrimmedCurve` is replaced by
    // its basis curve; the adaptor keeps the caller's range, which our `[0, 1]`
    // trim remap transmits as basis-curve parameters.
    if let Some((basis, b1, b2)) = curve.untrimmed_basis() {
        let den = b2 - b1;
        return box_curve(basis.as_ref(), b1 + u1 * den, b1 + u2 * den, tol);
    }
    if curve.is_line() {
        if let Some(lin) = lin_from_curve(curve) {
            return box_lin(&lin, u1, u2, tol);
        }
    }
    if let Some(circ) = curve.gp_circ() {
        return box_circ_range(&circ, u1, u2, tol);
    }
    if let Some(poles) = curve.bezier_poles() {
        return box_bezier(curve, poles, u1, u2, tol);
    }
    if let Some(poles) = curve.bspline_poles() {
        if let Some(knots) = curve.bspline_knots() {
            return box_bspline(curve, poles, knots, u1, u2, tol);
        }
    }
    box_other(curve, u1, u2, tol)
}
