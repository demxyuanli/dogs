//! 2D spline bounding-box helpers.
//!
//! Source: `GeomBndLib_SplineHelpers.pxx` (`FillBox`, `ReduceSplineBox`,
//! `ComputePoleIndexRange`, `BezierCurveBox`, `BSplineCurveBox`). PerformAreas
//! UV boxes call `GeomBndLib_BezierCurve2d::Box` / `BSplineCurve2d::Box` (pole
//! hull after optional Segment). These helpers are the FillBox + pole-range
//! path used when the full unique-knot span is kept.

use occt_core::bnd::BndBox2d;
use occt_core::bspl::knots as bspl_knots;
use occt_core::bspl::locate::hunt_occt;
use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::{Geom2dBezierCurve, Geom2dBSplineCurve};

use crate::geom_bnd_lib_elclib2d::adjust_periodic;
use crate::geom_bnd_lib_other_curve2d::{box_other_2d, fill_box2d};
use crate::geom_bnd_lib_sample2d::{compute_nb_samples2d, sample_kind_of};

const WEAKNESS: f64 = 1.5;

/// `PointOps<gp_Pnt2d>::Mid`.
pub fn mid2d(p1: &GpPnt2d, p2: &GpPnt2d) -> GpPnt2d {
    GpPnt2d::new(0.5 * (p1.x() + p2.x()), 0.5 * (p1.y() + p2.y()))
}

/// `GeomBndLib_SplineHelpers::ReduceSplineBox` 2D — intersect sampled box with
/// the poles' convex hull. `Bnd_Box2d::Get` returns `(xmin, ymin, xmax, ymax)`
/// in this port (same numbers as OCCT's `(xmin, xmax, ymin, ymax)` pair, just
/// a different tuple order).
pub fn reduce_spline_box2d(poles: &[GpPnt2d], orig: &BndBox2d) -> BndBox2d {
    let mut poles_box = BndBox2d::new();
    for p in poles {
        poles_box.add_point(p);
    }
    let Some((xmin, ymin, xmax, ymax)) = orig.get() else {
        return BndBox2d::new();
    };
    let mut reduced = BndBox2d::new();
    if !poles_box.is_void() {
        if let Some((pxmin, pymin, pxmax, pymax)) = poles_box.get() {
            reduced.update(
                xmin.max(pxmin),
                ymin.max(pymin),
                xmax.min(pxmax),
                ymax.min(pymax),
            );
        }
    } else {
        reduced.update(xmin, ymin, xmax, ymax);
    }
    reduced
}

/// `ReduceSplineBox` with a selected pole index range `[min_idx, max_idx]`
/// (1-based like OCCT). Indices beyond `poles.len()` wrap (periodic).
pub fn reduce_spline_box2d_range(
    poles: &[GpPnt2d],
    min_idx: i32,
    max_idx: i32,
    orig: &BndBox2d,
) -> BndBox2d {
    let mut poles_box = BndBox2d::new();
    let n = poles.len() as i32;
    if n <= 0 {
        return reduce_spline_box2d(&[], orig);
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
    let Some((xmin, ymin, xmax, ymax)) = orig.get() else {
        return BndBox2d::new();
    };
    let mut reduced = BndBox2d::new();
    if !poles_box.is_void() {
        if let Some((pxmin, pymin, pxmax, pymax)) = poles_box.get() {
            reduced.update(
                xmin.max(pxmin),
                ymin.max(pymin),
                xmax.min(pxmax),
                ymax.min(pymax),
            );
        }
    } else {
        reduced.update(xmin, ymin, xmax, ymax);
    }
    reduced
}

/// `GeomBndLib_SplineHelpers::ComputePoleIndexRange`.
///
/// `knots` / `mults` are the unique-knot tables (`Geom_BSplineCurve::Knots` /
/// `Multiplicities`). `hunt_occt` is the exact `BSplCLib::Hunt` port, 1-based
/// with `Lower() == 1`; the previous 0-based `bspl::knots::hunt` (last knot
/// `<= x`) plus `+ 1` was off by one whenever `min` / `the_max` sat exactly
/// on a knot, which widened or dropped a pole row.
pub fn compute_pole_index_range(
    knots: &[f64],
    mults: &[i32],
    degree: i32,
    min: f64,
    max_par: f64,
    max_pole_idx: i32,
    is_periodic: bool,
) -> (i32, i32) {
    if knots.is_empty() || mults.is_empty() {
        return (1, max_pole_idx.max(1));
    }
    let lower = 1i32;
    let upper = knots.len() as i32;
    let mut out_min = hunt_occt(knots, min);
    if out_min < lower {
        out_min = lower;
    }
    if out_min > upper {
        out_min = upper;
    }
    let mut out_max = hunt_occt(knots, max_par);
    out_max += 1;
    if out_max < lower {
        out_max = lower;
    }
    if out_max > upper {
        out_max = upper;
    }
    let multiplier = {
        let i = (out_max - 1) as usize;
        if i < mults.len() {
            mults[i]
        } else {
            0
        }
    };
    out_min = bspl_knots::pole_index(degree, out_min, is_periodic, mults) + 1;
    if out_min < 1 {
        out_min = 1;
    }
    out_max = bspl_knots::pole_index(degree, out_max, is_periodic, mults) + 1;
    out_max += degree - multiplier;
    if !is_periodic && out_max > max_pole_idx {
        out_max = max_pole_idx;
    }
    (out_min, out_max)
}

/// `GeomBndLib_SplineHelpers::BezierCurveBox`.
///
/// Sub-range: OCCT `Copy`+`Segment` then pole hull. Without Segment this
/// port uses `OtherCurve2d::Box`. Full range: `FillBox` with `N = Degree`
/// then `ReduceSplineBox`.
pub fn bezier_curve_box2d(
    geom: &Geom2dBezierCurve,
    u1: f64,
    u2: f64,
    tol: f64,
) -> BndBox2d {
    let first = 0.0;
    let last = 1.0;
    if u1 - first > PCONFUSION || last - u2 > PCONFUSION {
        return box_other_2d(geom, u1.max(first), u2.min(last), tol);
    }
    let mut b1 = BndBox2d::new();
    let n = geom.degree().max(1) as i32;
    let tol_val = fill_box2d(&mut b1, geom, u1, u2, n);
    b1.enlarge(WEAKNESS * tol_val);
    let mut box_ = reduce_spline_box2d(&geom.poles, &b1);
    box_.enlarge(tol);
    box_
}

fn bspline_poles(geom: &Geom2dBSplineCurve) -> Vec<GpPnt2d> {
    geom.xs
        .iter()
        .zip(geom.ys.iter())
        .map(|(&x, &y)| GpPnt2d::new(x, y))
        .collect()
}

/// Periodic branch of `BSplineCurveBox` (`SplineHelpers.pxx:279-329`).
///
/// OCCT `Copy`+`Segment` then FillBox per unique-knot span. This port has
/// no `Segment`; after `AdjustPeriodic` it falls back to `OtherCurve2d::Box`
/// (the same sampling `BndLib_Add2dCurve` uses for a non-adaptor curve).
pub fn bspline_curve_box2d_periodic(
    geom: &Geom2dBSplineCurve,
    u1: f64,
    u2: f64,
    tol: f64,
) -> BndBox2d {
    let mut u1_val = u1;
    let mut u2_val = u2;
    adjust_periodic(
        geom.first_parameter(),
        geom.last_parameter(),
        PCONFUSION,
        &mut u1_val,
        &mut u2_val,
    );
    box_other_2d(geom, u1_val, u2_val, tol)
}

/// Non-periodic `BSplineCurveBox` (`SplineHelpers.pxx:332-398`).
pub fn bspline_curve_box2d(
    geom: &Geom2dBSplineCurve,
    u1: f64,
    u2: f64,
    tol: f64,
) -> BndBox2d {
    let degree_val = geom.degree() as i32;
    let u1_val = u1.max(geom.first_parameter());
    let u2_val = u2.min(geom.last_parameter());
    if u1_val - geom.first_parameter() > PCONFUSION || geom.last_parameter() - u2_val > PCONFUSION {
        return box_other_2d(geom, u1_val, u2_val, tol);
    }
    let (uknots, umults) = bspl_knots::unique_knots_mults(&geom.knots);
    if uknots.len() < 2 {
        let mut box_ = BndBox2d::new();
        box_.add_point(&geom.d0(u1_val));
        box_.add_point(&geom.d0(u2_val));
        box_.enlarge(tol);
        return box_;
    }
    let lower = 0usize;
    let upper = uknots.len() - 1;
    let mut k_min = bspl_knots::hunt(&uknots, u1_val);
    k_min = k_min.clamp(lower, upper.saturating_sub(1));
    let mut k_max = bspl_knots::hunt(&uknots, u2_val);
    k_max = (k_max + 1).clamp(lower, upper);

    let mut b1 = BndBox2d::new();
    let mut tol_val = 0.0_f64;
    let mut first = u1_val;
    let n = degree_val.max(1);
    for k in (k_min + 1)..=k_max {
        let last = if k < k_max {
            u2_val.min(uknots[k])
        } else {
            u2_val
        };
        if last > first + PCONFUSION {
            tol_val = tol_val.max(fill_box2d(&mut b1, geom, first, last, n));
        }
        first = last;
        if first >= u2_val - PCONFUSION {
            break;
        }
    }
    if b1.is_void() {
        let mut box_ = BndBox2d::new();
        box_.add_point(&geom.d0(u1_val));
        box_.add_point(&geom.d0(u2_val));
        box_.enlarge(tol);
        return box_;
    }
    b1.enlarge(WEAKNESS * tol_val);
    let poles = bspline_poles(geom);
    let (pmin, pmax) = compute_pole_index_range(
        &uknots,
        &umults,
        degree_val,
        u1_val,
        u2_val,
        poles.len() as i32,
        false,
    );
    let mut box_ = reduce_spline_box2d_range(&poles, pmin, pmax, &b1);
    box_.enlarge(tol);
    box_
}

/// Sampling half of `CurveBoxOptimal` for 2D (`NDim = 2`) **without**
/// `AdjustExtrT` (PSO+Brent). PerformAreas never calls `BoxOptimal`; this is
/// the sample envelope that `BoxOptimal` starts from. The sample budget is
/// `GeomBndLib_SamplingHelpers::ComputeNbSamplesT` fed from the real curve
/// quantities (`CurveBoxOptimal`, `GeomBndLib_SplineHelpers.pxx:410-412`;
/// `ComputeNbSamples2d`, `GeomBndLib_SamplingHelpers.pxx:73-105`).
pub fn curve_box_sample2d(curve: &dyn Curve2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    let nb = compute_nb_samples2d(
        sample_kind_of(curve),
        curve.bezier_nb_poles().unwrap_or(0),
        curve.bspline_degree().unwrap_or(0),
        curve.bspline_nb_knots().unwrap_or(0),
        curve.first_parameter(),
        curve.last_parameter(),
        u1,
        u2,
    )
    .max(2);
    let du = (u2 - u1) / (nb - 1) as f64;
    let du2 = du / 2.0;
    let mut coord_min = [f64::MAX; 2];
    let mut coord_max = [f64::MIN; 2];
    for j in 0..nb {
        let u = u1 + j as f64 * du;
        let p = curve.d0(u);
        coord_min[0] = coord_min[0].min(p.x());
        coord_min[1] = coord_min[1].min(p.y());
        coord_max[0] = coord_max[0].max(p.x());
        coord_max[1] = coord_max[1].max(p.y());
        if j > 0 {
            let p_mid = curve.d0(u - du2);
            coord_min[0] = coord_min[0].min(p_mid.x());
            coord_min[1] = coord_min[1].min(p_mid.y());
            coord_max[0] = coord_max[0].max(p_mid.x());
            coord_max[1] = coord_max[1].max(p_mid.y());
        }
    }
    let mut bnd_box = BndBox2d::new();
    bnd_box.add_point(&GpPnt2d::new(coord_min[0], coord_min[1]));
    bnd_box.add_point(&GpPnt2d::new(coord_max[0], coord_max[1]));
    bnd_box.enlarge(tol.max(CONFUSION));
    bnd_box
}
