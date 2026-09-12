//! 2D spline bounding-box helpers.
//!
//! Source: `GeomBndLib_SplineHelpers.pxx` (`FillBox`, `ReduceSplineBox`,
//! `ComputePoleIndexRange`, `BezierCurveBox`, `BSplineCurveBox`). PerformAreas
//! UV boxes call `GeomBndLib_BezierCurve2d::Box` / `BSplineCurve2d::Box` (pole
//! hull after optional Segment). These helpers are the FillBox + pole-range
//! path used when the full unique-knot span is kept.

use occt_core::bnd::BndBox2d;
use occt_core::bspl::knots as bspl_knots;
use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::{Geom2dBezierCurve, Geom2dBSplineCurve};

use crate::geom_bnd_lib_elclib2d::adjust_periodic;
use crate::geom_bnd_lib_other2d::{box_other, fill_box2d};
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
/// `Multiplicities`). Hunt is 0-based here; OCCT's Hunt is 1-based with
/// `Lower() == 1`.
pub fn compute_pole_index_range(
    knots: &[f64],
    mults: &[i32],
    the_degree: i32,
    the_min: f64,
    the_max: f64,
    the_max_pole_idx: i32,
    the_is_periodic: bool,
) -> (i32, i32) {
    if knots.is_empty() || mults.is_empty() {
        return (1, the_max_pole_idx.max(1));
    }
    let lower = 1i32;
    let upper = knots.len() as i32;
    let mut out_min = bspl_knots::hunt(knots, the_min) as i32 + 1;
    if out_min < lower {
        out_min = lower;
    }
    if out_min > upper {
        out_min = upper;
    }
    let mut out_max = bspl_knots::hunt(knots, the_max) as i32 + 1;
    out_max += 1;
    if out_max < lower {
        out_max = lower;
    }
    if out_max > upper {
        out_max = upper;
    }
    let a_multiplier = {
        let i = (out_max - 1) as usize;
        if i < mults.len() {
            mults[i]
        } else {
            0
        }
    };
    out_min = bspl_knots::pole_index(the_degree, out_min, the_is_periodic, mults) + 1;
    if out_min < 1 {
        out_min = 1;
    }
    out_max = bspl_knots::pole_index(the_degree, out_max, the_is_periodic, mults) + 1;
    out_max += the_degree - a_multiplier;
    if !the_is_periodic && out_max > the_max_pole_idx {
        out_max = the_max_pole_idx;
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
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> BndBox2d {
    let first = 0.0;
    let last = 1.0;
    if the_u1 - first > PCONFUSION || last - the_u2 > PCONFUSION {
        return box_other(geom, the_u1.max(first), the_u2.min(last), the_tol);
    }
    let mut b1 = BndBox2d::new();
    let n = geom.degree().max(1) as i32;
    let a_tol = fill_box2d(&mut b1, geom, the_u1, the_u2, n);
    b1.enlarge(WEAKNESS * a_tol);
    let mut box_ = reduce_spline_box2d(&geom.poles, &b1);
    box_.enlarge(the_tol);
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
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> BndBox2d {
    let mut a_u1 = the_u1;
    let mut a_u2 = the_u2;
    adjust_periodic(
        geom.first_parameter(),
        geom.last_parameter(),
        PCONFUSION,
        &mut a_u1,
        &mut a_u2,
    );
    box_other(geom, a_u1, a_u2, the_tol)
}

/// Non-periodic `BSplineCurveBox` (`SplineHelpers.pxx:332-398`).
pub fn bspline_curve_box2d(
    geom: &Geom2dBSplineCurve,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> BndBox2d {
    let a_degree = geom.degree() as i32;
    let a_u1 = the_u1.max(geom.first_parameter());
    let a_u2 = the_u2.min(geom.last_parameter());
    if a_u1 - geom.first_parameter() > PCONFUSION || geom.last_parameter() - a_u2 > PCONFUSION {
        return box_other(geom, a_u1, a_u2, the_tol);
    }
    let (uknots, umults) = bspl_knots::unique_knots_mults(&geom.knots);
    if uknots.len() < 2 {
        let mut box_ = BndBox2d::new();
        box_.add_point(&geom.d0(a_u1));
        box_.add_point(&geom.d0(a_u2));
        box_.enlarge(the_tol);
        return box_;
    }
    let lower = 0usize;
    let upper = uknots.len() - 1;
    let mut a_k_min = bspl_knots::hunt(&uknots, a_u1);
    a_k_min = a_k_min.clamp(lower, upper.saturating_sub(1));
    let mut a_k_max = bspl_knots::hunt(&uknots, a_u2);
    a_k_max = (a_k_max + 1).clamp(lower, upper);

    let mut a_b1 = BndBox2d::new();
    let mut a_tol = 0.0_f64;
    let mut a_first = a_u1;
    let n = a_degree.max(1);
    for a_k in (a_k_min + 1)..=a_k_max {
        let a_last = if a_k < a_k_max {
            a_u2.min(uknots[a_k])
        } else {
            a_u2
        };
        if a_last > a_first + PCONFUSION {
            a_tol = a_tol.max(fill_box2d(&mut a_b1, geom, a_first, a_last, n));
        }
        a_first = a_last;
        if a_first >= a_u2 - PCONFUSION {
            break;
        }
    }
    if a_b1.is_void() {
        let mut box_ = BndBox2d::new();
        box_.add_point(&geom.d0(a_u1));
        box_.add_point(&geom.d0(a_u2));
        box_.enlarge(the_tol);
        return box_;
    }
    a_b1.enlarge(WEAKNESS * a_tol);
    let poles = bspline_poles(geom);
    let (pmin, pmax) = compute_pole_index_range(
        &uknots,
        &umults,
        a_degree,
        a_u1,
        a_u2,
        poles.len() as i32,
        false,
    );
    let mut box_ = reduce_spline_box2d_range(&poles, pmin, pmax, &a_b1);
    box_.enlarge(the_tol);
    box_
}

/// Sampling half of `CurveBoxOptimal` for 2D (`NDim = 2`) **without**
/// `AdjustExtrT` (PSO+Brent). PerformAreas never calls `BoxOptimal`; this is
/// the sample envelope that `BoxOptimal` starts from.
pub fn curve_box_sample2d(curve: &dyn Curve2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    let a_nb = compute_nb_samples2d(
        sample_kind_of(curve),
        8,
        3,
        4,
        curve.first_parameter(),
        curve.last_parameter(),
        the_u1,
        the_u2,
    )
    .max(2);
    let a_du = (the_u2 - the_u1) / (a_nb - 1) as f64;
    let a_du2 = a_du / 2.0;
    let mut coord_min = [f64::MAX; 2];
    let mut coord_max = [f64::MIN; 2];
    for j in 0..a_nb {
        let a_u = the_u1 + j as f64 * a_du;
        let a_p = curve.d0(a_u);
        coord_min[0] = coord_min[0].min(a_p.x());
        coord_min[1] = coord_min[1].min(a_p.y());
        coord_max[0] = coord_max[0].max(a_p.x());
        coord_max[1] = coord_max[1].max(a_p.y());
        if j > 0 {
            let a_p_mid = curve.d0(a_u - a_du2);
            coord_min[0] = coord_min[0].min(a_p_mid.x());
            coord_min[1] = coord_min[1].min(a_p_mid.y());
            coord_max[0] = coord_max[0].max(a_p_mid.x());
            coord_max[1] = coord_max[1].max(a_p_mid.y());
        }
    }
    let mut a_box = BndBox2d::new();
    a_box.add_point(&GpPnt2d::new(coord_min[0], coord_min[1]));
    a_box.add_point(&GpPnt2d::new(coord_max[0], coord_max[1]));
    a_box.enlarge(the_tol.max(CONFUSION));
    a_box
}
