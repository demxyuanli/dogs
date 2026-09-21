//! Port of `GeomBndLib_Curve.cxx` for 3D curves (`Box` only).
//!
//! Dispatch follows `GeomBndLib_Curve.cxx:101-172` (the handle constructor):
//! Line / Circle / Ellipse / Hyperbola / Parabola / Bezier / BSpline / Offset,
//! everything else `GeomBndLib_OtherCurve`.
//!
//! **Ported conic arms (batch 87, task T-12 first half)**: the Ellipse,
//! Hyperbola and Parabola boxes (`GeomBndLib_Ellipse.cxx:23-114`,
//! `GeomBndLib_Hyperbola.cxx:25-140`, `GeomBndLib_Parabola.cxx:23-94`) now run
//! their analytic extrema instead of the `OtherCurve` sampling path — the
//! `gp_Elips` / `gp_Hypr` / `gp_Parab` queries they need exist on
//! `occt_geom::Curve` since the A20/A29 work.
//!
//! **Ported (batch 90, task T-12 second half)**: the periodic arm of
//! `GeomBndLib_BSplineCurve::Box` (`BSplineCurve.cxx:43-50`):
//! `ElCLib::AdjustPeriodic` puts the range inside one period and the segmented
//! curve's knots are the original distinct knots rotated by whole periods, so
//! the span loop below reproduces OCCT's `Segment` + knot-span `FillBox`
//! without a `Segment` port — see [`box_bspline`].

use occt_core::bnd::BndBox;
use occt_core::bspl::knots as bspl_knots;
use occt_core::elib::clib::{ellipse_value, hyperbola_value, parabola_value};
use occt_core::gp::{GpElips, GpHypr, GpLin, GpParab, GpPnt};
use occt_core::precision::{epsilon, Precision, PCONFUSION};

use occt_geom::curve::Curve;

use crate::geom_bnd_lib_circle3d::box_circ_range;
use crate::geom_bnd_lib_elclib2d::{adjust_periodic, in_period};
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

/// `GeomBndLib_Ellipse::Box(theElips, theTol)` (`GeomBndLib_Ellipse.cxx:23-45`).
pub fn box_ellipse_full(the_elips: &GpElips, the_tol: f64) -> BndBox {
    let mut a_box = BndBox::new();
    let a_maj_r = the_elips.major_radius();
    let a_min_r = the_elips.minor_radius();
    let a_loc = the_elips.location();
    let a_o = a_loc.xyz();
    let a_xd = the_elips.pos.x_direction().xyz();
    let a_yd = the_elips.pos.y_direction().xyz();

    // Full ellipse: per-coordinate analytic extrema
    // `Amp = sqrt(Major²·Xd_k² + Minor²·Yd_k²)` (`cxx:38`).
    let mut a_min = [0.0_f64; 3];
    let mut a_max = [0.0_f64; 3];
    for k in 0..3 {
        let a_xk = a_xd.coord(k);
        let a_yk = a_yd.coord(k);
        let a_amp = (a_maj_r * a_maj_r * a_xk * a_xk + a_min_r * a_min_r * a_yk * a_yk).sqrt();
        a_min[k] = a_o.coord(k) - a_amp;
        a_max[k] = a_o.coord(k) + a_amp;
    }
    a_box.update(a_min[0], a_min[1], a_min[2], a_max[0], a_max[1], a_max[2]);
    a_box.enlarge(the_tol);
    a_box
}

/// `GeomBndLib_Ellipse::Box(theElips, theU1, theU2, theTol)`
/// (`GeomBndLib_Ellipse.cxx:49-114`).
pub fn box_ellipse_range(the_elips: &GpElips, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox {
    let a_period = 2.0 * std::f64::consts::PI - PCONFUSION;
    if the_u2 - the_u1 >= a_period {
        return box_ellipse_full(the_elips, the_tol);
    }

    let a_maj_r = the_elips.major_radius();
    let a_min_r = the_elips.minor_radius();
    let a_loc = the_elips.location();
    let a_o = a_loc.xyz();
    let a_xd = the_elips.pos.x_direction().xyz();
    let a_yd = the_elips.pos.y_direction().xyz();

    let mut a_box = BndBox::new();
    let mut a_u1 = the_u1;
    let mut a_u2 = the_u2;
    let a_tol = epsilon(1.0);
    adjust_periodic(0.0, 2.0 * std::f64::consts::PI, a_tol, &mut a_u1, &mut a_u2);

    // Arc endpoints (`cxx:71-72`).
    a_box.add_point(&ellipse_value(the_elips, a_u1));
    a_box.add_point(&ellipse_value(the_elips, a_u2));

    for k in 0..3 {
        let a_xk = a_xd.coord(k);
        let a_yk = a_yd.coord(k);

        let mut a_t_extr_min = if a_xk.abs() > occt_core::precision::RESOLUTION {
            let t = ((a_min_r * a_yk) / (a_maj_r * a_xk)).atan();
            in_period(t, 0.0, 2.0 * std::f64::consts::PI)
        } else {
            std::f64::consts::PI / 2.0
        };
        let mut a_t_extr_max = if a_t_extr_min <= std::f64::consts::PI {
            a_t_extr_min + std::f64::consts::PI
        } else {
            a_t_extr_min - std::f64::consts::PI
        };

        let a_val_min = a_maj_r * a_t_extr_min.cos() * a_xk
            + a_min_r * a_t_extr_min.sin() * a_yk
            + a_o.coord(k);
        let a_val_max = a_maj_r * a_t_extr_max.cos() * a_xk
            + a_min_r * a_t_extr_max.sin() * a_yk
            + a_o.coord(k);
        if a_val_min > a_val_max {
            std::mem::swap(&mut a_t_extr_min, &mut a_t_extr_max);
        }

        let a_tk = in_period(a_t_extr_min, a_u1, a_u1 + 2.0 * std::f64::consts::PI);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            a_box.add_point(&ellipse_value(the_elips, a_t_extr_min));
        }
        let a_tk = in_period(a_t_extr_max, a_u1, a_u1 + 2.0 * std::f64::consts::PI);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            a_box.add_point(&ellipse_value(the_elips, a_t_extr_max));
        }
    }

    a_box.enlarge(the_tol);
    a_box
}

/// `computeHyperbolaBox` (`GeomBndLib_Hyperbola.cxx:25-69`).
fn compute_hyperbola_box(the_hypr: &GpHypr, the_t1: f64, the_t2: f64, the_box: &mut BndBox) {
    let a_p1 = hyperbola_value(the_hypr, the_t1);
    let a_p2 = hyperbola_value(the_hypr, the_t2);
    the_box.add_point(&a_p1);
    the_box.add_point(&a_p2);

    if the_t1 * the_t2 < 0.0 {
        the_box.add_point(&hyperbola_value(the_hypr, 0.0));
    }

    let a_x_dir = the_hypr.pos.x_direction().xyz();
    let a_y_dir = the_hypr.pos.y_direction().xyz();
    let a_r_maj = the_hypr.major_radius;
    let a_r_min = the_hypr.minor_radius;
    let a_eps = epsilon(1.0);

    for i in 0..3 {
        let a_a = a_r_min * a_y_dir.coord(i);
        let a_b = a_r_maj * a_x_dir.coord(i);

        let a_abp = (a_a + a_b).abs();
        let a_bam = (a_b - a_a).abs();

        // A coordinate whose extremal equation degenerates has no interior
        // extremum (`cxx:55-58`).
        if a_abp < a_eps || a_bam < a_eps {
            continue;
        }

        let a_cf = a_bam / a_abp;
        let a_t3 = 0.5 * a_cf.ln();

        if a_t3 < the_t1 || a_t3 > the_t2 {
            continue;
        }
        the_box.add_point(&hyperbola_value(the_hypr, a_t3));
    }
}

/// `GeomBndLib_Hyperbola::Box(theHypr, theU1, theU2, theTol)`
/// (`GeomBndLib_Hyperbola.cxx:75-140`).
pub fn box_hyperbola_range(the_hypr: &GpHypr, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox {
    let mut a_box = BndBox::new();
    if Precision::is_negative_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            // `cxx:82` throws `Standard_Failure`; the port returns the void box
            // the caller's `IsVoid` check already handles.
            return a_box;
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
            a_box.open_zmax();
        } else {
            a_box.add_point(&hyperbola_value(the_hypr, the_u2));
        }
        a_box.open_xmin();
        a_box.open_ymin();
        a_box.open_zmin();
    } else if Precision::is_positive_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
            a_box.open_zmin();
        } else if Precision::is_positive_infinite(the_u2) {
            // `cxx:108` throws; see above.
            return a_box;
        } else {
            a_box.add_point(&hyperbola_value(the_hypr, the_u2));
        }
        a_box.open_xmax();
        a_box.open_ymax();
        a_box.open_zmax();
    } else {
        a_box.add_point(&hyperbola_value(the_hypr, the_u1));
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
            a_box.open_zmin();
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
            a_box.open_zmax();
        } else {
            compute_hyperbola_box(the_hypr, the_u1, the_u2, &mut a_box);
        }
    }
    a_box.enlarge(the_tol);
    a_box
}

/// `GeomBndLib_Parabola::Box(theParab, theU1, theU2, theTol)`
/// (`GeomBndLib_Parabola.cxx:23-94`).
pub fn box_parabola_range(the_parab: &GpParab, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox {
    let mut a_box = BndBox::new();
    if Precision::is_negative_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            // `cxx:33` throws; see `box_hyperbola_range`.
            return a_box;
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
            a_box.open_zmax();
        } else {
            a_box.add_point(&parabola_value(the_parab, the_u2));
        }
        a_box.open_xmin();
        a_box.open_ymin();
        a_box.open_zmin();
    } else if Precision::is_positive_infinite(the_u1) {
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
            a_box.open_zmin();
        } else if Precision::is_positive_infinite(the_u2) {
            // `cxx:59` throws; see above.
            return a_box;
        } else {
            a_box.add_point(&parabola_value(the_parab, the_u2));
        }
        a_box.open_xmax();
        a_box.open_ymax();
        a_box.open_zmax();
    } else {
        a_box.add_point(&parabola_value(the_parab, the_u1));
        if Precision::is_negative_infinite(the_u2) {
            a_box.open_xmin();
            a_box.open_ymin();
            a_box.open_zmin();
        } else if Precision::is_positive_infinite(the_u2) {
            a_box.open_xmax();
            a_box.open_ymax();
            a_box.open_zmax();
        } else {
            a_box.add_point(&parabola_value(the_parab, the_u2));
            if the_u1 * the_u2 < 0.0 {
                a_box.add_point(&parabola_value(the_parab, 0.0));
            }
        }
    }
    a_box.enlarge(the_tol);
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
/// Full range and sub-range share the knot-span `FillBox` loop (`cxx:85-105`)
/// and `ReduceSplineBox(myGeom->Poles(), ...)` (`cxx:110`, the ORIGINAL full pole
/// array, not a window). A strict sub-range is segmented by OCCT (`cxx:39-85`)
/// merely to re-parameterize/truncate the knot vector:
///
/// * non-periodic: `cxx:51-60` clamps the range to the curve's own bounds and
///   `Segment` keeps the parameterization, so clipping the spans to
///   `[a_u1, a_u2]` reproduces the segmented curve's samples;
/// * periodic (`cxx:43-50`, batch 90 / task T-12 second half):
///   `ElCLib::AdjustPeriodic` first moves the range inside one period, then
///   `Segment` rotates the knot vector by whole periods and truncates it, so the
///   segmented curve's distinct knots inside the range are the original distinct
///   knots shifted by multiples of the period — the span list below.
fn box_bspline(curve: &dyn Curve, poles: &[GpPnt], knots: &[f64], u1: f64, u2: f64, tol: f64) -> BndBox {
    let degree = curve.nurbs_degree().unwrap_or(1) as i32;
    let (a_u1, a_u2, cuts) = if curve.is_periodic() {
        let first = curve.first_parameter();
        let last = curve.last_parameter();
        let mut a_u1 = u1;
        let mut a_u2 = u2;
        adjust_periodic(first, last, PCONFUSION, &mut a_u1, &mut a_u2);
        let period = last - first;
        let (uknots, _umults) = bspl_knots::unique_knots_mults(knots);
        let mut cuts: Vec<f64> = Vec::new();
        if period > PCONFUSION && uknots.len() >= 2 {
            // Knots of the segmented curve: the original distinct knots rotated
            // by whole periods to cover `[a_u1, a_u2]`.
            let m_lo = ((a_u1 - uknots[uknots.len() - 1]) / period).floor() - 1.0;
            let m_hi = ((a_u2 - uknots[0]) / period).ceil() + 1.0;
            let mut m = m_lo;
            while m <= m_hi {
                for k in &uknots {
                    cuts.push(k + m * period);
                }
                m += 1.0;
            }
            cuts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        }
        (a_u1, a_u2, cuts)
    } else {
        let a_u1 = u1.max(curve.first_parameter());
        let a_u2 = u2.min(curve.last_parameter());
        let (uknots, _umults) = bspl_knots::unique_knots_mults(knots);
        (a_u1, a_u2, uknots)
    };

    if cuts.is_empty() && !curve.is_periodic() {
        let mut a_box = BndBox::new();
        a_box.add_point(&curve.d0(a_u1));
        a_box.add_point(&curve.d0(a_u2));
        a_box.enlarge(tol);
        return a_box;
    }

    // `cxx:92-104`: spans `[a_u1, k1] [k1, k2] ... [kn, a_u2]` over the
    // segmented curve's distinct knots; the interior knots outside the range are
    // dropped (the segment's end knots carry multiplicity `degree + 1`).
    let mut a_b1 = BndBox::new();
    let mut a_tol = 0.0_f64;
    let mut a_first = a_u1;
    let n = degree.max(1);
    for cut in cuts
        .iter()
        .copied()
        .filter(|k| *k > a_u1 + PCONFUSION && *k < a_u2 - PCONFUSION)
        .chain(std::iter::once(a_u2))
    {
        let a_last = cut.min(a_u2);
        if a_last > a_first + PCONFUSION {
            a_tol = a_tol.max(fill_box3d(&mut a_b1, curve, a_first, a_last, n));
        }
        a_first = a_last;
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
    // `GeomBndLib_Curve.cxx:98-147` dispatch order: Line, Circle, Ellipse,
    // Hyperbola, Parabola, Bezier, BSpline, Offset, Other. The three conic arms
    // were parked while `occt_geom::Curve` had no `gp_*` queries; those exist
    // since the A20/A29 work (task T-12 first half, batch 87).
    if let Some(elips) = curve.gp_ellipse() {
        return box_ellipse_range(&elips, u1, u2, tol);
    }
    if let Some(hypr) = curve.gp_hyperbola() {
        return box_hyperbola_range(&hypr, u1, u2, tol);
    }
    if let Some(parab) = curve.gp_parabola() {
        return box_parabola_range(&parab, u1, u2, tol);
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
