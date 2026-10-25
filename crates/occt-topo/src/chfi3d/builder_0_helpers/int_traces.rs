//! `ChFi3d_IntTraces` — intersect the p-curves of two neighbouring fillet
//! traces to sharpen an end.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:1222-1447`
//! (OCCT 8.0.0): the file-static `recadre` (`cxx:1222-1231`) and
//! `ChFi3d_IntTraces` (`cxx:1234-1447`).
//!
//! Translation notes:
//! - `p1` / `p2` are out-parameters in OCCT and part of the answer here; the
//!   `bool` result becomes `Option<(f64, f64)>`,
//! - `Geom2dAdaptor_Curve C1; C1.Load(curve)` is the adaptor's own domain, and
//!   `C1.Load(curve, first, last)` the explicit one; the port hands
//!   [`Geom2dIntGInter`] the `IntRes2d_Domain` those loads produce
//!   (`domain_over_range` reproduces `ComputeDomain` over an explicit range),
//! - the two-periodic comparison `PCurveOnFace() == PCurveOnFace()` is a handle
//!   identity test in OCCT, so it is `Arc::ptr_eq` over the two interferences'
//!   **original** p-curves (the ones read before any `TrimmedCurve`
//!   unwrapping),
//! - `gp_Pnt2d` / `IntRes2d_IntersectionPoint` are values here, so the loop
//!   rebinds them instead of assigning through references.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::intres2d::intersection::IntRes2dIntersectionPoint;
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::chfi_ds::surf_data::SurfData;

/// `recadre(p, ref, sens, first, last)` (`cxx:1222-1231`): shift `p` by one
/// period when that lands closer to `ref`.
fn recadre(p: f64, r#ref: f64, sens: i32, first: f64, last: f64) -> f64 {
    let pp = p + if sens > 0 { first - last } else { last - first };
    if (pp - r#ref).abs() < (p - r#ref).abs() {
        pp
    } else {
        p
    }
}

/// The p-curve of one interference and the domain the `Load` arm of
/// `cxx:1266-1314` gives it.
struct LoadedCurve {
    curve: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
}

/// `cxx:1266-1284` / `cxx:1297-1314`: unwrap a `Geom2d_TrimmedCurve` when the
/// p-curve is periodic, then record the domain the adaptor would use.
fn load(
    pcurve: Option<&Arc<dyn Curve2d>>,
    first: f64,
    last: f64,
    delta: f64,
    periodic: bool,
) -> Option<(LoadedCurve, f64, f64)> {
    let pcurve = pcurve?;
    let curve = if periodic {
        // `occ::down_cast<Geom2d_TrimmedCurve>(pcf1)` then
        // `pcf1 = tr1->BasisCurve()`.
        match pcurve.trimmed_basis() {
            Some(basis) => Arc::from(basis.clone_dyn()),
            None => Arc::clone(pcurve),
        }
    } else {
        Arc::clone(pcurve)
    };
    let first1 = curve.first_parameter();
    let last1 = curve.last_parameter();
    let (d0, d1) = if periodic {
        (first1, last1)
    } else {
        (first - delta, last + delta)
    };
    Some((
        LoadedCurve {
            curve,
            first: first1,
            last: last1,
        },
        d0,
        d1,
    ))
}

/// `ChFi3d_IntTraces(fd1, pref1, p1, jf1, sens1, fd2, pref2, p2, jf2, sens2,
/// RefP2d, Check2dDistance, enlarge)` (`cxx:1234-1447`).
///
/// `jf` selects the interference (`1` -> `S1`, anything else `S2`), `sens` the
/// direction criterion (`+1` / `-1`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi3d_int_traces(
    fd1: &SurfData,
    pref1: f64,
    jf1: i32,
    sens1: i32,
    fd2: &SurfData,
    pref2: f64,
    jf2: i32,
    sens2: i32,
    ref_p2d: &GpPnt2d,
    check_2d_distance: bool,
    enlarge: bool,
) -> Option<(f64, f64)> {
    // `cxx:1252-1256`.
    let mut first = fd1.interference(jf1).first_parameter();
    let mut last = fd1.interference(jf1).last_parameter();
    if (last - first) < PCONFUSION {
        return None;
    }
    // `cxx:1261-1264`.
    let mut delta = 0.0;
    if enlarge {
        delta = 0.1f64.min(0.05 * (last - first));
    }
    // `cxx:1266-1270`: read once, test for null.
    let pcf1 = fd1.interference(jf1).p_curve_on_face()?;
    let isper1 = pcf1.is_periodic();
    let (c1, c1_d0, c1_d1) = load(Some(pcf1), first, last, delta, isper1)?;

    // `cxx:1286-1294`.
    first = fd2.interference(jf2).first_parameter();
    last = fd2.interference(jf2).last_parameter();
    if (last - first) < PCONFUSION {
        return None;
    }
    if enlarge {
        delta = 0.1f64.min(0.05 * (last - first));
    }
    let pcf2 = fd2.interference(jf2).p_curve_on_face()?;
    let isper2 = pcf2.is_periodic();
    // `cxx:1313`: the non-periodic arm re-reads the interference's p-curve.
    let pcf2_loaded = fd2.interference(jf2).p_curve_on_face();
    let (c2, c2_d0, c2_d1) = load(pcf2_loaded, first, last, delta, isper2)?;

    // `cxx:1319-1328`: `TolConf = Tol = Precision::PIntersection()`.
    let tol = occt_core::precision::PINTERSECTION;
    let mut inter = Geom2dIntGInter::new();
    let d1 = Geom2dIntGInter::domain_over_range(c1.curve.as_ref(), c1_d0, c1_d1, tol);
    if Arc::ptr_eq(pcf1, pcf2) {
        // `Intersection.Perform(C1, TolConf, Tol)` (`cxx:1322`).
        inter.perform_curve(c1.curve.as_ref(), &d1, tol, tol);
    } else {
        // `Intersection.Perform(C1, C2, TolConf, Tol)` (`cxx:1326`).
        let d2 = Geom2dIntGInter::domain_over_range(c2.curve.as_ref(), c2_d0, c2_d1, tol);
        inter.perform(
            c1.curve.as_ref(),
            &d1,
            c2.curve.as_ref(),
            &d2,
            tol,
            tol,
        );
    }
    // `cxx:1329-1341`.
    if !inter.is_done() {
        return None;
    }
    if inter.nb_points() < 1 {
        return None;
    }

    // `cxx:1346-1357`.
    let first_int: IntRes2dIntersectionPoint = inter.point(1);
    let mut p2d = *first_int.value();
    let mut p1 = first_int.param_on_first();
    let mut p2 = first_int.param_on_second();
    if isper1 {
        p1 = recadre(p1, pref1, sens1, c1.first, c1.last);
    }
    if isper2 {
        p2 = recadre(p2, pref2, sens2, c2.first, c2.last);
    }

    // `cxx:1359-1433`.
    for i in 2..=inter.nb_points() {
        let int2d = inter.point(i);
        if isper1 {
            let pp1 = recadre(int2d.param_on_first(), pref1, sens1, c1.first, c1.last);
            if (pp1 - pref1).abs() < (p1 - pref1).abs() {
                p1 = pp1;
                p2 = int2d.param_on_second();
                p2d = *int2d.value();
            } else if check_2d_distance
                && ref_p2d.distance(int2d.value()) < ref_p2d.distance(&p2d)
            {
                let mut pp2 = int2d.param_on_second();
                if isper2 {
                    pp2 = recadre(pp2, pref2, sens2, c2.first, c2.last);
                }
                p1 = pp1;
                p2 = pp2;
                p2d = *int2d.value();
            }
        } else if isper2 {
            let pp2 = recadre(int2d.param_on_second(), pref2, sens2, c2.first, c2.last);
            if (pp2 - pref2).abs() < (p2 - pref2).abs() {
                p2 = pp2;
                p1 = int2d.param_on_first();
                p2d = *int2d.value();
            } else if check_2d_distance
                && ref_p2d.distance(int2d.value()) < ref_p2d.distance(&p2d)
            {
                let mut pp1 = int2d.param_on_first();
                if isper1 {
                    pp1 = recadre(pp1, pref1, sens1, c1.first, c1.last);
                }
                p1 = pp1;
                p2 = pp2;
                p2d = *int2d.value();
            }
        } else if ((int2d.param_on_first() - p1) * f64::from(sens1) < 0.0)
            && ((int2d.param_on_second() - p2) * f64::from(sens2) < 0.0)
        {
            p1 = int2d.param_on_first();
            p2 = int2d.param_on_second();
            p2d = *int2d.value();
        } else if ((int2d.param_on_first() - pref1).abs() < (p1 - pref1).abs())
            && ((int2d.param_on_second() - pref2).abs() < (p2 - pref2).abs())
        {
            p1 = int2d.param_on_first();
            p2 = int2d.param_on_second();
            p2d = *int2d.value();
        } else if check_2d_distance
            && ref_p2d.distance(int2d.value()) < ref_p2d.distance(&p2d)
        {
            p1 = int2d.param_on_first();
            p2 = int2d.param_on_second();
            p2d = *int2d.value();
        }
    }
    let _ = (&c1, &c2);
    Some((p1, p2))
}
