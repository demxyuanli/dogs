//! `GeomConvert::ConcatC1`, `GeomConvert::C0BSplineToC1BSplineCurve` and
//! `GeomConvert::C0BSplineToArrayOfC1BSplineCurve`, with the file-static
//! helpers they need.
//!
//! Source: `ModelingData/TKGeomBase/GeomConvert/GeomConvert.cxx` (OCCT 8.0.0):
//! * `GeomConvert_law_evaluator` `:463-492`;
//! * `MultNumandDenom` `:496-583`;
//! * `Pretreatment` `:587-606`;
//! * `NeedToBeTreated` `:612-632`;
//! * `Need2DegRepara` `:638-652`;
//! * `Indexmin` `:658-670`;
//! * `ReorderArrayOfG1Curves` `:676-728`;
//! * `GeomConvert_reparameterise_evaluator` `:734-760`;
//! * `ConcatG1` `:764-997`;
//! * `ConcatC1` 7-argument `:1022-1310` and the 6-argument overload
//!   `:1004-1019`;
//! * `C0BSplineToC1BSplineCurve` `:1313-1338`;
//! * `C0BSplineToArrayOfC1BSplineCurve` 3-argument `:1344-1354` and
//!   4-argument `:1356-1425`.
//!
//! The 2D twin of this file lives in `occt-geom2d::geom2d_convert_concat`. The
//! two `.cxx` files are not identical — the differences that matter are called
//! out inline where the port follows the 3D source:
//! * `C0BSplineToArrayOfC1BSplineCurve` counts the C1 spans from
//!   `FirstUKnotIndex()` (`cxx:1366`), the 2D version from `FirstUKnotIndex()+1`;
//! * the closed flag of that function is `BS->IsClosed() &&
//!   V1.IsParallel(V2, AngularTolerance)` (`cxx:1408-1413`), the 2D version
//!   tests the endpoint distance itself;
//! * `ConcatC1`'s closed branch calls `C.Add(Curve1, local_tolerance(j-1))`
//!   (`cxx:1232-1233`, i.e. `After = false`), the 2D version passes `true`;
//! * `Pretreatment` has no pole-count guard (`cxx:593-605`).
//!
//! Divergences, all in the failure channel only:
//! * `Standard_ConstructionError` / `Standard_DimensionError` and the `throw`
//!   of `Replacer` become `Err`.
//! * The two "uninitialised read" windows of the `.cxx` are bounded:
//!   `Pretreatment` skips a curve with fewer than two poles (OCCT reads
//!   `Weight(2)` past the end), and the `Segment` knot walk tests the stored
//!   multiplicity array bounds before the loop test that OCCT evaluates first
//!   (`cxx:1392`).

use crate::bspline_curve::GeomBSplineCurve;
use crate::convert_bspl::CompCurveToBSplineCurve;
use occt_core::convert::ParameterisationType;
use occt_geom2d::curve::Curve2d;
use crate::curve::Curve;
use crate::geom_lprop::{continuity, continuity_default};
use crate::hermit;
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::function_ops::{
    flat_knot_sequence, function_multiply, function_reparameterise, merge_bspline_knots,
};
use occt_core::bspl::knots::reparametrize;
use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;

/// `Geom2d_BSplineCurve::MaxDegree()` (`Geom2d_BSplineCurve.hxx`), read by
/// `ConcatC1`'s degree guard (`GeomConvert.cxx:1120`).
const MAX_DEGREE: usize = 25;

/// `GeomConvert_law_evaluator` (`GeomConvert.cxx:463-492`): the `Y` coordinate
/// of the 2D ancore curve when the derivative order is 0 and the parameter is
/// inside `[start_end[0], start_end[1]]`.
fn law_evaluator(
    ancore: &Geom2dBSplineCurve,
    der: i32,
    start_end: [f64; 2],
    param: f64,
) -> Result<f64, ()> {
    if der == 0 && param >= start_end[0] && param <= start_end[1] {
        Ok(ancore.d0(param).y())
    } else {
        Err(())
    }
}

/// `MultNumandDenom` (`GeomConvert.cxx:496-583`): the 3D curve whose
/// homogeneous poles are `BS`'s poles times the 2D law curve, then divided back
/// by the resulting weights. Used to turn a rational B-spline into one whose
/// weight function is the Hermite approximation.
pub fn mult_num_and_denom(
    a: &Geom2dBSplineCurve,
    bs: &GeomBSplineCurve,
) -> Result<GeomBSplineCurve, String> {
    let (bs_knots, bs_mults) = bs.distinct_knots_and_mults();
    if bs_knots.is_empty() {
        return Err("MultNumandDenom: empty BS knots".into());
    }
    let start_value = bs_knots[0];
    let end_value = bs_knots[bs_knots.len() - 1];
    // `cxx:510-518`: `Precision::PConfusion()`, clamped to a fifth of the span.
    let mut tolerance = PCONFUSION;
    if (end_value - start_value) / 5.0 < tolerance {
        tolerance = (end_value - start_value) / 5.0;
    }

    let (mut knots_val, mults) = a.distinct_knots_and_mults();
    reparametrize(bs.first_parameter(), bs.last_parameter(), &mut knots_val);
    let ancore = Geom2dBSplineCurve::from_flat(
        a.xs.clone(),
        a.ys.clone(),
        None,
        knot_sequence(&knots_val, &mults, a.degree() as i32),
        a.degree(),
        false,
    )
    .map_err(|e| format!("MultNumandDenom: ancore ({e})"))?;

    let (res_nb_poles, res_knots, res_mults) = merge_bspline_knots(
        tolerance,
        start_value,
        end_value,
        a.degree() as i32,
        &knots_val,
        &mults,
        bs.degree() as i32,
        &bs_knots,
        &bs_mults,
    );
    if res_nb_poles <= 0 {
        return Err("MultNumandDenom: MergeBSplineKnots".into());
    }
    let degree = (bs.degree() + a.degree()) as i32;
    let res_flat = flat_knot_sequence(&res_knots, &res_mults, degree);

    // `cxx:546-551`: homogeneous poles (`Pole * Weight`).
    let mut bs_num = Vec::with_capacity(bs.nb_poles() * 3);
    for i in 0..bs.nb_poles() {
        let w = bs.weight((i + 1) as i32);
        let p = bs.pole(i);
        bs_num.push(p.x() * w);
        bs_num.push(p.y() * w);
        bs_num.push(p.z() * w);
    }
    let bs_weights: Vec<f64> = (1..=bs.nb_poles()).map(|i| bs.weight(i as i32)).collect();

    let mut res_num = vec![0.0; res_nb_poles as usize * 3];
    let mut res_den = vec![0.0; res_nb_poles as usize];
    let eval = |der: i32, se: [f64; 2], p: f64| law_evaluator(&ancore, der, se, p);

    function_multiply(
        &eval,
        bs.degree() as i32,
        &bs.knots,
        &bs_num,
        3,
        &res_flat,
        degree,
        &mut res_num,
    )
    .map_err(|_| "MultNumandDenom: FunctionMultiply num".to_string())?;
    function_multiply(
        &eval,
        bs.degree() as i32,
        &bs.knots,
        &bs_weights,
        1,
        &res_flat,
        degree,
        &mut res_den,
    )
    .map_err(|_| "MultNumandDenom: FunctionMultiply den".to_string())?;

    // `cxx:577-582`: divide the homogeneous poles back by the weights.
    let mut poles = Vec::with_capacity(res_nb_poles as usize);
    for i in 0..res_nb_poles as usize {
        let d = res_den[i];
        if d == 0.0 {
            return Err("MultNumandDenom: zero denominator pole".into());
        }
        poles.push(GpPnt::new(
            res_num[i * 3] / d,
            res_num[i * 3 + 1] / d,
            res_num[i * 3 + 2] / d,
        ));
    }
    GeomBSplineCurve::rational(poles, res_den, res_flat, degree as usize)
        .map_err(|e| format!("MultNumandDenom: result curve ({e})"))
}

/// `Pretreatment` (`GeomConvert.cxx:587-606`): rescale the weights of a
/// rational curve whose four boundary weights share the first weight.
fn pretreatment(tab: &mut [GeomBSplineCurve]) {
    for curve in tab.iter_mut() {
        if !curve.is_rational() {
            continue;
        }
        let n = curve.nb_poles();
        // OCCT reads `Weight(2)` unconditionally (`cxx:596`); a one-pole curve
        // has no such weight, so the rescale is skipped here.
        if n < 2 {
            continue;
        }
        let a = curve.weight(1);
        if curve.weight(2) == a && curve.weight((n - 1) as i32) == a && curve.weight(n as i32) == a
        {
            for j in 1..=n {
                let w = curve.weight(j as i32);
                curve.set_weight(j - 1, w / a);
            }
        }
    }
}

/// `NeedToBeTreated` (`GeomConvert.cxx:612-632`).
fn need_to_be_treated(bs: &GeomBSplineCurve) -> bool {
    if !bs.is_rational() {
        return false;
    }
    let Some(weights) = bs.weights.as_deref() else {
        return false;
    };
    if !occt_core::bspl::rational::is_rational(weights) {
        return false;
    }
    let n = bs.nb_poles();
    if n < 2 {
        return false;
    }
    let out = |w: f64| w < 1.0 - CONFUSION || w > 1.0 + CONFUSION;
    out(bs.weight(1))
        || out(bs.weight(2))
        || out(bs.weight((n - 1) as i32))
        || out(bs.weight(n as i32))
}

/// `Need2DegRepara` (`GeomConvert.cxx:638-652`).
fn need_2deg_repara(tab: &[GeomBSplineCurve]) -> bool {
    let mut rapport = 1.0;
    for i in 0..tab.len().saturating_sub(1) {
        let (_, vec1) = tab[i + 1].d1(tab[i + 1].first_parameter());
        let (_, vec2) = tab[i].d1(tab[i].last_parameter());
        rapport = rapport * vec2.magnitude() / vec1.magnitude();
    }
    rapport > 1.0 + CONFUSION || rapport < 1.0 - CONFUSION
}

/// `Indexmin` (`GeomConvert.cxx:658-670`): the index of the lowest-degree
/// curve (the last one on ties).
fn indexmin(tab: &[GeomBSplineCurve]) -> usize {
    let mut index = 0usize;
    let mut degree = tab[0].degree();
    for (i, c) in tab.iter().enumerate() {
        if c.degree() <= degree {
            degree = c.degree();
            index = i;
        }
    }
    index
}

/// `ReorderArrayOfG1Curves` (`GeomConvert.cxx:676-728`): rotate the curve,
/// tolerance and G1 arrays so that the group boundary `StartIndex` becomes the
/// array start, and mark the closing junction with `ClosedTolerance`.
fn reorder_array_of_g1_curves(
    curves: &mut [GeomBSplineCurve],
    toler: &mut [f64],
    tab_g1: &mut [bool],
    start_index: usize,
    closed_tolerance: f64,
) {
    let n = curves.len();
    if n == 0 {
        return;
    }
    let bis_curves: Vec<GeomBSplineCurve> = curves.to_vec();
    let bis_toler: Vec<f64> = toler.to_vec();
    let bis_g1: Vec<bool> = tab_g1.to_vec();

    let pivot = n.saturating_sub(start_index + 2);
    for i in 0..=pivot {
        curves[i] = bis_curves[i + start_index + 1].clone();
        if i != pivot {
            toler[i] = bis_toler[i + start_index + 1];
            tab_g1[i] = bis_g1[i + start_index + 1];
        }
    }
    if pivot < toler.len() {
        toler[pivot] = closed_tolerance;
    }
    if pivot < tab_g1.len() {
        tab_g1[pivot] = true;
    }
    let base = n.saturating_sub(start_index + 1);
    for i in base..n {
        let src = i - base;
        curves[i] = bis_curves[src].clone();
        if i + 1 != n {
            toler[i] = bis_toler[src];
            tab_g1[i] = bis_g1[src];
        }
    }
}

/// `GeomConvert_reparameterise_evaluator::Evaluate` (`GeomConvert.cxx:742-756`):
/// `PLib::EvalPolynomial(u, der, 2, 1, coefficient[0], result)`.
fn reparameterise_evaluator(coefficient: &[f64; 3], der: i32, param: f64) -> Result<f64, ()> {
    if !(0..=2).contains(&der) {
        return Err(());
    }
    Ok(occt_core::bspl::plib::eval_polynomial_horner(coefficient, param))
}

/// `NeedToBeTreated` arm (`GeomConvert.cxx:1136-1142`, `:1281-1287`):
/// `MultNumandDenom(Hermit::Solution(curve), curve)`.
fn prepare_curve_for_concat(curve: &GeomBSplineCurve) -> Result<GeomBSplineCurve, String> {
    if need_to_be_treated(curve) {
        let law = hermit::solution_default(curve)?;
        mult_num_and_denom(&law, curve)
    } else {
        Ok(curve.clone())
    }
}

/// Closed `ConcatC1` last-segment degree-doubling reparameterisation
/// (`GeomConvert.cxx:1153-1230`).
fn closed_need2deg_repara(
    curve2: &GeomBSplineCurve,
    curve1: &GeomBSplineCurve,
    array0: &GeomBSplineCurve,
) -> Result<GeomBSplineCurve, String> {
    let (_, vec1) = curve2.d1(curve2.last_parameter());
    let (_, vec2) = curve1.d1(curve1.first_parameter());
    let m1 = vec1.magnitude();
    let m2 = vec2.magnitude();
    if m1 == 0.0 {
        return Err("GeomConvert Need2DegRepara: zero tangent".into());
    }
    let lambda = m2 / m1;

    let (_, vec2_end) = curve1.d1(curve1.last_parameter());
    let (_, vec1_start) = array0.d1(array0.first_parameter());
    let m2e = vec2_end.magnitude();
    let m1s = vec1_start.magnitude();
    if m2e == 0.0 {
        return Err("GeomConvert Need2DegRepara: zero tangent".into());
    }
    let lambda2 = m1s / m2e;

    let umin = curve1.first_parameter();
    let umax = curve1.last_parameter();
    let tmax = 2.0 * lambda * (umax - umin) / (1.0 + lambda * lambda2);
    let a = (lambda * lambda2 - 1.0) / (2.0 * lambda * tmax);
    let b = 1.0 / lambda;
    let c = umin;
    let poly = [c, b, a];

    let (mut knot_c1, mut knot_c1_mults) = curve1.distinct_knots_and_mults();
    if knot_c1.is_empty() {
        return Err("GeomConvert Need2DegRepara: empty knots".into());
    }
    let curve1_flat = knot_sequence(&knot_c1, &knot_c1_mults, curve1.degree() as i32);
    // `cxx:1172-1177`: remap the distinct knots through the inverse of
    // `a t^2 + b t + c`.
    knot_c1[0] = 0.0;
    for k in knot_c1.iter_mut().skip(1) {
        if a == 0.0 {
            return Err("GeomConvert Need2DegRepara: singular quadratic".into());
        }
        let disc = b * b - 4.0 * a * (c - *k);
        if disc < 0.0 {
            return Err("GeomConvert Need2DegRepara: negative discriminant".into());
        }
        *k = (-b + disc.sqrt()) / (2.0 * a);
    }
    // `cxx:1179-1183`: `KnotC1Mults(ii) = Degree + KnotC1Mults(ii)`.
    for m in &mut knot_c1_mults {
        *m = curve1.degree() as i32 + *m;
    }
    let new_degree = 2 * curve1.degree();
    let flat_knots = knot_sequence(&knot_c1, &knot_c1_mults, new_degree as i32);
    let num_new = flat_knots.len().saturating_sub(new_degree + 1);
    if num_new == 0 {
        return Err("GeomConvert Need2DegRepara: empty new poles".into());
    }

    // `cxx:1193-1199`: homogeneous poles.
    let mut poles_num = Vec::with_capacity(curve1.nb_poles() * 3);
    for i in 0..curve1.nb_poles() {
        let w = curve1.weight((i + 1) as i32);
        let p = curve1.pole(i);
        poles_num.push(p.x() * w);
        poles_num.push(p.y() * w);
        poles_num.push(p.z() * w);
    }
    let weights: Vec<f64> = (1..=curve1.nb_poles())
        .map(|i| curve1.weight(i as i32))
        .collect();

    let eval = |der: i32, _se: [f64; 2], p: f64| reparameterise_evaluator(&poly, der, p);

    let mut new_poles = vec![0.0; num_new * 3];
    let mut new_weights = vec![0.0; num_new];
    function_reparameterise(
        &eval,
        curve1.degree() as i32,
        &curve1_flat,
        &poles_num,
        3,
        &flat_knots,
        new_degree as i32,
        &mut new_poles,
    )
    .map_err(|_| "GeomConvert Need2DegRepara: FunctionReparameterise poles".to_string())?;
    function_reparameterise(
        &eval,
        curve1.degree() as i32,
        &curve1_flat,
        &weights,
        1,
        &flat_knots,
        new_degree as i32,
        &mut new_weights,
    )
    .map_err(|_| "GeomConvert Need2DegRepara: FunctionReparameterise weights".to_string())?;

    let mut poles = Vec::with_capacity(num_new);
    for i in 0..num_new {
        let d = new_weights[i];
        if d == 0.0 {
            return Err("GeomConvert Need2DegRepara: zero weight".into());
        }
        poles.push(GpPnt::new(
            new_poles[i * 3] / d,
            new_poles[i * 3 + 1] / d,
            new_poles[i * 3 + 2] / d,
        ));
    }
    GeomBSplineCurve::rational(poles, new_weights, flat_knots, new_degree)
        .map_err(|e| format!("GeomConvert Need2DegRepara: result curve ({e})"))
}

/// `GeomConvert::ConcatC1` (`GeomConvert.cxx:1004-1019`), the overload without
/// an angular tolerance.
pub fn concat_c1(
    curves: &mut [GeomBSplineCurve],
    toler: &[f64],
    closed_flag: &mut bool,
    closed_tolerance: f64,
) -> Result<(Vec<i32>, Vec<GeomBSplineCurve>), String> {
    concat_c1_angular(
        curves,
        toler,
        closed_flag,
        closed_tolerance,
        ANGULAR,
    )
}

/// `GeomConvert::ConcatC1` with `AngularTolerance` (`GeomConvert.cxx:1022-1310`).
///
/// `curves` is `ArrayOfCurves` (reordered in place by
/// `ReorderArrayOfG1Curves`, exactly as OCCT mutates the caller's array).
/// Returns `ArrayOfIndices` and `ArrayOfConcatenated`.
pub fn concat_c1_angular(
    curves: &mut [GeomBSplineCurve],
    toler: &[f64],
    closed_flag: &mut bool,
    closed_tolerance: f64,
    angular_tolerance: f64,
) -> Result<(Vec<i32>, Vec<GeomBSplineCurve>), String> {
    let nb_curve = curves.len();
    if nb_curve == 0 {
        return Ok((vec![0], Vec::new()));
    }
    let mut local_tolerance = toler.to_vec();
    let mut tab_g1 = vec![false; nb_curve.saturating_sub(1)];

    // `cxx:1051-1081`: junction continuity.
    let mut pre_last = 0.0;
    for i in 0..nb_curve {
        if i >= 1 {
            let first = curves[i].first_parameter();
            let cont = continuity(
                &curves[i - 1],
                &curves[i],
                pre_last,
                first,
                true,
                true,
                local_tolerance[i - 1],
                angular_tolerance,
            )?;
            if cont < Shape::C0 {
                return Err("GeomConvert curves not C0".into());
            }
            tab_g1[i - 1] = cont >= Shape::G1;
        }
        pre_last = curves[i].last_parameter();
    }

    // `cxx:1083-1095`: group count.
    let mut index = 0usize;
    let mut nb_group = 0usize;
    let mut nb_vertex_group0 = 0usize;
    while index <= nb_curve - 1 {
        let mut nb_vertex_g1 = 0usize;
        while (index + nb_vertex_g1) <= nb_curve.saturating_sub(2) && tab_g1[index + nb_vertex_g1]
        {
            nb_vertex_g1 += 1;
        }
        nb_group += 1;
        if index == 0 {
            nb_vertex_group0 = nb_vertex_g1;
        }
        index = index + 1 + nb_vertex_g1;
    }

    if *closed_flag && nb_group != 1 {
        nb_group -= 1;
        reorder_array_of_g1_curves(
            curves,
            &mut local_tolerance,
            &mut tab_g1,
            nb_vertex_group0,
            closed_tolerance,
        );
    }

    let mut array_of_indices = vec![0i32; nb_group + 1];
    let mut array_of_concatenated: Vec<GeomBSplineCurve> = Vec::with_capacity(nb_group);

    pretreatment(curves);
    let need_double = need_2deg_repara(curves);
    // `cxx:1117-1123`: only the `MaxDegree()/2` guard clears the closed flag.
    if nb_group == 1 && *closed_flag && need_double {
        let curve1 = &curves[nb_curve - 1];
        if curve1.degree() > MAX_DEGREE / 2 {
            *closed_flag = false;
        }
    }

    if nb_group == 1 && *closed_flag {
        // `cxx:1126-1247`: the closed single-group case.
        array_of_indices[0] = 0;
        array_of_indices[1] = 0;
        let indexmin = indexmin(curves);
        if indexmin != nb_curve - 1 {
            reorder_array_of_g1_curves(
                curves,
                &mut local_tolerance,
                &mut tab_g1,
                indexmin,
                closed_tolerance,
            );
        }
        let mut curve2: Option<GeomBSplineCurve> = None;
        for j in 0..nb_curve {
            let mut curve1 = prepare_curve_for_concat(&curves[j])?;
            if j == 0 {
                curve2 = Some(curve1);
            } else {
                if j == nb_curve - 1 && need_double {
                    let prev = curve2
                        .as_ref()
                        .ok_or_else(|| "GeomConvert Concatenation Error".to_string())?;
                    curve1 = closed_need2deg_repara(prev, &curve1, &curves[0])?;
                }
                let mut c = CompCurveToBSplineCurve::from_bspline(
                    curve2.take().unwrap(),
                    ParameterisationType::TgtThetaOver2,
                );
                // `cxx:1232-1233`: no `After` argument, i.e. `false`.
                let fusion = c
                    .add(&curve1, local_tolerance[j - 1], false, true, 0)
                    .map_err(|e| format!("GeomConvert Concatenation Error ({e:?})"))?;
                if !fusion {
                    return Err("GeomConvert Concatenation Error".into());
                }
                curve2 = c.into_curve();
            }
        }
        let mut curve2 = curve2.ok_or_else(|| "GeomConvert Concatenation Error".to_string())?;
        // `cxx:1240-1243`.
        curve2.set_periodic();
        let last_idx = curve2.last_u_knot_index();
        let mult = curve2.multiplicity(last_idx);
        curve2
            .remove_knot(last_idx, mult - 1, CONFUSION)
            .map_err(|e| format!("GeomConvert Concatenation Error ({e})"))?;
        array_of_concatenated.push(curve2);
    } else {
        // `cxx:1252-1306`: the open / multi-group case.
        let mut index = 0usize;
        let mut k = 0usize;
        for _ in 0..nb_group {
            let mut nb_vertex_g1 = 0usize;
            while (index + nb_vertex_g1) <= nb_curve.saturating_sub(2)
                && tab_g1[index + nb_vertex_g1]
            {
                nb_vertex_g1 += 1;
            }

            if !*closed_flag || nb_group == 1 {
                k += 1;
                array_of_indices[k - 1] = index as i32;
                if k == nb_group {
                    array_of_indices[k] = 0;
                }
            } else {
                k += 1;
                array_of_indices[k - 1] = (index + nb_vertex_group0 + 1) as i32;
                if k == nb_group {
                    array_of_indices[k] = (nb_vertex_group0 + 1) as i32;
                }
            }

            let mut group_curve: Option<GeomBSplineCurve> = None;
            for j in index..=(index + nb_vertex_g1) {
                let curve1 = prepare_curve_for_concat(&curves[j])?;
                if index == j {
                    group_curve = Some(curve1);
                } else {
                    let mut c = CompCurveToBSplineCurve::from_bspline(
                        group_curve
                            .take()
                            .ok_or_else(|| "GeomConvert Concatenation Error".to_string())?,
                        ParameterisationType::TgtThetaOver2,
                    );
                    // `cxx:1297`: `Add(Curve1, local_tolerance(j-1), true)`.
                    let fusion = c
                        .add(&curve1, local_tolerance[j - 1], true, true, 0)
                        .map_err(|e| format!("GeomConvert Concatenation Error ({e:?})"))?;
                    if !fusion {
                        return Err("GeomConvert Concatenation Error".into());
                    }
                    group_curve = c.into_curve();
                }
            }
            array_of_concatenated
                .push(group_curve.ok_or_else(|| "GeomConvert Concatenation Error".to_string())?);
            index = index + 1 + nb_vertex_g1;
        }
    }

    Ok((array_of_indices, array_of_concatenated))
}

/// `GeomConvert::ConcatG1` (`GeomConvert.cxx:764-997`): concatenates the curves
/// as far as G1 continuity holds; a junction that is only C0 starts a new curve.
///
/// `curves` is `ArrayOfCurves` and is mutated in place by `Pretreatment` and
/// `ReorderArrayOfG1Curves`, exactly as OCCT mutates the caller's array;
/// `toler` is `ArrayOfToler` (its length is `curves.len() - 1`).
/// `closed_flag` is `ClosedFlag` in/out — the `MaxDegree()/2` guard of
/// `cxx:849-857` clears it — and the answer is `ArrayOfConcatenated`.
///
/// The junction test is the 6-argument `GeomLProp::Continuity`
/// (`GeomLProp.cxx:207-215`), i.e. `Precision::Confusion()` /
/// `Precision::Angular()`. Unlike `ConcatC1` there is no
/// `NeedToBeTreated` / `MultNumandDenom` pass: the curves reach
/// `GeomConvert_CompCurveToBSplineCurve::Add` unchanged (`cxx:949-950`,
/// `:986-988`), and the open-group `Add` uses `ArrayOfToler`, not the possibly
/// reordered `local_tolerance` (`cxx:988`).
pub fn concat_g1(
    curves: &mut [GeomBSplineCurve],
    toler: &[f64],
    closed_flag: &mut bool,
    closed_tolerance: f64,
) -> Result<Vec<GeomBSplineCurve>, String> {
    let nb_curve = curves.len();
    if nb_curve == 0 {
        return Ok(Vec::new());
    }
    let mut local_tolerance = toler.to_vec();
    let mut tab_g1 = vec![false; nb_curve.saturating_sub(1)];

    // `cxx:791-815`: junction continuity.
    let mut pre_last = 0.0;
    for i in 0..nb_curve {
        if i >= 1 {
            let first = curves[i].first_parameter();
            let cont = continuity_default(&curves[i - 1], &curves[i], pre_last, first, true, true)?;
            if cont < Shape::C0 {
                return Err("GeomConvert curves not C0".into());
            }
            tab_g1[i - 1] = cont >= Shape::G1;
        }
        pre_last = curves[i].last_parameter();
    }

    // `cxx:817-829`: group count.
    let mut index = 0usize;
    let mut nb_group = 0usize;
    let mut nb_vertex_group0 = 0usize;
    while index <= nb_curve - 1 {
        let mut nb_vertex_g1 = 0usize;
        while index + nb_vertex_g1 < nb_curve - 1 && tab_g1[index + nb_vertex_g1] {
            nb_vertex_g1 += 1;
        }
        nb_group += 1;
        if index == 0 {
            nb_vertex_group0 = nb_vertex_g1;
        }
        index = index + 1 + nb_vertex_g1;
    }

    if *closed_flag && nb_group != 1 {
        nb_group -= 1;
        reorder_array_of_g1_curves(
            curves,
            &mut local_tolerance,
            &mut tab_g1,
            nb_vertex_group0,
            closed_tolerance,
        );
    }

    let mut array_of_concatenated: Vec<GeomBSplineCurve> = Vec::with_capacity(nb_group);

    pretreatment(curves);
    let need_double = need_2deg_repara(curves);
    // `cxx:849-857`: only the MaxDegree()/2 guard clears ClosedG1Flag.
    if nb_group == 1 && *closed_flag && need_double {
        let curve1 = &curves[nb_curve - 1];
        if curve1.degree() > MAX_DEGREE / 2 {
            *closed_flag = false;
        }
    }

    if nb_group == 1 && *closed_flag {
        // `cxx:859-964`: the closed single-group case.
        let indexmin = indexmin(curves);
        if indexmin != nb_curve - 1 {
            reorder_array_of_g1_curves(
                curves,
                &mut local_tolerance,
                &mut tab_g1,
                indexmin,
                closed_tolerance,
            );
        }
        let mut curve2: Option<GeomBSplineCurve> = None;
        for j in 0..nb_curve {
            let mut curve1 = curves[j].clone();
            if j == 0 {
                curve2 = Some(curve1);
            } else {
                if j == nb_curve - 1 && need_double {
                    let prev = curve2
                        .as_ref()
                        .ok_or_else(|| "GeomConvert Concatenation Error".to_string())?;
                    curve1 = closed_need2deg_repara(prev, &curve1, &curves[0])?;
                }
                let mut c = CompCurveToBSplineCurve::from_bspline(
                    curve2
                        .take()
                        .ok_or_else(|| "GeomConvert Concatenation Error".to_string())?,
                    ParameterisationType::TgtThetaOver2,
                );
                // `cxx:949-950`: `Add(Curve1, local_tolerance(j-1))`, `After` default.
                let fusion = c
                    .add(&curve1, local_tolerance[j - 1], false, true, 0)
                    .map_err(|e| format!("GeomConvert Concatenation Error ({e:?})"))?;
                if !fusion {
                    return Err("GeomConvert Concatenation Error".into());
                }
                curve2 = c.into_curve();
            }
        }
        let mut curve2 = curve2.ok_or_else(|| "GeomConvert Concatenation Error".to_string())?;
        // `cxx:956-961`.
        curve2.set_periodic();
        let last_idx = curve2.last_u_knot_index();
        let mult = curve2.multiplicity(last_idx);
        curve2
            .remove_knot(last_idx, mult - 1, CONFUSION)
            .map_err(|e| format!("GeomConvert Concatenation Error ({e})"))?;
        array_of_concatenated.push(curve2);
    } else {
        // `cxx:966-996`: one pass per G1 group.
        let mut index = 0usize;
        for _ in 0..nb_group {
            let mut nb_vertex_g1 = 0usize;
            while index + nb_vertex_g1 < nb_curve - 1 && tab_g1[index + nb_vertex_g1] {
                nb_vertex_g1 += 1;
            }
            let mut group_curve: Option<GeomBSplineCurve> = None;
            for j in index..=(index + nb_vertex_g1) {
                let curve1 = curves[j].clone();
                if index == j {
                    group_curve = Some(curve1);
                } else {
                    let mut c = CompCurveToBSplineCurve::from_bspline(
                        group_curve
                            .take()
                            .ok_or_else(|| "GeomConvert Concatenation Error".to_string())?,
                        ParameterisationType::TgtThetaOver2,
                    );
                    // `cxx:988`: `Add(Curve1, ArrayOfToler(j-1))`.
                    let fusion = c
                        .add(&curve1, toler[j - 1], false, true, 0)
                        .map_err(|e| format!("GeomConvert Concatenation Error ({e:?})"))?;
                    if !fusion {
                        return Err("GeomConvert Concatenation Error".into());
                    }
                    group_curve = c.into_curve();
                }
            }
            array_of_concatenated
                .push(group_curve.ok_or_else(|| "GeomConvert Concatenation Error".to_string())?);
            index = index + 1 + nb_vertex_g1;
        }
    }

    Ok(array_of_concatenated)
}

/// `GeomConvert::C0BSplineToArrayOfC1BSplineCurve(BS, tabBS, tolerance)`
/// (`GeomConvert.cxx:1344-1354`): `Precision::Angular()` default.
pub fn c0_bspline_to_array_of_c1_bspline_curve(
    bs: &GeomBSplineCurve,
    tolerance: f64,
) -> Result<Vec<GeomBSplineCurve>, String> {
    c0_bspline_to_array_of_c1_bspline_curve_angular(bs, ANGULAR, tolerance)
}

/// `GeomConvert::C0BSplineToArrayOfC1BSplineCurve(BS, tabBS,
/// AngularTolerance, tolerance)` (`GeomConvert.cxx:1356-1425`): split `BS` at
/// the knots of multiplicity `Degree` and concatenate every maximal
/// C1-continuous run.
pub fn c0_bspline_to_array_of_c1_bspline_curve_angular(
    bs: &GeomBSplineCurve,
    angular_tolerance: f64,
    tolerance: f64,
) -> Result<Vec<GeomBSplineCurve>, String> {
    let (bs_knots, bs_mults) = bs.distinct_knots_and_mults();
    let first_u = bs.first_u_knot_index();
    let last_u = bs.last_u_knot_index();
    let degree = bs.degree() as i32;

    // `cxx:1366-1374`: count of C1 spans, from `FirstUKnotIndex()`.
    let mut nbcurve_c1 = 1i32;
    for i in first_u..=(last_u - 1) {
        if bs_mults[(i - 1) as usize] == degree {
            nbcurve_c1 += 1;
        }
    }
    if nbcurve_c1 <= 1 {
        return Ok(vec![bs.clone()]);
    }

    let mut array_of_curves: Vec<GeomBSplineCurve> = Vec::with_capacity(nbcurve_c1 as usize);
    let array_of_toler = vec![tolerance; (nbcurve_c1 - 1) as usize];

    // `cxx:1386-1403`: one segment per C1 span.
    let mut u2 = bs.first_parameter();
    let mut j = first_u + 1;
    for _ in 0..nbcurve_c1 {
        let u1 = u2;
        while j < last_u && bs_mults[(j - 1) as usize] < degree {
            j += 1;
        }
        u2 = bs_knots[(j - 1) as usize];
        j += 1;
        let mut bsbis = bs.clone();
        bsbis
            .segment(u1, u2, PCONFUSION)
            .map_err(|e| format!("GeomConvert Segment: {e}"))?;
        array_of_curves.push(bsbis);
    }

    // `cxx:1408-1413`: the closed flag of the 3D version.
    let (_, v1) = bs.d1(bs.first_parameter());
    let (_, v2) = bs.d1(bs.last_parameter());
    let mut closed_flag = bs.is_closed() && v1.is_parallel_ang(&v2, angular_tolerance);

    let (_indices, concatenated) = concat_c1_angular(
        &mut array_of_curves,
        &array_of_toler,
        &mut closed_flag,
        tolerance,
        angular_tolerance,
    )?;
    Ok(concatenated)
}

/// `GeomConvert::C0BSplineToC1BSplineCurve(BS, tolerance, AngularTol)`
/// (`GeomConvert.cxx:1313-1338`).
///
/// The re-merge of the C1 pieces is unconditional in OCCT (`cxx:1336`), so the
/// call answers the concatenation even when the pieces collapsed to one curve.
pub fn c0_bspline_to_c1_bspline_curve(
    bs: &GeomBSplineCurve,
    tolerance: f64,
    angular_tolerance: f64,
) -> Result<GeomBSplineCurve, String> {
    let array_of_concatenated =
        c0_bspline_to_array_of_c1_bspline_curve_angular(bs, angular_tolerance, tolerance)?;
    let first = array_of_concatenated
        .first()
        .ok_or_else(|| "GeomConvert Concatenation Error".to_string())?;
    let mut c = CompCurveToBSplineCurve::from_bspline(
        first.clone(),
        ParameterisationType::TgtThetaOver2,
    );
    for curve in array_of_concatenated.iter().skip(1) {
        // `cxx:1330`: `Add(ArrayOfConcatenated->Value(i), tolerance)`.
        let fusion = c
            .add(curve, tolerance, false, true, 0)
            .map_err(|e| format!("GeomConvert Concatenation Error ({e:?})"))?;
        if !fusion {
            return Err("GeomConvert Concatenation Error".into());
        }
    }
    c.into_curve()
        .ok_or_else(|| "GeomConvert Concatenation Error".to_string())
}
