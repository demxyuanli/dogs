//! `MultNumandDenom` (`GeomConvert.cxx:676-731`).

use super::*;

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

