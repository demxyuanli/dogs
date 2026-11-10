//! Closed `ConcatC1` last-segment degree-doubling reparameterisation (`cxx:1153-1230`).

use super::*;

/// Closed `ConcatC1` last-segment degree-doubling reparameterisation
/// (`GeomConvert.cxx:1153-1230`).
pub(crate) fn closed_need2deg_repara(
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

