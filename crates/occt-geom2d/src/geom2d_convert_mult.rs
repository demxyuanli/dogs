//! `Geom2dConvert` helpers `MultNumandDenom` and `Geom2dConvert_law_evaluator`.
//! Source: `Geom2dConvert.cxx:453-567`.

use crate::bspline_curve::Geom2dBSplineCurve;
use crate::curve::Curve2d;
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::function_ops::{flat_knot_sequence, function_multiply, merge_bspline_knots};
use occt_core::bspl::knots::reparametrize;
/// `Geom2dConvert_law_evaluator` (`Geom2dConvert.cxx:453-484`): returns the
/// Y coordinate of the ancore curve at `parameter` when the request is 0 and
/// `parameter` lies in `[start_end[0], start_end[1]]`.
fn law_evaluator(ancore: &Geom2dBSplineCurve, der: i32, start_end: [f64; 2], param: f64) -> Result<f64, ()> {
    if der == 0 && param >= start_end[0] && param <= start_end[1] {
        Ok(ancore.d0(param).y())
    } else {
        Err(())
    }
}

/// `MultNumandDenom` (`Geom2dConvert.cxx:488-567`).
pub fn mult_num_and_denom(
    a: &Geom2dBSplineCurve,
    bs: &Geom2dBSplineCurve,
) -> Result<Geom2dBSplineCurve, &'static str> {
    let (bs_knots, bs_mults) = bs.distinct_knots_and_mults();
    if bs_knots.is_empty() {
        return Err("MultNumandDenom: empty BS knots");
    }
    let start_value = bs_knots[0];
    let end_value = bs_knots[bs_knots.len() - 1];
    let tolerance = 10.0 * occt_core::precision::epsilon(end_value.abs());

    let (mut a_knots, mults) = a.distinct_knots_and_mults();
    reparametrize(bs.first_parameter(), bs.last_parameter(), &mut a_knots);
    let ancore = Geom2dBSplineCurve::from_flat(
        a.xs.clone(),
        a.ys.clone(),
        None,
        knot_sequence(&a_knots, &mults, a.degree() as i32),
        a.degree(),
        false,
    )
    .map_err(|_| "MultNumandDenom: ancore")?;

    let (res_nb_poles, res_knots, res_mults) = merge_bspline_knots(
        tolerance,
        start_value,
        end_value,
        a.degree() as i32,
        &a_knots,
        &mults,
        bs.degree() as i32,
        &bs_knots,
        &bs_mults,
    );
    if res_nb_poles <= 0 {
        return Err("MultNumandDenom: MergeBSplineKnots");
    }
    let degree = (bs.degree() + a.degree()) as i32;
    let res_flat = flat_knot_sequence(&res_knots, &res_mults, degree);
    let bs_flat = bs.knots.clone();

    let mut bs_num = Vec::with_capacity(bs.nb_poles() * 2);
    for i in 0..bs.nb_poles() {
        let w = bs.weight(i + 1);
        bs_num.push(bs.xs[i] * w);
        bs_num.push(bs.ys[i] * w);
    }
    let bs_weights: Vec<f64> = (1..=bs.nb_poles()).map(|i| bs.weight(i)).collect();

    let mut res_num = vec![0.0; res_nb_poles as usize * 2];
    let mut res_den = vec![0.0; res_nb_poles as usize];
    let eval = |der: i32, se: [f64; 2], p: f64| law_evaluator(&ancore, der, se, p);

    function_multiply(
        &eval,
        bs.degree() as i32,
        &bs_flat,
        &bs_num,
        2,
        &res_flat,
        degree,
        &mut res_num,
    )
    .map_err(|_| "MultNumandDenom: FunctionMultiply num")?;
    function_multiply(
        &eval,
        bs.degree() as i32,
        &bs_flat,
        &bs_weights,
        1,
        &res_flat,
        degree,
        &mut res_den,
    )
    .map_err(|_| "MultNumandDenom: FunctionMultiply den")?;

    let mut xs = Vec::with_capacity(res_nb_poles as usize);
    let mut ys = Vec::with_capacity(res_nb_poles as usize);
    for i in 0..res_nb_poles as usize {
        let d = res_den[i];
        if d == 0.0 {
            return Err("MultNumandDenom: zero denominator pole");
        }
        xs.push(res_num[i * 2] / d);
        ys.push(res_num[i * 2 + 1] / d);
    }
    Geom2dBSplineCurve::rational(xs, ys, res_den, res_flat, degree as usize)
        .map_err(|_| "MultNumandDenom: result curve")
}
