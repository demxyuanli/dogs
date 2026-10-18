//! `Geom2dConvert::ConcatC1` and `C0BSplineToC1BSplineCurve`
//! (`Geom2dConvert.cxx:1149-1542`).
//!
//! Ported branches:
//! * junction `Continuity` (`cxx:748-852`);
//! * `Pretreatment` / `NeedToBeTreated` / `Need2DegRepara` / `Indexmin` /
//!   `ReorderArrayOfG1` (`cxx:571-712`);
//! * open / multi-group `ConcatC1` via `CompCurveToBSplineCurve2d::Add`
//!   (`cxx:1394-1453`);
//! * closed single-group path including Need2DegRepara
//!   (`cxx:1272-1391`);
//! * `C0BSplineToC1BSplineCurve` segment + concat + re-merge
//!   (`cxx:1459-1541`).
//!
use crate::bspline_curve::Geom2dBSplineCurve;
use crate::comp_curve_to_bspline::CompCurveToBSplineCurve2d;
use crate::curve::Curve2d;
use crate::geom2d_convert::DEFAULT_PARAMETERISATION;
use crate::geom2d_convert_mult::mult_num_and_denom;
use crate::hermit;
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::function_ops::function_reparameterise;
use occt_core::bspl::plib::eval_polynomial_horner;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};

/// `GeomAbsToInteger` (`Geom2dConvert.cxx:716-743`).
fn geom_abs_to_integer(gcont: u8) -> i32 {
    i32::from(gcont)
}

/// Local `Continuity` (`Geom2dConvert.cxx:748-852`).
/// Returns `Err` where OCCT throws `Standard_Failure("Courbes non jointives")`.
fn continuity_at(
    c1: &dyn Curve2d,
    c2: &dyn Curve2d,
    u1: f64,
    u2: f64,
    r1: bool,
    r2: bool,
    tl: f64,
    ta: f64,
) -> Result<Shape, &'static str> {
    let mut cont = Shape::C0;
    let mut cont1 = geom_abs_to_integer(c1.continuity());
    let mut cont2 = geom_abs_to_integer(c2.continuity());

    let basis1: &dyn Curve2d = c1.trimmed_basis().unwrap_or(c1);
    let basis2: &dyn Curve2d = c2.trimmed_basis().unwrap_or(c2);

    if let Some(bs) = basis1.bspline_copy2d() {
        let tolerance = bs.resolution(tl);
        let (index1, index2) = bs.locate_u(u1, tolerance);
        let nb = bs.nb_knots() as i32;
        if index1 > 1 && index2 < nb && index1 == index2 {
            let (_, umults) = bs.distinct_knots_and_mults();
            cont1 = bs.degree() as i32 - umults[(index1 - 1) as usize];
        } else {
            cont1 = 5;
        }
    }
    if let Some(bs) = basis2.bspline_copy2d() {
        let tolerance = bs.resolution(tl);
        let (index1, index2) = bs.locate_u(u2, tolerance);
        let nb = bs.nb_knots() as i32;
        if index1 > 1 && index2 < nb && index1 == index2 {
            let (_, umults) = bs.distinct_knots_and_mults();
            cont2 = bs.degree() as i32 - umults[(index1 - 1) as usize];
        } else {
            cont2 = 5;
        }
    }

    let (point1, mut d1) = c1.d1(u1);
    let (point2, mut d2) = c2.d1(u2);
    if point1.square_distance(&point2) > tl * tl {
        return Err("Geom2dConvert Continuity: curves not jointive");
    }
    if cont1 != 0 && cont2 != 0 {
        let tl2 = tl * tl;
        if d1.square_magnitude() >= tl2 && d2.square_magnitude() >= tl2 {
            if r1 {
                d1.set_coord(-d1.x(), -d1.y());
            }
            if r2 {
                d2.set_coord(-d2.x(), -d2.y());
            }
            let value = d1.dot(&d2);
            let m1 = d1.magnitude();
            let m2 = d2.magnitude();
            if m1 <= (m2 + tl) && m1 >= (m2 - tl) && value / (m1 * m2) >= 1.0 - ta * ta {
                cont = Shape::C1;
            } else if let (Ok(n1), Ok(n2)) = (d1.normalized(), d2.normalized()) {
                let value = n1.dot(&n2).abs();
                if value >= 1.0 - ta * ta {
                    cont = Shape::G1;
                }
            }
        }
    }
    Ok(cont)
}

/// `Pretreatment` (`Geom2dConvert.cxx:571-593`).
fn pretreatment(tab: &mut [Geom2dBSplineCurve]) {
    for curve in tab.iter_mut() {
        if !curve.is_rational() {
            continue;
        }
        let a = curve.weight(1);
        let n = curve.nb_poles();
        if n >= 2
            && curve.weight(2) == a
            && curve.weight(n - 1) == a
            && curve.weight(n) == a
        {
            for j in 1..=n {
                curve.set_weight(j, curve.weight(j) / a);
            }
        }
    }
}

/// `NeedToBeTreated` (`Geom2dConvert.cxx:597-617`).
fn need_to_be_treated(bs: &Geom2dBSplineCurve) -> bool {
    if !bs.is_rational() {
        return false;
    }
    let Some(weights) = bs.bspline_weights2d() else {
        return false;
    };
    if !occt_core::bspl::rational::is_rational(weights) {
        return false;
    }
    let n = bs.nb_poles();
    if n < 2 {
        return false;
    }
    let conf = CONFUSION;
    let out = |w: f64| w < (1.0 - conf) || w > (1.0 + conf);
    out(bs.weight(1))
        || out(bs.weight(2))
        || out(bs.weight(n - 1))
        || out(bs.weight(n))
}

/// `Need2DegRepara` (`Geom2dConvert.cxx:621-637`).
fn need_2deg_repara(tab: &[Geom2dBSplineCurve]) -> bool {
    let mut rapport = 1.0;
    for i in 0..tab.len().saturating_sub(1) {
        let (_, v1) = tab[i + 1].d1(tab[i + 1].first_parameter());
        let (_, v2) = tab[i].d1(tab[i].last_parameter());
        let m1 = v1.magnitude();
        let m2 = v2.magnitude();
        if m1 > 0.0 {
            rapport *= m2 / m1;
        }
    }
    rapport > (1.0 + CONFUSION) || rapport < (1.0 - CONFUSION)
}

/// `Indexmin` (`Geom2dConvert.cxx:641-655`).
fn indexmin(tab: &[Geom2dBSplineCurve]) -> usize {
    let mut index = 0;
    let mut degree = tab[0].degree();
    for (i, c) in tab.iter().enumerate() {
        if c.degree() <= degree {
            degree = c.degree();
            index = i;
        }
    }
    index
}

/// `ReorderArrayOfG1` (`Geom2dConvert.cxx:659-712`).
fn reorder_array_of_g1(
    curves: &mut [Geom2dBSplineCurve],
    toler: &mut [f64],
    tab_g1: &mut [bool],
    start_index: usize,
    closed_tolerance: f64,
) {
    let n = curves.len();
    if n == 0 {
        return;
    }
    let bis_curves: Vec<Geom2dBSplineCurve> = curves.to_vec();
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

/// `NeedToBeTreated` arm (`Geom2dConvert.cxx:1283-1285`, `:1428-1430`):
/// `MultNumandDenom(Hermit::Solution(curve), curve)`.
fn prepare_curve_for_concat(
    curve: &Geom2dBSplineCurve,
) -> Result<Geom2dBSplineCurve, &'static str> {
    if need_to_be_treated(curve) {
        let law = hermit::solution_default(curve)?;
        mult_num_and_denom(&law, curve)
    } else {
        Ok(curve.clone())
    }
}

/// Closed ConcatC1 last-segment degree-doubling reparameterisation
/// (`Geom2dConvert.cxx:1300-1376`).
fn closed_need2deg_repara(
    curve2: &Geom2dBSplineCurve,
    curve1: &Geom2dBSplineCurve,
    array0: &Geom2dBSplineCurve,
) -> Result<Geom2dBSplineCurve, &'static str> {
    let (_, vec1) = curve2.d1(curve2.last_parameter());
    let (_, vec2) = curve1.d1(curve1.first_parameter());
    let m1 = vec1.magnitude();
    let m2 = vec2.magnitude();
    if m1 == 0.0 {
        return Err("Geom2dConvert Need2DegRepara: zero tangent");
    }
    let lambda = m2 / m1;
    let (_, vec2_end) = curve1.d1(curve1.last_parameter());
    let (_, vec1_start) = array0.d1(array0.first_parameter());
    let m2e = vec2_end.magnitude();
    let m1s = vec1_start.magnitude();
    if m2e == 0.0 {
        return Err("Geom2dConvert Need2DegRepara: zero tangent");
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
    let curve1_flat = knot_sequence(&knot_c1, &knot_c1_mults, curve1.degree() as i32);
    if knot_c1.is_empty() {
        return Err("Geom2dConvert Need2DegRepara: empty knots");
    }
    // Remap distinct knots via the inverse of a t^2 + b t + c (`cxx:1320-1325`).
    knot_c1[0] = 0.0;
    for k in knot_c1.iter_mut().skip(1) {
        if a == 0.0 {
            return Err("Geom2dConvert Need2DegRepara: singular quadratic");
        }
        let disc = b * b - 4.0 * a * (c - *k);
        if disc < 0.0 {
            return Err("Geom2dConvert Need2DegRepara: negative discriminant");
        }
        *k = (-b + disc.sqrt()) / (2.0 * a);
    }
    for m in &mut knot_c1_mults {
        *m = curve1.degree() as i32 + *m;
    }
    let new_degree = 2 * curve1.degree();
    let flat_knots = knot_sequence(&knot_c1, &knot_c1_mults, new_degree as i32);
    let num_new = flat_knots.len().saturating_sub(new_degree + 1);
    if num_new == 0 {
        return Err("Geom2dConvert Need2DegRepara: empty new poles");
    }

    let mut poles_num = Vec::with_capacity(curve1.nb_poles() * 2);
    for i in 0..curve1.nb_poles() {
        let w = curve1.weight(i + 1);
        poles_num.push(curve1.xs[i] * w);
        poles_num.push(curve1.ys[i] * w);
    }
    let weights: Vec<f64> = (1..=curve1.nb_poles())
        .map(|i| curve1.weight(i))
        .collect();

    let eval = |_der: i32, _se: [f64; 2], p: f64| -> Result<f64, ()> {
        Ok(eval_polynomial_horner(&poly, p))
    };

    let mut new_poles = vec![0.0; num_new * 2];
    let mut new_weights = vec![0.0; num_new];
    function_reparameterise(
        &eval,
        curve1.degree() as i32,
        &curve1_flat,
        &poles_num,
        2,
        &flat_knots,
        new_degree as i32,
        &mut new_poles,
    )
    .map_err(|_| "Geom2dConvert Need2DegRepara: FunctionReparameterise poles")?;
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
    .map_err(|_| "Geom2dConvert Need2DegRepara: FunctionReparameterise weights")?;

    let mut xs = Vec::with_capacity(num_new);
    let mut ys = Vec::with_capacity(num_new);
    for i in 0..num_new {
        let d = new_weights[i];
        if d == 0.0 {
            return Err("Geom2dConvert Need2DegRepara: zero weight");
        }
        xs.push(new_poles[i * 2] / d);
        ys.push(new_poles[i * 2 + 1] / d);
    }
    Geom2dBSplineCurve::rational(xs, ys, new_weights, flat_knots, new_degree)
        .map_err(|_| "Geom2dConvert Need2DegRepara: result curve")
}

fn shape_ge(a: Shape, b: Shape) -> bool {
    a >= b
}

fn shape_lt(a: Shape, b: Shape) -> bool {
    a < b
}

/// `Geom2dConvert::ConcatC1` with angular tolerance
/// (`Geom2dConvert.cxx:1168-1455`).
pub fn concat_c1(
    array_of_curves: &mut [Geom2dBSplineCurve],
    array_of_toler: &[f64],
    closed_flag: &mut bool,
    closed_tolerance: f64,
    angular_tolerance: f64,
) -> Result<(Vec<i32>, Vec<Geom2dBSplineCurve>), &'static str> {
    let nb_curve = array_of_curves.len();
    if nb_curve == 0 {
        return Ok((vec![0], Vec::new()));
    }
    let mut local_tolerance = array_of_toler.to_vec();
    let mut tab_g1 = vec![false; nb_curve.saturating_sub(1)];
    let mut pre_last = 0.0;
    for i in 0..nb_curve {
        if i >= 1 {
            let first = array_of_curves[i].first_parameter();
            let cont = continuity_at(
                &array_of_curves[i - 1],
                &array_of_curves[i],
                pre_last,
                first,
                true,
                true,
                array_of_toler[i - 1],
                angular_tolerance,
            )?;
            if shape_lt(cont, Shape::C0) {
                return Err("Geom2dConvert curves not C0");
            }
            tab_g1[i - 1] = shape_ge(cont, Shape::G1);
        }
        pre_last = array_of_curves[i].last_parameter();
    }

    let mut index = 0usize;
    let mut nb_group = 0usize;
    let mut nb_vertex_group0 = 0usize;
    while index <= nb_curve - 1 {
        let mut nb_vertex_g1 = 0usize;
        while (index + nb_vertex_g1) <= nb_curve.saturating_sub(2)
            && tab_g1[index + nb_vertex_g1]
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
        reorder_array_of_g1(
            array_of_curves,
            &mut local_tolerance,
            &mut tab_g1,
            nb_vertex_group0,
            closed_tolerance,
        );
    }

    let mut array_of_indices = vec![0i32; nb_group + 1];
    let mut array_of_concatenated: Vec<Geom2dBSplineCurve> = Vec::with_capacity(nb_group);

    pretreatment(array_of_curves);
    let need_double = need_2deg_repara(array_of_curves);
    // `cxx:1263-1269`: only the MaxDegree/2 guard clears ClosedFlag.
    if nb_group == 1 && *closed_flag && need_double {
        let curve1 = &array_of_curves[nb_curve - 1];
        if curve1.degree() > 25 / 2 {
            *closed_flag = false;
        }
    }

    if nb_group == 1 && *closed_flag {
        array_of_indices[0] = 0;
        array_of_indices[1] = 0;
        let indexmin = indexmin(array_of_curves);
        if indexmin != nb_curve - 1 {
            reorder_array_of_g1(
                array_of_curves,
                &mut local_tolerance,
                &mut tab_g1,
                indexmin,
                closed_tolerance,
            );
        }
        let mut curve2: Option<Geom2dBSplineCurve> = None;
        for j in 0..nb_curve {
            let mut curve1 = prepare_curve_for_concat(&array_of_curves[j])?;
            if j == 0 {
                curve2 = Some(curve1);
            } else {
                // `cxx:1300-1376`: last segment of a closed Need2DegRepara group.
                if j == nb_curve - 1 && need_double {
                    let prev = curve2.as_ref().ok_or("Geom2dConvert Concatenation Error")?;
                    curve1 = closed_need2deg_repara(prev, &curve1, &array_of_curves[0])?;
                }
                let mut c = CompCurveToBSplineCurve2d::from_bspline(
                    curve2.take().unwrap(),
                    DEFAULT_PARAMETERISATION,
                );
                let fusion = c
                    .add(&curve1, local_tolerance[j - 1], true)
                    .map_err(|_| "Geom2dConvert Concatenation Error")?;
                if !fusion {
                    return Err("Geom2dConvert Concatenation Error");
                }
                curve2 = c.into_curve();
            }
        }
        let mut curve2 = curve2.ok_or("Geom2dConvert Concatenation Error")?;
        curve2.set_periodic();
        let last_idx = curve2.last_u_knot_index();
        let (_, umults) = curve2.distinct_knots_and_mults();
        let m = umults[(last_idx - 1) as usize];
        let _ = curve2.remove_knot(last_idx, m - 1, CONFUSION);
        array_of_concatenated.push(curve2);
    } else {
        let mut index = 0usize;
        let mut k = 0usize;
        for i in 0..nb_group {
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

            let mut group_curve: Option<Geom2dBSplineCurve> = None;
            for j in index..=(index + nb_vertex_g1) {
                let curve1 = prepare_curve_for_concat(&array_of_curves[j])?;
                if index == j {
                    group_curve = Some(curve1);
                } else {
                    let mut c = CompCurveToBSplineCurve2d::from_bspline(
                        group_curve.take().unwrap(),
                        DEFAULT_PARAMETERISATION,
                    );
                    let fusion = c
                        .add(&curve1, array_of_toler[j - 1], true)
                        .map_err(|_| "Geom2dConvert Concatenation Error")?;
                    if !fusion {
                        return Err("Geom2dConvert Concatenation Error");
                    }
                    group_curve = c.into_curve();
                }
            }
            array_of_concatenated.push(group_curve.ok_or("Geom2dConvert Concatenation Error")?);
            let _ = i;
            index = index + 1 + nb_vertex_g1;
        }
    }

    Ok((array_of_indices, array_of_concatenated))
}

/// `Geom2dConvert::C0BSplineToC1BSplineCurve` (`Geom2dConvert.cxx:1459-1542`).
pub fn c0_bspline_to_c1_bspline_curve(
    bs: &mut Geom2dBSplineCurve,
    tolerance: f64,
) -> Result<(), &'static str> {
    let (uknots, umults) = bs.distinct_knots_and_mults();
    let first_u = bs.first_u_knot_index();
    let last_u = bs.last_u_knot_index();
    let mut nbcurve_c1 = 1i32;
    for i in (first_u + 1)..=(last_u - 1) {
        if umults[(i - 1) as usize] == bs.degree() as i32 {
            nbcurve_c1 += 1;
        }
    }
    nbcurve_c1 = nbcurve_c1.min((bs.nb_knots() as i32) - 1);
    if nbcurve_c1 <= 1 {
        return Ok(());
    }

    let mut array_of_curves: Vec<Geom2dBSplineCurve> =
        Vec::with_capacity(nbcurve_c1 as usize);
    let array_of_toler = vec![tolerance; (nbcurve_c1 - 1) as usize];

    let mut u2 = bs.first_parameter();
    let mut j = first_u + 1;
    for _ in 0..nbcurve_c1 {
        let u1 = u2;
        while j < last_u && umults[(j - 1) as usize] < bs.degree() as i32 {
            j += 1;
        }
        u2 = uknots[(j - 1) as usize];
        j += 1;
        let mut bsbis = bs.clone();
        bsbis.segment(u1, u2, PCONFUSION)?;
        array_of_curves.push(bsbis);
    }

    const AN_ANGULAR_TOLER: f64 = 1.0e-7;
    let (point1, v1) = bs.d1(bs.first_parameter());
    let (point2, v2) = bs.d1(bs.last_parameter());
    let mut closed_flag = point1.square_distance(&point2) < tolerance * tolerance
        && v1.is_parallel(&v2, AN_ANGULAR_TOLER);

    let (_indices, concatenated) = concat_c1(
        &mut array_of_curves,
        &array_of_toler,
        &mut closed_flag,
        tolerance,
        ANGULAR,
    )?;

    if concatenated.is_empty() {
        return Ok(());
    }
    let mut c = CompCurveToBSplineCurve2d::from_bspline(
        concatenated[0].clone(),
        DEFAULT_PARAMETERISATION,
    );
    if concatenated.len() >= 2 {
        for curve in concatenated.iter().skip(1) {
            let fusion = c
                .add(curve, tolerance, true)
                .map_err(|_| "Geom2dConvert Concatenation Error")?;
            if !fusion {
                return Err("Geom2dConvert Concatenation Error");
            }
        }
    }
    if let Some(out) = c.into_curve() {
        *bs = out;
    }
    Ok(())
}
