//! `MaxDegree`, `law_evaluator` and the `ConcatC1` pre/post helpers.

use super::*;

/// `Geom2d_BSplineCurve::MaxDegree()` (`Geom2d_BSplineCurve.hxx`), read by
/// `ConcatC1`'s degree guard (`GeomConvert.cxx:1120`).
pub(crate) const MAX_DEGREE: usize = 25;

/// `GeomConvert_law_evaluator` (`GeomConvert.cxx:463-492`): the `Y` coordinate
/// of the 2D ancore curve when the derivative order is 0 and the parameter is
/// inside `[start_end[0], start_end[1]]`.
pub(crate) fn law_evaluator(
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


/// `Pretreatment` (`GeomConvert.cxx:587-606`): rescale the weights of a
/// rational curve whose four boundary weights share the first weight.
pub(crate) fn pretreatment(tab: &mut [GeomBSplineCurve]) {
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
pub(crate) fn need_to_be_treated(bs: &GeomBSplineCurve) -> bool {
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
pub(crate) fn need_2deg_repara(tab: &[GeomBSplineCurve]) -> bool {
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
pub(crate) fn indexmin(tab: &[GeomBSplineCurve]) -> usize {
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
pub(crate) fn reorder_array_of_g1_curves(
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
pub(crate) fn reparameterise_evaluator(coefficient: &[f64; 3], der: i32, param: f64) -> Result<f64, ()> {
    if !(0..=2).contains(&der) {
        return Err(());
    }
    Ok(occt_core::bspl::plib::eval_polynomial_horner(coefficient, param))
}

/// `NeedToBeTreated` arm (`GeomConvert.cxx:1136-1142`, `:1281-1287`):
/// `MultNumandDenom(Hermit::Solution(curve), curve)`.
pub(crate) fn prepare_curve_for_concat(curve: &GeomBSplineCurve) -> Result<GeomBSplineCurve, String> {
    if need_to_be_treated(curve) {
        let law = hermit::solution_default(curve)?;
        mult_num_and_denom(&law, curve)
    } else {
        Ok(curve.clone())
    }
}

