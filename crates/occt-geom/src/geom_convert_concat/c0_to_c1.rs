//! `C0BSplineToArrayOfC1BSplineCurve` and `C0BSplineToC1BSplineCurve` (`cxx:1303-1361`).

use super::*;

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
