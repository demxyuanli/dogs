//! `GeomConvert::ConcatG1` (`cxx:1029-1101`).

use super::*;

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

