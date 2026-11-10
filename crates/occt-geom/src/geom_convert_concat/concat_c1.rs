//! `GeomConvert::ConcatC1` both overloads (`cxx:1004-1027`).

use super::*;

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

