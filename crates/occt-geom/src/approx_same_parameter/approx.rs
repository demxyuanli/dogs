//! `IncreaseNbPoles`, the evaluator and `RebuildPCurve`.

use super::*;

/// `Approx_SameParameter::IncreaseNbPoles` (`cxx:857-991`).
pub(crate) fn increase_nb_poles(
    poles: &[f64],
    flat: &[f64],
    data: &mut DistData,
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    best_sq_tol: &mut f64,
) -> bool {
    let mut new_pc2d = [0.0; MAX_ARRAY];
    let mut new_pc3d = [0.0; MAX_ARRAY];
    let mut newcount = 0usize;
    for ii in 0..data.nb {
        new_pc2d[newcount] = data.pc2d[ii];
        new_pc3d[newcount] = data.pc3d[ii];
        newcount += 1;
        if data.nb - ii + newcount == MAX_ARRAY {
            continue;
        }
        let mid_c3d = 0.5 * (data.pc3d[ii] + data.pc3d[ii + 1]);
        let eval_result = eval_1d_d0(flat, poles, mid_c3d);
        if eval_result < data.pc2d[ii] || eval_result > data.pc2d[ii + 1] {
            let ucons = 0.5 * (data.pc2d[ii] + data.pc2d[ii + 1]);
            let uc3d = 0.5 * (data.pc3d[ii] + data.pc3d[ii + 1]);
            let pcons = cons_value(c2d, s, ucons);
            if let Some(curp) =
                project_point_on_curve(uc3d, &pcons, data.tol, c3d, data.c3d_f, data.c3d_l)
            {
                let dist_2 = c3d.d0(curp).square_distance(&pcons);
                if dist_2 > *best_sq_tol {
                    *best_sq_tol = dist_2;
                }
                if curp > data.pc3d[ii] + PCONFUSION && curp < data.pc3d[ii + 1] - PCONFUSION {
                    new_pc3d[newcount] = curp;
                    new_pc2d[newcount] = ucons;
                    newcount += 1;
                }
            }
        }
    }
    new_pc3d[newcount] = data.pc3d[data.nb];
    new_pc2d[newcount] = data.pc2d[data.nb];
    if data.nb != newcount && newcount < MAX_ARRAY - 1 {
        data.pc2d = new_pc2d;
        data.pc3d = new_pc3d;
        data.nb = newcount;
        return true;
    }

    newcount = 0;
    for n in 0..data.nb {
        new_pc3d[newcount] = data.pc3d[n];
        new_pc2d[newcount] = data.pc2d[n];
        newcount += 1;
        if data.nb - n + newcount == MAX_ARRAY {
            continue;
        }
        let ucons = 0.5 * (data.pc2d[n] + data.pc2d[n + 1]);
        let uc3d = 0.5 * (data.pc3d[n] + data.pc3d[n + 1]);
        let pcons = cons_value(c2d, s, ucons);
        if let Some(curp) =
            project_point_on_curve(uc3d, &pcons, data.tol, c3d, data.c3d_f, data.c3d_l)
        {
            let dist_2 = c3d.d0(curp).square_distance(&pcons);
            if dist_2 > *best_sq_tol {
                *best_sq_tol = dist_2;
            }
            if curp > data.pc3d[n] + PCONFUSION && curp < data.pc3d[n + 1] - PCONFUSION {
                new_pc3d[newcount] = curp;
                new_pc2d[newcount] = ucons;
                newcount += 1;
            }
        }
    }
    new_pc3d[newcount] = data.pc3d[data.nb];
    new_pc2d[newcount] = data.pc2d[data.nb];
    if data.nb != newcount {
        data.pc2d = new_pc2d;
        data.pc3d = new_pc3d;
        data.nb = newcount;
        return true;
    }
    false
}

pub(crate) fn approx_evaluator_eval<'a>(
    flat: &'a [f64],
    poles: &'a [f64],
    c2d: &'a dyn Curve2d,
) -> impl Fn(f64, f64, f64, i32, &mut [f64]) -> i32 + 'a {
    // `StartEnd` unused: same-parameter evaluator has no CurveOnSurface re-trim.
    move |_a: f64, _b: f64, param: f64, der: i32, out: &mut [f64]| -> i32 {
        if der == 0 {
            let tcons = eval_1d_d0(flat, poles, param);
            let uv = c2d.d0(tcons);
            if out.len() >= 2 {
                out[0] = uv.x();
                out[1] = uv.y();
            }
            0
        } else if der == 1 {
            let (tcons, dt) = eval_1d_d1(flat, poles, param);
            let (_, duv) = c2d.d1(tcons);
            if out.len() >= 2 {
                out[0] = duv.x() * dt;
                out[1] = duv.y() * dt;
            }
            0
        } else {
            1
        }
    }
}

pub(crate) fn rebuild_pcurve(
    data: &mut DistData,
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    tang_first: f64,
    tang_last: f64,
) -> (bool, f64, Option<Arc<dyn Curve2d>>) {
    let cont = if c2d.continuity() > 2 {
        Shape::C1
    } else if c2d.continuity() == 0 {
        Shape::C0
    } else {
        Shape::C1
    };
    let mut besttol2 = data.tol * data.tol;
    let mut tolsov = INFINITE;
    let mut done = false;
    let mut curve2d: Option<Arc<dyn Curve2d>> = None;
    let mut tol_reached = compute_tol_reached(c3d, c2d, s, data.c3d_f, data.c3d_l);
    let mut interpolok;
    let mut has_count_changed = false;
    // `cxx:400-480`: densify loop over IncreaseNbPoles until AdvApprox accepts.
    loop {
        let Some((poles, flat)) = interpolate_reparam(data, tang_first, tang_last) else {
            return (false, tol_reached, None);
        };
        let mut algtol = besttol2.sqrt();
        interpolok = check_reparam(
            &flat,
            &poles,
            data.nb + 1,
            &data.pc3d,
            c3d,
            c2d,
            s,
            &mut algtol,
            tolsov,
        );
        tolsov = algtol;
        if interpolok {
            let besttol = besttol2.sqrt();
            let tol_u = u_resolution(s, besttol);
            let tol_v = v_resolution(s, besttol);
            let eval = approx_evaluator_eval(&flat, &poles, c2d);
            if let Ok(approx) = ApproxAFunction1dPair::approx(
                data.c3d_f,
                data.c3d_l,
                cont,
                11,
                1000,
                tol_u,
                tol_v,
                &eval,
            ) {
                if approx.done || approx.has_result {
                    if let Some(new_c2d) = curve2d_from_approx(&approx) {
                        let new_tol =
                            compute_tol_reached(c3d, new_c2d.as_ref(), s, data.c3d_f, data.c3d_l);
                        const MULT: f64 = 250.0;
                        if new_tol < MULT * besttol {
                            done = true;
                            tol_reached = new_tol;
                            curve2d = Some(new_c2d);
                            break;
                        } else if data.nb < MAX_ARRAY - 1 {
                            interpolok = false;
                        } else {
                            break;
                        }
                    }
                }
            }
        }
        if !interpolok {
            has_count_changed =
                increase_nb_poles(&poles, &flat, data, c3d, c2d, s, &mut besttol2);
        }
        if interpolok || !has_count_changed {
            break;
        }
    }

    // `cxx:482-539`: post-loop fallback — accept AdvApprox if better than original.
    if !done {
        tol_reached = compute_tol_reached(c3d, c2d, s, data.c3d_f, data.c3d_l);
        let Some((poles, flat)) = interpolate_reparam(data, tang_first, tang_last) else {
            return (false, tol_reached, None);
        };
        let besttol = besttol2.sqrt();
        let tol_u = u_resolution(s, besttol);
        let tol_v = v_resolution(s, besttol);
        let eval = approx_evaluator_eval(&flat, &poles, c2d);
        let Ok(approx) = ApproxAFunction1dPair::approx(
            data.c3d_f,
            data.c3d_l,
            cont,
            11,
            40,
            tol_u,
            tol_v,
            &eval,
        ) else {
            return (false, tol_reached, None);
        };
        if !approx.done && !approx.has_result {
            return (false, tol_reached, None);
        }
        let Some(new_c2d) = curve2d_from_approx(&approx) else {
            return (false, tol_reached, None);
        };
        let approx_tol = compute_tol_reached(c3d, new_c2d.as_ref(), s, data.c3d_f, data.c3d_l);
        if approx_tol < tol_reached {
            tol_reached = approx_tol;
            curve2d = Some(new_c2d);
        }
        done = true;
    }
    (done, tol_reached, curve2d)
}

