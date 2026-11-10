//! `CheckSameParameter`, `ComputeTangents`, `Interpolate` and `Check`.

use super::*;

pub(crate) fn check_same_parameter(
    data: &mut DistData,
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    sq_dist: &mut f64,
) -> bool {
    let tol2 = data.tol * data.tol;
    let mut is_same = true;
    let pcons = cons_value(c2d, s, data.c2d_f);
    let pc3d = c3d.d0(data.c3d_f);
    let mut dmax2 = pcons.square_distance(&pc3d);
    let pcons = cons_value(c2d, s, data.c2d_l);
    let pc3d = c3d.d0(data.c3d_l);
    dmax2 = dmax2.max(pcons.square_distance(&pc3d));

    let mut count = 1usize;
    let mut previousp = data.c3d_f;
    let mut initp = 0.0;
    let bornesup = data.c3d_l - PCONFUSION;
    let mut is_proj_ok = false;
    for ii in 1..data.nb {
        let pcons = cons_value(c2d, s, data.pc2d[ii]);
        let pc3d_pt = c3d.d0(data.pc3d[ii]);
        let dist2 = pcons.square_distance(&pc3d_pt);
        let is_use = dist2 <= tol2 && data.pc3d[ii] > data.pc3d[count - 1] + PCONFUSION;
        if is_use {
            if dmax2 < dist2 {
                dmax2 = dist2;
            }
            initp = data.pc3d[ii];
            previousp = initp;
            data.pc3d[count] = data.pc3d[ii];
            data.pc2d[count] = data.pc2d[ii];
            count += 1;
            continue;
        }
        if !is_proj_ok {
            initp = data.pc3d[ii];
        }
        is_proj_ok = false;
        is_same = false;
        let mut curp = initp;
        if let Some(p) = project_point_on_curve(initp, &pcons, data.tol, c3d, data.c3d_f, data.c3d_l)
        {
            curp = p;
            is_proj_ok = true;
        }
        is_proj_ok = is_proj_ok && curp > previousp + PCONFUSION && curp < bornesup;
        if is_proj_ok {
            initp = curp;
            previousp = curp;
            data.pc3d[count] = curp;
            data.pc2d[count] = data.pc2d[ii];
            count += 1;
            continue;
        }
        // Whole parameter space search (`Extrema_ExtPC`, cxx:728-758).
        let extrema = crate::extrema_pc::point_curve_extrema_all(c3d, &pcons);
        if extrema.is_empty() {
            continue;
        }
        let mut best_u: Option<f64> = None;
        let mut best_d2 = f64::MAX;
        for e in &extrema {
            if e.u1 < data.c3d_f - PCONFUSION || e.u1 > data.c3d_l + PCONFUSION {
                continue;
            }
            let d2 = e.p2.square_distance(&pcons);
            if d2 < best_d2 {
                best_d2 = d2;
                best_u = Some(e.u1);
            }
        }
        if let Some(curp) = best_u {
            if curp > previousp + PCONFUSION && curp < bornesup {
                initp = curp;
                previousp = curp;
                data.pc3d[count] = curp;
                data.pc2d[count] = data.pc2d[ii];
                count += 1;
                is_proj_ok = true;
            }
        }
    }
    data.nb = count;
    data.pc2d[data.nb] = data.c2d_l;
    data.pc3d[data.nb] = data.c3d_l;
    *sq_dist = dmax2;
    is_same
}

pub(crate) fn compute_tangents(c3d: &dyn Curve, c2d: &dyn Curve2d, s: &dyn Surface) -> Option<(f64, f64)> {
    const SMALL: f64 = 1.0e-12;
    let param_first = c3d.first_parameter();
    let (_, vec_cons) = cons_d1(c2d, s, param_first);
    let (_, vec) = c3d.d1(param_first);
    let mag = vec_cons.magnitude();
    if mag <= SMALL {
        return None;
    }
    let first_t = vec.magnitude() / mag;
    let param_last = c3d.last_parameter();
    let (_, vec_cons) = cons_d1(c2d, s, param_last);
    let (_, vec) = c3d.d1(param_last);
    let mag = vec_cons.magnitude();
    if mag <= SMALL {
        return None;
    }
    let last_t = vec.magnitude() / mag;
    Some((first_t, last_t))
}

pub(crate) fn interpolate_reparam(data: &DistData, tang_first: f64, tang_last: f64) -> Option<(Vec<f64>, Vec<f64>)> {
    let num_poles = data.nb + 3;
    let num_knots = data.nb + 7;
    let mut poles = vec![0.0; num_poles];
    let mut flat = vec![0.0; num_knots];
    let mut contact = vec![0i32; num_poles];
    let mut parameters = vec![0.0; num_poles];
    contact[1] = 1;
    contact[num_poles - 2] = 1;
    for i in 0..4 {
        flat[i] = data.c3d_f;
        flat[num_poles + i] = data.c3d_l;
    }
    poles[0] = data.c2d_f;
    poles[num_poles - 1] = data.c2d_l;
    poles[1] = tang_first;
    poles[num_poles - 2] = tang_last;
    parameters[0] = data.c3d_f;
    parameters[1] = data.c3d_f;
    parameters[num_poles - 2] = data.c3d_l;
    parameters[num_poles - 1] = data.c3d_l;
    for ii in 3..=(num_poles - 2) {
        poles[ii - 1] = data.pc2d[ii - 2];
        parameters[ii - 1] = data.pc3d[ii - 2];
        flat[ii + 1] = data.pc3d[ii - 2];
    }
    interpolate_contact(3, &flat, &parameters, &contact, &mut poles, 1).ok()?;
    Some((poles, flat))
}

pub(crate) fn eval_1d_d0(flat: &[f64], poles: &[f64], u: f64) -> f64 {
    let pts: Vec<GpPnt> = poles.iter().map(|&x| GpPnt::new(x, 0.0, 0.0)).collect();
    eval_curve(&pts, flat, 3, u).x()
}

pub(crate) fn eval_1d_d1(flat: &[f64], poles: &[f64], u: f64) -> (f64, f64) {
    let pts: Vec<GpPnt> = poles.iter().map(|&x| GpPnt::new(x, 0.0, 0.0)).collect();
    let (p, d) = eval_curve_d1(&pts, flat, 3, u);
    (p.x(), d.x())
}

