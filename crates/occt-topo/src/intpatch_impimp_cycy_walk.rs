//! CyCyNoGeometric adaptive walker.
//! Source: `IntPatch_ImpImpIntersection.cxx` CyCyNoGeometric (~6573-7877).

use std::f64::consts::PI;

use occt_core::gp::{GpCylinder, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL, SQUARE_CONFUSION};

use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::{PatchPoint, WLine, WLineWay};

use super::bounds::{critical_points, cyl_cyl_monotonicity, URange};
use super::step::{fill_step_matrix, step_computing};
use super::wl::{add_boundary_point, add_point_into_wl, seek_additional_points};
use super::{
    compute_params, compute_u2, compute_v, inscribe_point, Coeffs, SurfDom, PERIOD,
};
use super::super::glines::PairOutcome;
use super::super::quad::{normale, ImplicitQuad};

#[derive(Clone, Copy, PartialEq, Eq)]
enum WlStatus {
    Absent,
    Exist,
    Broken,
}

fn u1_of(line: &WLine, i1: i32, reversed: bool) -> (f64, f64) {
    if reversed {
        line.point(i1).parameters_on_s2()
    } else {
        line.point(i1).parameters_on_s1()
    }
}

fn u2_of(line: &WLine, i1: i32, reversed: bool) -> (f64, f64) {
    if reversed {
        line.point(i1).parameters_on_s1()
    } else {
        line.point(i1).parameters_on_s2()
    }
}

fn adjust_u2_period(line: &WLine, reversed: bool, u2: &mut f64, period: f64) {
    let n = line.nb_pnts();
    if n < 1 {
        return;
    }
    let (u2p, _) = u2_of(line, n, reversed);
    let delta = *u2 - u2p;
    if 2.0 * delta.abs() > period {
        if delta > 0.0 {
            *u2 -= period;
        } else {
            *u2 += period;
        }
    }
}

fn new_wline() -> WLine {
    let mut w = WLine::new();
    w.creating_way = WLineWay::ImpImp;
    w
}

fn emit_wline(
    cyl1: &GpCylinder,
    cyl2: &GpCylinder,
    coeffs: &Coeffs,
    reversed: bool,
    w: &mut WLine,
    wl: i32,
    nb_points: i32,
    nb_max: i32,
    tol3d: f64,
    tol2d: f64,
    period: f64,
    lines: &mut Vec<GeomIntLine>,
    points: &mut Vec<PatchPoint>,
    added: &mut bool,
) {
    let n = w.nb_pnts();
    if n == 1 && !*added {
        let p = *w.point(1);
        let same = points.last().is_some_and(|q| q.p.square_distance(&p.p) <= CONFUSION * CONFUSION);
        if !same {
            points.push(PatchPoint::new(p.p, p.u1, p.u1, p.v1, p.u2, p.v2));
        }
        return;
    }
    if n <= 1 {
        *added = false;
        return;
    }
    let mut good = true;
    if n == 2 {
        if w.point(1).is_same(w.point(2), CONFUSION, -1.0) {
            good = false;
        }
    } else if n > 2 {
        let q1 = ImplicitQuad::Cylinder(cyl1.clone());
        let q2 = ImplicitQuad::Cylinder(cyl2.clone());
        let sq_tol = tol3d * tol3d;
        for j in 0..2 {
            loop {
                if w.nb_pnts() >= nb_max {
                    break;
                }
                let idx1 = if j == 1 { w.nb_pnts() - 1 } else { 2 };
                let idx2 = if j == 1 { w.nb_pnts() } else { 1 };
                let p1 = w.point(idx1).value();
                let p2 = w.point(idx2).value();
                let dir = GpVec::from_pnts(&p1, &p2);
                if dir.square_magnitude() < sq_tol {
                    break;
                }
                let n1 = normale(&q1, &p2);
                let n2 = normale(&q2, &p2);
                let tg = n1.crossed(&n2);
                if tg.square_magnitude() < SQUARE_CONFUSION {
                    break;
                }
                let mut ang = dir.angle(&tg);
                if ang > PI / 2.0 {
                    ang -= PI;
                }
                if ang.abs() > 0.25 {
                    let prev = w.nb_pnts();
                    seek_additional_points(
                        cyl1, cyl2, w, coeffs, wl, 3, idx1, idx2, tol2d, period, reversed,
                    );
                    if w.nb_pnts() == prev {
                        break;
                    }
                } else {
                    break;
                }
            }
        }
    }
    if good {
        *added = true;
        let n_end = w.nb_pnts();
        seek_additional_points(
            cyl1,
            cyl2,
            w,
            coeffs,
            wl,
            nb_points,
            1,
            n_end,
            tol2d,
            period,
            reversed,
        );
        w.ensure_end_vertices();
        lines.push(GeomIntLine::Walking(std::mem::replace(w, new_wline())));
    } else {
        *added = false;
    }
}

pub(crate) fn cy_cy_walk(
    cyl1: &GpCylinder,
    cyl2: &GpCylinder,
    coeffs: Coeffs,
    reversed: bool,
    dom: SurfDom,
    ranges: [URange; 2],
    dv1: f64,
    dv2: f64,
    is_good: bool,
    step_min: f64,
    step_max: f64,
    nb_points: i32,
    nb_min: i32,
    nb_max: i32,
) -> Result<PairOutcome, bool> {
    let mut lines = Vec::new();
    let mut points = Vec::new();
    let mut crit = critical_points(
        &coeffs,
        dom.u1f,
        dom.u1l,
        dom.u2f,
        dom.u2l,
        PERIOD,
        dom.tol2d,
    );
    let nb_wl = 2;
    for range in ranges {
        let Some((mut uf, ul)) = range.get_bounds() else {
            continue;
        };
        let is_delta_period = (ul - uf) == PERIOD;
        for u in &mut crit {
            let _ = inscribe_point(uf, ul, u, 0.0, PERIOD, false);
        }
        crit.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut added_into = [false; 2];
        while uf < ul {
            let mut u2 = [0.0; 2];
            let mut v1 = [0.0; 2];
            let mut v2 = [0.0; 2];
            let mut v1_prev = [0.0; 2];
            let mut v2_prev = [0.0; 2];
            let mut u_expect = [uf; 2];
            let mut status = [WlStatus::Absent; 2];
            let mut enabled = [true; 2];
            let mut wline = [new_wline(), new_wline()];
            let mut crit_delta = vec![0.0; crit.len()];
            for (i, c) in crit.iter().enumerate() {
                crit_delta[i] = uf - *c;
            }
            let mut u1 = uf;
            let u1_min = uf;
            let mut u1_prev = uf;
            let mut is_first = true;
            while u1 <= ul {
                for (i, c) in crit.iter().enumerate() {
                    if (u1 - *c) * crit_delta[i] < 0.0 {
                        u1 = *c;
                        for j in 0..nb_wl {
                            status[j] = WlStatus::Broken;
                            u_expect[j] = u1;
                        }
                        break;
                    }
                }
                if u1 == ul {
                    for i in 0..nb_wl {
                        status[i] = WlStatus::Broken;
                        u_expect[i] = u1;
                        enabled[i] = if is_delta_period {
                            !added_into[i]
                        } else {
                            (dom.tol2d >= (u_expect[i] - u1)) || status[i] == WlStatus::Absent
                        };
                    }
                } else {
                    for i in 0..nb_wl {
                        enabled[i] =
                            (dom.tol2d >= (u_expect[i] - u1)) || status[i] == WlStatus::Absent;
                    }
                }
                for i in 0..nb_wl {
                    let n = wline[i].nb_pnts();
                    let wl = i as i32;
                    if status[i] == WlStatus::Broken || status[i] == WlStatus::Absent {
                        let mut tol = dom.tol2d;
                        if let Some(uu) = compute_u2(u1, wl, &coeffs, Some(&mut tol)) {
                            u2[i] = uu;
                        }
                        let _ = inscribe_point(
                            dom.u2f,
                            dom.u2l,
                            &mut u2[i],
                            dom.tol2d,
                            PERIOD,
                            false,
                        );
                        tol = tol.max(dom.tol2d);
                        if u2[i].abs() <= tol {
                            u2[i] = 0.0;
                        } else if (u2[i] - PERIOD).abs() <= tol {
                            u2[i] = PERIOD;
                        } else if (u2[i] - dom.u2f).abs() <= tol {
                            u2[i] = dom.u2f;
                        } else if (u2[i] - dom.u2l).abs() <= tol {
                            u2[i] = dom.u2l;
                        }
                    } else if let Some(uu) = compute_u2(u1, wl, &coeffs, None) {
                        u2[i] = uu;
                        let _ = inscribe_point(
                            dom.u2f,
                            dom.u2l,
                            &mut u2[i],
                            dom.tol2d,
                            PERIOD,
                            false,
                        );
                    }
                    if n == 0 {
                        if (dom.u2f + PERIOD - dom.u2l) <= 2.0 * dom.tol2d
                            && ((u2[i] - dom.u2f).abs() < dom.tol2d
                                || (u2[i] - dom.u2l).abs() < dom.tol2d)
                        {
                            if let Some(inc) =
                                cyl_cyl_monotonicity(u1 + step_min, wl, &coeffs, PERIOD)
                            {
                                u2[i] = if inc { dom.u2f } else { dom.u2l };
                            }
                        }
                    } else if (dom.u2l - dom.u2f) >= PERIOD
                        && ((u2[i] - dom.u2f).abs() < dom.tol2d
                            || (u2[i] - dom.u2l).abs() < dom.tol2d)
                    {
                        let (u2prev, _) = u2_of(&wline[i], n, reversed);
                        if 2.0 * (u2prev - u2[i]).abs() > PERIOD {
                            if u2prev > u2[i] {
                                u2[i] += PERIOD;
                            } else {
                                u2[i] -= PERIOD;
                            }
                        }
                    }
                    let (vv1, vv2) = compute_v(u1, u2[i], &coeffs);
                    v1[i] = vv1;
                    v2[i] = vv2;
                    if is_first {
                        v1_prev[i] = v1[i];
                        v2_prev[i] = v2[i];
                    }
                }
                is_first = false;
                let mut is_broken = false;
                for i in 0..nb_wl {
                    let wl = i as i32;
                    if !enabled[i] {
                        let mut bound_x = false;
                        if (v1[i] - dom.v1f).abs() <= dom.tol2d
                            || (v1[i] - dom.v1f) * (v1_prev[i] - dom.v1f) < 0.0
                        {
                            bound_x = true;
                        } else if (v1[i] - dom.v1l).abs() <= dom.tol2d
                            || (v1[i] - dom.v1l) * (v1_prev[i] - dom.v1l) < 0.0
                        {
                            bound_x = true;
                        } else if (v2[i] - dom.v2f).abs() <= dom.tol2d
                            || (v2[i] - dom.v2f) * (v2_prev[i] - dom.v2f) < 0.0
                        {
                            bound_x = true;
                        } else if (v2[i] - dom.v2l).abs() <= dom.tol2d
                            || (v2[i] - dom.v2l) * (v2_prev[i] - dom.v2l) < 0.0
                        {
                            bound_x = true;
                        }
                        if status[i] == WlStatus::Broken {
                            is_broken = true;
                        }
                        if !bound_x {
                            continue;
                        }
                        u_expect[i] = u1;
                    }
                    let is_inscribe = (dom.u2f - u2[i]) <= dom.tol2d
                        && (u2[i] - dom.u2l) <= dom.tol2d
                        && (dom.v1f - v1[i]) <= dom.tol2d
                        && (v1[i] - dom.v1l) <= dom.tol2d
                        && (dom.v2f - v2[i]) <= dom.tol2d
                        && (v2[i] - dom.v2l) <= dom.tol2d;
                    // `IntPatch_ImpImpIntersection.cxx:7055-7058`:
                    // `(((aVSurf1f - aV1[i]) * (aVSurf1f - aV1Prev[i]) < RealSmall()) && ...)`.
                    let is_v_int = (((dom.v1f - v1[i]) * (dom.v1f - v1_prev[i]) < REAL_SMALL)
                        && ((dom.v1l - v1[i]) * (dom.v1l - v1_prev[i]) < REAL_SMALL))
                        || (((dom.v2f - v2[i]) * (dom.v2f - v2_prev[i]) < REAL_SMALL)
                            && ((dom.v2l - v2[i]) * (dom.v2l - v2_prev[i]) < REAL_SMALL));
                    let mut force = false;
                    if status[i] == WlStatus::Absent
                        && (dom.u2l - dom.u2f) >= PERIOD
                        && (u1 - dom.u1l).abs() < dom.tol2d
                    {
                        force = true;
                    }
                    let (found1, found2) = add_boundary_point(
                        cyl1,
                        cyl2,
                        &coeffs,
                        reversed,
                        &mut wline[i],
                        u1,
                        u1_prev,
                        u1_min,
                        u2[i],
                        v1[i],
                        v1_prev[i],
                        v2[i],
                        v2_prev[i],
                        wl,
                        force,
                        dom,
                    );
                    let prev_v_bound = !is_v_int
                        && ((v1_prev[i] - dom.v1f).abs() <= dom.tol2d
                            || (v1_prev[i] - dom.v1l).abs() <= dom.tol2d
                            || (v2_prev[i] - dom.v2f).abs() <= dom.tol2d
                            || (v2_prev[i] - dom.v2l).abs() <= dom.tol2d);
                    v1_prev[i] = v1[i];
                    v2_prev[i] = v2[i];
                    if status[i] == WlStatus::Exist && (found1 || found2) && !prev_v_bound {
                        status[i] = WlStatus::Broken;
                    } else if is_inscribe {
                        if status[i] == WlStatus::Absent && (found1 || found2) {
                            status[i] = WlStatus::Exist;
                        }
                        if status[i] != WlStatus::Broken
                            || wline[i].nb_pnts() >= 1
                            || u1 == ul
                        {
                            adjust_u2_period(&wline[i], reversed, &mut u2[i], PERIOD);
                            if add_point_into_wl(
                                cyl1,
                                cyl2,
                                &coeffs,
                                reversed,
                                true,
                                u1,
                                v1[i],
                                u2[i],
                                v2[i],
                                dom,
                                &mut wline[i],
                                wl,
                                force,
                                false,
                            ) {
                                if status[i] == WlStatus::Absent {
                                    status[i] = WlStatus::Exist;
                                }
                            } else if !found1 && !found2 && status[i] == WlStatus::Exist {
                                status[i] = WlStatus::Broken;
                            }
                        }
                    } else if status[i] == WlStatus::Exist {
                        status[i] = WlStatus::Broken;
                    }
                    if status[i] == WlStatus::Broken {
                        is_broken = true;
                    }
                }
                if is_broken {
                    uf = u1;
                    let mut is_added = true;
                    for i in 0..nb_wl {
                        if enabled[i] {
                            continue;
                        }
                        is_added = false;
                        let (f1, f2) = add_boundary_point(
                            cyl1,
                            cyl2,
                            &coeffs,
                            reversed,
                            &mut wline[i],
                            u1,
                            u1_prev,
                            u1_min,
                            u2[i],
                            v1[i],
                            v1_prev[i],
                            v2[i],
                            v2_prev[i],
                            i as i32,
                            false,
                            dom,
                        );
                        if f1 || f2 {
                            is_added = true;
                        }
                        adjust_u2_period(&wline[i], reversed, &mut u2[i], PERIOD);
                        if add_point_into_wl(
                            cyl1,
                            cyl2,
                            &coeffs,
                            reversed,
                            true,
                            u1,
                            v1[i],
                            u2[i],
                            v2[i],
                            dom,
                            &mut wline[i],
                            i as i32,
                            false,
                            false,
                        ) {
                            is_added = true;
                        }
                    }
                    if !is_added {
                        let mut umax = -f64::MAX;
                        let mut changed = false;
                        for i in 0..nb_wl {
                            if status[i] == WlStatus::Absent || wline[i].nb_pnts() == 0 {
                                continue;
                            }
                            let n = wline[i].nb_pnts();
                            let (u1c, _) = u1_of(&wline[i], n, reversed);
                            umax = umax.max(u1c);
                            changed = true;
                        }
                        if !changed {
                            break;
                        }
                        for i in 0..nb_wl {
                            if enabled[i] {
                                continue;
                            }
                            if let Some((uu2, vv1, vv2)) =
                                compute_params(umax, i as i32, &coeffs)
                            {
                                u2[i] = uu2;
                                v1[i] = vv1;
                                v2[i] = vv2;
                                let _ = add_point_into_wl(
                                    cyl1,
                                    cyl2,
                                    &coeffs,
                                    reversed,
                                    true,
                                    umax,
                                    v1[i],
                                    u2[i],
                                    v2[i],
                                    dom,
                                    &mut wline[i],
                                    i as i32,
                                    false,
                                    false,
                                );
                            }
                        }
                    }
                    break;
                }
                {
                    let delta_v1 = dv1 / nb_points as f64;
                    let delta_v2 = dv2 / nb_points as f64;
                    let mut min_uexp = f64::MAX;
                    for i in 0..nb_wl {
                        if dom.tol2d < (u_expect[i] - u1) {
                            continue;
                        }
                        if status[i] == WlStatus::Absent || is_good {
                            u_expect[i] += step_max;
                            min_uexp = min_uexp.min(u_expect[i]);
                            continue;
                        }
                        let mut step_tmp = step_max;
                        let m = fill_step_matrix(&coeffs, u1, u2[i]);
                        if let Some(s) = step_computing(&m, v1[i], v2[i], delta_v1, delta_v2) {
                            step_tmp = s.clamp(step_min, step_max);
                            u_expect[i] = u1 + step_tmp;
                        } else {
                            u_expect[i] += step_max;
                        }
                        min_uexp = min_uexp.min(u_expect[i]);
                    }
                    u1_prev = u1;
                    u1 = min_uexp;
                }
                if PCONFUSION >= (ul - u1) {
                    u1 = ul;
                }
                uf = u1;
                for i in 0..nb_wl {
                    if wline[i].nb_pnts() != 1 {
                        added_into[i] = false;
                    }
                    if u1 == ul {
                        u_expect[i] = ul;
                    }
                }
            }
            for i in 0..nb_wl {
                emit_wline(
                    cyl1,
                    cyl2,
                    &coeffs,
                    reversed,
                    &mut wline[i],
                    i as i32,
                    nb_points,
                    nb_max,
                    dom.tol3d,
                    dom.tol2d,
                    PERIOD,
                    &mut lines,
                    &mut points,
                    &mut added_into[i],
                );
            }
        }
    }
    prune_isolated(&mut points, &lines, dom.tol3d);
    densify_isolated(
        cyl1,
        cyl2,
        &coeffs,
        reversed,
        dom,
        step_min,
        step_max,
        nb_min,
        &mut lines,
        &mut points,
    );
    if lines.is_empty() && points.is_empty() {
        Ok(PairOutcome::Empty)
    } else {
        Ok(PairOutcome::Result { lines, points })
    }
}

fn prune_isolated(points: &mut Vec<PatchPoint>, lines: &[GeomIntLine], tol3d: f64) {
    let mut i = 0;
    while i < points.len() {
        let mut drop = false;
        for line in lines {
            if let GeomIntLine::Walking(w) = line {
                if w.nb_pnts() < 1 {
                    continue;
                }
                let a = w.point(1);
                let b = w.point(w.nb_pnts());
                if points[i].p.square_distance(&a.p) <= tol3d * tol3d
                    || points[i].p.square_distance(&b.p) <= tol3d * tol3d
                {
                    drop = true;
                    break;
                }
            }
        }
        if drop {
            points.remove(i);
        } else {
            i += 1;
        }
    }
}

fn densify_isolated(
    cyl1: &GpCylinder,
    cyl2: &GpCylinder,
    coeffs: &Coeffs,
    reversed: bool,
    dom: SurfDom,
    step_min: f64,
    step_max: f64,
    nb_min: i32,
    lines: &mut Vec<GeomIntLine>,
    points: &mut Vec<PatchPoint>,
) {
    let mut i = 0;
    while i < points.len() {
        let pt = points[i];
        let (u1, v1, u2, v2) = (pt.u1, pt.v1, pt.u2, pt.v2);
        let (uf0, ul0, cur_u2) = if reversed {
            (u2 - step_max, u2 + step_max, u1)
        } else {
            (u1 - step_max, u1 + step_max, u2)
        };
        let uf0 = uf0.max(dom.u1f);
        let ul0 = ul0.min(dom.u1l);
        let umid = 0.5 * (uf0 + ul0);
        let mut index = 0;
        let mut best = f64::MAX;
        for k in 0..2 {
            if let Some(u2t) = compute_u2(umid, k, coeffs, None) {
                let mut du = (u2t - cur_u2).abs() % PERIOD;
                du = du.min((du - PERIOD).abs());
                if du < best {
                    best = du;
                    index = k;
                }
            }
        }
        let mut added = [
            if reversed { u2 } else { u1 },
            if reversed { u2 } else { u1 },
        ];
        let mut w = new_wline();
        for par in 0..2 {
            let (mut lo, mut hi) = if par == 0 {
                (uf0, umid)
            } else {
                (umid, ul0)
            };
            while (hi - lo).abs() > step_min {
                let uc = 0.5 * (lo + hi);
                let mut ok = false;
                if let Some((uu2, mut vv1, mut vv2)) = compute_params(uc, index, coeffs) {
                    if (vv1 - dom.v1f).abs() <= dom.tol2d {
                        vv1 = dom.v1f;
                    }
                    if (vv1 - dom.v1l).abs() <= dom.tol2d {
                        vv1 = dom.v1l;
                    }
                    if (vv2 - dom.v2f).abs() <= dom.tol2d {
                        vv2 = dom.v2f;
                    }
                    if (vv2 - dom.v2l).abs() <= dom.tol2d {
                        vv2 = dom.v2l;
                    }
                    ok = add_point_into_wl(
                        cyl1, cyl2, coeffs, reversed, true, uc, vv1, uu2, vv2, dom, &mut w, index,
                        false, true,
                    );
                }
                if ok {
                    added[0] = added[0].min(uc);
                    added[1] = added[1].max(uc);
                    if par == 0 {
                        hi = uc;
                    } else {
                        lo = uc;
                    }
                } else if par == 0 {
                    lo = uc;
                } else {
                    hi = uc;
                }
            }
        }
        if added[1] - added[0] > step_min {
            for par in 0..2 {
                if let Some((uu2, vv1, vv2)) = compute_params(added[par], index, coeffs) {
                    let _ = add_point_into_wl(
                        cyl1,
                        cyl2,
                        coeffs,
                        reversed,
                        true,
                        added[par],
                        vv1,
                        uu2,
                        vv2,
                        dom,
                        &mut w,
                        index,
                        false,
                        false,
                    );
                }
            }
            let n_end = w.nb_pnts();
            seek_additional_points(
                cyl1,
                cyl2,
                &mut w,
                coeffs,
                index,
                nb_min,
                1,
                n_end,
                dom.tol2d,
                PERIOD,
                reversed,
            );
            w.ensure_end_vertices();
            lines.push(GeomIntLine::Walking(w));
            points.remove(i);
        } else {
            i += 1;
        }
        let _ = (v1, v2);
    }
}
