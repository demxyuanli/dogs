use super::prelude::*;
use super::*;

/// Newton on F = 0 from a seed, clamped to the parameter box. Returns the
/// converged stationary point, or `None` if the iterate diverged.
pub(super) fn solve_point_surface(
    s: &dyn Surface,
    p: &GpPnt,
    u0: f64,
    v0: f64,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
) -> Option<(f64, f64)> {
    let cl = |x: f64, a: f64, b: f64| x.clamp(a.min(b), a.max(b));
    let (mut u, mut v) = (cl(u0, ua, ub), cl(v0, va, vb));
    let hu = ((ub - ua).abs() * 1e-6).max(1e-9);
    let hv = ((vb - va).abs() * 1e-6).max(1e-9);
    let base = ps_f(s, p, u, v).iter().map(|x| x.abs()).sum::<f64>();
    let mut converged = false;
    for _ in 0..40 {
        let f = ps_f(s, p, u, v);
        let j = ps_jac(s, p, u, v, hu, hv);
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() < 1e-300 {
            break;
        }
        let du = (-f[0] * j[1][1] + f[1] * j[0][1]) / det;
        let dv = (-j[0][0] * f[1] + j[1][0] * f[0]) / det;
        let (un, vn) = (cl(u + du, ua, ub), cl(v + dv, va, vb));
        if (un - u).abs() < 1e-10 * (1.0 + u.abs())
            && (vn - v).abs() < 1e-10 * (1.0 + v.abs())
        {
            u = un;
            v = vn;
            converged = true;
            break;
        }
        if (un - u).abs() + (vn - v).abs() > 1e6 * (1.0 + u.abs() + v.abs()) {
            return None; // diverged
        }
        u = un;
        v = vn;
    }
    // A converged step locates a stationary point of the (finite-difference)
    // system even when the residual cannot reach machine zero at a boundary
    // (one-sided derivatives); accept it. Otherwise require a genuine root.
    if converged {
        return Some((u, v));
    }
    let resid = ps_f(s, p, u, v).iter().map(|x| x.abs()).sum::<f64>();
    if resid < base * 1e-2 + 1e-6 {
        Some((u, v))
    } else {
        None
    }
}

/// All local extrema of |S(u,v)-P| via grid seeding + Newton, deduplicated and
/// sorted by distance. Replaces the sampling path of `Extrema_ExtPS` for
/// non-analytic surfaces.
pub(crate) fn point_surface_newton_all(s: &dyn Surface, p: &GpPnt) -> Vec<ExtremaPair> {
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let (nu, nv) = (24, 24);
    let mut d2 = vec![vec![0.0; nv + 1]; nu + 1];
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = s.d0(u, v).square_distance(p);
            d2[i][j] = if d.is_finite() { d } else { f64::INFINITY };
        }
    }
    let mut seeds: Vec<(f64, f64)> = Vec::new();
    // Local extrema of the sampled squared distance (both minima and maxima of
    // the distance are stationary points of F = 0).
    for i in 1..nu {
        for j in 1..nv {
            let d = d2[i][j];
            let (du_m, du_p, dv_m, dv_p) = (d2[i - 1][j], d2[i + 1][j], d2[i][j - 1], d2[i][j + 1]);
            if d <= du_m && d <= du_p && d <= dv_m && d <= dv_p {
                seeds.push((u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64));
            } else if d >= du_m && d >= du_p && d >= dv_m && d >= dv_p {
                seeds.push((u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64));
            }
        }
    }
    // Global grid min/max + corners (boundary extrema for trimmed surfaces).
    let (mut gmin, mut gmax) = ((0usize, 0usize), (0usize, 0usize));
    for i in 0..=nu {
        for j in 0..=nv {
            if d2[i][j] < d2[gmin.0][gmin.1] {
                gmin = (i, j);
            }
            if d2[i][j] > d2[gmax.0][gmax.1] {
                gmax = (i, j);
            }
        }
    }
    seeds.push((u0 + (u1 - u0) * gmin.0 as f64 / nu as f64, v0 + (v1 - v0) * gmin.1 as f64 / nv as f64));
    seeds.push((u0 + (u1 - u0) * gmax.0 as f64 / nu as f64, v0 + (v1 - v0) * gmax.1 as f64 / nv as f64));
    seeds.push((u0, v0));
    seeds.push((u0, v1));
    seeds.push((u1, v0));
    seeds.push((u1, v1));

    let mut out: Vec<ExtremaPair> = Vec::new();
    for (u, v) in seeds {
        if let Some((uu, vv)) = solve_point_surface(s, p, u, v, u0, u1, v0, v1) {
            let q = s.d0(uu, vv);
            if !q.x().is_finite() {
                continue;
            }
            if !out.iter().any(|e: &ExtremaPair| e.p2.distance(&q) < 1e-6) {
                out.push(ps_pair(p, uu, vv, q));
            }
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

/// Fallback: coarse grid min + coordinate descent (degenerate surfaces).
pub(super) fn fallback_point_surface(s: &dyn Surface, p: &GpPnt) -> ExtremaPair {
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let n = 40;
    let mut best = (0.0, 0.0, f64::INFINITY);
    for i in 0..=n {
        for j in 0..=n {
            let u = u0 + (u1 - u0) * i as f64 / n as f64;
            let v = v0 + (v1 - v0) * j as f64 / n as f64;
            let d = s.d0(u, v).square_distance(p);
            if d < best.2 {
                best = (u, v, d);
            }
        }
    }
    let (mut u, mut v) = (best.0, best.1);
    let dist2 = |u: f64, v: f64| s.d0(u, v).square_distance(p);
    let mut step = (v1 - v0).abs().max(1.0) / 4.0;
    for _ in 0..256 {
        let cur = dist2(u, v);
        let cand = [(u + step, v), (u - step, v), (u, v + step), (u, v - step)];
        let mut bi = usize::MAX;
        let mut bb = cur;
        for (k, (cu, cv)) in cand.iter().enumerate() {
            let dd = dist2(*cu, *cv);
            if dd < bb {
                bi = k;
                bb = dd;
            }
        }
        if bi == usize::MAX {
            step *= 0.5;
            if step < 1e-12 {
                break;
            }
        } else {
            match bi {
                0 => u += step,
                1 => u -= step,
                2 => v += step,
                _ => v -= step,
            }
            u = u.clamp(u0.min(u1), u0.max(u1));
            v = v.clamp(v0.min(v1), v0.max(v1));
        }
    }
    let q = s.d0(u, v);
    ps_pair(p, u, v, q)
}

// ---------------------------------------------------------------------------
// General Newton path (curve-surface, port of `Extrema_GenExtCS` +
// `Extrema_FuncExtCS`, numeric 3×3 Jacobian).
// ---------------------------------------------------------------------------

/// F(t, u, v) = ((C-S)·C′, (C-S)·Su, (C-S)·Sv).
pub(super) fn cs_f(c: &dyn Curve, s: &dyn Surface, t: f64, u: f64, v: f64) -> Option<[f64; 3]> {
    let (pc, dc) = c.d1(t);
    let (ps, su, sv) = s.d1(u, v);
    if !(pc.x().is_finite() && ps.x().is_finite()) {
        return None;
    }
    let w = GpVec::from_pnts(&ps, &pc);
    Some([w.dot(&dc), w.dot(&su), w.dot(&sv)])
}

/// Newton on F = 0 (3×3) from a seed, clamped to the parameter box.
pub(super) fn solve_curve_surface(
    c: &dyn Curve,
    s: &dyn Surface,
    t0: f64,
    u0: f64,
    v0: f64,
    ta: f64,
    tb: f64,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
) -> Option<(f64, f64, f64)> {
    let cl = |x: f64, a: f64, b: f64| x.clamp(a.min(b), a.max(b));
    let mut x = [cl(t0, ta, tb), cl(u0, ua, ub), cl(v0, va, vb)];
    let hs = [
        ((tb - ta).abs() * 1e-6).max(1e-9),
        ((ub - ua).abs() * 1e-6).max(1e-9),
        ((vb - va).abs() * 1e-6).max(1e-9),
    ];
    let base = cs_f(c, s, x[0], x[1], x[2])?.iter().map(|x| x.abs()).sum::<f64>();
    let mut converged = false;
    for _ in 0..40 {
        let f = cs_f(c, s, x[0], x[1], x[2])?;
        let mut j = [[0.0; 3]; 3];
        for k in 0..3 {
            let h = hs[k];
            let mut xp = x;
            xp[k] += h;
            let mut xm = x;
            xm[k] -= h;
            let fp = cs_f(c, s, xp[0], xp[1], xp[2])?;
            let fm = cs_f(c, s, xm[0], xm[1], xm[2])?;
            for r in 0..3 {
                j[r][k] = (fp[r] - fm[r]) / (2.0 * h);
            }
        }
        let d = match solve3(&j, &f) {
            Some(d) => d,
            None => break,
        };
        let xn = [
            cl(x[0] - d[0], ta, tb),
            cl(x[1] - d[1], ua, ub),
            cl(x[2] - d[2], va, vb),
        ];
        if (xn[0] - x[0]).abs() < 1e-10 * (1.0 + x[0].abs())
            && (xn[1] - x[1]).abs() < 1e-10 * (1.0 + x[1].abs())
            && (xn[2] - x[2]).abs() < 1e-10 * (1.0 + x[2].abs())
        {
            x = xn;
            converged = true;
            break;
        }
        if (xn[0] - x[0]).abs() + (xn[1] - x[1]).abs() + (xn[2] - x[2]).abs()
            > 1e6 * (1.0 + x[0].abs() + x[1].abs() + x[2].abs())
        {
            return None;
        }
        x = xn;
    }
    if converged {
        return Some((x[0], x[1], x[2]));
    }
    let resid = cs_f(c, s, x[0], x[1], x[2])?.iter().map(|x| x.abs()).sum::<f64>();
    if resid < base * 1e-2 + 1e-6 {
        Some((x[0], x[1], x[2]))
    } else {
        None
    }
}

/// All local extrema of the curve-surface distance via grid seeding + Newton.
pub(crate) fn curve_surface_newton_all(c: &dyn Curve, s: &dyn Surface) -> Vec<ExtremaPair> {
    let (t0, t1) = curve_bound(c);
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let (nt, nu, nv) = (15, 11, 11);
    let mut d2 = vec![vec![vec![0.0; nv + 1]; nu + 1]; nt + 1];
    for i in 0..=nt {
        for j in 0..=nu {
            for k in 0..=nv {
                let t = t0 + (t1 - t0) * i as f64 / nt as f64;
                let u = u0 + (u1 - u0) * j as f64 / nu as f64;
                let v = v0 + (v1 - v0) * k as f64 / nv as f64;
                let d = c.d0(t).square_distance(&s.d0(u, v));
                d2[i][j][k] = if d.is_finite() { d } else { f64::INFINITY };
            }
        }
    }
    let at = |i: usize, j: usize, k: usize| (t0 + (t1 - t0) * i as f64 / nt as f64, u0 + (u1 - u0) * j as f64 / nu as f64, v0 + (v1 - v0) * k as f64 / nv as f64);
    let mut seeds: Vec<(f64, f64, f64)> = Vec::new();
    for i in 1..nt {
        for j in 1..nu {
            for k in 1..nv {
                let d = d2[i][j][k];
                let lo = d2[i - 1][j][k] >= d
                    && d2[i + 1][j][k] >= d
                    && d2[i][j - 1][k] >= d
                    && d2[i][j + 1][k] >= d
                    && d2[i][j][k - 1] >= d
                    && d2[i][j][k + 1] >= d;
                let hi = d2[i - 1][j][k] <= d
                    && d2[i + 1][j][k] <= d
                    && d2[i][j - 1][k] <= d
                    && d2[i][j + 1][k] <= d
                    && d2[i][j][k - 1] <= d
                    && d2[i][j][k + 1] <= d;
                if lo || hi {
                    seeds.push(at(i, j, k));
                }
            }
        }
    }
    let (mut gmin, mut gmax) = ((0usize, 0usize, 0usize), (0usize, 0usize, 0usize));
    for i in 0..=nt {
        for j in 0..=nu {
            for k in 0..=nv {
                if d2[i][j][k] < d2[gmin.0][gmin.1][gmin.2] {
                    gmin = (i, j, k);
                }
                if d2[i][j][k] > d2[gmax.0][gmax.1][gmax.2] {
                    gmax = (i, j, k);
                }
            }
        }
    }
    seeds.push(at(gmin.0, gmin.1, gmin.2));
    seeds.push(at(gmax.0, gmax.1, gmax.2));

    let mut out: Vec<ExtremaPair> = Vec::new();
    for (t, u, v) in seeds {
        if let Some((tt, uu, vv)) =
            solve_curve_surface(c, s, t, u, v, t0, t1, u0, u1, v0, v1)
        {
            let pc = c.d0(tt);
            let ps = s.d0(uu, vv);
            if !pc.x().is_finite() || !ps.x().is_finite() {
                continue;
            }
            if !out
                .iter()
                .any(|e: &ExtremaPair| e.p1.distance(&pc) < 1e-6 && e.p2.distance(&ps) < 1e-6)
            {
                out.push(cs_pair(pc, tt, ps, uu, vv, pc.distance(&ps)));
            }
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

/// Fallback: grid min (degenerate curve-surface pairs).
pub(super) fn fallback_curve_surface(c: &dyn Curve, s: &dyn Surface) -> ExtremaPair {
    let (t0, t1) = curve_bound(c);
    let (u0, u1) = surf_bound_u(s);
    let (v0, v1) = surf_bound_v(s);
    let (nt, nu, nv) = (20, 16, 16);
    let mut best = (0.0, 0.0, 0.0, f64::INFINITY);
    for i in 0..=nt {
        for j in 0..=nu {
            for k in 0..=nv {
                let t = t0 + (t1 - t0) * i as f64 / nt as f64;
                let u = u0 + (u1 - u0) * j as f64 / nu as f64;
                let v = v0 + (v1 - v0) * k as f64 / nv as f64;
                let d = c.d0(t).square_distance(&s.d0(u, v));
                if d < best.3 {
                    best = (t, u, v, d);
                }
            }
        }
    }
    let pc = c.d0(best.0);
    let ps = s.d0(best.1, best.2);
    cs_pair(pc, best.0, ps, best.1, best.2, pc.distance(&ps))
}

// ---------------------------------------------------------------------------
// Public dispatch.
// ---------------------------------------------------------------------------

/// All local extrema of the point-surface distance, deduplicated and sorted.
/// Planes and spheres are classified and solved analytically (exact); every
/// other surface goes through the grid + Newton path.
pub fn point_surface_extrema_all(s: &dyn Surface, p: &GpPnt) -> Vec<ExtremaPair> {
    if let Some(pl) = classify_plane(s) {
        return vec![point_plane_extrema(&pl, p)];
    }
    if let Some(sp) = classify_sphere(s) {
        let mut v = point_sphere_extrema(&sp, p);
        v.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
        return v;
    }
    point_surface_newton_all(s, p)
}

/// Minimum distance from `p` to `s` (with the closest point and parameters).
pub fn point_surface_extrema(s: &dyn Surface, p: &GpPnt) -> ExtremaPair {
    match point_surface_extrema_all(s, p).into_iter().next() {
        Some(e) => e,
        None => fallback_point_surface(s, p),
    }
}

/// All local extrema of the curve-surface distance, deduplicated and sorted.
/// Lines against planes/spheres are classified and solved analytically; all
/// other pairs go through the grid + Newton path.
pub fn curve_surface_extrema_all(c: &dyn Curve, s: &dyn Surface) -> Vec<ExtremaPair> {
    if let Some(l) = reconstruct_line(c) {
        if let Some(sp) = classify_sphere(s) {
            return line_sphere_extrema(&l, &sp);
        }
        if let Some(pl) = classify_plane(s) {
            return line_plane_extrema(&l, &pl);
        }
    }
    curve_surface_newton_all(c, s)
}

/// Minimum distance between curve `c` and surface `s`.
pub fn curve_surface_extrema(c: &dyn Curve, s: &dyn Surface) -> ExtremaPair {
    match curve_surface_extrema_all(c, s).into_iter().next() {
        Some(e) => e,
        None => fallback_curve_surface(c, s),
    }
}
