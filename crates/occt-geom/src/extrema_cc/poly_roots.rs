use super::prelude::*;
use super::*;

pub(super) const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Polynomial root helpers.
//
// These are the ported `math_DirectPolynomialRoots`
// (`occt-core/src/math_direct_poly_roots.rs`) consumed as the OCCT
// `Value(1 .. NbSolutions)` list. The local closed-form quadratic/cubic, the
// Ferrari quartic and their `poly4`/`polish4` helpers were a third copy of the
// same invented solver and are gone.
// ---------------------------------------------------------------------------

pub(super) use occt_core::math_direct_poly_roots::{cubic_roots, quadratic_roots, quartic_roots};

/// Push `u + 2πk` for every integer k that lands in [uinf, usup].
pub(super) fn push_periodic(roots: &mut Vec<f64>, u: f64, uinf: f64, usup: f64) {
    if !(uinf.is_finite() && usup.is_finite()) {
        return;
    }
    let period = 2.0 * PI;
    let k0 = ((uinf - u) / period).ceil() as i64;
    let k1 = ((usup - u) / period).floor() as i64;
    for k in k0..=k1 {
        roots.push(u + k as f64 * period);
    }
}

/// Evaluate `cc·cos² + 2·sc·cos·sin + c·cos + s·sin + cte` and its derivative.
pub(super) fn trig_poly(cc: f64, sc: f64, c: f64, s: f64, cte: f64, u: f64) -> (f64, f64) {
    let (su, cu) = u.sin_cos();
    let g = cc * cu * cu + 2.0 * sc * cu * su + c * cu + s * su + cte;
    let gp = -2.0 * cc * cu * su + 2.0 * sc * (cu * cu - su * su) - c * su + s * cu;
    (g, gp)
}

/// Newton-polish a root of the full trigonometric polynomial.
pub(super) fn polish_trig_full(cc: f64, sc: f64, c: f64, s: f64, cte: f64, mut u: f64) -> f64 {
    for _ in 0..8 {
        let (g, gp) = trig_poly(cc, sc, c, s, cte, u);
        if gp.abs() < 1e-300 {
            break;
        }
        let un = u - g / gp;
        if !un.is_finite() || (un - u).abs() < 1e-13 * (1.0 + u.abs()) {
            break;
        }
        u = un;
    }
    u
}

/// All roots of `cc·cos² + 2·sc·cos·sin + c·cos + s·sin + cte = 0` in
/// [uinf, usup].  Substituting t = tan(u/2) yields the quartic
/// (cc−c+cte)·t⁴ + (−4sc+2s)·t³ + (−2cc+2cte)·t² + (4sc+2s)·t + (cc+c+cte) = 0;
/// u = π (t = ∞) is a root iff cc−c+cte = 0.
pub(super) fn trig_roots_full(cc: f64, sc: f64, c: f64, s: f64, cte: f64, uinf: f64, usup: f64) -> Vec<f64> {
    let eps = 1e-13;
    let a = cc - c + cte;
    let b = -4.0 * sc + 2.0 * s;
    let c2 = -2.0 * cc + 2.0 * cte;
    let d = 4.0 * sc + 2.0 * s;
    let e = cc + c + cte;
    let mut roots: Vec<f64> = Vec::new();
    for t in quartic_roots(a, b, c2, d, e) {
        push_periodic(&mut roots, 2.0 * t.atan(), uinf, usup);
    }
    if a.abs() < eps {
        push_periodic(&mut roots, PI, uinf, usup);
    }
    for r in roots.iter_mut() {
        *r = polish_trig_full(cc, sc, c, s, cte, *r);
    }
    let amax = cc.abs().max(sc.abs()).max(c.abs()).max(s.abs()).max(cte.abs());
    let tol = 1e-8 * (1.0 + amax);
    roots.retain(|u| *u >= uinf - 1e-8 && *u <= usup + 1e-8 && trig_poly(cc, sc, c, s, cte, *u).0.abs() < tol);
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    roots.dedup_by(|x, y| (*x - *y).abs() < 1e-8);
    roots
}

// ---------------------------------------------------------------------------
// Analytic solvers (port of `Extrema_ExtElC`).
// ---------------------------------------------------------------------------

pub(super) fn pair_cc(p1: GpPnt, u1: f64, p2: GpPnt, u2: f64) -> ExtremaPair {
    ExtremaPair { p1, p2, distance: p1.distance(&p2), u1, v1: None, u2, v2: None }
}

pub(super) fn circle_param(c: &GpCirc, p: &GpPnt) -> f64 {
    let o = c.location();
    let v = GpVec::from_pnts(&o, p);
    let xd = *c.x_axis().direction();
    let yd = *c.y_axis().direction();
    // C(u) = O + r·cos·xd + r·sin·yd → u = atan2(v·yd, v·xd).
    v.dot(&GpVec::from_xyz(yd.xyz())).atan2(v.dot(&GpVec::from_xyz(xd.xyz())))
}

pub(super) fn line_value(l: &GpLin, u: f64) -> GpPnt {
    clib::line_value(l, u)
}

pub(super) fn circle_value(c: &GpCirc, u: f64) -> GpPnt {
    clib::circle_value(c, u)
}

/// Exact extrema between two (unbounded) straight lines.
///
/// Port of `Extrema_ExtElC(const gp_Lin&, const gp_Lin&, AngTol)`. Skew lines
/// have a single closest pair; parallel lines have constant distance.
pub fn line_line_extrema(l1: &GpLin, l2: &GpLin) -> Vec<ExtremaPair> {
    let d1 = l1.direction();
    let d2 = l2.direction();
    let cosa = d1.dot(&d2);
    let sqsina = 1.0 - cosa * cosa;
    if sqsina < RESOLUTION || d1.is_parallel(&d2) {
        // Parallel: distance is the distance from any point of C1 to C2.
        let p1 = l1.location();
        let u2 = GpVec::from_pnts(&l2.location(), &p1).dot(&GpVec::from_xyz(d2.xyz()));
        let p2 = line_value(l2, u2);
        return vec![pair_cc(p1, 0.0, p2, u2)];
    }
    let l1l2 = GpVec::from_pnts(&l1.location(), &l2.location());
    let d1l = GpVec::from_xyz(d1.xyz()).dot(&l1l2);
    let d2l = GpVec::from_xyz(d2.xyz()).dot(&l1l2);
    let u1 = (d1l - cosa * d2l) / sqsina;
    let u2 = (cosa * d1l - d2l) / sqsina;
    let p1 = line_value(l1, u1);
    let p2 = line_value(l2, u2);
    vec![pair_cc(p1, u1, p2, u2)]
}

/// Extrema between a straight line and a circle.
///
/// Port of `Extrema_ExtElC(const gp_Lin&, const gp_Circ&, Tol)`. Lines parallel
/// to the circle plane are solved in the plane (2D line–circle plus
/// intersections); otherwise the full 3D trigonometric condition
/// `A1·cos² + 2·A2·cos·sin + A3·cos + A4·sin + A5 = 0` is solved.
pub fn line_circle_extrema(l: &GpLin, c: &GpCirc) -> Vec<ExtremaPair> {
    let dir_c = *c.axis().direction();
    let dir_l = l.direction();
    let o1 = l.location();
    let o2 = c.location();
    let r = c.radius();

    // Planar (line parallel to the circle plane): solve in the circle frame.
    if dir_c.dot(&dir_l).abs() <= ANGULAR {
        return planar_line_circle(l, c);
    }

    // Full 3D case. Work in the circle's reference frame.
    let x2 = *c.x_axis().direction();
    let y2 = *c.y_axis().direction();
    let dx = dir_l.dot(&x2);
    let dy = dir_l.dot(&y2);
    let _dz = dir_l.dot(&dir_c);

    // O2O1 in the circle frame.
    let o2o1 = GpVec::from_pnts(&o2, &o1);
    let vx = o2o1.dot(&GpVec::from_xyz(x2.xyz()));
    let vy = o2o1.dot(&GpVec::from_xyz(y2.xyz()));

    // V = (D·O2O1)·D − O2O1 in the circle frame.
    let do2o1 = dx * vx + dy * vy;
    let vx = do2o1 * dx - vx;
    let vy = do2o1 * dy - vy;

    // Coefficients of the trigonometric equation (divided by R for A3/A4).
    let a5 = r * dx * dy;
    let a1 = -2.0 * a5;
    let a2 = 0.5 * r * (dx * dx - dy * dy);
    let a3 = vy;
    let a4 = -vx;

    // Infinite solutions: line is the circle axis (constant distance R).
    if (a1.abs() + a2.abs() + a3.abs() + a4.abs()) < 1e-10 && a5.abs() < 1e-10 {
        let p2 = circle_value(c, 0.0);
        let u1 = GpVec::from_pnts(&o1, &p2).dot(&GpVec::from_xyz(dir_l.xyz()));
        let p1 = line_value(l, u1);
        return vec![pair_cc(p1, u1, p2, 0.0)];
    }

    let us = trig_roots_full(a1, a2, a3, a4, a5, 0.0, 2.0 * PI);
    let mut out = Vec::with_capacity(us.len());
    for u2 in us {
        let p2 = circle_value(c, u2);
        let u1 = GpVec::from_pnts(&o1, &p2).dot(&GpVec::from_xyz(dir_l.xyz()));
        let p1 = line_value(l, u1);
        out.push(pair_cc(p1, u1, p2, u2));
    }
    dedupe_sort(out)
}

/// Planar line–circle (port of `Extrema_ExtElC::PlanarLineCircleExtrema`).
///
/// The line is projected onto the circle plane; the 2D tangent extrema (radius
/// perpendicular to the projected direction) and the 2D intersections are both
/// mapped back to the actual 3D line/circle parameters (the "intersections" of
/// the projection are candidate closest pairs in 3D, not literal crossings).
pub(super) fn planar_line_circle(l: &GpLin, c: &GpCirc) -> Vec<ExtremaPair> {
    let o1 = l.location();
    let o2 = c.location();
    let r = c.radius();
    let x2 = *c.x_axis().direction();
    let y2 = *c.y_axis().direction();

    // Project the line onto the circle plane.
    let dl = GpVec::from_xyz(l.direction().xyz());
    let p2dx = GpVec::from_xyz(x2.xyz());
    let p2dy = GpVec::from_xyz(y2.xyz());
    let ldx = dl.dot(&p2dx);
    let ldy = dl.dot(&p2dy);
    let o2o1 = GpVec::from_pnts(&o2, &o1);
    let loc2d_x = o2o1.dot(&p2dx);
    let loc2d_y = o2o1.dot(&p2dy);

    let mut out = Vec::new();

    // 2D tangent extrema: radius perpendicular to the projected line direction.
    let mut tetas: Vec<f64> = Vec::new();
    if ldy.abs() <= 1e-15 {
        tetas.push(PI / 2.0);
    } else {
        tetas.push((-ldx / ldy).atan());
    }
    tetas.push(tetas[0] + PI);
    if tetas[0] < 0.0 {
        tetas[0] += 2.0 * PI;
    }
    for teta in tetas {
        let p2 = circle_value(c, teta);
        let u1 = GpVec::from_pnts(&o1, &p2).dot(&dl);
        let p1 = line_value(l, u1);
        out.push(pair_cc(p1, u1, p2, teta));
    }

    // 2D projected intersections: |(loc2d_x, loc2d_y) + t·(ldx, ldy)|² = r².
    let a = ldx * ldx + ldy * ldy;
    let b = 2.0 * (loc2d_x * ldx + loc2d_y * ldy);
    let c0 = loc2d_x * loc2d_x + loc2d_y * loc2d_y - r * r;
    let disc = b * b - 4.0 * a * c0;
    if disc >= 0.0 && a > 1e-300 {
        let sq = disc.sqrt();
        for t in [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
            let p1 = line_value(l, t);
            let x2d = loc2d_x + t * ldx;
            let y2d = loc2d_y + t * ldy;
            let u2 = y2d.atan2(x2d);
            let p2 = circle_value(c, u2);
            out.push(pair_cc(p1, t, p2, u2));
        }
    }
    dedupe_sort(out)
}

/// Extrema between two circles.
///
/// Port of `Extrema_ExtElC(const gp_Circ&, const gp_Circ&)`. Only same-plane
/// circles are handled exactly (coaxial constant-distance, the four collinear
/// critical pairs, and intersections). Skew circles return an empty vector and
/// must fall back to the general Newton path.
pub fn circle_circle_extrema(c1: &GpCirc, c2: &GpCirc) -> Vec<ExtremaPair> {
    let oc1 = c1.location();
    let oc2 = c2.location();
    let dc1 = *c1.axis().direction();
    let dc2 = *c2.axis().direction();
    let tol_d2 = CONFUSION * CONFUSION;

    // Same-plane check.
    let plane_dist2 = {
        let v = GpVec::from_pnts(&oc1, &oc2);
        let n = GpVec::from_xyz(dc1.xyz());
        let d = v.dot(&n);
        d * d
    };
    if !(dc1.is_parallel(&dc2) && plane_dist2 < tol_d2) {
        return Vec::new(); // skew circles
    }

    let r1 = c1.radius();
    let r2 = c2.radius();

    // Coaxial (concentric): constant distance |r1 − r2|.
    let d12 = oc1.distance(&oc2);
    if d12 < CONFUSION {
        let u = GpVec::from_xyz(c1.x_axis().direction().xyz());
        let p1 = oc1.translated_vec(&u.multiplied_scalar(r1));
        let p2 = oc2.translated_vec(&u.multiplied_scalar(r2));
        return vec![pair_cc(p1, circle_param(c1, &p1), p2, circle_param(c2, &p2))];
    }

    let dir12 = GpVec::from_pnts(&oc1, &oc2).multiplied_scalar(1.0 / d12);
    let mut out = Vec::with_capacity(6);

    // Four collinear critical pairs.
    let p11 = oc1.translated_vec(&dir12.multiplied_scalar(-r1));
    let p12 = oc1.translated_vec(&dir12.multiplied_scalar(r1));
    let p21 = oc2.translated_vec(&dir12.multiplied_scalar(-r2));
    let p22 = oc2.translated_vec(&dir12.multiplied_scalar(r2));
    for (pa, pb) in [(p11, p21), (p11, p22), (p12, p21), (p12, p22)] {
        out.push(pair_cc(pa, circle_param(c1, &pa), pb, circle_param(c2, &pb)));
    }

    // Intersections (zero distance) when the circles cross or touch.
    let b_out = d12 > (r1 + r2 + CONFUSION);
    let b_in = d12 < (r1 - r2).abs() - CONFUSION;
    if !b_out && !b_in {
        let alpha = 0.5 * (r1 * r1 - r2 * r2 + d12 * d12) / d12;
        let val = r1 * r1 - alpha * alpha;
        let beta = val.abs().sqrt();
        let pt = oc1.translated_vec(&dir12.multiplied_scalar(alpha));
        let dlt = GpVec::from_xyz(dc1.xyz()).crossed(&dir12);
        let pl1 = pt.translated_vec(&dlt.multiplied_scalar(beta));
        let pl2 = pt.translated_vec(&dlt.multiplied_scalar(-beta));
        if pl1.square_distance(&pl2) > tol_d2 {
            out.push(pair_cc(pl1, circle_param(c1, &pl1), pl1, circle_param(c2, &pl1)));
            out.push(pair_cc(pl2, circle_param(c1, &pl2), pl2, circle_param(c2, &pl2)));
        } else {
            out.push(pair_cc(pl1, circle_param(c1, &pl1), pl1, circle_param(c2, &pl1)));
        }
    }
    dedupe_sort(out)
}

