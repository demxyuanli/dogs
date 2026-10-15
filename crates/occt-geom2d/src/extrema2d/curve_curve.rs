use super::prelude::*;
use super::*;

/// `Extrema_ExtPC2d` solution turned into the port's `Extrema2d` pair carrier:
/// `u1`/`u2` both hold the curve parameter, `p1` the query point.
fn point_solution2d(p: &GpPnt2d, s: &ExtPc2dSolution) -> Extrema2d {
    Extrema2d { p1: *p, p2: s.point, distance: p.distance(&s.point), u1: s.u, u2: s.u }
}

/// All extremal distances between a point and a 2D curve, in the order
/// `Extrema_ExtPC2d` reports them.
///
/// Faithful route: [`ExtremaExtPC2d`] = `Extrema_GGExtPC<Adaptor2d_Curve2d, ...>`
/// (`Extrema_ExtPC2d.hxx:31-38`, body `Extrema_GGExtPC.hxx`). The previous body
/// was a self-invented `span/0.1` seed grid (`analytic_solvers::build_samples2d`)
/// with bracket sign changes plus Newton (`newton_point_curve_all2d`); both are
/// deleted. OCCT neither sorts nor dedups this list, and it only adds a range
/// endpoint when the query point already coincides with it
/// (`Extrema_GGExtPC.hxx:473-502`), hence the empty list on `!IsDone()`.
pub fn point_curve_extrema2d_all(c: &dyn Curve2d, p: &GpPnt2d) -> Vec<Extrema2d> {
    let ext = ExtremaExtPC2d::new(c, p);
    if !ext.is_done() {
        return Vec::new();
    }
    ext.solutions().iter().map(|s| point_solution2d(p, s)).collect()
}

/// Closest point of `p` on `c` as `Geom2dAPI_ProjectPointOnCurve::Perform`
/// reports it (`Geom2dAPI_ProjectPointOnCurve.cxx:59-86`): `Extrema_ExtPC2d`,
/// done only when `IsDone() && NbExt() > 0`, then the smallest
/// `SquareDistance(i)` with the first index winning ties. `None` is the
/// not-done/no-extremum case.
pub fn point_curve_extrema2d(c: &dyn Curve2d, p: &GpPnt2d) -> Option<Extrema2d> {
    let ext = ExtremaExtPC2d::new(c, p);
    if !ext.is_done() || ext.nb_ext() == 0 {
        return None;
    }
    let mut dist2_min = ext.square_distance(1);
    let mut index = 1;
    for i in 2..=ext.nb_ext() {
        let dist2 = ext.square_distance(i);
        if dist2 < dist2_min {
            dist2_min = dist2;
            index = i;
        }
    }
    Some(point_solution2d(p, &ext.solutions()[index - 1]))
}

/// Minimum distance between two curves (all extrema, first = closest).
///
/// Delegates to `curve_curve_extrema2d_all`. `samples` is kept for signature
/// compatibility but no longer drives the computation.
pub fn curve_curve_extrema2d(c1: &dyn Curve2d, c2: &dyn Curve2d, _samples: usize) -> Vec<Extrema2d> {
    curve_curve_extrema2d_all(c1, c2)
}

// ---------------------------------------------------------------------------
// Curve–curve extrema (port of `Extrema_ExtCC2d`, `Extrema_ExtElC2d`).
// ---------------------------------------------------------------------------

pub(super) fn pair2d_cc(p1: GpPnt2d, u1: f64, p2: GpPnt2d, u2: f64) -> Extrema2d {
    Extrema2d { p1, p2, distance: p1.distance(&p2), u1, u2 }
}

pub(super) fn line2d_value(l: &GpLin2d, u: f64) -> GpPnt2d {
    GpPnt2d::new(l.pos.loc.x() + u * l.pos.vdir.x, l.pos.loc.y() + u * l.pos.vdir.y)
}

pub(super) fn circle_param2d(c: &GpCirc2d, p: &GpPnt2d) -> f64 {
    let v = GpVec2d::new(p.x() - c.location().x(), p.y() - c.location().y());
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    // C(u) = O + r·cos·xd + r·sin·yd → u = atan2(v·yd, v·xd).
    v.dot(&GpVec2d::new(yd.x, yd.y)).atan2(v.dot(&GpVec2d::new(xd.x, xd.y)))
}

/// Exact extrema between two 2D lines (port of `Extrema_ExtElC2d`).
/// Parallel lines have constant distance; otherwise they intersect (distance 0).
pub fn line_line_extrema2d(l1: &GpLin2d, l2: &GpLin2d) -> Vec<Extrema2d> {
    let d1 = GpVec2d::from_dir2d(&l1.pos.vdir);
    let d2 = GpVec2d::from_dir2d(&l2.pos.vdir);
    if d1.is_parallel(&d2, ANGULAR) {
        let p1 = l1.pos.loc;
        let u2 = (p1.x() - l2.pos.loc.x()) * d2.x() + (p1.y() - l2.pos.loc.y()) * d2.y();
        let p2 = line2d_value(l2, u2);
        return vec![pair2d_cc(p1, 0.0, p2, u2)];
    }
    // Cramer's rule for the intersection.
    let a_p1p2 = GpVec2d::new(l2.pos.loc.x() - l1.pos.loc.x(), l2.pos.loc.y() - l1.pos.loc.y());
    let delim = 1.0 / d1.crossed(&d2);
    let param1 = a_p1p2.crossed(&d2) * delim;
    let param2 = -(d1.crossed(&a_p1p2)) * delim;
    let p1 = line2d_value(l1, param1);
    let p2 = line2d_value(l2, param2);
    vec![pair2d_cc(p1, param1, p2, param2)]
}

/// Exact extrema between a 2D line and a 2D circle (port of
/// `Extrema_ExtElC2d(const gp_Lin2d&, const gp_Circ2d&)`), plus the
/// zero-distance intersections.
pub fn line_circle_extrema2d(l: &GpLin2d, c: &GpCirc2d) -> Vec<Extrema2d> {
    let d = l.pos.vdir;
    let x2 = c.pos.vxdir;
    let y2 = c.pos.vydir;
    let dx = d.dot(&x2);
    let dy = d.dot(&y2);
    let o1 = l.pos.loc;

    // Radius perpendicular to the line direction.
    let mut tetas: Vec<f64> = Vec::new();
    if dy.abs() <= 1e-15 {
        tetas.push(PI / 2.0);
    } else {
        tetas.push((-dx / dy).atan());
    }
    tetas.push(tetas[0] + PI);
    if tetas[0] < 0.0 {
        tetas[0] += 2.0 * PI;
    }

    let mut out = Vec::new();
    for teta in tetas {
        let p2 = circle2d_val(c, teta);
        let u1 = (p2.x() - o1.x()) * d.x + (p2.y() - o1.y()) * d.y;
        let p1 = line2d_value(l, u1);
        out.push(pair2d_cc(p1, u1, p2, teta));
    }

    // Intersections: |O1 + t·D − O2|² = r².
    let o2 = c.location();
    let r = c.radius();
    let v12 = GpVec2d::new(o1.x() - o2.x(), o1.y() - o2.y());
    let a = d.x * d.x + d.y * d.y;
    let b = 2.0 * (v12.x() * d.x + v12.y() * d.y);
    let c0 = v12.square_magnitude() - r * r;
    let disc = b * b - 4.0 * a * c0;
    if disc >= 0.0 && a > 1e-300 {
        let sq = disc.sqrt();
        for t in [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
            let p1 = line2d_value(l, t);
            let u2 = circle_param2d(c, &p1);
            out.push(pair2d_cc(p1, t, p1, u2));
        }
    }
    dedupe_sort_cc2d(out)
}

/// Exact extrema between two 2D circles (port of
/// `Extrema_ExtElC2d(const gp_Circ2d&, const gp_Circ2d&)`), plus the
/// zero-distance intersections for crossing circles.
pub fn circle_circle_extrema2d(c1: &GpCirc2d, c2: &GpCirc2d) -> Vec<Extrema2d> {
    let o1 = c1.location();
    let o2 = c2.location();
    let r1 = c1.radius();
    let r2 = c2.radius();

    // Concentric: constant distance |r1 − r2|.
    let d12 = o1.distance(&o2);
    if d12 < CONFUSION {
        let u = GpVec2d::new(c1.pos.vxdir.x, c1.pos.vxdir.y);
        let p1 = o1.translated_vec(&u.multiplied_scalar(r1));
        let p2 = o2.translated_vec(&u.multiplied_scalar(r2));
        return vec![pair2d_cc(p1, circle_param2d(c1, &p1), p2, circle_param2d(c2, &p2))];
    }

    let dir = GpVec2d::new(o2.x() - o1.x(), o2.y() - o1.y()).multiplied_scalar(1.0 / d12);
    let mut out = Vec::with_capacity(6);

    // Four collinear critical pairs.
    let p11 = o1.translated_vec(&dir.multiplied_scalar(-r1));
    let p12 = o1.translated_vec(&dir.multiplied_scalar(r1));
    let p21 = o2.translated_vec(&dir.multiplied_scalar(-r2));
    let p22 = o2.translated_vec(&dir.multiplied_scalar(r2));
    for (pa, pb) in [(p11, p21), (p11, p22), (p12, p21), (p12, p22)] {
        out.push(pair2d_cc(pa, circle_param2d(c1, &pa), pb, circle_param2d(c2, &pb)));
    }

    // Intersections (zero distance).
    let b_out = d12 > (r1 + r2 + CONFUSION);
    let b_in = d12 < (r1 - r2).abs() - CONFUSION;
    if !b_out && !b_in {
        let alpha = 0.5 * (r1 * r1 - r2 * r2 + d12 * d12) / d12;
        let val = r1 * r1 - alpha * alpha;
        let beta = val.abs().sqrt();
        let pt = o1.translated_vec(&dir.multiplied_scalar(alpha));
        // Perpendicular direction in the plane.
        let dlt = GpVec2d::new(-dir.y(), dir.x());
        let pl1 = pt.translated_vec(&dlt.multiplied_scalar(beta));
        let pl2 = pt.translated_vec(&dlt.multiplied_scalar(-beta));
        if pl1.square_distance(&pl2) > CONFUSION * CONFUSION {
            out.push(pair2d_cc(pl1, circle_param2d(c1, &pl1), pl1, circle_param2d(c2, &pl1)));
            out.push(pair2d_cc(pl2, circle_param2d(c1, &pl2), pl2, circle_param2d(c2, &pl2)));
        } else {
            out.push(pair2d_cc(pl1, circle_param2d(c1, &pl1), pl1, circle_param2d(c2, &pl1)));
        }
    }
    dedupe_sort_cc2d(out)
}

/// Deduplicate curve–curve extrema by parameter proximity, sorted by distance.
pub(super) fn dedupe_sort_cc2d(v: Vec<Extrema2d>) -> Vec<Extrema2d> {
    let mut out: Vec<Extrema2d> = Vec::new();
    for e in v {
        if let Some(o) = out
            .iter_mut()
            .find(|o| (o.u1 - e.u1).abs() < 1e-5 && (o.u2 - e.u2).abs() < 1e-5)
        {
            if e.distance < o.distance {
                *o = e;
            }
        } else {
            out.push(e);
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

/// Map `u` into `[lo, lo+period)` for a periodic curve (OCCT `ElCLib::InPeriod`).
pub(super) fn in_period2d(u: f64, lo: f64, period: f64) -> f64 {
    if period <= 0.0 || !lo.is_finite() {
        u
    } else {
        lo + (u - lo).rem_euclid(period)
    }
}

/// All extremal distances between two planar curves, in the order
/// `Extrema_ExtCC2d` reports them.
///
/// Faithful route: [`ExtremaExtCC2d`] = `Extrema_ExtCC2d` (`Extrema_ExtCC2d.hxx`,
/// `Extrema_ExtCC2d.cxx:104-145` for the `Extrema_ExtElC2d` dispatch and
/// `:435-451` for the general arm) with its engine `Extrema_ECC2d` =
/// `Extrema_GGenExtCC`. The previous body re-implemented that dispatch locally
/// and fell back to the self-invented seed grid plus Newton
/// (`newton_curve_curve_all2d`), which is deleted.
pub fn curve_curve_extrema2d_all(c1: &dyn Curve2d, c2: &dyn Curve2d) -> Vec<Extrema2d> {
    let ext = ExtremaExtCC2d::new(c1, c2);
    (1..=ext.nb_ext())
        .map(|n| {
            let (u1, p1, u2, p2) = ext.points(n);
            pair2d_cc(p1, u1, p2, u2)
        })
        .collect()
}

/// Local extremum of the distance between two planar curves from seed
/// parameters. Port of `Extrema_LocateExtCC2d`
/// (`Extrema_LocateExtCC2d.cxx:26-82`), which runs
/// `Extrema_LocECC2d` = `Extrema_GenLocateExtCC` -> `Extrema_CCLocFOfLocECC2d`
/// (`Extrema_GFuncExtCC` 2D) under `math_FunctionSetRoot`. The previous body was
/// a self-invented Newton loop (`refine_curve_curve2d_newton` /
/// `f_and_jac2d`), both deleted.
pub fn locate_extcc2d(c1: &dyn Curve2d, c2: &dyn Curve2d, u0: f64, v0: f64) -> Option<Extrema2d> {
    let loc = ExtremaLocateExtCC2d::new(c1, c2, u0, v0);
    if !loc.is_done() {
        return None;
    }
    let (p1, p2) = loc.point();
    Some(pair2d_cc(p1.value(), p1.parameter(), p2.value(), p2.parameter()))
}

/// All curve–curve intersection points (zero-distance extrema).
pub fn curve_curve_intersections2d(c1: &dyn Curve2d, c2: &dyn Curve2d, tol: f64) -> Vec<GpPnt2d> {
    let (a1, b1) = bound(c1);
    let (a2, b2) = bound(c2);
    let n = 32;
    let mut pts: Vec<GpPnt2d> = Vec::new();
    for i in 0..=n {
        let u = a1 + (b1 - a1) * i as f64 / n as f64;
        for j in 0..=n {
            let v = a2 + (b2 - a2) * j as f64 / n as f64;
            let p = c1.d0(u);
            let q = c2.d0(v);
            if p.distance(&q) < tol {
                if pts.iter().all(|x| x.distance(&p) > 1e-6) {
                    pts.push(p);
                }
            }
        }
    }
    pts
}

/// Distance from a point to the closest point on a polyline.
pub fn point_polyline_extrema2d(poly: &[GpPnt2d], p: &GpPnt2d) -> Extrema2d {
    let mut best: Option<Extrema2d> = None;
    for w in poly.windows(2) {
        let (a, b) = (w[0], w[1]);
        let abx = b.x() - a.x();
        let aby = b.y() - a.y();
        let len2 = abx * abx + aby * aby;
        let t = if len2 < 1e-30 {
            0.0
        } else {
            (((p.x() - a.x()) * abx + (p.y() - a.y()) * aby) / len2).clamp(0.0, 1.0)
        };
        let q = GpPnt2d::new(a.x() + t * abx, a.y() + t * aby);
        let d = p.distance(&q);
        if best.as_ref().map_or(true, |e: &Extrema2d| d < e.distance) {
            best = Some(Extrema2d { p1: *p, p2: q, distance: d, u1: t, u2: t });
        }
    }
    best.unwrap()
}

/// Tangent direction of a 2D curve at `u` (normalized).
pub fn tangent2d(c: &dyn Curve2d, u: f64) -> GpVec2d {
    let (_, d) = c.d1(u);
    let m = d.magnitude();
    if m > 1e-12 {
        GpVec2d::new(d.x() / m, d.y() / m)
    } else {
        GpVec2d::new(1.0, 0.0)
    }
}
