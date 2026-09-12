//! `GeomInt_IntSS::BuildPCurves`, `TreatRLine`, `TrimILineOnSurfBoundaries`.

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt2d;
use occt_core::precision::{PCONFUSION, Precision, RESOLUTION};
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::Geom2dBSplineCurve;

use crate::geom_int::adjust_periodic;
use crate::geom_int::types::RLine;
use crate::geom_int::intss_bspline::make_bspline;

/// `GeomInt_IntSS::BuildPCurves` with explicit UV box.
pub fn build_pcurves_box(
    the_first: f64,
    the_last: f64,
    the_umin: f64,
    the_umax: f64,
    the_vmin: f64,
    the_vmax: f64,
    the_tol: &mut f64,
    the_surface: &dyn Surface,
    the_curve: &dyn Curve,
) -> Option<Arc<dyn Curve2d>> {
    if the_last - the_first > 2.0e-9 {
        let n = 32usize.max(2);
        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);
        for i in 0..n {
            let t = the_first + (the_last - the_first) * (i as f64 / (n - 1) as f64);
            let p = the_curve.d0(t);
            let proj = project_point_on_surface(the_surface, &p, 0.0)?;
            xs.push(proj.u);
            ys.push(proj.v);
        }
        let mut knots = Vec::with_capacity(n + 2);
        knots.push(the_first);
        for i in 0..n {
            let t = the_first + (the_last - the_first) * (i as f64 / (n - 1) as f64);
            knots.push(t);
        }
        knots.push(the_last);
        let c2d = Geom2dBSplineCurve::new(xs, ys, knots, 1).ok()?;
        let mut out: Arc<dyn Curve2d> = Arc::new(c2d);
        if the_surface.is_u_periodic() {
            recadre_u_periodic(
                the_first,
                the_last,
                the_umin,
                the_umax,
                the_surface,
                &mut out,
            );
        }
        let _ = (the_vmin, the_vmax);
        Some(out)
    } else if (the_last - the_first) > epsilon_abs(the_first) {
        let p3d1 = the_curve.d0(the_first);
        let p3d2 = the_curve.d0(the_last);
        let p2d1 = project_point_on_surface(the_surface, &p3d1, 0.0)?;
        let p2d2 = project_point_on_surface(the_surface, &p3d2, 0.0)?;
        let a = GpPnt2d::new(p2d1.u, p2d1.v);
        let b = GpPnt2d::new(p2d2.u, p2d2.v);
        if a.distance(&b) <= RESOLUTION {
            return None;
        }
        let xs = vec![a.x(), b.x()];
        let ys = vec![a.y(), b.y()];
        let knots = vec![the_first, the_first, the_last, the_last];
        let c2d = Geom2dBSplineCurve::new(xs, ys, knots, 1).ok()?;
        let pmid = the_curve.d0(0.5 * (the_first + the_last));
        let pmidc = GpPnt2d::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()));
        let pc = the_surface.d0(pmidc.x(), pmidc.y());
        *the_tol = the_tol.max(pmid.distance(&pc));
        Some(Arc::new(c2d))
    } else {
        None
    }
}

fn epsilon_abs(x: f64) -> f64 {
    x.abs() * f64::EPSILON
}

/// `AdjustUPeriodic` (`GeomInt_IntSS_1.cxx:55`).
pub fn adjust_u_periodic_curve2d(a_s: &dyn Surface, a_c2d: &mut Arc<dyn Curve2d>) {
    if !a_s.is_u_periodic() {
        return;
    }
    let (umin, umax) = a_s.u_range();
    let f = a_c2d.first_parameter();
    let l = a_c2d.last_parameter();
    recadre_u_periodic(f, l, umin, umax, a_s, a_c2d);
}

fn recadre_u_periodic(
    the_first: f64,
    the_last: f64,
    the_umin: f64,
    the_umax: f64,
    the_surface: &dyn Surface,
    the_curve2d: &mut Arc<dyn Curve2d>,
) {
    let a_eps = PCONFUSION;
    let period = {
        let (a, b) = the_surface.u_range();
        if a.is_finite() && b.is_finite() {
            (b - a).abs()
        } else {
            2.0 * std::f64::consts::PI
        }
    };
    let a_tm = 0.5 * (the_first + the_last);
    let pm = the_curve2d.d0(a_tm);
    let u0 = pm.x();
    let (u0x, du) = adjust_periodic(u0, the_umin, the_umax, period, a_eps);
    if (u0x - u0).abs() > 0.0 || du.abs() > 0.0 {
        let n = 8;
        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);
        for i in 0..n {
            let t = the_first + (the_last - the_first) * (i as f64 / (n - 1) as f64);
            let p = the_curve2d.d0(t);
            xs.push(p.x() + du);
            ys.push(p.y());
        }
        let mut knots = Vec::with_capacity(n + 2);
        knots.push(the_first);
        for i in 0..n {
            knots.push(the_first + (the_last - the_first) * (i as f64 / (n - 1) as f64));
        }
        knots.push(the_last);
        if let Ok(c) = Geom2dBSplineCurve::new(xs, ys, knots, 1) {
            *the_curve2d = Arc::new(c);
        }
    }
}

/// `GeomInt_IntSS::BuildPCurves(f, l, Tol, S, C, C2d)`.
pub fn build_pcurves(
    f: f64,
    l: f64,
    tol: &mut f64,
    s: &dyn Surface,
    c: &dyn Curve,
) -> Option<Arc<dyn Curve2d>> {
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    build_pcurves_box(f, l, umin, umax, vmin, vmax, tol, s, c)
}

/// `GeomInt_IntSS::TreatRLine`.
pub fn treat_rline(
    the_rl: &RLine,
    the_hs1: &dyn Surface,
    the_hs2: &dyn Surface,
) -> (
    Option<Arc<dyn Curve>>,
    Option<Arc<dyn Curve2d>>,
    Option<Arc<dyn Curve2d>>,
    f64,
) {
    let Some(c_orig) = the_rl.c2d.clone() else {
        return (None, None, None, 0.0);
    };
    let (a_gahs, on_s1) = if the_rl.arc_on_s1 {
        (the_hs1, true)
    } else if the_rl.arc_on_s2 {
        (the_hs2, false)
    } else {
        return (None, None, None, 0.0);
    };
    let tf = the_rl.param_f.max(c_orig.first_parameter());
    let tl = the_rl.param_l.min(c_orig.last_parameter());
    if is_degenerated(a_gahs, c_orig.as_ref(), tf, tl) {
        return (None, None, None, 0.0);
    }
    let trimmed: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(c_orig.clone(), tf, tl));
    let mut the_c2d1 = if on_s1 { Some(trimmed.clone()) } else { None };
    let mut the_c2d2 = if on_s1 { None } else { Some(trimmed) };
    let mut wl = crate::int_tools_wline::WLine::new();
    let n = 32;
    for i in 0..n {
        let t = tf + (tl - tf) * (i as f64 / (n - 1) as f64);
        let uv = c_orig.d0(t);
        let p = a_gahs.d0(uv.x(), uv.y());
        wl.add(crate::int_tools_wline::PntOn2S {
            p,
            u1: uv.x(),
            v1: uv.y(),
            u2: uv.x(),
            v2: uv.y(),
        });
    }
    let the_c3d = make_bspline(&wl, 1, wl.nb_pnts());
    let mut the_tol_reached = 0.0;
    if let Some(c3) = &the_c3d {
        let mut a_tol = occt_core::precision::CONFUSION;
        let f = c3.first_parameter();
        let l = c3.last_parameter();
        if on_s1 {
            the_c2d2 = build_pcurves(f, l, &mut a_tol, the_hs2, c3.as_ref());
        } else {
            the_c2d1 = build_pcurves(f, l, &mut a_tol, the_hs1, c3.as_ref());
        }
        the_tol_reached = a_tol;
        return (Some(c3.clone()), the_c2d1, the_c2d2, the_tol_reached);
    }
    (None, the_c2d1, the_c2d2, the_tol_reached)
}

fn is_degenerated(s: &dyn Surface, c: &dyn Curve2d, tf: f64, tl: f64) -> bool {
    let p0 = c.d0(tf);
    let p1 = c.d0(tl);
    let a = s.d0(p0.x(), p0.y());
    let b = s.d0(p1.x(), p1.y());
    if a.square_distance(&b) > occt_core::precision::SQUARE_CONFUSION {
        return false;
    }
    let pm = c.d0(0.5 * (tf + tl));
    let m = s.d0(pm.x(), pm.y());
    a.square_distance(&m) <= occt_core::precision::SQUARE_CONFUSION
}

/// `GeomInt_IntSS::TrimILineOnSurfBoundaries`.
pub fn trim_iline_on_surf_boundaries(
    the_c2d1: Option<&dyn Curve2d>,
    the_c2d2: Option<&dyn Curve2d>,
    the_bound1: &BndBox2d,
    the_bound2: &BndBox2d,
    the_array: &mut Vec<f64>,
) {
    let Some((u1f, v1f, u1l, v1l)) = the_bound1.get() else {
        return;
    };
    let Some((u2f, v2f, u2l, v2l)) = the_bound2.get() else {
        return;
    };
    let an_int_tol = 10.0 * occt_core::precision::CONFUSION;
    let b1 = bounds_as_segments(u1f, v1f, u1l, v1l);
    let b2 = bounds_as_segments(u2f, v2f, u2l, v2l);
    intersect_curve_and_boundary(the_c2d1, &b1, an_int_tol, the_array);
    intersect_curve_and_boundary(the_c2d2, &b2, an_int_tol, the_array);
    the_array.sort_by(|a, b| a.total_cmp(b));
}

fn bounds_as_segments(uf: f64, vf: f64, ul: f64, vl: f64) -> Vec<(GpPnt2d, GpPnt2d)> {
    let mut out = Vec::new();
    let dy = vl - vf;
    let dx = ul - uf;
    if dy.abs() > f64::MIN_POSITIVE {
        if !Precision::is_infinite(uf) {
            out.push((GpPnt2d::new(uf, vf), GpPnt2d::new(uf, vl)));
        }
        if !Precision::is_infinite(ul) {
            out.push((GpPnt2d::new(ul, vf), GpPnt2d::new(ul, vl)));
        }
    }
    if dx.abs() > f64::MIN_POSITIVE {
        if !Precision::is_infinite(vf) {
            out.push((GpPnt2d::new(uf, vf), GpPnt2d::new(ul, vf)));
        }
        if !Precision::is_infinite(vl) {
            out.push((GpPnt2d::new(uf, vl), GpPnt2d::new(ul, vl)));
        }
    }
    out
}

fn intersect_curve_and_boundary(
    c: Option<&dyn Curve2d>,
    bounds: &[(GpPnt2d, GpPnt2d)],
    tol: f64,
    out: &mut Vec<f64>,
) {
    let Some(c) = c else {
        return;
    };
    let f = c.first_parameter();
    let l = c.last_parameter();
    if !f.is_finite() || !l.is_finite() {
        return;
    }
    let n = 48;
    for i in 0..n {
        let t0 = f + (l - f) * (i as f64 / n as f64);
        let t1 = f + (l - f) * ((i + 1) as f64 / n as f64);
        let p0 = c.d0(t0);
        let p1 = c.d0(t1);
        for &(q0, q1) in bounds {
            if let Some(t) = seg_seg_param(p0, p1, q0, q1, t0, t1, tol) {
                out.push(t);
            }
        }
    }
}

fn seg_seg_param(
    a0: GpPnt2d,
    a1: GpPnt2d,
    b0: GpPnt2d,
    b1: GpPnt2d,
    t0: f64,
    t1: f64,
    _tol: f64,
) -> Option<f64> {
    let dxa = a1.x() - a0.x();
    let dya = a1.y() - a0.y();
    let dxb = b1.x() - b0.x();
    let dyb = b1.y() - b0.y();
    let den = dxa * dyb - dya * dxb;
    if den.abs() < 1e-30 {
        return None;
    }
    let dx = b0.x() - a0.x();
    let dy = b0.y() - a0.y();
    let s = (dx * dyb - dy * dxb) / den;
    let u = (dx * dya - dy * dxa) / den;
    if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&u) {
        Some(t0 + s * (t1 - t0))
    } else {
        None
    }
}
