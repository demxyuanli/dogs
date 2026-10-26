//! Boundary helpers of `ProjLib_CompProjectedCurve`: `ExactBound`,
//! `DichExactBound` and `InitialPoint`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx:307-544`
//! (OCCT 8.0.0).
//!
//! `ExactBound` walks the four surface borders ordered by their distance from
//! the current `(u, v)` and solves the normal-projection system on each of them
//! (`myFix` 2 for a `u` border, 3 for a `v` border). `DichExactBound` is the
//! bisection fallback: it walks `t` towards `NotSol` keeping the last parameter
//! where the projection is still solvable. `InitialPoint` asks
//! `Extrema_ExtPS` for the point-surface extrema and keeps the nearest one that
//! is *also* a normal projection.
//!
//! Translation notes:
//! - `Precision::Infinite()` (`cxx:331-332`) is [`occt_core::precision::INFINITE`].
//! - the `StrictInside` argument is `true` in both `InitialPoint` call sites
//!   (`cxx:524`) and absent (default `false`) elsewhere (`cxx:389`, `cxx:406`);
//!   the port passes the same values.
//! - `Extrema_ExtPS(P, S, FirstU, LastU, FirstV, LastV, TolU, TolV,
//!   Extrema_ExtFlag_MIN)` (`cxx:488-499`) is [`ExtPs::initialize`] +
//!   [`ExtPs::set_flag`] + [`ExtPs::perform`].
//! - `aExtPS.SquareDistance(i)` / `Point(i)` (`cxx:510-540`) are
//!   [`ExtPs::square_distance`] / [`ExtPs::point`], both 1-based.
//! - `d1` raises `Standard_ConstructionError` on a singular Hessian
//!   (`cxx:158-160`); `ExactBound` and `DichExactBound` do not catch it, so the
//!   port panics the same way instead of reporting `false`.

use occt_core::gp::GpPnt;
use occt_core::precision::{INFINITE, RESOLUTION};

use crate::extrema_surf::{ExtPs, ExtremaExtFlag};
use crate::projlib::prj_resolve::PrjResolve;
use crate::{Curve, Surface};

use super::jets::d1;
use super::FUNC_TOL;

/// `ExactBound(Sol, NotSol, Tol, TolU, TolV, Curve, Surface)`
/// (`cxx:307-420`).
///
/// `Sol` is `(t, u, v)`; on success it is replaced by the boundary solution.
pub(crate) fn exact_bound(
    sol: &mut GpPnt,
    not_sol: f64,
    tol: f64,
    tol_u: f64,
    tol_v: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
) -> bool {
    let u0 = sol.y();
    let v0 = sol.z();
    let (first_u, last_u) = surface.u_range();
    let (first_v, last_v) = surface.v_range();
    let d2d = d1(sol.x(), u0, v0, curve, surface)
        .expect("ExactBound: d1 raised Standard_ConstructionError (cxx:158)");

    // Estimate which border has more opportunity to be intersected
    // (`cxx:326-349`).
    let (ru1, ru2, rv1, rv2);
    if d2d.x().abs() < RESOLUTION {
        ru1 = INFINITE;
        ru2 = INFINITE;
        rv1 = v0 - first_v;
        rv2 = last_v - v0;
    } else if d2d.y().abs() < RESOLUTION {
        ru1 = u0 - first_u;
        ru2 = last_u - u0;
        rv1 = INFINITE;
        rv2 = INFINITE;
    } else {
        ru1 = pnt2d_distance(u0, v0, first_u, v0 + (first_u - u0) * d2d.y() / d2d.x());
        ru2 = pnt2d_distance(u0, v0, last_u, v0 + (last_u - u0) * d2d.y() / d2d.x());
        rv1 = pnt2d_distance(u0, v0, u0 + (first_v - v0) * d2d.x() / d2d.y(), first_v);
        rv2 = pnt2d_distance(u0, v0, u0 + (last_v - v0) * d2d.x() / d2d.y(), last_v);
    }

    // `Seq` (`cxx:350-354`): `(border parameter, estimated distance, kind)`.
    let mut seq = [
        (first_u, ru1, 2i32),
        (last_u, ru2, 2),
        (first_v, rv1, 3),
        (last_v, rv2, 3),
    ];
    // Bubble sort on the estimated distance, largest first (`cxx:355-368`);
    // the loop below pops from the end, so the nearest border is tried first.
    for i in 0..3 {
        for j in 0..3 - i {
            if seq[j].1 < seq[j + 1].1 {
                seq.swap(j, j + 1);
            }
        }
    }

    let t = sol.x();
    let t1 = t.min(not_sol);
    let t2 = t.max(not_sol);

    for k in (0..4).rev() {
        let p_x = seq[k].0;
        if seq[k].2 == 2 {
            let mut prj = PrjResolve::new(curve, surface, 2).expect("myFix in 1..=3");
            prj.perform(
                t,
                p_x,
                v0,
                (tol, tol_v),
                (t1, first_v),
                (t2, last_v),
                FUNC_TOL,
                false,
            );
            if !prj.is_done() {
                continue;
            }
            let (su, sv) = prj.solution().expect("IsDone");
            *sol = GpPnt::new(su, p_x, sv);
            return true;
        } else {
            let mut prj = PrjResolve::new(curve, surface, 3).expect("myFix in 1..=3");
            prj.perform(
                t,
                u0,
                p_x,
                (tol, tol_u),
                (t1, first_u),
                (t2, last_u),
                FUNC_TOL,
                false,
            );
            if !prj.is_done() {
                continue;
            }
            let (su, sv) = prj.solution().expect("IsDone");
            *sol = GpPnt::new(su, sv, p_x);
            return true;
        }
    }
    false
}

/// `DichExactBound(Sol, NotSol, Tol, TolU, TolV, Curve, Surface)`
/// (`cxx:423-472`).
pub(crate) fn dich_exact_bound(
    sol: &mut GpPnt,
    not_sol: f64,
    tol: f64,
    tol_u: f64,
    tol_v: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
) {
    let mut u0 = sol.y();
    let mut v0 = sol.z();
    let mut prj = PrjResolve::new(curve, surface, 1).expect("myFix in 1..=3");
    let mut a_not_sol = not_sol;
    while (sol.x() - a_not_sol).abs() > tol {
        let t = (sol.x() + a_not_sol) / 2.0;
        prj.perform(
            t,
            u0,
            v0,
            (tol_u, tol_v),
            (surface.u_range().0, surface.v_range().0),
            (surface.u_range().1, surface.v_range().1),
            FUNC_TOL,
            true,
        );
        if prj.is_done() {
            let (su, sv) = prj.solution().expect("IsDone");
            *sol = GpPnt::new(t, su, sv);
            u0 = su;
            v0 = sv;
        } else {
            a_not_sol = t;
        }
    }
}

/// `InitialPoint(Point, t, C, S, TolU, TolV, U, V, MaxDist)` (`cxx:475-543`).
///
/// `Some((u, v))` is OCCT's `true` with `U` / `V` written back.
pub(crate) fn initial_point(
    point: &GpPnt,
    t: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
    tol_u: f64,
    tol_v: f64,
    max_dist: f64,
) -> Option<(f64, f64)> {
    let mut prj = PrjResolve::new(curve, surface, 1).expect("myFix in 1..=3");
    let mut ext = ExtPs::new();
    ext.initialize(
        surface,
        surface.u_range().0,
        surface.u_range().1,
        surface.v_range().0,
        surface.v_range().1,
        tol_u,
        tol_v,
    );
    ext.set_flag(ExtremaExtFlag::Min);
    ext.perform(point);

    let mut a_max_dist = max_dist;
    if a_max_dist > 0.0 {
        a_max_dist *= a_max_dist;
    }

    let mut argmin = 0usize;
    if ext.is_done() && ext.nb_ext() > 0 {
        let nend = ext.nb_ext();
        for i in 1..=nend {
            if a_max_dist > 0.0 && a_max_dist < ext.square_distance(i) {
                continue;
            }
            let (par_u, par_v, _) = ext.point(i);
            prj.perform(
                t,
                par_u,
                par_v,
                (tol_u, tol_v),
                (surface.u_range().0, surface.v_range().0),
                (surface.u_range().1, surface.v_range().1),
                FUNC_TOL,
                true,
            );
            if prj.is_done() && (argmin == 0 || ext.square_distance(i) < ext.square_distance(argmin)) {
                argmin = i;
            }
        }
    }
    if argmin == 0 {
        return None;
    }
    let (u, v, _) = ext.point(argmin);
    Some((u, v))
}

/// `gp_Pnt2d(U0, V0).Distance(gp_Pnt2d(U1, V1))`.
fn pnt2d_distance(u0: f64, v0: f64, u1: f64, v1: f64) -> f64 {
    let (du, dv) = (u1 - u0, v1 - v0);
    (du * du + dv * dv).sqrt()
}
