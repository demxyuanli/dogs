//! `IntPatch_TheSOnBounds` / `IntPatch_ArcFunction` on UV-box restriction arcs.
//! Source: `IntStart_SearchOnBoundaries.gxx` BoundedArc, `IntPatch_ArcFunction.cxx`.

use occt_core::gp::GpPnt;

use occt_geom::Surface;

use crate::geom_int::{RestrictionArc, TopolTool};

use super::quad::{distance, ImplicitQuad};

/// `IntPatch_ThePathPointOfTheSOnBounds`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PathPoint {
    pub p: GpPnt,
    pub u: f64,
    pub v: f64,
    pub param_on_arc: f64,
    pub arc: RestrictionArc,
}

/// Result of `TheSOnBounds::Perform`.
pub(crate) struct BoundsSol {
    pub done: bool,
    pub all_arc_solution: bool,
    pub points: Vec<PathPoint>,
    pub segments: Vec<BoundSegment>,
}

/// `IntPatch_TheSegmentOfTheSOnBounds`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BoundSegment {
    pub arc: RestrictionArc,
    pub first: Option<PathPoint>,
    pub last: Option<PathPoint>,
}

/// `IntPatch_TheSOnBounds::Perform(F, Domain, TolBoundary, TolTangency)`.
pub(crate) fn search_on_bounds(
    surf: &dyn Surface,
    quad: &ImplicitQuad,
    domain: &TopolTool,
    tol_boundary: f64,
    tol_tangency: f64,
) -> BoundsSol {
    let arcs = domain.restriction_arcs();
    if arcs.is_empty() {
        return BoundsSol {
            done: true,
            all_arc_solution: false,
            points: Vec::new(),
            segments: Vec::new(),
        };
    }
    let mut all = true;
    let mut points = Vec::new();
    let mut segments = Vec::new();
    for arc in &arcs {
        let mut arc_sol = false;
        bounded_arc(
            surf,
            quad,
            arc,
            tol_boundary,
            tol_tangency,
            &mut points,
            &mut segments,
            &mut arc_sol,
        );
        all = all && arc_sol;
    }
    BoundsSol {
        done: true,
        all_arc_solution: all,
        points,
        segments,
    }
}

fn bounded_arc(
    surf: &dyn Surface,
    quad: &ImplicitQuad,
    arc: &RestrictionArc,
    tol_boundary: f64,
    tol_tangency: f64,
    points: &mut Vec<PathPoint>,
    segments: &mut Vec<BoundSegment>,
    arc_sol: &mut bool,
) {
    let pdeb = arc.first;
    let pfin = arc.last;
    if !pdeb.is_finite() || !pfin.is_finite() || (pfin - pdeb).abs() < 1e-14 {
        *arc_sol = false;
        return;
    }
    let mut n_tol = tol_tangency;
    if (pfin - pdeb) < (tol_tangency * 10.0) {
        n_tol = (pfin - pdeb) * 0.1;
    }
    let mut vals: Vec<(f64, f64, GpPnt, f64, f64)> = Vec::with_capacity(101);
    let mut max_abs: f64 = 0.0;
    let nb_echant = 100;
    for i in 0..=nb_echant {
        let t = pdeb + (pfin - pdeb) * i as f64 / nb_echant as f64;
        let uv = arc.value(t);
        let p = surf.d0(uv.x(), uv.y());
        let f = distance(quad, &p);
        max_abs = max_abs.max(f.abs());
        vals.push((t, f, p, uv.x(), uv.y()));
    }
    *arc_sol = max_abs <= tol_boundary;
    if *arc_sol {
        let (t0, _, p0, u0, v0) = vals[0];
        let (t1, _, p1, u1, v1) = vals[vals.len() - 1];
        let pf = PathPoint {
            p: p0,
            u: u0,
            v: v0,
            param_on_arc: t0,
            arc: *arc,
        };
        let pl = PathPoint {
            p: p1,
            u: u1,
            v: v1,
            param_on_arc: t1,
            arc: *arc,
        };
        push_point(points, p0, u0, v0, t0, n_tol, *arc);
        push_point(points, p1, u1, v1, t1, n_tol, *arc);
        segments.push(BoundSegment {
            arc: *arc,
            first: Some(pf),
            last: Some(pl),
        });
        return;
    }

    for i in 0..vals.len() {
        let (t, f, p, u, v) = vals[i];
        if f.abs() <= n_tol {
            push_point(points, p, u, v, t, n_tol, *arc);
        }
        if i + 1 < vals.len() {
            let (t2, f2, _, _, _) = vals[i + 1];
            if f * f2 < 0.0 {
                if let Some((tr, pr, ur, vr)) =
                    root_on_interval(surf, quad, arc, t, t2, f, f2, n_tol)
                {
                    push_point(points, pr, ur, vr, tr, n_tol, *arc);
                }
            }
        }
    }
}

fn push_point(
    points: &mut Vec<PathPoint>,
    p: GpPnt,
    u: f64,
    v: f64,
    t: f64,
    tol: f64,
    arc: RestrictionArc,
) {
    if points.iter().any(|q| q.p.square_distance(&p) <= tol * tol) {
        return;
    }
    points.push(PathPoint {
        p,
        u,
        v,
        param_on_arc: t,
        arc,
    });
}

fn root_on_interval(
    surf: &dyn Surface,
    quad: &ImplicitQuad,
    arc: &RestrictionArc,
    mut a: f64,
    mut b: f64,
    mut fa: f64,
    mut fb: f64,
    tol: f64,
) -> Option<(f64, GpPnt, f64, f64)> {
    for _ in 0..40 {
        let m = 0.5 * (a + b);
        let uv = arc.value(m);
        let p = surf.d0(uv.x(), uv.y());
        let fm = distance(quad, &p);
        if fm.abs() <= tol || (b - a).abs() <= 1e-10 {
            return Some((m, p, uv.x(), uv.y()));
        }
        if fa * fm <= 0.0 {
            b = m;
            fb = fm;
        } else {
            a = m;
            fa = fm;
        }
        let _ = fb;
    }
    let m = 0.5 * (a + b);
    let uv = arc.value(m);
    let p = surf.d0(uv.x(), uv.y());
    Some((m, p, uv.x(), uv.y()))
}
