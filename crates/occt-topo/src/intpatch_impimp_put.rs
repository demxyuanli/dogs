//! `PutPointsOnLine` / `FindLine` / `GLine::ComputeVertexParameters`.
//! Source: `IntPatch_ImpImpIntersection.cxx` PutPointsOnLine, FindLine.

use occt_core::elib::clib;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::PCONFUSION;
use occt_geom::Surface;

use crate::geom_int::{surface_parameters, GLineKind, GeomIntLine};
use crate::int_tools_wline::PatchPoint;

use super::quad::ImplicitQuad;
use super::sonb::PathPoint;

/// `PutPointsOnLine`: snap restriction zeros onto existing GLines as vertices.
pub(crate) fn put_points_on_line(
    s1: &dyn Surface,
    s2: &dyn Surface,
    list: &[PathPoint],
    slin: &mut [GeomIntLine],
    on_first: bool,
    _q_dom: &ImplicitQuad,
    _q_other: &ImplicitQuad,
    tol: f64,
) {
    if slin.is_empty() || list.is_empty() {
        return;
    }
    let sq_tol = tol * tol;
    for pt in list {
        if let Some((idx, para)) = find_line(&pt.p, slin, sq_tol) {
            let Some(g) = slin[idx].as_gline_mut() else {
                continue;
            };
            let (u_dom, v_dom) = (pt.u, pt.v);
            let other_s = if on_first { s2 } else { s1 };
            let (u_oth, v_oth) = surface_parameters(other_s, &pt.p).unwrap_or((0.0, 0.0));
            let (u1, v1, u2, v2) = if on_first {
                (u_dom, v_dom, u_oth, v_oth)
            } else {
                (u_oth, v_oth, u_dom, v_dom)
            };
            let mut vp = PatchPoint::new(pt.p, para, u1, v1, u2, v2);
            if on_first {
                vp.on_dom_s1 = true;
            } else {
                vp.on_dom_s2 = true;
            }
            g.add_vertex(vp);
        }
    }
}

fn find_line(p: &GpPnt, slin: &[GeomIntLine], sq_tol: f64) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64, f64)> = None;
    for (i, line) in slin.iter().enumerate() {
        let Some(g) = line.as_gline() else {
            continue;
        };
        let para = parameter_on_gline(&g.kind, p);
        let q = value_on_gline(&g.kind, para);
        let sq = p.square_distance(&q);
        if sq < sq_tol {
            if best.map(|b| sq < b.2).unwrap_or(true) {
                best = Some((i, para, sq));
            }
        }
    }
    best.map(|(i, para, _)| (i, para))
}

/// Squared 3D distance from `p` to a GLine (`ProcessRLine` / `IsRLineGood`).
pub(crate) fn square_distance_to_gline(kind: &GLineKind, p: &GpPnt) -> f64 {
    let para = parameter_on_gline(kind, p);
    p.square_distance(&value_on_gline(kind, para))
}

fn parameter_on_gline(kind: &GLineKind, p: &GpPnt) -> f64 {
    match kind {
        GLineKind::Lin(l) => {
            let dir = GpVec::from_xyz(l.direction().xyz());
            GpVec::from_pnts(&l.location(), p).dot(&dir)
        }
        GLineKind::Circ(c) => {
            let pos = c.position();
            let w = GpVec::from_pnts(&c.location(), p);
            let x = w.dot(&GpVec::from_xyz(pos.x_direction().xyz()));
            let y = w.dot(&GpVec::from_xyz(pos.y_direction().xyz()));
            y.atan2(x)
        }
        GLineKind::Elips(e) => {
            let pos = e.position();
            let w = GpVec::from_pnts(&e.location(), p);
            let x = w.dot(&GpVec::from_xyz(pos.x_direction().xyz())) / e.major_radius;
            let y = w.dot(&GpVec::from_xyz(pos.y_direction().xyz())) / e.minor_radius;
            y.atan2(x)
        }
        GLineKind::Parab(pb) => {
            let pos = pb.position();
            let w = GpVec::from_pnts(&pb.location(), p);
            w.dot(&GpVec::from_xyz(pos.y_direction().xyz()))
        }
        GLineKind::Hypr(h) => {
            let pos = h.position();
            let w = GpVec::from_pnts(&h.location(), p);
            let y = w.dot(&GpVec::from_xyz(pos.y_direction().xyz()));
            let b = h.minor_radius;
            if b.abs() < 1e-14 {
                0.0
            } else {
                (y / b).asinh()
            }
        }
    }
}

fn value_on_gline(kind: &GLineKind, t: f64) -> GpPnt {
    match kind {
        GLineKind::Lin(l) => clib::line_value(l, t),
        GLineKind::Circ(c) => clib::circle_value(c, t),
        GLineKind::Elips(e) => clib::ellipse_value(e, t),
        GLineKind::Parab(p) => clib::parabola_value(p, t),
        GLineKind::Hypr(h) => clib::hyperbola_value(h, t),
    }
}

/// Sort / dedup GLine and RLine vertices (`ComputeVertexParameters` core).
pub(crate) fn compute_vertex_parameters(slin: &mut [GeomIntLine], _tol: f64) {
    let a_tol_pc = 1000.0 * PCONFUSION;
    for line in slin.iter_mut() {
        let vertices = match line {
            GeomIntLine::Geometric(g) => &mut g.vertices,
            GeomIntLine::Restriction(r) => &mut r.vertices,
            GeomIntLine::Analytic(a) => &mut a.vertices,
            GeomIntLine::Walking(_) => continue,
        };
        vertices.sort_by(|a, b| a.param_on_line.total_cmp(&b.param_on_line));
        let mut i = 1;
        while i < vertices.len() {
            if (vertices[i].param_on_line - vertices[i - 1].param_on_line).abs() <= a_tol_pc {
                vertices.remove(i);
            } else {
                i += 1;
            }
        }
    }
}
