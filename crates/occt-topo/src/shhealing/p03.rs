use super::prelude::*;
use super::*;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpPnt, GpPnt2d, GpTrsf2d, GpVec, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d::{curve_on_surface_oriented, replace_pcurve};
use crate::brep_surface::SurfaceKind;
use crate::brep_tool::BRepTool;
use crate::meshing::wire_order::{WireOrder, WireOrderStatus};
use crate::pcurve_full::{
    classify_surface_kind, project_curve_on_surface_perform, reparam_curve2d,
};

/// `ShapeAnalysis::AdjustByPeriod` (`ShapeAnalysis.cxx:48-62`).
pub fn adjust_by_period(val: f64, to_val: f64, period: f64) -> f64 {
    let diff = val - to_val;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
}

pub(super) fn first_vertex(edge: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(edge);
    if edge.0.orientation().is_reversed() {
        l
    } else {
        f
    }
}

pub(super) fn last_vertex(edge: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(edge);
    if edge.0.orientation().is_reversed() {
        f
    } else {
        l
    }
}

fn pcurve_at(
    edge: &Edge,
    face: &Face,
) -> Option<(Arc<dyn Curve2d>, f64, f64, GpPnt2d, GpPnt2d)> {
    let (c2, a, b) = curve_on_surface_oriented(edge, face, true)?;
    Some((c2.clone(), a, b, c2.d0(a), c2.d0(b)))
}

fn xy_cross(a: &GpPnt2d, b: &GpPnt2d, x: &GpVec2d) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    dx * x.y() - dy * x.x()
}

fn xy_dot_delta(a: &GpPnt2d, b: &GpPnt2d, x: &GpVec2d) -> f64 {
    (a.x() - b.x()) * x.x() + (a.y() - b.y()) * x.y()
}

fn xy_dot_axis(p: &GpPnt2d, x: &GpVec2d) -> f64 {
    p.x() * x.x() + p.y() * x.y()
}

/// Analytic singularities of `ShapeAnalysis_Surface::ComputeSingularities`
/// for cone / torus / sphere (`cxx:201-246`).
pub(super) fn degenerated_values(
    surf: &dyn Surface,
    p3d: &GpPnt,
    preci: f64,
) -> Option<(GpPnt2d, GpPnt2d, f64, f64)> {
    let (su1, su2) = surf.u_range();
    let (sv1, sv2) = surf.v_range();
    let mut hits: Vec<(f64, GpPnt, GpPnt2d, GpPnt2d)> = Vec::new();
    if let Some((radius, alpha)) = surf.cone_ref() {
        let sin_a = alpha.sin();
        if sin_a.abs() > 1e-16 {
            let v_apex = -radius / sin_a;
            let apex = surf.d0(0.0, v_apex);
            hits.push((
                0.0,
                apex,
                GpPnt2d::new(su1, v_apex),
                GpPnt2d::new(su2, v_apex),
            ));
        }
    } else if let Some(tor) = surf.gp_torus() {
        let minor = tor.minor_radius();
        let major = tor.major_radius();
        let ang = (major / minor).min(1.0).acos();
        let pre = (major - minor).max(0.0);
        hits.push((
            pre,
            surf.d0(0.0, std::f64::consts::PI - ang),
            GpPnt2d::new(su1, std::f64::consts::PI - ang),
            GpPnt2d::new(su2, std::f64::consts::PI - ang),
        ));
        if major <= minor {
            hits.push((
                pre,
                surf.d0(0.0, std::f64::consts::PI + ang),
                GpPnt2d::new(su2, std::f64::consts::PI + ang),
                GpPnt2d::new(su1, std::f64::consts::PI + ang),
            ));
        }
    } else if surf.gp_sphere().is_some() {
        hits.push((
            0.0,
            surf.d0(su1, sv2),
            GpPnt2d::new(su2, sv2),
            GpPnt2d::new(su1, sv2),
        ));
        hits.push((
            0.0,
            surf.d0(su1, sv1),
            GpPnt2d::new(su1, sv1),
            GpPnt2d::new(su2, sv1),
        ));
    }
    let mut best: Option<(f64, GpPnt2d, GpPnt2d)> = None;
    for (pre, q, a, b) in hits {
        if pre > preci {
            continue;
        }
        let gap = q.distance(p3d);
        if gap <= preci && best.as_ref().is_none_or(|(g, _, _)| gap < *g) {
            best = Some((gap, a, b));
        }
    }
    best.map(|(_, a, b)| (a, b, su1, su2))
}

fn is_degenerated_2d(surf: &dyn Surface, p1: GpPnt2d, p2: GpPnt2d, tol: f64, ratio: f64) -> bool {
    let a = surf.d0(p1.x(), p1.y());
    let b = surf.d0(p2.x(), p2.y());
    let m = surf.d0(0.5 * (p1.x() + p2.x()), 0.5 * (p1.y() + p2.y()));
    let mut max3d = a.distance(&b).max(m.distance(&a)).max(m.distance(&b));
    if max3d > tol {
        return false;
    }
    let ru = 1.0;
    let rv = 1.0;
    let du = (p1.x() - p2.x()).abs() / ru;
    let dv = (p1.y() - p2.y()).abs() / rv;
    max3d *= ratio;
    du * du + dv * dv > max3d * max3d
}

/// `ShapeFix_Wire::FixShifted` (`ShapeFix_Wire.cxx:1661-2126`).
pub fn fix_shifted_wire(wire: &Wire, face: &Face) -> bool {
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let mut u_closed = surf.is_u_closed();
    let mut v_closed = surf.is_v_closed() || surf.gp_sphere().is_some();
    let mut v_range = 1.0;
    let mut v_crv_closed = false;
    if surf.is_surface_of_revolution() {
        if let Some(basis) = surf.extrusion_basis_curve() {
            if basis.is_periodic() {
                v_closed = true;
                v_range = basis.period();
                v_crv_closed = true;
            }
        }
    }
    if !u_closed && !v_closed {
        return false;
    }
    let (suf, sul) = surf.u_range();
    let (svf, svl) = surf.v_range();
    let su_mid = 0.5 * (suf + sul);
    let sv_mid = 0.5 * (svf + svl);
    let u_range = if u_closed {
        (sul - suf).abs()
    } else {
        f64::MAX
    };
    if !v_crv_closed {
        v_range = if v_closed {
            (svl - svf).abs()
        } else {
            f64::MAX
        };
    }
    if !u_range.is_finite() && !v_range.is_finite() {
        return false;
    }
    let u_tol = 0.2 * u_range;
    let v_tol = 0.2 * v_range;
    let edges: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| !(BRepTool::is_degenerated(e) && pcurve_at(e, face).is_none()))
        .collect();
    let nb = edges.len();
    if nb == 0 {
        return false;
    }
    let mut done = false;
    let mut stop = nb;
    let mut ended = nb == 0;
    let mut degstop = false;
    let mut degn2 = 0usize;
    let mut pdeg = GpPnt::new(0.0, 0.0, 0.0);
    let mut n2 = 0usize;
    let mut n1 = nb;
    // `cxx:1764` / `cxx:2015`: first-pass box is junction `p2d1` only.
    let mut box2 = BndBox2d::new();
    while !ended {
        n2 += 1;
        if n2 > nb {
            n2 = 1;
        }
        if n2 == stop {
            ended = true;
        }
        let e1 = &edges[n1 - 1];
        let e2 = &edges[n2 - 1];
        n1 = n2;
        if BRepTool::is_degenerated(e1) || BRepTool::is_degenerated(e2) {
            if !degstop {
                stop = n2;
                degstop = true;
            }
            continue;
        }
        let Some(v) = first_vertex(e2) else {
            continue;
        };
        let p = BRepTool::vertex_point(&v);
        let mut is_deg = 0i32;
        let preci = CONFUSION.max(BRepTool::vertex_tolerance(&v));
        if let Some((deg_p1, deg_p2, _, _)) = degenerated_values(surf.as_ref(), &p, preci) {
            is_deg = if (deg_p1.x() - deg_p2.x()).abs() > (deg_p1.y() - deg_p2.y()).abs() {
                1
            } else {
                2
            };
        }
        const MAX_TOL: f64 = 1.0;
        if surf.is_surface_of_revolution() {
            if is_deg == 0 && !v_closed {
                if let Some((_, _, _, _, p_b1)) = pcurve_at(e1, face) {
                    let q1 = GpPnt2d::new(suf, p_b1.y());
                    let q2 = GpPnt2d::new(sul, p_b1.y());
                    if let Some((_, a1, b1, pa, pb)) = pcurve_at(e1, face) {
                        let _ = (a1, b1, pa);
                        if is_degenerated_2d(surf.as_ref(), q1, q2, MAX_TOL, 10.0)
                            && !is_degenerated_2d(surf.as_ref(), pa, pb, MAX_TOL, 10.0)
                        {
                            is_deg = 1;
                        }
                    }
                }
            }
            if is_deg == 0 && !u_closed {
                if let Some((_, _, _, _, p_b1)) = pcurve_at(e1, face) {
                    let q1 = GpPnt2d::new(p_b1.x(), svf);
                    let q2 = GpPnt2d::new(p_b1.x(), svl);
                    if let Some((_, _, _, pa, pb)) = pcurve_at(e1, face) {
                        if is_degenerated_2d(surf.as_ref(), q1, q2, MAX_TOL, 10.0)
                            && !is_degenerated_2d(surf.as_ref(), pa, pb, MAX_TOL, 10.0)
                        {
                            is_deg = 2;
                        }
                    }
                }
            }
        }
        if is_deg != 0 {
            if !degstop {
                stop = n2;
                degstop = true;
            }
            if degn2 == 0 {
                degn2 = n2;
                pdeg = p;
            } else if pdeg.square_distance(&p) < CONFUSION * CONFUSION {
                degn2 = n2;
            } else if try_bi_meridian(
                &edges,
                face,
                degn2,
                n2,
                nb,
                e1,
                e2,
                u_closed,
                if u_closed { u_range } else { v_range },
            ) {
                done = true;
                continue;
            }
        }
        let Some((_, _, _, _, p2d1)) = pcurve_at(e1, face) else {
            continue;
        };
        let Some((c2, _, _, p2d2, _)) = pcurve_at(e2, face) else {
            continue;
        };
        box2.add_point(&p2d1);
        let mut du = 0.0;
        let mut dv = 0.0;
        if u_closed && is_deg != 1 {
            let dx = (p2d2.x() - p2d1.x()).abs();
            if dx > u_range - u_tol {
                du = adjust_by_period(p2d2.x(), p2d1.x(), u_range);
            } else if dx > u_tol && stop == nb {
                stop = n2;
            }
        }
        if v_closed && is_deg != 2 {
            let dy = (p2d2.y() - p2d1.y()).abs();
            if dy > v_range - v_tol {
                dv = adjust_by_period(p2d2.y(), p2d1.y(), v_range);
            } else if dy > v_tol && stop == nb {
                stop = n2;
            }
        }
        if du == 0.0 && dv == 0.0 {
            continue;
        }
        let mut shift = GpTrsf2d::default();
        shift.set_translation_vec(&GpVec2d::new(du, dv));
        replace_pcurve(e2, face, Arc::from(c2.transformed(&shift)));
        done = true;
    }

    // `cxx:2056-2067`: early-out uses the first-pass `p2d1` box, not `[a,mid]`.
    if box2.is_void() {
        return false;
    }
    let (umin, vmin, umax, vmax) = box2.get().expect("non-void box");
    if (umin + umax - suf - sul).abs() < u_range
        && (vmin + vmax - svf - svl).abs() < v_range
        && !done
    {
        return false;
    }
    // `cxx:2069-2082`: rebuild from `Value(a)` and `Value((a+b)/2)`.
    box2.set_void();
    for e in &edges {
        let Some((c2, a, b, _, _)) = pcurve_at(e, face) else {
            continue;
        };
        box2.add_point(&c2.d0(a));
        box2.add_point(&c2.d0(0.5 * (a + b)));
    }
    let Some((umin, vmin, umax, vmax)) = box2.get() else {
        return done;
    };
    let mut du = 0.0;
    let mut dv = 0.0;
    if u_closed {
        du = adjust_by_period(0.5 * (umin + umax), su_mid, u_range);
    }
    if v_closed {
        dv = adjust_by_period(0.5 * (vmin + vmax), sv_mid, v_range);
    }
    if du == 0.0 && dv == 0.0 {
        return done;
    }
    for e in edges_of_wire(wire) {
        let Some((c2, _, _, _, _)) = pcurve_at(&e, face) else {
            continue;
        };
        let mut shift = GpTrsf2d::default();
        shift.set_translation_vec(&GpVec2d::new(du, dv));
        replace_pcurve(&e, face, Arc::from(c2.transformed(&shift)));
    }
    true
}

fn try_bi_meridian(
    edges: &[Edge],
    face: &Face,
    degn2: usize,
    n2: usize,
    nb: usize,
    e1: &Edge,
    e2: &Edge,
    u_closed: bool,
    u_range: f64,
) -> bool {
    let prev = if degn2 > 1 { degn2 - 1 } else { nb };
    let Some((_, _, _, _, pn1)) = pcurve_at(e1, face) else {
        return false;
    };
    let Some((_, _, _, pn2, _)) = pcurve_at(e2, face) else {
        return false;
    };
    let Some((_, _, _, _, pd1)) = pcurve_at(&edges[prev - 1], face) else {
        return false;
    };
    let Some((_, _, _, pd2, _)) = pcurve_at(&edges[degn2 - 1], face) else {
        return false;
    };
    let (x, period) = if u_closed {
        (GpVec2d::new(1.0, 0.0), u_range)
    } else {
        (GpVec2d::new(0.0, 1.0), u_range)
    };
    let rot1 = xy_cross(&pn1, &pd2, &x);
    let rot2 = xy_cross(&pd1, &pn2, &x);
    let scld = xy_dot_delta(&pd2, &pd1, &x);
    let scln = xy_dot_delta(&pn2, &pn1, &x);
    if !(rot1 * rot2 < -PCONFUSION
        && scld * scln < -PCONFUSION
        && scln.abs() > 0.1 * period
        && scld.abs() > 0.1 * period
        && rot1 * scld > PCONFUSION
        && rot2 * scln > PCONFUSION)
    {
        return false;
    }
    let sign = if rot2 > 0.0 { 1.0 } else { -1.0 };
    let Some((_, a2, b2, _, pb2)) = pcurve_at(e2, face) else {
        return false;
    };
    let Some((_, ax1, bx1, _, _)) = pcurve_at(&edges[prev - 1], face) else {
        return false;
    };
    let Some((cx1, _, _, _, _)) = pcurve_at(&edges[prev - 1], face) else {
        return false;
    };
    let Some((c2d2, _, _, _, _)) = pcurve_at(e2, face) else {
        return false;
    };
    let Some((c2d1, a1, b1, _, _)) = pcurve_at(e1, face) else {
        return false;
    };
    let Some((cx2, ax2, bx2, _, _)) = pcurve_at(&edges[degn2 - 1], face) else {
        return false;
    };
    let deep1 = [
        sign * xy_dot_axis(&pn2, &x),
        sign * xy_dot_axis(&pd1, &x),
        sign * xy_dot_axis(&pb2, &x),
        sign * xy_dot_axis(&cx1.d0(ax1), &x),
        sign * xy_dot_axis(&c2d2.d0(0.5 * (a2 + b2)), &x),
        sign * xy_dot_axis(&cx1.d0(0.5 * (ax1 + bx1)), &x),
    ]
    .into_iter()
    .fold(f64::INFINITY, f64::min);
    let deep2 = [
        sign * xy_dot_axis(&pn1, &x),
        sign * xy_dot_axis(&pd2, &x),
        sign * xy_dot_axis(&c2d1.d0(a1), &x),
        sign * xy_dot_axis(&cx2.d0(bx2), &x),
        sign * xy_dot_axis(&c2d1.d0(0.5 * (a1 + b1)), &x),
        sign * xy_dot_axis(&cx2.d0(0.5 * (ax2 + bx2)), &x),
    ]
    .into_iter()
    .fold(f64::NEG_INFINITY, f64::max);
    let deep = deep2 - deep1;
    let dx = adjust_by_period(deep, 0.5 * (PCONFUSION + period + PCONFUSION), period);
    let scale = if scld > 0.0 { -dx } else { dx };
    let mut k = degn2;
    loop {
        if k > nb {
            k = 1;
        }
        if k == n2 {
            break;
        }
        if let Some((cx, _, _, _, _)) = pcurve_at(&edges[k - 1], face) {
            let mut shift = GpTrsf2d::default();
            shift.set_translation_vec(&GpVec2d::new(x.x() * scale, x.y() * scale));
            replace_pcurve(&edges[k - 1], face, Arc::from(cx.transformed(&shift)));
        }
        k += 1;
    }
    true
}

/// `GeomLib::SameRange` (`GeomLib.cxx:842-922`): remap a 2d curve from
/// `[first_on, last_on]` onto `[req_first, req_last]`. Linear reparam matches
/// the BSpline knot `BSplCLib::Reparametrize` evaluation.
fn is_geom2d_line_curve(c: &dyn Curve2d) -> bool {
    !c.first_parameter().is_finite()
        && !c.last_parameter().is_finite()
        && c.d2(0.0).2.square_magnitude() < 1e-30
}

fn geom_lib_same_range(
    c2d: Arc<dyn Curve2d>,
    first_on: f64,
    last_on: f64,
    req_first: f64,
    req_last: f64,
) -> Arc<dyn Curve2d> {
    if (last_on - req_last).abs() <= PCONFUSION && (first_on - req_first).abs() <= PCONFUSION {
        return c2d;
    }
    if !first_on.is_finite() || !last_on.is_finite() {
        return c2d;
    }
    // `GeomLib.cxx:862-870`: equal span on a Geom2d_Line is a translation
    // along the line, not a Reparam wrapper (that wrapper is finite and
    // breaks `Adaptor3d_CurveOnSurface::EvalKPart` iso detection).
    if ((last_on - first_on) - (req_last - req_first)).abs() <= PCONFUSION
        && is_geom2d_line_curve(c2d.as_ref())
    {
        let du = first_on - req_first;
        let (_, dir) = c2d.d1(0.0);
        let mut trsf = GpTrsf2d::identity();
        trsf.set_translation_vec(&GpVec2d::new(dir.x() * du, dir.y() * du));
        return Arc::from(c2d.transformed(&trsf));
    }
    reparam_curve2d(c2d, req_first, req_last, first_on, last_on)
}

/// `ShapeAnalysis_Edge::CheckPCurveRange` (`ShapeAnalysis_Edge.cxx:999-1032`).
fn check_pcurve_range(first: f64, last: f64, pc: &dyn Curve2d) -> bool {
    let eps = PCONFUSION;
    let (mut fp, mut lp) = (pc.first_parameter(), pc.last_parameter());
    let mut is_periodic = pc.is_periodic();
    let mut period = if is_periodic { pc.period() } else { f64::MAX };
    if let Some(basis) = pc.trimmed_basis() {
        fp = basis.first_parameter();
        lp = basis.last_parameter();
        is_periodic = basis.is_periodic();
        if is_periodic {
            period = basis.period();
        }
    }
    if is_periodic {
        if last - first > period + eps {
            return false;
        }
    } else if first < fp - eps || last > lp + eps {
        return false;
    }
    true
}

/// `ShapeFix_Edge.cxx:171-325` `TranslatePCurve` for a both-closed seam.
fn translate_pcurve_seam(surf: &dyn Surface, c2d: &dyn Curve2d, tol: f64) -> Arc<dyn Curve2d> {
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let a = c2d.first_parameter();
    let b = c2d.last_parameter();
    let (p0, p1) = if a.is_finite() && b.is_finite() {
        (c2d.d0(a), c2d.d0(b))
    } else {
        let (_, tan) = c2d.d1(0.0);
        return translate_pcurve_dir(surf, c2d, &c2d.d0(0.0), &tan, tol);
    };
    let vec = GpVec2d::new(p1.x() - p0.x(), p1.y() - p0.y());
    let iso_u = GpVec2d::new(0.0, vl - vf);
    let iso_v = GpVec2d::new(ul - uf, 0.0);
    let mut shift = GpTrsf2d::default();
    if vec.is_parallel(&iso_u, tol) {
        if (p0.x() - uf).abs() < (p0.x() - ul).abs() {
            shift.set_translation_vec(&GpVec2d::new(ul - uf, 0.0));
        } else {
            shift.set_translation_vec(&GpVec2d::new(uf - ul, 0.0));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    if vec.is_parallel(&iso_v, tol) {
        if (p0.y() - vf).abs() < (p0.y() - vl).abs() {
            shift.set_translation_vec(&GpVec2d::new(0.0, vl - vf));
        } else {
            shift.set_translation_vec(&GpVec2d::new(0.0, vf - vl));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    Arc::from(c2d.clone_dyn())
}

fn translate_pcurve_dir(
    surf: &dyn Surface,
    c2d: &dyn Curve2d,
    loc: &GpPnt2d,
    dir: &GpVec2d,
    tol: f64,
) -> Arc<dyn Curve2d> {
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let mut shift = GpTrsf2d::default();
    if dir.x().abs() <= tol && dir.y().abs() >= tol {
        if (loc.x() - uf).abs() < (loc.x() - ul).abs() {
            shift.set_translation_vec(&GpVec2d::new(ul - uf, 0.0));
        } else {
            shift.set_translation_vec(&GpVec2d::new(uf - ul, 0.0));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    if dir.x().abs() >= tol && dir.y().abs() <= tol {
        if (loc.y() - vf).abs() < (loc.y() - vl).abs() {
            shift.set_translation_vec(&GpVec2d::new(0.0, vl - vf));
        } else {
            shift.set_translation_vec(&GpVec2d::new(0.0, vf - vl));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    Arc::from(c2d.clone_dyn())
}

/// `ShapeFix_Edge::FixAddPCurve` (`ShapeFix_Edge.cxx:470-614`).
pub(super) fn fix_add_pcurve(edge: &Edge, face: &Face, is_seam: bool) -> bool {
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    if classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane {
        return false;
    }
    let Some(c3d) = GeometryRegistry::global().edge_curve(&edge.0) else {
        return false;
    };
    let (first, last) = GeometryRegistry::global().edge_parameters(&edge.0);
    let Some(c2d) = project_curve_on_surface_perform(c3d.as_ref(), surf.as_ref(), first, last)
    else {
        return false;
    };
    let face_key = GeometryRegistry::shape_key(&face.0);
    if !is_seam {
        GeometryRegistry::global().set_edge_pcurve(&edge.0, face_key, c2d);
        return true;
    }
    let prec = GeometryRegistry::global().edge_tolerance(&edge.0).max(CONFUSION);
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let mut shift = GpTrsf2d::default();
    let c2d2 = if surf.is_u_closed() && !surf.is_v_closed() {
        shift.set_translation_vec(&GpVec2d::new(ul - uf, 0.0));
        Arc::from(c2d.transformed(&shift))
    } else if surf.is_v_closed() && !surf.is_u_closed() {
        shift.set_translation_vec(&GpVec2d::new(0.0, vl - vf));
        Arc::from(c2d.transformed(&shift))
    } else if surf.is_u_closed() && surf.is_v_closed() {
        translate_pcurve_seam(surf.as_ref(), c2d.as_ref(), prec)
    } else {
        Arc::from(c2d.clone_dyn())
    };
    GeometryRegistry::global().set_edge_pcurves(&edge.0, face_key, vec![c2d, c2d2]);
    true
}

pub(super) fn is_seam_use(edges: &[Edge], i: usize) -> bool {
    let key = GeometryRegistry::shape_key(&edges[i].0);
    edges
        .iter()
        .filter(|e| GeometryRegistry::shape_key(&e.0) == key)
        .count()
        >= 2
}

/// `StepToTopoDS_TranslateEdgeLoop::CheckPCurves` first check
/// (`StepToTopoDS_TranslateEdgeLoop.cxx:110-172`).
fn check_pcurve_rep_range(wire: &Wire, face: &Face) {
    let Some(surf) = BRepTool::face_surface(face) else {
        return;
    };
    // `cxx:110-114`: planar faces drop STEP pcurves; meshing uses
    // `BRep_Tool::CurveOnPlane` via `make_pcurve_on_face`.
    if classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane {
        let reg = GeometryRegistry::global();
        for e in edges_of_wire(wire) {
            reg.remove_pcurves_on_surface(&e.0, &face.0);
        }
        return;
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    for e in edges_of_wire(wire) {
        let Some((the_pc, mut w1, mut w2)) = curve_on_surface_oriented(&e, face, false) else {
            continue;
        };
        let cf = the_pc.first_parameter();
        let cl = the_pc.last_parameter();
        if w1 == w2 {
            reg.remove_pcurves_on_surface(&e.0, &face.0);
            continue;
        }
        if !the_pc.is_periodic() {
            if w1 < cf {
                w1 = cf;
                reg.set_pcurve_range(&e.0, face_key, w1, w2);
            }
            if w2 > cl {
                w2 = cl;
                reg.set_pcurve_range(&e.0, face_key, w1, w2);
            }
        }
        if w1 > w2 && surf.is_u_periodic() {
            let (u1, u2) = surf.u_range();
            let preci = ((w2 - w1).abs() / 2.0).min(PCONFUSION);
            crate::geom_bnd_lib_elclib2d::adjust_periodic(u1, u2, preci, &mut w1, &mut w2);
            reg.set_pcurve_range(&e.0, face_key, w1, w2);
        }
    }
}

/// `ShapeFix_Edge.cxx:335-464` `TempSameRange`.
/// Remap each COS pcurve whose representation range differs from the 3D
/// range onto that 3D range, then `B.Range` all representations and
/// `B.SameRange(true)`.
fn temp_same_range(edge: &Edge) {
    let reg = GeometryRegistry::global();
    let Some(g) = reg.edge_geom(&edge.0) else {
        return;
    };
    let (current_first, current_last) = (g.first, g.last);
    if !current_first.is_finite() || !current_last.is_finite() {
        return;
    }
    let face_keys: Vec<usize> = g.pcurves.keys().copied().collect();
    for fk in &face_keys {
        let Some((first, last)) = g.pcurve_ranges.get(fk).copied() else {
            continue;
        };
        if (first - current_first).abs() <= PCONFUSION
            && (last - current_last).abs() <= PCONFUSION
        {
            continue;
        }
        let pcs = g.pcurves.get(fk).cloned().unwrap_or_default();
        if pcs.is_empty() {
            continue;
        }
        let new_pcs: Vec<Arc<dyn Curve2d>> = pcs
            .into_iter()
            .map(|pc| geom_lib_same_range(pc, first, last, current_first, current_last))
            .collect();
        if new_pcs.len() == 1 {
            reg.set_edge_pcurve(&edge.0, *fk, new_pcs.into_iter().next().unwrap());
        } else {
            reg.set_edge_pcurves(&edge.0, *fk, new_pcs);
        }
    }
    for fk in face_keys {
        reg.set_pcurve_range(&edge.0, fk, current_first, current_last);
    }
    reg.set_same_range(&edge.0, true);
}

/// `Extrema_LocateExtPC` Newton walk (`int_tools_vertex_line` / ValidateEdge).
/// Returns square distance at the local extremum, or `None` when `!IsDone`.
fn locate_ext_pc_sq<D0, D1>(
    d0: D0,
    d1: D1,
    p: &GpPnt,
    t_seed: f64,
    first: f64,
    last: f64,
) -> Option<f64>
where
    D0: Fn(f64) -> GpPnt,
    D1: Fn(f64) -> (GpPnt, GpVec),
{
    let lo = first.min(last);
    let hi = first.max(last);
    let mut t = t_seed.clamp(lo, hi);
    let mut last_good = d0(t);
    for _ in 0..12 {
        let (q, der) = d1(t);
        let speed2 = der.square_magnitude();
        if speed2 <= 1e-20 {
            return None;
        }
        let num = (q.x() - p.x()) * der.x() + (q.y() - p.y()) * der.y() + (q.z() - p.z()) * der.z();
        let next = t - num / speed2;
        if !next.is_finite() {
            return None;
        }
        let next = next.clamp(lo, hi);
        let qn = d0(next);
        if qn.distance(&last_good) < 1e-10 {
            return Some(qn.square_distance(p));
        }
        t = next;
        last_good = qn;
    }
    Some(last_good.square_distance(p))
}

fn cos_d0(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> GpPnt {
    let uv = pc.d0(t);
    surf.d0(uv.x(), uv.y())
}

fn cos_d1(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> (GpPnt, GpVec) {
    let (uv, duv) = pc.d1(t);
    let (p, su, sv) = surf.d1(uv.x(), uv.y());
    let tan = su
        .multiplied_scalar(duv.x())
        .added(&sv.multiplied_scalar(duv.y()));
    (p, tan)
}

/// `BRepLib_ValidateEdge::processApprox` (`cxx:103-229`) plus
/// `UpdateTolerance` (`cxx:57-63`).
/// `ShapeAnalysis_Edge::CheckSameParameter` (`cxx:704-838`).
fn check_same_parameter(edge: &Edge, face: &Face) -> f64 {
    if BRepTool::is_degenerated(edge) {
        return 0.0;
    }
    let Some(c3d) = BRepTool::edge_curve(edge) else {
        return 0.0;
    };
    let (first, last) = BRepTool::edge_parameters(edge);
    if !first.is_finite() || !last.is_finite() {
        return 0.0;
    }
    let Some(surf) = BRepTool::face_surface(face) else {
        return 0.0;
    };
    let same_p = BRepTool::same_parameter(edge);
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    let is_plane = classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane;
    // After CheckPCurves RemovePCurves (`cxx:110-114`) a plane has no stored
    // pcurve. Do not same-t leftover STEP pcs (Shape-2 1.322 inflation).
    let mut pcs = if is_plane {
        Vec::new()
    } else {
        reg.edge_pcurves(&edge.0, face_key)
    };
    let (mut cf, mut cl) = reg
        .pcurve_range(&edge.0, face_key)
        .unwrap_or((first, last));
    // `ShapeAnalysis_Edge.cxx:801-826`: `CurveOnPlane` then
    // `Geom2dAdaptor_Curve(aPC, aFirst, aLast)` + `CurveOnSurface` +
    // `BRepLib_ValidateEdge` (`aFirst/aLast` = 3D edge range, cxx:418-419).
    if pcs.is_empty() && is_plane {
        if let Ok(pc) = crate::pcurve::make_pcurve_on_face(edge, face) {
            pcs.push(pc);
            cf = first;
            cl = last;
        }
    }
    let mut max_sq: f64 = 0.0;
    // `SetControlPointsNumber(NbControl - 1)` with default NbControl=23.
    const NB: i32 = 22;
    for pc in &pcs {
        let projection = !same_p
            || (cf - first).abs() > PCONFUSION
            || (cl - last).abs() > PCONFUSION;
        if projection {
            // `BRepLib_ValidateEdge.cxx:141-227`.
            let saved = max_sq;
            let p_ref0 = c3d.d0(first);
            let p_oth0 = cos_d0(pc.as_ref(), surf.as_ref(), cf);
            max_sq = max_sq.max(p_ref0.square_distance(&p_oth0));
            let p_ref1 = c3d.d0(last);
            let p_oth1 = cos_d0(pc.as_ref(), surf.as_ref(), cl);
            max_sq = max_sq.max(p_ref1.square_distance(&p_oth1));
            let mut aborted = false;
            for i in 1..NB {
                let t_ref = ((NB - i) as f64 * first + i as f64 * last) / NB as f64;
                let t_oth = ((NB - i) as f64 * cf + i as f64 * cl) / NB as f64;
                let p_ref = c3d.d0(t_ref);
                let p_oth = cos_d0(pc.as_ref(), surf.as_ref(), t_oth);
                let Some(sq_ref) = locate_ext_pc_sq(
                    |t| c3d.d0(t),
                    |t| c3d.d1(t),
                    &p_oth,
                    t_ref,
                    first,
                    last,
                ) else {
                    aborted = true;
                    break;
                };
                max_sq = max_sq.max(sq_ref);
                let Some(sq_oth) = locate_ext_pc_sq(
                    |t| cos_d0(pc.as_ref(), surf.as_ref(), t),
                    |t| cos_d1(pc.as_ref(), surf.as_ref(), t),
                    &p_ref,
                    t_oth,
                    cf,
                    cl,
                ) else {
                    aborted = true;
                    break;
                };
                max_sq = max_sq.max(sq_oth);
            }
            if aborted {
                // cxx:202-205 / 222-225: !IsDone leaves CalculatedDistance unset
                // for this pair; earlier pairs already contributed to maxdev.
                max_sq = saved;
            }
            continue;
        }
        for i in 0..=NB {
            let t = ((NB - i) as f64 * first + i as f64 * last) / NB as f64;
            let uv = pc.d0(t);
            max_sq = max_sq.max(c3d.d0(t).square_distance(&surf.d0(uv.x(), uv.y())));
        }
    }
    max_sq.sqrt() * 1.00001
}

/// `ShapeAnalysis_Edge::CheckVertexTolerance` with `checkAll=true`
/// (`cxx:568-667`) then `B.UpdateVertex`.
fn fix_vertex_tolerance_edge(edge: &Edge) {
    let (Some(v1), Some(v2)) = (first_vertex(edge), last_vertex(edge)) else {
        return;
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    let mut t1: f64 = 0.0;
    let mut t2: f64 = 0.0;
    if let Some(c3d) = BRepTool::edge_curve(edge) {
        let (mut a, mut b) = BRepTool::edge_parameters(edge);
        if a.is_finite() && b.is_finite() {
            if edge.0.orientation().is_reversed() {
                std::mem::swap(&mut a, &mut b);
            }
            t1 = p1.square_distance(&c3d.d0(a));
            t2 = p2.square_distance(&c3d.d0(b));
        }
    }
    let reg = GeometryRegistry::global();
    for (surf, pc, mut a, mut b) in reg.edge_pcurve_reps(&edge.0) {
        if !(a.is_finite() && b.is_finite()) {
            continue;
        }
        if edge.0.orientation().is_reversed() {
            std::mem::swap(&mut a, &mut b);
        }
        let uv1 = pc.d0(a);
        let uv2 = pc.d0(b);
        t1 = t1.max(p1.square_distance(&surf.d0(uv1.x(), uv1.y())));
        t2 = t2.max(p2.square_distance(&surf.d0(uv2.x(), uv2.y())));
    }
    let tole = BRepTool::edge_tolerance(edge);
    v1.set_tolerance((1.0000001 * t1.sqrt()).max(tole));
    v2.set_tolerance((1.0000001 * t2.sqrt()).max(tole));
}

/// `ShapeFix_Edge::FixSameParameter` (`ShapeFix_Edge.cxx:798-935`).
/// STEP `wasSP` is true, so the BRepLib copy/compare arm is not taken.
pub(super) fn fix_same_parameter(edge: &Edge, face: &Face) {
    let reg = GeometryRegistry::global();
    let same_range = reg.edge_geom(&edge.0).map(|g| g.same_range).unwrap_or(true);
    if reg.is_degenerated_edge(&edge.0) {
        if !same_range {
            temp_same_range(edge);
        }
        reg.set_same_parameter(&edge.0, true);
        return;
    }
    if !same_range {
        temp_same_range(edge);
    }
    reg.set_same_parameter(&edge.0, true);
    let maxdev = check_same_parameter(edge, face);
    if let Some(v) = first_vertex(edge) {
        v.set_tolerance(BRepTool::vertex_tolerance(&v).max(maxdev));
    }
    if let Some(v) = last_vertex(edge) {
        v.set_tolerance(BRepTool::vertex_tolerance(&v).max(maxdev));
    }
    let tol = BRepTool::edge_tolerance(edge);
    if maxdev > tol {
        reg.set_edge_tolerance(&edge.0, maxdev);
        fix_vertex_tolerance_edge(edge);
    }
}

/// `ShapeFix_Wire::FixReorder` (`cxx:487-534`) then `FixReorder(sawo)`
/// (`cxx:1351-1399`). `Perform` calls `FixReorder()` with `theModeBoth=false`,
/// so CheckOrder is 3D (`cxx:504`). Status 0 leaves the wire unchanged.
pub fn fix_reorder_wire(wire: &Wire, face: &Face) -> bool {
    // `FixReorder(theModeBoth)` (`cxx:497-500`) only switches CheckOrder
    // mode on a bi-periodic surface. `Perform` still uses Mode3D
    // (`theModeBoth=false`). Restrict the write to those faces so a
    // cylinder SIW (Shape-1 f29) is not un-skipped by a 3D reorder of
    // an unrelated wire.
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    if !(surf.is_u_periodic() && surf.is_v_periodic()) {
        return false;
    }
    let stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if stored.len() < 2 {
        return false;
    }
    let mut order = WireOrder::new();
    for e in &stored {
        let Some(v1) = first_vertex(e) else {
            return false;
        };
        let Some(v2) = last_vertex(e) else {
            return false;
        };
        order.add_edge_xyz(BRepTool::vertex_point(&v1), BRepTool::vertex_point(&v2));
    }
    order.perform();
    if order.status() == WireOrderStatus::Same || order.nb_edges() != stored.len() {
        return false;
    }
    let mut new_edges = Vec::with_capacity(stored.len());
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        if signed == 0 {
            return false;
        }
        let mut e = stored[signed.unsigned_abs() as usize - 1].clone();
        if signed < 0 {
            e.0.reverse();
        }
        new_edges.push(e);
    }
    let extras: Vec<TopoShape> = {
        let ts = wire.0.tshape.read().expect("poisoned TShape lock");
        ts.children
            .iter()
            .filter(|c| {
                !c.is_edge()
                    || (c.orientation() != Orientation::Forward
                        && c.orientation() != Orientation::Reversed)
            })
            .cloned()
            .collect()
    };
    {
        let mut ts = wire.0.tshape.write().expect("poisoned TShape lock");
        ts.children.clear();
        for e in &new_edges {
            ts.children.push(e.0.clone());
        }
        ts.children.extend(extras);
    }
    true
}

/// STEP post-pass after pcurve association.
/// Source: `TranslateEdgeLoop::CheckPCurves` first check, then
/// `ShapeFix_Wire::Perform` FixReorder (`cxx:317-325` / `cxx:487`) then
/// `FromSTEP.FixShape.FixShiftedMode` (default on) / `cxx:939` FixShifted,
/// then SameRange / CheckPCurveRange / FixAddPCurve / FixSameParameter
/// (`ShapeFix_Wire.cxx:938-995`).
pub fn check_pcurves_and_shift(wire: &Wire, face: &Face) {
    check_pcurve_rep_range(wire, face);
    // `ShapeFix_Wire::Perform` (`cxx:317-325`): FixReorder before FixEdgeCurves.
    // `cxx:368-374` keeps FixShiftedMode when reorder did not FAIL.
    let _ = fix_reorder_wire(wire, face);
    let _ = fix_shifted_wire(wire, face);
    let edges = edges_of_wire(wire);
    let reg = GeometryRegistry::global();
    for (i, e) in edges.iter().enumerate() {
        if let Some((c2d, fp2d, lp2d)) = curve_on_surface_oriented(e, face, false) {
            let (first, last) = BRepTool::edge_parameters(e);
            if first.is_finite() && last.is_finite() {
                let (cfp, clp) = (c2d.first_parameter(), c2d.last_parameter());
                if (first - fp2d).abs() > PCONFUSION || (last - lp2d).abs() > PCONFUSION {
                    reg.set_same_range(&e.0, false);
                } else if !check_pcurve_range(first, last, c2d.as_ref()) {
                    // `GeomLib.cxx:862`: same parametrisation length uses SameRange.
                    // Different length stays on CheckPCurveRange -> RemovePCurve +
                    // FixAddPCurve (`ShapeFix_Wire.cxx:973-983`).
                    if cfp.is_finite()
                        && clp.is_finite()
                        && !c2d.is_periodic()
                        && ((clp - cfp) - (last - first)).abs() <= PCONFUSION
                    {
                        let remapped = geom_lib_same_range(c2d, cfp, clp, first, last);
                        replace_pcurve(e, face, remapped);
                        reg.set_same_range(&e.0, true);
                    } else {
                        let is_seam = is_seam_use(&edges, i);
                        reg.remove_pcurves_on_surface(&e.0, &face.0);
                        let _ = fix_add_pcurve(e, face, is_seam);
                    }
                }
            }
        }
        fix_same_parameter(e, face);
    }
}
