use super::prelude::*;
use super::*;

use occt_core::elib::{clib, slib};
use occt_core::gp::{GpCirc, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION, Precision};
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d::{curve_on_surface_oriented, replace_pcurve};
use crate::brep_surface::SurfaceKind;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::pcurve_full::classify_surface_kind;
use crate::topo_tools_full::{edge_vertices, edges_of_wire, is_same};

/// When the 3D range is a proper subset of the pcurve geometry domain, set
/// the COS representation back to that domain so `FixEdgeCurves` sees
/// `First != fp2d` and `TempSameRange` remaps (`ShapeFix_Wire.cxx:967-971`).
pub fn restore_pcurve_geom_range(wire: &Wire, face: &Face) {
    let Some(surf) = BRepTool::face_surface(face) else {
        return;
    };
    if classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane {
        return;
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    for e in edges_of_wire(wire) {
        let Some((pc, _, _)) = curve_on_surface_oriented(&e, face, false) else {
            continue;
        };
        if pc.is_periodic() {
            continue;
        }
        let (cf, cl) = (pc.first_parameter(), pc.last_parameter());
        if !cf.is_finite() || !cl.is_finite() {
            continue;
        }
        let (first, last) = BRepTool::edge_parameters(&e);
        if !first.is_finite() || !last.is_finite() {
            continue;
        }
        if first + PCONFUSION < cf || last > cl + PCONFUSION {
            continue;
        }
        if (first - cf).abs() <= PCONFUSION && (last - cl).abs() <= PCONFUSION {
            continue;
        }
        if last - first <= PCONFUSION {
            continue;
        }
        reg.set_pcurve_range(&e.0, face_key, cf, cl);
    }
}

/// `TranslateEdgeLoop.cxx:844-868`: `ShapeFix_EdgeProjAux::Compute` then
/// `B.Range(edge, face, First, Last)` or `RemovePCurve`.
pub fn project_wire_pcurve_ranges(wire: &Wire, face: &Face, preci: f64) {
    let Some(surf) = BRepTool::face_surface(face) else {
        return;
    };
    if classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane {
        return;
    }
    // `TranslateEdgeLoop.cxx:844`: EdgeProjAux runs only when
    // `!aTool.ComputePCurve()`. That flag becomes true if any loop edge
    // lacked a STEP pcurve (`cxx:570/590/670/681`).
    if edges_of_wire(wire)
        .iter()
        .any(|e| curve_on_surface_oriented(e, face, false).is_none())
    {
        return;
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    for e in edges_of_wire(wire) {
        match edge_proj_aux_compute(&e, face, preci) {
            Some((first, last, _need_rev)) if (first - last).abs() >= PCONFUSION => {
                // `TranslateEdgeLoop.cxx:853-859`: whenever `IsFirstDone` and
                // `IsLastDone` and the window is non-degenerate, the loop
                // writes `B.Range(edge, Face, FirstParam, LastParam)` — the
                // pcurve range only (`BRep_Builder.hxx:311-315`); the 3d range
                // of the edge is untouched. The window is written verbatim
                // even when it sits a whole period away from the 3d range
                // (Torus V-seam: Init2d `(-pi,pi)` vs 3d `(pi,3pi)`): OCCT does
                // the same, and the resulting counts match `data/occ-torus.obj`
                // exactly (1369/2592), plus linkrods 3474/5025 vs occ 3494/5078
                // (was 4313/6386 before this write) and Shape-2 3091/4764 vs occ
                // 3098/4778 (was 2815/4034).
                //
                // Reverse of the pcurve (`ShapeFix_EdgeProjAux.cxx:455-457` and
                // the `ReversedParameter` arm `cxx:679-683`) stays deferred:
                // applying it collapses the Torus face to a single parallel
                // (74 v / 72 f, all vertices on r=8, z=0) because our mesher
                // mis-handles a reversed pcurve whose stored window is
                // period-shifted from the 3d range. The Init2d control flow
                // above is arm-for-arm identical to cxx (verified against the
                // `cxx:275-400` Line-cut, `cxx:449-471` end tests, `cxx:486-533`
                // period arm and `cxx:619-683` UpdateParam2d), so the missing
                // piece is downstream (pcurve-range consumers / mesh side), not
                // an Init2d branch. No invented predicate is added here.
                reg.set_pcurve_range(&e.0, face_key, first, last);
            }
            Some(_) => {}
            None => {
                let has_pc = curve_on_surface_oriented(&e, face, false).is_some();
                if has_pc {
                    reg.remove_pcurves_on_surface(&e.0, &face.0);
                }
            }
        }
    }
}

/// `TranslateEdgeLoop.cxx:175` `XSAlgo_ShapeProcessor::CheckPCurve` per edge.
pub fn check_pcurves_xsalgo(wire: &Wire, face: &Face, preci: f64) {
    let Some(surf) = BRepTool::face_surface(face) else {
        return;
    };
    if classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane {
        return;
    }
    let edges = edges_of_wire(wire);
    for (i, e) in edges.iter().enumerate() {
        let _ = check_pcurve(e, face, preci, is_seam_use(&edges, i));
    }
}

/// `ShapeFix_EdgeProjAux::Compute` (`cxx:76-103`) via `Init2d`.
pub fn edge_proj_aux_params(edge: &Edge, face: &Face, preci: f64) -> Option<(f64, f64)> {
    let (u1, u2, _) = edge_proj_aux_compute(edge, face, preci)?;
    Some((u1, u2))
}

/// Returns `(first, last, need_reverse)`. Reverse is not applied here; the
/// caller applies `replace_pcurve` only when accepting the Init2d window.
fn edge_proj_aux_compute(edge: &Edge, face: &Face, preci: f64) -> Option<(f64, f64, bool)> {
    edge_proj_aux_compute_inner(edge, face, preci)
}

fn edge_proj_aux_compute_inner(edge: &Edge, face: &Face, preci: f64) -> Option<(f64, f64, bool)> {
    let (u1, u2, need_rev) = init2d(edge, face, preci)?;
    if u1 >= u2 {
        Some((u2, u1, need_rev))
    } else {
        Some((u1, u2, need_rev))
    }
}

fn is_geom2d_line(c: &dyn Curve2d) -> bool {
    !c.first_parameter().is_finite()
        && !c.last_parameter().is_finite()
        && c.d2(0.0).2.square_magnitude() < 1e-30
}

fn reverse_circ_dir(circ: &mut GpCirc) {
    let mut ax = circ.position();
    ax.set_direction(ax.direction().reversed());
    circ.set_position(&ax);
}

/// `Adaptor3d_CurveOnSurface::EvalKPart` (`cxx:1579-1728`) when the pcurve
/// is a 2d line on a sphere / cylinder / torus iso.
fn cons_circle(pc: &dyn Curve2d, surf: &dyn Surface) -> Option<GpCirc> {
    if !is_geom2d_line(pc) {
        return None;
    }
    let loc = pc.d0(0.0);
    let dir = pc.d1(0.0).1;
    let (dx, dy) = (dir.x(), dir.y());
    if dy.abs() <= ANGULAR && dx.abs() > ANGULAR {
        let v = loc.y();
        let (mut circ, zaxis) = if let Some(t) = surf.gp_torus() {
            (
                slib::torus_v_iso(t.position(), t.major_radius(), t.minor_radius(), v),
                t.position().axis().clone(),
            )
        } else if let Some(s) = surf.gp_sphere() {
            (
                slib::sphere_v_iso(s.position(), s.radius(), v),
                s.position().axis().clone(),
            )
        } else if let Some(c) = surf.gp_cylinder() {
            (
                slib::cylinder_v_iso(&c.position(), c.radius(), v),
                c.position().axis().clone(),
            )
        } else {
            return None;
        };
        circ.rotate(&zaxis, loc.x());
        if dx < 0.0 {
            reverse_circ_dir(&mut circ);
        }
        Some(circ)
    } else if dx.abs() <= ANGULAR && dy.abs() > ANGULAR {
        if let Some(t) = surf.gp_torus() {
            let mut circ =
                slib::torus_u_iso(t.position(), t.major_radius(), t.minor_radius(), loc.x());
            let axis = circ.position().axis().clone();
            circ.rotate(&axis, loc.y());
            if dy < 0.0 {
                reverse_circ_dir(&mut circ);
            }
            Some(circ)
        } else if let Some(s) = surf.gp_sphere() {
            let mut circ = slib::sphere_u_iso(s.position(), s.radius(), 0.0);
            let yrev = s
                .position()
                .x_direction()
                .crossed(&s.position().direction())
                .ok()?;
            let axe_y = occt_core::gp::GpAx1::new(s.position().location(), yrev);
            circ.rotate(&axe_y, loc.y());
            circ.rotate(s.position().axis(), loc.x());
            if dy < 0.0 {
                reverse_circ_dir(&mut circ);
            }
            Some(circ)
        } else {
            None
        }
    } else {
        None
    }
}

/// `Geom2d_Curve::ReversedParameter` (`Geom2d_Curve.hxx`; the concrete
/// overrides are `FirstParameter() + LastParameter() - U`, e.g.
/// `Geom2d_Circle::ReversedParameter`). Shared with
/// `ShapeFix_Wire::CopyReversePcurves` (`ShapeFix_Wire.cxx:4154-4156`).
pub(super) fn reversed_parameter(c: &dyn Curve2d, u: f64) -> f64 {
    let (cf, cl) = (c.first_parameter(), c.last_parameter());
    if cf.is_finite() && cl.is_finite() {
        cf + cl - u
    } else {
        -u
    }
}

fn cut_infinite_line(c: &dyn Curve2d, surf: &dyn Surface) -> (f64, f64, bool, bool) {
    let (uf, ul) = surf.u_range();
    let (mut vf, mut vl) = surf.v_range();
    // `ShapeFix_EdgeProjAux.cxx:281-296`: the bounds are clamped to +/-23 only
    // when the basis curve is a hyperbola. UNPORTED: the `Curve` trait has no
    // `GeomAbs_Hyperbola` query, so the extrusion arm (`cxx:281-288`) is not
    // ported and the revolution arm below uses a fully infinite basis range as
    // its stand-in for the hyperbola test.
    if surf.is_surface_of_revolution() {
        if let Some(basis) = surf.revolution_basis_curve() {
            let (bf, bl) = (basis.first_parameter(), basis.last_parameter());
            if !bf.is_finite() && !bl.is_finite() {
                vf = vf.max(-23.0);
                vl = vl.min(23.0);
            }
        }
    }
    let pnt = c.d0(0.0);
    let dir = c.d1(0.0).1;
    let dx = dir.x();
    let dy = dir.y();
    let mut par_u = false;
    let mut par_v = false;
    if uf.is_finite() && ul.is_finite() && vf.is_finite() && vl.is_finite() {
        let (cfi, cli) = if dy == 0.0 && dx != 0.0 {
            par_u = true;
            ((uf - pnt.x()) / dx, (ul - pnt.x()) / dx)
        } else if dx == 0.0 && dy != 0.0 {
            par_v = true;
            ((vf - pnt.y()) / dy, (vl - pnt.y()) / dy)
        } else if dx != 0.0 && dy != 0.0 {
            let xfi = (uf - pnt.x()) / dx;
            let xli = (ul - pnt.x()) / dx;
            let yfi = (vf - pnt.y()) / dy;
            let yli = (vl - pnt.y()) / dy;
            if dx * dy > 0.0 {
                (
                    if (xli - xfi).abs() < (xli - yfi).abs() {
                        xfi
                    } else {
                        yfi
                    },
                    if (xfi - xli).abs() < (xfi - yli).abs() {
                        xli
                    } else {
                        yli
                    },
                )
            } else {
                (
                    if (xli - xfi).abs() < (xli - yli).abs() {
                        xfi
                    } else {
                        yli
                    },
                    if (yli - xli).abs() < (yli - yfi).abs() {
                        xli
                    } else {
                        yfi
                    },
                )
            }
        } else {
            return (-10000.0, 10000.0, false, false);
        };
        return if cfi < cli {
            (cfi, cli, par_u, par_v)
        } else {
            (cli, cfi, par_u, par_v)
        };
    }
    if uf.is_finite() && ul.is_finite() {
        if dx != 0.0 {
            if dy == 0.0 {
                par_u = true;
            }
            let cfi = (uf - pnt.x()) / dx;
            let cli = (ul - pnt.x()) / dx;
            return if cfi < cli {
                (cfi, cli, par_u, par_v)
            } else {
                (cli, cfi, par_u, par_v)
            };
        }
        return (-10000.0, 10000.0, par_u, par_v);
    }
    (-10000.0, 10000.0, par_u, par_v)
}

fn cons_d0(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> GpPnt {
    let uv = pc.d0(t);
    surf.d0(uv.x(), uv.y())
}

fn cons_d1(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> (GpPnt, GpVec) {
    let (uv, duv) = pc.d1(t);
    let (p, su, sv) = surf.d1(uv.x(), uv.y());
    let tan = su
        .multiplied_scalar(duv.x())
        .added(&sv.multiplied_scalar(duv.y()));
    (p, tan)
}

/// `ShapeAnalysis_Curve::Project` on `Adaptor3d_CurveOnSurface` with
/// `AdjustToEnds=false` (`ShapeAnalysis_Curve.cxx:205-261`).
fn project_cons(
    pc: &dyn Curve2d,
    surf: &dyn Surface,
    pt: &GpPnt,
    cf: f64,
    cl: f64,
) -> Option<(f64, f64)> {
    if !cf.is_finite() || !cl.is_finite() {
        return project_act(pc, surf, pt, cf, cl);
    }
    let low = cons_d0(pc, surf, cf);
    let high = cons_d0(pc, surf, cl);
    let dl = low.distance(pt);
    let dh = high.distance(pt);
    if dl <= CONFUSION {
        return Some((cf, dl));
    }
    if dh <= CONFUSION {
        return Some((cl, dh));
    }
    let (param, dist) = project_act(pc, surf, pt, cf, cl)?;
    if dist < dl + CONFUSION && dist < dh + CONFUSION {
        return Some((param, dist));
    }
    if dl < dh {
        Some((cf, dl))
    } else {
        Some((cl, dh))
    }
}

/// `ShapeAnalysis_Curve::ProjectAct` extrema scan on CurveOnSurface.
fn project_act(
    pc: &dyn Curve2d,
    surf: &dyn Surface,
    pt: &GpPnt,
    cf: f64,
    cl: f64,
) -> Option<(f64, f64)> {
    let (lo, hi) = if cf.is_finite() && cl.is_finite() && cl > cf {
        (cf, cl)
    } else {
        (-10000.0, 10000.0)
    };
    if let Some(circ) = cons_circle(pc, surf) {
        let pos = circ.position();
        let mut w = clib::circle_parameter(&pos, pt);
        if lo.is_finite() && hi.is_finite() && hi > lo && (w < lo || w > hi) {
            w += crate::shhealing::adjust_by_period(w, 0.5 * (lo + hi), hi - lo);
        }
        let q = clib::circle_value(&circ, w);
        return Some((w, pt.distance(&q)));
    }
    let mut best_t = lo;
    let mut best_d = cons_d0(pc, surf, lo).distance(pt);
    const N: i32 = 64;
    for i in 1..=N {
        let t = lo + (hi - lo) * i as f64 / N as f64;
        let d = cons_d0(pc, surf, t).distance(pt);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    let mut t = best_t;
    for _ in 0..16 {
        let (q, der) = cons_d1(pc, surf, t);
        let speed2 = der.square_magnitude();
        if speed2 <= 1e-20 {
            break;
        }
        let num = (q.x() - pt.x()) * der.x() + (q.y() - pt.y()) * der.y() + (q.z() - pt.z()) * der.z();
        let next = (t - num / speed2).clamp(lo, hi);
        if !next.is_finite() || (next - t).abs() < 1e-12 {
            t = next;
            break;
        }
        let d_next = cons_d0(pc, surf, next).distance(pt);
        if d_next > best_d {
            t = best_t;
            break;
        }
        t = next;
        best_d = d_next;
        best_t = next;
    }
    let dist = cons_d0(pc, surf, t).distance(pt);
    if Precision::is_infinite(dist) {
        return None;
    }
    Some((t, dist))
}

fn update_param2d(c: &dyn Curve2d, first: &mut f64, last: &mut f64) -> bool {
    if *first < *last {
        return false;
    }
    let cf = c.first_parameter();
    let cl = c.last_parameter();
    let preci2d = PCONFUSION;
    if c.is_periodic() && cf.is_finite() && cl.is_finite() {
        crate::geom_bnd_lib_elclib2d::adjust_periodic(cf, cl, preci2d, first, last);
        return false;
    }
    let closed = cf.is_finite()
        && cl.is_finite()
        && c.d0(cf).distance(&c.d0(cl)) <= preci2d;
    if closed {
        if (*first - cl).abs() <= preci2d {
            *first = cf;
            return false;
        }
        if (*last - cf).abs() <= preci2d {
            *last = cl;
            return false;
        }
        return false;
    }
    true
}

/// `ShapeFix_EdgeProjAux::Init2d` (`cxx:208-536`).
/// Returns `(first, last, need_reverse)` without mutating the edge pcurve.
fn init2d(edge: &Edge, face: &Face, preci: f64) -> Option<(f64, f64, bool)> {
    let surf = BRepTool::face_surface(face)?;
    let (pc, _, _) = curve_on_surface_oriented(edge, face, false)?;
    let (pt1, pt2) = if let Some(c3d) = BRepTool::edge_curve(edge) {
        let (a, b) = BRepTool::edge_parameters(edge);
        (c3d.d0(a), c3d.d0(b))
    } else {
        let (v1, v2) = edge_vertices(edge);
        (BRepTool::vertex_point(&v1?), BRepTool::vertex_point(&v2?))
    };
    let (v1, v2) = edge_vertices(edge);
    if let (Some(v1), Some(v2)) = (v1.as_ref(), v2.as_ref()) {
        if is_same(&v1.0, &v2.0) {
            let dg = degenerated_values(surf.as_ref(), &pt1, preci);
            if let Some((a_pt1, a_pt2, firstpar, lastpar)) = dg
            {
                if is_geom2d_line(pc.as_ref())
                    && a_pt1.distance(&pc.d0(firstpar)) <= preci
                    && a_pt2.distance(&pc.d0(lastpar)) <= preci
                {
                    return Some((firstpar, lastpar, false));
                }
            }
        }
    }

    let mut cf = pc.first_parameter();
    let mut cl = pc.last_parameter();
    let mut par_u = false;
    let mut par_v = false;
    if !cf.is_finite() || !cl.is_finite() {
        if is_geom2d_line(pc.as_ref()) {
            let cut = cut_infinite_line(pc.as_ref(), surf.as_ref());
            cf = cut.0;
            cl = cut.1;
            par_u = cut.2;
            par_v = cut.3;
        } else {
            // `cxx:386-399` BSpline Reparametrize on infinite non-line is
            // unported (`BSplCLib::Reparametrize` on a Geom2d_BSpline copy).
            cf = -10000.0;
            cl = 10000.0;
        }
    }

    // `cxx:421-439`: `Project` returns the distance only to guard
    // `Precision::IsInfinite(dist)`; `project_cons` already maps that case to
    // `None`, so the value itself has no further use.
    let (w1, _d1) = project_cons(pc.as_ref(), surf.as_ref(), &pt1, cf, cl)?;
    let (w2, _d2) = project_cons(pc.as_ref(), surf.as_ref(), &pt2, cf, cl)?;
    let mut my_first = w1;
    let mut my_last = w2;
    if (w1 - w2).abs() < PCONFUSION && !surf.is_u_periodic() && !surf.is_v_periodic() {
        return Some((my_first, my_last, false));
    }
    if my_first == cf && my_last == cl {
        return Some((my_first, my_last, false));
    }
    if my_first == cl && my_last == cf {
        // `cxx:455-457` Reverse — deferred until `B.Range` accepts this window.
        my_first = reversed_parameter(pc.as_ref(), cl);
        my_last = reversed_parameter(pc.as_ref(), cf);
        return Some((my_first, my_last, true));
    }
    if cons_d0(pc.as_ref(), surf.as_ref(), cf)
        .distance(&cons_d0(pc.as_ref(), surf.as_ref(), cl))
        < CONFUSION
    {
        if (my_first - cf).abs() < PCONFUSION && (my_last - cf).abs() < PCONFUSION {
            my_last = cl;
        } else if (my_first - cl).abs() < PCONFUSION && (my_last - cl).abs() < PCONFUSION {
            my_first = cf;
        }
    }
    if par_u || par_v {
        let (uf, ul) = surf.u_range();
        let (vf, vl) = surf.v_range();
        let period = if par_u { ul - uf } else { vl - vf };
        // `cxx:486-489` AdjustToPeriod(w, 0, period) == AdjustByPeriod(w, period/2, period).
        my_first += adjust_by_period(my_first, 0.5 * period, period);
        my_last += adjust_by_period(my_last, 0.5 * period, period);
        if let Some(c3d) = BRepTool::edge_curve(edge) {
            let (a, b) = BRepTool::edge_parameters(edge);
            let mid = c3d.d0((a + b) / 2.0);
            if let Some((mut wmid, _)) = project_cons(pc.as_ref(), surf.as_ref(), &mid, cf, cl) {
                wmid += adjust_by_period(wmid, 0.5 * period, period);
                // `ShapeFix_EdgeProjAux.cxx:500-533` (equal ends use `w1 >= w2`).
                // Torus V-seam LINE Point(t)=(0,t): the cut window is the surface
                // V bounds `cf=0, cl=2pi`, both vertices project to `pi` (cxx:449
                // is skipped because `pi != 0/2pi`), the COS ends coincide
                // (`cxx:462`), and `Project(mid)` returns `uMin=cf=0` because the
                // mid point equals the low bound within `Confusion`
                // (`ShapeAnalysis_Curve.cxx:231-236`), so cxx:503-504 yields
                // `(-pi, pi)`. Measured: `cf=0 cl=2pi w1=w2=pi wmid=0` -> window
                // `(-pi,pi)`; u-seam: `cf=-2pi cl=-0 w1=w2=-2pi wmid=pi` ->
                // `(-2pi,0)`. Both windows are written verbatim (see the caller).
                if my_first >= my_last {
                    if my_last > wmid {
                        my_first -= period;
                    } else if my_first > wmid {
                        let _ = update_param2d(pc.as_ref(), &mut my_first, &mut my_last);
                    } else {
                        my_last += period;
                    }
                } else if my_first > wmid {
                    my_last -= period;
                    let _ = update_param2d(pc.as_ref(), &mut my_first, &mut my_last);
                } else if my_last < wmid {
                    my_first += period;
                    let _ = update_param2d(pc.as_ref(), &mut my_first, &mut my_last);
                }
            }
        } else {
            let _ = update_param2d(pc.as_ref(), &mut my_first, &mut my_last);
            return Some((my_first, my_last, false));
        }
    }
    // `cxx:535` UpdateParam2d — Reverse deferred to caller (`cxx:673-682`).
    let do_reverse = update_param2d(pc.as_ref(), &mut my_first, &mut my_last);
    if do_reverse {
        let tmp1 = my_first;
        let tmp2 = my_last;
        my_first = reversed_parameter(pc.as_ref(), tmp1);
        my_last = reversed_parameter(pc.as_ref(), tmp2);
        // Reverse deferred — caller applies only when accepting this window.
        return Some((my_first, my_last, true));
    }
    Some((my_first, my_last, false))
}

fn make_edge_on_curve(edge: &Edge) -> Option<Edge> {
    let c3d = BRepTool::edge_curve(edge)?;
    let (first, last) = BRepTool::edge_parameters(edge);
    if !first.is_finite() || !last.is_finite() {
        return None;
    }
    let b = TopoBuilder::new();
    let mut tmp = b.make_edge(c3d.clone(), first, last);
    let v1 = b.make_vertex(c3d.d0(first), 0.0);
    let v2 = b.make_vertex(c3d.d0(last), 0.0);
    b.add_edge_vertices(&mut tmp, &v1, &v2);
    Some(tmp)
}

/// `ShapeAnalysis_Edge::PCurve(TopoDS::Edge(theEdge.Reversed()), theFace, C,
/// ..., false)`: the curve `BRep_Tool::CurveOnSurface` returns for the
/// *flipped* orientation of `edge`, i.e. `PCurve2` instead of `PCurve` when the
/// representation is a curve on a closed surface (`BRep_Tool.cxx:354-361`). On a
/// seam (closed surface) that is the other curve of the pair; when the edge
/// carries a single pcurve the rule collapses to that same curve, which is what
/// `XSAlgo_ShapeProcessor.cxx:422-427` / `:473-478` compare against the forward
/// one before falling back to a copy.
///
/// `TopoDS::Edge(theEdge.Reversed())` is `TopoDS_Shape::Reversed()`
/// (`TopoDS_Shape.lxx`), which *flips* `FORWARD <-> REVERSED`; it is not
/// `Oriented(TopAbs_REVERSED)`, which *sets* the orientation and would leave an
/// already-REVERSED edge unchanged.
fn pcurve_on_reversed(edge: &Edge, face: &Face) -> Option<Arc<dyn Curve2d>> {
    let mut flipped = edge.0.clone();
    flipped.reverse();
    curve_on_surface_oriented(&Edge(flipped), face, false).map(|(pc, _, _)| pc)
}

/// `XSAlgo_ShapeProcessor::CheckPCurve` (`XSAlgo_ShapeProcessor.cxx:344-505`).
/// See `CHECK_PCURVE_REPROJECT` for the parked tail of that range.
const CHECK_PCURVE_REPROJECT: bool = true;

fn check_pcurve(edge: &Edge, face: &Face, preci: f64, is_seam: bool) -> bool {
    let Some((c2d, p1, p2)) = curve_on_surface_oriented(edge, face, false) else {
        return false;
    };
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let uv1 = c2d.d0(p1);
    let uv2 = c2d.d0(p2);
    let (u1, u2) = surf.u_range();
    let (v1, v2) = surf.v_range();
    let span_x = (uv1.x() - uv2.x()).abs();
    let span_y = (uv1.y() - uv2.y()).abs();
    if span_x / 8.0 > (u2 / 6.0 - u1 / 6.0) || span_y / 8.0 > (v2 / 6.0 - v1 / 6.0) {
        GeometryRegistry::global().remove_pcurves_on_surface(&edge.0, &face.0);
        return false;
    }
    let c3d_pt1;
    let c3d_pt2;
    let (a, b) = BRepTool::edge_parameters(edge);
    if let Some(c3d) = BRepTool::edge_curve(edge) {
        c3d_pt1 = c3d.d0(a);
        c3d_pt2 = c3d.d0(b);
    } else {
        let (v1, v2) = (first_vertex(edge), last_vertex(edge));
        let (Some(v1), Some(v2)) = (v1, v2) else {
            return false;
        };
        c3d_pt1 = BRepTool::vertex_point(&v1);
        c3d_pt2 = BRepTool::vertex_point(&v2);
    }
    let s1 = surf.d0(uv1.x(), uv1.y());
    let s2 = surf.d0(uv2.x(), uv2.y());
    let d1 = c3d_pt1.distance(&s1);
    let d2 = c3d_pt2.distance(&s2);
    if d1 > preci || d2 > preci {
        GeometryRegistry::global().remove_pcurves_on_surface(&edge.0, &face.0);
        return false;
    }
    // PARKED (t336): `XSAlgo_ShapeProcessor.cxx:403-505` -- the temporary-edge
    // `ShapeFix_Edge::FixSameParameter`, the `FixAddPCurve` re-projection and
    // the copy back to the original edge (pcurve, vertices and the `aTolerance`
    // write at `cxx:491`). Reason for parking is the re-projection quality, not
    // this control flow: the only port of the re-projection, `fix_add_pcurve` ->
    // `pcurve_full::project_curve_on_surface_perform`
    // (`pcurve_full/common.rs:417-460`), builds a degree-1 B-spline with one pole
    // per sample instead of the tolerance-driven `Approx_BSplineApprox`
    // reduction of `ShapeConstruct_ProjectCurveOnSurface::ApproxPCurve` (no
    // port of `Approx_BSplineApprox` exists), and `cxx:435-439` forces
    // SameRange(false) so `cxx:457` always prefers that re-projection.
    //
    // Measured (t336) with this gate ON and the whole tail unconditional
    // (`export_data_obj`, T312 also on): all 16 model vertex/face counts are
    // byte-identical to this gate OFF, at +35 s for the run, so the tail is
    // currently inert for counts. It therefore stays OFF, but no longer on the
    // earlier `Offset` 708/876 -> 30264/57316 evidence: that blow-up was
    // `pcurve_on_reversed` (see its doc comment), and `Offset` now measures
    // 712/892 with the T312 gate on either way. The `cxx:392-401` gate above is
    // the part that explains linkrods FACE 35.
    if !CHECK_PCURVE_REPROJECT {
        return true;
    }
    let Some(mut tmp) = make_edge_on_curve(edge) else {
        return false;
    };
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    if is_seam {
        // `XSAlgo_ShapeProcessor.cxx:415-429`: `aSeamPCurve` is
        // `PCurve(theEdge.Reversed(), theFace, ..., false)`; a null handle (no
        // second pcurve) or a handle equal to `aCurve2D` makes it a copy of
        // `aCurve2D`.
        let seam = pcurve_on_reversed(edge, face)
            .filter(|pc| !Arc::ptr_eq(pc, &c2d))
            .unwrap_or_else(|| Arc::from(c2d.clone_dyn()));
        reg.set_edge_pcurves(&tmp.0, face_key, vec![c2d.clone(), seam]);
    } else {
        reg.set_edge_pcurve(&tmp.0, face_key, c2d.clone());
    }
    reg.set_pcurve_range(&tmp.0, face_key, p1, p2);
    reg.set_same_range(&tmp.0, false);
    fix_same_parameter(&tmp, face);
    let mut a_tol = BRepTool::edge_tolerance(&tmp);
    let mut same_range = reg.edge_geom(&tmp.0).map(|g| g.same_range).unwrap_or(true);
    let mut same_param = BRepTool::same_parameter(&tmp);
    if a_tol > 1.0_f64.min(2.0 * preci) || !same_range {
        if let Some(pr) = make_edge_on_curve(edge) {
            // `XSAlgo` drives its own `ShapeFix_Edge`, so this cache is fresh.
            let mut cache = crate::pcurve_full::ProjectorCache::new();
            let _ = fix_add_pcurve(&pr, face, is_seam, preci, &mut cache);
            fix_same_parameter(&pr, face);
            let tol_pr = BRepTool::edge_tolerance(&pr);
            if tol_pr < a_tol || !same_range {
                same_range = reg.edge_geom(&pr.0).map(|g| g.same_range).unwrap_or(true);
                same_param = BRepTool::same_parameter(&pr);
                a_tol = tol_pr;
                tmp = pr;
            }
        }
    }
    let Some((out_pc, out_f, out_l)) = curve_on_surface_oriented(&tmp, face, false) else {
        return false;
    };
    // `cxx:465-487`: take the corrected pcurves off the temporary edge and put
    // them back on the original one.
    if is_seam {
        // `aSeamPCurve` is `PCurve(aTmpEdge.Reversed(), theFace, ..., false)`
        // (`cxx:471-477`); a null handle or a handle equal to `aCurve2D` makes
        // it a copy of `aCurve2D`. It is written before `aCurve2D` when the
        // original edge is REVERSED (`cxx:479-486`).
        let seam = pcurve_on_reversed(&tmp, face)
            .filter(|pc| !Arc::ptr_eq(pc, &out_pc))
            .unwrap_or_else(|| Arc::from(out_pc.clone_dyn()));
        if edge.0.orientation().is_reversed() {
            reg.set_edge_pcurves(&edge.0, face_key, vec![seam, out_pc]);
        } else {
            reg.set_edge_pcurves(&edge.0, face_key, vec![out_pc, seam]);
        }
    } else {
        replace_pcurve(edge, face, out_pc);
    }
    // `cxx:491` `aBuilder.UpdateEdge(..., aTolerance)`: `BRep_Builder::UpdateEdge`
    // (`BRep_Builder.cxx:668-671`, two-curve overload `:715-719`) ends in
    // `BRep_TEdge::UpdateTolerance(Tol)`, which only raises
    // (`BRep_TEdge.lxx:33-37`: `if (T > myTolerance) myTolerance = T;`), so the
    // write cannot lower the edge tolerance.
    let tol_cur = BRepTool::edge_tolerance(edge);
    reg.set_edge_tolerance(&edge.0, tol_cur.max(a_tol));
    // `cxx:492-493` `aBuilder.UpdateVertex(aVertex1/aVertex2, aTolerance)`; the
    // vertex write is `ShapeFix_ShapeTolerance`-free raise-only in practice
    // because `UpdateVertex` also routes through `BRep_TVertex::UpdateTolerance`.
    if let Some(v) = first_vertex(edge) {
        v.set_tolerance(BRepTool::vertex_tolerance(&v).max(a_tol));
    }
    if let Some(v) = last_vertex(edge) {
        v.set_tolerance(BRepTool::vertex_tolerance(&v).max(a_tol));
    }
    // `cxx:494-503`.
    reg.set_pcurve_range(&edge.0, face_key, out_f, out_l);
    if reg.edge_geom(&edge.0).map(|g| g.same_range).unwrap_or(true) {
        reg.set_same_range(&edge.0, same_range);
    }
    if BRepTool::same_parameter(edge) {
        reg.set_same_parameter(&edge.0, same_param);
    }
    true
}
