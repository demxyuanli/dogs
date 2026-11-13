//! `IntPatch_TheIWalking` (`IntWalk_IWalking.gxx`) — start-driven walk of F=0.

use occt_core::gp::{GpDir2d, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Surface;

use crate::int_tools_wline::{u_resolution, v_resolution, PntOn2S};
use crate::int_tools_wline::{WLine, WLineWay};
use crate::intpatch::impimp::{normale, ImplicitQuad};
use crate::intsurf_transition::IntSurfTypeTrans;
use crate::geom_int::surface_parameters;

use super::search_inside::InteriorPoint;
use super::surf_func::{function_set_root, SurfFunction};

#[path = "intpatch_iw_line_of_iwalking.rs"]
mod iw_line_of_iwalking;

/// Path start after `ComputeTangency` (`IntSurf_PathPoint`).
#[derive(Clone, Copy)]
pub(crate) struct WalkStart {
    pub p: GpPnt,
    pub u: f64,
    pub v: f64,
    pub d3d: GpVec,
    pub d2d: GpDir2d,
    pub tangent: bool,
}

/// `IntPatch_TheIWLineOfTheIWalking` (the `IntWalk_TheIWLine` of this
/// instantiation, `IntPatch_TheIWLineOfTheIWalking.hxx:248`).
pub(crate) use self::iw_line_of_iwalking::IwLineOfIWalking;

/// `IntWalk_IWalking::Perform(seqpdep, seqpins, Func, PSurf, reversed)`.
pub(crate) fn iwalking_perform(
    path: &[WalkStart],
    interior: &[InteriorPoint],
    func: &mut SurfFunction<'_>,
    prm: &dyn Surface,
    reversed: bool,
    pas: f64,
    fleche: f64,
) -> Vec<IwLineOfIWalking> {
    let (mut um, mut umax) = prm.u_range();
    let (mut vm, mut vmax) = prm.v_range();
    if umax < um {
        std::mem::swap(&mut um, &mut umax);
    }
    if vmax < vm {
        std::mem::swap(&mut vm, &mut vmax);
    }
    if !um.is_finite() || !umax.is_finite() || !vm.is_finite() || !vmax.is_finite() {
        return Vec::new();
    }

    let step_u = pas * (umax - um);
    let step_v = pas * (vmax - vm);
    let tol_u = u_resolution(prm, CONFUSION);
    let tol_v = v_resolution(prm, CONFUSION);

    let mut used_path = vec![false; path.len()];
    let mut used_int = vec![false; interior.len()];
    let mut lines: Vec<IwLineOfIWalking> = Vec::new();

    for (i, st) in path.iter().enumerate() {
        if st.tangent || used_path[i] {
            continue;
        }
        if is_on_lines(&st.p, &lines, reversed) {
            used_path[i] = true;
            continue;
        }
        let Some(mut line) = walk_open(
            st.u,
            st.v,
            st.d2d,
            st.d3d,
            func,
            prm,
            reversed,
            pas,
            fleche,
            um,
            umax,
            vm,
            vmax,
            tol_u,
            tol_v,
        ) else {
            continue;
        };
        line.add_status_first_at(false, true, (i + 1) as i32);
        line.set_tangency_at_beginning(st.tangent);
        used_path[i] = true;
        mark_hits(&line, path, interior, &mut used_path, &mut used_int, reversed);
        bind_last(&mut line, path, Some(i));
        lines.push(line);
    }

    for (i, ip) in interior.iter().enumerate() {
        if used_int[i] {
            continue;
        }
        if is_tangent_ext(
            func, ip.u, ip.v, step_u, step_v, um, umax, vm, vmax,
        ) {
            continue;
        }
        if is_on_lines(&ip.p, &lines, reversed) {
            used_int[i] = true;
            continue;
        }
        let Some(mut line) = walk_closed(
            ip.u,
            ip.v,
            ip.d2d,
            ip.d3d,
            func,
            prm,
            reversed,
            pas,
            fleche,
            um,
            umax,
            vm,
            vmax,
            tol_u,
            tol_v,
        ) else {
            continue;
        };
        used_int[i] = true;
        mark_hits(&line, path, interior, &mut used_path, &mut used_int, reversed);
        bind_first(&mut line, path);
        bind_last(&mut line, path, None);
        lines.push(line);
    }

    lines.into_iter().filter(|l| l.points.len() >= 2).collect()
}

/// Fill quadric UV on walking points and wrap as `WLine` (`ImpPrm` post-walk).
pub(crate) fn lines_to_wlines(
    raw: &[IwLineOfIWalking],
    quad: &ImplicitQuad,
    quad_surf: &dyn Surface,
    prm: &dyn Surface,
    reversed: bool,
) -> Vec<(WLine, IwLineOfIWalking)> {
    let (vmin, vmax) = quad_surf.v_range();
    let recadre_u = matches!(
        quad,
        ImplicitQuad::Cylinder(_) | ImplicitQuad::Cone(_) | ImplicitQuad::Sphere(_)
    );
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;
    let mut out = Vec::new();
    for line in raw {
        if line.points.len() < 2 {
            continue;
        }
        // `IntPatch_ImpPrmIntersection.cxx:877-912`: the transitions of the
        // new `IntPatch_WLine` (`:1073`) come from the walking tangent vector
        // and the two normals at that point.
        let (tgline, trans1, trans2) = line_transitions(line, quad, prm, reversed);
        let mut points = Vec::with_capacity(line.points.len());
        let mut an_u1 = 0.0;
        let mut an_u2 = 0.0;
        for (k, src) in line.points.iter().enumerate() {
            let (up, vp) = if reversed {
                src.parameters_on_s1()
            } else {
                src.parameters_on_s2()
            };
            let mut uq = 0.0;
            let mut vq = 0.0;
            if let Some((u, mut v)) = surface_parameters(quad_surf, &src.p) {
                if vmin.is_finite() && v < vmin && vmin - v < 1.0e-14 {
                    v = vmin;
                }
                if vmax.is_finite() && v > vmax && v - vmax < 1.0e-14 {
                    v = vmax;
                }
                uq = u;
                vq = v;
            }
            let mut u_prm = up;
            let v_prm = vp;
            if k == 0 {
                an_u1 = uq;
                an_u2 = u_prm;
            } else if recadre_u {
                let mut cf = 0.0;
                if (uq - an_u1) > 1.5 * std::f64::consts::PI {
                    while (uq - an_u1) > 1.5 * std::f64::consts::PI + cf * two_pi {
                        cf += 1.0;
                    }
                    uq -= cf * two_pi;
                } else {
                    while (uq - an_u1) < -1.5 * std::f64::consts::PI - cf * two_pi {
                        cf += 1.0;
                    }
                    uq += cf * two_pi;
                }
                while u_prm < an_u2 - 1.5 * std::f64::consts::PI {
                    u_prm += two_pi;
                }
                while u_prm > an_u2 + 1.5 * std::f64::consts::PI {
                    u_prm -= two_pi;
                }
                an_u1 = uq;
                an_u2 = u_prm;
            }
            let pnt = if reversed {
                PntOn2S {
                    p: src.p,
                    u1: u_prm,
                    v1: v_prm,
                    u2: uq,
                    v2: vq,
                }
            } else {
                PntOn2S {
                    p: src.p,
                    u1: uq,
                    v1: vq,
                    u2: u_prm,
                    v2: v_prm,
                }
            };
            points.push(pnt);
        }
        let n = points.len();
        let _ = tgline;
        let wl = WLine {
            base: crate::intpatch_line::IntPatchLine::with_transitions(false, trans1, trans2),
            points,
            vertices: Vec::new(),
            has_first_point: false,
            has_last_point: false,
            first_index: 0,
            last_index: 0,
            arc_on_s1: None,
            arc_on_s2: None,
            u1_period: 0.0,
            v1_period: 0.0,
            u2_period: 0.0,
            v2_period: 0.0,
            creating_way: WLineWay::ImpPrm,
            purger_allowed: true,
        };
        if n >= 2 {
            out.push((wl, line.clone()));
        }
    }
    out
}

/// `IntPatch_ImpPrmIntersection.cxx:877-912`: the walking tangent vector and
/// the normals at point `k` give the line transitions.
///
/// `k` starts as `TangentVector(Index)` (`:877-879`); when the stored index is
/// outside `1..NbPoints` OCCT falls back to the middle point (`:880-886`)
/// while keeping the stored vector, so an unset tangent leaves the two
/// normals to decide (`DotCross` of a null vector is `0.` →
/// `IntSurf_In` / `IntSurf_Out`).
fn line_transitions(
    line: &IwLineOfIWalking,
    quad: &ImplicitQuad,
    prm: &dyn Surface,
    reversed: bool,
) -> (GpVec, IntSurfTypeTrans, IntSurfTypeTrans) {
    let (tgline, index) = line.tangent_vector();
    let nbpts = line.nb_points();
    let k = if index >= 1 && index <= nbpts {
        index
    } else {
        nbpts >> 1
    };
    let Some(p) = line.points.get(k.max(1).saturating_sub(1) as usize) else {
        return (tgline, IntSurfTypeTrans::In, IntSurfTypeTrans::Out);
    };
    let valpt = p.p;
    let (norm1, norm2) = if !reversed {
        let (u2, v2) = p.parameters_on_s2();
        let n1 = normale(quad, &valpt);
        let (_pt, d1u, d1v) = prm.d1(u2, v2);
        (n1, d1u.crossed(&d1v))
    } else {
        let (u1, v1) = p.parameters_on_s1();
        let n2 = normale(quad, &valpt);
        let (_pt, d1u, d1v) = prm.d1(u1, v1);
        (d1u.crossed(&d1v), n2)
    };
    if tgline.dot_cross(&norm2, &norm1) > 0.0 {
        (tgline, IntSurfTypeTrans::Out, IntSurfTypeTrans::In)
    } else {
        (tgline, IntSurfTypeTrans::In, IntSurfTypeTrans::Out)
    }
}

fn walk_open(
    u0: f64,
    v0: f64,
    d2d: GpDir2d,
    d3d: GpVec,
    func: &mut SurfFunction<'_>,
    _prm: &dyn Surface,
    reversed: bool,
    pas: f64,
    fleche: f64,
    um: f64,
    umax: f64,
    vm: f64,
    vmax: f64,
    tol_u: f64,
    tol_v: f64,
) -> Option<IwLineOfIWalking> {
    walk_dir(
        u0, v0, d2d, d3d, 1, false, func, _prm, reversed, pas, fleche, um, umax, vm, vmax, tol_u,
        tol_v,
    )
}

fn walk_closed(
    u0: f64,
    v0: f64,
    d2d: GpDir2d,
    d3d: GpVec,
    func: &mut SurfFunction<'_>,
    _prm: &dyn Surface,
    reversed: bool,
    pas: f64,
    fleche: f64,
    um: f64,
    umax: f64,
    vm: f64,
    vmax: f64,
    tol_u: f64,
    tol_v: f64,
) -> Option<IwLineOfIWalking> {
    let mut fwd = walk_dir(
        u0, v0, d2d, d3d, 1, true, func, _prm, reversed, pas, fleche, um, umax, vm, vmax, tol_u,
        tol_v,
    )?;
    if fwd.points.len() >= 3 {
        let last = fwd.points[fwd.points.len() - 1];
        let (ul, vl) = prm_uv(&last, reversed);
        if (ul - u0).abs() <= tol_u * 20.0 && (vl - v0).abs() <= tol_v * 20.0 {
            return Some(fwd);
        }
    }
    if let Some(mut back) = walk_dir(
        u0, v0, d2d, d3d, -1, false, func, _prm, reversed, pas, fleche, um, umax, vm, vmax, tol_u,
        tol_v,
    ) {
        back.reverse();
        if back.points.len() > 1 {
            back.points.pop();
        }
        back.points.extend(fwd.points);
        // `IntWalk_IWalking::OpenLine` (`gxx:2984-2985`): after reversing the
        // line the tangent vector is stored reversed at the new last point.
        let n = back.nb_points();
        back.set_tangent_vector(d3d.reversed(), n);
        fwd = back;
    }
    if fwd.points.len() >= 2 {
        Some(fwd)
    } else {
        None
    }
}

fn walk_dir(
    u0: f64,
    v0: f64,
    mut d2d: GpDir2d,
    mut d3d: GpVec,
    step_sign: i32,
    stop_closed: bool,
    func: &mut SurfFunction<'_>,
    _prm: &dyn Surface,
    reversed: bool,
    pas: f64,
    fleche: f64,
    um: f64,
    umax: f64,
    vm: f64,
    vmax: f64,
    tol_u: f64,
    tol_v: f64,
) -> Option<IwLineOfIWalking> {
    let _ = func.values(u0, v0);
    if !func.is_tangent() {
        d2d = func.direction2d();
        d3d = func.direction3d();
    }
    let p0 = func.point();
    let mut line = IwLineOfIWalking::new();
    line.add_point(make_pnt(p0, u0, v0, reversed));
    // `IntWalk_IWalking.gxx:1509-1528` (`ComputeOpenLine`) and `:2119`
    // (`ComputeCloseLine`): the direction that produced the first point is
    // stored as the tangent of point 1; a negative step walks the reversed
    // direction.
    line.set_tangent_vector(if step_sign < 0 { d3d.reversed() } else { d3d }, 1);
    let mut pasc = first_step(d2d, pas, um, umax, vm, vmax, tol_u, tol_v);
    if pasc <= 0.0 {
        return None;
    }
    let mut u = u0;
    let mut v = v0;
    let mut nb_div = 0;
    // `aNbBadRootIter` (`gxx:1569`, OCCT 8.0.1): consecutive Root > Tolerance.
    let mut a_nb_bad_root_iter = 0i32;
    let mut prev_d3d = d3d;
    for _ in 0..8000 {
        let mut binf = (um, vm);
        let mut bsup = (umax, vmax);
        let mut uvap = (u, v);
        let cadre = cadrage(
            &mut binf,
            &mut bsup,
            &mut uvap,
            &mut pasc,
            d2d,
            step_sign,
        );
        let Some((un, vn)) = function_set_root(func, uvap, binf, bsup, (tol_u, tol_v)) else {
            pasc *= 0.5;
            nb_div += 1;
            if nb_div > 12 || pas_too_small(pasc, d2d, tol_u, tol_v) {
                break;
            }
            continue;
        };
        if func.root().abs() > func.tolerance() {
            // `IntWalk_IWalking.gxx:1596-1601` (OCCT 8.0.1): stop after
            // `MaxDivision` (60) consecutive bad-root subdivisions.
            pasc *= 0.5;
            a_nb_bad_root_iter += 1;
            if pas_too_small(pasc, d2d, tol_u, tol_v) || a_nb_bad_root_iter > 60 {
                break;
            }
            continue;
        }
        a_nb_bad_root_iter = 0;
        let pn = func.point();
        if line.points.iter().any(|q| q.p.square_distance(&pn) <= CONFUSION * CONFUSION)
        {
            break;
        }
        let new_d3d = if func.is_tangent() {
            d3d
        } else {
            func.direction3d()
        };
        if fleche > 0.0 && line.points.len() >= 1 {
            let last_p = line.points[line.points.len() - 1].p;
            let norme = last_p.square_distance(&pn);
            let fleche_c = prev_d3d
                .normalized()
                .subtracted(&new_d3d.normalized())
                .square_magnitude()
                * norme
                / 64.0;
            if fleche_c > fleche * fleche {
                pasc *= 0.5;
                nb_div += 1;
                if nb_div > 12 || pas_too_small(pasc, d2d, tol_u, tol_v) {
                    break;
                }
                continue;
            }
            // gxx:2799-2815 (`TestDeflection`) raises the step when the current
            // deflection is small. Here `pasc` is unconditionally recomputed by
            // `first_step` at the end of the iteration, so that store is dead in
            // this port; kept as a comment to preserve the OCCT branch.
        }
        nb_div = 0;
        line.points.push(make_pnt(pn, un, vn, reversed));
        u = un;
        v = vn;
        if !func.is_tangent() {
            d2d = func.direction2d();
            d3d = func.direction3d();
        }
        prev_d3d = d3d;
        if stop_closed && line.points.len() >= 4 {
            if (u - u0).abs() <= tol_u * 20.0 && (v - v0).abs() <= tol_v * 20.0 {
                line.points.push(make_pnt(p0, u0, v0, reversed));
                break;
            }
        }
        if cadre {
            // `IntWalk_IWalking::TestArretCadre` (`gxx:1290-1314`): when the
            // walk stops on the natural bound of the surface, the tangent
            // vector is recorded at the last point, guarded by
            // `Indextg > NBP` and `j > 3 && j <= NBP + 1` (`gxx:1294-1306`).
            let nbp = line.nb_points();
            if line.indextg > nbp {
                let j = nbp + 1;
                if j > 3 && j <= nbp + 1 {
                    line.set_tangent_vector(prev_d3d, j - 1);
                }
            }
            break;
        }
        pasc = first_step(d2d, pas, um, umax, vm, vmax, tol_u, tol_v);
    }
    if line.points.len() >= 2 {
        Some(line)
    } else {
        None
    }
}

fn first_step(
    d2d: GpDir2d,
    pas: f64,
    um: f64,
    umax: f64,
    vm: f64,
    vmax: f64,
    tol_u: f64,
    tol_v: f64,
) -> f64 {
    let dx = d2d.x().abs();
    let dy = d2d.y().abs();
    if dx < tol_u {
        if dy <= 0.0 {
            return 0.0;
        }
        pas * (vmax - vm) / dy
    } else if dy < tol_v {
        pas * (umax - um) / dx
    } else {
        pas * ((umax - um) / dx).min((vmax - vm) / dy)
    }
}

fn pas_too_small(pasc: f64, d2d: GpDir2d, tol_u: f64, tol_v: f64) -> bool {
    pasc.abs() * d2d.x().abs() <= tol_u && pasc.abs() * d2d.y().abs() <= tol_v
}

/// `IntWalk_IWalking::Cadrage`.
fn cadrage(
    binf: &mut (f64, f64),
    bsup: &mut (f64, f64),
    uvap: &mut (f64, f64),
    step: &mut f64,
    d2d: GpDir2d,
    step_sign: i32,
) -> bool {
    let duvx = d2d.x();
    let duvy = d2d.y();
    let sg = step_sign as f64;
    let u1 = uvap.0 + *step * duvx * sg;
    let v1 = uvap.1 + *step * duvy * sg;
    let infu = u1 <= binf.0 + PCONFUSION;
    let supu = u1 >= bsup.0 - PCONFUSION;
    let infv = v1 <= binf.1 + PCONFUSION;
    let supv = v1 >= bsup.1 - PCONFUSION;
    if !infu && !supu && !infv && !supv {
        uvap.0 = u1;
        uvap.1 = v1;
        return false;
    }
    if (infu || supu) && (infv || supv) {
        let step_u = if infu {
            if duvx.abs() > 0.0 {
                (binf.0 - uvap.0).abs() / duvx.abs()
            } else {
                *step
            }
        } else if duvx.abs() > 0.0 {
            (bsup.0 - uvap.0).abs() / duvx.abs()
        } else {
            *step
        };
        let step_v = if infv {
            if duvy.abs() > 0.0 {
                (binf.1 - uvap.1).abs() / duvy.abs()
            } else {
                *step
            }
        } else if duvy.abs() > 0.0 {
            (bsup.1 - uvap.1).abs() / duvy.abs()
        } else {
            *step
        };
        if step_u <= step_v {
            *step = step_u;
            if infu {
                uvap.0 = binf.0;
                bsup.0 = binf.0;
            } else {
                uvap.0 = bsup.0;
                binf.0 = bsup.0;
            }
            uvap.1 += *step * duvy * sg;
        } else {
            *step = step_v;
            if infv {
                uvap.1 = binf.1;
                bsup.1 = binf.1;
            } else {
                uvap.1 = bsup.1;
                binf.1 = bsup.1;
            }
            uvap.0 += *step * duvx * sg;
        }
        return true;
    }
    if infu {
        if duvx.abs() > 0.0 {
            let a = (binf.0 - uvap.0).abs() / duvx.abs();
            if a < *step {
                *step = a;
            }
        }
        bsup.0 = binf.0;
        uvap.0 = binf.0;
        uvap.1 += *step * duvy * sg;
        return true;
    }
    if supu {
        if duvx.abs() > 0.0 {
            let a = (bsup.0 - uvap.0).abs() / duvx.abs();
            if a < *step {
                *step = a;
            }
        }
        binf.0 = bsup.0;
        uvap.0 = bsup.0;
        uvap.1 += *step * duvy * sg;
        return true;
    }
    if infv {
        if duvy.abs() > 0.0 {
            let a = (binf.1 - uvap.1).abs() / duvy.abs();
            if a < *step {
                *step = a;
            }
        }
        bsup.1 = binf.1;
        uvap.1 = binf.1;
        uvap.0 += *step * duvx * sg;
        return true;
    }
    if duvy.abs() > 0.0 {
        let a = (bsup.1 - uvap.1).abs() / duvy.abs();
        if a < *step {
            *step = a;
        }
    }
    binf.1 = bsup.1;
    uvap.1 = bsup.1;
    uvap.0 += *step * duvx * sg;
    true
}

fn is_tangent_ext(
    func: &mut SurfFunction<'_>,
    u: f64,
    v: f64,
    step_u: f64,
    step_v: f64,
    uinf: f64,
    usup: f64,
    vinf: f64,
    vsup: f64,
) -> bool {
    let tol = func.tolerance();
    let pu = [(u + step_u).min(usup), (u - step_u).max(uinf), u, u];
    let pv = [v, v, (v + step_v).min(vsup), (v - step_v).max(vinf)];
    for i in 0..4 {
        let _ = func.value(pu[i], pv[i]);
        if func.root().abs() > tol {
            return false;
        }
    }
    true
}

fn make_pnt(p: GpPnt, u: f64, v: f64, reversed: bool) -> PntOn2S {
    if reversed {
        PntOn2S {
            p,
            u1: u,
            v1: v,
            u2: 0.0,
            v2: 0.0,
        }
    } else {
        PntOn2S {
            p,
            u1: 0.0,
            v1: 0.0,
            u2: u,
            v2: v,
        }
    }
}

fn prm_uv(p: &PntOn2S, reversed: bool) -> (f64, f64) {
    if reversed {
        p.parameters_on_s1()
    } else {
        p.parameters_on_s2()
    }
}

fn bind_first(line: &mut IwLineOfIWalking, path: &[WalkStart]) {
    if line.has_first {
        return;
    }
    let Some(first) = line.points.first() else {
        return;
    };
    let eps2 = CONFUSION * CONFUSION * 100.0;
    for (j, st) in path.iter().enumerate() {
        if first.p.square_distance(&st.p) <= eps2 {
            line.add_status_first_at(false, true, (j + 1) as i32);
            line.set_tangency_at_beginning(st.tangent);
            return;
        }
    }
}

fn bind_last(line: &mut IwLineOfIWalking, path: &[WalkStart], skip: Option<usize>) {
    let Some(last) = line.points.last() else {
        return;
    };
    let eps2 = CONFUSION * CONFUSION * 100.0;
    for (j, st) in path.iter().enumerate() {
        if Some(j) == skip {
            continue;
        }
        if last.p.square_distance(&st.p) <= eps2 {
            line.add_status_last_at(true, (j + 1) as i32);
            line.set_tangency_at_end(st.tangent);
            return;
        }
    }
}

fn is_on_lines(p: &GpPnt, lines: &[IwLineOfIWalking], _reversed: bool) -> bool {
    let eps2 = CONFUSION * CONFUSION;
    lines
        .iter()
        .any(|l| l.points.iter().any(|q| q.p.square_distance(p) <= eps2))
}

fn mark_hits(
    line: &IwLineOfIWalking,
    path: &[WalkStart],
    interior: &[InteriorPoint],
    used_path: &mut [bool],
    used_int: &mut [bool],
    _reversed: bool,
) {
    let eps2 = CONFUSION * CONFUSION * 100.0;
    for (i, st) in path.iter().enumerate() {
        if line.points.iter().any(|q| q.p.square_distance(&st.p) <= eps2) {
            used_path[i] = true;
        }
    }
    for (i, ip) in interior.iter().enumerate() {
        if line.points.iter().any(|q| q.p.square_distance(&ip.p) <= eps2) {
            used_int[i] = true;
        }
    }
}
