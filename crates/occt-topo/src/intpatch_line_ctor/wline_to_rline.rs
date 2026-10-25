//! The walking-line to restriction-line statics of `IntPatch_LineConstructor`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_LineConstructor.cxx`
//! - `TestWLineIsARLine` (`cxx:949-1005`)
//! - `TestIfWLineIsRestriction` (`cxx:1008-1071`)
//! - `ProjectOnArc` (`cxx:1073-1094`)
//! - `TestWLineToRLine` (`cxx:1095-1432`)
//!
//! Translation notes:
//! - `aDomain->Init()/More()/Next()/Value()` is walked by [`domain_arcs`],
//!   which materialises one `Arc<Curve2dAdaptor>` per domain arc in iteration
//!   order. The `Arc` address is the port's stand-in for OCCT's handle pointer,
//!   so two vertices resolve to "the same arc" exactly when they land on the
//!   same index, and `arc1 == arc` of `cxx:1258` becomes an address compare
//!   (`int_tools_wline::patch_point::arc_id`).
//! - `IntPatch_RLine::ArcOnS1()/ArcOnS2()` hand back the stored
//!   `Adaptor2d_Curve2d` handle. The port's [`RLine`] keeps the loaded 2D
//!   curve (`c2d`) and its edge (`arc_edge`); `FindParameter(C2d, P, 1e-7)`
//!   only reads `FirstParameter`/`LastParameter`/`D0`, which
//!   `Curve2dAdaptor` forwards to that curve, so `c2d` is the adaptor body.
//! - `RealLast()` (`cxx:1187`) is `DBL_MAX`, kept as `f64::MAX`.
//! - `slin.Value(lastwline)` (`cxx:1107`) is read through a handle while the
//!   sequence still takes the new `RLine`s; the port clones the walking line so
//!   the borrow does not overlap the `push`.
//! - OCCT dereferences `*arc` without a null check in `TestWLineIsARLine`
//!   (`cxx:991`) and `ProjectOnArc` (`cxx:1085`). A null arc cannot occur from
//!   `Perform` (only arcs that resolved or are linked are visited) but the port
//!   answers "not on the arc" instead of panicking, as it does in
//!   [`super::helpers::get_vertex_tolerance`].

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::find_parameter;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_int::{GeomIntLine, RLine};
use crate::int_tools_wline::{ArcId, TransType, WLine};
use crate::intpatch_line::IntPatchLine;
use crate::intsurf_transition::IntSurfSituation;

use super::helpers::compute_parametric_tolerance;

/// One domain arc, materialised once so its address is a stable identity.
fn domain_arcs(d: &mut TopolTool) -> Vec<Arc<Curve2dAdaptor>> {
    let mut arcs = Vec::new();
    d.init();
    while d.more() {
        if let Some(c) = d.value() {
            arcs.push(Arc::new(c.clone()));
        }
        d.next();
    }
    arcs
}

/// [`ArcId`] of a materialised domain arc.
fn arc_token(a: &Arc<Curve2dAdaptor>) -> usize {
    Arc::as_ptr(a) as usize
}

/// `ProjectOnArc(u, v, arc, surf, TolArc, par, dist)` (`cxx:1073-1094`):
/// `Some((par, dist))` when the projected 2D point lies within the parametric
/// tolerance.
fn project_on_arc(
    u: f64,
    v: f64,
    arc: &Curve2dAdaptor,
    surf: &dyn Surface,
    tol_arc: f64,
) -> Option<(f64, f64)> {
    let (_, ad1u, ad1v) = surf.d1(u, v);
    let tol2d = compute_parametric_tolerance(tol_arc, &ad1u, &ad1v);
    let curve: &Arc<dyn Curve2d> = arc.curve()?;
    let p = GpPnt2d::new(u, v);
    let par = find_parameter(curve.as_ref(), &p, 1e-7);
    let dist = curve.d0(par).distance(&p);
    if dist <= tol2d {
        Some((par, dist))
    } else {
        None
    }
}

/// `TestWLineIsARLine(slinref, wlin, tol2d)` (`cxx:949-1005`).
fn test_wline_is_arline(slinref: &[GeomIntLine], wlin: &WLine, tol2d: f64) -> bool {
    let nbpnt = wlin.nb_pnts();
    let indicepnt = nbpnt / 2;
    if indicepnt < 1 {
        return false;
    }
    let p_on_2s = wlin.point(indicepnt);
    let p_on_2s1 = wlin.point(indicepnt + 1);
    for line in slinref {
        let Some(rlin) = line.as_rline() else {
            continue;
        };
        for on_first in [true, false] {
            // `(onFirst && rlin->IsArcOnS1()) || (!onFirst && rlin->IsArcOnS2())`.
            if (on_first && !rlin.arc_on_s1) || (!on_first && !rlin.arc_on_s2) {
                continue;
            }
            let Some(c2d) = rlin.c2d.as_ref() else {
                continue;
            };
            let (mut u, mut v) = if on_first {
                p_on_2s.parameters_on_s1()
            } else {
                p_on_2s.parameters_on_s2()
            };
            let (u1, v1) = if on_first {
                p_on_2s1.parameters_on_s1()
            } else {
                p_on_2s1.parameters_on_s2()
            };
            if indicepnt == 1 {
                u = (u + u1) * 0.5;
                v = (v + v1) * 0.5;
            }
            let p2d = GpPnt2d::new(u, v);
            let par = find_parameter(c2d.as_ref(), &p2d, 1e-7);
            if c2d.d0(par).distance(&p2d) < tol2d {
                return true;
            }
        }
    }
    false
}

/// `TestIfWLineIsRestriction(slinref, wlin, S1, D1, S2, D2, TolArc)`
/// (`cxx:1008-1071`).
#[allow(clippy::too_many_arguments)]
fn test_if_wline_is_restriction(
    slinref: &[GeomIntLine],
    wlin: &WLine,
    s1: &dyn Surface,
    d1: &mut TopolTool,
    s2: &dyn Surface,
    d2: &mut TopolTool,
    tol_arc: f64,
) -> bool {
    let nb_pnts = wlin.nb_pnts();
    let mut allon1 = 0;
    let mut allon2 = 0;
    let mut tol2d1 = 0.0f64;
    let mut tol2d2 = 0.0f64;
    for i in 1..=nb_pnts {
        let (u1, v1, u2, v2) = wlin.point(i).parameters();

        let (_, ad1u, ad1v) = s1.d1(u1, v1);
        let tol = compute_parametric_tolerance(tol_arc, &ad1u, &ad1v);
        if tol > tol2d1 {
            tol2d1 = tol;
        }
        if allon1 + 1 == i && d1.is_the_point_on(GpPnt2d::new(u1, v1), tol, true) {
            allon1 += 1;
        }

        let (_, ad1u, ad1v) = s2.d1(u2, v2);
        let tol = compute_parametric_tolerance(tol_arc, &ad1u, &ad1v);
        if tol > tol2d2 {
            tol2d2 = tol;
        }
        if allon2 + 1 == i && d2.is_the_point_on(GpPnt2d::new(u2, v2), tol, true) {
            allon2 += 1;
        }

        if allon1 != i && allon2 != i {
            break;
        }
    }
    if allon1 == nb_pnts || allon2 == nb_pnts {
        let tol2d = if tol2d1 > tol2d2 { tol2d1 } else { tol2d2 };
        return test_wline_is_arline(slinref, wlin, tol2d);
    }
    false
}

/// The resolved arc / parameter of one walking-line vertex (`arcsResolved(iv)`
/// / `paramsResolved(iv)`, `cxx:1176-1196`).
#[derive(Clone, Copy, Default)]
struct Resolved {
    arc: Option<usize>,
    param: f64,
}

/// `TestWLineToRLine(slinref, slin, Surf1, Dom1, Surf2, Dom2, TolArc)`
/// (`cxx:1095-1432`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn test_wline_to_rline(
    slinref: &[GeomIntLine],
    slin: &mut Vec<GeomIntLine>,
    s1: &dyn Surface,
    d1: &mut TopolTool,
    s2: &dyn Surface,
    d2: &mut TopolTool,
    tol_arc: f64,
) {
    let Some(last) = slin.len().checked_sub(1) else {
        return;
    };
    let Some(wline) = slin[last].as_wline().cloned() else {
        return;
    };

    let nbvtx = wline.nb_vertex();
    if nbvtx < 2 {
        return;
    }
    let param_min_on_line = wline.vertex(1).parameter_on_line() as i32;
    let param_max_on_line = wline.vertex(nbvtx).parameter_on_line() as i32;
    if param_min_on_line >= param_max_on_line {
        return;
    }
    let mid_ind = (param_max_on_line + param_min_on_line) / 2;

    let mut indices_v1: Vec<i32> = Vec::new();
    let mut indices_v2: Vec<i32> = Vec::new();
    for iv in 1..=nbvtx {
        let plin = wline.vertex(iv).parameter_on_line() as i32;
        if plin == param_min_on_line {
            indices_v1.push(iv);
        } else if plin == param_max_on_line {
            indices_v2.push(iv);
        }
    }

    let mut is_rline = false;

    for on_first in [true, false] {
        if (on_first && !wline.has_arc_on_s1()) || (!on_first && !wline.has_arc_on_s2()) {
            continue;
        }
        let (surf, a_domain): (&dyn Surface, &mut TopolTool) = if on_first {
            (s1, d1)
        } else {
            (s2, d2)
        };
        let arcs = domain_arcs(a_domain);

        // -- resolve arcs for vertices not having a link to an arc
        //    (`cxx:1176-1196`).
        let mut resolved = vec![Resolved::default(); (nbvtx + 1) as usize];
        for iv in 1..=nbvtx {
            let vtx = wline.vertex(iv);
            let on_dom = if on_first {
                vtx.is_on_dom_s1()
            } else {
                vtx.is_on_dom_s2()
            };
            if on_dom {
                continue;
            }
            let ip = vtx.parameter_on_line() as i32;
            let (utst, vtst) = if on_first {
                wline.point(ip).parameters_on_s1()
            } else {
                wline.point(ip).parameters_on_s2()
            };
            let mut distmin = f64::MAX;
            for (index, arc) in arcs.iter().enumerate() {
                if let Some((par, dist)) = project_on_arc(utst, vtst, arc, surf, tol_arc) {
                    if dist < distmin {
                        resolved[iv as usize] = Resolved {
                            arc: Some(index),
                            param: par,
                        };
                        distmin = dist;
                    }
                }
            }
        }

        // -- the arc a vertex is linked with (`cxx:1202-1237`): the vertex'
        //    own arc when it is on the domain, otherwise the resolved one.
        let vertex_arc = |iv: i32| -> ArcId {
            let vtx = wline.vertex(iv);
            let on_dom = if on_first {
                vtx.is_on_dom_s1()
            } else {
                vtx.is_on_dom_s2()
            };
            if on_dom {
                vtx.arc_id(on_first)
            } else {
                resolved[iv as usize]
                    .arc
                    .and_then(|index| arcs.get(index))
                    .map(arc_token)
            }
        };
        // `arcsResolved(iv).IsNull() ? ParameterOnArc : paramsResolved(iv)`
        // (`cxx:1281-1284`).
        let vertex_param = |iv: i32| -> f64 {
            if resolved[iv as usize].arc.is_some() {
                resolved[iv as usize].param
            } else {
                wline.vertex(iv).param_on_arc(on_first)
            }
        };

        // -- `mapArcsV1` / `mapArcs` (`cxx:1201-1237`): `IndexedMap` keeps the
        //    first-seen order without duplicates, a null handle is skipped.
        let mut map_arcs_v1: Vec<usize> = Vec::new();
        for &iv in &indices_v1 {
            if let Some(token) = vertex_arc(iv) {
                if !map_arcs_v1.contains(&token) {
                    map_arcs_v1.push(token);
                }
            }
        }
        let mut map_arcs: Vec<usize> = Vec::new();
        for &iv in &indices_v2 {
            if let Some(token) = vertex_arc(iv) {
                if map_arcs_v1.contains(&token) && !map_arcs.contains(&token) {
                    map_arcs.push(token);
                }
            }
        }

        // -- for each common arc (`cxx:1239-1420`)
        for &token in &map_arcs {
            let Some(arc) = arcs.iter().find(|a| arc_token(a) == token).cloned() else {
                continue;
            };
            let arc = arc.as_ref();

            let mut iv1 = 0;
            for &iv in &indices_v1 {
                if vertex_arc(iv) == Some(token) {
                    iv1 = iv;
                    break;
                }
            }
            let mut iv2 = 0;
            for &iv in &indices_v2 {
                if vertex_arc(iv) == Some(token) {
                    iv2 = iv;
                    break;
                }
            }
            if iv1 == 0 || iv2 == 0 {
                continue;
            }

            let par1 = vertex_param(iv1);
            let par2 = vertex_param(iv2);

            // -- check that the middle point is on the arc (`cxx:1296-1310`)
            let (mut utst, mut vtst) = if on_first {
                wline.point(mid_ind).parameters_on_s1()
            } else {
                wline.point(mid_ind).parameters_on_s2()
            };
            if mid_ind == param_min_on_line {
                let (u1, v1) = if on_first {
                    wline.point(mid_ind + 1).parameters_on_s1()
                } else {
                    wline.point(mid_ind + 1).parameters_on_s2()
                };
                utst = (utst + u1) * 0.5;
                vtst = (vtst + v1) * 0.5;
            }
            if project_on_arc(utst, vtst, arc, surf, tol_arc).is_none() {
                continue;
            }

            // -- encoding of the WLine as an RLine (`cxx:1311-1320`):
            //    `new IntPatch_RLine(true, IntSurf_Unknown, IntSurf_Unknown)`
            //    then `SetArcOnS1` / `SetArcOnS2`.
            let mut rlig = RLine {
                base: IntPatchLine::with_situations(
                    true,
                    IntSurfSituation::Unknown,
                    IntSurfSituation::Unknown,
                ),
                points: Vec::new(),
                vertices: Vec::new(),
                has_first_point: false,
                has_last_point: false,
                first_index: 0,
                last_index: 0,
                arc_on_s1: on_first,
                arc_on_s2: false,
                c2d: arc.curve().cloned(),
                param_f: arc.first_parameter(),
                param_l: arc.last_parameter(),
                uv_arc: None,
                arc_edge: Some(arc.edge().clone()),
            };

            // `cxx:1329-1357`: forward or reversed sample sequence.
            let (ivmin, ivmax, parmin, parmax, reverse, ind_vmin, ind_vmax) = if par1 < par2 {
                let points = (param_min_on_line..=param_max_on_line)
                    .map(|i| *wline.point(i))
                    .collect();
                rlig.add_points(points);
                (
                    iv1,
                    iv2,
                    par1,
                    par2,
                    false,
                    indices_v1.clone(),
                    indices_v2.clone(),
                )
            } else {
                let points = (param_min_on_line..=param_max_on_line)
                    .rev()
                    .map(|i| *wline.point(i))
                    .collect();
                rlig.add_points(points);
                (
                    iv2,
                    iv1,
                    par2,
                    par1,
                    true,
                    indices_v2.clone(),
                    indices_v1.clone(),
                )
            };

            // `cxx:1358-1418`: the two end vertices, then the remaining ones
            // of each side, all carrying the arc parameter of their side.
            let undecided = TransType::Undecided;
            let mut vtx_first = wline.vertex(ivmin).clone();
            vtx_first.set_parameter(parmin);
            if resolved[ivmin as usize].arc.is_some() {
                vtx_first.set_arc(on_first, Some(token), parmin, undecided);
            }
            if reverse {
                vtx_first.reverse_transition();
            }
            rlig.add_vertex(vtx_first);

            for &iv in &ind_vmin {
                if iv != ivmin {
                    let mut vtx = wline.vertex(iv).clone();
                    vtx.set_parameter(parmin);
                    if resolved[iv as usize].arc.is_some() {
                        vtx.set_arc(on_first, Some(token), parmin, undecided);
                    }
                    if reverse {
                        vtx.reverse_transition();
                    }
                    rlig.add_vertex(vtx);
                }
            }
            for &iv in &ind_vmax {
                if iv != ivmax {
                    let mut vtx = wline.vertex(iv).clone();
                    vtx.set_parameter(parmax);
                    if resolved[iv as usize].arc.is_some() {
                        vtx.set_arc(on_first, Some(token), parmax, undecided);
                    }
                    if reverse {
                        vtx.reverse_transition();
                    }
                    rlig.add_vertex(vtx);
                }
            }
            let mut vtx_last = wline.vertex(ivmax).clone();
            vtx_last.set_parameter(parmax);
            if resolved[ivmax as usize].arc.is_some() {
                vtx_last.set_arc(on_first, Some(token), parmax, undecided);
            }
            if reverse {
                vtx_last.reverse_transition();
            }
            rlig.add_vertex(vtx_last);

            // `cxx:1416-1417`: the last-point index is the plain vertex count
            // of both sides, not the number of vertices actually appended.
            rlig.set_first_point(1);
            rlig.set_last_point(indices_v1.len() as i32 + indices_v2.len() as i32);
            slin.push(GeomIntLine::Restriction(rlig));
            is_rline = true;
        }
    }

    if is_rline
        || test_if_wline_is_restriction(slinref, &wline, s1, d1, s2, d2, tol_arc)
    {
        slin.remove(last);
    }
}
