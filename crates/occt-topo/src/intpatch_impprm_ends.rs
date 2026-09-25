//! ImpPrm walking-line end vertices: Destination / Recadre / MakeTransition.
//! Source: `IntPatch_ImpPrmIntersection.cxx` Recadre ~473, Perform ~1080–1290;
//! `IntSurf.cxx` MakeTransition.

use occt_core::gp::{GpDir, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::Surface;

use crate::brep_surface::SurfaceKind;
use crate::int_tools_wline::{PatchPoint, TransType, WLine};
use crate::intpatch::impimp::{ImplicitQuad, PathPoint};

use super::iwalking::{IwLine, WalkStart};

/// `IntSurf::MakeTransition` — TLine / TArc from mixed product with the surface normal.
pub(crate) fn make_transition(tg_first: &GpVec, tg_second: &GpVec, normale: &GpDir) -> (TransType, TransType) {
    let pvect = tg_second.crossed(tg_first);
    let n_second = tg_second.magnitude();
    let n_first = tg_first.magnitude();
    let n_ang = n_second * n_first * ANGULAR;
    if n_first <= CONFUSION {
        (TransType::Undecided, TransType::Undecided)
    } else if n_second <= CONFUSION || pvect.magnitude() <= n_ang {
        (TransType::Touch, TransType::Touch)
    } else {
        let mut yu = pvect.dot(&GpVec::from_xyz(normale.xyz()));
        yu /= n_second * n_first;
        if yu > 0.0001 {
            (TransType::In, TransType::Out)
        } else if yu < -0.0001 {
            (TransType::Out, TransType::In)
        } else {
            (TransType::Undecided, TransType::Undecided)
        }
    }
}

/// `Recadre` — wrap periodic UV to the neighbour walking sample.
pub(crate) fn recadre(
    kind1: SurfaceKind,
    kind2: SurfaceKind,
    pt: &mut PatchPoint,
    neighbour: &crate::int_tools_wline::PntOn2S,
    mut u1: f64,
    mut v1: f64,
    mut u2: f64,
    mut v2: f64,
) {
    let (u1p, v1p) = neighbour.parameters_on_s1();
    let (u2p, v2p) = neighbour.parameters_on_s2();
    wrap_surface(kind1, &mut u1, &mut v1, u1p, v1p);
    wrap_surface(kind2, &mut u2, &mut v2, u2p, v2p);
    pt.u1 = u1;
    pt.v1 = v1;
    pt.u2 = u2;
    pt.v2 = v2;
}

fn wrap_surface(kind: SurfaceKind, u: &mut f64, v: &mut f64, up: f64, vp: f64) {
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;
    if kind == SurfaceKind::Torus {
        while *v < vp - 1.5 * std::f64::consts::PI {
            *v += two_pi;
        }
        while *v > vp + 1.5 * std::f64::consts::PI {
            *v -= two_pi;
        }
    }
    match kind {
        SurfaceKind::Torus | SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere => {
            while *u < up - 1.5 * std::f64::consts::PI {
                *u += two_pi;
            }
            while *u > up + 1.5 * std::f64::consts::PI {
                *u -= two_pi;
            }
        }
        _ => {}
    }
}

/// Bind first/last vertices (`ImpPrmIntersection.cxx` ~1080–1290).
pub(crate) fn attach_wline_ends(
    w: &mut WLine,
    iw: &IwLine,
    seqpdep: &[WalkStart],
    rst: &[PathPoint],
    dest: &[i32],
    reversed: bool,
    kind1: SurfaceKind,
    kind2: SurfaceKind,
    quad: &ImplicitQuad,
    prm: &dyn Surface,
    quad_surf: &dyn Surface,
    tol_arc: f64,
) {
    w.vertices.clear();
    w.has_first_point = false;
    w.has_last_point = false;
    let nbpts = w.points.len();
    if nbpts < 2 {
        return;
    }
    let (vmin, vmax) = quad_surf.v_range();
    const TOL_V: f64 = 1.0e-14;

    if iw.has_first && !iw.tgt_begin {
        let indfirst = iw.first_index;
        if indfirst >= 1 {
            let idx = (indfirst as usize).saturating_sub(1);
            if idx < seqpdep.len() {
                attach_one_end(
                    w,
                    rst,
                    dest,
                    seqpdep[idx],
                    indfirst,
                    reversed,
                    kind1,
                    kind2,
                    quad,
                    prm,
                    vmin,
                    vmax,
                    TOL_V,
                    tol_arc,
                    1.0,
                    1,
                    true,
                );
            }
        }
    } else if iw.tgt_begin {
        push_sample_vertex(w, 0, 1.0, true, true);
    } else {
        push_sample_vertex(w, 0, 1.0, false, true);
    }

    if iw.has_last && !iw.tgt_end {
        let indlast = iw.last_index;
        if indlast >= 1 {
            let idx = (indlast as usize).saturating_sub(1);
            if idx < seqpdep.len() {
                let neigh = if nbpts > 1 { nbpts - 1 } else { 1 };
                attach_one_end(
                    w,
                    rst,
                    dest,
                    seqpdep[idx],
                    indlast,
                    reversed,
                    kind1,
                    kind2,
                    quad,
                    prm,
                    vmin,
                    vmax,
                    TOL_V,
                    tol_arc,
                    nbpts as f64,
                    neigh,
                    false,
                );
            }
        }
    } else if iw.tgt_end {
        push_sample_vertex(w, nbpts - 1, nbpts as f64, true, false);
    } else {
        push_sample_vertex(w, nbpts - 1, nbpts as f64, false, false);
    }
}

fn push_sample_vertex(w: &mut WLine, i: usize, param: f64, tangent: bool, first: bool) {
    let a = w.points[i];
    let pt = PatchPoint::new(a.p, param, a.u1, a.v1, a.u2, a.v2);
    let _ = tangent;
    w.vertices.push(pt);
    if first {
        w.has_first_point = true;
    } else {
        w.has_last_point = true;
    }
}

fn attach_one_end(
    w: &mut WLine,
    rst: &[PathPoint],
    dest: &[i32],
    ppoint: WalkStart,
    ind: i32,
    reversed: bool,
    kind1: SurfaceKind,
    kind2: SurfaceKind,
    quad: &ImplicitQuad,
    prm: &dyn Surface,
    vmin: f64,
    vmax: f64,
    tol_v: f64,
    tol_arc: f64,
    param_on_line: f64,
    neighbour_1based: usize,
    first: bool,
) {
    let tgline = if first {
        ppoint.d3d
    } else {
        ppoint.d3d.reversed()
    };
    let _ = tgline;
    // `Multiplicity() == sequv.Length()-1`; UV-box PathPoints are IsNew (one UV).
    let mut themult: i32 = 0;
    for i in (0..rst.len()).rev() {
        if dest.get(i).copied().unwrap_or(0) != ind {
            continue;
        }
        let rp = rst[i];
        let (mut u1, mut v1, mut u2, mut v2, d1u, d1v) = if !reversed {
            let (uq, mut vq) = quad.parameters(&ppoint.p);
            clamp_v(&mut vq, vmin, vmax, tol_v);
            let (up, vp) = (ppoint.u, ppoint.v);
            let (_, du, dv) = prm.d1(up, vp);
            (uq, vq, up, vp, du, dv)
        } else {
            let (uq, mut vq) = quad.parameters(&ppoint.p);
            clamp_v(&mut vq, vmin, vmax, tol_v);
            let (up, vp) = (ppoint.u, ppoint.v);
            let (_, du, dv) = prm.d1(up, vp);
            (up, vp, uq, vq, du, dv)
        };
        let vec_n = d1u.crossed(&d1v);
        let mut pt = PatchPoint::new(ppoint.p, param_on_line, u1, v1, u2, v2);
        pt.tolerance = tol_arc;
        let neigh_i = neighbour_1based.saturating_sub(1).min(w.points.len().saturating_sub(1));
        recadre(
            kind1,
            kind2,
            &mut pt,
            &w.points[neigh_i],
            u1,
            v1,
            u2,
            v2,
        );
        u1 = pt.u1;
        v1 = pt.v1;
        u2 = pt.u2;
        v2 = pt.v2;
        let _ = (u1, v1, u2, v2);
        let (_p2d, d2d) = rp.arc.d1(rp.param_on_arc);
        let tgrst = d1u
            .multiplied_scalar(d2d.x())
            .added(&d1v.multiplied_scalar(d2d.y()));
        let (t_line, t_arc) = if vec_n.square_magnitude() > 1e-13 {
            match GpDir::from_vec(&vec_n) {
                Ok(dir) => make_transition(&tgline, &tgrst, &dir),
                Err(_) => (TransType::Undecided, TransType::Undecided),
            }
        } else {
            (TransType::Undecided, TransType::Undecided)
        };
        if reversed {
            pt.on_dom_s1 = true;
            pt.trans1 = t_arc;
            pt.trans2 = t_line;
        } else {
            pt.on_dom_s2 = true;
            pt.trans2 = t_arc;
            pt.trans1 = t_line;
        }
        w.vertices.push(pt);
        if themult == 0 {
            if first {
                w.has_first_point = true;
            } else {
                w.has_last_point = true;
            }
        }
        themult -= 1;
    }
}

fn clamp_v(v: &mut f64, vmin: f64, vmax: f64, tol_v: f64) {
    if vmin.is_finite() && *v < vmin && vmin - *v < tol_v {
        *v = vmin;
    }
    if vmax.is_finite() && *v > vmax && *v - vmax < tol_v {
        *v = vmax;
    }
}
