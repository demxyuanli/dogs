//! `GeomInt_IntSS::MakeCurve`. Source: `GeomInt_IntSS_1.cxx:275`.

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::{PCONFUSION, Precision};
use occt_geom::{
    Curve, GeomCircle, GeomEllipse, GeomHyperbola, GeomLine, GeomParabola, GeomTrimmedCurve,
    Surface,
};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::fclass2d::FaceState;
use crate::geom_int::intss::IntSS;
use crate::geom_int::intss_bspline::{make_bspline, make_bspline2d};
use crate::geom_int::intss_pcurve::{
    adjust_u_periodic_curve2d, build_pcurves, treat_rline, trim_iline_on_surf_boundaries,
};
use crate::geom_int::line_tool;
use crate::geom_int::quadric::surface_parameters;
use crate::geom_int::types::{GLineKind, GeomIntLine, IntPatchIType};
use crate::geom_int::wl_approx::{ParametrizationType, WlApprox};
use crate::int_tools_wline::WLine;
use occt_geom::GeomBSplineCurve;

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

pub(crate) fn make_curve(
    iss: &mut IntSS,
    line: &GeomIntLine,
    approx: bool,
    approx_s1: bool,
    approx_s2: bool,
    tol: f64,
) {
    let mut tolpc = tol;
    let my_approx = approx;
    let my_approx1 = approx_s1;
    let my_approx2 = approx_s2;
    let my_tol_approx = 1.0e-7;
    let typl = line.arc_type();

    iss.l_construct.perform(line);
    if !iss.l_construct.is_done() || iss.l_construct.nb_parts() <= 0 {
        return;
    }

    match typl {
        IntPatchIType::Lin | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
            make_lin_parab_hypr(iss, line, typl, my_approx1, my_approx2, &mut tolpc);
        }
        IntPatchIType::Circle | IntPatchIType::Ellipse => {
            make_circle_ellipse(iss, line, typl, my_approx1, my_approx2, &mut tolpc);
        }
        IntPatchIType::Analytic => {
            if let Some(al) = line.as_aline() {
                if let Some((s1, s2)) = surface_arcs(iss) {
                    let conv = crate::intpatch::aline_to_wline::ALineToWLine::new(
                        s1.as_ref(),
                        s2.as_ref(),
                        200,
                    );
                    for wl in conv.make_wline(al) {
                        make_walking(
                            iss,
                            &wl,
                            my_approx,
                            my_approx1,
                            my_approx2,
                            my_tol_approx,
                        );
                    }
                }
            }
        }
        IntPatchIType::Walking => {
            if let Some(wl) = line.as_wline() {
                make_walking(
                    iss,
                    wl,
                    my_approx,
                    my_approx1,
                    my_approx2,
                    my_tol_approx,
                );
            }
        }
        IntPatchIType::Restriction => {
            if let Some(rl) = line.as_rline() {
                make_restriction(iss, rl, my_approx1, my_approx2);
            }
        }
    }
}

fn surface_arcs(iss: &IntSS) -> Option<(Arc<dyn Surface>, Arc<dyn Surface>)> {
    Some((iss.hs1.clone()?, iss.hs2.clone()?))
}

fn gline_curve(typl: IntPatchIType, line: &GeomIntLine) -> Option<Arc<dyn Curve>> {
    let g = line.as_gline()?;
    Some(match (typl, &g.kind) {
        (IntPatchIType::Lin, GLineKind::Lin(l)) => Arc::new(GeomLine::new(l.clone())),
        (IntPatchIType::Parabola, GLineKind::Parab(p)) => Arc::new(GeomParabola::new(p.clone())),
        (IntPatchIType::Hyperbola, GLineKind::Hypr(h)) => Arc::new(GeomHyperbola::new(h.clone())),
        (IntPatchIType::Circle, GLineKind::Circ(c)) => Arc::new(GeomCircle::new(c.clone())),
        (IntPatchIType::Ellipse, GLineKind::Elips(e)) => Arc::new(GeomEllipse::new(e.clone())),
        _ => return None,
    })
}

fn maybe_pcurve(
    iss: &mut IntSS,
    want: bool,
    fprm: f64,
    lprm: f64,
    s: &dyn Surface,
    c3: &dyn Curve,
    tolpc: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    if !want {
        return None;
    }
    let c2d = build_pcurves(fprm, lprm, tolpc, s, c3)?;
    iss.bump_tol_2d(*tolpc);
    Some(c2d)
}

fn maybe_trimmed_pcurve(
    iss: &mut IntSS,
    want: bool,
    fprm: f64,
    lprm: f64,
    s: &dyn Surface,
    c3: &dyn Curve,
    tolpc: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    let c2d = maybe_pcurve(iss, want, fprm, lprm, s, c3, tolpc)?;
    Some(Arc::new(Geom2dTrimmedCurve::new(c2d, fprm, lprm)))
}

fn make_lin_parab_hypr(
    iss: &mut IntSS,
    line: &GeomIntLine,
    typl: IntPatchIType,
    my_approx1: bool,
    my_approx2: bool,
    tolpc: &mut f64,
) {
    let Some(newc) = gline_curve(typl, line) else {
        return;
    };
    let Some((s1, s2)) = surface_arcs(iss) else {
        return;
    };
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    let a_nb_parts = iss.l_construct.nb_parts();
    for i in 1..=a_nb_parts {
        let Some((fprm, lprm)) = iss.l_construct.part(i) else {
            continue;
        };
        if !Precision::is_negative_infinite(fprm) && !Precision::is_positive_infinite(lprm) {
            let a_ct3d: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(newc.clone(), fprm, lprm));
            let c1 = maybe_trimmed_pcurve(iss, my_approx1, fprm, lprm, s1, newc.as_ref(), tolpc);
            let c2 = maybe_trimmed_pcurve(iss, my_approx2, fprm, lprm, s2, newc.as_ref(), tolpc);
            iss.append_line(a_ct3d, c1, c2);
            continue;
        }
        let typ_s1 = classify_surface(s1);
        let typ_s2 = classify_surface(s2);
        if matches!(typ_s1, SurfaceKind::Other) || matches!(typ_s2, SurfaceKind::Other) {
            iss.append_line(newc.clone(), None, None);
            continue;
        }
        let b_fnit = Precision::is_negative_infinite(fprm);
        let b_lpit = Precision::is_positive_infinite(lprm);
        let d_t = 100.0;
        let a_test_prm = if b_fnit && !b_lpit {
            lprm - d_t
        } else if !b_fnit && b_lpit {
            fprm + d_t
        } else {
            0.0
        };
        let ptref = newc.d0(a_test_prm);
        if classify_both(iss, &ptref) {
            iss.append_line(newc.clone(), None, None);
        }
    }
}

fn classify_both(iss: &IntSS, ptref: &GpPnt) -> bool {
    let Some((s1, s2)) = surface_arcs(iss) else {
        return false;
    };
    let Some((d1, d2)) = iss.l_construct.domain_clones() else {
        return false;
    };
    let tol_x = Precision::CONFUSION;
    let Some((u1, v1)) = surface_parameters(s1.as_ref(), ptref) else {
        return false;
    };
    let Some((u2, v2)) = surface_parameters(s2.as_ref(), ptref) else {
        return false;
    };
    d1.classify(GpPnt2d::new(u1, v1), tol_x) != FaceState::Out
        && d2.classify(GpPnt2d::new(u2, v2), tol_x) != FaceState::Out
}

fn make_circle_ellipse(
    iss: &mut IntSS,
    line: &GeomIntLine,
    typl: IntPatchIType,
    my_approx1: bool,
    my_approx2: bool,
    tolpc: &mut f64,
) {
    let Some(newc) = gline_curve(typl, line) else {
        return;
    };
    let Some((s1, s2)) = surface_arcs(iss) else {
        return;
    };
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    let a_real_eps = f64::EPSILON;
    let a_period = TWO_PI;
    let a_nb_parts = iss.l_construct.nb_parts();
    for i in 1..=a_nb_parts {
        let Some((fprm, lprm)) = iss.l_construct.part(i) else {
            continue;
        };
        if fprm.abs() > a_real_eps || (lprm - a_period).abs() > a_real_eps {
            append_trimmed_conic(iss, newc.clone(), fprm, lprm, s1, s2, my_approx1, my_approx2, tolpc);
            continue;
        }
        if a_nb_parts == 1 && fprm.abs() < a_real_eps && (lprm - TWO_PI).abs() < a_real_eps {
            append_trimmed_conic(iss, newc.clone(), fprm, lprm, s1, s2, my_approx1, my_approx2, tolpc);
            break;
        }
        let a_two_pi_div_17 = TWO_PI / 17.0;
        for j in 0..=17 {
            let ptref = newc.d0(j as f64 * a_two_pi_div_17);
            if classify_both(iss, &ptref) {
                let c1 = maybe_pcurve(iss, my_approx1, fprm, lprm, s1, newc.as_ref(), tolpc);
                let c2 = maybe_pcurve(iss, my_approx2, fprm, lprm, s2, newc.as_ref(), tolpc);
                iss.append_line(newc.clone(), c1, c2);
                break;
            }
        }
    }
}

fn append_trimmed_conic(
    iss: &mut IntSS,
    newc: Arc<dyn Curve>,
    fprm: f64,
    lprm: f64,
    s1: &dyn Surface,
    s2: &dyn Surface,
    my_approx1: bool,
    my_approx2: bool,
    tolpc: &mut f64,
) {
    let a_tc3d = GeomTrimmedCurve::new(newc.clone(), fprm, lprm);
    let f = a_tc3d.first_parameter();
    let l = a_tc3d.last_parameter();
    let c1 = maybe_pcurve(iss, my_approx1, fprm, lprm, s1, newc.as_ref(), tolpc);
    let c2 = maybe_pcurve(iss, my_approx2, fprm, lprm, s2, newc.as_ref(), tolpc);
    let _ = (f, l);
    iss.append_line(Arc::new(a_tc3d), c1, c2);
}

fn make_walking(
    iss: &mut IntSS,
    wl: &WLine,
    my_approx: bool,
    my_approx1: bool,
    my_approx2: bool,
    my_tol_approx: f64,
) {
    let Some((s1, s2)) = surface_arcs(iss) else {
        return;
    };
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    if !my_approx {
        let a_nb_parts = iss.l_construct.nb_parts();
        for i in 1..=a_nb_parts {
            let Some((fprm, lprm)) = iss.l_construct.part(i) else {
                continue;
            };
            let ifprm = fprm as i32;
            let ilprm = lprm as i32;
            let a_h1 = if my_approx1 {
                make_bspline2d(wl, ifprm, ilprm, true)
            } else {
                None
            };
            let a_h2 = if my_approx2 {
                make_bspline2d(wl, ifprm, ilprm, false)
            } else {
                None
            };
            if let Some(a_bsp) = make_bspline(wl, ifprm, ilprm) {
                iss.append_line(a_bsp, a_h1, a_h2);
            }
        }
        return;
    }

    let mut theapp3d = WlApprox::new();
    let tol2d = my_tol_approx;
    let a_tol_ss = 2.0e-7;
    theapp3d.set_parameters(
        my_tol_approx,
        tol2d,
        4,
        8,
        0,
        30,
        !iss.same_surfaces,
        ParametrizationType::ChordLength,
    );
    let mut a_seq_of_l: Vec<WLine> = Vec::new();
    let b_is_decomposited = line_tool::decomposition_of_wline(
        wl,
        s1,
        s2,
        a_tol_ss,
        &iss.l_construct,
        &mut a_seq_of_l,
    );
    let a_nb_parts = iss.l_construct.nb_parts();
    let a_nb_seq_of_l = a_seq_of_l.len() as i32;
    let nbiter = if b_is_decomposited {
        a_nb_seq_of_l
    } else {
        a_nb_parts
    };
    for i in 1..=nbiter {
        let (piece, ifprm, ilprm) = if b_is_decomposited {
            let w = &a_seq_of_l[(i as usize).saturating_sub(1)];
            (w, 1, w.nb_pnts())
        } else {
            let Some((fprm, lprm)) = iss.l_construct.part(i) else {
                continue;
            };
            (wl, fprm as i32, lprm as i32)
        };
        walking_approx_one(
            iss,
            piece,
            ifprm,
            ilprm,
            my_approx,
            my_approx1,
            my_approx2,
            my_tol_approx,
            &mut theapp3d,
        );
    }
}

fn walking_approx_one(
    iss: &mut IntSS,
    wl: &WLine,
    ifprm: i32,
    ilprm: i32,
    my_approx: bool,
    my_approx1: bool,
    my_approx2: bool,
    my_tol_approx: f64,
    theapp3d: &mut WlApprox,
) {
    let Some((s1, s2)) = surface_arcs(iss) else {
        return;
    };
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    let mut an_approx = my_approx;
    let mut an_approx1 = my_approx1;
    let mut an_approx2 = my_approx2;
    let typs1 = classify_surface(s1);
    let typs2 = classify_surface(s2);
    if typs1 == SurfaceKind::Plane {
        an_approx = false;
        an_approx1 = true;
    } else if typs2 == SurfaceKind::Plane {
        an_approx = false;
        an_approx2 = true;
    }
    let _ = (an_approx, an_approx1, an_approx2);
    let tol2d = my_tol_approx;
    theapp3d.set_parameters(
        my_tol_approx,
        tol2d,
        4,
        8,
        0,
        30,
        !iss.same_surfaces,
        ParametrizationType::ChordLength,
    );
    if typs1 == SurfaceKind::Plane {
        theapp3d.perform(s1, s2, wl, false, true, my_approx2, ifprm, ilprm);
    } else if typs2 == SurfaceKind::Plane {
        theapp3d.perform(s1, s2, wl, false, my_approx1, true, ifprm, ilprm);
    } else {
        theapp3d.perform(s1, s2, wl, true, my_approx1, my_approx2, ifprm, ilprm);
    }

    if !theapp3d.is_done() {
        let a_h1 = if my_approx1 {
            make_bspline2d(wl, ifprm, ilprm, true)
        } else {
            None
        };
        let a_h2 = if my_approx2 {
            make_bspline2d(wl, ifprm, ilprm, false)
        } else {
            None
        };
        if let Some(a_bsp) = make_bspline(wl, ifprm, ilprm) {
            iss.append_line(a_bsp, a_h1, a_h2);
        }
        return;
    }

    if my_approx1 || my_approx2 || typs1 == SurfaceKind::Plane || typs2 == SurfaceKind::Plane {
        iss.bump_tol_2d(theapp3d.tol_reached_2d());
    }
    if typs1 == SurfaceKind::Plane || typs2 == SurfaceKind::Plane {
        iss.tol_reached_3d = iss.tol_reached_2d;
    } else {
        iss.bump_tol_3d(theapp3d.tol_reached_3d());
    }

    let a_nb_multi = theapp3d.nb_multi_curves();
    for j in 1..=a_nb_multi {
        let Some(mbspc) = theapp3d.value(j).cloned() else {
            continue;
        };
        if typs1 == SurfaceKind::Plane {
            append_plane_lift(iss, s1, s2, &mbspc, true, my_approx1, my_approx2);
        } else if typs2 == SurfaceKind::Plane {
            append_plane_lift(iss, s1, s2, &mbspc, false, my_approx1, my_approx2);
        } else {
            append_walking_done(iss, s1, s2, &mbspc, my_approx1, my_approx2);
        }
    }
}

fn bspline_from_poles(poles: Vec<occt_core::gp::GpPnt>, knots: &[f64]) -> Option<Arc<dyn Curve>> {
    if poles.len() < 2 {
        return None;
    }
    let k = if knots.len() == poles.len() + 2 {
        knots.to_vec()
    } else {
        degree1_knots(poles.len())
    };
    GeomBSplineCurve::new(poles, k, 1)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve>)
}

fn bspline2d_from_poles(
    pts: &[occt_core::gp::GpPnt2d],
    knots: &[f64],
) -> Option<Arc<dyn Curve2d>> {
    if pts.len() < 2 {
        return None;
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y()).collect();
    let k = if knots.len() == pts.len() + 2 {
        knots.to_vec()
    } else {
        degree1_knots(pts.len())
    };
    occt_geom2d::Geom2dBSplineCurve::new(xs, ys, k, 1)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve2d>)
}

fn degree1_knots(n: usize) -> Vec<f64> {
    let mut k = Vec::with_capacity(n + 2);
    k.push(0.0);
    for i in 0..n {
        k.push(i as f64);
    }
    k.push((n - 1) as f64);
    k
}

fn append_plane_lift(
    iss: &mut IntSS,
    s1: &dyn Surface,
    s2: &dyn Surface,
    mbspc: &crate::geom_int::wl_approx::MultiBSpCurve,
    plane_is_s1: bool,
    my_approx1: bool,
    my_approx2: bool,
) {
    let (uvs, plane) = if plane_is_s1 {
        (&mbspc.poles2d_s1, s1)
    } else {
        (&mbspc.poles2d_s2, s2)
    };
    if uvs.is_empty() {
        return;
    }
    let mut poles = Vec::with_capacity(uvs.len());
    for uv in uvs {
        poles.push(plane.d0(uv.x(), uv.y()));
    }
    if poles.len() == 2 && poles[0].square_distance(&poles[poles.len() - 1]) < 2.0 * f64::EPSILON {
        return;
    }
    let Some(bs) = bspline_from_poles(poles, &mbspc.knots) else {
        return;
    };
    let mut c1 = if my_approx1 {
        if plane_is_s1 {
            bspline2d_from_poles(&mbspc.poles2d_s1, &mbspc.knots)
        } else {
            mbspc.curve2d_s1.clone()
        }
    } else {
        None
    };
    let mut c2 = if my_approx2 {
        if plane_is_s1 {
            mbspc.curve2d_s2.clone()
        } else {
            bspline2d_from_poles(&mbspc.poles2d_s2, &mbspc.knots)
        }
    } else {
        None
    };
    if let Some(c) = c1.as_mut() {
        adjust_u_periodic_curve2d(s1, c);
    }
    if let Some(c) = c2.as_mut() {
        adjust_u_periodic_curve2d(s2, c);
    }
    iss.append_line(bs, c1, c2);
}

fn append_walking_done(
    iss: &mut IntSS,
    s1: &dyn Surface,
    s2: &dyn Surface,
    mbspc: &crate::geom_int::wl_approx::MultiBSpCurve,
    my_approx1: bool,
    my_approx2: bool,
) {
    let Some(mut bs) = mbspc.curve3d.clone() else {
        return;
    };
    if mbspc.poles3d.len() == 2 {
        let a = mbspc.poles3d[0];
        let b = mbspc.poles3d[1];
        let a_dist = a.coord.square_modulus().max(b.coord.square_modulus());
        let eps = f64::EPSILON * a_dist.max(1.0);
        if a.square_distance(&b) < 2.0 * eps {
            return;
        }
    }
    let _ = &mut bs;
    let mut c1 = if my_approx1 {
        mbspc.curve2d_s1.clone()
    } else {
        None
    };
    let mut c2 = if my_approx2 {
        mbspc.curve2d_s2.clone()
    } else {
        None
    };
    if let Some(c) = c1.as_mut() {
        adjust_u_periodic_curve2d(s1, c);
    }
    if let Some(c) = c2.as_mut() {
        adjust_u_periodic_curve2d(s2, c);
    }
    iss.append_line(bs, c1, c2);
}

fn make_restriction(
    iss: &mut IntSS,
    rl: &crate::geom_int::types::RLine,
    my_approx1: bool,
    my_approx2: bool,
) {
    let Some((s1, s2)) = surface_arcs(iss) else {
        return;
    };
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    let (a_c3d, a_c2d1, a_c2d2, a_tol_reached) = treat_rline(rl, s1, s2);
    let Some(a_c3d) = a_c3d else {
        return;
    };
    iss.bump_tol_3d(a_tol_reached);
    let mut box1 = BndBox2d::new();
    let mut box2 = BndBox2d::new();
    let (u1f, u1l) = s1.u_range();
    let (v1f, v1l) = s1.v_range();
    let (u2f, u2l) = s2.u_range();
    let (v2f, v2l) = s2.v_range();
    box1.update_point(u1f, v1f);
    box1.update_point(u1l, v1l);
    box2.update_point(u2f, v2f);
    box2.update_point(u2l, v2l);
    let mut an_array = vec![a_c3d.first_parameter(), a_c3d.last_parameter()];
    trim_iline_on_surf_boundaries(
        a_c2d1.as_deref(),
        a_c2d2.as_deref(),
        &box1,
        &box2,
        &mut an_array,
    );
    let n = an_array.len().saturating_sub(1);
    for an_ind in 0..n {
        let a_par_f = an_array[an_ind];
        let a_par_l = an_array[an_ind + 1];
        if a_par_l - a_par_f <= PCONFUSION {
            continue;
        }
        let a_par = 0.5 * (a_par_f + a_par_l);
        let mut a_curv2d1 = None;
        let mut a_curv2d2 = None;
        if let Some(c) = &a_c2d1 {
            let a_pt = c.d0(a_par);
            if box1.is_out(&a_pt) {
                continue;
            }
            if my_approx1 {
                a_curv2d1 = Some(Arc::new(Geom2dTrimmedCurve::new(c.clone(), a_par_f, a_par_l))
                    as Arc<dyn Curve2d>);
            }
        }
        if let Some(c) = &a_c2d2 {
            let a_pt = c.d0(a_par);
            if box2.is_out(&a_pt) {
                continue;
            }
            if my_approx2 {
                a_curv2d2 = Some(Arc::new(Geom2dTrimmedCurve::new(c.clone(), a_par_f, a_par_l))
                    as Arc<dyn Curve2d>);
            }
        }
        let a_curv3d: Arc<dyn Curve> =
            Arc::new(GeomTrimmedCurve::new(a_c3d.clone(), a_par_f, a_par_l));
        iss.append_line(a_curv3d, a_curv2d1, a_curv2d2);
    }
}
