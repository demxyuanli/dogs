//! `IntPatch_LineConstructor::Perform` (`cxx:1433-2085`): cut one raw
//! intersection line at the `VPoint`s that fall on the two faces' domains and
//! collect the surviving pieces in `slin`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_LineConstructor.cxx`
//! - Analytic arm (`cxx:1451-1514`)
//! - Walking arm (`cxx:1515-1689`)
//! - Geometric arm (`cxx:1690-1867`)
//! - Restriction arm (`cxx:1868-2084`)
//!
//! Translation notes:
//! - `Tol` (`cxx:1447`) starts at `Precision::PConfusion() * 100`, is left
//!   untouched by the Analytic and Walking arms, and is overwritten by the
//!   Geometric and Restriction arms. There it survives across loop iterations
//!   because OCCT assigns to the function-scope variable; the port keeps one
//!   mutable `tol` per arm to preserve that.
//! - `Parameters(HS1, HS2, P, U1, V1, U2, V2)` (`cxx:156-210`) raises
//!   `Standard_ConstructionError` when a surface is not a quadric
//!   (`cxx:184`, `:204`). The port's [`super::helpers::parameters`] answers
//!   `Err`, and the arms skip the segment, leaving `TabIndex` at 0 exactly like
//!   a segment that never entered `AddLine`.
//! - `myDom1->Classify(..., Tol, false)` / `IsThePointOn(..., Tol)` mutate the
//!   lazy `myFClass2d` (`BRepTopAdaptor_TopolTool.cxx:176-202`), so the domains
//!   are taken by `&mut`.
//! - The Geometric arm's `myDom1->Init(); if (myDom2->More())` (`cxx:1742-1746`)
//!   checks the *second* iterator's state and then reads the *first*'s arc.
//!   That is what OCCT runs, so the port keeps the two calls on `d2` / `d1`
//!   rather than "fixing" the pair.
//! - `TabIndex` (`cxx:1459`) is `new int[nbvtx + 2]` with 1-based access; the
//!   port allocates `nbvtx + 2` zeroes and indexes with the same 1-based index.

use occt_core::elib::clib;
use occt_core::gp::{GpPnt2d, GpVec};
use occt_core::precision::PCONFUSION;
use occt_geom::Surface;

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::fclass2d::FaceState;
use crate::geom_int::{GLineKind, GeomIntLine, IntPatchIType};

use super::add_line::add_line;
use super::helpers::{
    compute_parametric_tolerance, is_segment_small, local_last_parameter,
    local_first_parameter, parameters, recadre,
};
use super::wline_to_rline::test_wline_to_rline;
use super::LineConstructor;

/// `Perform(SL, L, S1, D1, S2, D2, TolArc)` (`cxx:1433-2085`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform(
    this: &mut LineConstructor,
    sl: &[GeomIntLine],
    l: &GeomIntLine,
    s1: &dyn Surface,
    d1: &mut TopolTool,
    s2: &dyn Surface,
    d2: &mut TopolTool,
    tol_arc: f64,
) {
    // `double Tol = Precision::PConfusion()*100.;` (`cxx:1447`).
    let tol = PCONFUSION * 100.0;
    let typl = l.arc_type();
    if typl == IntPatchIType::Analytic {
        perform_analytic(this, l, s1, d1, s2, d2, tol);
        return;
    }
    if typl == IntPatchIType::Walking {
        perform_walking(this, sl, l, s1, d1, s2, d2, tol, tol_arc);
        return;
    }
    if typl != IntPatchIType::Restriction {
        perform_geometric(this, l, s1, d1, s2, d2, tol);
        return;
    }
    perform_restriction(this, l, s1, d1, s2, d2, tol);
}

/// `Parameters` + `Recadre` of a 3D point (`cxx:1474-1477` and friends).
fn uv_on_both(
    s1: &dyn Surface,
    s2: &dyn Surface,
    p: &occt_core::gp::GpPnt,
) -> Option<(f64, f64, f64, f64)> {
    let (u1, v1, u2, v2) = parameters(s1, s2, p).ok()?;
    let mut uv = (u1, v1, u2, v2);
    recadre(s1, s2, &mut uv);
    Some(uv)
}

/// `if (in1 == TopAbs_OUT || in2 == TopAbs_OUT) {} else { ... }` — the two
/// `TopAbs_UNKNOWN` reclassifications of `cxx:1758-1766` and `:1816-1824`
/// belong to the Geometric arm only.
fn both_not_out(a: FaceState, b: FaceState) -> bool {
    a != FaceState::Out && b != FaceState::Out
}

/// Analytic arm (`cxx:1451-1514`).
fn perform_analytic(
    this: &mut LineConstructor,
    l: &GeomIntLine,
    s1: &dyn Surface,
    d1: &mut TopolTool,
    s2: &dyn Surface,
    d2: &mut TopolTool,
    tol: f64,
) {
    this.slin_mut().clear();
    let Some(aline) = l.as_aline() else {
        return;
    };
    let nbvtx = aline.nb_vertex();
    let mut tab_index = vec![0i32; (nbvtx + 2) as usize];
    let mut numline = 0;

    for i in 1..nbvtx {
        let firstp = aline.vertex(i).parameter_on_line();
        let lastp = aline.vertex(i + 1).parameter_on_line();
        if firstp == lastp {
            continue;
        }
        let pmid = (firstp + lastp) * 0.5;
        let p_mid = aline.value(pmid);
        let Some((u1, v1, u2, v2)) = uv_on_both(s1, s2, &p_mid) else {
            continue;
        };
        let in1 = d1.classify(GpPnt2d::new(u1, v1), tol, false);
        let in2 = if in1 != FaceState::Out {
            d2.classify(GpPnt2d::new(u2, v2), tol, false)
        } else {
            FaceState::Out
        };
        if both_not_out(in1, in2) {
            numline += 1;
            tab_index[i as usize] = numline;
            tab_index[(i + 1) as usize] = numline;
            add_line(l, i, i + 1, &mut tab_index, this.slin_mut());
        }
    }

    // `cxx:1496-1509`: external edge-edge vertices, which are restrictions on
    // both edges and so belong to no line yet.
    for i in 1..=nbvtx {
        if tab_index[i as usize] != 0 {
            continue;
        }
        let vtx = aline.vertex(i);
        if vtx.is_on_dom_s1() && vtx.is_on_dom_s2() {
            numline += 1;
            tab_index[i as usize] = numline;
            add_line(l, i, i, &mut tab_index, this.slin_mut());
        }
    }
}

/// Walking arm (`cxx:1515-1689`).
#[allow(clippy::too_many_arguments)]
fn perform_walking(
    this: &mut LineConstructor,
    sl: &[GeomIntLine],
    l: &GeomIntLine,
    s1: &dyn Surface,
    d1: &mut TopolTool,
    s2: &dyn Surface,
    d2: &mut TopolTool,
    tol: f64,
    tol_arc: f64,
) {
    this.slin_mut().clear();
    let Some(wline) = l.as_wline() else {
        return;
    };
    let nbvtx = wline.nb_vertex();
    let mut tab_index = vec![0i32; (nbvtx + 2) as usize];
    let mut numline = 0;

    for i in 1..nbvtx {
        let firstp = wline.vertex(i).parameter_on_line();
        let lastp = wline.vertex(i + 1).parameter_on_line();
        if firstp == lastp || is_segment_small(wline, i, i + 1) {
            continue;
        }
        let int_lastp = lastp as i32;
        let int_firstp = firstp as i32;
        let mut pmid = ((firstp + lastp) / 2.0) as i32;
        if pmid == int_lastp {
            pmid = int_firstp;
        }
        let (mut u1, mut v1, mut u2, mut v2) = wline.point(pmid).parameters();
        let mut uv = (u1, v1, u2, v2);
        recadre(s1, s2, &mut uv);
        u1 = uv.0;
        v1 = uv.1;
        u2 = uv.2;
        v2 = uv.3;

        // `cxx:1549-1570`: the tolerance comes from the surface derivative,
        // not from `Tol`.
        let (_, ad1u, ad1v) = s1.d1(u1, v1);
        let mut a_tolerance = compute_parametric_tolerance(tol_arc, &ad1u, &ad1v);
        let mut in1 = d1.classify(GpPnt2d::new(u1, v1), a_tolerance, false);
        let mut in2 = FaceState::Out;
        if in1 != FaceState::Out {
            let (_, ad1u, ad1v) = s2.d1(u2, v2);
            a_tolerance = compute_parametric_tolerance(tol_arc, &ad1u, &ad1v);
            in2 = d2.classify(GpPnt2d::new(u2, v2), a_tolerance, false);
        }

        // `cxx:1579-1630`: `(lastp - firstp) == 1` can mean a bad
        // parametrisation; classify the "virtual" point between the two
        // vertices instead.
        if (int_lastp - int_firstp).abs() == 1 {
            let mut vf = wline.vertex(i).parameters();
            recadre(s1, s2, &mut vf);
            let mut vl = wline.vertex(i + 1).parameters();
            recadre(s1, s2, &mut vl);
            let (vfu1, vfv1, _, _) = vf;
            let (vlu1, vlv1, _, _) = vl;
            if in1 != FaceState::In {
                let pv_f = GpPnt2d::new(vfu1, vfv1);
                let pv_l = GpPnt2d::new(vlu1, vlv1);
                let p_pm = GpPnt2d::new(u1, v1);
                let (du, dv) = if pv_f.distance(&p_pm) > pv_l.distance(&p_pm) {
                    ((vfu1 + u1) * 0.5, (vfv1 + v1) * 0.5)
                } else {
                    ((vlu1 + u1) * 0.5, (vlv1 + v1) * 0.5)
                };
                let (_, ad1u, ad1v) = s1.d1(du, dv);
                a_tolerance = compute_parametric_tolerance(tol_arc, &ad1u, &ad1v);
                in1 = d1.classify(GpPnt2d::new(du, dv), a_tolerance, false);
            }
            if in2 != FaceState::In {
                let (vfu2, vfv2) = (vf.2, vf.3);
                let (vlu2, vlv2) = (vl.2, vl.3);
                let pv_f = GpPnt2d::new(vfu2, vfv2);
                let pv_l = GpPnt2d::new(vlu2, vlv2);
                let p_pm = GpPnt2d::new(u2, v2);
                let (du, dv) = if pv_f.distance(&p_pm) > pv_l.distance(&p_pm) {
                    ((vfu2 + u2) * 0.5, (vfv2 + v2) * 0.5)
                } else {
                    ((vlu2 + u2) * 0.5, (vlv2 + v2) * 0.5)
                };
                let (_, ad1u, ad1v) = s2.d1(du, dv);
                a_tolerance = compute_parametric_tolerance(tol_arc, &ad1u, &ad1v);
                in2 = d2.classify(GpPnt2d::new(du, dv), a_tolerance, false);
            }
        }

        if !both_not_out(in1, in2) {
            continue;
        }

        // `cxx:1635-1656`: two samples whose UV on the first surface coincide
        // mean the segment is far too small to keep.
        let mut uv_a = wline.point(firstp as i32).parameters();
        recadre(s1, s2, &mut uv_a);
        let mut uv_b = wline.point(lastp as i32).parameters();
        recadre(s1, s2, &mut uv_b);
        let mut dd12_u = (uv_a.0 - uv_b.0).abs();
        let mut dd12_v = (uv_a.1 - uv_b.1).abs();
        let mut lignetrop_petite = false;
        if dd12_u + dd12_v < 1e-12 {
            dd12_u = (u1 - uv_b.0).abs();
            dd12_v = (v1 - uv_b.1).abs();
            if dd12_u + dd12_v < 1e-12 {
                lignetrop_petite = true;
            }
        }
        if lignetrop_petite {
            continue;
        }

        numline += 1;
        tab_index[i as usize] = numline;
        tab_index[(i + 1) as usize] = numline;
        add_line(l, i, i + 1, &mut tab_index, this.slin_mut());
        // `TestWLineToRLine` inspects the last entry of `slin` (`cxx:1666`).
        test_wline_to_rline(sl, this.slin_mut(), s1, d1, s2, d2, tol_arc);
    }

    // `cxx:1673-1686`: external edge-edge vertices.
    for i in 1..=nbvtx {
        if tab_index[i as usize] != 0 {
            continue;
        }
        let vtx = wline.vertex(i);
        if vtx.is_on_dom_s1() && vtx.is_on_dom_s2() {
            numline += 1;
            tab_index[i as usize] = numline;
            add_line(l, i, i, &mut tab_index, this.slin_mut());
        }
    }
    let _ = tol;
}

/// Geometric arm (`cxx:1690-1867`).
fn perform_geometric(
    this: &mut LineConstructor,
    l: &GeomIntLine,
    s1: &dyn Surface,
    d1: &mut TopolTool,
    s2: &dyn Surface,
    d2: &mut TopolTool,
    mut tol: f64,
) {
    this.slin_mut().clear();
    let Some(gline) = l.as_gline() else {
        return;
    };
    let typl = gline.ityp();
    let nbvtx = gline.nb_vertex();
    let mut tab_index = vec![0i32; (nbvtx + 2) as usize];
    let mut numline = 0;

    let value_at = |t: f64| -> occt_core::gp::GpPnt {
        match &gline.kind {
            GLineKind::Lin(x) => clib::line_value(x, t),
            GLineKind::Circ(c) => clib::circle_value(c, t),
            GLineKind::Elips(e) => clib::ellipse_value(e, t),
            GLineKind::Hypr(h) => clib::hyperbola_value(h, t),
            GLineKind::Parab(p) => clib::parabola_value(p, t),
        }
    };

    let mut intrvtested = false;
    for i in 1..nbvtx {
        let firstp = gline.vertex(i).parameter_on_line();
        let lastp = gline.vertex(i + 1).parameter_on_line();
        if (firstp - lastp).abs() <= PCONFUSION {
            continue;
        }
        intrvtested = true;
        let pmid = (firstp + lastp) * 0.5;
        let p_mid = value_at(pmid);
        let Some((u1, v1, u2, v2)) = uv_on_both(s1, s2, &p_mid) else {
            continue;
        };

        // `cxx:1740-1756`.
        d1.init();
        if d2.more() {
            let (_, a, b) = s1.d1(u1, v1);
            if let Some(c) = d1.value() {
                tol = compute_parametric_tolerance(TopolTool::tol3d_of_curve(c), &a, &b);
            }
        }
        let mut in1 = d1.classify(GpPnt2d::new(u1, v1), tol, false);

        d2.init();
        if in1 != FaceState::Out && d2.more() {
            let (_, a, b) = s2.d1(u2, v2);
            if let Some(c) = d2.value() {
                tol = compute_parametric_tolerance(TopolTool::tol3d_of_curve(c), &a, &b);
            }
        }
        let mut in2 = if in1 != FaceState::Out {
            d2.classify(GpPnt2d::new(u2, v2), tol, false)
        } else {
            FaceState::Out
        };

        // `cxx:1758-1766`.
        if in1 == FaceState::Unknown {
            in1 = FaceState::Out;
        }
        if in2 == FaceState::Unknown {
            in2 = FaceState::Out;
        }

        if both_not_out(in1, in2) {
            numline += 1;
            tab_index[i as usize] = numline;
            tab_index[(i + 1) as usize] = numline;
            add_line(l, i, i + 1, &mut tab_index, this.slin_mut());
        }
    }

    // `cxx:1780-1838`: the closing interval `[last vertex, first vertex + 2π]`
    // of a periodic conic.
    if typl == IntPatchIType::Circle || typl == IntPatchIType::Ellipse {
        let firstp = gline.vertex(nbvtx).parameter_on_line();
        let lastp = std::f64::consts::PI * 2.0 + gline.vertex(1).parameter_on_line();
        let cadrinf = local_first_parameter(l);
        let cadrsup = local_last_parameter(l);
        let mut acadr = (firstp + lastp) * 0.5;
        while acadr < cadrinf {
            acadr += std::f64::consts::PI * 2.0;
        }
        while acadr > cadrsup {
            acadr -= std::f64::consts::PI * 2.0;
        }
        if acadr >= cadrinf && acadr <= cadrsup && (firstp - lastp).abs() > PCONFUSION {
            intrvtested = true;
            let pmid = (firstp + lastp) * 0.5;
            let p_mid = value_at(pmid);
            if let Some((u1, v1, u2, v2)) = uv_on_both(s1, s2, &p_mid) {
                let mut in1 = d1.classify(GpPnt2d::new(u1, v1), tol, false);
                let mut in2 = if in1 != FaceState::Out {
                    d2.classify(GpPnt2d::new(u2, v2), tol, false)
                } else {
                    FaceState::Out
                };
                if in1 == FaceState::Unknown {
                    in1 = FaceState::Out;
                }
                if in2 == FaceState::Unknown {
                    in2 = FaceState::Out;
                }
                if both_not_out(in1, in2) {
                    numline += 1;
                    tab_index[nbvtx as usize] = numline;
                    tab_index[1] = numline;
                    add_line(l, nbvtx, 1, &mut tab_index, this.slin_mut());
                }
            }
        }
    }
    // `cxx:1839-1848` keeps the line a priori when no interval was tested; the
    // `AddLine` body of that branch is commented out in OCCT.
    let _ = intrvtested;

    // `cxx:1852-1865`: external edge-edge vertices.
    for i in 1..=nbvtx {
        if tab_index[i as usize] != 0 {
            continue;
        }
        let vtx = gline.vertex(i);
        if vtx.is_on_dom_s1() && vtx.is_on_dom_s2() {
            numline += 1;
            tab_index[i as usize] = numline;
            add_line(l, i, i, &mut tab_index, this.slin_mut());
        }
    }
}

/// Restriction arm (`cxx:1868-2084`).
fn perform_restriction(
    this: &mut LineConstructor,
    l: &GeomIntLine,
    s1: &dyn Surface,
    d1: &mut TopolTool,
    s2: &dyn Surface,
    d2: &mut TopolTool,
    mut tol: f64,
) {
    this.slin_mut().clear();
    let Some(rline) = l.as_rline() else {
        return;
    };
    let nbvtx = rline.nb_vertex();
    let rest_on_s1 = rline.arc_on_s1;
    let rest_on_s2 = rline.arc_on_s2;
    let mut tab_index = vec![0i32; (nbvtx + 2) as usize];

    // `cxx:1892` / `:2006`: both arms gate the tolerance on `myDom2->More()`,
    // the second one even though it just called `myDom1->Init()`. Only
    // `Init()` and `Value()` are used below, so `More()` is constant here and
    // is evaluated once, exactly as OCCT's repeated calls would answer.
    d2.init();
    let dom2_more = d2.more();

    for i in 1..nbvtx {
        let vtx1 = rline.vertex(i).clone();
        let vtx2 = rline.vertex(i + 1).clone();

        if rest_on_s1 && rest_on_s2 {
            add_line(l, i, i + 1, &mut tab_index, this.slin_mut());
            continue;
        }
        // `cxx:1893-1988` is the `RestOnS1` arm and `:1988-2082` the
        // `RestOnS2` one; only the roles of the two surfaces swap.
        let on_first = rest_on_s1;
        let Some(c2d) = rline.c2d.as_ref() else {
            continue;
        };

        let u0 = vtx1.parameter_on_line();
        let u1_line = vtx2.parameter_on_line();
        if (u1_line - u0).abs() <= PCONFUSION {
            continue;
        }
        let u = (999.0 * u0 + u1_line) * 0.001;
        let p0 = vtx1.p;
        let px2d = c2d.d0(u);
        let host = if on_first { s1 } else { s2 };
        let other = if on_first { s2 } else { s1 };
        let px = host.d0(px2d.x(), px2d.y());
        let p0px = GpVec::from_pnts(&p0, &px);

        let (su1, sv1, su2, sv2) = vtx1.parameters();
        let (ou, ov) = if on_first { (su2, sv2) } else { (su1, sv1) };

        let (_, d1u, d1v) = other.d1(ou, ov);
        let other_dom: &mut TopolTool = if on_first { &mut *d2 } else { &mut *d1 };
        other_dom.init();
        if dom2_more {
            if let Some(c) = other_dom.value() {
                tol = compute_parametric_tolerance(TopolTool::tol3d_of_curve(c), &d1u, &d1v);
            }
        }

        // `cxx:1918-1928`: the classification of both end vertices' other-side
        // parameters gates the whole interval.
        let mut bornin = other_dom.classify(GpPnt2d::new(ou, ov), tol, false);
        if bornin != FaceState::Out {
            let (tu1, tv1, tu2, tv2) = vtx2.parameters();
            let (tu, tv) = if on_first { (tu2, tv2) } else { (tu1, tv1) };
            bornin = other_dom.classify(GpPnt2d::new(tu, tv), tol, false);
        }
        if bornin == FaceState::Out {
            continue;
        }

        // `cxx:1931-1953`: least-squares projection of `P0Px` onto the tangent
        // plane of the other surface.
        let d1u_d1u = d1u.dot(&d1u);
        let d1v_d1v = d1v.dot(&d1v);
        let d1u_d1v = d1u.dot(&d1v);
        let tg_d1u = p0px.dot(&d1u);
        let tg_d1v = p0px.dot(&d1v);
        let dis = d1u_d1u * d1v_d1v - d1u_d1v * d1u_d1v;

        let mut deltau = 1e-10;
        let mut deltav = 1e-10;
        if dis < -1e-10 || dis > 1e-10 {
            deltau = (tg_d1u * d1v_d1v - tg_d1v * d1u_d1v) / dis;
            deltav = (tg_d1v * d1u_d1u - tg_d1u * d1u_d1v) / dis;
        }
        let (nu, nv) = (ou + deltau, ov + deltav);

        // `cxx:1955-1978`: accept when one of the five probes lands inside.
        let mut in2 = other_dom.classify(GpPnt2d::new(nu, nv), tol, false);
        deltau *= 0.05;
        deltav *= 0.05;
        if in2 == FaceState::Out {
            in2 = other_dom.classify(GpPnt2d::new(nu + deltau, nv), tol, false);
        }
        if in2 == FaceState::Out {
            in2 = other_dom.classify(GpPnt2d::new(nu - deltau, nv), tol, false);
        }
        if in2 == FaceState::Out {
            in2 = other_dom.classify(GpPnt2d::new(nu, nv + deltav), tol, false);
        }
        if in2 == FaceState::Out {
            in2 = other_dom.classify(GpPnt2d::new(nu, nv - deltav), tol, false);
        }
        if in2 != FaceState::Out {
            add_line(l, i, i + 1, &mut tab_index, this.slin_mut());
        }
    }
}
