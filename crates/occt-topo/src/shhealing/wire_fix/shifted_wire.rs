use super::*;

/// `ShapeFix_Wire::FixShifted` (`ShapeFix_Wire.cxx:1661-2126`).
pub fn fix_shifted_wire(wire: &Wire, face: &Face) -> bool {
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    // `ShapeFix_Wire.cxx:1671-1673`: `surf->IsUClosed(Precision())` /
    // `IsVClosed(Precision())` are the `ShapeAnalysis_Surface` versions, not the
    // bare `Geom_Surface` flags, with `Precision() = ShapeFix_Root::myPrecision`
    // (`ShapeFix_Root.lxx:34`, `ShapeFix_Root.cxx:26` = `Precision::Confusion()`).
    // The extra `Geom_SphericalSurface` arm is OCCT's own, because a sphere is
    // closed in V without being V-periodic.
    let u_closed = crate::pcurve_full::sa_is_u_closed(surf.as_ref(), CONFUSION);
    let mut v_closed =
        crate::pcurve_full::sa_is_v_closed(surf.as_ref(), CONFUSION) || surf.gp_sphere().is_some();
    let mut v_range = 1.0;
    let mut v_crv_closed = false;
    if surf.is_surface_of_revolution() {
        // `ShapeFix_Wire.cxx:1682-1707`: for a `Geom_SurfaceOfRevolution` the
        // basis curve decides `vclosed` when it is periodic (CTS18546-2: a 2d
        // contour shifted by 2*PI against a V range of length 2*PI). The
        // `Geom_TrimmedCurve` unwrap is `cxx:1692-1696`; the
        // `Geom_OffsetCurve` unwrap (`cxx:1688-1691`) is UNPORTED, this port
        // has no 3d offset curve type.
        let mut basis = surf.revolution_basis_curve();
        loop {
            let Some(cur) = basis.as_ref() else {
                break;
            };
            if !cur.is_geom_trimmed() {
                break;
            }
            let cur = cur.clone();
            let Some((base, _, _)) = cur.untrimmed_basis() else {
                break;
            };
            basis = Some(base);
        }
        if let Some(basis) = basis {
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
        const MAX_TOL: f64 = SHAPE_FIX_MAX_TOLERANCE;
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
        if du != 0.0 || dv != 0.0 {
            let mut shift = GpTrsf2d::default();
            shift.set_translation_vec(&GpVec2d::new(du, dv));
            replace_pcurve(e2, face, Arc::from(c2.transformed(&shift)));
            done = true;
        }
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

pub(in crate::shhealing) fn try_bi_meridian(
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
