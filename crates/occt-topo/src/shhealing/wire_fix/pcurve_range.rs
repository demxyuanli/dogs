use super::*;

/// `GeomLib::SameRange` (`GeomLib.cxx:842-922`): remap a 2d curve from
/// `[first_on, last_on]` onto `[req_first, req_last]`. Linear reparam matches
/// the BSpline knot `BSplCLib::Reparametrize` evaluation.
pub(in crate::shhealing) fn is_geom2d_line_curve(c: &dyn Curve2d) -> bool {
    !c.first_parameter().is_finite()
        && !c.last_parameter().is_finite()
        && c.d2(0.0).2.square_magnitude() < 1e-30
}

pub(in crate::shhealing) fn same_range_is_line(c: &dyn Curve2d) -> bool {
    is_geom2d_line_curve(c) || c.is_line()
}

/// `Geom2dConvert::CurveToBSplineCurve` Line arm (`Geom2dConvert.cxx:211-225`)
/// then `BSplCLib::Reparametrize` (`GeomLib.cxx:962-968`).
pub(in crate::shhealing) fn line_trim_to_reparam_bspline(
    c2d: Arc<dyn Curve2d>,
    udeb: f64,
    ufin: f64,
    req_first: f64,
    req_last: f64,
) -> Option<Arc<dyn Curve2d>> {
    if (ufin - udeb).abs() <= PCONFUSION {
        return None;
    }
    // Intermediate TrimmedCurve so Start/End and First/Last match Convert.
    let tc = Geom2dTrimmedCurve::new(c2d, udeb, ufin);
    let k0 = tc.first_parameter();
    let k1 = tc.last_parameter();
    if (k1 - k0).abs() <= PCONFUSION {
        return None;
    }
    let p0 = tc.d0(k0);
    let p1 = tc.d0(k1);
    // Flat knots for Mults=[2,2], Degree=1 unique Knots=[k0,k1].
    let mut knots = vec![k0, k0, k1, k1];
    // `BSplCLib::Reparametrize(RequestedFirst, RequestedLast, Knots)`.
    let u_first = req_first.min(req_last);
    let u_last = req_first.max(req_last);
    reparameterize(&mut knots, k0, k1, u_first, u_last);
    let bs = Geom2dBSplineCurve::new(vec![p0.x(), p1.x()], vec![p0.y(), p1.y()], knots, 1).ok()?;
    Some(Arc::new(bs))
}

/// `GeomLib.cxx:908-921` / `:960-968`: build `Geom2dTrimmedCurve(CurvePtr, u1,
/// u2)`, convert with `Geom2dConvert::CurveToBSplineCurve` (default
/// parameterisation), then `BSplCLib::Reparametrize(req_first, req_last)` over
/// the distinct knots followed by `SetKnots`. Covers the non-Line/Circle/
/// Trimmed kinds (BSpline, Bezier, Ellipse, Hyperbola, Parabola, Offset).
pub(in crate::shhealing) fn trim_to_reparam_bspline(
    c2d: Arc<dyn Curve2d>,
    u1: f64,
    u2: f64,
    req_first: f64,
    req_last: f64,
) -> Option<Arc<dyn Curve2d>> {
    if !u1.is_finite() || !u2.is_finite() {
        return None;
    }
    let tc = Geom2dTrimmedCurve::new(c2d, u1, u2);
    let mut bs = occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
        &tc,
        occt_geom2d::geom2d_convert::DEFAULT_PARAMETERISATION,
    )?;
    let (mut distinct, _) = bs.distinct_knots_and_mults();
    occt_core::bspl::knots::reparametrize(req_first, req_last, &mut distinct);
    bs.set_knots(&distinct).ok()?;
    Some(Arc::new(bs))
}

/// `GeomLib.cxx:929-958`: resolve the `[Udeb, Ufin]` window used for the
/// intermediate `Geom2d_TrimmedCurve` of the unequal-span branch.
pub(in crate::shhealing) fn unequal_trim_window(c2d: &dyn Curve2d, first_on: f64, last_on: f64) -> (f64, f64) {
    let mut a_check: &dyn Curve2d = c2d;
    if let Some(basis) = c2d.trimmed_basis() {
        a_check = basis;
    }
    if a_check.is_periodic() {
        // `cxx:934-944`: periodic basis — trim to [FirstOn, LastOn] when span
        // is non-degenerate, else the curve's own parameter window.
        if (last_on - first_on).abs() > PCONFUSION {
            (first_on, last_on)
        } else {
            (c2d.first_parameter(), c2d.last_parameter())
        }
    } else {
        // `cxx:946-958`.
        let cf = c2d.first_parameter();
        let cl = c2d.last_parameter();
        let udeb = if cf.is_finite() { cf.max(first_on) } else { first_on };
        let ufin = if cl.is_finite() { cl.min(last_on) } else { last_on };
        if (ufin - udeb).abs() > PCONFUSION {
            (udeb, ufin)
        } else {
            (cf, cl)
        }
    }
}

/// `GeomLib.cxx:924-969` unequal-span segment for a Line: TrimmedCurve bounds
/// then CurveToBSplineCurve + knot Reparametrize (not a ReparamCurve2d wrapper).
pub(in crate::shhealing) fn same_range_unequal_line(
    c2d: Arc<dyn Curve2d>,
    first_on: f64,
    last_on: f64,
    req_first: f64,
    req_last: f64,
) -> Option<Arc<dyn Curve2d>> {
    let (udeb, ufin) = unequal_trim_window(c2d.as_ref(), first_on, last_on);
    if !udeb.is_finite() || !ufin.is_finite() {
        return None;
    }
    line_trim_to_reparam_bspline(c2d, udeb, ufin, req_first, req_last)
}

/// `GeomLib.cxx:924-969` unequal-span segment for a non-Line kind (same
/// TrimmedCurve + CurveToBSplineCurve + Reparametrize path as `:908-921`).
pub(in crate::shhealing) fn same_range_unequal_other(
    c2d: Arc<dyn Curve2d>,
    first_on: f64,
    last_on: f64,
    req_first: f64,
    req_last: f64,
) -> Option<Arc<dyn Curve2d>> {
    let (udeb, ufin) = unequal_trim_window(c2d.as_ref(), first_on, last_on);
    if !udeb.is_finite() || !ufin.is_finite() {
        return None;
    }
    trim_to_reparam_bspline(c2d, udeb, ufin, req_first, req_last)
}

/// `GeomLib::SameRange` (`GeomLib.cxx:842-969`) for equal-span Line/reparam
/// and unequal-span Line Trimmed+Convert.
pub(crate) fn geom_lib_same_range(
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
    let equal_span =
        ((last_on - first_on) - (req_last - req_first)).abs() <= PCONFUSION;
    if equal_span {
        // `GeomLib.cxx:864-870` `IsKind(Geom2d_Line)` — exact infinite Line only
        // (not Trimmed/Reparam wrappers that forward `is_line`).
        if is_geom2d_line_curve(c2d.as_ref()) {
            let du = first_on - req_first;
            let (_, dir) = c2d.d1(0.0);
            let mut trsf = GpTrsf2d::identity();
            trsf.set_translation_vec(&GpVec2d::new(dir.x() * du, dir.y() * du));
            return Arc::from(c2d.transformed(&trsf));
        }
        // `GeomLib.cxx:871-888` `IsKind(Geom2d_Circle)` — exact Geom2d_Circle
        // only. Rotate a copy about its location so `FirstOnCurve` maps onto
        // `RequestedFirst`; the rotation sense follows `Circ2d().IsDirect()`.
        if c2d.is_geom2d_circle() {
            if let Some(circ) = c2d.gp_circ2d() {
                let p = circ.location();
                let du = if circ.is_direct() {
                    first_on - req_first
                } else {
                    req_first - first_on
                };
                let mut trsf = GpTrsf2d::identity();
                trsf.set_rotation(&p, du);
                return Arc::from(c2d.transformed(&trsf));
            }
        }
        // `GeomLib.cxx:890-900`: TrimmedCurve recurses on the basis then
        // re-trims to RequestedFirst/Last (period-shifted Init2d windows).
        if c2d.trimmed_basis().is_some() {
            let basis = Arc::from(c2d.trimmed_basis().unwrap().clone_dyn());
            let new_basis =
                geom_lib_same_range(basis, first_on, last_on, req_first, req_last);
            return Arc::new(Geom2dTrimmedCurve::new(new_basis, req_first, req_last));
        }
        // `GeomLib.cxx:908-921` equal-span generic arm: TrimmedCurve over
        // [FirstOn, LastOn] -> CurveToBSplineCurve -> knot Reparametrize (the
        // condition is OCCT's verbatim `|LastOn-FirstOn| > PConfusion ||
        // |RequestedLast+RequestedFirst| > PConfusion`). Kinds: BSpline, Bezier,
        // Ellipse, Hyperbola, Parabola, Offset.
        if (last_on - first_on).abs() > PCONFUSION
            || (req_last + req_first).abs() > PCONFUSION
        {
            if let Some(bs) =
                trim_to_reparam_bspline(c2d.clone(), first_on, last_on, req_first, req_last)
            {
                return bs;
            }
        }
        return reparam_curve2d(c2d, req_first, req_last, first_on, last_on);
    }
    // `GeomLib.cxx:924-969`: unequal span — Line uses Trimmed+CurveToBSpline
    // + Reparametrize; every other kind takes the same path via
    // `same_range_unequal_other` (`:960-968`).
    if same_range_is_line(c2d.as_ref()) {
        if let Some(bs) =
            same_range_unequal_line(c2d.clone(), first_on, last_on, req_first, req_last)
        {
            return bs;
        }
    } else if let Some(bs) =
        same_range_unequal_other(c2d.clone(), first_on, last_on, req_first, req_last)
    {
        return bs;
    }
    reparam_curve2d(c2d, req_first, req_last, first_on, last_on)
}

/// `ShapeAnalysis_Edge::CheckPCurveRange` (`ShapeAnalysis_Edge.cxx:999-1032`).
pub(in crate::shhealing) fn check_pcurve_range(first: f64, last: f64, pc: &dyn Curve2d) -> bool {
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
pub(in crate::shhealing) fn translate_pcurve_seam(surf: &dyn Surface, c2d: &dyn Curve2d, tol: f64) -> Arc<dyn Curve2d> {
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

pub(in crate::shhealing) fn translate_pcurve_dir(
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

/// `ShapeFix_Edge::FixAddPCurve` (`ShapeFix_Edge.cxx:470-614`). `prec` is the
/// caller's tolerance; `ShapeFix_Edge.cxx:499` falls back to the edge tolerance
/// when it is not positive, and that value becomes `myPreci` of the projector
/// (`cxx:533`).
pub(in crate::shhealing) fn fix_add_pcurve(
    edge: &Edge,
    face: &Face,
    is_seam: bool,
    prec: f64,
    cache: &mut crate::pcurve_full::ProjectorCache,
) -> bool {
    let face_key = GeometryRegistry::shape_key(&face.0);
    // `ShapeFix_Edge.cxx:475-485`: nothing to do when the edge already carries
    // the pcurve (`ShapeAnalysis_Edge::HasPCurve`) or both seam pcurves
    // (`ShapeAnalysis_Edge::IsSeam`).
    let existing = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
    if (!is_seam && !existing.is_empty()) || (is_seam && existing.len() >= 2) {
        return false;
    }
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    // `ShapeFix_Edge.cxx:487-489`: a pcurve on a plane is not computed. The
    // C++ tests `surf->IsKind(STANDARD_TYPE(Geom_Plane))`, i.e. the *concrete*
    // `Geom_Plane`; a trimmed / offset surface built on a plane is not one and
    // goes to the projector below (`ShapeConstruct_ProjectCurveOnSurface`,
    // whose `projectAnalytic` resolves the wrapper's basis plane).
    if surf.gp_pln().is_some() {
        return false;
    }
    let Some(c3d) = GeometryRegistry::global().edge_curve(&edge.0) else {
        return false;
    };
    let (first, last) = GeometryRegistry::global().edge_parameters(&edge.0);
    // `ShapeFix_Edge.cxx:499`: `preci = (prec > 0. ? prec : BRep_Tool::Tolerance(edge))`.
    let preci = if prec > 0.0 {
        prec
    } else {
        BRepTool::edge_tolerance(edge)
    };
    // `ShapeFix_Edge.cxx:521-531`: `TolFirst` / `TolLast` are the tolerances of
    // the edge's end vertices, `-1` when a vertex is null.
    let tol_first = first_vertex(edge).map_or(-1.0, |v| BRepTool::vertex_tolerance(&v));
    let tol_last = last_vertex(edge).map_or(-1.0, |v| BRepTool::vertex_tolerance(&v));
    let Some(c2d) = project_curve_on_surface_perform(
        c3d.as_ref(),
        surf.as_ref(),
        first,
        last,
        preci,
        tol_first,
        tol_last,
        cache,
    ) else {
        return false;
    };
    if !is_seam {
        GeometryRegistry::global().set_edge_pcurve(&edge.0, face_key, c2d);
        return true;
    }
    let prec = GeometryRegistry::global().edge_tolerance(&edge.0).max(CONFUSION);
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let mut shift = GpTrsf2d::default();
    // `ShapeFix_Edge.cxx:563-579`: which of the two seam pcurves moves by one
    // period is decided by `sas->IsUClosed(prec)` / `IsVClosed(prec)` (the same
    // `preci` as `cxx:499`), falling back to the double-closed
    // `TranslatePCurve` for `sas->IsUClosed()` / `IsVClosed()` with the default
    // `preci = -1`, i.e. `Precision::Confusion()` (`cxx:578`).
    let u_closed = crate::pcurve_full::sa_is_u_closed(surf.as_ref(), preci);
    let v_closed = crate::pcurve_full::sa_is_v_closed(surf.as_ref(), preci);
    let c2d2 = if u_closed && !v_closed {
        shift.set_translation_vec(&GpVec2d::new(ul - uf, 0.0));
        Arc::from(c2d.transformed(&shift))
    } else if v_closed && !u_closed {
        shift.set_translation_vec(&GpVec2d::new(0.0, vl - vf));
        Arc::from(c2d.transformed(&shift))
    } else if crate::pcurve_full::sa_is_u_closed(surf.as_ref(), -1.0)
        && crate::pcurve_full::sa_is_v_closed(surf.as_ref(), -1.0)
    {
        translate_pcurve_seam(surf.as_ref(), c2d.as_ref(), prec)
    } else {
        Arc::from(c2d.clone_dyn())
    };
    GeometryRegistry::global().set_edge_pcurves(&edge.0, face_key, vec![c2d, c2d2]);
    true
}

pub(in crate::shhealing) fn is_seam_use(edges: &[Edge], i: usize) -> bool {
    let key = GeometryRegistry::shape_key(&edges[i].0);
    edges
        .iter()
        .filter(|e| GeometryRegistry::shape_key(&e.0) == key)
        .count()
        >= 2
}

/// `StepToTopoDS_TranslateEdgeLoop::CheckPCurves`
/// (`StepToTopoDS_TranslateEdgeLoop.cxx:106-177`): the planar early return,
/// then the per-edge 2D parameter checks (`cxx:126-172`) and the
/// `XSAlgo_ShapeProcessor::CheckPCurve` advanced check (`cxx:175`).
pub(in crate::shhealing) fn check_pcurve_rep_range(wire: &Wire, face: &Face, preci: f64) {
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
            let preci_p = ((w2 - w1).abs() / 2.0).min(PCONFUSION);
            crate::geom_bnd_lib_elclib2d::adjust_periodic(u1, u2, preci_p, &mut w1, &mut w2);
            reg.set_pcurve_range(&e.0, face_key, w1, w2);
        }
    }
    // `cxx:175` ends the `cxx:120-176` loop body with
    // `XSAlgo_ShapeProcessor::CheckPCurve(myEdge, aFace, preci, sbwd->IsSeam(i))`
    // (`XSAlgo_ShapeProcessor.cxx:344-505`), reached from
    // `StepToTopoDS_TranslateEdgeLoop.cxx:875` `CheckPCurves` and therefore from
    // `step/pcurve_ranges.rs` where this function is invoked. That call is the missing
    // piece behind linkrods FACE 35: the gate at `cxx:392-401`
    // (`aDist11/aDist22 > thePrecision` -> `RemovePCurve` + return) is what
    // leaves the 2.172 edge's tolerance at the imported 1.2432146551310009e-07
    // instead of the 3.415648804e-03 our `fix_same_parameter` writes. Measured
    // on that edge (instrumented `check_pcurve`): d1=1.230581e-04,
    // d2=8.621578e-05. `preci` is `StepToTopoDS_Root::Precision()` set from
    // `STEPControl_ActorRead` (`STEPControl_ActorRead.cxx:2370-2384`): the
    // file's `UNCERTAINTY_MEASURE_WITH_UNIT` length measure times the length
    // factor (linkrods.step declares 2.E-005, so preci = 2e-05), else
    // `read.precision.val` (`Interface_StaticStandards.cxx:39`, default
    // 1.e-03); `step/wire_fix::step_precision` parses that.
    //
    // ENABLED (t336). The two blockers recorded here against enabling it are
    // both gone:
    //  - the `Offset` blow-up (708/876 -> 30264/57316) was not caused by this
    //    gate: it was `pcurve_on_reversed` using `Oriented(TopAbs_REVERSED)`
    //    (`TopoDS_Shape::Oriented`, which *sets* the orientation) where OCCT
    //    uses `TopoDS::Edge(theEdge.Reversed())` (`TopoDS_Shape::Reversed()`,
    //    which *flips* FORWARD <-> REVERSED). After that fix the gate on/off
    //    measures `Offset` 712/892 both ways.
    //  - the 53 dropped linkrods pcurves were a symptom of the same wrong
    //    pcurve selection on reversed edges, not of `preci`.
    // Measured on/off over all 16 models (`export_data_obj`, preci from the
    // STEP file): only `linkrods` 3461/5007 -> 3485/5060 (OCCT reference
    // 3494/5078) and `Shape-1` 3327/4304 -> 3291/4230 (OCCT 3343/4336) move;
    // every other model, including all seven hard holds, is unchanged. Cost:
    // the full run goes 44.8 s -> 279 s.
    //
    // The re-projection tail `cxx:403-505` is separately gated by
    // `CHECK_PCURVE_REPROJECT` in `shhealing/pcurve_ranges.rs`; measured inert for all 16
    // models, and still parked for the reason recorded there.
    if T312_WIRE_XSALGO_CHECKPCURVE {
        check_pcurves_xsalgo(wire, face, preci);
    }
}
