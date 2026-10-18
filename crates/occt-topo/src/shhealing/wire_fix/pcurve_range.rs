use super::*;
use crate::pcurve_full::trimmed_reparam_bspline;

/// `IsKind(STANDARD_TYPE(Geom2d_Line))` (`GeomLib.cxx:864`) as the port can
/// still test it: an exact `Geom2d_Line` is unbounded on both sides with a
/// vanishing second derivative, and the narrower test deliberately excludes the
/// wrappers that forward `is_line` (`Geom2d_TrimmedCurve`, `ReparamCurve2d`),
/// which OCCT would resolve through its own `Geom2d_TrimmedCurve` arm
/// (`cxx:890-901`) rather than the line translate.
pub(in crate::shhealing) fn is_geom2d_line_curve(c: &dyn Curve2d) -> bool {
    !c.first_parameter().is_finite()
        && !c.last_parameter().is_finite()
        && c.d2(0.0).2.square_magnitude() < 1e-30
}

/// `GeomLib::SameRange` (`GeomLib.cxx:842-970`).
///
/// `tolerance` is the first `GeomLib::SameRange` argument (`cxx:842`): the
/// early-out (`cxx:854-859`) and the span-preserving test (`cxx:862`) use
/// it; the interior `Precision::PConfusion()` guards stay as in the cxx. The
/// conversion itself is the single faithful port in
/// [`crate::pcurve_full::trimmed_reparam_bspline`] (the `cxx:912-921` and
/// `cxx:962-968` tails); the unequal-span arm has no kind dispatch in the cxx,
/// so a line, a conic and a B-spline all leave through the same call.
pub(crate) fn geom_lib_same_range(
    tolerance: f64,
    c2d: Arc<dyn Curve2d>,
    first_on: f64,
    last_on: f64,
    req_first: f64,
    req_last: f64,
) -> Arc<dyn Curve2d> {
    // `cxx:854-859`: the ranges already agree; OCCT hands back the same curve.
    if (last_on - req_last).abs() <= tolerance && (first_on - req_first).abs() <= tolerance {
        return c2d;
    }
    // Port-local guard: `BRepLib::SameRange` feeds finite `BRep_GCurve` ranges
    // (`BRepLib.cxx:236-254`) and an unset range would only reach OCCT's
    // `Geom2d_TrimmedCurve` as an infinite bound.
    if !first_on.is_finite() || !last_on.is_finite() {
        return c2d;
    }
    let equal_span =
        ((last_on - first_on) - (req_last - req_first)).abs() <= tolerance;
    if equal_span {
        // `cxx:864-871` `IsKind(Geom2d_Line)`: translate a copy of the line by
        // `dU * Direction` with `dU = FirstOnCurve - RequestedFirst`.
        if is_geom2d_line_curve(c2d.as_ref()) {
            let du = first_on - req_first;
            let (_, dir) = c2d.d1(0.0);
            let mut trsf = GpTrsf2d::identity();
            trsf.set_translation_vec(&GpVec2d::new(dir.x() * du, dir.y() * du));
            return Arc::from(c2d.transformed(&trsf));
        }
        // `cxx:872-889` `IsKind(Geom2d_Circle)`: rotate a copy about its
        // location so `FirstOnCurve` maps onto `RequestedFirst`; the rotation
        // sense follows `Circ2d().IsDirect()`.
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
        // `cxx:890-901`: recurse on `BasisCurve()` with the same arguments, then
        // re-trim the result to `[RequestedFirst, RequestedLast]`.
        if let Some(basis) = c2d.trimmed_basis() {
            let new_basis = geom_lib_same_range(
                tolerance,
                Arc::from(basis.clone_dyn()),
                first_on,
                last_on,
                req_first,
                req_last,
            );
            return Arc::new(Geom2dTrimmedCurve::new(new_basis, req_first, req_last));
        }
        // `cxx:908-922`: the span-preserving generic arm. The guard is OCCT's
        // verbatim `|LastOnCurve - FirstOnCurve| > PConfusion ||
        // |RequestedLast + RequestedFirst| > PConfusion`; the conversion runs on
        // the unclamped `[FirstOnCurve, LastOnCurve]` window.
        if (last_on - first_on).abs() > PCONFUSION
            || (req_last + req_first).abs() > PCONFUSION
        {
            if let Some(bs) =
                trimmed_reparam_bspline(c2d.clone(), first_on, last_on, req_first, req_last)
            {
                return Arc::new(bs);
            }
        }
        // `cxx:908-909` is false (so the arm above is skipped) when both spans
        // and the request sit inside `PConfusion`, and OCCT then assigns nothing,
        // leaving the out-parameter null; `BRepLib::SameRange`
        // stores that null handle as the pcurve (`BRepLib.cxx:236-254`), which the
        // port cannot represent, so the linear stand-in stays instead.
        return reparam_curve2d(c2d, req_first, req_last, first_on, last_on);
    }
    // `cxx:924-969`: unequal span. One path for every kind -- the window
    // resolution (`:934-959`) plus the conversion tail (`:962-968`).
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
    if std::env::var("PCFIX").is_ok() {
        eprintln!(
            "PCFIX edge={:p} face={:p} key={face_key} seam={} n={}",
            Arc::as_ptr(&edge.0.tshape),
            Arc::as_ptr(&face.0.tshape),
            is_seam,
            existing.len()
        );
    }
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
