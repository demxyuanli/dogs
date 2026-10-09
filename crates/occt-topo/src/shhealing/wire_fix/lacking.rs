use super::*;

/// `ShapeAnalysis_Wire::CheckLacking(num, Tolerance, p2d1, p2d2)`
/// (`ShapeAnalysis_Wire.cxx:1711-1793`).
///
/// `None` covers both `FAIL` (missing vertex / pcurve, `cxx:1732-1750`) and the
/// "no gap" exit `myMax2d < tol2d * tol2d` (`cxx:1776-1780`) - the two states
/// `FixLacking` treats the same way (`cxx:3629-3635`).
pub(in crate::shhealing) struct LackingCheck {
    p2d1: GpPnt2d,
    p2d2: GpPnt2d,
    /// `MaxDistance2d()` (`cxx:1782`).
    max2d: f64,
    /// `MaxDistance3d()` (`cxx:1783`).
    max3d: f64,
    /// `LastCheckStatus(ShapeExtend_DONE2)` (`cxx:1786-1791`).
    done2: bool,
}

pub(in crate::shhealing) fn check_lacking(wire: &Wire, face: &Face, num: usize, tolerance: f64) -> Option<LackingCheck> {
    let nb = wire_edges_nb(wire);
    if nb < 1 {
        return None;
    }
    // `cxx:1723-1726`.
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let surf = BRepTool::face_surface(face)?;
    // `cxx:1729-1741`: `LastVertex(E1)` / `FirstVertex(E2)` through
    // `BRepTools::Compare`.
    let (v1, v2) = (last_vertex(&e1)?, first_vertex(&e2)?);
    if !vertices_coincide(&v1, &v2) {
        return None;
    }
    // `cxx:1746-1767`: `sae.PCurve(..., orient=true)` orders the range as
    // FirstVertex -> LastVertex, so `b` is the junction parameter of E1 and `a`
    // the junction parameter of E2 (`ShapeAnalysis_Edge.cxx:218-227`).
    let (c2d1, _a1, b1) = curve_on_surface_oriented(&e1, face, true)?;
    let (p2d1, mut v1t) = c2d1.d1(b1);
    if e1.0.orientation().is_reversed() {
        v1t.reverse();
    }
    let (c2d2, a2, _b2) = curve_on_surface_oriented(&e2, face, true)?;
    let (p2d2, mut v2t) = c2d2.d1(a2);
    if e2.0.orientation().is_reversed() {
        v2t.reverse();
    }
    // `cxx:1768-1780`.
    let v12 = GpVec2d::new(p2d2.x() - p2d1.x(), p2d2.y() - p2d1.y());
    let mut max2d = v12.square_magnitude();
    let mut tol = BRepTool::vertex_tolerance(&v1).max(BRepTool::vertex_tolerance(&v2));
    if tolerance > occt_core::precision::REAL_SMALL && tolerance < tol {
        tol = tolerance; // `cxx:1773`
    }
    let tol2d = 2.0
        * occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol)
            .max(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol));
    if max2d < tol2d * tol2d {
        return None;
    }
    max2d = max2d.sqrt();
    let max3d = tol * max2d / tol2d.max(occt_core::precision::REAL_SMALL);
    // `cxx:1786-1791`: `myMax2d < PConfusion`, or the 2d gap points against the
    // wire direction at either end ("back-going" zigzag).
    let done2 = max2d < PCONFUSION
        || (v1t.square_magnitude() > occt_core::precision::REAL_SMALL
            && v12.angle(&v1t).abs() > 0.9 * std::f64::consts::PI)
        || (v2t.square_magnitude() > occt_core::precision::REAL_SMALL
            && v12.angle(&v2t).abs() > 0.9 * std::f64::consts::PI);
    Some(LackingCheck {
        p2d1,
        p2d2,
        max2d,
        max3d,
        done2,
    })
}

/// `TryNewPCurve` (`ShapeFix_Wire.cxx:2215-2246`): deviation of `c2d` against
/// the edge's own 3D curve, measured by `ShapeFix_Edge::FixSameParameter`
/// (`cxx:2241`) on a temporary edge.
pub(in crate::shhealing) fn try_new_pcurve(
    edge: &Edge,
    face: &Face,
    c2d: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
) -> Option<(Arc<dyn Curve2d>, f64, f64, f64)> {
    let curve = BRepTool::edge_curve(edge)?; // `cxx:2223-2227`
    let (f, l) = BRepTool::edge_parameters(edge);
    // `cxx:2230-2236`: `BRepBuilderAPI_MakeEdge(crv, f, l)` + `SetRange3d`.
    let builder = TopoBuilder::new();
    let tmp = builder.make_edge(curve, f, l);
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    // `cxx:2238-2243`: `BRepBuilderAPI_MakeEdge` leaves the edge tolerance at
    // `BRep_TEdge`'s constructor value `RealEpsilon()`
    // (`BRep_TEdge.cxx:32-38`), and `B.UpdateEdge(edge, c2d, face, 0.)`
    // (`cxx:2239`) is `TE->UpdateTolerance(0.)`, a max, so it stays there.
    // `B.Range(edge, face, first, last)` (`cxx:2240`), `B.SameRange(edge,
    // false)` (`cxx:2242`) and no `SameParameter` flag (`cxx:2243`).
    reg.set_edge_pcurve(&tmp.0, face_key, c2d);
    reg.set_pcurve_range(&tmp.0, face_key, first, last);
    reg.set_edge_tolerance(&tmp.0, occt_core::precision::COMPUTATIONAL);
    reg.set_same_range(&tmp.0, false);
    // `cxx:2244`: `sfe->FixSameParameter(edge, face)`.
    fix_same_parameter(&tmp, face);
    // `cxx:2245-2246`.
    let (new_c2d, nf, nl) = curve_on_surface_oriented(&tmp, face, false)?;
    Some((new_c2d, nf, nl, BRepTool::edge_tolerance(&tmp)))
}

/// The four in-place outputs of `TryBendingPCurve`
/// (`ShapeFix_Wire.cxx:3534-3538`): `c2d`, `first`, `last`, `tol`. The C++
/// function assigns only the ones its control flow reaches, so a failure after
/// `c2d = bs` (`cxx:3593`) leaves `c2d` assigned at the new candidate while
/// `first`/`last` keep the range `sae.PCurve` wrote (`cxx:3542`) and `tol` keeps
/// its previous value. `None` models a null handle.
pub(in crate::shhealing) struct BendParam {
    curve: Option<Arc<dyn Curve2d>>,
    first: f64,
    last: f64,
    tol: f64,
}

impl BendParam {
    /// `double bendtol1 = 0., bendtol2 = 0.; ... bendf1 = 0., bendl1 = 0. ...`
    /// (`cxx:3664-3666`).
    fn new() -> Self {
        BendParam {
            curve: None,
            first: 0.0,
            last: 0.0,
            tol: 0.0,
        }
    }
}

/// `TryBendingPCurve(E, face, p2d, end, c2d, first, last, tol)`
/// (`ShapeFix_Wire.cxx:3532-3613`). `out` is mutated in place, exactly as the
/// C++ reference parameters are.
pub(in crate::shhealing) fn try_bending_pcurve(
    edge: &Edge,
    face: &Face,
    p2d: GpPnt2d,
    end: bool,
    out: &mut BendParam,
) -> bool {
    // `cxx:3542-3545`: `sae.PCurve(E, face, c2d, first, last, false)`.
    let Some((c2d, first, last)) = curve_on_surface_oriented(edge, face, false) else {
        return false;
    };
    out.curve = Some(c2d.clone());
    out.first = first;
    out.last = last;
    // `cxx:3550-3586`: `c2d->IsKind(Geom2d_BSplineCurve)` (with a seam, the
    // second pcurve) takes the copy arm; every other kind goes through
    // `Geom2dTrimmedCurve` + `Geom2dConvert::CurveToBSplineCurve`
    // (`cxx:3556-3559`). This port does not expose a knot vector, so the
    // `SetPole` edit cannot be reproduced on an imported B-spline pcurve; only
    // the Line arm (`Geom2dConvert.cxx:211-225`, mirrored by
    // `line_trim_to_reparam_bspline` above) is ported. Any other kind reports
    // "no bend candidate" - the same outcome OCCT gives when its own
    // `CurveToBSplineCurve` throws (`cxx:3605-3612`).
    if !is_geom2d_line_curve(c2d.as_ref()) {
        return false;
    }
    // `cxx:3555-3559`: the trimmed-line B-spline. `line_trim_to_reparam_bspline`
    // builds the same curve for `GeomLib::SameRange`; here the knots stay the
    // trimmed curve's own, because `Segment` is never reached (below).
    let tc = Geom2dTrimmedCurve::new(c2d.clone(), first, last);
    let (k0, k1) = (tc.first_parameter(), tc.last_parameter());
    let (u0, u1) = (k0, k1);
    if (u1 - u0).abs() <= PCONFUSION {
        return false;
    }
    let (p_at_u0, p_at_u1) = if first <= last {
        (tc.d0(first), tc.d0(last))
    } else {
        (tc.d0(last), tc.d0(first))
    };
    // `cxx:3567-3591`: `par = (end ? last : first)`; the pole test succeeds
    // because the flat knot vector has multiplicity 2 > degree 1 at both ends
    // (`Multiplicity(1) > Degree()` at `cxx:3569`, `Multiplicity(NbKnots) >
    // Degree()` at `cxx:3573`), so `SetPole` always applies at `par` and
    // `Segment` at `cxx:3579` is unreachable for this arm.
    let par = if end { last } else { first };
    let (np0, np1) = if (par - u0).abs() <= (par - u1).abs() {
        (p2d, p_at_u1) // `SetPole(1, p2d)`
    } else {
        (p_at_u0, p2d) // `SetPole(NbPoles, p2d)`
    };
    let Ok(bs) = Geom2dBSplineCurve::new(
        vec![np0.x(), np1.x()],
        vec![np0.y(), np1.y()],
        vec![u0, u0, u1, u1],
        1,
    ) else {
        return false;
    };
    let candidate: Arc<dyn Curve2d> = Arc::new(bs);
    // `cxx:3593`: `c2d = bs` happens before the tolerance check.
    out.curve = Some(candidate.clone());
    // `cxx:3596-3601`.
    match try_new_pcurve(edge, face, candidate, first, last) {
        Some((c, f, l, tol)) => {
            out.curve = Some(c);
            out.first = f;
            out.last = l;
            out.tol = tol;
            true
        }
        None => false,
    }
}

/// `ShapeFix_Wire::FixLacking(const int num, const bool force)`
/// (`ShapeFix_Wire.cxx:3617-3975`).
///
/// NOT PORTED: the `Context()->Replace` calls at `cxx:3881-3899` only update
/// the `ShapeBuild_ReShape` map (this pass runs without a context), and the two
/// `BRep_Tool::Degenerated` arms at `cxx:3819-3825` have empty bodies in OCCT
/// itself.
pub(in crate::shhealing) fn fix_lacking_one(
    wire: &mut Wire,
    face: &Face,
    num: usize,
    force: bool,
    preci: f64,
    max_tol: f64,
) -> bool {
    // `ShapeAnalysis_Wire::IsReady()`: `!myWire.IsNull() && NbEdges() > 0`
    // (`ShapeAnalysis_Wire.hxx`), checked at `cxx:3620-3623`.
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    // `cxx:3627-3635`: `CheckLacking(num, force ? Precision() : 0.)`, then bail
    // unless it reported `DONE`.
    let tolerance = if force { preci } else { 0.0 };
    let Some(check) = check_lacking(wire, face, num, tolerance) else {
        return false;
    };
    // `cxx:3640-3647`.
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let (Some(v1), Some(v2)) = (last_vertex(&e1), first_vertex(&e2)) else {
        return false;
    };
    // `cxx:3652`.
    let tol = BRepTool::vertex_tolerance(&v1).max(BRepTool::vertex_tolerance(&v2));
    let dist2d = check.max2d;
    let inctol = check.max3d;
    let p2d1 = check.p2d1;
    let p2d2 = check.p2d2;
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    // `BRep_Tool::IsClosed(E, face)` (`BRep_Tool.cxx:795-841`): a plane is
    // never closed; otherwise the edge carries two pcurves on this face.
    let is_closed = |e: &Edge| crate::brep_tool::BRepTool::is_closed_edge_face(e, face);
    let tol_e1 = BRepTool::edge_tolerance(&e1);
    let tol_e2 = BRepTool::edge_tolerance(&e2);
    let mut tol1 = CONFUSION; // `cxx:3659-3660`
    let mut tol2 = CONFUSION;
    let mut p3d1: Option<GpPnt> = None;
    let mut p3d2: Option<GpPnt> = None;

    // `cxx:3662-3713`: bend speculation. `myGeomMode` is true
    // (`ShapeFix_Wire.cxx:169`).
    let mut bp1 = BendParam::new();
    let mut bp2 = BendParam::new();
    if !is_closed(&e1) && !is_closed(&e2) {
        let mid = GpPnt2d::new(0.5 * (p2d1.x() + p2d2.x()), 0.5 * (p2d1.y() + p2d2.y()));
        let end1 = e1.0.orientation() == Orientation::Forward;
        let end2 = e2.0.orientation() == Orientation::Reversed;
        let mut ok1 = try_bending_pcurve(&e1, face, mid, end1, &mut bp1);
        let mut ok2 = try_bending_pcurve(&e2, face, mid, end2, &mut bp2);
        if ok1 && !ok2 {
            // `cxx:3686-3697`.
            bp2.tol = tol_e2;
            ok1 = try_bending_pcurve(&e1, face, p2d2, end1, &mut bp1);
        } else if !ok1 && ok2 {
            // `cxx:3698-3709`. OCCT passes `E2.Orientation() == TopAbs_FORWARD`
            // for this retry while the first call at `cxx:3681` tested
            // REVERSED; reproduced as written.
            bp1.tol = tol_e1;
            let end2_retry = e2.0.orientation() == Orientation::Forward;
            ok2 = try_bending_pcurve(&e2, face, p2d1, end2_retry, &mut bp2);
        }
        if !ok1 && !ok2 {
            bp1.curve = None; // `cxx:3710-3712`: `bendc1.Nullify()`
        }
    }

    // `cxx:3716-3817`: selector of solutions.
    let mut do_increase = false;
    let mut do_add_long = false;
    let mut do_add_closed = false;
    let mut do_add_degen = false;
    let mut do_bend = false;
    let bendtol1 = bp1.tol;
    let bendtol2 = bp2.tol;
    let bend_ready = bp1.curve.is_some() && bp2.curve.is_some();
    if bend_ready
        && ((bendtol1 < tol_e1 && bendtol2 < tol_e2)
            || (inctol < preci && bendtol1 < inctol && bendtol2 < inctol))
    {
        do_bend = true; // `cxx:3727-3732`
    } else if inctol < preci {
        do_increase = true; // `cxx:3735-3738`
    } else if !reg.is_degenerated_edge(&e2.0) && !reg.is_degenerated_edge(&e1.0) {
        // `cxx:3744-3745`: `myTopoMode` is forced true for the face path by
        // `ShapeFix_Shape.cxx:200`, so this block always runs.
        {
            // `cxx:3748-3762`: `sae.Curve3d(E1, c3d, a, b, true)` gives `b` at
            // LastVertex and `a` at FirstVertex.
            let Some(c1) = BRepTool::edge_curve(&e1) else {
                return false; // `cxx:3749-3753` FAIL1
            };
            let (mut a1, mut b1p) = BRepTool::edge_parameters(&e1);
            if e1.0.orientation().is_reversed() {
                std::mem::swap(&mut a1, &mut b1p);
            }
            let q1 = c1.d0(b1p);
            let dist2d3d1 = q1.distance(&surf.d0(p2d1.x(), p2d1.y()));
            let Some(c2) = BRepTool::edge_curve(&e2) else {
                return false; // `cxx:3756-3760` FAIL1
            };
            let (mut a2, mut b2p) = BRepTool::edge_parameters(&e2);
            if e2.0.orientation().is_reversed() {
                std::mem::swap(&mut a2, &mut b2p);
            }
            let q2 = c2.d0(a2);
            let dist2d3d2 = q2.distance(&surf.d0(p2d2.x(), p2d2.y()));
            // `cxx:3764-3769`.
            tol1 = tol_e1.max(dist2d3d1);
            tol2 = tol_e2.max(dist2d3d2);
            let tol0 = tol1 + tol2;
            let dist3d2 = q1.square_distance(&q2);
            p3d1 = Some(q1);
            p3d2 = Some(q2);
            if !check.done2
                && dist3d2 > 1.25 * tol0 * tol0
                && (force || dist3d2 > preci * preci || inctol > max_tol)
            {
                do_add_long = true; // `cxx:3771-3777`
            }
        }
        // `cxx:3781-3791`.
        if !do_add_long
            && inctol < max_tol
            && !is_degenerated_2d(surf.as_ref(), p2d1, p2d2, 2.0 * tol, 10.0)
        {
            if bend_ready && bendtol1 < inctol && bendtol2 < inctol {
                do_bend = true;
            } else {
                do_increase = true;
            }
        } else if !do_add_long {
            // `cxx:3799-3815`.
            let p1 = BRepTool::vertex_point(&v1);
            let p2 = BRepTool::vertex_point(&v2);
            let pv = GpPnt::new(
                0.5 * (p1.x() + p2.x()),
                0.5 * (p1.y() + p2.y()),
                0.5 * (p1.z() + p2.z()),
            );
            let pm = surf.d0(0.5 * (p2d1.x() + p2d2.x()), 0.5 * (p2d1.y() + p2d2.y()));
            let dist = pv.distance(&pm);
            if dist <= tol {
                do_add_degen = true;
            } else {
                // `cxx:3806`: `myTopoMode`, which is true here. The
                // `dist <= MaxTolerance()` fallback at `cxx:3810-3814` belongs
                // to the `!myTopoMode` arm and is unreachable.
                do_add_closed = true;
            }
        }
    }

    // `cxx:3830-3921`: add the new edge. `DONE2` is set at `cxx:3920` on all
    // three arms (`DONE3` at `cxx:3910` for the degenerated edge, `DONE4` at
    // `cxx:3916` for the closed one), so every arm reports a fix.
    let mut done = false;
    if do_add_long || do_add_degen || do_add_closed {
        let builder = TopoBuilder::new();
        // `cxx:3835-3849`.
        let (new_v1, new_v2) = if do_add_long {
            let v1n = builder.make_vertex(p3d1.expect("doAddLong sets p3d1"), 0.0);
            let v2n = builder.make_vertex(p3d2.expect("doAddLong sets p3d2"), 0.0);
            v1n.set_tolerance(1.001 * tol1); // `cxx:3846`
            v2n.set_tolerance(1.001 * tol2); // `cxx:3847`
            (v1n, v2n)
        } else {
            (v1.clone(), v2.clone())
        };
        // `cxx:3852-3862`: `B.MakeEdge`, `B.Degenerated(edge, true)` before the
        // curve, `B.UpdateEdge(edge, theLine2d, face, Precision::Confusion())`,
        // `B.Range(edge, face, 0, dist2d)`.
        let v12 = GpVec2d::new(p2d2.x() - p2d1.x(), p2d2.y() - p2d1.y());
        let Ok(dir2d) = GpDir2d::from_vec2d(&v12) else {
            return false;
        };
        let line2d: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2d1, dir2d));
        // `ShapeBuild_Edge::BuildCurve3d` (`cxx:3866`) approximates the 3D image
        // of the pcurve on the face. This port's `EdgeGeom` always needs a 3D
        // curve, and the pcurve image is the curve `BRepAdaptor_Curve` resolves
        // once the same edge exists - the adapter `fix_degenerated` uses for its
        // own curve-less edge (`ShapeFix_Wire.cxx:2165-2168`).
        let curve3d: Arc<dyn Curve> = Arc::new(crate::meshing::edge_discret::CurveOnSurface::new(
            line2d.clone(),
            surf.clone(),
            0.0,
            dist2d,
        ));
        let mut new_edge = builder.make_edge(curve3d, 0.0, dist2d);
        reg.set_edge_pcurve(&new_edge.0, face_key, line2d);
        reg.set_pcurve_range(&new_edge.0, face_key, 0.0, dist2d);
        reg.set_edge_tolerance(&new_edge.0, CONFUSION);
        if do_add_degen {
            reg.set_degenerated(&new_edge.0, true);
        }
        builder.add_edge_vertices(&mut new_edge, &new_v1, &new_v2);
        // `cxx:3873-3901`: `doAddLong` re-points the two adjacent edges at the
        // new vertices.
        if do_add_long {
            let first_arg = if n1 == n2 { Some(&new_v2) } else { None };
            let ne1 = copy_replace_vertices_with(&e1, first_arg, Some(&new_v1));
            wire_set_edge_composed(wire, n1, &ne1);
            if n1 != n2 {
                let ne2 = copy_replace_vertices_with(&e2, Some(&new_v2), None);
                wire_set_edge_composed(wire, n2, &ne2);
            }
        }
        // `cxx:3905-3920`: `DONE2` is encoded on every arm of this branch.
        done = true;
        // `ShapeExtend_WireData::Add(edge, n2)` (`cxx:3919`) inserts before
        // `n2`; the stored children keep the wire's own orientation composed
        // (`ShapeExtend_WireData.cxx:114-121`), so undo it here as
        // `wire_set_edge_composed` does for `Set`.
        let mut ins = new_edge.clone();
        if wire.0.orientation() == Orientation::Reversed {
            ins.0.reverse();
        }
        wire_insert_edge_before(wire, n2, &ins);
    } else if inctol > tol && inctol < max_tol {
        // `cxx:3924-3933`.
        if bend_ready && bendtol1 < inctol && bendtol2 < inctol {
            do_bend = true;
        } else {
            do_increase = true;
        }
    }

    // `cxx:3937-3956`: bend the pcurves. `B.UpdateEdge(E, bendc, face, bendtol)`
    // writes the pcurve, the COS range and the edge tolerance (`TE->
    // UpdateTolerance(Tol)` is a max, `BRep_Builder.cxx:655-671` /
    // `BRep_TEdge.hxx`).
    if do_bend {
        if let Some(bc1) = bp1.curve.as_ref() {
            reg.set_edge_pcurve(&e1.0, face_key, bc1.clone());
            reg.set_pcurve_range(&e1.0, face_key, bp1.first, bp1.last);
            reg.set_edge_tolerance(&e1.0, tol_e1.max(bendtol1));
        }
        if let Some(bc2) = bp2.curve.as_ref() {
            reg.set_edge_pcurve(&e2.0, face_key, bc2.clone());
            reg.set_pcurve_range(&e2.0, face_key, bp2.first, bp2.last);
            reg.set_edge_tolerance(&e2.0, tol_e2.max(bendtol2));
        }
        // `cxx:3943-3946`: `B.UpdateVertex` on all four vertices.
        for (v, t) in [
            (first_vertex(&e1), bendtol1),
            (last_vertex(&e1), bendtol1),
            (first_vertex(&e2), bendtol2),
            (last_vertex(&e2), bendtol2),
        ] {
            if let Some(v) = v {
                v.set_tolerance(t);
            }
        }
        // `cxx:3948-3951`: re-check the bent edges for self-intersections.
        let _ = fix_self_intersecting_edge(wire, face, n1, preci);
        let _ = fix_self_intersecting_edge(wire, face, n2, preci);
        let _ = fix_intersecting_edges(wire, face, n2, preci);
        done = true;
    }

    // `cxx:3959-3963`.
    if do_increase {
        v1.set_tolerance(1.001 * inctol);
        v2.set_tolerance(1.001 * inctol);
        done = true;
    }

    done
}

/// `ShapeFix_Wire::FixLacking()` (`ShapeFix_Wire.cxx:1284-1295`) - `Perform`
/// calls it at `cxx:448` under `NeedFix(myFixLackingMode, ReorderOK)` with
/// `FromSTEP.FixShape.FixLackingMode = -1` (`STEPControl_Controller.cxx:237`),
/// so the gate is `ReorderOK`. Wired in `check_pcurves_and_shift` below.
///
/// The pass is translated branch for branch, including the `doBend` tail's
/// `FixSelfIntersectingEdge` / `FixIntersectingEdges` calls
/// (`ShapeFix_Wire.cxx:3948-3951`).
/// Previously it could not be wired: the in-memory cylinder used by
/// `step::tests::cylinder_roundtrip` had a lateral wire whose two rings ran the
/// same way (sum(du) = +4pi), so the closure junction showed
/// `p2d1 = (2*pi, 0)` / `p2d2 = (-2*pi, 0)` - both mapping to the same 3D
/// vertex `(2, 0, 0)`. `ShapeAnalysis_Wire::CheckLacking`
/// (`ShapeAnalysis_Wire.cxx:1767-1790`) has no period normalization, so it
/// reported DONE1 with `myMax2d = 4*pi` / `myMax3d = 2*pi` and `FixLacking`
/// selected the degenerated-edge branch (`cxx:3782-3806`). That malformed input
/// wire is fixed in `primitives.rs` (`BRepPrimCylinder::make_cylinder` now
/// traverses the two rings oppositely, as OCCT's swept lateral face does), so
/// every junction of the cylinder's lateral loop is gap-free and this pass is a
/// no-op there.

pub(in crate::shhealing) fn fix_lacking_all(wire: &mut Wire, face: &Face, force: bool, preci: f64, max_tol: f64) -> bool {
    let mut done = false;
    // `cxx:1291-1295`: `myClosedMode` is true (`ShapeFix_Wire.cxx:171`), so
    // `start = 1`; `NbEdges()` is re-read every iteration, so an inserted edge
    // extends the loop.
    let mut i = 1usize;
    while i <= wire_edges_nb(wire) {
        done |= fix_lacking_one(wire, face, i, force, preci, max_tol);
        i += 1;
    }
    done
}

