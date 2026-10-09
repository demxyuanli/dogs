use super::*;

// ---------------------------------------------------------------------------
// ShapeFix_Wire::FixNotchedEdges / ShapeAnalysis_Wire::CheckNotchedEdges
// ---------------------------------------------------------------------------

/// `Geom_Curve::IsClosed()`: `Geom_Circle::IsClosed()` (`Geom_Circle.cxx:73-76`)
/// and `Geom_Ellipse` (base `Geom_Curve::IsClosed()` is `IsPeriodic()`) are the
/// periodic flag; `Geom_TrimmedCurve::IsClosed()` (`Geom_TrimmedCurve.cxx:155-168`)
/// is `IsPeriodic()` (whole periods of a periodic basis, see
/// `occt-geom/trimmed.rs:159-172`) or coincident ends; `Geom_BSplineCurve`
/// (`Geom_BSplineCurve_1.cxx:146-149`) and `Geom_OffsetCurve`
/// (`Geom_OffsetCurve.cxx:436-442`) compare the ends against
/// `Precision::Computational()`, standing in for `Geom_BezierCurve`'s
/// construction flag too.
///
/// NOTE this is NOT `GeomAdaptor_Curve::IsClosed()`
/// (`GeomAdaptor_Curve.cxx:576-584`), which only compares the *window* ends and
/// is what `shape_analysis_curve::project_act` needs (`cxx:340`).
pub(in crate::shhealing) fn geom_curve_is_closed(c: &dyn Curve) -> bool {
    if c.is_periodic() {
        return true;
    }
    let (f, l) = (c.first_parameter(), c.last_parameter());
    f.is_finite()
        && l.is_finite()
        && c.d0(f).square_distance(&c.d0(l)) <= occt_core::precision::COMPUTATIONAL
}

/// `ShapeAnalysis_Edge::IsClosed3d(edge)` (`ShapeAnalysis_Edge.cxx:129-142`):
/// the edge carries a 3D curve that is closed, and its first and last vertices
/// are the same `TShape` (`FirstVertex` / `LastVertex` with the default
/// `CumOri = false`, i.e. the stored FORWARD / REVERSED child).
pub(in crate::shhealing) fn edge_is_closed_3d(edge: &Edge) -> bool {
    let Some(c3d) = BRepTool::edge_curve(edge) else {
        return false;
    };
    if !geom_curve_is_closed(c3d.as_ref()) {
        return false;
    }
    match (edge_vertices(edge).0, edge_vertices(edge).1) {
        (Some(v1), Some(v2)) => is_same(&v1.0, &v2.0),
        _ => false,
    }
}

/// `ShapeAnalysis_Wire::ProjectInside` (`ShapeAnalysis_Wire.cxx:1837-1860`).
pub(in crate::shhealing) fn project_inside(
    ad: &crate::meshing::edge_discret::CurveOnSurface,
    pnt: &GpPnt,
    preci: f64,
    adjust_to_ends: bool,
) -> Projection {
    let proj = project_adaptor(ad, pnt, preci, adjust_to_ends);
    let u_first = ad.first_parameter();
    let u_last = ad.last_parameter();
    if proj.param < u_first {
        let p = ad.d0(u_first);
        return Projection { distance: p.distance(pnt), point: p, param: u_first };
    }
    if proj.param > u_last {
        let p = ad.d0(u_last);
        return Projection { distance: p.distance(pnt), point: p, param: u_last };
    }
    proj
}

/// `ShapeAnalysis_Wire::CheckNotchedEdges(num, shortNum, param, Tolerance)`
/// (`ShapeAnalysis_Wire.cxx:1864-2000`). `Some` carries the `shortNum` and
/// `param` out-parameters; the `ShapeExtend_FAIL*` bits it leaves in `myStatus`
/// (`cxx:1890`, `:1896`, `:1906`, `:1919`) are not modelled, as no caller of
/// this port reads them.
pub(in crate::shhealing) struct NotchedCheck {
    short_num: usize,
    param: f64,
}

pub(in crate::shhealing) fn check_notched_edges(
    wire: &Wire,
    face: &Face,
    num: usize,
    tolerance: f64,
) -> Option<NotchedCheck> {
    // `cxx:1869-1878`: `IsReady()` and the `num == 0` meaning "the last edge".
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    if nb == 0 {
        return None;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let e1 = edges.get(n1.checked_sub(1)?)?.clone();
    let e2 = edges.get(n2.checked_sub(1)?)?.clone();
    if BRepTool::is_degenerated(&e1) || BRepTool::is_degenerated(&e2) {
        return None; // `cxx:1880-1883`
    }

    // `cxx:1885-1898`: `V1 = sae.LastVertex(E1)`, `V2 = sae.FirstVertex(E2)`,
    // then `BRepTools::Compare(V1, V2)`. Both accessors are orientation-aware
    // (`ShapeAnalysis_Edge.cxx:228-258`): on a REVERSED edge the first vertex
    // is the raw *last* one and vice versa. Reading the raw `edge_vertices`
    // here swapped the pair for every REVERSED edge, so the notched pair was
    // never recognised (F113 of `data/occ/a3n00.stp` is exactly that).
    let (Some(v1), Some(v2)) = (last_vertex(&e1), first_vertex(&e2)) else {
        return None;
    };
    if !vertices_coincide(&v1, &v2) {
        return None;
    }

    // `cxx:1900-1932`: the end tangents of both pcurves.
    let (c2d1, a1, b1) = curve_on_surface_oriented(&e1, face, false)?;
    let (c2d2, a2, b2) = curve_on_surface_oriented(&e2, face, false)?;
    let (p2d1, tan1) = if e1.0.orientation().is_reversed() {
        c2d1.d1(a1)
    } else {
        let (p, d) = c2d1.d1(b1);
        (p, d.reversed())
    };
    let (p2d2, tan2) = if e2.0.orientation().is_reversed() {
        let (p, d) = c2d2.d1(b2);
        (p, d.reversed())
    } else {
        c2d2.d1(a2)
    };
    // `cxx:1934-1937`: `gp::Resolution()` is `RealSmall()`.
    if tan2.magnitude() < REAL_SMALL || tan1.magnitude() < REAL_SMALL {
        return None;
    }
    // `cxx:1939-1942`.
    if tan2.angle(&tan1).abs() > 0.1 || p2d1.distance(&p2d2) > tolerance {
        return None;
    }

    // `cxx:1944-1952`: both pcurves are read as 3D curves lying in the XY plane
    // (`Geom_Plane(gp_Pln())`), which is what `Adaptor3d_CurveOnSurface` over a
    // plane resolves.
    let plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
    let ad1 = crate::meshing::edge_discret::CurveOnSurface::new(c2d1.clone(), plane.clone(), a1, b1);
    let ad2 = crate::meshing::edge_discret::CurveOnSurface::new(c2d2.clone(), plane.clone(), a2, b2);

    // `cxx:1954-1968`.
    let pt2 = if e2.0.orientation().is_forward() { c2d2.d0(b2) } else { c2d2.d0(a2) };
    let pt1 = if e1.0.orientation().is_forward() { c2d1.d0(a1) } else { c2d1.d0(b1) };
    let start2 = GpPnt::new(pt2.x(), pt2.y(), 0.0);
    let start1 = GpPnt::new(pt1.x(), pt1.y(), 0.0);
    let proj1 = project_inside(&ad1, &start2, tolerance, false);
    let proj2 = project_inside(&ad2, &start1, tolerance, false);
    if proj1.distance > tolerance && proj2.distance > tolerance {
        return None;
    }

    // `cxx:1970-1986`: the longer edge is the one the other is projected onto.
    let (long_ad, short_ad, len_p, first_p, short_num, param) = if proj1.distance < proj2.distance
    {
        (&ad1, &ad2, b2 - a2, a2, n2, proj1.param)
    } else {
        (&ad2, &ad1, b1 - a1, a1, n1, proj2.param)
    };

    // `cxx:1988-1998`: 22 interior samples of the short edge must stay on the
    // long one. The last `Project` default is `AdjustToEnds = true` (`cxx:208`).
    let step = len_p / 23.0;
    let mut sample = first_p;
    for _ in 1..23 {
        let p = short_ad.d0(sample);
        if project_adaptor(long_ad, &p, tolerance, true).distance > tolerance {
            return None;
        }
        sample += step;
    }
    Some(NotchedCheck { short_num, param })
}

/// `ShapeFix_Wire::CopyReversePcurves` (`ShapeFix_Wire.cxx:4108-4176`),
/// file-static in the cxx: append to `to` every `CurveOnSurface` representation
/// of `from` that `to` does not already carry for the same surface *and*
/// location (`cxx:4125-4143`), reversing the first pcurve and its range when
/// `reverse` (`cxx:4152-4160`). The representation's location participates in
/// the identity test, and locations are identity throughout this port, so the
/// face key alone decides it.
pub(in crate::shhealing) fn copy_reverse_pcurves(to: &Edge, from: &Edge, reverse: bool) {
    let reg = GeometryRegistry::global();
    let Some(from_geom) = reg.edge_geom(&from.0) else {
        return;
    };
    let Some(mut to_geom) = reg.edge_geom(&to.0) else {
        return;
    };
    for (face_key, from_pcurves) in &from_geom.pcurves {
        if to_geom.pcurves.contains_key(face_key) {
            continue;
        }
        let (fp0, lp0) = from_geom
            .pcurve_ranges
            .get(face_key)
            .copied()
            .unwrap_or((from_geom.first, from_geom.last));
        let mut range = (fp0, lp0);
        let mut copied: Vec<Arc<dyn Curve2d>> = Vec::with_capacity(from_pcurves.len());
        for (i, pc) in from_pcurves.iter().enumerate() {
            if i == 0 && reverse {
                // `cxx:4152-4160`: `ReversedParameter` on the old range, then
                // `Reverse`, then the swap.
                let fp = super::pcurve_ranges::reversed_parameter(pc.as_ref(), fp0);
                let lp = super::pcurve_ranges::reversed_parameter(pc.as_ref(), lp0);
                let mut c = pc.clone_dyn();
                c.reverse();
                copied.push(Arc::from(c));
                range = (lp, fp);
            } else {
                // `cxx:4168-4171`: `PCurve2` is copied as-is.
                copied.push(Arc::from(pc.clone_dyn()));
            }
        }
        to_geom.pcurves.insert(*face_key, copied);
        to_geom.pcurve_ranges.insert(*face_key, range);
    }
    reg.set_edge(&to.0, to_geom);
}

/// `ShapeFix_Wire::FixDummySeam(num)` (`ShapeFix_Wire.cxx:4213-4289`).
///
/// The `Context()` branches (`cxx:4241-4256`, `:4264-4272`) are the port's
/// in-place wire mutations instead of a `ShapeBuild_ReShape` record, exactly as
/// in the other wire passes; `toRemove` is hard-coded false at `cxx:4227`, so
/// the `Context()->Remove` arm is dead.
pub(in crate::shhealing) fn fix_dummy_seam(wire: &mut Wire, num: usize) {
    let nb = wire_edges_nb(wire);
    if nb < 2 || num == 0 || num > nb {
        return;
    }
    let num1 = if num == nb { 1 } else { num + 1 };
    let edges = edges_of_wire(wire);
    let (Some(e1), Some(e2)) = (edges.get(num - 1).cloned(), edges.get(num1 - 1).cloned()) else {
        return;
    };
    // `cxx:4219-4221`: `V1 = sae.FirstVertex(E1)`, `V2 = sae.LastVertex(E2)`,
    // `Vm = sbv.CombineVertex(V1, V2, 1.0001)` (`CumOri = false`).
    // `ShapeAnalysis_Edge::FirstVertex/LastVertex` are the cumulated
    // (`CumOri = true`) accessors, i.e. orientation-aware, so the port must use
    // `first_vertex` / `last_vertex`, not the raw stored children.
    let (Some(v1), Some(v2)) = (first_vertex(&e1), last_vertex(&e2)) else {
        return;
    };
    let vm = combine_vertex(&v1, &v2, 1.0001);

    // `cxx:4230-4233`: `Vs = sae.FirstVertex(E2)`, replaced by `Vm` when it is
    // already one of the two merged vertices.
    let mut vs = first_vertex(&e2);
    if let Some(vs_ref) = vs.as_ref() {
        if is_same(&vs_ref.0, &v1.0) || is_same(&vs_ref.0, &v2.0) {
            vs = Some(vm.clone());
        }
    }
    let new_edge = copy_replace_vertices_with(&e2, vs.as_ref(), Some(&vm));
    copy_reverse_pcurves(&new_edge, &e1, e1.0.orientation() == e2.0.orientation());
    let reg = GeometryRegistry::global();
    reg.set_same_range(&new_edge.0, false); // `cxx:4237-4238`
    reg.set_same_parameter(&new_edge.0, false);

    // `cxx:4258-4272`: the neighbouring edges are re-pointed at `Vm`.
    let next = if num1 == nb { 1 } else { num1 + 1 };
    let prev = if num > 1 { num - 1 } else { nb };
    if let Some(prev_e) = edges.get(prev - 1).cloned() {
        let tmp = copy_replace_vertices_with(&prev_e, None, Some(&vm));
        wire_set_edge_composed(wire, prev, &tmp);
    }
    if let Some(next_e) = edges.get(next - 1).cloned() {
        let tmp = copy_replace_vertices_with(&next_e, Some(&vm), None);
        wire_set_edge_composed(wire, next, &tmp);
    }

    // `cxx:4275-4288`: both notch edges leave the wire, the higher index first.
    let (n1, n2) = if num < num1 { (num, num1) } else { (num1, num) };
    wire_remove_edge(wire, n2);
    wire_remove_edge(wire, n1);
}

/// `ShapeFix_Wire::FixNotchedEdges()` (`ShapeFix_Wire.cxx:3977-4310`).
///
/// `Perform` calls it at `cxx:400` under
/// `myFixTailMode <= 0 && NeedFix(myFixNotchedEdgesMode, ReorderOK)`.
/// `FromSTEP.FixShape.FixTailMode` is "0" (`STEPControl_Controller.cxx:249`) so
/// the first half holds, and `FixNotchedEdgesMode` is "-1"
/// (`STEPControl_Controller.cxx:248`), so `NeedFix(-1, ReorderOK)`
/// (`ShapeFix_Root.lxx:101-104`) is ReorderOK.
///
/// `WireData()`/`UpdateWire()` (`cxx:3989-3991`, `:4083-4086`) are the port's
/// in-memory wire; `myStatusNotches` (`cxx:4102`) is a status no caller of this
/// port reads.
///
/// Wiring this pass moved no count: a temporary probe on `check_notched_edges`'s
/// `Some` return (since removed) fired 0 times on a full `export_data_obj` run,
/// i.e. `ShapeAnalysis_Wire::CheckNotchedEdges` finds no notched pair in this
/// corpus, exactly as it does in OCCT.
pub(in crate::shhealing) fn fix_notched_edges(wire: &mut Wire, face: &Face, min_tol: f64, max_tol: f64) -> bool {
    let mut done = false;
    let mut i = 1usize;
    // `cxx:3994`: `i <= NbEdges() && NbEdges() > 2`, both re-read every pass.
    while i <= wire_edges_nb(wire) && wire_edges_nb(wire) > 2 {
        let Some(check) = check_notched_edges(wire, face, i, min_tol) else {
            i += 1;
            continue;
        };
        // `cxx:3999-4007`.
        let nb = wire_edges_nb(wire);
        let n2 = if i > 0 { i } else { nb };
        let n1 = if n2 > 1 { n2 - 1 } else { nb };
        let is_remove_first = n1 == check.short_num;
        let to_split = if n2 == check.short_num { n1 } else { n2 };
        let edges = edges_of_wire(wire);
        let Some(split_e) = edges.get(to_split - 1).cloned() else {
            i += 1;
            continue;
        };
        let Some((c2d, a, b)) = curve_on_surface_oriented(&split_e, face, true) else {
            i += 1;
            continue;
        };
        let orient = split_e.0.orientation();

        // `cxx:4011-4018`: the split point falls on an end of the edge; for a
        // closed edge it may fall on the other end (issue #0029780).
        let on_end = (check.param - if is_remove_first { b } else { a }).abs() <= PCONFUSION
            || (edge_is_closed_3d(&split_e)
                && (check.param - if is_remove_first { a } else { b }).abs() <= PCONFUSION);
        if on_end {
            fix_dummy_seam(wire, n1);
            // `cxx:4019-4021`: the seam edge left the list; the cxx's `i--`
            // cancels the loop's `i++`, so `i` stays put.
        } else {
            // `cxx:4027-4029`: the split point is already the edge's own end.
            if ((if is_remove_first { a } else { b }) - check.param).abs() < PCONFUSION {
                i += 1;
                continue;
            }
            // `cxx:4031-4035`.
            let mut transfer = TransferParametersProj::new();
            transfer.set_max_tolerance(max_tol);
            transfer.init(&split_e, face);
            let (first, last) = if a < b { (a, b) } else { (b, a) };

            // `cxx:4046-4050`: `B.MakeVertex(Vnew,
            // Analyzer()->Surface()->Value(c2d->Value(param)),
            // Precision::Confusion())`.
            let Some(surface) = BRepTool::face_surface(face) else {
                i += 1;
                continue;
            };
            let uv = c2d.d0(check.param);
            let vnew = TopoBuilder::new().make_vertex(surface.d0(uv.x(), uv.y()), CONFUSION);

            // `cxx:4051-4058`: first half, from the edge start to the notch.
            let mut we = split_e.clone();
            we.0.set_orientation(Orientation::Forward);
            let (Some(fv), Some(lv)) = (first_vertex(&we), last_vertex(&we)) else {
                i += 1;
                continue;
            };
            let reg = GeometryRegistry::global();
            let mut new_e1 = copy_replace_vertices_with(&we, Some(&fv), Some(&vnew));
            copy_pcurves(&new_e1, &we);
            transfer.transfer_range(&mut new_e1, first, check.param, true);
            reg.set_same_range(&new_e1.0, false);
            reg.set_same_parameter(&new_e1.0, false);
            // `cxx:4059-4067`: second half, from the notch to the edge end.
            let mut new_e2 = copy_replace_vertices_with(&we, Some(&vnew), Some(&lv));
            copy_pcurves(&new_e2, &we);
            transfer.transfer_range(&mut new_e2, check.param, last, true);
            reg.set_same_range(&new_e2.0, false);
            reg.set_same_parameter(&new_e2.0, false);

            // `cxx:4078-4084`: both halves take the split edge's orientation,
            // and a REVERSED edge swaps them.
            new_e1.0.set_orientation(orient);
            new_e2.0.set_orientation(orient);
            if orient == Orientation::Reversed {
                std::mem::swap(&mut new_e1, &mut new_e2);
            }
            // `cxx:4087-4088`: `Set(newE1, toSplit)` then `Add(newE2, ...)`
            // (0 = append).
            let nb_now = wire_edges_nb(wire);
            wire_set_edge_composed(wire, to_split, &new_e1);
            let at = if to_split == nb_now { 0 } else { to_split + 1 };
            let mut ins = new_e2.clone();
            if wire.0.orientation() == Orientation::Reversed {
                ins.0.reverse();
            }
            wire_insert_edge_before(wire, at, &ins);
            // `cxx:4089`: `FixDummySeam(isRemoveLast ? NbEdges() : toRemove)`.
            let is_remove_last = n1 == nb_now && n2 == 1;
            let target = if is_remove_last { wire_edges_nb(wire) } else { check.short_num };
            fix_dummy_seam(wire, target);
        }
        // `cxx:4099`: `DONE1` on both arms. The cxx's `i--` after each fix
        // cancels the loop's `i++`, so `i` is left unchanged here; the `Set` /
        // `Add` above shift the list under it.
        done = true;
    }
    done
}

/// STEP post-pass after pcurve association.
/// Source: `TranslateEdgeLoop::CheckPCurves` first check, then
/// `ShapeFix_Wire::Perform` FixReorder (`cxx:317-325` / `cxx:487`) then
/// `FromSTEP.FixShape.FixShiftedMode` (default on) / `cxx:939` FixShifted,
/// then SameRange / CheckPCurveRange / FixAddPCurve / FixSameParameter
/// (`ShapeFix_Wire.cxx:938-995`), then `ShapeFix_Wire::FixDegenerated`
/// (`cxx:386-392`, `fix_degenerated_all` below).
///
/// `TranslateEdgeLoop.cxx:105-112`: for a plane face `CheckPCurves` calls
/// `RemovePCurves(aWire, aFace)` and returns, so the pcurve-range healing does
/// not run on a plane's edges and their tolerances keep the imported values
/// (linkrods FACE 35 edge 4 keeps 1.2432146551310009e-07). That removal is
/// ported in `check_pcurve_rep_range` below (`cxx:110-114`), together with the
/// non-planar `XSAlgo_ShapeProcessor::CheckPCurve` call at `cxx:175`. The
/// `ShapeFix_Wire` passes that follow are independent of the face type.
pub fn check_pcurves_and_shift(wire: &mut Wire, face: &Face, preci: f64, fix_lacking: bool) {
    check_pcurve_rep_range(wire, face, preci);
    // `ShapeFix_Wire::Perform` (`cxx:317-325`): FixReorder before FixEdgeCurves.
    let mut reorder_ok = fix_reorder_wire(wire, face);
    // `ShapeFix_Wire::Perform` (`cxx:333-347`): FixSmall under
    // `NeedFix(myFixSmallMode, myTopoMode)`. `FromSTEP.FixShape.FixSmallMode` is
    // "-1" (`STEPControl_Controller.cxx:233`) and `myTopoMode` is forced true by
    // `ShapeFix_Shape.cxx:200` for the face path, so the pass runs, with
    // `lockvtx = !myTopoMode || !ReorderOK` and `precsmall = MinTolerance()`.
    if fix_small_all(wire, face, preci, SHAPE_FIX_MIN_TOLERANCE, !reorder_ok) {
        // `cxx:339-343`: after a removal the wire may reorder better; the retry
        // is gated by `NeedFix(myFixReorderMode = -1, !ReorderOK)` evaluated with
        // the ReorderOK from before FixSmall.
        if !reorder_ok {
            reorder_ok = fix_reorder_wire(wire, face);
        }
    }
    // `ShapeFix_Wire::Perform` (`cxx:353-357`): FixConnected after FixSmall.
    if reorder_ok {
        let _ = fix_connected_all(wire, preci, SHAPE_FIX_MAX_TOLERANCE);
    }
    let edges = edges_of_wire(wire);
    // `ShapeFix_Wire::FixEdgeCurves` (`cxx:652-670`): FixAddPCurve runs for
    // every edge of the wire. `FromSTEP.exec.op = FixShape`
    // (`STEPControl_Controller.cxx:201`) routes every imported face through
    // `ShapeFix_Shape` -> `ShapeFix_Face` -> `ShapeFix_Wire::Perform`, so a
    // face whose STEP edges carry a 3D curve only (no `SURFACE_CURVE` /
    // `PCURVE`, as in ATU01038) is projected onto the surface here. With no
    // pcurve stored the face keeps the natural surface bounds
    // (`BRepTools::AddUVBounds`, `BRepTools.cxx:139-153`), which are infinite
    // for a cylinder/cone. The `cxx:673-863` sub-branch (edge over a
    // singularity) needs a `ShapeBuild_ReShape` context and stays unported.
    // `ShapeFix_Edge::myProjector` lives on the `ShapeFix_Edge` that
    // `ShapeFix_Wire` keeps for the whole wire (`ShapeFix_Wire.cxx:280-300`), so
    // one `myCache` covers every edge of this face's wire.
    let mut cache = crate::pcurve_full::ProjectorCache::new();
    for (i, e) in edges.iter().enumerate() {
        let is_seam = is_seam_use(&edges, i);
        let _ = fix_add_pcurve(e, face, is_seam, preci, &mut cache);
    }
    // `cxx:368-374` turns FixShifted off when the reorder reported FAIL.
    if reorder_ok {
        let _ = fix_shifted_wire(wire, face);
    }
    let reg = GeometryRegistry::global();
    for (i, e) in edges.iter().enumerate() {
        if let Some((c2d, fp2d, lp2d)) = curve_on_surface_oriented(e, face, false) {
            let (first, last) = BRepTool::edge_parameters(e);
            if first.is_finite() && last.is_finite() {
                let (cfp, clp) = (c2d.first_parameter(), c2d.last_parameter());
                if (first - fp2d).abs() > PCONFUSION || (last - lp2d).abs() > PCONFUSION {
                    reg.set_same_range(&e.0, false);
                } else if !check_pcurve_range(first, last, c2d.as_ref()) {
                    // `ShapeFix_Wire.cxx:970-985`: when the pcurve cannot cover the
                    // edge's 3D range OCCT does exactly two things —
                    // `ShapeBuild_Edge().RemovePCurve(E, S, L)` followed by
                    // `myFixEdge->FixAddPCurve(E, face, isSeam, myAnalyzer->Surface(),
                    // Precision())`, i.e. it *re-projects* the pcurve from the 3D
                    // curve. There is no range-remapping branch here, so neither is
                    // there one in the port: a stored pcurve whose 2D points fall
                    // outside the surface domain (e.g. the spherical pcurve of
                    // `data/Offset.step` `#1290`, whose line sits at v = -3*pi/4,
                    // outside v in [-pi/2, pi/2]) must be replaced by the projection,
                    // not re-parameterised.
                    let is_seam = is_seam_use(&edges, i);
                    reg.remove_pcurves_on_surface(&e.0, &face.0);
                    let _ = fix_add_pcurve(e, face, is_seam, preci, &mut cache);
                }
            }
        }
        fix_same_parameter(e, face);
    }
    // `ShapeFix_Wire::Perform` (`cxx:386-392`): FixDegenerated runs right after
    // FixEdgeCurves (the loop above) and before FixNotchedEdges / FixTails /
    // FixSelfIntersection / FixLacking / FixVertexTolerance.
    fix_degenerated_all(wire, face, preci);
    // `ShapeFix_Wire::Perform` (`cxx:400`): FixNotchedEdges under
    // `myFixTailMode <= 0 && NeedFix(myFixNotchedEdgesMode, ReorderOK)`.
    // `FromSTEP.FixShape.FixTailMode` is "0" (`STEPControl_Controller.cxx:249`)
    // and `FixNotchedEdgesMode` is "-1" (`STEPControl_Controller.cxx:248`), so
    // the gate is ReorderOK. `cxx:402-406`: a fix is followed by FixShifted.
    if reorder_ok && fix_notched_edges(wire, face, SHAPE_FIX_MIN_TOLERANCE, SHAPE_FIX_MAX_TOLERANCE)
    {
        let _ = fix_shifted_wire(wire, face);
    }
    // `ShapeFix_Wire::Perform` (`cxx:412`): FixTails under
    // `myFixTailMode != 0`. `FromSTEP.FixShape.FixTailMode` is "0"
    // (`STEPControl_Controller.cxx:249`), so OCCT skips this pass in the STEP
    // read-in and porting it here would deviate. Intentionally absent.
    // `ShapeFix_Wire::Perform` (`cxx:428`): FixSelfIntersection under
    // `NeedFix(myFixSelfIntersectionMode, myClosedMode)`. `FixSelfIntersectionMode`
    // is "-1" (`STEPControl_Controller.cxx:238`), so `NeedFix(-1, myClosedMode)`
    // (`ShapeFix_Root.lxx:101-104`) is `myClosedMode` (`ShapeFix_Wire.cxx:171`)
    // and the pass is requested, including its non-adjacent tail
    // (`ShapeFix_Wire.cxx:1200-1223`) with
    // `FromSTEP.FixShape.FixNonAdjacentIntersectingEdgesMode` = "-1"
    // (`STEPControl_Controller.cxx:254-257`) and the `maxtol = 1.0` default
    // (`ShapeFix_IntersectionTool.hxx:46`).
    let _ = fix_self_intersection(wire, face, preci);
    // `ShapeFix_Wire::Perform` (`cxx:448`): FixLacking under
    // `NeedFix(myFixLackingMode, ReorderOK)`; `NeedFix(flag, need)` is
    // `flag < 0 ? need : flag > 0` (`ShapeFix_Root.lxx:101-104`), and
    // `FromSTEP.FixShape.FixLackingMode` is "-1"
    // (`STEPControl_Controller.cxx:237`), so OCCT runs the pass exactly when
    // ReorderOK is true.
    // `ShapeFix_Face::Perform` (`ShapeFix_Face.cxx:369-374`) *disables*
    // `FixLacking` for its first wire round and restores it (`cxx:457`) only
    // for the post-`FixMissingSeam` second round (`cxx:509+`). This reader-time
    // pass is that first round (`TranslateEdgeLoop::CheckPCurves` +
    // `ShapeFix_Face::Perform.cxx:365-480`), so `FixLacking` must not run
    // here: on a periodic face whose bounds are single closed seam edges it
    // inserted a spurious back-going edge, which made `FixMissingSeam`'s
    // `check_wire` sum to 0 and never merge the two wires (see
    // specs/_a3n00_gap_analysis.md 9.95). The second-round call belongs after
    // `FixMissingSeam`, which is not wired into the reader yet.
    if reorder_ok && fix_lacking {
        let _ = fix_lacking_all(wire, face, false, preci, SHAPE_FIX_MAX_TOLERANCE);
    }
}
