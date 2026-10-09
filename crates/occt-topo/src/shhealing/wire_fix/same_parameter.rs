use super::*;

/// True when `check_pcurve_rep_range` runs the `cxx:175`
/// `XSAlgo_ShapeProcessor::CheckPCurve` call. ENABLED (t336); see the call site
/// for the on/off measurements that replaced the earlier block on it.
pub(in crate::shhealing) const T312_WIRE_XSALGO_CHECKPCURVE: bool = true;

/// `Extrema_LocateExtPC` Newton walk (`int_tools_vertex_line` / ValidateEdge).
/// Returns square distance at the local extremum, or `None` when `!IsDone`.
pub(in crate::shhealing) fn locate_ext_pc_sq<D0, D1>(
    d0: D0,
    d1: D1,
    p: &GpPnt,
    t_seed: f64,
    first: f64,
    last: f64,
) -> Option<f64>
where
    D0: Fn(f64) -> GpPnt,
    D1: Fn(f64) -> (GpPnt, GpVec),
{
    let lo = first.min(last);
    let hi = first.max(last);
    let mut t = t_seed.clamp(lo, hi);
    let mut last_good = d0(t);
    for _ in 0..12 {
        let (q, der) = d1(t);
        let speed2 = der.square_magnitude();
        if speed2 <= 1e-20 {
            return None;
        }
        let num = (q.x() - p.x()) * der.x() + (q.y() - p.y()) * der.y() + (q.z() - p.z()) * der.z();
        let next = t - num / speed2;
        if !next.is_finite() {
            return None;
        }
        let next = next.clamp(lo, hi);
        let qn = d0(next);
        if qn.distance(&last_good) < 1e-10 {
            return Some(qn.square_distance(p));
        }
        t = next;
        last_good = qn;
    }
    Some(last_good.square_distance(p))
}

pub(in crate::shhealing) fn cos_d0(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> GpPnt {
    let uv = pc.d0(t);
    surf.d0(uv.x(), uv.y())
}

pub(in crate::shhealing) fn cos_d1(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> (GpPnt, GpVec) {
    let (uv, duv) = pc.d1(t);
    let (p, su, sv) = surf.d1(uv.x(), uv.y());
    let tan = su
        .multiplied_scalar(duv.x())
        .added(&sv.multiplied_scalar(duv.y()));
    (p, tan)
}

/// `BRepLib_ValidateEdge::processApprox` (`cxx:103-229`) plus
/// `UpdateTolerance` (`cxx:57-63`).
/// `ShapeAnalysis_Edge::CheckSameParameter` (`cxx:704-838`).
pub(in crate::shhealing) fn check_same_parameter(edge: &Edge, face: &Face) -> f64 {
    if BRepTool::is_degenerated(edge) {
        return 0.0;
    }
    let Some(c3d) = BRepTool::edge_curve(edge) else {
        return 0.0;
    };
    let (first, last) = BRepTool::edge_parameters(edge);
    if !first.is_finite() || !last.is_finite() {
        return 0.0;
    }
    let Some(surf) = BRepTool::face_surface(face) else {
        return 0.0;
    };
    let same_p = BRepTool::same_parameter(edge);
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    let is_plane = classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane;
    // After CheckPCurves RemovePCurves (`cxx:110-114`) a plane has no stored
    // pcurve. Do not same-t leftover STEP pcs (Shape-2 1.322 inflation).
    let mut pcs = if is_plane {
        Vec::new()
    } else {
        reg.edge_pcurves(&edge.0, face_key)
    };
    let (mut cf, mut cl) = reg
        .pcurve_range(&edge.0, face_key)
        .unwrap_or((first, last));
    // `ShapeAnalysis_Edge.cxx:801-826`: `CurveOnPlane` then
    // `Geom2dAdaptor_Curve(aPC, aFirst, aLast)` + `CurveOnSurface` +
    // `BRepLib_ValidateEdge` (`aFirst/aLast` = 3D edge range, cxx:418-419).
    if pcs.is_empty() && is_plane {
        if let Ok(pc) = crate::pcurve::make_pcurve_on_face(edge, face) {
            pcs.push(pc);
            cf = first;
            cl = last;
        }
    }
    let mut max_sq: f64 = 0.0;
    // `SetControlPointsNumber(NbControl - 1)` with default NbControl=23.
    const NB: i32 = 22;
    for pc in &pcs {
        let projection = !same_p
            || (cf - first).abs() > PCONFUSION
            || (cl - last).abs() > PCONFUSION;

        if projection {
            // `BRepLib_ValidateEdge.cxx:141-227`.
            let saved = max_sq;
            let p_ref0 = c3d.d0(first);
            let p_oth0 = cos_d0(pc.as_ref(), surf.as_ref(), cf);
            max_sq = max_sq.max(p_ref0.square_distance(&p_oth0));
            let p_ref1 = c3d.d0(last);
            let p_oth1 = cos_d0(pc.as_ref(), surf.as_ref(), cl);
            max_sq = max_sq.max(p_ref1.square_distance(&p_oth1));
            let mut aborted = false;
            for i in 1..NB {
                let t_ref = ((NB - i) as f64 * first + i as f64 * last) / NB as f64;
                let t_oth = ((NB - i) as f64 * cf + i as f64 * cl) / NB as f64;
                let p_ref = c3d.d0(t_ref);
                let p_oth = cos_d0(pc.as_ref(), surf.as_ref(), t_oth);
                let Some(sq_ref) = locate_ext_pc_sq(
                    |t| c3d.d0(t),
                    |t| c3d.d1(t),
                    &p_oth,
                    t_ref,
                    first,
                    last,
                ) else {
                    aborted = true;
                    break;
                };
                max_sq = max_sq.max(sq_ref);
                let Some(sq_oth) = locate_ext_pc_sq(
                    |t| cos_d0(pc.as_ref(), surf.as_ref(), t),
                    |t| cos_d1(pc.as_ref(), surf.as_ref(), t),
                    &p_ref,
                    t_oth,
                    cf,
                    cl,
                ) else {
                    aborted = true;
                    break;
                };
                max_sq = max_sq.max(sq_oth);
            }
            if aborted {
                // cxx:202-205 / 222-225: !IsDone leaves CalculatedDistance unset
                // for this pair; earlier pairs already contributed to maxdev.
                max_sq = saved;
            }
            continue;
        }
        for i in 0..=NB {
            let t = ((NB - i) as f64 * first + i as f64 * last) / NB as f64;
            let uv = pc.d0(t);
            let sq = c3d.d0(t).square_distance(&surf.d0(uv.x(), uv.y()));
            max_sq = max_sq.max(sq);
        }
    }
    let maxdev = max_sq.sqrt() * 1.00001;
    maxdev
}

/// `ShapeAnalysis_Edge::CheckVertexTolerance` with `checkAll=true`
/// (`cxx:568-667`) then `B.UpdateVertex`.
pub(in crate::shhealing) fn fix_vertex_tolerance_edge(edge: &Edge) {
    let (Some(v1), Some(v2)) = (first_vertex(edge), last_vertex(edge)) else {
        return;
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    let mut t1: f64 = 0.0;
    let mut t2: f64 = 0.0;
    if let Some(c3d) = BRepTool::edge_curve(edge) {
        let (mut a, mut b) = BRepTool::edge_parameters(edge);
        if a.is_finite() && b.is_finite() {
            if edge.0.orientation().is_reversed() {
                std::mem::swap(&mut a, &mut b);
            }
            t1 = p1.square_distance(&c3d.d0(a));
            t2 = p2.square_distance(&c3d.d0(b));
        }
    }
    let reg = GeometryRegistry::global();
    for (surf, pc, mut a, mut b) in reg.edge_pcurve_reps(&edge.0) {
        if !(a.is_finite() && b.is_finite()) {
            continue;
        }
        if edge.0.orientation().is_reversed() {
            std::mem::swap(&mut a, &mut b);
        }
        let uv1 = pc.d0(a);
        let uv2 = pc.d0(b);
        t1 = t1.max(p1.square_distance(&surf.d0(uv1.x(), uv1.y())));
        t2 = t2.max(p2.square_distance(&surf.d0(uv2.x(), uv2.y())));
    }
    let tole = BRepTool::edge_tolerance(edge);
    v1.set_tolerance((1.0000001 * t1.sqrt()).max(tole));
    v2.set_tolerance((1.0000001 * t2.sqrt()).max(tole));
}

/// `ShapeBuild_Edge::CopyReplaceVertices` (`ShapeBuild_Edge.cxx:59-140`): a new
/// edge over the same 3D curve and location, carrying the source vertices
/// (`cxx:99-119`, first FORWARD and last REVERSED) and the source
/// tolerance/flags/ranges/pcurves (`BRep_TEdge::EmptyCopy`, `BRep_TEdge.cxx:104-129`,
/// appends a `Copy()` of every `BRep_GCurve` including its pcurve).
/// `EmptyCopied` also keeps the source orientation (`TopoDS_Shape.hxx:294-302`)
/// and `CopyRanges` (`cxx:125`) re-ranges every representation.
pub fn copy_replace_vertices(edge: &Edge) -> Edge {
    copy_replace_vertices_with(edge, None, None)
}

/// `ShapeBuild_Edge::CopyReplaceVertices(edge, V1, V2)`
/// (`ShapeBuild_Edge.cxx:59-140`). A null argument keeps the edge's own
/// endpoint (`cxx:65-96`), which is the `TopExp` FORWARD / REVERSED vertex
/// (`cxx:78-93`) and, in this port, `first_vertex` / `last_vertex`.
pub fn copy_replace_vertices_with(edge: &Edge, v1: Option<&Vertex>, v2: Option<&Vertex>) -> Edge {
    let reg = GeometryRegistry::global();
    let Some(geom) = reg.edge_geom(&edge.0) else {
        return edge.clone();
    };
    let builder = TopoBuilder::new();
    let mut copy = builder.make_edge(geom.curve.clone(), geom.first, geom.last);
    // `TopoDS_Shape::EmptyCopied` (`TopoDS_Shape.hxx:294-302`) keeps the source
    // orientation, and `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:104-129`) keeps
    // the tolerance, the `SameParameter` / `SameRange` / `Degenerated` flags and
    // a `Copy()` of every curve representation - the 3D curve and the pcurves.
    // Registering the source `geom` reproduces all of that (see `EdgeGeom`).
    //
    // Order matters: `ShapeBuild_Edge.cxx:103` runs `edge.EmptyCopied()`, which
    // keeps the source orientation and location (`TopoDS_Shape.hxx:297-302`),
    // *before* the two `B.Add(E, V)` calls (`cxx:107-114`).
    // `TopoDS_Builder::Add` (`TopoDS_Builder.cxx:74-91`) composes E's own
    // orientation into each child: for a REVERSED E the stored children are V1
    // with REVERSED and V2 with FORWARD storage orientation.
    // `TopExp::Vertices(E, V1, V2, true)` (`TopExp.cxx:214-253`) iterates with
    // `CumOri = true` (`TopoDS_Iterator.cxx:26-83`, `Compose` at :73-82), so the
    // cumulated orientation is V1 FORWARD / V2 REVERSED again and
    // `TopExp::FirstVertex` is V1, `TopExp::LastVertex` is V2 - independent of
    // E's orientation. Setting the orientation and the location first reproduces
    // exactly that; `first_vertex` / `last_vertex` below are the cumulated
    // accessors.
    copy.0.set_orientation(edge.0.orientation());
    copy.0.set_location(edge.0.location());
    let fv = v1.cloned().or_else(|| first_vertex(edge));
    let lv = v2.cloned().or_else(|| last_vertex(edge));
    if let (Some(v1), Some(v2)) = (fv, lv) {
        builder.add_edge_vertices(&mut copy, &v1, &v2);
    }
    reg.set_edge(&copy.0, geom);
    copy
}

/// `BRepTools::Compare(V1, V2)` (`BRepTools.cxx:527-545`): same vertex, or a
/// gap no larger than either vertex tolerance.
pub(crate) fn vertices_coincide(a: &Vertex, b: &Vertex) -> bool {
    if is_same(&a.0, &b.0) {
        return true;
    }
    let l = BRepTool::vertex_point(a).distance(&BRepTool::vertex_point(b));
    l <= BRepTool::vertex_tolerance(a) || l <= BRepTool::vertex_tolerance(b)
}

/// `ShapeBuild_Vertex::CombineVertex(V1, V2, tolFactor)`
/// (`ShapeBuild_Vertex.cxx:26-71`).
pub fn combine_vertex(v1: &Vertex, v2: &Vertex, tol_factor: f64) -> Vertex {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    let tol1 = BRepTool::vertex_tolerance(v1);
    let tol2 = BRepTool::vertex_tolerance(v2);
    let dist = p1.distance(&p2);
    let (pos, tol) = if dist + tol2 <= tol1 {
        (p1, tol1)
    } else if dist + tol1 <= tol2 {
        (p2, tol2)
    } else {
        let tol = 0.5 * (dist + tol1 + tol2);
        // `ShapeBuild_Vertex.cxx:64`: guarded against a zero distance.
        let s = if dist > 0.0 { (tol2 - tol1) / dist } else { 0.0 };
        (
            GpPnt::new(
                0.5 * ((1.0 - s) * p1.x() + (1.0 + s) * p2.x()),
                0.5 * ((1.0 - s) * p1.y() + (1.0 + s) * p2.y()),
                0.5 * ((1.0 - s) * p1.z() + (1.0 + s) * p2.z()),
            ),
            tol,
        )
    };
    TopoBuilder::new().make_vertex(pos, tol_factor * tol)
}

/// `ShapeBuild_Edge::CopyPCurves` (`ShapeBuild_Edge.cxx:360-413`): copy every
/// `CurveOnSurface` representation of `from` onto `to` (replacing the
/// representation on the same surface, `cxx:376-398`) together with its
/// pcurve `Copy()` and range (`cxx:400-411`).
pub fn copy_pcurves(to: &Edge, from: &Edge) {
    // `ShapeBuild_Edge::CopyPCurves` (`ShapeBuild_Edge.cxx:360-413`) walks the
    // `BRep_TEdge` CurveRepresentation list, which is independent of the 3D
    // curve. Copying through `edge_geom` dropped the pcurves of curve-less
    // edges (the `ShapeFix_ComposeShell::SplitByLine` seam edges,
    // `ShapeFix_ComposeShell.cxx:2061-2070`), so dispatch produced a face
    // whose seam edges had no pcurve and never meshed.
    GeometryRegistry::global().copy_edge_pcurve_slots(&to.0, &from.0);
}

/// `ShapeFix_Edge::FixSameParameter` (`ShapeFix_Edge.cxx:798-935`).
///
/// Both arms are ported: `wasSP` (`cxx:872-928`) and `!wasSP`
/// (`cxx:839-858` copy + `BRepLib::SameParameter(copy, tol)` with its
/// `Approx_SameParameter` re-projection at `BRepLib.cxx:1631`, then the
/// pick-best block `cxx:889-912` with `ShapeBuild_Edge::CopyPCurves`).
///
/// UNPORTED, inside the `!wasSP` arm: the cxx widens the deviation check to
/// *all* pcurves by passing an empty face (`cxx:877-881`) and
/// `BRepLib::SameParameter(edge, tol)` walks every `CurveOnSurface`
/// representation (`BRepLib.cxx:1263-1291`); our `check_same_parameter` and
/// `same_parameter_edge` are face based, so they cover the call-site face only.
///
/// Also unmodelled in both arms: the `ShapeExtend` status writes of the method
/// (`cxx:853-855` FAIL2, `cxx:883-886` FAIL1, `cxx:894` DONE3, `cxx:902-905`,
/// `cxx:910` DONE5, `cxx:931-934` DONE2); this port has no `myStatus`
/// accumulator.
pub(in crate::shhealing) fn fix_same_parameter(edge: &Edge, face: &Face) {
    let reg = GeometryRegistry::global();
    let same_range = reg.edge_geom(&edge.0).map(|g| g.same_range).unwrap_or(true);
    if reg.is_degenerated_edge(&edge.0) {
        if !same_range {
            crate::shhealing::temp_same_range(edge);
        }
        reg.set_same_parameter(&edge.0, true);
        return;
    }
    let tol = BRepTool::edge_tolerance(edge);
    let (first_3d, last_3d) = BRepTool::edge_parameters(edge);
    // `cxx:820-823`: `V1`/`V2` and their tolerances are read *before*
    // `BRepLib::SameParameter` runs. The copy made by the `!wasSP` arm carries
    // these same two vertices (`ShapeBuild_Edge.cxx:105-115`), and that call
    // widens them (`BRepLib.cxx:1220-1232`), so the restore below needs the
    // values captured here.
    let v1 = first_vertex(edge);
    let v2 = last_vertex(edge);
    let tol_fv = v1.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
    let tol_lv = v2.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
    // `cxx:826`: `bool wasSP = BRep_Tool::SameParameter(edge), SP = false`.
    let was_sp = BRepTool::same_parameter(edge);
    let mut copy_edge: Option<Edge> = None;
    let mut sp = false;
    if !same_range {
        crate::shhealing::temp_same_range(edge);
    }
    if !was_sp {
        // `cxx:839-858`: heal on a copy of the edge, keeping the original.
        let copy = copy_replace_vertices(edge); // cxx:841
        copy_pcurves(&copy, edge); // cxx:841 `Copy(edge, false)`
        reg.set_same_parameter(&copy.0, false); // cxx:842
        // `cxx:844-849`: `Copy` may shift a periodic 3D range; enforce the
        // original range so the copy's pcurves are not rewritten under it.
        if let Some(mut geom) = reg.edge_geom(&copy.0) {
            geom.first = first_3d;
            geom.last = last_3d;
            reg.set_edge(&copy.0, geom);
        }
        // `cxx:850`: the two argument overload passes the default
        // `tolerance = 0` (`ShapeFix_Edge.hxx:193`), so `tol` is selected.
        crate::brep_lib_same_parameter::same_parameter_edge_inplace(&copy, face, tol);
        sp = BRepTool::same_parameter(&copy); // cxx:851
        copy_edge = Some(copy);
    }
    reg.set_same_parameter(&edge.0, true); // cxx:872
    // `cxx:883`: deviation on the original pcurves.
    let mut maxdev = check_same_parameter(edge, face);
    if sp {
        if let Some(copy) = copy_edge.as_ref() {
            // `cxx:889-912`: compare and select the best variant.
            let mut brl_tol = BRepTool::edge_tolerance(copy); // cxx:891
            let brl_dev = check_same_parameter(copy, face); // cxx:892
            if brl_tol < brl_dev {
                brl_tol = brl_dev; // cxx:895-898
            }
            if brl_tol < maxdev {
                copy_pcurves(edge, copy); // cxx:906
                maxdev = brl_tol; // cxx:907
                reg.set_edge_tolerance(&edge.0, brl_tol); // cxx:908
            }
        }
    }
    // `cxx:915-922`: restore the vertex tolerances captured at `cxx:820-823`,
    // because `BRepLib::SameParameter` may have widened them through the copy.
    // `ShapeFix_ShapeTolerance::SetTolerance` writes the value *exactly*
    // (`ShapeFix_ShapeTolerance.cxx:142-164`) and ignores a non-positive
    // precision (`cxx:149-152`).
    for (v, tol_v) in [(v1.as_ref(), tol_fv), (v2.as_ref(), tol_lv)] {
        if let Some(v) = v {
            let value = maxdev.max(tol_v);
            if value > 0.0 {
                v.set_tolerance(value);
            }
        }
    }
    if maxdev > tol {
        reg.set_edge_tolerance(&edge.0, maxdev);
        fix_vertex_tolerance_edge(edge);
    }
}
