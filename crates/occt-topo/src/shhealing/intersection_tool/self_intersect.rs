use super::*;

/// `ShapeFix_IntersectionTool::FixSelfIntersectWire`
/// (`cxx:1028-1830`).
///
/// `max_tol` is the tool's `myMaxTol` (`ShapeFix_IntersectionTool.hxx:67`
/// default `1.0`; `ShapeFix_Wire.cxx:1204` constructs the tool without a
/// `maxtol` argument). The tool's `myPreci` is stored only (`cxx:56`) and is
/// never read in this method, so it has no parameter here.
pub fn fix_self_intersect_wire(
    wire: &mut Wire,
    face: &Face,
    max_tol: f64,
    ctx: &mut dyn ReShape,
) -> FixSelfIntersectWireStatus {
    // `cxx:1035-1038`: `if (myContext.IsNull() || face.IsNull()) return false;`
    // The port always hands in a context, and a null face is not representable.
    let mut max_tol_vert = 0.0f64;
    // `cxx:1043-1050`.
    {
        let sf = ctx.apply(&face.0);
        for v in vertices_of(&sf) {
            max_tol_vert = max_tol_vert.max(BRepTool::vertex_tolerance(&v));
        }
        max_tol_vert = max_tol_vert.min(max_tol);
    }

    // `cxx:1055-1059`: step 2 box map over the wire's edges. (`cxx:1051`
    // "step 1" has no code: adjacent-edge intersection is the `FixIntersectingEdges`
    // pair driven by `ShapeFix_Wire::FixSelfIntersection` itself.)
    let mut boxes: Boxes = HashMap::new();
    create_boxes_2d(wire, face, &mut boxes);
    // OCCT hands the same `boxes` map (mutable, in/out) to every helper below;
    // shadow it with a mutable borrow so the calls read like the C++ source.
    let boxes = &mut boxes;

    // `cxx:1061-1064`.
    let mut nb_split = 0i64;
    let mut nb_cut = 0i64;
    let mut nb_replaced = 0i64;
    let mut nb_removed = 0i64;

    let status = |nb_split: i64, nb_cut: i64, nb_removed: i64, done: bool| FixSelfIntersectWireStatus {
        nb_split,
        nb_cut,
        nb_removed,
        done,
    };

    // `cxx:1065-1068`: `for (int num1 = 1; num1 < sewd->NbEdges() && NbSplit < 30; num1++)`.
    let mut num1 = 1i64;
    while num1 < wire_len(wire) && nb_split < 30 {
        // `cxx:1069`: `for (int num2 = num1 + 2; num2 <= sewd->NbEdges() && NbSplit < 30; num2++)`.
        let mut num2 = num1 + 2;
        loop {
            if !(num2 <= wire_len(wire) && nb_split < 30) {
                break;
            }
            // C++ `continue` advances `num2` by one; `num2--; continue;` nets
            // zero. Each `break 'body` states which one applies through `step2`.
            let mut step2 = 1i64;
            let mut restart1 = false;

            'body: {
                // `cxx:1071-1074`.
                if num1 == 1 && num2 == wire_len(wire) {
                    break 'body;
                }
                let Some(edge1) = wire_edge_at(wire, num1) else {
                    break 'body;
                };
                let Some(edge2) = wire_edge_at(wire, num2) else {
                    break 'body;
                };
                // `cxx:1078-1082`.
                if is_same(&edge1.0, &edge2.0) {
                    break 'body;
                }
                if BRepTool::is_degenerated(&edge1) || BRepTool::is_degenerated(&edge2) {
                    break 'body;
                }
                // `cxx:1085-1092`.
                let (Some(b1_box), Some(b2_box)) =
                    (boxes.get(&skey(&edge1)).copied(), boxes.get(&skey(&edge2)).copied())
                else {
                    break 'body;
                };
                if b1_box.is_out_box(&b2_box) {
                    break 'body;
                }

                // `cxx:1097-1103`: `sae.PCurve` failure aborts the whole call.
                let Some((crv1, a1, b1)) = curve_on_surface_oriented(&edge1, face, false) else {
                    return status(nb_split, nb_cut, nb_removed, false);
                };
                let Some((crv2, a2, b2)) = curve_on_surface_oriented(&edge2, face, false) else {
                    return status(nb_split, nb_cut, nb_removed, false);
                };

                // `cxx:1105-1112`.
                let d1 = IntRes2dDomain::bounded(&crv1.d0(a1), a1, TOLINT, &crv1.d0(b1), b1, TOLINT);
                let d2 = IntRes2dDomain::bounded(&crv2.d0(a2), a2, TOLINT, &crv2.d0(b2), b2, TOLINT);
                let mut inter = Geom2dIntGInter::new();
                inter.perform(crv1.as_ref(), &d1, crv2.as_ref(), &d2, TOLINT, TOLINT);
                if !inter.is_done() {
                    break 'body;
                }

                if inter.nb_points() > 0 && inter.nb_points() < 3 {
                    // `cxx:1115-1119`: intersection is a point.
                    let ip = select_int_pnt(&inter);
                    let tr1 = *ip.transition_of_first();
                    let tr2 = *ip.transition_of_second();

                    // `cxx:1120-1243`.
                    if tr1.position_on_curve() == IntRes2dPosition::Middle
                        && tr2.position_on_curve() == IntRes2dPosition::Middle
                    {
                        let param1 = ip.param_on_first();
                        let param2 = ip.param_on_second();
                        let Some(surf) = BRepTool::face_surface(face) else {
                            return status(nb_split, nb_cut, nb_removed, false);
                        };
                        let pi1 = point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), param1);
                        let pi2 = point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), param2);
                        let mut tol_v = 0.0f64;
                        let mut v_sel: Option<Vertex> = None;

                        // `cxx:1129-1163`: analysis for edge1.
                        let modif_e1 = analyze_point_edge(
                            &edge1, &pi1, a1, b1, param1, max_tol_vert, face,
                            &mut nb_cut, &mut tol_v, &mut v_sel,
                        );
                        // `cxx:1165-1197`: analysis for edge2.
                        let modif_e2 = analyze_point_edge(
                            &edge2, &pi2, a2, b2, param2, max_tol_vert, face,
                            &mut nb_cut, &mut tol_v, &mut v_sel,
                        );

                        // `cxx:1199-1208`.
                        if modif_e1 && !modif_e2 {
                            if let Some(v) = &v_sel {
                                if split_edge1(wire, ctx, face, num2, param2, v, tol_v, boxes)
                                {
                                    nb_split += 1;
                                    step2 = 0; // `num2--; continue;`
                                    break 'body;
                                }
                            }
                        }
                        // `cxx:1209-1217`.
                        if !modif_e1 && modif_e2 {
                            if let Some(v) = &v_sel {
                                if split_edge1(wire, ctx, face, num1, param1, v, tol_v, boxes)
                                {
                                    nb_split += 1;
                                    num1 -= 1; // `num1--; break;`
                                    restart1 = true;
                                    break 'body;
                                }
                            }
                        }
                        // `cxx:1218-1239`.
                        if !modif_e1 && !modif_e2 {
                            let p0 = GpPnt::new(
                                0.5 * (pi1.x() + pi2.x()),
                                0.5 * (pi1.y() + pi2.y()),
                                0.5 * (pi1.z() + pi2.z()),
                            );
                            tol_v = ((pi1.distance(&pi2) / 2.0) * 1.00001).max(CONFUSION);
                            let v = TopoBuilder::new().make_vertex(p0, tol_v);
                            max_tol_vert = max_tol_vert.max(tol_v);
                            let is_edge_split2 =
                                split_edge1(wire, ctx, face, num2, param2, &v, tol_v, boxes);
                            if is_edge_split2 {
                                nb_split += 1;
                                step2 = 0; // `num2--`
                            }
                            if split_edge1(wire, ctx, face, num1, param1, &v, tol_v, boxes) {
                                nb_split += 1;
                                num1 -= 1;
                                restart1 = true;
                                break 'body;
                            }
                            if is_edge_split2 {
                                break 'body; // `continue` with the `num2--` already applied
                            }
                        }
                    }

                    // `cxx:1244-1259`.
                    if tr1.position_on_curve() == IntRes2dPosition::Middle
                        && tr2.position_on_curve() != IntRes2dPosition::Middle
                    {
                        let param1 = ip.param_on_first();
                        if find_vert_and_split_edge(
                            param1,
                            &edge1,
                            &edge2,
                            crv1.as_ref(),
                            &mut max_tol_vert,
                            &mut num1,
                            wire,
                            face,
                            boxes,
                            false,
                            ctx,
                        ) {
                            nb_split += 1;
                            // `num1` was decremented inside, `break` then
                            // `num1++` nets zero.
                            restart1 = true;
                            break 'body;
                        }
                    }
                    // `cxx:1260-1276`.
                    if tr1.position_on_curve() != IntRes2dPosition::Middle
                        && tr2.position_on_curve() == IntRes2dPosition::Middle
                    {
                        let param2 = ip.param_on_second();
                        if find_vert_and_split_edge(
                            param2,
                            &edge2,
                            &edge1,
                            crv2.as_ref(),
                            &mut max_tol_vert,
                            &mut num2,
                            wire,
                            face,
                            boxes,
                            false,
                            ctx,
                        ) {
                            nb_split += 1;
                            // `num2` was decremented inside, `num2++` nets zero.
                            break 'body;
                        }
                    }
                    // `cxx:1278-1285`.
                    if tr1.position_on_curve() != IntRes2dPosition::Middle
                        && tr2.position_on_curve() != IntRes2dPosition::Middle
                    {
                        if union_vertexes(wire, ctx, &edge1, &edge2, num2, boxes, &b2_box) {
                            nb_replaced += 1;
                        }
                    }
                }

                // `cxx:1288-1299`: intersection is a segment.
                if inter.nb_segments() == 1 {
                    let seg = inter.segment(1);
                    if seg.has_first_point() && seg.has_last_point() {
                        let ipf = seg.first_point().clone();
                        let ipl = seg.last_point().clone();
                        let p11 = ipf.param_on_first();
                        let p21 = ipf.param_on_second();
                        let p12 = ipl.param_on_first();
                        let p22 = ipl.param_on_second();
                        let Some(surf) = BRepTool::face_surface(face) else {
                            return status(nb_split, nb_cut, nb_removed, false);
                        };
                        let pnt11 = point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), p11);
                        let pnt12 = point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), p12);
                        let pnt21 = point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), p21);
                        let pnt22 = point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), p22);
                        // `cxx:1311-1314`.
                        if pnt11.distance(&pnt21) > max_tol_vert
                            || pnt12.distance(&pnt22) > max_tol_vert
                        {
                            break 'body;
                        }

                        let mut is_modified1 = false;
                        let mut is_modified2 = false;
                        let mut new_v: Option<Vertex> = None;
                        let mut newtol = 0.0f64;

                        // `cxx:1316-1397`: analysis for edge1.
                        if let (Some(v1), Some(v2)) =
                            (first_vertex(&edge1), last_vertex(&edge1))
                        {
                            let pv1 = vertex_position(&v1);
                            let pv2 = vertex_position(&v2);
                            let dist1 = pnt11.distance(&pv1);
                            let dist2 = pnt12.distance(&pv1);
                            let maxdist = dist1.max(dist2);
                            let pdist = if edge1.0.orientation() == Orientation::Reversed {
                                (b1 - p11).abs().max((b1 - p12).abs())
                            } else {
                                (a1 - p11).abs().max((a1 - p12).abs())
                            };
                            if maxdist < max_tol_vert || pdist < (b1 - a1).abs() * 0.01 {
                                newtol = maxdist;
                                new_v = Some(v1.clone());
                                is_modified1 = true;
                            }
                            let dist1 = pnt11.distance(&pv2);
                            let dist2 = pnt12.distance(&pv2);
                            let maxdist = dist1.max(dist2);
                            let pdist = if edge1.0.orientation() == Orientation::Reversed {
                                (a1 - p11).abs().max((a1 - p12).abs())
                            } else {
                                (b1 - p11).abs().max((b1 - p12).abs())
                            };
                            if maxdist < max_tol_vert || pdist < (b1 - a1).abs() * 0.01 {
                                if (is_modified1 && maxdist < newtol) || !is_modified1 {
                                    newtol = maxdist;
                                    new_v = Some(v2.clone());
                                    is_modified1 = true;
                                }
                            }
                            if is_modified1 {
                                // `cxx:1362-1396`: cut edge1 and widen NewV.
                                let dista = (a1 - p11).abs() + (a1 - p12).abs();
                                let distb = (b1 - p11).abs() + (b1 - p12).abs();
                                let pend = if dista > distb { a1 } else { b1 };
                                let cut = if (pend - p11).abs() > (pend - p12).abs() {
                                    p12
                                } else {
                                    p11
                                };
                                let mut is_cut_line = false;
                                if cut_edge(&edge1, pend, cut, face, &mut is_cut_line) {
                                    nb_cut += 1;
                                }
                                if let Some(nv) = &new_v {
                                    if newtol > BRepTool::vertex_tolerance(nv) {
                                        update_vertex_tolerance(nv, newtol);
                                    } else {
                                        newtol = BRepTool::vertex_tolerance(nv);
                                    }
                                }
                            }
                        }

                        // `cxx:1399-1483`: analysis for edge2.
                        if let (Some(v12), Some(v22)) =
                            (first_vertex(&edge2), last_vertex(&edge2))
                        {
                            let pv12 = vertex_position(&v12);
                            let pv22 = vertex_position(&v22);
                            let dist1 = pnt21.distance(&pv12);
                            let dist2 = pnt22.distance(&pv12);
                            let maxdist = dist1.max(dist2);
                            let pdist = if edge2.0.orientation() == Orientation::Reversed {
                                (b2 - p21).abs().max((b2 - p22).abs())
                            } else {
                                (a2 - p21).abs().max((a2 - p22).abs())
                            };
                            if maxdist < max_tol_vert || pdist < (b2 - a2).abs() * 0.01 {
                                newtol = maxdist;
                                new_v = Some(v12.clone());
                                is_modified2 = true;
                            }
                            let dist1 = pnt21.distance(&pv22);
                            let dist2 = pnt22.distance(&pv22);
                            let maxdist = dist1.max(dist2);
                            let pdist = if edge2.0.orientation() == Orientation::Reversed {
                                (a2 - p21).abs().max((a2 - p22).abs())
                            } else {
                                (b2 - p21).abs().max((b2 - p22).abs())
                            };
                            if maxdist < max_tol_vert || pdist < (b2 - a2).abs() * 0.01 {
                                if (is_modified2 && maxdist < newtol) || !is_modified2 {
                                    newtol = maxdist;
                                    new_v = Some(v22.clone());
                                    is_modified2 = true;
                                }
                            }
                            if is_modified2 {
                                // `cxx:1445-1482`: cut edge2 and widen NewV.
                                let dista = (a2 - p21).abs() + (a2 - p22).abs();
                                let distb = (b2 - p21).abs() + (b2 - p22).abs();
                                let pend = if dista > distb { a2 } else { b2 };
                                let cut = if (pend - p21).abs() > (pend - p22).abs() {
                                    p22
                                } else {
                                    p21
                                };
                                let mut is_cut_line = false;
                                if cut_edge(&edge2, pend, cut, face, &mut is_cut_line) {
                                    nb_cut += 1;
                                }
                                if let Some(nv) = &new_v {
                                    if newtol > BRepTool::vertex_tolerance(nv) {
                                        update_vertex_tolerance(nv, newtol);
                                    } else {
                                        newtol = BRepTool::vertex_tolerance(nv);
                                    }
                                }
                            }
                        }

                        // `cxx:1486-1497`.
                        if is_modified1 && !is_modified2 {
                            if let Some(nv) = &new_v {
                                if split_edge2(
                                    wire, ctx, face, num2, p21, p22, nv, newtol, boxes,
                                ) {
                                    nb_split += 1;
                                    step2 = 0; // `num2--; continue;`
                                    break 'body;
                                }
                            }
                        }
                        // `cxx:1494-1509`.
                        if !is_modified1 && is_modified2 {
                            if let Some(nv) = &new_v {
                                if split_edge2(
                                    wire, ctx, face, num1, p11, p12, nv, newtol, boxes,
                                ) {
                                    nb_split += 1;
                                    num1 -= 1; // `num1--; break;`
                                    restart1 = true;
                                    break 'body;
                                }
                            }
                        }
                        // `cxx:1501-1818`.
                        if !is_modified1 && !is_modified2 {
                            let param1 = 0.5 * (p11 + p12);
                            let param2 = 0.5 * (p21 + p22);
                            let pnt10 = point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), param1);
                            let pnt20 = point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), param2);
                            let p0 = GpPnt::new(
                                0.5 * (pnt10.x() + pnt20.x()),
                                0.5 * (pnt10.y() + pnt20.y()),
                                0.5 * (pnt10.z() + pnt20.z()),
                            );
                            let dist1 = pnt11.distance(&p0).max(pnt12.distance(&p0));
                            let dist2 = pnt21.distance(&p0).max(pnt22.distance(&p0));
                            let mut tol_v = dist1.max(dist2);
                            tol_v = tol_v.max(pnt10.distance(&pnt20)) * 1.00001;
                            let fix_segment = true;
                            if tol_v < max_tol_vert {
                                // `cxx:1515-1529`: create a new vertex and split
                                // each intersecting edge into two.
                                let nv = TopoBuilder::new().make_vertex(pnt10, tol_v);
                                if split_edge2(
                                    wire, ctx, face, num2, p21, p22, &nv, tol_v, boxes,
                                ) {
                                    nb_split += 1;
                                    step2 = 0;
                                }
                                if split_edge2(
                                    wire, ctx, face, num1, p11, p12, &nv, tol_v, boxes,
                                ) {
                                    nb_split += 1;
                                    num1 -= 1;
                                    restart1 = true;
                                    break 'body;
                                }
                            } else if fix_segment {
                                // `cxx:1536-1817`: the segment is too far to be
                                // a single vertex - split each edge on three
                                // edges and remove the middle (segment) pair.
                                let p01 = GpPnt::new(
                                    0.5 * (pnt11.x() + pnt21.x()),
                                    0.5 * (pnt11.y() + pnt21.y()),
                                    0.5 * (pnt11.z() + pnt21.z()),
                                );
                                let p02 = GpPnt::new(
                                    0.5 * (pnt12.x() + pnt22.x()),
                                    0.5 * (pnt12.y() + pnt22.y()),
                                    0.5 * (pnt12.z() + pnt22.z()),
                                );
                                let mut tol_v1 = pnt11.distance(&p01).max(pnt21.distance(&p01));
                                tol_v1 = tol_v1.max(CONFUSION) * 1.00001;
                                let mut tol_v2 = pnt12.distance(&p02).max(pnt22.distance(&p02));
                                tol_v2 = tol_v2.max(CONFUSION) * 1.00001;
                                if tol_v1 > max_tol_vert || tol_v2 > max_tol_vert {
                                    break 'body;
                                }

                                // `cxx:1573-1608`: analysis of P01 / P02 against
                                // the two endpoints of edge1.
                                let (Some(v1), Some(v2)) =
                                    (first_vertex(&edge1), last_vertex(&edge1))
                                else {
                                    break 'body;
                                };
                                let pv1 = vertex_position(&v1);
                                let pv2 = vertex_position(&v2);
                                let mut new_v1: Option<Vertex> = None;
                                let mut new_v2: Option<Vertex> = None;
                                let mut akey1 = 0i32;
                                let mut akey2 = 0i32;
                                let new_tolerance = tol_v1.max(BRepTool::vertex_tolerance(&v1));
                                if p01.distance(&pv1) < new_tolerance {
                                    let mut nv = TopoBuilder::new().make_vertex(pv1, new_tolerance);
                                    nv.0.set_orientation(v1.0.orientation());
                                    new_v1 = Some(nv);
                                    akey1 += 1;
                                }
                                let new_tolerance = tol_v1.max(BRepTool::vertex_tolerance(&v2));
                                if p01.distance(&pv2) < new_tolerance {
                                    let mut nv = TopoBuilder::new().make_vertex(pv2, new_tolerance);
                                    nv.0.set_orientation(v2.0.orientation());
                                    new_v1 = Some(nv);
                                    akey1 += 1;
                                }
                                let new_tolerance = tol_v2.max(BRepTool::vertex_tolerance(&v1));
                                if p02.distance(&pv1) < new_tolerance {
                                    let mut nv = TopoBuilder::new().make_vertex(pv1, new_tolerance);
                                    nv.0.set_orientation(v1.0.orientation());
                                    new_v2 = Some(nv);
                                    akey2 += 1;
                                }
                                let new_tolerance = tol_v2.max(BRepTool::vertex_tolerance(&v2));
                                if p02.distance(&pv2) < new_tolerance {
                                    let mut nv = TopoBuilder::new().make_vertex(pv2, new_tolerance);
                                    nv.0.set_orientation(v2.0.orientation());
                                    new_v2 = Some(nv);
                                    akey2 += 1;
                                }
                                if akey1 > 1 || akey2 > 1 {
                                    break 'body;
                                }
                                // `cxx:1611-1620`: prepare vertices.
                                if akey1 == 0 {
                                    new_v1 = Some(TopoBuilder::new().make_vertex(p01, tol_v1));
                                }
                                if akey2 == 0 {
                                    new_v2 = Some(TopoBuilder::new().make_vertex(p02, tol_v2));
                                }
                                let (Some(nv1), Some(nv2)) = (new_v1, new_v2) else {
                                    break 'body;
                                };

                                // `cxx:1621-1669`: split edge1.
                                let mut dnum1 = 0i64;
                                let mut numseg1 = num1;
                                if akey1 == 0 && akey2 > 0 {
                                    if split_edge1(
                                        wire, ctx, face, num1, p11, &nv1, tol_v1, boxes,
                                    ) {
                                        nb_split += 1;
                                        dnum1 = 1;
                                        numseg1 = num1 + 1;
                                    }
                                }
                                if akey1 > 0 && akey2 == 0 {
                                    if split_edge1(
                                        wire, ctx, face, num1, p12, &nv2, tol_v2, boxes,
                                    ) {
                                        nb_split += 1;
                                        dnum1 = 1;
                                        numseg1 = num1;
                                    }
                                }
                                if akey1 == 0 && akey2 == 0 {
                                    if split_edge1(
                                        wire, ctx, face, num1, p11, &nv1, tol_v1, boxes,
                                    ) {
                                        nb_split += 1;
                                        dnum1 = 1;
                                    }
                                    if let Some(tmp_e) = wire_edge_at(wire, num1) {
                                        if let Some((_, a, b)) =
                                            curve_on_surface_oriented(&tmp_e, face, false)
                                        {
                                            if (a - p12) * (b - p12) > 0.0 {
                                                // p12 external for [a, b] => split
                                                // the next edge.
                                                if split_edge1(
                                                    wire, ctx, face, num1 + 1, p12, &nv2, tol_v2,
                                                    boxes,
                                                ) {
                                                    nb_split += 1;
                                                    dnum1 += 1;
                                                    numseg1 = num1 + 1;
                                                }
                                            } else if split_edge1(
                                                wire, ctx, face, num1, p12, &nv2, tol_v2, boxes,
                                            ) {
                                                nb_split += 1;
                                                dnum1 += 1;
                                                numseg1 = num1 + 1;
                                            }
                                        }
                                    }
                                }

                                // `cxx:1672-1765`: split edge2 / replace its
                                // vertices. `edge2` follows every rebuild.
                                let mut cur_edge2 = edge2.clone();
                                let mut akey1 = 0i32;
                                let mut akey2 = 0i32;
                                let (Some(mut v12), Some(mut v22)) =
                                    (first_vertex(&cur_edge2), last_vertex(&cur_edge2))
                                else {
                                    break 'body;
                                };
                                let pv12 = vertex_position(&v12);
                                let pv22 = vertex_position(&v22);

                                if p01.distance(&pv12) < tol_v1 {
                                    tol_v1 += p01.distance(&pv12);
                                    update_vertex_tolerance(&nv1, tol_v1);
                                    v12 = replace_vertex_keep_orientation(ctx, &v12, &nv1);
                                    nb_replaced += 1;
                                    let new_e = copy_replace_vertices_with(&cur_edge2, Some(&nv1), Some(&v22));
                                    ctx.replace(&cur_edge2.0, &new_e.0);
                                    wire_set_edge_composed(wire, (num2 + dnum1) as usize, &new_e);
                                    boxes.insert(skey(&new_e), b2_box);
                                    cur_edge2 = new_e;
                                    akey1 = 1;
                                }
                                if p01.distance(&pv22) < tol_v1 {
                                    tol_v1 += p01.distance(&pv22);
                                    update_vertex_tolerance(&nv1, tol_v1);
                                    v22 = replace_vertex_keep_orientation(ctx, &v22, &nv1);
                                    nb_replaced += 1;
                                    let new_e = copy_replace_vertices_with(&cur_edge2, Some(&v12), Some(&nv1));
                                    ctx.replace(&cur_edge2.0, &new_e.0);
                                    wire_set_edge_composed(wire, (num2 + dnum1) as usize, &new_e);
                                    boxes.insert(skey(&new_e), b2_box);
                                    cur_edge2 = new_e;
                                    akey1 = 2;
                                }
                                if p02.distance(&pv12) < tol_v2 {
                                    tol_v2 += p02.distance(&pv12);
                                    update_vertex_tolerance(&nv2, tol_v2);
                                    v12 = replace_vertex_keep_orientation(ctx, &v12, &nv2);
                                    nb_replaced += 1;
                                    let new_e = copy_replace_vertices_with(&cur_edge2, Some(&nv2), Some(&v22));
                                    ctx.replace(&cur_edge2.0, &new_e.0);
                                    wire_set_edge_composed(wire, (num2 + dnum1) as usize, &new_e);
                                    boxes.insert(skey(&new_e), b2_box);
                                    cur_edge2 = new_e;
                                    akey2 = 1;
                                }
                                if p02.distance(&pv22) < tol_v2 {
                                    tol_v2 += p02.distance(&pv22);
                                    update_vertex_tolerance(&nv2, tol_v2);
                                    v22 = replace_vertex_keep_orientation(ctx, &v22, &nv2);
                                    nb_replaced += 1;
                                    let new_e = copy_replace_vertices_with(&cur_edge2, Some(&v12), Some(&nv2));
                                    ctx.replace(&cur_edge2.0, &new_e.0);
                                    wire_set_edge_composed(wire, (num2 + dnum1) as usize, &new_e);
                                    boxes.insert(skey(&new_e), b2_box);
                                    cur_edge2 = new_e;
                                    akey2 = 2;
                                }

                                // `cxx:1766-1814`: split edge2.
                                let mut dnum2 = 0i64;
                                let mut numseg2 = num2 + dnum1;
                                if akey1 == 0 && akey2 > 0 {
                                    if split_edge1(
                                        wire, ctx, face, num2 + dnum1, p21, &nv1, tol_v1, boxes,
                                    ) {
                                        nb_split += 1;
                                        dnum2 = 1;
                                        numseg2 = num2 + dnum1;
                                    }
                                }
                                if akey1 > 0 && akey2 == 0 {
                                    if split_edge1(
                                        wire, ctx, face, num2 + dnum1, p22, &nv2, tol_v2, boxes,
                                    ) {
                                        nb_split += 1;
                                        dnum2 = 1;
                                        numseg2 = num2 + dnum1 + 1;
                                    }
                                }
                                if akey1 == 0 && akey2 == 0 {
                                    if split_edge1(
                                        wire, ctx, face, num2 + dnum1, p21, &nv1, tol_v1, boxes,
                                    ) {
                                        nb_split += 1;
                                        dnum2 = 1;
                                    }
                                    if let Some(tmp_e) = wire_edge_at(wire, num2 + dnum1) {
                                        if let Some((_, a, b)) =
                                            curve_on_surface_oriented(&tmp_e, face, false)
                                        {
                                            if (a - p22) * (b - p22) > 0.0 {
                                                // p22 external for [a, b] => split
                                                // the next edge.
                                                if split_edge1(
                                                    wire,
                                                    ctx,
                                                    face,
                                                    num2 + dnum1 + dnum2,
                                                    p22,
                                                    &nv2,
                                                    tol_v2,
                                                    boxes,
                                                ) {
                                                    nb_split += 1;
                                                    numseg2 = num2 + dnum1 + dnum2;
                                                    dnum2 += 1;
                                                }
                                            } else if split_edge1(
                                                wire,
                                                ctx,
                                                face,
                                                num2 + dnum1,
                                                p22,
                                                &nv2,
                                                tol_v2,
                                                boxes,
                                            ) {
                                                nb_split += 1;
                                                dnum2 += 1;
                                                numseg2 = num2 + dnum1 + 1;
                                            }
                                        }
                                    }
                                }

                                // `cxx:1818-1819`: remove the segment pair.
                                if numseg2 >= 1 {
                                    wire_remove_edge(wire, numseg2 as usize);
                                }
                                if numseg1 >= 1 {
                                    wire_remove_edge(wire, numseg1 as usize);
                                }
                                nb_removed += 2;
                            }
                        }
                    }
                }
            }

            if restart1 {
                // Each `restart1` site already applied the C++ `num1--`, which
                // the outer `num1++` cancels: the same `num1` is re-scanned.
                break;
            }
            num2 += step2;
        }
        num1 += 1;
    }

    // `cxx:1828-1829`.
    let is_done = nb_split != 0 || nb_cut != 0 || nb_replaced != 0 || nb_removed != 0;
    status(nb_split, nb_cut, nb_removed, is_done)
}

/// Endpoint analysis of one edge of a point intersection
/// (`cxx:1129-1163` for edge1, `cxx:1165-1197` for edge2). Updates `tol_v` and
/// `v_sel` in place and returns `modif`.
fn analyze_point_edge(
    edge: &Edge,
    pi: &GpPnt,
    a: f64,
    b: f64,
    param: f64,
    max_tol_vert: f64,
    face: &Face,
    nb_cut: &mut i64,
    tol_v: &mut f64,
    v_sel: &mut Option<Vertex>,
) -> bool {
    let mut modif = false;
    if let (Some(vf), Some(vl)) = (first_vertex(edge), last_vertex(edge)) {
        let pvf = vertex_position(&vf);
        let pvl = vertex_position(&vl);
        let dist1 = pi.distance(&pvf);
        let dist2 = pi.distance(&pvl);
        let distmin = dist1.min(dist2);
        if dist1 != dist2 && distmin < max_tol_vert {
            if dist1 < dist2 {
                *tol_v = (dist1 * 1.00001).max(BRepTool::vertex_tolerance(&vf));
                update_vertex_tolerance(&vf, *tol_v);
                *v_sel = Some(vf);
            } else {
                *tol_v = (dist2 * 1.00001).max(BRepTool::vertex_tolerance(&vl));
                update_vertex_tolerance(&vl, *tol_v);
                *v_sel = Some(vl);
            }
            let dista = (a - param).abs();
            let distb = (b - param).abs();
            let mut is_cut_line = false;
            modif = cut_edge(edge, if dista > distb { a } else { b }, param, face, &mut is_cut_line);
            if modif {
                *nb_cut += 1;
            }
            // not needed split edge, if one of parts is too small
            modif = modif || distmin < CONFUSION;
        }
    }
    modif
}
