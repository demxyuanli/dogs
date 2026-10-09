use super::*;

/// `ShapeFix_IntersectionTool::UnionVertexes` (`cxx:495-880`).
///
/// Picks the closest of the four endpoint pairs, welds `edge1`'s vertex onto
/// `edge2` (a fresh `CopyReplaceVertices`) and then replaces the dropped vertex
/// in the two neighbouring edges. `false` when no pair is within tolerance.
pub(super) fn union_vertexes(
    wire: &mut Wire,
    ctx: &mut dyn ReShape,
    edge1: &Edge,
    edge2: &Edge,
    num2: i64,
    boxes: &mut Boxes,
    b2: &BndBox2d,
) -> bool {
    // `cxx:509-520`.
    let (Some(v1f), Some(v1l), Some(v2f), Some(v2l)) = (
        first_vertex(edge1),
        last_vertex(edge1),
        first_vertex(edge2),
        last_vertex(edge2),
    ) else {
        return false;
    };
    let pv1f = vertex_position(&v1f);
    let pv1l = vertex_position(&v1l);
    let pv2f = vertex_position(&v2f);
    let pv2l = vertex_position(&v2l);
    let d11 = pv1f.distance(&pv2f);
    let d12 = pv1f.distance(&pv2l);
    let d21 = pv1l.distance(&pv2f);
    let d22 = pv1l.distance(&pv2l);

    if d11 < d12 && d11 < d21 && d11 < d22 {
        // `cxx:521-528`: union vertexes V1F and V2F.
        let tolv = BRepTool::vertex_tolerance(&v1f).max(BRepTool::vertex_tolerance(&v2f));
        if is_same(&v2f.0, &v1f.0) || d11 >= tolv {
            return false;
        }
        update_vertex_tolerance(&v1f, tolv);
        let new_e = copy_replace_vertices_with(edge2, Some(&v1f), Some(&v2l));
        ctx.replace(&edge2.0, &new_e.0);
        wire_set_edge_composed(wire, num2 as usize, &new_e);
        boxes.insert(skey(&new_e), *b2);
        replace_in_neighbours(wire, ctx, boxes, num2, &v1f, &v2f);
        true
    } else if d12 < d21 && d12 < d22 {
        // `cxx:610-617`: union vertexes V1F and V2L.
        let tolv = BRepTool::vertex_tolerance(&v1f).max(BRepTool::vertex_tolerance(&v2l));
        if is_same(&v2l.0, &v1f.0) || d12 >= tolv {
            return false;
        }
        update_vertex_tolerance(&v1f, tolv);
        let new_e = copy_replace_vertices_with(edge2, Some(&v2f), Some(&v1f));
        ctx.replace(&edge2.0, &new_e.0);
        wire_set_edge_composed(wire, num2 as usize, &new_e);
        boxes.insert(skey(&new_e), *b2);
        replace_in_neighbours(wire, ctx, boxes, num2, &v1f, &v2l);
        true
    } else if d21 < d22 {
        // `cxx:700-707`: union vertexes V1L and V2F.
        let tolv = BRepTool::vertex_tolerance(&v1l).max(BRepTool::vertex_tolerance(&v2f));
        if is_same(&v2f.0, &v1l.0) || d21 >= tolv {
            return false;
        }
        update_vertex_tolerance(&v1l, tolv);
        let new_e = copy_replace_vertices_with(edge2, Some(&v1l), Some(&v2l));
        ctx.replace(&edge2.0, &new_e.0);
        wire_set_edge_composed(wire, num2 as usize, &new_e);
        boxes.insert(skey(&new_e), *b2);
        replace_in_neighbours(wire, ctx, boxes, num2, &v1l, &v2f);
        true
    } else {
        // `cxx:789-796`: union vertexes V1L and V2L.
        let tolv = BRepTool::vertex_tolerance(&v1l).max(BRepTool::vertex_tolerance(&v2l));
        if is_same(&v2l.0, &v1l.0) || d22 >= tolv {
            return false;
        }
        update_vertex_tolerance(&v1l, tolv);
        let new_e = copy_replace_vertices_with(edge2, Some(&v2f), Some(&v1l));
        ctx.replace(&edge2.0, &new_e.0);
        wire_set_edge_composed(wire, num2 as usize, &new_e);
        boxes.insert(skey(&new_e), *b2);
        replace_in_neighbours(wire, ctx, boxes, num2, &v1l, &v2l);
        true
    }
}

/// `CreateBoxes2d` (`cxx:884-920`): the 2D pcurve boxes of every wire edge plus
/// their total box.
pub(super) fn create_boxes_2d(wire: &Wire, face: &Face, boxes: &mut Boxes) -> BndBox2d {
    let mut total = BndBox2d::new();
    for e in edges_of_wire(wire) {
        let Some((c2d, cf, cl)) = curve_on_surface_oriented(&e, face, false) else {
            continue;
        };
        let box2d = curve2d_box(c2d.as_ref(), cf, cl);
        boxes.insert(skey(&e), box2d);
        total.add_box(&box2d);
    }
    total
}

/// `SelectIntPnt` (`cxx:923-962`): `Inter.Point(1)`, replaced by `Point(2)`
/// when the second point has more `Middle` positions.
pub(super) fn select_int_pnt(inter: &Geom2dIntGInter) -> IntRes2dIntersectionPoint {
    let mut ip = inter.point(1);
    if inter.nb_points() == 2 {
        let status1 = int_pnt_status(&ip);
        let ip2 = inter.point(2);
        let status2 = int_pnt_status(&ip2);
        if status2 > status1 {
            ip = ip2;
        }
    }
    ip
}

pub(super) fn int_pnt_status(ip: &IntRes2dIntersectionPoint) -> i32 {
    let mut status = 0;
    if ip.transition_of_first().position_on_curve() == IntRes2dPosition::Middle {
        status += 1;
    }
    if ip.transition_of_second().position_on_curve() == IntRes2dPosition::Middle {
        status += 2;
    }
    status
}

/// `ShapeFix_IntersectionTool::FindVertAndSplitEdge` (`cxx:967-1025`).
///
/// Takes the closer endpoint of `edge2` to the point of `edge1` at `param1` and
/// splits `edge1` there. `num1` is decremented on success so the caller's
/// `for (num1++)` re-visits the same edge, matching `cxx:1019`.
#[allow(clippy::too_many_arguments)]
pub(super) fn find_vert_and_split_edge(
    param1: f64,
    edge1: &Edge,
    edge2: &Edge,
    crv1: &dyn Curve2d,
    max_tol_vert: &mut f64,
    num1: &mut i64,
    wire: &mut Wire,
    face: &Face,
    boxes: &mut Boxes,
    a_tmp_key: bool,
    ctx: &mut dyn ReShape,
) -> bool {
    // `cxx:981-983`.
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let pi1 = point_on_edge(edge1, surf.as_ref(), crv1, param1);

    // `cxx:988-991`.
    let (Some(v1), Some(v2)) = (first_vertex(edge2), last_vertex(edge2)) else {
        return false;
    };
    let pv1 = vertex_position(&v1);
    let pv2 = vertex_position(&v2);
    let v11 = first_vertex(edge1);
    let v12 = last_vertex(edge1);

    // `cxx:992-1011`.
    let mut need_split = true;
    let tol1 = ((pi1.distance(&pv1) / 2.0) * 1.00001).max(BRepTool::vertex_tolerance(&v1));
    let tol2 = ((pi1.distance(&pv2) / 2.0) * 1.00001).max(BRepTool::vertex_tolerance(&v2));
    let (v, tol_v) = if pi1.distance(&pv1) < pi1.distance(&pv2) {
        if v11.as_ref().is_some_and(|x| is_same(&x.0, &v1.0))
            || v12.as_ref().is_some_and(|x| is_same(&x.0, &v1.0))
        {
            need_split = false;
        }
        (v1, tol1)
    } else {
        if v11.as_ref().is_some_and(|x| is_same(&x.0, &v2.0))
            || v12.as_ref().is_some_and(|x| is_same(&x.0, &v2.0))
        {
            need_split = false;
        }
        (v2, tol2)
    };

    // `cxx:1013-1024`. OCCT passes `tolV` in the `preci` slot of `SplitEdge1`.
    if (need_split || a_tmp_key) && split_edge1(wire, ctx, face, *num1, param1, &v, tol_v, boxes) {
        update_vertex_tolerance(&v, tol_v);
        *max_tol_vert = max_tol_vert.max(tol_v);
        *num1 -= 1;
        return true;
    }
    false
}

