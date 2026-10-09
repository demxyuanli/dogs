use super::*;

/// `ShapeFix_IntersectionTool::SplitEdge` (`cxx:85-189`).
///
/// Cuts `edge` at the pcurve parameter `param` using the existing vertex
/// `vert`, returning the two pieces `[first, param]` / `[param, last]` in
/// wire order. `None` is the OCCT `return false`.
pub(super) fn split_edge(
    edge: &Edge,
    param: f64,
    vert: &Vertex,
    face: &Face,
    preci: f64,
) -> Option<(Edge, Edge)> {
    // `cxx:99-103`.
    let v1 = first_vertex(edge)?;
    let v2 = last_vertex(edge)?;
    if is_same(&v1.0, &vert.0) || is_same(&v2.0, &vert.0) {
        return None;
    }

    // `cxx:105-110`: `sae.PCurve(edge, face, c2d, a, b, true)`.
    let (c2d, a, b) = curve_on_surface_oriented(edge, face, true)?;
    if (a - param).abs() < 0.01 * preci || (b - param).abs() < 0.01 * preci {
        return None;
    }

    // `cxx:112-137`: check the distance between the edge and the new vertex.
    let p1 = if BRepTool::same_parameter(edge) && !BRepTool::is_degenerated(edge) {
        let c3d = BRepTool::edge_curve_world(edge)?;
        c3d.d0(param)
    } else {
        // `cxx:130-136`: `BRep_Tool::Surface(face, L)` then `P1.Transformed(L)`.
        let surf = BRepTool::face_surface_world(face)?;
        let p2d = c2d.d0(param);
        surf.d0(p2d.x(), p2d.y())
    };
    let p2 = vertex_position(vert);
    if p1.distance(&p2) > preci {
        // `cxx:139-143`: the cxx keeps going after widening the tolerance.
        update_vertex_tolerance(vert, p1.distance(&p2));
    }

    // `cxx:145-150`: `ShapeAnalysis_TransferParametersProj`, `SetMaxTolerance`
    // then `Init`.
    let mut transfer = TransferParametersProj::new();
    transfer.set_max_tolerance(preci);
    transfer.init(edge, face);

    // `cxx:151-160`.
    let (first, last) = if a < b { (a, b) } else { (b, a) };

    // `cxx:162-176`.
    let orient = edge.0.orientation();
    let mut we = edge.clone();
    we.0.set_orientation(Orientation::Forward);

    let vert_rev = Vertex(vert.0.oriented(Orientation::Reversed));
    let fv = first_vertex(&we)?;
    let mut new_e1 = copy_replace_vertices_with(&we, Some(&fv), Some(&vert_rev));
    copy_pcurves(&new_e1, &we);
    transfer.transfer_range(&mut new_e1, first, param, true);
    let reg = GeometryRegistry::global();
    reg.set_same_range(&new_e1.0, false);
    reg.set_same_parameter(&new_e1.0, false);

    let vert_fwd = Vertex(vert.0.oriented(Orientation::Forward));
    let lv = last_vertex(&we)?;
    let mut new_e2 = copy_replace_vertices_with(&we, Some(&vert_fwd), Some(&lv));
    copy_pcurves(&new_e2, &we);
    transfer.transfer_range(&mut new_e2, param, last, true);
    reg.set_same_range(&new_e2.0, false);
    reg.set_same_parameter(&new_e2.0, false);

    // `cxx:180-188`.
    new_e1.0.set_orientation(orient);
    new_e2.0.set_orientation(orient);
    if orient == Orientation::Reversed {
        std::mem::swap(&mut new_e1, &mut new_e2);
    }

    Some((new_e1, new_e2))
}

/// `ShapeFix_IntersectionTool::CutEdge` (`cxx:194-268`).
///
/// Distinct from `ShapeFix_SplitTool::CutEdge` (`split_tool::cut_edge`): the
/// same-parameter tail here has **no** `ShapeAnalysis_Curve::ValidateRange` and
/// no `ShapeFix_Edge::FixSameParameter` - it only writes the range
/// (`cxx:258-259`).
pub(super) fn cut_edge(edge: &Edge, pend: f64, cut: f64, face: &Face, is_cut_line: &mut bool) -> bool {
    // `cxx:210-212`.
    if (cut - pend).abs() < 10.0 * PCONFUSION {
        return false;
    }
    let a_range = (cut - pend).abs();
    let reg = GeometryRegistry::global();
    // `cxx:213-214`: `BRep_Tool::Range(edge, a, b)`.
    let (a, b) = reg.edge_parameters(&edge.0);
    *is_cut_line = false; // `cxx:219`
    if a_range < 10.0 * PCONFUSION {
        return false;
    }

    // `cxx:225-246`: pcurve trimmed-of-line arm.
    if !reg.same_parameter(&edge.0) {
        if let Some((crv, fp, lp)) = curve_on_surface_oriented(edge, face, false) {
            // `cxx:221`: only a `Geom2d_TrimmedCurve` pcurve enters the line arm;
            // any other pcurve falls through to `return false` (`cxx:245`).
            let Some(basis) = crv.trimmed_basis() else {
                return false;
            };
            if basis.is_line() {
                // `cxx:239`: all representations get the 2D cut interval.
                reg.set_ranges_all(&edge.0, pend.min(cut), pend.max(cut));
                if (pend - lp).abs() < PCONFUSION {
                    // `cxx:229-233`: cut from the beginning.
                    let cut3d = (cut - fp) * (b - a) / (lp - fp);
                    reg.set_edge_range(&edge.0, a + cut3d, b);
                    *is_cut_line = true;
                } else if (pend - fp).abs() < PCONFUSION {
                    // `cxx:234-238`: cut from the end.
                    let cut3d = (lp - cut) * (b - a) / (lp - fp);
                    reg.set_edge_range(&edge.0, a, b - cut3d);
                    *is_cut_line = true;
                }
            }
            return true; // `cxx:242-245`
        }
        return false;
    }

    // `cxx:249-255`.
    if ((a - b).abs() - a_range).abs() < PCONFUSION {
        return false;
    }
    if a_range < 10.0 * PCONFUSION {
        return false;
    }
    // `cxx:258-259`: `B.Range(edge, min, max)`, every representation.
    reg.set_ranges_all(&edge.0, pend.min(cut), pend.max(cut));
    true
}

/// `ShapeFix_IntersectionTool::SplitEdge1` (`cxx:270-356`).
pub(super) fn split_edge1(
    wire: &mut Wire,
    ctx: &mut dyn ReShape,
    face: &Face,
    num: i64,
    param: f64,
    vert: &Vertex,
    preci: f64,
    boxes: &mut Boxes,
) -> bool {
    // `cxx:279`.
    if num < 1 || num > wire_len(wire) {
        return false;
    }
    let Some(edge) = wire_edge_at(wire, num) else {
        return false;
    };
    let Some((new_e1, new_e2)) = split_edge(&edge, param, vert, face, preci) else {
        return false;
    };

    // `cxx:288-301`: change context.
    let mut wd = Wire::new();
    for e in [&new_e1, &new_e2] {
        wire_add_edge_composed(&mut wd, 0, e);
    }
    wd.0.set_orientation(Orientation::Forward);
    ctx.replace(&edge.0, &wd.0);
    // `BRepTools::Update(E)` (`cxx:298`) only invalidates the edge's cached
    // triangulation - the port carries none, so it is a no-op.

    // `cxx:304-311`: change sewd.
    wire_set_edge_composed(wire, num as usize, &new_e1);
    if num == wire_len(wire) {
        wire_add_edge_composed(wire, 0, &new_e2);
    } else {
        wire_add_edge_composed(wire, (num + 1) as usize, &new_e2);
    }

    // `cxx:314-354`: change boxes.
    boxes.remove(&skey(&edge));
    for e in [&new_e1, &new_e2] {
        if let Some((c2d, cf, cl)) = curve_on_surface_oriented(e, face, false) {
            boxes.insert(skey(e), curve2d_box(c2d.as_ref(), cf, cl));
        }
    }

    true
}

/// `ShapeFix_IntersectionTool::SplitEdge2` (`cxx:363-490`).
///
/// Splits at the middle of `(param1, param2)` then cuts the removed segment
/// away from both halves.
#[allow(clippy::too_many_arguments)]
pub(super) fn split_edge2(
    wire: &mut Wire,
    ctx: &mut dyn ReShape,
    face: &Face,
    num: i64,
    param1: f64,
    param2: f64,
    vert: &Vertex,
    preci: f64,
    boxes: &mut Boxes,
) -> bool {
    if num < 1 || num > wire_len(wire) {
        return false;
    }
    let Some(edge) = wire_edge_at(wire, num) else {
        return false;
    };
    let param = 0.5 * (param1 + param2);
    let Some((new_e1, new_e2)) = split_edge(&edge, param, vert, face, preci) else {
        return false;
    };

    // `cxx:383-420`: cut new edges by param1 and param2.
    let mut is_cut_line = false;
    if let Some((crv1, fp1, lp1)) = curve_on_surface_oriented(&new_e1, face, false) {
        let _ = &crv1;
        if let Some((_crv2, fp2, lp2)) = curve_on_surface_oriented(&new_e2, face, false) {
            if lp1 == param {
                if (lp1 - fp1) * (lp1 - param1) > 0.0 {
                    cut_edge(&new_e1, fp1, param1, face, &mut is_cut_line);
                    cut_edge(&new_e2, lp2, param2, face, &mut is_cut_line);
                } else {
                    cut_edge(&new_e1, fp1, param2, face, &mut is_cut_line);
                    cut_edge(&new_e2, lp2, param1, face, &mut is_cut_line);
                }
            } else if (fp1 - lp1) * (fp1 - param1) > 0.0 {
                cut_edge(&new_e1, lp1, param1, face, &mut is_cut_line);
                cut_edge(&new_e2, fp2, param2, face, &mut is_cut_line);
            } else {
                cut_edge(&new_e1, lp1, param2, face, &mut is_cut_line);
                cut_edge(&new_e2, fp2, param1, face, &mut is_cut_line);
            }
        }
    }

    // `cxx:422-435`: change context.
    let mut wd = Wire::new();
    for e in [&new_e1, &new_e2] {
        wire_add_edge_composed(&mut wd, 0, e);
    }
    wd.0.set_orientation(Orientation::Forward);
    ctx.replace(&edge.0, &wd.0);

    // `cxx:437-444`: change sewd.
    wire_set_edge_composed(wire, num as usize, &new_e1);
    if num == wire_len(wire) {
        wire_add_edge_composed(wire, 0, &new_e2);
    } else {
        wire_add_edge_composed(wire, (num + 1) as usize, &new_e2);
    }

    // `cxx:447-487`: change boxes.
    boxes.remove(&skey(&edge));
    for e in [&new_e1, &new_e2] {
        if let Some((c2d, cf, cl)) = curve_on_surface_oriented(e, face, false) {
            boxes.insert(skey(e), curve2d_box(c2d.as_ref(), cf, cl));
        }
    }

    true
}

/// The four `if (V21*.IsSame(V2*))` blocks every `UnionVertexes` arm repeats
/// (`cxx:539-608`, `:629-698`, `:718-787`, `:808-877`).
///
/// `keep` is the vertex of `edge1` the pair is welded onto, `drop` the
/// vertex of `edge2` being replaced in the two neighbouring edges of the wire
/// (`num2 - 1` and `num2 + 1`, wrapping).
pub(super) fn replace_in_neighbours(
    wire: &mut Wire,
    ctx: &mut dyn ReShape,
    boxes: &mut Boxes,
    num2: i64,
    keep: &Vertex,
    drop: &Vertex,
) {
    let nb = wire_len(wire);
    let num21 = if num2 > 1 { num2 - 1 } else { nb };
    let num22 = if num2 < nb { num2 + 1 } else { 1 };
    let (Some(edge21), Some(edge22)) = (wire_edge_at(wire, num21), wire_edge_at(wire, num22)) else {
        return;
    };
    let v21f = first_vertex(&edge21);
    let v21l = last_vertex(&edge21);
    let v22f = first_vertex(&edge22);
    let v22l = last_vertex(&edge22);

    let subst = |wire: &mut Wire,
                     ctx: &mut dyn ReShape,
                     boxes: &mut Boxes,
                     old: &Edge,
                     new_e: &Edge,
                     at: i64| {
        if let Some(b) = boxes.get(&skey(old)).copied() {
            boxes.insert(skey(new_e), b);
        }
        ctx.replace(&old.0, &new_e.0);
        wire_set_edge_composed(wire, at as usize, new_e);
    };

    // `cxx:564-585`.
    if v21f.as_ref().is_some_and(|v| is_same(&v.0, &drop.0)) {
        let new_e = copy_replace_vertices_with(&edge21, Some(keep), v21l.as_ref());
        subst(wire, ctx, boxes, &edge21, &new_e, num21);
    }
    if v21l.as_ref().is_some_and(|v| is_same(&v.0, &drop.0)) {
        let new_e = copy_replace_vertices_with(&edge21, v21f.as_ref(), Some(keep));
        subst(wire, ctx, boxes, &edge21, &new_e, num21);
    }
    // `cxx:586-605`.
    if v22f.as_ref().is_some_and(|v| is_same(&v.0, &drop.0)) {
        let new_e = copy_replace_vertices_with(&edge22, Some(keep), v22l.as_ref());
        subst(wire, ctx, boxes, &edge22, &new_e, num22);
    }
    if v22l.as_ref().is_some_and(|v| is_same(&v.0, &drop.0)) {
        let new_e = copy_replace_vertices_with(&edge22, v22f.as_ref(), Some(keep));
        subst(wire, ctx, boxes, &edge22, &new_e, num22);
    }
}

