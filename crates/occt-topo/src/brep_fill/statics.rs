use super::prelude::*;


pub(super) fn shared_location(edge1: &Edge, edge2: &Edge) -> (bool, TopLocLocation) {
    let l1 = edge1.0.location().clone();
    let l2 = edge2.0.location().clone();
    if l1.is_equal(&l2) {
        (true, l1)
    } else {
        (false, TopLocLocation::identity())
    }
}


/// Trim or copy, then `Transform` when the edges do not share a location,
/// then `Reverse` when the edge is reversed (`BRepFill.cxx:241-290`).
pub(super) fn prepare_curve(edge: &Edge, same_loc: bool) -> Result<(Arc<dyn Curve>, Vertex, Vertex), String> {
    let reg = GeometryRegistry::global();
    let curve = reg
        .edge_curve(&edge.0)
        .ok_or_else(|| "BRep_Tool::Curve".to_string())?;
    let (first, last) = reg.edge_parameters(&edge.0);
    let trimmed = (first - curve.first_parameter()).abs() > PCONFUSION
        || (last - curve.last_parameter()).abs() > PCONFUSION;
    let mut owned = if trimmed {
        Box::new(occt_geom::geom_fill::ParameterTrim::new(curve, first, last)) as Box<dyn Curve>
    } else {
        curve.clone_dyn()
    };
    if !same_loc {
        owned.transform(&edge.0.location().transformation());
    }
    let (stored_first, stored_last) = edge_vertices(edge);
    let (v_first, v_last) = if edge.0.orientation() == Orientation::Reversed {
        owned.reverse();
        (stored_last, stored_first)
    } else {
        (stored_first, stored_last)
    };
    let v_first = v_first.ok_or_else(|| "TopExp::Vertices".to_string())?;
    let v_last = v_last.ok_or_else(|| "TopExp::Vertices".to_string())?;
    Ok((Arc::from(owned), v_first, v_last))
}


pub(super) fn iso_edge(
    surf: &dyn Surface,
    u: f64,
    v0: f64,
    v1: f64,
    v_first: &Vertex,
    v_second: &Vertex,
) -> Result<Edge, String> {
    let iso = surf
        .u_iso_curve(u)
        .ok_or_else(|| "Geom_Surface::UIso".to_string())?;
    let tol = BRepTool::vertex_tolerance(v_first).max(BRepTool::vertex_tolerance(v_second));
    let mut edge = if iso.d0(v0).distance(&iso.d0(v1)) > tol {
        TopoBuilder::new().make_edge(iso, v0, v1)
    } else {
        TopoBuilder::new().make_degenerated_edge(&BRepTool::vertex_point(v_first))
    };
    GeometryRegistry::global().raise_edge_tolerance(&edge.0, CONFUSION);
    TopoBuilder::new().add_edge_vertices(&mut edge, v_first, v_second);
    GeometryRegistry::global().set_ranges_all(&edge.0, v0, v1);
    Ok(edge)
}


pub(super) fn bind_profile_pcurve(
    reg: &GeometryRegistry,
    edge: &Edge,
    face_key: usize,
    _origin_u: f64,
    v: f64,
    u0: f64,
    u1: f64,
) {
    // `BRepFill.cxx:369-389`: a reversed edge stores the negative-X line and
    // the negated U range.
    if edge.0.orientation() == Orientation::Reversed {
        reg.set_edge_pcurve(&edge.0, face_key, uv_line(0.0, v, -1.0, 0.0));
        reg.set_pcurve_range(&edge.0, face_key, -u1, -u0);
    } else {
        reg.set_edge_pcurve(&edge.0, face_key, uv_line(0.0, v, 1.0, 0.0));
        reg.set_pcurve_range(&edge.0, face_key, u0, u1);
    }
}


pub(super) fn uv_line(x: f64, y: f64, dx: f64, dy: f64) -> Arc<dyn Curve2d> {
    let dir = GpDir2d::new(dx, dy).unwrap_or(GpDir2d { x: dx, y: dy });
    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(x, y), dir))
}


pub(super) fn clear_same_flags(edge: &Edge) {
    let reg = GeometryRegistry::global();
    reg.set_same_parameter(&edge.0, false);
    reg.set_same_range(&edge.0, false);
}


/// `BRep_Tool::IsClosed(Edge)` (`BRep_Tool.cxx:1750-1755`).
pub(super) fn edge_is_closed(edge: &Edge) -> bool {
    match edge_vertices(edge) {
        (Some(first), Some(last)) => is_same(&first.0, &last.0),
        _ => false,
    }
}


/// `BRep_Tool::IsClosed(Shell)` (`BRep_Tool.cxx:1709-1728`).
pub(super) fn shell_is_closed(shell: &TopoShape) -> bool {
    let mut odd: Vec<TopoShape> = Vec::new();
    let mut has_bound = false;
    for face in faces_of(shell) {
        for wire in wires_of_face(&face) {
            for edge in edges_of_wire(&wire) {
                if BRepTool::is_degenerated(&edge)
                    || matches!(
                        edge.0.orientation(),
                        Orientation::Internal | Orientation::External
                    )
                {
                    continue;
                }
                has_bound = true;
                toggle_same(&mut odd, &edge.0);
            }
        }
    }
    has_bound && odd.is_empty()
}


/// `TopTools_Map` toggle of `BRepLib_Fill`'s edge pairing: `IsSame` already
/// covers `TShape` and location (`TopoDS_Shape.hxx:262-267`).
pub(super) fn toggle_same(odd: &mut Vec<TopoShape>, shape: &TopoShape) {
    if let Some(i) = odd.iter().position(|kept| is_same(kept, shape)) {
        odd.remove(i);
    } else {
        odd.push(shape.clone());
    }
}
