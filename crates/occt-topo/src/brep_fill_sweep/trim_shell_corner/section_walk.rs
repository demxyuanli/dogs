//! The section walk: `CheckAndOrientEdges`, `FindVertex`, `FindNextVertex`,
//! `GetPave`, `FindNextEdge` and `FillGap`.
//!
//! Moved here verbatim from the single-file port; the `BRepFill_TrimShellCorner.cxx`
//! branch structure is unchanged.

use super::*;
use super::helpers::*;

/// `CheckAndOrientEdges` (`cxx:1705-1827`): orient the section edges into one
/// chain and report whether their pcurve ends reach both gap endpoints.
///
/// `the_ordered_list` is the chain from `FindNextEdge`; `true` means both the
/// first and the last 2D point of the chain land on the gap's endpoints within
/// the surface resolution of the adjacent vertex tolerance.
pub(super) fn check_and_orient_edges(
    the_ordered_list: &[TopoShape],
    the_first_point: &GpPnt2d,
    the_last_point: &GpPnt2d,
    the_face: &TopoShape,
    the_oriented_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:1712-1715`.
    if the_ordered_list.is_empty() {
        return true;
    }
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(the_face);

    // `cxx:1718-1721`.
    let mut a_e_prev = Edge(the_ordered_list[0].clone());

    // `cxx:1723-1735`.
    let Some(a_curve) = reg.edge_pcurve(&a_e_prev.0, key) else {
        return false;
    };
    let (pf, pl) = reg.pcurve_range(&a_e_prev.0, key).unwrap_or((0.0, 0.0));
    let (vf, vl) = edge_vertices(&a_e_prev);
    let a_tolerance1 = match &vf {
        Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
        None => surface_tolerance(the_face, CONFUSION),
    };
    let a_tolerance2 = match &vl {
        Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
        None => surface_tolerance(the_face, CONFUSION),
    };

    // `cxx:1737-1777`.
    let mut ap = a_curve.d0(pf);
    let mut first_found = false;
    let mut last_found = false;
    if pnt2d_distance(&ap, the_first_point) < a_tolerance1 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Forward));
        }
        first_found = true;
    } else if pnt2d_distance(&ap, the_last_point) < a_tolerance1 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Reversed));
        }
        last_found = true;
    }
    ap = a_curve.d0(pl);
    if pnt2d_distance(&ap, the_last_point) < a_tolerance2 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Forward));
        }
        last_found = true;
    } else if pnt2d_distance(&ap, the_first_point) < a_tolerance2 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Reversed));
        }
        first_found = true;
    }

    // `cxx:1779-1782`.
    if !the_oriented_list.is_empty() {
        a_e_prev = Edge(the_oriented_list[the_oriented_list.len() - 1].clone());
    }

    // `cxx:1784-1823`.
    for a_shape in the_ordered_list.iter().skip(1) {
        let a_e = Edge(a_shape.clone());
        let (a_v11, a_v12) = edge_vertices_cum_ori(&a_e_prev);
        let (a_v21, a_v22) = edge_vertices(&a_e);

        // `cxx:1788-1790`.
        let closed = match (&a_v12, &a_v21, &a_v11, &a_v22) {
            (Some(v12), Some(v21), _, _) if is_same(&v12.0, &v21.0) => true,
            (_, _, Some(v11), Some(v22)) if is_same(&v11.0, &v22.0) => true,
            _ => false,
        };
        let an_ori = if closed {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        the_oriented_list.push(oriented_to(&a_e.0, an_ori));
        a_e_prev = Edge(the_oriented_list[the_oriented_list.len() - 1].clone());

        // `cxx:1794-1801`.
        let a_tolerance1 = match &a_v21 {
            Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
            None => surface_tolerance(the_face, CONFUSION),
        };
        let a_tolerance2 = match &a_v22 {
            Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
            None => surface_tolerance(the_face, CONFUSION),
        };
        // `cxx:1802-1822`.
        let Some(a_curve) = reg.edge_pcurve(&a_e.0, key) else {
            return false;
        };
        let (pf, pl) = reg.pcurve_range(&a_e.0, key).unwrap_or((0.0, 0.0));
        let ap = a_curve.d0(pf);
        if pnt2d_distance(&ap, the_first_point) < a_tolerance1 {
            first_found = true;
        } else if pnt2d_distance(&ap, the_last_point) < a_tolerance1 {
            last_found = true;
        }
        let ap = a_curve.d0(pl);
        if pnt2d_distance(&ap, the_first_point) < a_tolerance2 {
            first_found = true;
        } else if pnt2d_distance(&ap, the_last_point) < a_tolerance2 {
            last_found = true;
        }
    }

    first_found && last_found
}

/// `FindVertex` (`cxx:1829-1887`): the vertex and pave of a `u`-iso edge at
/// rank `the_rank`, taken from its first `myHistMap` image.
pub(super) fn find_vertex(
    the_edge: &TopoShape,
    the_rank: usize,
    the_ds: &BopdsDS,
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_vertex: &mut Option<TopoShape>,
    the_pave: &mut BopdsPave,
) -> bool {
    // `cxx:1839-1848`.
    let Some(lst) = the_hist_map.get(&shape_map_key(the_edge)) else {
        return false;
    };
    let Some(first) = lst.first() else {
        return false;
    };

    // `cxx:1851-1866`.
    let a_new_edge = Edge(first.clone());
    let (f, l) = GeometryRegistry::global().edge_parameters(&a_new_edge.0);
    let v = if the_rank == 0 {
        the_pave.set_parameter(l);
        last_vertex(&a_new_edge)
    } else {
        the_pave.set_parameter(f);
        first_vertex(&a_new_edge)
    };
    *the_vertex = v.map(|v| v.0);
    let vertex = match the_vertex.as_ref() {
        Some(v) => v,
        None => return true, // `cxx:2001-2004` guards the null vertex.
    };

    // `cxx:1867-1881`: the index in the DS, searching the new shapes when the
    // vertex is not one of the source shapes.
    let an_index = match the_ds.index(vertex) {
        Some(i) => i,
        None => {
            let mut found = usize::MAX;
            for i in the_ds.nb_source_shapes()..the_ds.nb_shapes() {
                if let Some(s) = the_ds.shape(i) {
                    if is_same(s, vertex) {
                        found = i;
                        break;
                    }
                }
            }
            found
        }
    };
    if an_index == usize::MAX {
        // The reference stores `-1`; no caller can then read a shape back from
        // `theDS->Shape(-1)`, so the walk reports failure where OCCT reaches a
        // null vertex. UNPORTED: `BOPDS_DS::Shape(-1)` has no Rust equivalent.
        the_pave.set_index(usize::MAX);
    } else {
        the_pave.set_index(an_index);
    }
    true
}

/// `FindNextVertex` (`cxx:1889-1931`): the first pave after `the_prev_pave` on
/// the bound edge whose shape is a *new* shape of the DS.
pub(super) fn find_next_vertex(
    the_edge_index: usize,
    the_prev_pave: &BopdsPave,
    the_ds: &BopdsDS,
    the_next_vertex: &mut Option<TopoShape>,
    the_pave: &mut BopdsPave,
) -> bool {
    // `cxx:1900-1902`: `BOPDS_Pave anullpave;` is the default-constructed
    // (index 0, parameter 0) pave.
    let null_pave = BopdsPave::new(0, 0.0);
    let mut b_take_pave = the_prev_pave.is_equal(&null_pave);
    let mut b_found = false;

    // `cxx:1904-1922`.
    for a_tmp_pave in the_ds.paves(the_edge_index) {
        if b_take_pave {
            if the_ds.is_new_shape(a_tmp_pave.index()) {
                *the_next_vertex = the_ds.shape(a_tmp_pave.index()).cloned();
                *the_pave = a_tmp_pave;
                b_found = true;
                break;
            }
        } else if a_tmp_pave.is_equal(the_prev_pave) {
            b_take_pave = true;
        }
    }
    b_found
}

/// `GetPave` (`cxx:1933-1958`): the first or the last pave of the edge.
pub(super) fn get_pave(
    the_edge_index: usize,
    is_first: bool,
    the_ds: &BopdsDS,
    the_pave: &mut BopdsPave,
) -> bool {
    let a_lp = the_ds.paves(the_edge_index);
    if a_lp.is_empty() {
        return false;
    }
    *the_pave = if is_first {
        a_lp[0]
    } else {
        a_lp[a_lp.len() - 1]
    };
    true
}

/// `FindNextEdge` (`cxx:1616-1703`): the edge chain that walks from
/// `the_first_vertex` to `the_last_vertex` inside `the_map_ve`.
///
/// The recursion carries its own copy of `the_map_to_avoid` (`cxx:1626-1627`),
/// so a dead end unwinds with the caller's list cleared and its cursor
/// restored (`cxx:1669-1672`).
#[allow(clippy::too_many_arguments)]
pub(super) fn find_next_edge(
    the_first_vertex: &TopoShape,
    the_last_vertex: &TopoShape,
    the_map_ve: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_map_to_avoid: &[ShapeMapKey],
    the_ordered_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:1626-1629`.
    let mut a_cur_vertex = the_first_vertex.clone();
    let mut a_map_to_avoid: Vec<ShapeMapKey> = the_map_to_avoid.to_vec();
    let mut a_list_of_edge: Vec<TopoShape> = Vec::new();
    let extent = the_map_ve.len();

    // `cxx:1630-1698`.
    for _ in 0..extent {
        let Some(lste) = the_map_ve.get(&shape_map_key(&a_cur_vertex)) else {
            break; // `cxx:1632-1635`.
        };
        let mut b_efound = false;
        for an_edge in lste.clone() {
            let a_save_cur_vertex = a_cur_vertex.clone();
            if a_map_to_avoid.contains(&shape_map_key(&an_edge)) {
                continue;
            }
            // `cxx:1647-1659`: move the cursor to the far end.
            let (v1, v2) = edge_vertices(&Edge(an_edge.clone()));
            if let Some(v1) = &v1 {
                if !is_same(&a_cur_vertex, &v1.0) {
                    a_cur_vertex = v1.0.clone();
                } else if let Some(v2) = &v2 {
                    if !is_same(&a_cur_vertex, &v2.0) {
                        a_cur_vertex = v2.0.clone();
                    }
                }
            } else if let Some(v2) = &v2 {
                if !is_same(&a_cur_vertex, &v2.0) {
                    a_cur_vertex = v2.0.clone();
                }
            }
            a_map_to_avoid.push(shape_map_key(&an_edge));
            b_efound = true;
            a_list_of_edge.push(an_edge);

            // `cxx:1663-1680`.
            if !is_same(&a_cur_vertex, the_last_vertex) {
                let mut a_list_tmp: Vec<TopoShape> = Vec::new();
                if !find_next_edge(
                    &a_cur_vertex,
                    the_last_vertex,
                    the_map_ve,
                    &a_map_to_avoid,
                    &mut a_list_tmp,
                ) {
                    a_list_of_edge.clear();
                    a_cur_vertex = a_save_cur_vertex;
                    continue;
                }
                a_list_of_edge.extend(a_list_tmp);
                the_ordered_list.extend(a_list_of_edge);
                return true;
            }
            break;
        }

        // `cxx:1683-1686`.
        if is_same(&a_cur_vertex, the_last_vertex) {
            break;
        }
        // `cxx:1688-1691`.
        if !b_efound {
            return false;
        }
    }

    // `cxx:1694-1697`.
    if is_same(&a_cur_vertex, the_last_vertex) {
        the_ordered_list.extend(a_list_of_edge);
        return true;
    }
    false
}

/// `FillGap` (`cxx:1577-1613`): the oriented section edges that bridge the 2D
/// points `the_first_point` / `the_last_point` between the two vertices.
pub(super) fn fill_gap(
    the_first_vertex: &TopoShape,
    the_last_vertex: &TopoShape,
    the_first_point: &GpPnt2d,
    the_last_point: &GpPnt2d,
    the_face: &TopoShape,
    the_section_edges: &[TopoShape],
    the_ordered_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:1587-1593`: the `(VERTEX, EDGE)` ancestors of the section compound.
    let comp = TopoBuilder::new().make_compound_of(the_section_edges);
    let a_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
    if a_map.is_empty() {
        return false;
    }
    // `cxx:1595-1598`.
    if !a_map.contains_key(&shape_map_key(the_first_vertex))
        || !a_map.contains_key(&shape_map_key(the_last_vertex))
    {
        return false;
    }

    // `cxx:1599-1612`.
    let mut a_list_of_edge: Vec<TopoShape> = Vec::new();
    let a_map_to_avoid: Vec<ShapeMapKey> = Vec::new();
    if find_next_edge(
        the_first_vertex,
        the_last_vertex,
        &a_map,
        &a_map_to_avoid,
        &mut a_list_of_edge,
    ) && !a_list_of_edge.is_empty()
    {
        return check_and_orient_edges(
            &a_list_of_edge,
            the_first_point,
            the_last_point,
            the_face,
            the_ordered_list,
        );
    }
    false
}

