//! `ChooseSection` and its geometry leaves.
//!
//! Moved here verbatim from the single-file port; the `BRepFill_TrimShellCorner.cxx`
//! branch structure is unchanged.

use super::*;
use super::helpers::*;

/// `ComputeAveragePlaneAndMaxDeviation` (`cxx:2620-2658`).
///
/// `IsSingular` reports a degenerate point set (all points on a line), in
/// which case the reference returns `-1` (`cxx:2645-2648`) and leaves
/// `thePlane` at whatever `GeomLib::AxeOfInertia` produced.
pub(super) fn compute_average_plane_and_max_deviation(
    the_wire: &TopoShape,
    the_plane: &mut GpPln,
    is_singular: &mut bool,
) -> f64 {
    // `cxx:2623-2624`.
    const N: usize = 40;
    *is_singular = false;
    let reg = GeometryRegistry::global();
    let mut pnts: Vec<GpPnt> = Vec::new();
    // `cxx:2628-2639`: `GCPnts_UniformAbscissa Distribution(aCurve, N + 1)`,
    // then `for (i = 1; i <= N; i++)`.
    for e in shape_children(the_wire) {
        let Some(curve) = reg.edge_curve(&e) else {
            continue;
        };
        let Ok(params) = uniform_abscissa(curve.as_ref(), N + 1) else {
            continue;
        };
        for i in 1..=N {
            if let Some(&par) = params.get(i) {
                pnts.push(curve.d0(par));
            }
        }
    }

    // `cxx:2641-2648`: `GeomLib::AxeOfInertia(Pnts, Axe, IsSingular)` with the
    // declaration's default `Tol = 1.0e-7` (`GeomLib.hxx`).
    let (axe, singular) = geom_lib_axe_of_inertia(&pnts, 1.0e-7);
    *is_singular = singular;
    if singular {
        return -1.0;
    }

    // `cxx:2649`.
    *the_plane = GpPln::new(axe.to_ax3());
    // `cxx:2650-2657`.
    let mut max_deviation = 0.0f64;
    for p in &pnts {
        let dist = the_plane.signed_distance_pnt(p).abs();
        if dist > max_deviation {
            max_deviation = dist;
        }
    }
    max_deviation
}

/// `UpdateSectionEdge` (`cxx:2661-2705`): replace one endpoint of
/// `the_edge` by `the_vertex` and shorten the edge range to `the_param`.
///
/// `F_Edge = theEdge` shares the TShape, so every write lands on `the_edge`
/// exactly as in the reference; `the_vertex`'s orientation is set in place
/// (`cxx:2700`).
pub(super) fn update_section_edge(
    the_edge: &TopoShape,
    the_const_vertex: &TopoShape,
    the_vertex: &mut TopoShape,
    the_param: f64,
) {
    let reg = GeometryRegistry::global();
    // `cxx:2666-2667`: `F_Edge = theEdge; F_Edge.Orientation(TopAbs_FORWARD)`.
    let f_edge = forward(the_edge);
    let (v1, v2) = edge_vertices(&Edge(f_edge.clone()));
    // `cxx:2671-2684`: the vertex opposite `theConstVertex`.
    let (or_of_vertex, another_vertex) = match (&v1, &v2) {
        (Some(v1), Some(v2)) if is_same(&v1.0, the_const_vertex) => {
            (Orientation::Reversed, Some(v2.0.clone()))
        }
        (Some(v1), Some(_)) => (Orientation::Forward, Some(v1.0.clone())),
        _ => (Orientation::Forward, None),
    };

    // `cxx:2687-2698`: `BRep_Tool::Range` + `B.Begin`/`Range`.
    let (mut fpar, mut lpar) = reg.edge_parameters(&f_edge);
    if or_of_vertex == Orientation::Forward {
        fpar = the_param;
    } else {
        lpar = the_param;
    }
    reg.set_ranges_all(&f_edge, fpar, lpar);

    // `cxx:2699-2703`: `F_Edge.Free(true)`, `B.Remove(F_Edge, AnotherVertex)`,
    // `theVertex.Orientation(OrOfVertex)`, `B.Add(F_Edge, theVertex)`.
    f_edge.set_free(true);
    let builder = TopoBuilder::new();
    let mut target = f_edge.clone();
    if let Some(another) = another_vertex {
        builder.remove(&mut target, &another);
    }
    the_vertex.set_orientation(or_of_vertex);
    builder.add(&mut target, the_vertex);
}

/// `FindEdgeCloseToBisectorPlane` (`cxx:2708-2777`): the edge at
/// `the_vertex` whose direction / plane is closest to the bisector axis, taken
/// out of `the_comp`.
///
/// `the_comp` holds the compound's edges (the reference mutates a
/// `TopoDS_Compound`); `None` is OCCT's null edge.
pub(super) fn find_edge_close_to_bisector_plane(
    the_vertex: &TopoShape,
    the_comp: &mut Vec<TopoShape>,
    the_axis: &GpAx1,
) -> Option<Edge> {
    // `cxx:2711-2713`: `MapShapesAndAncestors(theComp, VERTEX, EDGE)`.
    let comp = TopoBuilder::new().make_compound_of(the_comp);
    let ve_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
    // `cxx:2717-2720`.
    let edges = ve_map.get(&shape_map_key(the_vertex))?;

    // `cxx:2725-2729`.
    let min_edge = if edges.len() == 1 {
        Edge(edges[0].clone())
    } else {
        // `cxx:2731-2774`.
        let reg = GeometryRegistry::global();
        let mut min_angle = f64::MAX;
        let mut best: Option<Edge> = None;
        for an_edge in edges {
            let wire = Wire(an_edge.clone());
            let mut a_pln = GpPln::new(GpAx2::standard().to_ax3());
            let mut issing = false;
            compute_average_plane_and_max_deviation(&wire.0, &mut a_pln, &mut issing);
            let an_angle = if issing {
                // `cxx:2744-2758`: the edge is a line segment.
                let Some(curve) = reg.edge_curve(an_edge) else {
                    continue;
                };
                let (first_par, last_par) = reg.edge_parameters(an_edge);
                let first_pnt = curve.d0(first_par);
                let last_pnt = curve.d0(last_par);
                let edge_vec = GpVec::from_pnts(&first_pnt, &last_pnt);
                let Ok(edge_dir) = GpDir::from_vec(&edge_vec) else {
                    continue;
                };
                let edge_axis = GpAx1::new(first_pnt, edge_dir);
                let mut angle = edge_axis.direction().angle(the_axis.direction());
                if angle > std::f64::consts::FRAC_PI_2 {
                    angle = std::f64::consts::PI - angle;
                }
                std::f64::consts::FRAC_PI_2 - angle
            } else {
                // `cxx:2759-2766`.
                let mut angle = a_pln.axis().direction().angle(the_axis.direction());
                if angle > std::f64::consts::FRAC_PI_2 {
                    angle = std::f64::consts::PI - angle;
                }
                angle
            };
            if an_angle < min_angle {
                min_angle = an_angle;
                best = Some(Edge(an_edge.clone()));
            }
        }
        best?
    };

    // `cxx:2775-2776`: `B.Remove(theComp, MinEdge)`.
    comp_remove(the_comp, &min_edge.0);
    Some(min_edge)
}

/// `FindMiddleEdges` (`cxx:2779-2819`): walk the section compound from
/// `the_vertex1` to `the_vertex2`, always taking the edge closest to the
/// bisector plane.
pub(super) fn find_middle_edges(
    the_vertex1: &TopoShape,
    the_vertex2: &TopoShape,
    the_axis: &GpAx1,
    the_comp: &mut Vec<TopoShape>,
    the_e_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:2786-2794`.
    let comp = TopoBuilder::new().make_compound_of(the_comp);
    let ve_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
    if ve_map.is_empty() {
        return false;
    }
    if !ve_map.contains_key(&shape_map_key(the_vertex1))
        || !ve_map.contains_key(&shape_map_key(the_vertex2))
    {
        return false;
    }

    // `cxx:2796-2818`.
    let mut cur_vertex = the_vertex1.clone();
    loop {
        let Some(cur_edge) = find_edge_close_to_bisector_plane(&cur_vertex, the_comp, the_axis)
        else {
            return false;
        };
        let (v1, v2) = edge_vertices(&cur_edge);
        let Some(v1) = v1 else {
            return false;
        };
        let Some(v2) = v2 else {
            return false;
        };
        cur_vertex = if is_same(&v1.0, &cur_vertex) {
            v2.0.clone()
        } else {
            v1.0.clone()
        };
        the_e_list.push(cur_edge.0.clone());
        if is_same(&cur_vertex, the_vertex2) {
            return true;
        }
    }
}

/// `FindCommonVertex(theFirstEdge, theLastEdge, ...)` (`cxx:2821-2855`): the
/// shared vertex of two section edges, going through
/// `TopExp::CommonVertex` unless both ends coincide.
pub(super) fn find_common_vertex_of_edges(
    the_first_edge: &Edge,
    the_last_edge: &Edge,
    the_first_vertex: Option<&TopoShape>,
    the_last_vertex: Option<&TopoShape>,
) -> Option<TopoShape> {
    // `cxx:2827-2831`: `!theFirstVertex.IsSame(theLastVertex)`.
    let same_ends = match (the_first_vertex, the_last_vertex) {
        (Some(a), Some(b)) => is_same(a, b),
        (None, None) => true,
        _ => false,
    };
    if !same_ends {
        return common_vertex(the_first_edge, the_last_edge);
    }

    // `cxx:2833-2853`.
    let first = the_first_vertex?;
    let (Some(v1), Some(v2)) = edge_vertices(the_first_edge) else {
        return None;
    };
    let (Some(v3), Some(v4)) = edge_vertices(the_last_edge) else {
        return None;
    };
    if is_same(&v1.0, first) {
        if is_same(&v2.0, &v3.0) || is_same(&v2.0, &v4.0) {
            return Some(v2.0.clone());
        }
    } else if is_same(&v1.0, &v3.0) || is_same(&v1.0, &v4.0) {
        return Some(v1.0.clone());
    }
    None
}

/// `ChooseSection` (`cxx:1020-1285`): the continuous section wire closest to
/// the bisector plane, with its average plane.
///
/// `None` is the reference's `false`; the third field is `IsSingular`.
pub(super) fn choose_section(
    the_transition: TransitionStyle,
    the_axe_of_bis_plane: &GpAx2,
    the_comp: &[TopoShape],
    the_first_vertex: Option<&TopoShape>,
    the_last_vertex: Option<&TopoShape>,
) -> Option<(TopoShape, GpPln, bool)> {
    // `cxx:1028`.
    let mut is_singular = false;
    let axis = the_axe_of_bis_plane.axis();

    // `cxx:1032-1137`: with `myTransition == BRepFill_Right` and both end
    // vertices known, the section wire runs from the first to the last one.
    if the_transition == TransitionStyle::Right
        && the_first_vertex.is_some()
        && the_last_vertex.is_some()
    {
        let first_v = the_first_vertex?;
        let last_v = the_last_vertex?;

        // `cxx:1039-1045`.
        let mut old_comp: Vec<TopoShape> = the_comp.to_vec();

        // `cxx:1047-1052`.
        let first_edge = find_edge_close_to_bisector_plane(first_v, &mut old_comp, axis)?;

        // `cxx:1054-1059`: the compound must not run dry before the last edge.
        if old_comp.is_empty() {
            if let Some(v) = the_comp.first() {
                old_comp.push(v.clone());
            }
        }

        // `cxx:1060-1064`.
        let last_edge = find_edge_close_to_bisector_plane(last_v, &mut old_comp, axis)?;

        // `cxx:1072`.
        let mut new_wire: Vec<TopoShape> = vec![first_edge.0.clone()];

        // `cxx:1074-1132`.
        if !is_same(&first_edge.0, &last_edge.0) {
            let common =
                find_common_vertex_of_edges(&first_edge, &last_edge, Some(first_v), Some(last_v));
            if common.is_some() {
                // `cxx:1079-1082`.
                new_wire.push(last_edge.0.clone());
            } else {
                // `cxx:1084-1089`: the far end of each edge.
                let (fv1, fv2) = edge_vertices(&first_edge);
                let (fv1, fv2) = (fv1?, fv2?);
                let vertex1 = if is_same(&fv1.0, first_v) {
                    fv2.0.clone()
                } else {
                    fv1.0.clone()
                };
                let (lv1, lv2) = edge_vertices(&last_edge);
                let (lv1, lv2) = (lv1?, lv2?);
                let vertex2 = if is_same(&lv1.0, last_v) {
                    lv2.0.clone()
                } else {
                    lv1.0.clone()
                };

                // `cxx:1090-1101`.
                if find_middle_edges(&vertex1, &vertex2, axis, &mut old_comp, &mut new_wire) {
                    new_wire.push(last_edge.0.clone());
                } else {
                    // `cxx:1103-1131`: trim both edges at their closest points
                    // and join them with a vertex at the midpoint.
                    let sols = extrema_of_edges(&first_edge, &last_edge);
                    if !sols.is_empty() {
                        let mut imin = 0usize;
                        for i in 1..sols.len() {
                            let d_i = sols[i].0.square_distance(&sols[i].1);
                            let d_min = sols[imin].0.square_distance(&sols[imin].1);
                            if d_i < d_min {
                                imin = i;
                            }
                        }
                        let (p_on_first, p_on_last, par_first, par_last) = sols[imin].clone();
                        let a_min_dist = p_on_first.distance(&p_on_last);
                        let mid = GpPnt::new(
                            (p_on_first.x() + p_on_last.x()) * 0.5,
                            (p_on_first.y() + p_on_last.y()) * 0.5,
                            (p_on_first.z() + p_on_last.z()) * 0.5,
                        );
                        // `cxx:1121-1122`.
                        let mut common_vertex = TopoBuilder::new().make_vertex(mid, CONFUSION);
                        common_vertex.set_tolerance(1.001 * a_min_dist / 2.0);
                        // `cxx:1124-1125`.
                        update_section_edge(
                            &first_edge.0,
                            first_v,
                            &mut common_vertex.0,
                            par_first,
                        );
                        update_section_edge(
                            &last_edge.0,
                            last_v,
                            &mut common_vertex.0,
                            par_last,
                        );
                        // `cxx:1127-1129`.
                        new_wire.push(last_edge.0.clone());
                    }
                }
            }
        }

        // `cxx:1134-1136`.
        let edges: Vec<Edge> = new_wire.iter().map(|s| Edge(s.clone())).collect();
        let wire = wire_from_edges(&edges)?;
        return Some((
            wire.0,
            GpPln::new(the_axe_of_bis_plane.to_ax3()),
            is_singular,
        ));
    }

    // `cxx:1139-1148`: the general case walks the section compound edge by edge.
    let mut old_comp: Vec<TopoShape> = the_comp.to_vec();
    let mut an_error = false;
    let mut wseq: Vec<TopoShape> = Vec::new();

    loop {
        // `cxx:1156-1160`.
        let Some(first_shape) = old_comp.first().cloned() else {
            break;
        };
        // `cxx:1161-1163`.
        let first_edge = Edge(first_shape.clone());
        let Some(mut new_wire) = wire_from_edges(&[first_edge.clone()]) else {
            an_error = true;
            break;
        };
        comp_remove(&mut old_comp, &first_shape);
        // `cxx:1164-1168`.
        if new_wire.0.is_closed() {
            wseq.push(new_wire.0);
            continue;
        }

        loop {
            // `cxx:1171-1177`.
            let (ext0, ext1) = wire_vertices(&new_wire);
            let (Some(ext0), Some(ext1)) = (ext0, ext1) else {
                an_error = true;
                break;
            };
            // `cxx:1178-1190`.
            let comp = TopoBuilder::new().make_compound_of(&old_comp);
            let ve_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
            let mut vedges: [Vec<TopoShape>; 2] = [Vec::new(), Vec::new()];
            for (j, ext) in [&ext0, &ext1].into_iter().enumerate() {
                if let Some(list) = ve_map.get(&shape_map_key(&ext.0)) {
                    vedges[j] = list.clone();
                }
            }
            // `cxx:1191-1196`.
            if vedges[0].is_empty() && vedges[1].is_empty() {
                break;
            }
            // `cxx:1197-1207`: a single continuation is appended directly.
            let mut modified = false;
            for j in 0..2 {
                if vedges[j].len() == 1 {
                    let an_edge = Edge(vedges[j][0].clone());
                    let Some(w) = wire_with_edge(&new_wire.0, &an_edge) else {
                        an_error = true;
                        break;
                    };
                    new_wire = w;
                    comp_remove(&mut old_comp, &an_edge.0);
                    modified = true;
                }
            }
            if an_error {
                break;
            }
            // `cxx:1208-1233`: several continuations - take the one whose
            // average plane deviates least.
            if !modified {
                let ind = if vedges[0].is_empty() { 1 } else { 0 };
                let edges = vedges[ind].clone();
                let mut the_ind = 0usize;
                let mut min_deviation = f64::MAX;
                for (j, cand) in edges.iter().enumerate() {
                    let Some(a_wire) = wire_with_edge(&new_wire.0, &Edge(cand.clone())) else {
                        break;
                    };
                    let mut a_pln = GpPln::new(GpAx2::standard().to_ax3());
                    let mut issing = false;
                    let deviation =
                        compute_average_plane_and_max_deviation(&a_wire.0, &mut a_pln, &mut issing);
                    if deviation < min_deviation {
                        min_deviation = deviation;
                        the_ind = j;
                    }
                }
                let Some(w) = wire_with_edge(&new_wire.0, &Edge(edges[the_ind].clone())) else {
                    an_error = true;
                    break;
                };
                new_wire = w;
                // `cxx:1232-1233`.
                comp_remove(&mut old_comp, &edges[the_ind]);
            }
            // `cxx:1234-1237`.
            if new_wire.0.is_closed() {
                break;
            }
        }
        // `cxx:1238`.
        wseq.push(new_wire.0);
        if an_error {
            break;
        }
    }

    // `cxx:1246-1248`.
    let mut min_angle = f64::MAX;
    if an_error || !old_comp.is_empty() {
        return None;
    }
    // `cxx:1249-1254`: one wire is the answer.
    if wseq.len() == 1 {
        let res_wire = wseq[0].clone();
        let mut res_plane = GpPln::new(GpAx2::standard().to_ax3());
        compute_average_plane_and_max_deviation(&res_wire, &mut res_plane, &mut is_singular);
        return Some((res_wire, res_plane, is_singular));
    }
    // `cxx:1255-1281`: keep the wire whose average plane is closest to the
    // bisector plane.
    let mut res_wire: Option<TopoShape> = None;
    let mut res_plane = GpPln::new(GpAx2::standard().to_ax3());
    for a_wire in &wseq {
        let mut a_pln = GpPln::new(GpAx2::standard().to_ax3());
        let mut issing = false;
        compute_average_plane_and_max_deviation(a_wire, &mut a_pln, &mut issing);
        if issing {
            continue;
        }
        // `cxx:1266-1273`.
        let mut angle = a_pln.axis().direction().angle(axis.direction());
        if angle > std::f64::consts::FRAC_PI_2 {
            angle = std::f64::consts::PI - angle;
        }
        if angle < min_angle {
            min_angle = angle;
            res_wire = Some(a_wire.clone());
            res_plane = a_pln;
        }
    }
    // `cxx:1249-1283`: with an empty `Wseq` the reference still answers `true`
    // and hands back its untouched (null) wire; the port reports the empty
    // wire so `MakeFacesSec` sees the same empty section compound.
    Some((
        res_wire.unwrap_or_else(|| Wire::new().0),
        res_plane,
        is_singular,
    ))
}

/// `StoreVedgeInHistMap` (`cxx:1393-1407`): the `v`-iso edge of the face
/// starts with a default (empty) entry in the map, then the intersection of
/// the two consecutive faces is appended to it.
pub(super) fn store_vedge_in_hist_map(
    the_v_edges: &[Option<TopoShape>],
    the_index: usize,
    the_new_vshape: &TopoShape,
    the_hist_map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) {
    let Some(Some(a_v_edge)) = the_v_edges.get(the_index) else {
        return;
    };
    the_hist_map
        .entry(shape_map_key(a_v_edge))
        .or_default()
        .push(the_new_vshape.clone());
}

/// `FindFreeVertices(theShape, theVerticesToAvoid, theListOfVertex)`
/// (`cxx:1409-1435`): the vertices of `theShape` outside
/// `the_vertices_to_avoid` whose ancestor count is below two.
///
/// `TopExp::MapShapesAndAncestors(theShape, TopAbs_VERTEX, TopAbs_EDGE, aMap)`
/// is the non-unique variant: every `(vertex, edge)` occurrence is one list
/// entry, so an edge whose two endpoints are the same vertex contributes two.
pub(super) fn find_free_vertices(
    the_shape: &TopoShape,
    the_vertices_to_avoid: &[TopoShape],
    the_list_of_vertex: &mut Vec<TopoShape>,
) {
    the_list_of_vertex.clear();
    // Insertion-ordered ancestor counts, matching the `IndexedDataMap` walk of
    // the reference.
    let mut counts: Vec<(TopoShape, usize)> = Vec::new();
    let bump = |v: &TopoShape, counts: &mut Vec<(TopoShape, usize)>| {
        let key = shape_map_key(v);
        if let Some(slot) = counts
            .iter_mut()
            .find(|(s, _)| shape_map_key(s) == key)
        {
            slot.1 += 1;
        } else {
            counts.push((v.clone(), 1));
        }
    };
    for e in explorer_shapes(the_shape, ShapeType::Edge) {
        for v in explorer_shapes(&e, ShapeType::Vertex) {
            bump(&v, &mut counts);
        }
    }
    for (v, n) in counts {
        if the_vertices_to_avoid.iter().any(|a| is_same(a, &v)) {
            continue;
        }
        if n < 2 {
            the_list_of_vertex.push(v);
        }
    }
}

/// `FindCommonVertex(theDS, theEIndex1, theEIndex2, theCrossDirection,
/// theCommonVertex, theParamOnE1, theParamOnE2)` (`cxx:1440-1508`).
///
/// Walks the `InterfEE` records of the pair, keeps the ones whose
/// `IntTools_CommonPrt::Type()` is `TopAbs_VERTEX`, takes the vertex parameters
/// through `IntTools_Tools::VertexParameters` and accepts the record only when
/// the cross product of the two edge tangents follows `theCrossDirection`
/// (`cxx:1494-1501`).
pub(super) fn find_common_vertex(
    the_ds: &BopdsDS,
    the_e_index1: usize,
    the_e_index2: usize,
    the_cross_direction: &GpVec,
    the_common_vertex: &mut Option<TopoShape>,
    the_param_on_e1: &mut f64,
    the_param_on_e2: &mut f64,
) -> bool {
    let Some(e1s) = the_ds.shape(the_e_index1).cloned() else {
        return false;
    };
    let Some(e2s) = the_ds.shape(the_e_index2).cloned() else {
        return false;
    };
    let (Some(c1), Some(c2)) = (
        BRepTool::edge_curve(&Edge(e1s)),
        BRepTool::edge_curve(&Edge(e2s)),
    ) else {
        return false;
    };

    for ee in the_ds.interf_ee() {
        if !((the_e_index1 == ee.index1() && the_e_index2 == ee.index2())
            || (the_e_index1 == ee.index2() && the_e_index2 == ee.index1()))
        {
            continue;
        }
        // `cxx:1469-1472`: records without a new vertex are skipped.
        let Some(n_new) = ee.new_index() else {
            continue;
        };
        let Some(cp) = ee.common_part() else {
            continue;
        };
        if cp.part_type() != CommonPartType::Vertex {
            continue;
        }
        let Some(vs) = the_ds.shape(n_new).cloned() else {
            continue;
        };

        // `cxx:1478-1484`: `IntTools_Tools::VertexParameters(aCP, t1, t2)`,
        // with the argument order following `Index1()`.
        let range1 = cp.range;
        let range2 = cp.ranges2().first().copied().unwrap_or_default();
        let vp1 = cp.vertex_parameter1().unwrap_or(f64::NAN);
        let vp2 = cp.vertex_parameter2().unwrap_or(f64::NAN);
        let (t1, t2) = if the_e_index1 == ee.index1() {
            crate::int_tools_curve_box::curve_box_vertex_parameters(range1, vp1, range2, vp2)
        } else {
            crate::int_tools_curve_box::curve_box_vertex_parameters(range2, vp2, range1, vp1)
        };

        // `cxx:1488-1501`: the intersection direction must follow the path.
        let dir1 = c1.d1(t1).1;
        let dir2 = c2.d1(t2).1;
        let cross = dir1.crossed(&dir2);
        if cross.dot(the_cross_direction) > SQUARE_CONFUSION {
            *the_common_vertex = Some(vs);
            *the_param_on_e1 = t1;
            *the_param_on_e2 = t2;
            return true;
        }
    }
    false
}


