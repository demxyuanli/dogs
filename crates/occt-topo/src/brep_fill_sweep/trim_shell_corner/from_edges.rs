//! `FindFromUEdge` / `FindFromVEdge` and the two section-pair filters.
//!
//! Moved here verbatim from the single-file port; the `BRepFill_TrimShellCorner.cxx`
//! branch structure is unchanged.

use super::*;
use super::helpers::*;
use super::section_walk::*;

/// `FindFromUEdge` outputs (`cxx:1962-2070`).
pub(super) struct FromUEdge {
    /// `theSecEdgesNew`.
    pub(super) comp_of_sec_edges: Vec<TopoShape>,
    /// `theListOfWireEdges`.
    pub(super) list_of_wire_edges: Vec<TopoShape>,
    /// `theFoundPave`.
    pub(super) found_pave: BopdsPave,
    /// `isOnUEdge`.
    pub(super) is_on_u_edge: bool,
}

/// `FindFromUEdge` (`cxx:1962-2070`): the section edges that bridge the first
/// `u`-iso edge's end to the first reachable pave of the bound edge.
///
/// `None` is the reference's `false`. The `cxx:1998-2000` early `return true`
/// leaves the outputs at the caller's values, i.e. the incoming section
/// compound and an empty wire-edge list, with `theFoundPave` reset to `(0, 0)`.
#[allow(clippy::too_many_arguments)]
pub(super) fn find_from_u_edge(
    the_ue1_old: &TopoShape,
    the_ue2_old: &TopoShape,
    the_ue1_new: &TopoShape,
    the_ue2_new: &TopoShape,
    the_face: &TopoShape,
    the_sec_edges: &[TopoShape],
    the_rank: usize,
    the_bound_edge: &TopoShape,
    the_bound_edge_index: usize,
    the_ds: &BopdsDS,
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_sec_edges_in: Vec<TopoShape>,
    the_list_in: Vec<TopoShape>,
) -> Option<FromUEdge> {
    // `cxx:1980-1982`.
    let mut found_pave = BopdsPave::new(0, 0.0);
    let mut is_on_u_edge = true;

    // `cxx:1984-1988`.
    let mut a_face_f = the_face.clone();
    a_face_f.set_orientation(Orientation::Forward);
    let mut comp_of_sec_edges = the_sec_edges_in;
    let mut list_of_wire_edges = the_list_in;

    // `cxx:1992-1996`.
    let mut a_pave1 = BopdsPave::new(0, 0.0);
    let mut a_pave2 = BopdsPave::new(0, 0.0);
    let mut a_prev_vertex: Option<TopoShape> = None;
    let mut a_next_vertex: Option<TopoShape> = None;
    let mut atmp_pave = BopdsPave::new(0, 0.0);

    // `cxx:1998-2000`: an unrecorded `u`-iso edge leaves every output alone.
    if !find_vertex(
        the_ue1_old,
        the_rank,
        the_ds,
        the_hist_map,
        &mut a_prev_vertex,
        &mut atmp_pave,
    ) {
        return Some(FromUEdge {
            comp_of_sec_edges,
            list_of_wire_edges,
            found_pave,
            is_on_u_edge,
        });
    }
    // `cxx:2002-2005`.
    let Some(a_prev) = a_prev_vertex.clone() else {
        return None;
    };

    // `cxx:2007-2010`.
    let a_first_v = a_prev.clone();
    let mut b_sec_found = false;
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(&a_face_f);
    let Some(a_c1) = reg.edge_pcurve(the_ue1_new, key) else {
        return None;
    };
    let (c1_f, c1_l) = reg.pcurve_range(the_ue1_new, key).unwrap_or((0.0, 0.0));
    let mut p1 = if the_rank == 0 {
        a_c1.d0(c1_l)
    } else {
        a_c1.d0(c1_f)
    };

    // `cxx:2013-2017`.
    let nbpave = the_ds.paves(the_bound_edge_index).len();
    let mut pit = 0usize;

    // `cxx:2018-2039`: walk the bound edge's paves.
    let Some(a_c2) = reg.edge_pcurve(the_bound_edge, key) else {
        return None;
    };
    while find_next_vertex(
        the_bound_edge_index,
        &a_pave1,
        the_ds,
        &mut a_next_vertex,
        &mut a_pave2,
    ) && pit < nbpave
    {
        let Some(a_last_v) = a_next_vertex.clone() else {
            return None;
        };
        let p2 = a_c2.d0(a_pave2.parameter());
        let mut a_ordered_list: Vec<TopoShape> = Vec::new();

        if fill_gap(
            &a_first_v,
            &a_last_v,
            &p1,
            &p2,
            &a_face_f,
            &comp_of_sec_edges,
            &mut a_ordered_list,
        ) {
            comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
            list_of_wire_edges.extend(a_ordered_list);
            found_pave = a_pave2;
            is_on_u_edge = false;
            b_sec_found = true;
            break;
        }
        if let Some(v) = a_next_vertex.clone() {
            a_prev_vertex = Some(v);
        }
        a_pave1 = a_pave2;
        pit += 1;
    }

    // `cxx:2041-2060`: fall back to the second `u`-iso edge.
    if !b_sec_found {
        let mut a_next2: Option<TopoShape> = None;
        let mut a_pave2b = BopdsPave::new(0, 0.0);
        if find_vertex(
            the_ue2_old,
            the_rank,
            the_ds,
            the_hist_map,
            &mut a_next2,
            &mut a_pave2b,
        ) {
            if let Some(a_last_v) = a_next2.clone() {
                let Some(a_c3) = reg.edge_pcurve(the_ue2_new, key) else {
                    return None;
                };
                let p2 = a_c3.d0(a_pave2b.parameter());
                let mut a_ordered_list: Vec<TopoShape> = Vec::new();
                if fill_gap(
                    &a_first_v,
                    &a_last_v,
                    &p1,
                    &p2,
                    &a_face_f,
                    &comp_of_sec_edges,
                    &mut a_ordered_list,
                ) {
                    comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                    list_of_wire_edges.extend(a_ordered_list);
                    found_pave = a_pave2b;
                    b_sec_found = true;
                    is_on_u_edge = true;
                }
            }
        }
    }

    // `cxx:2062-2068`: without a section the outputs stay as they came in.
    if !b_sec_found {
        comp_of_sec_edges = the_sec_edges.to_vec();
        list_of_wire_edges.clear();
    }
    Some(FromUEdge {
        comp_of_sec_edges,
        list_of_wire_edges,
        found_pave,
        is_on_u_edge,
    })
}

/// `FindFromVEdge` outputs (`cxx:2073-2513`).
pub(super) struct FromVEdge {
    /// `theListOfWireEdges`.
    pub(super) list_of_wire_edges: Vec<TopoShape>,
    /// `isSectionFound`.
    pub(super) is_section_found: bool,
}

/// `FindFromVEdge` (`cxx:2073-2513`): the section edges that complete the wire
/// from the pave found by `FindFromUEdge` down to the bound edge's last pave.
///
/// The walk exists twice in the reference: a "by pairs non continuously" pass
/// (`cxx:2144-2226`) that collects candidate `(pave1, pave2, flag)` triples and
/// a "by pairs continuously" pass (`cxx:2228-2511`) that consumes them. Both
/// are reproduced in order; `anArrayOfListOfSec` stays 1-based as in OCCT.
#[allow(clippy::too_many_arguments)]
pub(super) fn find_from_v_edge(
    the_prev_pave: &BopdsPave,
    is_on_u_edge: bool,
    the_ue1_old: &TopoShape,
    the_ue2_old: &TopoShape,
    the_face: &TopoShape,
    the_sec_edges: &[TopoShape],
    the_rank: usize,
    the_bound_edge: &TopoShape,
    the_bound_edge_index: usize,
    the_ds: &BopdsDS,
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> Option<FromVEdge> {
    // `cxx:2090-2097`.
    let mut list_of_wire_edges: Vec<TopoShape> = Vec::new();
    let mut is_section_found = false;
    let mut a_face_f = the_face.clone();
    a_face_f.set_orientation(Orientation::Forward);
    let mut comp_of_sec_edges: Vec<TopoShape> = the_sec_edges.to_vec();

    // `cxx:2102-2113`: when the previous pave sits on the bound edge, a `u`-iso
    // end matching it means the walk is already closed.
    if is_on_u_edge {
        let mut atmp_vertex: Option<TopoShape> = None;
        let mut a_pave_of_e2 = BopdsPave::new(0, 0.0);
        if find_vertex(
            the_ue2_old,
            the_rank,
            the_ds,
            the_hist_map,
            &mut atmp_vertex,
            &mut a_pave_of_e2,
        ) && the_prev_pave.is_equal(&a_pave_of_e2)
        {
            return Some(FromVEdge {
                list_of_wire_edges,
                is_section_found,
            });
        }
    }

    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(&a_face_f);
    let Some(a_c2) = reg.edge_pcurve(the_bound_edge, key) else {
        return None;
    };

    // `cxx:2123-2136`.
    let mut a_pave1 = *the_prev_pave;
    if is_on_u_edge {
        let mut atmp_pave = BopdsPave::new(0, 0.0);
        if !get_pave(the_bound_edge_index, true, the_ds, &mut atmp_pave) {
            return None;
        }
        a_pave1 = atmp_pave;
    }
    let mut p1 = a_c2.d0(a_pave1.parameter());
    let mut a_prev_vertex: Option<TopoShape> = the_ds.shape(a_pave1.index()).cloned();

    // `cxx:2138-2142`.
    let nbpave = the_ds.paves(the_bound_edge_index).len();
    let mut pit = 0usize;
    // `NCollection_Array1<NCollection_List<TopoDS_Shape>>(1, nbpave)`.
    let mut an_array_of_list_of_sec: Vec<Vec<TopoShape>> = vec![Vec::new(); nbpave + 1];

    // `by pairs non continuously. begin` (`cxx:2144-2226`).
    let mut a_first_pave = a_pave1;
    let mut a_first_vertex = a_prev_vertex.clone();
    let mut apfirst = p1;
    let mut a_first_paves: Vec<BopdsPave> = Vec::new();
    let mut a_last_paves: Vec<BopdsPave> = Vec::new();
    let mut a_list_of_flags: Vec<i32> = Vec::new();
    let mut apaircounter = 1usize;

    for _k in 0..nbpave {
        a_pave1 = a_first_pave;
        p1 = apfirst;
        a_prev_vertex = a_first_vertex.clone();
        let mut bfound = false;
        pit = 0;
        let mut a_pave2 = BopdsPave::new(0, 0.0);
        let mut a_next_vertex: Option<TopoShape> = None;

        while find_next_vertex(
            the_bound_edge_index,
            &a_pave1,
            the_ds,
            &mut a_next_vertex,
            &mut a_pave2,
        ) && pit < nbpave
        {
            let Some(a_first_v) = a_prev_vertex.clone() else {
                return None;
            };
            let Some(a_last_v) = a_next_vertex.clone() else {
                return None;
            };
            let p2 = a_c2.d0(a_pave2.parameter());
            let mut a_ordered_list: Vec<TopoShape> = Vec::new();

            if fill_gap(
                &a_first_v,
                &a_last_v,
                &p1,
                &p2,
                &a_face_f,
                &comp_of_sec_edges,
                &mut a_ordered_list,
            ) {
                comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                if apaircounter <= nbpave {
                    an_array_of_list_of_sec[apaircounter].extend(a_ordered_list);
                }
                apaircounter += 1;
                a_first_paves.push(a_first_pave);
                a_last_paves.push(a_pave2);
                a_list_of_flags.push(1);
                a_first_pave = a_pave2;
                a_first_vertex = a_next_vertex.clone();
                apfirst = p2;
                a_prev_vertex = a_next_vertex.clone();
                is_section_found = true;
                bfound = true;
            }
            a_pave1 = a_pave2;
            pit += 1;
        }

        // `cxx:2190-2210`: the closing pair against the second `u`-iso edge.
        let mut a_pave2v = BopdsPave::new(0, 0.0);
        let mut a_next2: Option<TopoShape> = None;
        if find_vertex(
            the_ue2_old,
            the_rank,
            the_ds,
            the_hist_map,
            &mut a_next2,
            &mut a_pave2v,
        ) {
            if let (Some(a_first_v), Some(a_last_v)) =
                (a_prev_vertex.clone(), a_next2.clone())
            {
                let Some(a_c3) = reg.edge_pcurve(the_ue2_old, key) else {
                    return None;
                };
                let p2 = a_c3.d0(a_pave2v.parameter());
                let mut a_ordered_list: Vec<TopoShape> = Vec::new();
                if fill_gap(
                    &a_first_v,
                    &a_last_v,
                    &p1,
                    &p2,
                    &a_face_f,
                    &comp_of_sec_edges,
                    &mut a_ordered_list,
                ) {
                    comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                    if apaircounter <= nbpave {
                        an_array_of_list_of_sec[apaircounter].extend(a_ordered_list);
                    }
                    apaircounter += 1;
                    a_first_paves.push(a_first_pave);
                    a_last_paves.push(a_pave2v);
                    a_list_of_flags.push(0);
                    is_section_found = true;
                    break;
                }
            }
        }

        // `cxx:2212-2224`: advance the anchor pave.
        if !bfound {
            let mut a_next3: Option<TopoShape> = None;
            let mut a_pave2c = BopdsPave::new(0, 0.0);
            if !find_next_vertex(
                the_bound_edge_index,
                &a_first_pave,
                the_ds,
                &mut a_next3,
                &mut a_pave2c,
            ) {
                break;
            }
            a_first_pave = a_pave2c;
            apfirst = a_c2.d0(a_pave2c.parameter());
            a_first_vertex = a_next3;
        }
    }
    // `by pairs non continuously. end`

    // `by pairs continuously. begin` (`cxx:2228-2511`).
    a_pave1 = *the_prev_pave;
    if is_on_u_edge {
        let mut atmp_pave = BopdsPave::new(0, 0.0);
        if !get_pave(the_bound_edge_index, true, the_ds, &mut atmp_pave) {
            return None;
        }
        a_pave1 = atmp_pave;
    }
    p1 = a_c2.d0(a_pave1.parameter());
    a_prev_vertex = the_ds.shape(a_pave1.index()).cloned();

    pit = 0;
    let mut a_next_vertex: Option<TopoShape> = None;
    let mut a_pave2 = BopdsPave::new(0, 0.0);
    while find_next_vertex(
        the_bound_edge_index,
        &a_pave1,
        the_ds,
        &mut a_next_vertex,
        &mut a_pave2,
    ) && pit < nbpave
    {
        let Some(a_first_v) = a_prev_vertex.clone() else {
            return None;
        };
        let Some(a_last_v) = a_next_vertex.clone() else {
            return None;
        };
        let p2 = a_c2.d0(a_pave2.parameter());

        // `cxx:2250-2318`.
        let (bisinside, apbindex) = pair_inside_continuous(
            &a_pave1,
            &a_pave2,
            &p2,
            &a_first_paves,
            &a_last_paves,
            &a_list_of_flags,
            the_ue2_old,
            &a_face_f,
            key,
            the_ds,
        )?;

        if !bisinside {
            let mut a_ordered_list: Vec<TopoShape> = Vec::new();
            if fill_gap(
                &a_first_v,
                &a_last_v,
                &p1,
                &p2,
                &a_face_f,
                &comp_of_sec_edges,
                &mut a_ordered_list,
            ) {
                comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                list_of_wire_edges.extend(a_ordered_list);
                is_section_found = true;
            } else {
                // `cxx:2333-2353`: fall back to the split edge of the pave pair.
                if let Some(a_e_split) = find_split_edge(
                    the_bound_edge_index,
                    &a_pave1,
                    &a_pave2,
                    the_ds,
                ) {
                    list_of_wire_edges.push(a_e_split);
                }
            }
        } else if apbindex > 0 {
            // `cxx:2357-2362`.
            if let Some(a_list_of_sec) = an_array_of_list_of_sec.get(apbindex) {
                list_of_wire_edges.extend(a_list_of_sec.clone());
            }
        }
        a_pave1 = a_pave2;
        a_prev_vertex = a_next_vertex.clone();
        p1 = p2;
        pit += 1;
    }

    // `cxx:2368-2478`: the closing step against the second `u`-iso edge.
    let mut a_next2: Option<TopoShape> = None;
    let mut a_pave2v = BopdsPave::new(0, 0.0);
    if find_vertex(
        the_ue2_old,
        the_rank,
        the_ds,
        the_hist_map,
        &mut a_next2,
        &mut a_pave2v,
    ) {
        if let (Some(a_first_v), Some(a_last_v)) = (a_prev_vertex.clone(), a_next2.clone()) {
            let Some(a_c3) = reg.edge_pcurve(the_ue2_old, key) else {
                return None;
            };
            let p2 = a_c3.d0(a_pave2v.parameter());

            // `cxx:2377-2426`: the `flag == 0` arm here is unconditional
            // (`cxx:2410-2414`), unlike the in-loop variant above.
            let (bisinside, apbindex) = pair_inside_closing(
                &a_pave1,
                &a_pave2v,
                &a_first_paves,
                &a_last_paves,
                &a_list_of_flags,
            );

            if !bisinside {
                let mut a_ordered_list: Vec<TopoShape> = Vec::new();
                if fill_gap(
                    &a_first_v,
                    &a_last_v,
                    &p1,
                    &p2,
                    &a_face_f,
                    &comp_of_sec_edges,
                    &mut a_ordered_list,
                ) {
                    comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                    list_of_wire_edges.extend(a_ordered_list);
                    is_section_found = true;
                } else {
                    // `cxx:2442-2468`: the last pave replaces `aPave2`.
                    let mut a_last_pave = BopdsPave::new(0, 0.0);
                    if !get_pave(the_bound_edge_index, false, the_ds, &mut a_last_pave) {
                        return None;
                    }
                    if let Some(a_e_split) =
                        find_split_edge(the_bound_edge_index, &a_pave1, &a_last_pave, the_ds)
                    {
                        list_of_wire_edges.push(a_e_split);
                    }
                }
            } else if apbindex > 0 {
                if let Some(a_list_of_sec) = an_array_of_list_of_sec.get(apbindex) {
                    list_of_wire_edges.extend(a_list_of_sec.clone());
                }
            }
        }
    } else {
        // `cxx:2479-2509`.
        let mut a_last_pave = BopdsPave::new(0, 0.0);
        if !get_pave(the_bound_edge_index, false, the_ds, &mut a_last_pave) {
            return None;
        }
        if let Some(a_e_split) =
            find_split_edge(the_bound_edge_index, &a_pave1, &a_last_pave, the_ds)
        {
            list_of_wire_edges.push(a_e_split);
        }
    }

    Some(FromVEdge {
        list_of_wire_edges,
        is_section_found,
    })
}

/// `FindFromVEdge`'s in-loop "is this pair already covered" test
/// (`cxx:2254-2315`).
///
/// Returns `(bisinside, apbindex)`; `None` is an unreachable null pcurve, which
/// the reference would dereference.
#[allow(clippy::too_many_arguments)]
pub(super) fn pair_inside_continuous(
    the_pave1: &BopdsPave,
    the_pave2: &BopdsPave,
    p2: &GpPnt2d,
    the_first_paves: &[BopdsPave],
    the_last_paves: &[BopdsPave],
    the_flags: &[i32],
    the_ue2_old: &TopoShape,
    the_face: &TopoShape,
    key: usize,
    the_ds: &BopdsDS,
) -> Option<(bool, usize)> {
    let mut bisinside = false;
    let mut apbindex = 0usize;
    let count = the_first_paves
        .len()
        .min(the_last_paves.len())
        .min(the_flags.len());

    // `cxx:2258-2262`: `apbcounter` starts at 1.
    for apbcounter in 1..=count {
        let a_pit1 = &the_first_paves[apbcounter - 1];
        let a_pit2 = &the_last_paves[apbcounter - 1];
        let a_flag = the_flags[apbcounter - 1];

        // `cxx:2267-2275`.
        let bfin = if the_pave1.is_equal(a_pit1) {
            true
        } else {
            the_pave1.parameter() > a_pit1.parameter()
        };

        // `cxx:2277-2305`.
        let blin = if a_flag != 0 {
            if the_pave2.is_equal(a_pit2) {
                true
            } else {
                the_pave2.parameter() < a_pit2.parameter()
            }
        } else if the_pave2.index() == a_pit2.index() && the_pave2.index() > 0 {
            let reg = GeometryRegistry::global();
            let pc = reg.edge_pcurve(the_ue2_old, key)?;
            let p3 = pc.d0(a_pit2.parameter());
            let a_v = the_ds.shape(the_pave2.index())?;
            let a_tolerance = surface_tolerance(the_face, BRepTool::vertex_tolerance(&Vertex(a_v.clone())));
            pnt2d_distance(p2, &p3) < a_tolerance
        } else {
            false
        };

        if bfin && blin {
            apbindex = apbcounter;
            bisinside = true;
            break;
        }
    }
    Some((bisinside, apbindex))
}

/// `FindFromVEdge`'s closing "is this pair already covered" test
/// (`cxx:2388-2425`): same shape as [`pair_inside_continuous`] except that the
/// `flag == 0` arm is unconditional (`cxx:2410-2414`).
pub(super) fn pair_inside_closing(
    the_pave1: &BopdsPave,
    the_pave2: &BopdsPave,
    the_first_paves: &[BopdsPave],
    the_last_paves: &[BopdsPave],
    the_flags: &[i32],
) -> (bool, usize) {
    let mut bisinside = false;
    let mut apbindex = 0usize;
    let count = the_first_paves
        .len()
        .min(the_last_paves.len())
        .min(the_flags.len());
    for apbcounter in 1..=count {
        let a_pit1 = &the_first_paves[apbcounter - 1];
        let a_pit2 = &the_last_paves[apbcounter - 1];
        let a_flag = the_flags[apbcounter - 1];

        let bfin = if the_pave1.is_equal(a_pit1) {
            true
        } else {
            the_pave1.parameter() > a_pit1.parameter()
        };
        let blin = if a_flag != 0 {
            if the_pave2.is_equal(a_pit2) {
                true
            } else {
                the_pave2.parameter() < a_pit2.parameter()
            }
        } else {
            true
        };
        if bfin && blin {
            apbindex = apbcounter;
            bisinside = true;
            break;
        }
    }
    (bisinside, apbindex)
}

/// The `aESplit` lookup of `FindFromVEdge` (`cxx:2333-2353`, `cxx:2450-2463`):
/// the split edge of the pave block of the bound edge whose ends are exactly
/// `(the_pave1, the_pave2)`.
pub(super) fn find_split_edge(
    the_bound_edge_index: usize,
    the_pave1: &BopdsPave,
    the_pave2: &BopdsPave,
    the_ds: &BopdsDS,
) -> Option<TopoShape> {
    for pb in the_ds.pave_blocks(the_bound_edge_index) {
        // `cxx:2338-2339`: `OriginalEdge() == theBoundEdgeIndex` plus both ends.
        if pb.original_edge() != the_bound_edge_index
            || !pb.pave1().is_equal(the_pave1)
            || !pb.pave2().is_equal(the_pave2)
        {
            continue;
        }
        // `cxx:2341-2346`: `aPB1->Edge() > 0`.
        if pb.has_edge() && pb.edge() > 0 {
            if let Some(s) = the_ds.shape(pb.edge()) {
                return Some(s.clone());
            }
        }
    }
    None
}
