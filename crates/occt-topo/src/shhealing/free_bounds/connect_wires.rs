//! `ConnectEdgesToWires` / `ConnectWiresToWires`
//! (`ShapeAnalysis_FreeBounds.cxx:135-453`).

use super::*;

/// `ShapeAnalysis_FreeBounds::ConnectEdgesToWires(edges, toler, shared)`
/// (`ShapeAnalysis_FreeBounds.cxx:135-163`).
pub(crate) fn connect_edges_to_wires(edges: &mut Vec<Edge>, toler: f64, shared: bool) -> Vec<Wire> {
    let builder = crate::builder::TopoBuilder::new();
    // `cxx:139-150`: one wire per edge.
    let mut iwires: Vec<Wire> = edges
        .iter()
        .map(|e| builder.make_wire(std::slice::from_ref(e)))
        .collect();

    // `cxx:152-153`.
    let mut vertices: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
    let wires = connect_wires_to_wires(&mut iwires, toler, shared, &mut vertices);

    // `cxx:155-161`: an edge whose wire was reversed is reversed back.
    for (i, edge) in edges.iter_mut().enumerate() {
        if iwires[i].0.orientation() == Orientation::Reversed {
            edge.0.reverse();
        }
    }

    wires
}

/// `ShapeAnalysis_FreeBounds::ConnectWiresToWires(iwires, toler, shared,
/// vertices)` (`ShapeAnalysis_FreeBounds.cxx:432-440`) through
/// `connectWiresToWiresImpl` (`cxx:201-427`).
pub(crate) fn connect_wires_to_wires(
    iwires: &mut Vec<Wire>,
    toler: f64,
    shared: bool,
    vertices: &mut HashMap<ShapeMapKey, TopoShape>,
) -> Vec<Wire> {
    if iwires.is_empty() {
        return Vec::new(); // `cxx:207-210`.
    }

    // `cxx:211-218`: `arrwires`.
    let mut arrwires: Vec<Wire> = iwires.clone();
    let mut owires: Vec<Wire> = Vec::new();
    let tolerance = toler.max(CONFUSION); // `cxx:220`.

    // `cxx:222-231`.
    let mut sewd = WireData::from_wire(&arrwires[0], false, true);
    let mut is_used_manifold_mode = true;
    if sewd.nb_edges() < 1 && sewd.nb_non_manifold_edges() > 0 {
        is_used_manifold_mode = false;
        sewd = WireData::from_wire(&arrwires[0], true, is_used_manifold_mode);
    }

    // `cxx:233-236`.
    let mut saw = WireAnalyzer::new();
    saw.set_precision(tolerance);

    // `cxx:238-257`.
    let mut tree = UbTree::new();
    let mut filler = UbTreeFiller::new();
    let mut sel = ShapeAnalysisBoxBndTreeSelector::new(shared);
    sel.load_list(1);
    for inb_w in 2..=arrwires.len() {
        let (v1, v2) = find_bounds(&arrwires[inb_w - 1].0);
        let (Some(v1), Some(v2)) = (v1, v2) else {
            continue;
        };
        let p1 = BRepTool::vertex_point(&v1);
        let p2 = BRepTool::vertex_point(&v2);
        let mut a_box = BndBox::new();
        a_box.add_point(&p1);
        a_box.add_point(&p2);
        a_box.set_gap(tolerance);
        filler.add(inb_w, a_box);
    }
    filler.fill(&mut tree);

    // `cxx:260-420`.
    let mut done = false;
    while !done {
        let mut found = false;
        let mut tail = false;
        let mut direct = false;
        let mut lwire: isize = 0;
        sel.set_stop();

        // `cxx:268-271`.
        let nb_edges = sewd.nb_edges();
        if nb_edges == 0 {
            // OCCT reads `sewd->Edge(1)` out of range here; the port stops.
            break;
        }
        let Some(first_edge) = sewd.edge(1) else {
            break;
        };
        let Some(last_edge) = sewd.edge(nb_edges as isize) else {
            break;
        };
        let vf = super::super::first_vertex(&first_edge);
        let vl = super::super::last_vertex(&last_edge);

        // `cxx:273-278`.
        let mut fv_box = BndBox::new();
        let mut lv_box = BndBox::new();
        if let Some(vf) = vf.as_ref() {
            let pf = BRepTool::vertex_point(vf);
            fv_box.add_point(&pf);
            fv_box.set_gap(tolerance);
        }
        if let Some(vl) = vl.as_ref() {
            let pl = BRepTool::vertex_point(vl);
            lv_box.add_point(&pl);
            lv_box.set_gap(tolerance);
        }
        sel.define_boxes(fv_box, lv_box);

        // `cxx:281-290`.
        if shared {
            sel.define_vertexes(vf.clone(), vl.clone());
        } else {
            let pf = vf
                .as_ref()
                .map(BRepTool::vertex_point)
                .unwrap_or_else(GpPnt::zero);
            let pl = vl
                .as_ref()
                .map(BRepTool::vertex_point)
                .unwrap_or_else(GpPnt::zero);
            sel.define_pnt(pf, pl);
            sel.set_tolerance(tolerance);
        }

        // `cxx:292-300`.
        let nsel = tree.select(&mut sel, &arrwires);
        if nsel != 0 && !sel.last_check_status(SE_FAIL) {
            found = true;
            lwire = sel.get_nb() as isize;
            tail = sel.last_check_status(SE_DONE1)
                || sel.last_check_status(super::super::box_bnd_tree::SE_DONE2);
            direct = sel.last_check_status(SE_DONE1)
                || sel.last_check_status(super::super::box_bnd_tree::SE_DONE3);
            sel.load_list(lwire as usize);
        }

        if found {
            // `cxx:303-317`.
            if !direct {
                arrwires[lwire as usize - 1].0.reverse();
            }
            let Some(acurwd) = Some(WireData::from_wire(
                &arrwires[lwire as usize - 1],
                true,
                is_used_manifold_mode,
            )) else {
                continue;
            };
            if acurwd.nb_edges() == 0 {
                continue; // `cxx:313-315`.
            }
            sewd.add_wiredata(&acurwd, if tail { 0 } else { 1 }); // `cxx:317`.
        } else {
            // `cxx:318-359`: heal the junctions of the wire being built.
            for i in 1..=sewd.nb_edges() {
                if !saw.check_connected(&sewd, i) {
                    continue;
                }
                let n2 = i;
                let n1 = if n2 > 1 { n2 - 1 } else { sewd.nb_edges() };
                let (Some(e1), Some(e2)) = (sewd.edge(n1 as isize), sewd.edge(n2 as isize)) else {
                    continue;
                };
                let vprev = super::super::last_vertex(&e1);
                let vfol = super::super::first_vertex(&e2);
                let (Some(vprev), Some(vfol)) = (vprev, vfol) else {
                    continue;
                };
                let v = if saw.last_check_status(SE_DONE1) {
                    vprev.clone() // `cxx:333-335`.
                } else {
                    combine_vertex(&vprev, &vfol, 1.0001) // `cxx:337-341`.
                };
                vertices.insert(shape_map_key(&vprev.0), v.0.clone());
                vertices.insert(shape_map_key(&vfol.0), v.0.clone());

                if sewd.nb_edges() < 2 {
                    // `cxx:346-349`.
                    let new_edge = copy_replace_vertices_with(&e2, Some(&v), Some(&v));
                    sewd.set(&new_edge, n2);
                } else {
                    // `cxx:350-357`.
                    let new_edge = copy_replace_vertices_with(&e2, Some(&v), None);
                    sewd.set(&new_edge, n2);
                    if !saw.last_check_status(SE_DONE1) {
                        let new_edge1 = copy_replace_vertices_with(&e1, None, Some(&v));
                        sewd.set(&new_edge1, n1);
                    }
                }
            }

            // `cxx:361-393`.
            let wire = sewd.wire();
            if is_used_manifold_mode {
                if !saw.check_connected(&sewd, 1) && saw.last_check_status(SE_OK) {
                    wire.0.set_closed(true);
                }
            } else {
                let mut vmap: Vec<TopoShape> = Vec::new();
                for child in cumulated_children(&wire.0) {
                    if child.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    for v in stored_children(&child) {
                        if v.shape_type() != ShapeType::Vertex {
                            continue;
                        }
                        if v.orientation() == Orientation::Forward
                            || v.orientation() == Orientation::Reversed
                        {
                            // `NCollection_Map::Add` toggling.
                            let key = shape_map_key(&v);
                            match vmap.iter().position(|s| shape_map_key(s) == key) {
                                Some(pos) => {
                                    vmap.remove(pos);
                                }
                                None => vmap.push(v),
                            }
                        }
                    }
                }
                if vmap.is_empty() {
                    wire.0.set_closed(true);
                }
            }

            // `cxx:395-420`.
            owires.push(wire);
            sewd.clear();
            sewd.set_manifold_mode(is_used_manifold_mode);

            lwire = -1;
            for i in 1..=arrwires.len() {
                if !sel.cont_wire(i) {
                    lwire = i as isize;
                    sewd.add_wire(&arrwires[i - 1], 0);
                    sel.load_list(i);
                    if sewd.nb_edges() > 0 {
                        break;
                    }
                    sewd.clear();
                }
            }
            if lwire == -1 {
                done = true;
            }
        }
    }

    // `cxx:422-426`.
    for i in 1..=iwires.len() {
        iwires[i - 1] = arrwires[i - 1].clone();
    }
    owires
}
