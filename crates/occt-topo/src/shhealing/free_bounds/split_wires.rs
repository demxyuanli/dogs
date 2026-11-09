//! `SplitWire` / `SplitWires` (`ShapeAnalysis_FreeBounds.cxx:455-606`).

use super::*;

/// The file-local `SplitWire` (`ShapeAnalysis_FreeBounds.cxx:455-587`).
fn split_wire(
    wire: &Wire,
    toler: f64,
    shared: bool,
    closed: &mut Vec<Wire>,
    open: &mut Vec<Wire>,
) {
    let tolerance = toler.max(CONFUSION); // `cxx:462`.
    let builder = crate::builder::TopoBuilder::new();

    let sewd = WireData::from_wire(wire, false, true); // `cxx:468`.
    let nbedges = sewd.nb_edges(); // `cxx:469`.

    // `cxx:473-477`: `ces` and the per-edge status array.
    let mut ces: Vec<usize> = Vec::new();
    let mut statuses = vec![0i32; nbedges + 1];

    // `cxx:480-574`: build the closed wires.
    for i in 1..=nbedges {
        if statuses[i] != 0 {
            continue;
        }
        ces.push(i);
        statuses[i] = 1;
        let mut search_backward = true;

        loop {
            if search_backward {
                search_backward = false;
                // `cxx:499-536`: look for a connection inside the chain.
                let last = *ces.last().unwrap();
                let edge = sewd.edge(last as isize);
                let mut found = false;
                let mut match_index = 0usize;
                if let Some(edge) = edge.as_ref() {
                    let lvertex = super::super::last_vertex(edge);
                    if let Some(lvertex) = lvertex.as_ref() {
                        let lpoint = BRepTool::vertex_point(lvertex);
                        let mut j = ces.len() as isize;
                        while j >= 1 && !found {
                            let fv = sewd
                                .edge(ces[j as usize - 1] as isize)
                                .and_then(|e| super::super::first_vertex(&e));
                            if let Some(fv) = fv.as_ref() {
                                let fp = BRepTool::vertex_point(fv);
                                if (shared && is_same(&lvertex.0, &fv.0))
                                    || (!shared && lpoint.distance(&fp) <= tolerance)
                                {
                                    found = true;
                                    match_index = j as usize;
                                }
                            }
                            j -= 1;
                        }
                    }
                }

                if found {
                    // `cxx:518-534`: `j` is the matching slot.
                    let start = match_index;
                    let mut wire1 = Wire::new();
                    for cesindex in start..=ces.len() {
                        builder.add_edge(&mut wire1, &sewd.edge(ces[cesindex - 1] as isize).unwrap());
                        statuses[ces[cesindex - 1]] = 2;
                    }
                    wire1.0.set_closed(true);
                    closed.push(wire1);
                    ces.truncate(start - 1);
                    if ces.is_empty() {
                        break;
                    }
                }
            }

            // `cxx:537-559`: look for a connection among the free edges.
            let mut found = false;
            let mut match_index = 0usize;
            let last = *ces.last().unwrap();
            let edge = sewd.edge(last as isize);
            if let Some(edge) = edge.as_ref() {
                if let Some(lvertex) = super::super::last_vertex(edge) {
                    let lpoint = BRepTool::vertex_point(&lvertex);
                    for j in 1..=nbedges {
                        if statuses[j] != 0 {
                            continue;
                        }
                        let fv = sewd.edge(j as isize).and_then(|e| super::super::first_vertex(&e));
                        if let Some(fv) = fv.as_ref() {
                            let fp = BRepTool::vertex_point(fv);
                            if (shared && is_same(&lvertex.0, &fv.0))
                                || (!shared && lpoint.distance(&fp) <= tolerance)
                            {
                                found = true;
                                match_index = j;
                                break;
                            }
                        }
                    }
                }
            }

            if found {
                // `cxx:554-562`.
                ces.push(match_index);
                statuses[match_index] = 1;
                search_backward = true;
                continue;
            }

            // `cxx:563-571`.
            let last = *ces.last().unwrap();
            statuses[last] = 3;
            ces.pop();
            if ces.is_empty() {
                break;
            }
        }
    }

    // `cxx:576-587`: build the open wires out of the remaining edges.
    let mut edges: Vec<Edge> = Vec::new();
    for i in 1..=nbedges {
        if statuses[i] != 2 {
            if let Some(e) = sewd.edge(i as isize) {
                edges.push(e);
            }
        }
    }
    let open_wires = connect_edges_to_wires(&mut edges, toler, shared);
    open.extend(open_wires);
}

/// `ShapeAnalysis_FreeBounds::SplitWires(wires, toler, shared, closed, open)`
/// (`cxx:589-606`).
pub(crate) fn split_wires_impl(
    wires: &[TopoShape],
    toler: f64,
    shared: bool,
) -> (Vec<TopoShape>, Vec<TopoShape>) {
    let mut closed: Vec<TopoShape> = Vec::new();
    let mut open: Vec<TopoShape> = Vec::new();
    for value in wires {
        let Some(wire) = Wire::wrap(value.clone()) else {
            continue;
        };
        let mut tmp_closed: Vec<Wire> = Vec::new();
        let mut tmp_open: Vec<Wire> = Vec::new();
        split_wire(&wire, toler, shared, &mut tmp_closed, &mut tmp_open);
        closed.extend(tmp_closed.into_iter().map(|w| w.0));
        open.extend(tmp_open.into_iter().map(|w| w.0));
    }
    (closed, open)
}
