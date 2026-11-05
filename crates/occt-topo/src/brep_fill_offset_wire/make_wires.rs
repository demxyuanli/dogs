//! `BRepFill_OffsetWire::MakeWires` (`cxx:1405-1579`): chains the parallel
//! edges into wires by walking the vertex-to-edge map, then stores the result in
//! `myShape` (a wire, or a compound when several wires come out).

use crate::brep_fill_offset_wire::vertices::{check_small_param_on_edge, is_small_closed_edge};
use crate::brep_fill_offset_wire::{BRepFillOffsetWire, IndexedShapeMap};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same};

impl BRepFillOffsetWire {
    /// `MakeWires()` (`cxx:1405-1579`).
    pub(crate) fn make_wires(&mut self) {
        // `cxx:1413`: `MVE` — vertex to the list of edges incident to it.
        let mut mve = IndexedShapeMap::new();

        // `cxx:1422-1448`.
        for i in 1..=self.my_map.extent() {
            let edges = self.my_map.value(i).cloned().unwrap_or_default();
            for item in edges {
                let e = Edge(item);
                let (v1, v2) = edge_vertices(&e);
                // `cxx:1422-1424`: `TopExp::Vertices` leaves two null vertices
                // for an edge built without them, `V1.IsSame(V2)` holds for
                // those, and `IsSmallClosedEdge` then raises
                // `Standard_NullObject` on its first statement,
                // `BRep_Tool::Pnt(aVertex)` (`cxx:2633`).
                if v1.is_none() && v2.is_none() {
                    panic!("BRep_Tool:: TopoDS_Vertex hasn't gp_Pnt");
                }
                let (Some(v1), Some(v2)) = (v1, v2) else {
                    continue;
                };
                if is_same(&v1.0, &v2.0) && is_small_closed_edge(&e, &v1) {
                    continue;
                }
                if !check_small_param_on_edge(&e) {
                    continue;
                }
                if !mve.contains(&v1.0) {
                    mve.add(&v1.0, Vec::new());
                }
                if let Some(list) = mve.change_from_key(&v1.0) {
                    list.push(e.0.clone());
                }
                if !mve.contains(&v2.0) {
                    mve.add(&v2.0, Vec::new());
                }
                if let Some(list) = mve.change_from_key(&v2.0) {
                    list.push(e.0.clone());
                }
            }
        }

        let builder = TopoBuilder::new();
        let mut the_wires: Vec<TopoShape> = Vec::new();

        // `cxx:1466-1558`.
        while !mve.is_empty() {
            let mut wire_edges: Vec<Edge> = Vec::new();

            // `cxx:1470-1481`: start from a vertex with a single edge if any.
            let mut start = 0usize;
            for i in 1..=mve.extent() {
                if mve.value(i).map_or(0, |l| l.len()) == 1 {
                    start = i;
                    break;
                }
            }
            if start == 0 {
                start = 1;
            }

            let mut cv = Vertex(
                mve.find_key(start)
                    .expect("MakeWires: empty map slot")
                    .clone(),
            );
            let vf = cv.clone();
            let mut ce = Edge(
                mve.value(start)
                    .and_then(|l| l.first())
                    .expect("MakeWires: empty incident list")
                    .clone(),
            );
            let mut end = false;
            if let Some(list) = mve.change_from_key(&cv.0) {
                if !list.is_empty() {
                    list.remove(0);
                }
            }
            if self.my_is_open_result && mve.find_from_key(&cv.0).map_or(true, |l| l.is_empty()) {
                mve.remove_key(&cv.0);
            }

            let mut is_closed = false;
            while !end {
                // `cxx:1505-1519`.
                let (v1, v2) = edge_vertices(&ce);
                cv = match (v1, v2) {
                    (Some(a), Some(b)) => {
                        if !is_same(&cv.0, &a.0) {
                            a
                        } else {
                            b
                        }
                    }
                    _ => cv.clone(),
                };

                wire_edges.push(Edge(ce.0.clone()));

                if is_same(&vf.0, &cv.0)
                    || mve.find_from_key(&cv.0).map_or(true, |l| l.is_empty())
                {
                    is_closed = is_same(&vf.0, &cv.0);
                    end = true;
                    mve.remove_key(&vf.0);
                }

                if !end {
                    // `cxx:1529-1535`: drop the edge just walked.
                    if let Some(list) = mve.change_from_key(&cv.0) {
                        if let Some(pos) = list.iter().position(|x| is_same(x, &ce.0)) {
                            list.remove(pos);
                        }
                    }
                    let next = mve
                        .find_from_key(&cv.0)
                        .and_then(|l| l.first())
                        .cloned();
                    if let Some(next) = next {
                        ce = Edge(next);
                        if let Some(list) = mve.change_from_key(&cv.0) {
                            if !list.is_empty() {
                                list.remove(0);
                            }
                        }
                    } else if self.my_is_open_result {
                        end = true;
                    }

                    if mve.find_from_key(&cv.0).map_or(true, |l| l.is_empty()) {
                        mve.remove_key(&cv.0);
                    }
                }
            }

            let wire = builder.make_wire(&wire_edges);
            wire.0.set_closed(is_closed);
            the_wires.push(wire.0);
        }

        // `cxx:1560-1576`.
        if the_wires.len() == 1 {
            self.my_shape = Some(the_wires.remove(0));
        } else {
            let compound = builder.make_compound_of(&the_wires);
            self.my_shape = Some(compound.0);
        }
    }
}
