//! `LocOpe_BuildWires` (`LocOpe_BuildWires.hxx:25-52`, `.cxx:44-279`).
//!
//! Chains a loose list of edges into wires along the boundary of the shape that
//! `LocOpe_WiresOnShape` describes. A vertex counts as a border vertex when it
//! has a single incident edge in the list, or when `PW` reports it on the
//! original shape (`LocOpe_BuildWires.cxx:94-110`). Starting from the first
//! remaining edge at a border vertex, the walk follows unused incident edges
//! until it either comes back to a visited vertex (closed wire) or reaches
//! another border vertex (open wire).
//!
//! Statement mapping:
//! - `theMapVE` (`NCollection_IndexedDataMap<TopoDS_Shape,
//!   NCollection_List<TopoDS_Shape>>`) is [`AncestorMap`]; the loop that drops
//!   consumed edges from a vertex's list (`cxx:206-222`) uses
//!   [`AncestorMap::change_ancestors`].
//! - `mapE` (`NCollection_IndexedMap`) keeps edge occurrences with their walk
//!   orientation; membership (`mapE.Contains`, `cxx:270` / `:293`) is the
//!   orientation-insensitive `TopTools_ShapeMapHasher` test, so a [`ShapeSet`]
//!   stands in for it next to the ordered vector.
//! - `mapV` / `Bords` (`NCollection_Map`) are [`ShapeSet`]s.

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeSet;
use crate::locope::wires_on_shape::WiresOnShape;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{self, AncestorMap};

/// `LocOpe_BuildWires` (`hxx:28-52`).
#[derive(Debug, Clone, Default)]
pub struct BuildWires {
    my_done: bool,
    my_res: Vec<TopoShape>,
}

impl BuildWires {
    /// `LocOpe_BuildWires()` (`cxx:44-46`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `LocOpe_BuildWires(L, PW)` (`cxx:52-55`): the workload constructor.
    pub fn with(&mut self, ledges: &[TopoShape], pw: &WiresOnShape) {
        let _ = self.perform(ledges, pw);
    }

    /// `Perform(L, PW)` (`cxx:63-229`).
    ///
    /// The reference raises `Standard_ConstructionError` when the walk reaches a
    /// vertex whose incident edges are all used (`cxx:180-183`); the port
    /// reports it as the `Err` string and leaves `my_done` `false`, which is
    /// what `IsDone()` / `Result()` read.
    pub fn perform(&mut self, ledges: &[TopoShape], pw: &WiresOnShape) -> Result<(), String> {
        self.my_done = false;
        self.my_res.clear();

        let builder = TopoBuilder::new();

        // `cxx:72-84`: compound of the distinct edges, all FORWARD (the comment
        // at `cxx:79` notes the orientation matters for `TopExp::Vertices`).
        let mut the_map = ShapeSet::new();
        let mut edges: Vec<TopoShape> = Vec::new();
        for edg in ledges {
            if the_map.contains(edg) {
                continue;
            }
            the_map.add(edg);
            if edg.shape_type() == ShapeType::Edge {
                edges.push(edg.oriented(Orientation::Forward));
            }
        }
        let compound = builder.make_compound_of(&edges);

        // `cxx:90`: `TopExp::MapShapesAndAncestors(C, VERTEX, EDGE, theMapVE)`.
        let mut the_map_ve =
            AncestorMap::map_shapes_and_ancestors(&compound.0, ShapeType::Vertex, ShapeType::Edge);

        // `cxx:94-110`: border vertices.
        let mut bords = ShapeSet::new();
        for i in 0..the_map_ve.length() {
            let vtx = the_map_ve.find_key(i).clone();
            if the_map_ve.ancestors(i).len() == 1 {
                bords.add(&vtx);
                continue;
            }
            let Some(v) = Vertex::wrap(vtx.clone()) else {
                continue;
            };
            if pw.on_vertex(&v).is_some() || pw.on_edge_of_vertex(&v).is_some() {
                bords.add(&vtx);
            }
        }

        // `cxx:112-228`: one wire per walk.
        let extent = the_map_ve.length();
        loop {
            let i = find_first_edge(&the_map_ve, &bords, extent);
            if i >= extent {
                break;
            }

            // `cxx:114-116`.
            let mut map_e: Vec<TopoShape> = Vec::new();
            let mut map_e_keys = ShapeSet::new();
            let mut map_v = ShapeSet::new();

            let first = the_map_ve.ancestors(i)[0].clone();
            let Some(edgf) = Edge::wrap(first) else {
                return Err("LocOpe_BuildWires: ancestor is not an edge".to_string());
            };
            // `TopExp::Vertices(edgf, VF, VL)` (`cxx:121`), `CumOri = false`.
            let (Some(vf), Some(vl)) = topo_tools::edge_vertices(&edgf) else {
                return Err("LocOpe_BuildWires: edge without two vertices".to_string());
            };
            let mut vf = vf;
            let mut vl = vl;

            // `cxx:123-133`.
            if bords.contains(&vl.0) && !bords.contains(&vf.0) {
                push_edge(&mut map_e, &mut map_e_keys, edgf.0.oriented(Orientation::Reversed));
                std::mem::swap(&mut vf, &mut vl);
            } else {
                push_edge(&mut map_e, &mut map_e_keys, edgf.0.oriented(Orientation::Forward));
            }
            map_v.add(&vf.0);

            // `cxx:136-158`: follow unused incident edges.
            while !(map_v.contains(&vl.0) || bords.contains(&vl.0)) {
                let Some(ind) = the_map_ve.find_index(&vl.0) else {
                    return Err("LocOpe_BuildWires: vertex is not in the map".to_string());
                };
                // `cxx:139-147`: first incident edge not already in the walk.
                let mut next: Option<TopoShape> = None;
                for anc in the_map_ve.ancestors(ind) {
                    if !map_e_keys.contains(anc) {
                        next = Some(anc.clone());
                        break;
                    }
                }
                let Some(next) = next else {
                    // `cxx:180-183`: `throw Standard_ConstructionError()`.
                    return Err("LocOpe_BuildWires: no chaining edge".to_string());
                };
                let Some(the_edge) = Edge::wrap(next) else {
                    return Err("LocOpe_BuildWires: ancestor is not an edge".to_string());
                };
                // `cxx:150`: `TopExp::Vertices(theEdge, Vf, Vl)`.
                let (Some(edge_f), Some(edge_l)) = topo_tools::edge_vertices(&the_edge) else {
                    return Err("LocOpe_BuildWires: edge without two vertices".to_string());
                };
                map_v.add(&vl.0);
                // `cxx:151-157`.
                if topo_tools::is_same(&edge_f.0, &vl.0) {
                    push_edge(
                        &mut map_e,
                        &mut map_e_keys,
                        the_edge.0.oriented(Orientation::Forward),
                    );
                    vl = edge_l;
                } else {
                    push_edge(
                        &mut map_e,
                        &mut map_e_keys,
                        the_edge.0.oriented(Orientation::Reversed),
                    );
                    vl = edge_f;
                }
            }

            // `cxx:160-201`: turn `mapE` into a wire.
            let mut wire = Wire::new();
            if map_v.contains(&vl.0) {
                // `cxx:164-186`: came back to a visited vertex - trim the tail
                // so the wire starts and ends at `VL`.
                let mut j = 0usize;
                while j < map_e.len() {
                    let Some(edg) = Edge::wrap(map_e[j].clone()) else {
                        return Err("LocOpe_BuildWires: walk holds a non-edge".to_string());
                    };
                    // `cxx:169-175`: `TopExp::FirstVertex` / `LastVertex` with
                    // `CumOri = false`, chosen by the walk orientation.
                    let vf = if map_e[j].orientation() == Orientation::Forward {
                        topo_tools::edge_vertices(&edg).0
                    } else {
                        topo_tools::edge_vertices(&edg).1
                    };
                    let Some(vf) = vf else {
                        return Err("LocOpe_BuildWires: edge without first vertex".to_string());
                    };
                    if topo_tools::is_same(&vf.0, &vl.0) {
                        break;
                    }
                    map_v.remove(&vf.0);
                    j += 1;
                }
                let kept: Vec<Edge> = map_e[j..]
                    .iter()
                    .filter_map(|s| Edge::wrap(s.clone()))
                    .collect();
                for e in &kept {
                    builder.add_edge(&mut wire, e);
                }
                // `cxx:186`: `newWire.Closed(true)`.
                wire.0.set_closed(true);
            } else {
                // `cxx:187-194`: left the chain through a border vertex.
                map_v.add(&vl.0);
                let kept: Vec<Edge> = map_e
                    .iter()
                    .filter_map(|s| Edge::wrap(s.clone()))
                    .collect();
                for e in &kept {
                    builder.add_edge(&mut wire, e);
                }
                // `cxx:194`: `newWire.Closed(false)`.
                wire.0.set_closed(false);
            }
            self.my_res.push(wire.0);

            // `cxx:204-222`: the walk's vertices become borders and its edges
            // leave the incident lists.
            for vtx in map_v.iter() {
                bords.add(vtx);
                let Some(ind) = the_map_ve.find_index(vtx) else {
                    continue;
                };
                the_map_ve
                    .change_ancestors(ind)
                    .retain(|anc| !map_e_keys.contains(anc));
            }
        }

        self.my_done = true;
        Ok(())
    }

    /// `IsDone()` (`cxx:232-235`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Result()` (`cxx:239-246`); `None` stands for the `StdFail_NotDone` raise.
    pub fn result(&self) -> Option<&[TopoShape]> {
        if self.my_done {
            Some(&self.my_res)
        } else {
            None
        }
    }
}

/// `mapE.Add(edge)` (`NCollection_IndexedMap::Add`): append unless the shape is
/// already present.
fn push_edge(map_e: &mut Vec<TopoShape>, keys: &mut ShapeSet, edge: TopoShape) {
    if keys.contains(&edge) {
        return;
    }
    keys.add(&edge);
    map_e.push(edge);
}

/// `FindFirstEdge(theMapVE, Bords)` (`cxx:250-279`), answered as a 0-based index;
/// `extent` is OCCT's `Extent() + 1`.
fn find_first_edge(map: &AncestorMap, bord: &ShapeSet, extent: usize) -> usize {
    // `cxx:256-260`: first vertex with any incident edge.
    let mut i = 0usize;
    while i < extent && map.ancestors(i).is_empty() {
        i += 1;
    }
    if i >= extent {
        return extent;
    }

    // `cxx:262-276`: prefer one of those that is a border vertex.
    let goodi = i;
    while i < extent {
        if !map.ancestors(i).is_empty() && bord.contains(map.find_key(i)) {
            return i;
        }
        i += 1;
    }
    goodi
}
