//! `TopOpeBRepBuild_Builder1` — the per-shape state machinery.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1_1.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `PerformShapeWithStates()` | `183-360` |
//! | `PerformShapeWithStates(anObj, aReference)` | `363-464` |
//! | `PerformFacesWithStates(anObj, aFaces, aSplF)` | `467-552` |
//!
//! The state of every sub-shape of `anObj` against `aReference` ends up in
//! `myDataStructure->ChangeDS().ChangeMapOfShapeWithState(anObj)`:
//! - a shape the DS does not hold has no interference, so its state comes from
//!   one of its vertices (`Tools::FindStateThroughVertex`),
//! - a shape the DS holds is split, and the state is propagated from the split
//!   parts to the rest (`Tools::PropagateState`).
//!
//! OCCT keeps the two maps inside `TopOpeBRepDS_DataStructure` and holds the DS
//! through a handle, so it can mutate a map and read the DS in the same loop.
//! The port owns the DS, so each walk takes the map out of the DS, works on the
//! local copy and puts it back: the insertion order, and therefore the
//! `IndexedDataMap` indexing, is preserved exactly.
//!
//! `NCollection_DataMap::Bind` overwrites, while the port's
//! [`IndexedShapeData`] `add` keeps the first value, so the `Bind` sites go
//! through [`bind_state`].

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::edge_vertices;
use crate::topopebrep::build::tools::{
    find_state_through_vertex, map_shapes_of, propagate_state, propagate_state_for_wires,
};
use crate::topopebrep::ds::data_structure::{
    IndexedShapeData, IndexedShapeSet, ShapeWithStateMap,
};
use crate::topopebrep::ds::ShapeWithState;
use crate::Edge;

use super::super::Builder;
use super::{
    processed_edges_clear, unk_state_ver_clear, used_vertex_add, used_vertex_clear,
    THE_UNK_STATE_VER,
};

impl Builder {
    /// `PerformShapeWithStates()` (`cxx:183-360`): seeds the two maps with
    /// both arguments and with every rejected face as `OUT`, then walks each
    /// argument against the other.
    pub(super) fn perform_shape_with_states(&mut self) {
        // `cxx:185-187`
        used_vertex_clear();
        unk_state_ver_clear();

        let (mut obj, mut tool) = match self.my_data_structure.as_mut() {
            Some(hds) => {
                let ds = hds.change_ds();
                // `cxx:187-188`
                let obj = std::mem::replace(
                    ds.change_map_of_shape_with_state_obj(),
                    IndexedShapeData::new(),
                );
                let tool = std::mem::replace(
                    ds.change_map_of_shape_with_state_tool(),
                    IndexedShapeData::new(),
                );
                (obj, tool)
            }
            None => return,
        };

        // `cxx:190-204`: every new edge of every section curve.
        let nb_c = self
            .my_data_structure
            .as_ref()
            .map(|h| h.ds().nb_curves())
            .unwrap_or(0);
        for i in 1..=nb_c {
            let edges = self.change_new_edges(i).clone();
            for e in &edges {
                self.add_edge_vertices(e);
            }
        }

        // `cxx:206-252`: every split part of every section edge, `ON` first.
        let n = self
            .my_data_structure
            .as_ref()
            .map(|h| h.ds().nb_section_edges())
            .unwrap_or(0);
        for i in 1..=n {
            let e = match self.my_data_structure.as_ref() {
                Some(hds) => hds.ds().section_edge(i, true),
                None => None,
            };
            let Some(e) = e else {
                continue;
            };
            if e.is_null() {
                continue;
            }
            for st in [State::On, State::In, State::Out] {
                let parts = self.splits(&e, st).to_vec();
                for part in &parts {
                    self.add_edge_vertices(part);
                }
            }
        }

        // `cxx:256-266`: 1) add both arguments to facilitate the search.
        if let Some(s1) = self.my_shape1.clone() {
            obj.add(&s1, ShapeWithState::new());
        }
        if let Some(s2) = self.my_shape2.clone() {
            tool.add(&s2, ShapeWithState::new());
        }

        // `cxx:272-337`: 2) add all rejected shapes as `OUT`.
        let mut rejected_state = ShapeWithState::new();
        rejected_state.set_is_splitted(false);
        rejected_state.set_state(State::Out);

        let reject_obj = self.rejected_shapes(true);
        add_rejected_faces(&reject_obj, &rejected_state, &mut obj);
        let reject_tool = self.rejected_shapes(false);
        add_rejected_faces(&reject_tool, &rejected_state, &mut tool);

        // Put the maps back before the two walks, which look them up by
        // `Contains` (`cxx:369-374`).
        self.store_shape_with_state_maps(obj, tool);

        // `cxx:336-339`
        if let Some(s1) = self.my_shape1.clone() {
            if let Some(s2) = self.my_shape2.clone() {
                self.perform_shape_with_states_pair(&s1, &s2);
                processed_edges_clear();
                self.perform_shape_with_states_pair(&s2, &s1);
                processed_edges_clear();
            }
        }
        // `cxx:357-359`: the ON phase (`PerformOn2D`) is commented out in OCCT.
    }

    /// `PerformShapeWithStates(anObj, aReference)` (`cxx:363-464`).
    pub(super) fn perform_shape_with_states_pair(
        &mut self,
        an_obj: &TopoShape,
        a_reference: &TopoShape,
    ) {
        // `cxx:365`
        self
            .my_shape_classifier
            .borrow_mut()
            .set_reference(a_reference);

        // `cxx:366-375`
        let Some((which, mut map)) = self.take_shape_with_state_map(an_obj) else {
            return;
        };

        // `cxx:380-385`
        let shells = map_shapes_of(an_obj, ShapeType::Shell);
        for i in 1..=shells.extent() {
            let Some(shell) = shells.find_key(i).cloned() else {
                continue;
            };

            // `cxx:392-394`
            if map.contains(&shell) {
                continue;
            } else if !self.has_shape(&shell) {
                // `cxx:395-404`: A.1, the shell has no interference.
                let avoid = THE_UNK_STATE_VER.with(|m| m.borrow().clone());
                find_state_through_vertex(
                    &shell,
                    &mut self.my_shape_classifier.borrow_mut(),
                    &mut map,
                    &avoid,
                );
                continue;
            }

            // `cxx:406-462`: A.2, split it into faces to rest and interfered
            // faces.
            let mut faces_with_interf = IndexedShapeSet::new();
            let mut faces_to_rest = IndexedShapeSet::new();
            let mut spl_faces_state = IndexedShapeData::<State>::new();

            let faces = map_shapes_of(&shell, ShapeType::Face);
            for j in 1..=faces.extent() {
                let Some(face) = faces.find_key(j).cloned() else {
                    continue;
                };
                if map.contains(&face) {
                    // `cxx:420-438`: a known face means its edges are known
                    // too; carry their states into `aSplFacesState`.
                    let edges = map_shapes_of(&face, ShapeType::Edge);
                    for k in 1..=edges.extent() {
                        let Some(edge) = edges.find_key(k) else {
                            continue;
                        };
                        if let Some(sws) = map.find_from_key(edge) {
                            let st = sws.state();
                            bind_state(&mut spl_faces_state, edge, st);
                        }
                    }
                    continue;
                } else if self.has_shape(&face) {
                    faces_with_interf.add(&face);
                } else {
                    faces_to_rest.add(&face);
                }
            }

            // `cxx:448`
            self.perform_faces_with_states(
                an_obj,
                &faces_with_interf,
                &mut spl_faces_state,
                &mut map,
            );

            // `cxx:451-461`: OCCT declares a *fresh empty* `anUnkStateEdge`
            // here, so nothing is avoided; `theUnkStateVer` is only consulted
            // by the `StatusEdgesToSplit` call site.
            let avoid = IndexedShapeSet::new();
            propagate_state(
                &spl_faces_state,
                &faces_to_rest,
                ShapeType::Edge,
                ShapeType::Face,
                &mut self.my_shape_classifier.borrow_mut(),
                &mut map,
                &avoid,
            );
            propagate_state_for_wires(&faces_to_rest, &mut map);
        }

        self.store_shape_with_state_map(which, map);
    }

    /// `PerformFacesWithStates(anObj, aFaces, aSplF)` (`cxx:467-552`).
    ///
    /// The map is threaded in from `PerformShapeWithStates(anObj, aReference)`,
    /// which is the OCCT reference `aMapOfShapeWithState`; OCCT re-fetches it
    /// with `ChangeMapOfShapeWithState(anObj, aFlag)` (`cxx:472-480`) but the
    /// `aFlag` gate has already been passed there.
    pub(super) fn perform_faces_with_states(
        &mut self,
        an_obj: &TopoShape,
        a_faces_with_interferences_map: &IndexedShapeSet,
        a_spl_faces_state: &mut IndexedShapeData<State>,
        map: &mut IndexedShapeData<ShapeWithState>,
    ) {
        // `cxx:485-551`
        let n_f = a_faces_with_interferences_map.extent();
        for i in 1..=n_f {
            let Some(face) = a_faces_with_interferences_map.find_key(i).cloned() else {
                continue;
            };

            let wires = map_shapes_of(&face, ShapeType::Wire);
            for j in 1..=wires.extent() {
                let Some(wire) = wires.find_key(j).cloned() else {
                    continue;
                };

                if !self.has_shape(&wire) {
                    // `cxx:501-509`
                    let avoid = THE_UNK_STATE_VER.with(|m| m.borrow().clone());
                    find_state_through_vertex(
                        &wire,
                        &mut self.my_shape_classifier.borrow_mut(),
                        map,
                        &avoid,
                    );
                    continue;
                }

                // `cxx:511-549`: the wire has interferences.
                let mut edges_to_split = IndexedShapeSet::new();
                let mut edges_to_rest = IndexedShapeSet::new();
                let edges = map_shapes_of(&wire, ShapeType::Edge);
                for k in 1..=edges.extent() {
                    let Some(edge) = edges.find_key(k).cloned() else {
                        continue;
                    };
                    if self.has_shape(&edge) {
                        edges_to_split.add(&edge);
                    } else {
                        edges_to_rest.add(&edge);
                    }
                }

                // `cxx:534`
                self.status_edges_to_split(an_obj, &edges_to_split, &edges_to_rest, map);

                // `cxx:538-548`: the states of the rest edges are known now.
                for k in 1..=edges_to_rest.extent() {
                    let Some(edge) = edges_to_rest.find_key(k) else {
                        continue;
                    };
                    if let Some(sws) = map.find_from_key(edge) {
                        let st = sws.state();
                        bind_state(a_spl_faces_state, edge, st);
                    }
                }
            }
        }
    }

    /// `TopExp::Vertices(TopoDS::Edge(E), Vf, Vl)` followed by two
    /// `theUsedVertexMap.Add` calls (`cxx:199-203` and the other sites).
    ///
    /// OCCT binds the two vertex handles even when they are null; the port's
    /// map key cannot stand for a null shape, so a null vertex is skipped. The
    /// map is only ever queried with a real vertex (`theUsedVertexMap
    /// .Contains(V)` at `cxx:932`), so the answer is the same.
    fn add_edge_vertices(&self, e: &TopoShape) {
        let (vf, vl) = edge_vertices(&Edge(e.clone()));
        if let Some(vf) = vf {
            used_vertex_add(&vf);
        }
        if let Some(vl) = vl {
            used_vertex_add(&vl);
        }
    }

    /// `myDataStructure->HasShape(S)` (default `FindKeep = true`).
    fn has_shape(&self, s: &TopoShape) -> bool {
        match self.my_data_structure.as_ref() {
            Some(hds) => hds.has_shape(s, true),
            None => false,
        }
    }

    /// The rejected-shape map of one operand (`cxx:269-272`): `true` is
    /// `ChangeMapOfRejectedShapesObj()`, `false` the tool map.
    fn rejected_shapes(&mut self, obj: bool) -> IndexedShapeSet {
        let Some(hds) = self.my_data_structure.as_mut() else {
            return IndexedShapeSet::new();
        };
        let ds = hds.change_ds();
        if obj {
            ds.change_map_of_rejected_shapes_obj().clone()
        } else {
            ds.change_map_of_rejected_shapes_tool().clone()
        }
    }

    /// `aDS.ChangeMapOfShapeWithState(anObj, aFlag)` (`cxx:367-375`): the map
    /// is taken out of the DS so the walk can also read the DS. `None` is
    /// `aFlag == false`.
    pub(super) fn take_shape_with_state_map(
        &mut self,
        an_obj: &TopoShape,
    ) -> Option<(ShapeWithStateMap, IndexedShapeData<ShapeWithState>)> {
        let hds = self.my_data_structure.as_mut()?;
        let ds = hds.change_ds();
        let (view, which) = ds.change_map_of_shape_with_state(an_obj)?;
        Some((which, std::mem::replace(view, IndexedShapeData::new())))
    }

    /// Put one map back where `take_shape_with_state_map` found it.
    pub(super) fn store_shape_with_state_map(
        &mut self,
        which: ShapeWithStateMap,
        map: IndexedShapeData<ShapeWithState>,
    ) {
        let Some(hds) = self.my_data_structure.as_mut() else {
            return;
        };
        let ds = hds.change_ds();
        match which {
            ShapeWithStateMap::Obj => *ds.change_map_of_shape_with_state_obj() = map,
            ShapeWithStateMap::Tool => *ds.change_map_of_shape_with_state_tool() = map,
        }
    }

    /// Put both maps back (`cxx:184-188` emptied them).
    fn store_shape_with_state_maps(
        &mut self,
        obj: IndexedShapeData<ShapeWithState>,
        tool: IndexedShapeData<ShapeWithState>,
    ) {
        let Some(hds) = self.my_data_structure.as_mut() else {
            return;
        };
        let ds = hds.change_ds();
        *ds.change_map_of_shape_with_state_obj() = obj;
        *ds.change_map_of_shape_with_state_tool() = tool;
    }
}

/// `NCollection_DataMap::Bind`: overwrite an existing binding, append a new
/// one at the end (so the key order of an existing key is kept).
pub(super) fn bind_state(map: &mut IndexedShapeData<State>, key: &TopoShape, st: State) {
    match map.change_from_key(key) {
        Some(v) => *v = st,
        None => {
            map.add(key, st);
        }
    }
}

/// `cxx:280-337`: every edge, wire and face of every rejected face gets the
/// `OUT` state.
fn add_rejected_faces(
    rejected: &IndexedShapeSet,
    state: &ShapeWithState,
    map: &mut IndexedShapeData<ShapeWithState>,
) {
    for i in 1..=rejected.extent() {
        let Some(a_face) = rejected.find_key(i).cloned() else {
            continue;
        };
        if a_face.shape_type() != ShapeType::Face {
            continue;
        }
        let wires = map_shapes_of(&a_face, ShapeType::Wire);
        for iw in 1..=wires.extent() {
            let Some(a_wire) = wires.find_key(iw).cloned() else {
                continue;
            };
            let edges = map_shapes_of(&a_wire, ShapeType::Edge);
            for j in 1..=edges.extent() {
                if let Some(an_edge) = edges.find_key(j) {
                    map.add(an_edge, state.clone());
                }
            }
            map.add(&a_wire, state.clone());
        }
        map.add(&a_face, state.clone());
    }
}
