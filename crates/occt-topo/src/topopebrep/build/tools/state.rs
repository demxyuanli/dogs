//! `TopOpeBRepBuild_Tools` — the state propagation family.
//!
//! Sources: `TopOpeBRepBuild_Tools.cxx:60-92` (`FindState`), `:94-252`
//! (`PropagateState`), `:254-290` (`FindState2`), `:292-333` (`FindState1`),
//! `:335-387` (`FindStateThroughVertex`), `:389-410` (`SpreadStateToChild`),
//! `:412-457` (`PropagateStateForWires`) — OCCT 8.0.0.
//!
//! Map conventions of the port (`ds/data_structure/index_map.rs`):
//! `NCollection_DataMap` / `IndexedDataMap` -> [`IndexedShapeData`] (insertion
//! order kept, keys are `IsSame`), `NCollection_Map` / `IndexedMap` ->
//! [`IndexedShapeSet`], `NCollection_List<TopoDS_Shape>` -> `Vec<TopoShape>`.

use std::collections::HashMap;

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors_into, shape_map_key, ShapeMapKey};
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};
use crate::topopebrep::ds::shape_with_state::ShapeWithState;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;

/// `TopExp::MapShapes(S, T, M)` (`TopExp.cxx:32-44`): one `M.Add` per
/// `TopExp_Explorer(S, T)` hit, so first-occurrence order is kept and repeats
/// are dropped through `IsSame`.
pub fn map_shapes_of(shape: &TopoShape, t: ShapeType) -> IndexedShapeSet {
    let mut m = IndexedShapeSet::new();
    for s in explorer_shapes(shape, t) {
        m.add(&s);
    }
    m
}

/// `NCollection_DataMap::Bind(S, V)`: a fresh key appends, an existing key
/// overwrites in place (the index of the key never changes).
fn bind(map: &mut IndexedShapeData<State>, s: &TopoShape, state: State) {
    match map.change_from_key(s) {
        Some(v) => *v = state,
        None => {
            map.add(s, state);
        }
    }
}

/// `aMapSubshAnc.FindFromKey(aSubsh)`, i.e. the ancestor list of one sub-shape.
fn ancestors_of<'a>(
    map: &'a HashMap<ShapeMapKey, Vec<TopoShape>>,
    s: &TopoShape,
) -> &'a [TopoShape] {
    map.get(&shape_map_key(s))
        .map(|v| v.as_slice())
        .unwrap_or(&[])
}

/// `FindState(aSubsh, aState, aSubshEnum, aMapSubshAnc, aMapProcessedSubsh,
/// aMapSS)` (`cxx:60-92`).
pub fn find_state(
    a_subsh: &TopoShape,
    a_state: State,
    a_subsh_enum: ShapeType,
    a_map_subsh_anc: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    a_map_processed_subsh: &mut IndexedShapeSet,
    a_map_ss: &mut IndexedShapeData<State>,
) {
    let list: Vec<TopoShape> = ancestors_of(a_map_subsh_anc, a_subsh).to_vec();
    for a_s in &list {
        let subsh_map = map_shapes_of(a_s, a_subsh_enum);
        for i in 1..=subsh_map.extent() {
            let Some(a_ss) = subsh_map.find_key(i).cloned() else {
                continue;
            };
            if !a_map_processed_subsh.contains(&a_ss) {
                a_map_processed_subsh.add(&a_ss);
                bind(a_map_ss, &a_ss, a_state);
                find_state(
                    &a_ss,
                    a_state,
                    a_subsh_enum,
                    a_map_subsh_anc,
                    a_map_processed_subsh,
                    a_map_ss,
                );
            }
        }
    }
}

/// `FindState2(aSubsh, aState, aMapSubshAnc, aMapProcessedSubsh, aMapSS)`
/// (`cxx:254-290`) — the `aSubshEnum == TopAbs_VERTEX` arm of
/// `PropagateState`: the ancestor shape itself is bound, then its vertices.
pub fn find_state2(
    a_subsh: &TopoShape,
    a_state: State,
    a_map_subsh_anc: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    a_map_processed_subsh: &mut IndexedShapeSet,
    a_map_ss: &mut IndexedShapeData<State>,
) {
    let list: Vec<TopoShape> = ancestors_of(a_map_subsh_anc, a_subsh).to_vec();
    for a_s in &list {
        bind(a_map_ss, a_s, a_state);
        let subsh_map = map_shapes_of(a_s, ShapeType::Vertex);
        for i in 1..=subsh_map.extent() {
            let Some(a_ss) = subsh_map.find_key(i).cloned() else {
                continue;
            };
            if !a_map_processed_subsh.contains(&a_ss) {
                a_map_processed_subsh.add(&a_ss);
                bind(a_map_ss, &a_ss, a_state);
                find_state2(
                    &a_ss,
                    a_state,
                    a_map_subsh_anc,
                    a_map_processed_subsh,
                    a_map_ss,
                );
            }
        }
    }
}

/// `FindState1(aSubsh, aState, aMapSubshAnc, aMapProcessedSubsh, aMapSS)`
/// (`cxx:292-333`) — the `aSubshEnum == TopAbs_EDGE` arm of `PropagateState`:
/// the ancestor face, its wires, then its edges are bound.
pub fn find_state1(
    a_subsh: &TopoShape,
    a_state: State,
    a_map_subsh_anc: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    a_map_processed_subsh: &mut IndexedShapeSet,
    a_map_ss: &mut IndexedShapeData<State>,
) {
    let list: Vec<TopoShape> = ancestors_of(a_map_subsh_anc, a_subsh).to_vec();
    for a_s in &list {
        bind(a_map_ss, a_s, a_state);
        let wire_map = map_shapes_of(a_s, ShapeType::Wire);
        for j in 1..=wire_map.extent() {
            if let Some(w) = wire_map.find_key(j).cloned() {
                bind(a_map_ss, &w, a_state);
            }
        }
        let subsh_map = map_shapes_of(a_s, ShapeType::Edge);
        for i in 1..=subsh_map.extent() {
            let Some(a_ss) = subsh_map.find_key(i).cloned() else {
                continue;
            };
            if !a_map_processed_subsh.contains(&a_ss) {
                a_map_processed_subsh.add(&a_ss);
                bind(a_map_ss, &a_ss, a_state);
                find_state1(
                    &a_ss,
                    a_state,
                    a_map_subsh_anc,
                    a_map_processed_subsh,
                    a_map_ss,
                );
            }
        }
    }
}

/// `SpreadStateToChild(aShape, aState, aMapOfShapeWithState)` (`cxx:389-410`):
/// binds every face, wire and edge of `aShape` to `aState`.
pub fn spread_state_to_child(
    a_shape: &TopoShape,
    a_state: State,
    a_map_of_shape_with_state: &mut IndexedShapeData<ShapeWithState>,
) {
    let mut child_map = IndexedShapeSet::new();
    for t in [ShapeType::Face, ShapeType::Wire, ShapeType::Edge] {
        for s in explorer_shapes(a_shape, t) {
            child_map.add(&s);
        }
    }

    let mut a_shape_with_state = ShapeWithState::new();
    a_shape_with_state.set_state(a_state);
    a_shape_with_state.set_is_splitted(false);

    for i in 1..=child_map.extent() {
        if let Some(c) = child_map.find_key(i) {
            a_map_of_shape_with_state.add(c, a_shape_with_state.clone());
        }
    }
}

/// `FindStateThroughVertex(aShape, aShapeClassifier, aMapOfShapeWithState,
/// anAvoidSubshMap)` (`cxx:335-387`). Returns `TopAbs_UNKNOWN` when every
/// vertex *and* every edge of `aShape` is avoided (`cxx:369-374`).
pub fn find_state_through_vertex(
    a_shape: &TopoShape,
    a_shape_classifier: &mut ShapeClassifier,
    a_map_of_shape_with_state: &mut IndexedShapeData<ShapeWithState>,
    an_avoid_subsh_map: &IndexedShapeSet,
) -> State {
    let mut a_subsh = first_avoided_by(a_shape, ShapeType::Vertex, an_avoid_subsh_map);
    if a_subsh.is_none() {
        a_subsh = first_avoided_by(a_shape, ShapeType::Edge, an_avoid_subsh_map);
        if a_subsh.is_none() {
            return State::Unknown;
        }
    }
    let a_subsh = a_subsh.expect("checked above");

    // `aShapeClassifier.StateShapeReference(aSubsh, aNullShape)` (`cxx:378`).
    let a_state = a_shape_classifier
        .state_shape_reference_null_avoid(&a_subsh)
        .unwrap_or(State::Unknown);

    let mut a_shape_with_state = ShapeWithState::new();
    a_shape_with_state.set_state(a_state);
    a_shape_with_state.set_is_splitted(false);
    a_map_of_shape_with_state.add(a_shape, a_shape_with_state);
    spread_state_to_child(a_shape, a_state, a_map_of_shape_with_state);
    a_state
}

/// The first `t` sub-shape of `shape` that `avoid` does not contain
/// (`cxx:346-357`, `:359-367`).
fn first_avoided_by(
    shape: &TopoShape,
    t: ShapeType,
    avoid: &IndexedShapeSet,
) -> Option<TopoShape> {
    let subsh_map = map_shapes_of(shape, t);
    for i in 1..=subsh_map.extent() {
        let s = subsh_map.find_key(i)?;
        if !avoid.contains(s) {
            return Some(s.clone());
        }
    }
    None
}

/// `PropagateState(aSplShapesState, aShapesToRestMap, aSubshEnum, aShapeEnum,
/// aShapeClassifier, aMapOfShapeWithState, anAvoidSubshMap)` (`cxx:94-252`).
pub fn propagate_state(
    a_spl_shapes_state: &IndexedShapeData<State>,
    a_shapes_to_rest_map: &IndexedShapeSet,
    a_subsh_enum: ShapeType,
    a_shape_enum: ShapeType,
    a_shape_classifier: &mut ShapeClassifier,
    a_map_of_shape_with_state: &mut IndexedShapeData<ShapeWithState>,
    an_avoid_subsh_map: &IndexedShapeSet,
) {
    // 1. states of the sub-shapes of every split shape (`cxx:114-123`).
    let mut a_map_ss: IndexedShapeData<State> = IndexedShapeData::new();
    for k in 1..=a_spl_shapes_state.extent() {
        let (Some(a_shape), Some(&a_state)) = (
            a_spl_shapes_state.find_key(k),
            a_spl_shapes_state.find_from_index(k),
        ) else {
            continue;
        };
        let a_subshapes = map_shapes_of(a_shape, a_subsh_enum);
        for j in 1..=a_subshapes.extent() {
            let Some(sub) = a_subshapes.find_key(j).cloned() else {
                continue;
            };
            if !an_avoid_subsh_map.contains(&sub) {
                bind(&mut a_map_ss, &sub, a_state);
            }
        }
    }

    let a_map_ss1 = a_map_ss.clone();

    // 2. shapes and ancestors of the shapes to rest (`cxx:127-134`).
    let mut a_map_subsh_anc: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    for j in 1..=a_shapes_to_rest_map.extent() {
        if let Some(s) = a_shapes_to_rest_map.find_key(j) {
            map_shapes_and_ancestors_into(s, a_subsh_enum, a_shape_enum, &mut a_map_subsh_anc);
        }
    }

    // 3. propagate over the connexity graph (`cxx:137-148`).
    let mut a_processed_subshapes = IndexedShapeSet::new();
    for k in 1..=a_map_ss1.extent() {
        let (Some(a_subsh), Some(&a_state)) =
            (a_map_ss1.find_key(k), a_map_ss1.find_from_index(k))
        else {
            continue;
        };
        if a_map_subsh_anc.contains_key(&shape_map_key(a_subsh)) {
            a_processed_subshapes.add(a_subsh);
            find_state(
                a_subsh,
                a_state,
                a_subsh_enum,
                &a_map_subsh_anc,
                &mut a_processed_subshapes,
                &mut a_map_ss,
            );
        }
    }

    // 4. bind the resting shapes whose first sub-shape is already known
    //    (`cxx:151-180`).
    let mut a_non_passed_shapes = IndexedShapeSet::new();
    for j in 1..=a_shapes_to_rest_map.extent() {
        let Some(a_s) = a_shapes_to_rest_map.find_key(j).cloned() else {
            continue;
        };
        let a_subsh_map = map_shapes_of(&a_s, a_subsh_enum);
        let Some(a_subsh) = a_subsh_map.find_key(1).cloned() else {
            // `aSubshMap(1)` on an empty map: OCCT reads out of range here.
            a_non_passed_shapes.add(&a_s);
            continue;
        };
        match a_map_ss.find_from_key(&a_subsh) {
            Some(&st) => {
                let mut a_state = st;
                if a_state == State::On {
                    a_state = a_shape_classifier
                        .state_shape_reference_null_avoid(&a_s)
                        .unwrap_or(State::Unknown);
                }
                let mut a_shape_with_state = ShapeWithState::new();
                a_shape_with_state.set_state(a_state);
                a_shape_with_state.set_is_splitted(false);
                a_map_of_shape_with_state.add(&a_s, a_shape_with_state);
            }
            None => {
                a_non_passed_shapes.add(&a_s);
            }
        }
    }

    // 5. shapes that could not be reached by propagation (`cxx:182-252`).
    if a_non_passed_shapes.extent() == 0 {
        return;
    }

    a_map_subsh_anc.clear();
    let mut non_passed: Vec<TopoShape> = Vec::new();
    for i in 1..=a_non_passed_shapes.extent() {
        if let Some(s) = a_non_passed_shapes.find_key(i) {
            non_passed.push(s.clone());
            map_shapes_and_ancestors_into(s, a_subsh_enum, a_shape_enum, &mut a_map_subsh_anc);
        }
    }

    a_map_ss.clear();
    for a_non_passed_shape in &non_passed {
        if a_map_ss.contains(a_non_passed_shape) {
            continue;
        }
        let a_state = find_state_through_vertex(
            a_non_passed_shape,
            a_shape_classifier,
            a_map_of_shape_with_state,
            an_avoid_subsh_map,
        );
        bind(&mut a_map_ss, a_non_passed_shape, a_state);

        // First sub-shape (`cxx:210-224`).
        let a_tmp_map = map_shapes_of(a_non_passed_shape, a_subsh_enum);
        let mut a_first_subsh: Option<TopoShape> = None;
        for j in 1..=a_tmp_map.extent() {
            if let Some(s) = a_tmp_map.find_key(j) {
                if !an_avoid_subsh_map.contains(s) {
                    a_first_subsh = Some(s.clone());
                    break;
                }
            }
        }
        let Some(a_first_subsh) = a_first_subsh else {
            continue;
        };
        bind(&mut a_map_ss, &a_first_subsh, a_state);

        let mut a_map_processed_subsh = IndexedShapeSet::new();
        if a_subsh_enum == ShapeType::Edge {
            find_state1(
                &a_first_subsh,
                a_state,
                &a_map_subsh_anc,
                &mut a_map_processed_subsh,
                &mut a_map_ss,
            );
        } else {
            find_state2(
                &a_first_subsh,
                a_state,
                &a_map_subsh_anc,
                &mut a_map_processed_subsh,
                &mut a_map_ss,
            );
        }
    }

    // Fill aShapeWithState (`cxx:240-251`); vertices are not recorded.
    let mut a_shape_with_state = ShapeWithState::new();
    a_shape_with_state.set_is_splitted(false);
    for k in 1..=a_map_ss.extent() {
        let (Some(key), Some(&st)) = (a_map_ss.find_key(k), a_map_ss.find_from_index(k)) else {
            continue;
        };
        a_shape_with_state.set_state(st);
        if key.shape_type() != ShapeType::Vertex {
            a_map_of_shape_with_state.add(key, a_shape_with_state.clone());
        }
    }
}

/// `PropagateStateForWires(aFacesToRestMap, aMapOfShapeWithState)`
/// (`cxx:412-457`): gives every wire (and every still-unknown edge) of a face
/// that already carries a state the state of that face.
pub fn propagate_state_for_wires(
    a_faces_to_rest_map: &IndexedShapeSet,
    a_map_of_shape_with_state: &mut IndexedShapeData<ShapeWithState>,
) {
    for i in 1..=a_faces_to_rest_map.extent() {
        let Some(a_f) = a_faces_to_rest_map.find_key(i).cloned() else {
            continue;
        };
        let Some(a_sws) = a_map_of_shape_with_state.find_from_key(&a_f).cloned() else {
            continue;
        };
        let a_st = a_sws.state();

        let a_wire_map = map_shapes_of(&a_f, ShapeType::Wire);
        for j in 1..=a_wire_map.extent() {
            let Some(a_w) = a_wire_map.find_key(j).cloned() else {
                continue;
            };
            let mut a_wire_sws = ShapeWithState::new();
            a_wire_sws.set_state(a_st);
            a_wire_sws.set_is_splitted(false);
            a_map_of_shape_with_state.add(&a_w, a_wire_sws);

            let a_edge_map = map_shapes_of(&a_w, ShapeType::Edge);
            for k in 1..=a_edge_map.extent() {
                let Some(a_e) = a_edge_map.find_key(k) else {
                    continue;
                };
                if !a_map_of_shape_with_state.contains(a_e) {
                    let mut an_edge_sws = ShapeWithState::new();
                    an_edge_sws.set_state(a_st);
                    an_edge_sws.set_is_splitted(false);
                    a_map_of_shape_with_state.add(a_e, an_edge_sws);
                }
            }
        }
    }
}
