//! `TopOpeBRepBuild_Builder1` — splitting the edges of a face into their
//! state parts.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1_1.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `StatusEdgesToSplit(anObj, anEdgesToSplitMap, anEdgesToRestMap)` | `556-855` |
//! | `SplitEdge(anEdge, aLNew, aDataMapOfShapeState)` | `859-1115` |
//! | `ClassifyEdgeToSolidByOnePoint(E, SC)` | `1117-1144` |
//!
//! `StatusEdgesToSplit` gives every interfered edge of a wire its split parts
//! and their states, then hands the untouched edges to `Tools::PropagateState`
//! so they get a state too.
//!
//! `SplitEdge` cuts one edge at the paves `GFillPointTopologyPVS` collected.
//! Two paves that land on the same parameter are merged first
//! (`cxx:898-1002`): the surviving vertex is the one that is same-domain,
//! already used, or from a 3-D interference, in that order, and the dropped
//! vertex's orientation is remembered so the second pass can restore it
//! (`cxx:1004-1052`). The parameter comparison uses the edge resolution at
//! `max(tolEdge, tolV1, tolV2)` and `Precision::PConfusion`, exactly as OCCT
//! writes it.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::g_tool::GTool;
use crate::topopebrep::build::loop_set::LoopSetSource;
use crate::topopebrep::build::loop_shape::LoopHandle;
use crate::topopebrep::build::pave_set::PaveSet;
use crate::topopebrep::build::tools::propagate_state;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::list_of_shape_on1_state::ListOfShapeOn1State;
use crate::topopebrep::ds::ShapeWithState;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;
use crate::{Edge, Vertex};

use super::super::grid_ee::{global_faces2d, set_global_faces2d};
use super::shape_states::bind_state;
use super::super::Builder;
use super::{
    processed_edges_add, processed_edges_contains, unk_state_ver_add, used_vertex_add,
    used_vertex_contains, PAR_T, THE_UNK_STATE_VER,
};

/// `TopAbs_Orientation` to the integer OCCT stores in `aVerOriMap`
/// (`cxx:894`).
fn ori_to_i(o: Orientation) -> i32 {
    match o {
        Orientation::Forward => 0,
        Orientation::Reversed => 1,
        Orientation::Internal => 2,
        Orientation::External => 3,
    }
}

/// The inverse of [`ori_to_i`].
fn ori_from_i(v: i32) -> Orientation {
    match v {
        1 => Orientation::Reversed,
        2 => Orientation::Internal,
        3 => Orientation::External,
        _ => Orientation::Forward,
    }
}

/// The `(vertex, parameter, interference type)` of a `Pave` handle.
fn pave_of(handle: &LoopHandle) -> Option<(TopoShape, f64, Kind)> {
    let p = match handle.read() {
        Ok(p) => p,
        Err(p) => p.into_inner(),
    };
    let pave = p.as_pave()?;
    Some((pave.vertex.clone(), pave.parameter(), pave.interference_type()))
}

impl Builder {
    /// `StatusEdgesToSplit(anObj, anEdgesToSplitMap, anEdgesToRestMap)`
    /// (`cxx:556-855`).
    ///
    /// `map` is the `aMapOfShapeWithState = aDS.ChangeMapOfShapeWithState(
    /// anObj, aFlag)` of `cxx:565-566`; the caller already holds it (see
    /// `PerformFacesWithStates`).
    pub(super) fn status_edges_to_split(
        &mut self,
        _an_obj: &TopoShape,
        an_edges_to_split_map: &IndexedShapeSet,
        an_edges_to_rest_map: &IndexedShapeSet,
        map: &mut IndexedShapeData<ShapeWithState>,
    ) {
        // `cxx:573-577`
        let n_e = an_edges_to_split_map.extent();
        if n_e == 0 {
            // `cxx:575-578` returns without touching the map.
            return;
        }

        // `cxx:579-580`
        let mut a_spl_edges_state = IndexedShapeData::<State>::new();
        let sc = Rc::clone(&self.my_shape_classifier);

        for i in 1..=n_e {
            let Some(an_edge) = an_edges_to_split_map.find_key(i).cloned() else {
                continue;
            };

            // `cxx:588-618`: an edge already processed in the other direction
            // only contributes the states of its known split parts.
            if processed_edges_contains(&an_edge) {
                if let Some(sws) = map.find_from_key(&an_edge).cloned() {
                    if sws.is_splitted() {
                        for st in [State::On, State::Out, State::In] {
                            for v in sws.part(st) {
                                bind_state(&mut a_spl_edges_state, v, st);
                            }
                        }
                    }
                }
                continue;
            }

            // `cxx:620-622`
            processed_edges_add(&an_edge);

            let mut a_shape_with_state = ShapeWithState::new();

            // `cxx:626-663`: an edge from same-domain faces was split already.
            let is_split_on = self.is_split(&an_edge, State::On);
            if is_split_on {
                for st in [State::On, State::In, State::Out] {
                    let parts = self.splits(&an_edge, st).to_vec();
                    for part in &parts {
                        let mut ne = part.clone();
                        ne.set_orientation(an_edge.orientation());
                        a_shape_with_state.add_part(ne.clone(), st);
                        bind_state(&mut a_spl_edges_state, &ne, st);
                    }
                }
                a_shape_with_state.set_is_splitted(true);
                map.add(&an_edge, a_shape_with_state);
                continue;
            }

            // `cxx:668-676`: attempt to split the edge.
            let mut a_data_map_of_shape_state = IndexedShapeData::<State>::new();
            let mut a_l_new: Vec<TopoShape> = Vec::new();

            let old_state = global_faces2d();
            set_global_faces2d(true);
            self.split_edge_pvs(&an_edge, &mut a_l_new, &mut a_data_map_of_shape_state);
            set_global_faces2d(old_state);

            if a_l_new.is_empty() {
                // `cxx:678-761`: no new edge means the whole edge is `IN`,
                // and the Fuse operation dropped it.
                let (vf0, vl0) = edge_vertices(&Edge(an_edge.clone()));
                let (Some(vf0), Some(vl0)) = (vf0, vl0) else {
                    map.add(&an_edge, a_shape_with_state);
                    continue;
                };
                let (mut vf0, mut vl0): (TopoShape, TopoShape) = (vf0.into(), vl0.into());

                // `cxx:690-698`: `myDataStructure->HasSameDomain(Vf/Vl/aNewEdge)`.
                let (has_sdv1, has_sdv2, has_sd_edge) = {
                    let Some(hds) = self.my_data_structure.as_ref() else {
                        return;
                    };
                    (
                        hds.has_same_domain(&vf0, true),
                        hds.has_same_domain(&vl0, true),
                        hds.has_same_domain(&an_edge, true),
                    )
                };

                // An edge with same-domain partners must have come through
                // `SplitSectionEdges`; nothing is done for it here.
                let (has_sdv1, has_sdv2) = if has_sd_edge {
                    (false, false)
                } else {
                    (has_sdv1, has_sdv2)
                };

                let mut a_new_edge = an_edge.clone();
                if has_sdv1 || has_sdv2 {
                    let mut eor = an_edge.clone();
                    eor.set_orientation(Orientation::Forward);
                    let par_f = BRepTool::parameter_of_vertex(&eor, &vf0);
                    let par_l = BRepTool::parameter_of_vertex(&eor, &vl0);
                    let mut ne = self.my_build_tool.copy_edge(&eor);

                    if has_sdv1 {
                        if let Some(s) = self.same_domain_reference_shape(&vf0) {
                            vf0 = s;
                        }
                        vf0.set_orientation(Orientation::Forward);
                    }
                    if has_sdv2 {
                        if let Some(s) = self.same_domain_reference_shape(&vl0) {
                            vl0 = s;
                        }
                        vl0.set_orientation(Orientation::Reversed);
                    }

                    self.my_build_tool.add_edge_vertex(&mut ne, &vf0);
                    self.my_build_tool.parameter(&ne, &vf0, par_f);
                    self.my_build_tool.add_edge_vertex(&mut ne, &vl0);
                    self.my_build_tool.parameter(&ne, &vl0, par_l);

                    ne.set_orientation(an_edge.orientation());
                    a_new_edge = ne;
                }

                // `cxx:731-739`
                let a_state = classify_edge_to_solid_by_one_point(&Edge(a_new_edge.clone()), &sc);
                a_shape_with_state.set_is_splitted(true);
                a_shape_with_state.add_part(a_new_edge.clone(), a_state);
                bind_state(&mut a_spl_edges_state, &a_new_edge, a_state);

                let (vf, vl) = edge_vertices(&Edge(a_new_edge.clone()));
                if let Some(vf) = &vf {
                    used_vertex_add(&vf.0);
                }
                if let Some(vl) = &vl {
                    used_vertex_add(&vl.0);
                }

                // `cxx:741-760`: a vertex whose own state disagrees with the
                // edge is marked unknown so the propagation does not use it.
                if !BRepTool::is_degenerated(&Edge(a_new_edge.clone())) {
                    note_unknown_vertex_state(&sc, vf.as_ref(), vl.as_ref(), a_state);
                }
            } else {
                // `cxx:762-808`: the usual case, the edge was split.
                for a_s in a_l_new.clone() {
                    let mut a_state = a_data_map_of_shape_state
                        .find_from_key(&a_s)
                        .copied()
                        .unwrap_or(State::Unknown);
                    let isdegen = BRepTool::is_degenerated(&Edge(a_s.clone()));

                    // `cxx:770-781`: an `IN` state from the split is not proof;
                    // classify the piece again.
                    if a_state == State::In && !isdegen {
                        a_state = classify_edge_to_solid_by_one_point(&Edge(a_s.clone()), &sc);
                    }

                    a_shape_with_state.add_part(a_s.clone(), a_state);
                    a_shape_with_state.set_is_splitted(true);
                    bind_state(&mut a_spl_edges_state, &a_s, a_state);

                    let (vf, vl) = edge_vertices(&Edge(a_s.clone()));
                    if let Some(vf) = &vf {
                        used_vertex_add(&vf.0);
                    }
                    if let Some(vl) = &vl {
                        used_vertex_add(&vl.0);
                    }
                    if !isdegen {
                        note_unknown_vertex_state(&sc, vf.as_ref(), vl.as_ref(), a_state);
                    }
                }
            }

            // `cxx:810-822`: keep the ON parts as splits of this edge and
            // register the edge as a section edge.
            let esp_on = a_shape_with_state.part(State::On).to_vec();
            if !is_split_on && !esp_on.is_empty() {
                let mut on_spl = ListOfShapeOn1State::new();
                on_spl.change_list_on_state().clear();
                on_spl.change_list_on_state().extend(esp_on);
                on_spl.split(true);
                self.my_split_on.insert(shape_map_key(&an_edge), on_spl);
                if let Some(hds) = self.my_data_structure.as_mut() {
                    hds.change_ds().add_section_edge(&an_edge);
                }
            }

            map.add(&an_edge, a_shape_with_state);
        }

        // `cxx:826-843`: the untouched edges only contribute the state they
        // already carry.
        let n_e = an_edges_to_rest_map.extent();
        for i in 1..=n_e {
            let Some(an_edge) = an_edges_to_rest_map.find_key(i) else {
                continue;
            };
            if let Some(sws) = map.find_from_key(an_edge).cloned() {
                if !sws.is_splitted() {
                    let a_state = sws.state();
                    bind_state(&mut a_spl_edges_state, an_edge, a_state);
                    continue;
                }
            }
        }

        // `cxx:845-855`
        if n_e != 0 {
            let avoid = THE_UNK_STATE_VER.with(|m| m.borrow().clone());
            propagate_state(
                &a_spl_edges_state,
                an_edges_to_rest_map,
                ShapeType::Vertex,
                ShapeType::Edge,
                &mut sc.borrow_mut(),
                map,
                &avoid,
            );
        }
    }

    /// `myDataStructure->Shape(SameDomainReference(V))`.
    fn same_domain_reference_shape(&self, v: &TopoShape) -> Option<TopoShape> {
        let hds = self.my_data_structure.as_ref()?;
        let iref = hds.same_domain_reference(v);
        hds.shape(iref, true)
    }

    /// `SplitEdge(anEdge, aLNew, aDataMapOfShapeState)` (`cxx:859-1114`).
    ///
    /// OCCT spells the second pass as a `NCollection_List` iterator whose
    /// `More()` is tested before reading the next item; the port walks
    /// `a_pvs_list` by index so the "append the last pave and break" shape of
    /// `cxx:898-910` is kept.
    pub(super) fn split_edge_pvs(
        &mut self,
        an_edge: &TopoShape,
        a_l_new: &mut Vec<TopoShape>,
        a_data_map_of_shape_state: &mut IndexedShapeData<State>,
    ) {
        // `cxx:866-869`: the paves are only meaningful on a FORWARD edge.
        let mut edge_f = an_edge.clone();
        edge_f.set_orientation(Orientation::Forward);

        // `cxx:871-879`
        let mut pvs = PaveSet::new(&edge_f);
        let g1 = GTool::g_fus_same(ShapeType::Face, ShapeType::Face);
        self.my_edge_reference = Some(edge_f.clone());
        self.fill_point_topology_pvs(&edge_f, &g1, &mut pvs);

        pvs.init_loop();

        // `cxx:888-892`
        let mut a_pvs_list: Vec<LoopHandle> = Vec::new();
        let mut a_ver_ori_map: HashMap<ShapeMapKey, Vec<i32>> = HashMap::new();

        let edge_obj = Edge(an_edge.clone());
        let tol_edge = BRepTool::edge_tolerance(&edge_obj);
        let a_curve_adaptor = BRepTool::edge_curve(&edge_obj);

        // `cxx:898-1002`: merge paves with equal parameters.
        while pvs.more_loop() {
            let Some(pave1) = pvs.loop_() else {
                break;
            };
            let Some((a_v1, a_par1, int_type1)) = pave_of(&pave1) else {
                pvs.next_loop();
                continue;
            };
            let an_ori1 = a_v1.orientation();

            pvs.next_loop();
            if !pvs.more_loop() {
                a_pvs_list.push(pave1);
                break;
            }

            let Some(pave2) = pvs.loop_() else {
                break;
            };
            let Some((a_v2, a_par2, _int_type2)) = pave_of(&pave2) else {
                a_pvs_list.push(pave1);
                continue;
            };
            let an_ori2 = a_v2.orientation();

            let tol_v1 = BRepTool::vertex_tolerance(&Vertex(a_v1.clone()));
            let tol_v2 = BRepTool::vertex_tolerance(&Vertex(a_v2.clone()));
            let tol_max = tol_edge.max(tol_v1.max(tol_v2));
            let resol = match &a_curve_adaptor {
                Some(c) => c.resolution(tol_max),
                None => 0.0,
            };
            let delta = (a_par1 - a_par2).abs();

            if delta < resol {
                let dist = BRepTool::vertex_point(&Vertex(a_v1.clone()))
                    .distance(&BRepTool::vertex_point(&Vertex(a_v2.clone())));
                if dist < tol_max || delta < PCONFUSION {
                    let int_3d1 = int_type1 == Kind::Face;

                    let (Some(hds), used_v1, used_v2) = (
                        self.my_data_structure.as_ref(),
                        used_vertex_contains(&a_v1),
                        used_vertex_contains(&a_v2),
                    ) else {
                        return;
                    };
                    let has_sdv1 = hds.has_same_domain(&a_v1, true);
                    let has_sdv2 = hds.has_same_domain(&a_v2, true);

                    // `cxx:936-957`, in that order.
                    let take_first = if has_sdv1 {
                        true
                    } else if has_sdv2 {
                        false
                    } else if used_v1 {
                        true
                    } else if used_v2 {
                        false
                    } else {
                        int_3d1
                    };

                    let (mut a_ver, has_sdv, an_ori_opp) = if take_first {
                        (a_v1.clone(), has_sdv1, an_ori2)
                    } else {
                        (a_v2.clone(), has_sdv2, an_ori1)
                    };

                    // `cxx:978-996`: remember the dropped orientation.
                    if an_ori1 != an_ori2 {
                        let key = shape_map_key(&a_ver);
                        let list = a_ver_ori_map.entry(key).or_default();
                        list.push(i32::from(take_first));
                        list.push(ori_to_i(an_ori_opp));

                        if has_sdv {
                            let iref = hds.same_domain_reference(&a_ver);
                            if let Some(s) = hds.shape(iref, true) {
                                a_ver = s;
                            }
                        }
                        unk_state_ver_add(&a_ver);
                    }

                    a_pvs_list.push(if take_first { pave1 } else { pave2 });

                    pvs.next_loop();
                    continue;
                }
            }
            a_pvs_list.push(pave1);
        }

        // `cxx:1004-1052`: build the pieces between consecutive paves.
        let mut it = 0usize;
        while it < a_pvs_list.len() {
            let pave1 = a_pvs_list[it].clone();
            let Some((mut a_v1, a_par1, _)) = pave_of(&pave1) else {
                it += 1;
                continue;
            };
            // `cxx:1010-1011`: the orientation is read before the local copy
            // is forced FORWARD.
            let mut an_or1 = a_v1.orientation();
            a_v1.set_orientation(Orientation::Forward);
            if let Some(list) = a_ver_ori_map.get_mut(&shape_map_key(&a_v1)) {
                if list.len() >= 2 {
                    if list[0] != 0 {
                        an_or1 = ori_from_i(list[1]);
                    }
                    list.drain(0..2);
                }
            }

            it += 1;
            if it >= a_pvs_list.len() {
                break;
            }

            let pave2 = a_pvs_list[it].clone();
            let Some((mut a_v2, a_par2, _)) = pave_of(&pave2) else {
                continue;
            };
            // `cxx:1039-1040`
            let mut an_or2 = a_v2.orientation();
            a_v2.set_orientation(Orientation::Reversed);
            if let Some(list) = a_ver_ori_map.get_mut(&shape_map_key(&a_v2)) {
                if !list.is_empty() {
                    if list[0] == 0 {
                        an_or2 = ori_from_i(list[1]);
                    }
                }
            }

            // `cxx:1055-1058`: avoid creating an edge with an invalid range.
            if a_par1 > a_par2 {
                continue;
            }

            // `cxx:1060-1076`
            let (has_sdv1, has_sdv2) = match self.my_data_structure.as_ref() {
                Some(hds) => (
                    hds.has_same_domain(&a_v1, true),
                    hds.has_same_domain(&a_v2, true),
                ),
                None => return,
            };
            if has_sdv1 {
                if let Some(s) = self.same_domain_reference_shape(&a_v1) {
                    a_v1 = s;
                }
                a_v1.set_orientation(Orientation::Forward);
            }
            if has_sdv2 {
                if let Some(s) = self.same_domain_reference_shape(&a_v2) {
                    a_v2 = s;
                }
                a_v2.set_orientation(Orientation::Reversed);
            }

            // `cxx:1078-1085`
            let mut a_new_edge = self.my_build_tool.copy_edge(&edge_f);
            self.my_build_tool.add_edge_vertex(&mut a_new_edge, &a_v1);
            self
                .my_build_tool
                .parameter(&a_new_edge, &a_v1, a_par1);
            self.my_build_tool.add_edge_vertex(&mut a_new_edge, &a_v2);
            self
                .my_build_tool
                .parameter(&a_new_edge, &a_v2, a_par2);

            // `cxx:1087-1106`: the piece is `IN` unless the two ends point
            // inward.
            let mut a_state = State::In;
            if (an_or1 == Orientation::Forward && an_or2 == Orientation::Reversed)
                || (an_or1 == Orientation::Forward && an_or2 == Orientation::Internal)
                || (an_or1 == Orientation::Internal && an_or2 == Orientation::Reversed)
                || (an_or1 == Orientation::Internal && an_or2 == Orientation::Internal)
            {
                a_state = State::Out;
            }

            // `cxx:1109-1112`
            a_new_edge.set_orientation(an_edge.orientation());
            a_l_new.push(a_new_edge.clone());
            bind_state(a_data_map_of_shape_state, &a_new_edge, a_state);
        }
        // `cxx:1114`: `GEDBUMakeEdges` is commented out in OCCT.
    }
}

/// `cxx:741-760` / `cxx:792-807`: a vertex whose classified state disagrees
/// with the edge, and is not `ON`, is added to `theUnkStateVer`.
fn note_unknown_vertex_state(
    sc: &Rc<RefCell<ShapeClassifier>>,
    vf: Option<&Vertex>,
    vl: Option<&Vertex>,
    a_state: State,
) {
    if let Some(vf) = vf {
        sc.borrow_mut()
            .state_p3d_reference(&BRepTool::vertex_point(vf))
            .ok();
        let sv = sc.borrow().state();
        if sv != a_state && sv != State::On {
            unk_state_ver_add(&vf.0);
        }
        if let Some(vl) = vl {
            if !is_same(&vf.0, &vl.0) {
                sc.borrow_mut()
                    .state_p3d_reference(&BRepTool::vertex_point(vl))
                    .ok();
                let sv = sc.borrow().state();
                if sv != a_state && sv != State::On {
                    unk_state_ver_add(&vl.0);
                }
            }
        }
    }
}

/// `static TopAbs_State ClassifyEdgeToSolidByOnePoint(E, SC)`
/// (`cxx:1118-1143`).
fn classify_edge_to_solid_by_one_point(e: &Edge, sc: &Rc<RefCell<ShapeClassifier>>) -> State {
    let (f2, l2) = BRepTool::edge_parameters(e);
    let a_p3d = match BRepTool::edge_curve(e) {
        Some(c) => {
            // `cxx:1136-1137`
            let par = f2 * PAR_T + (1.0 - PAR_T) * l2;
            c.d0(par)
        }
        None => {
            // `cxx:1123-1134`: a degenerated edge has no curve; its first
            // vertex answers instead.
            let (fv, _) = edge_vertices(e);
            match fv {
                Some(fv) => BRepTool::vertex_point(&fv),
                None => return State::Unknown,
            }
        }
    };
    sc.borrow_mut().state_p3d_reference(&a_p3d).ok();
    sc.borrow().state()
}

