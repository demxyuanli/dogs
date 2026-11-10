//! `TopOpeBRepBuild_GridFF` — the edge half of the `GTopo` grid.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridFF.cxx`
//! (OCCT 8.0.0) and `TopOpeBRepBuild_makeedges.cxx`:
//!
//! | OCCT | lines |
//! |---|---|
//! | `FUN_samgeomori(BDS, iref, ifil, samgeomori)` | `803-814` |
//! | `GSplitEdgeWES(EOR, LSclass, G1, WES)` | `826-1207` |
//! | `GMergeEdgeWES(EOR, G1, WES)` | `1217-1347` |
//! | `GSplitEdge(EOR, G1, LSclass)` | `1350-1485` |
//! | `ClassifyEdgeToSolidByOnePoint(E, Ref)` | `1493-1521` |
//! | `GPVSMakeEdges(EF, PVS, LOE)` | `makeedges.cxx:41-76` |
//! | `GEDBUMakeEdges(EF, EDBU, LOE)` | `makeedges.cxx:79-228` |
//!
//! `GLOBAL_classifysplitedge` (`GridFF.cxx:89`) is set by
//! `TopOpeBRepBuild_ffsfs.cxx:420`, `:492` and `:580`, a class the port does
//! not carry yet; the flag therefore stays `false` here and the `pls` choice
//! (`cxx:1439-1462`) always falls through to the `toclass` arm. The setter is
//! exposed so that `ffsfs` can drive it once it lands.
//!
//! UNPORTED: the `OCCT_DEBUG` blocks (`cxx:893-906`, `:963-971`, `:1007-1020`,
//! `:1052-1071`, `:1078-1083`, `:1121-1131`, `:1133-1147`, `:1200-1204`,
//! `:1221-1232`, `:1288-1299`, `:1320-1329`, `:1337-1343`, `:1363-1372`,
//! `:1379-1380`, `:1384-1389`, `makeedges.cxx:46-51`, `:84-100`) — trace,
//! dump and the `deb*` hooks only.

use std::cell::Cell;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, is_equal};
use crate::topopebrep::build::edge_builder::EdgeBuilder;
use crate::topopebrep::build::pave_classifier::PaveClassifier;
use crate::topopebrep::build::pave_set::PaveSet;
use crate::topopebrep::build::loop_set::LoopSetSource;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::complete_ds::{
    fun_ds_ori_e_in_f, ORI_E_IN_F_FOR_REV_OPPO,
    ORI_E_IN_F_ON_OPPO_SHA, ORI_E_IN_F_ON_SAME_SHA,
};
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::tool::{curves_same_oriented_point, fun_tool_bounds};

use super::super::g_topo::GTopo;
use super::Builder;

thread_local! {
    /// `GLOBAL_classifysplitedge` (`GridFF.cxx:89`).
    static GLOBAL_CLASSIFY_SPLIT_EDGE: Cell<bool> = const { Cell::new(false) };
}

/// Current `GLOBAL_classifysplitedge`.
pub fn global_classify_split_edge() -> bool {
    GLOBAL_CLASSIFY_SPLIT_EDGE.with(|c| c.get())
}

/// Assign `GLOBAL_classifysplitedge`.
pub fn set_global_classify_split_edge(v: bool) {
    GLOBAL_CLASSIFY_SPLIT_EDGE.with(|c| c.set(v));
}

/// `PAR_T` of `ClassifyEdgeToSolidByOnePoint` (`cxx:1495`): `10.*e^-PI`.
const PAR_T: f64 = 0.43213918;

impl Builder {
    /// `GFillEdgeWES(EOR, LSclass, G1, WES)` (`cxx:761-800`).
    ///
    /// `cxx:777-796` computes `GToSplit` / `GToMerge` for the trace block only,
    /// so the port does not evaluate them.
    pub fn g_fill_edge_wes(
        &mut self,
        eor: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        // `cxx:793-794`
        self.g_merge_edge_wes(eor, g1, wes);

        // `cxx:796-799`
        let mut gse = g1.clone();
        gse.change_config(Config::UnshGeometry, Config::UnshGeometry);
        self.g_split_edge_wes(eor, ls_class, &gse, wes);
    }

    /// `GSplitEdgeWES(EOR, LSclass, G1, WES)` (`cxx:826-1207`).
    pub fn g_split_edge_wes(
        &mut self,
        eor: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();
        let ori_e = eor.orientation();
        let new_ori_e = Builder::orient(ori_e, rev_ori1);

        // `cxx:838-857`
        let mut o_ein_f = Orientation::External;
        let oinref;
        let (hsdm, hsdm_e) = match self.my_data_structure.as_ref() {
            Some(hds) => {
                let hsdm = match self.my_face_to_fill.clone() {
                    Some(f) => hds.has_same_domain(&f, true),
                    None => false,
                };
                (hsdm, hds.has_same_domain(eor, true))
            }
            None => (false, false),
        };

        if hsdm {
            let (oref_face, ofill_face) = (
                self.my_face_reference.clone(),
                self.my_face_to_fill.clone(),
            );
            if let (Some(fref), Some(ffill)) = (oref_face, ofill_face) {
                let bds_ref = oinref_for(self, eor, &fref);
                oinref = bds_ref.0;
                o_ein_f = bds_ref.1;

                // `cxx:849-876`
                let (o_infill, o_einfill) = fun_ds_ori_e_in_f(
                    match self.my_data_structure.as_ref() {
                        Some(hds) => hds.ds(),
                        None => return,
                    },
                    eor,
                    &ffill,
                );
                if oinref == o_infill {
                    let mut reverse = false;
                    let (ifil, iref) = match self.my_data_structure.as_ref() {
                        Some(hds) => (
                            hds.shape_index(&ffill, true),
                            hds.shape_index(&fref, true),
                        ),
                        None => (0, 0),
                    };
                    if iref != ifil {
                        let oref = fref.orientation();
                        let same_geom_ori = {
                            let bds = match self.my_data_structure.as_ref() {
                                Some(hds) => hds.ds(),
                                None => return,
                            };
                            fun_samgeomori(bds, iref, ifil)
                        };
                        reverse = !same_geom_ori;
                        if oref == Orientation::Reversed {
                            reverse = !reverse;
                        }
                    }
                    let o_einfill_to_ref = if reverse {
                        o_einfill.complement()
                    } else {
                        o_einfill
                    };
                    let same = o_ein_f == o_einfill_to_ref;
                    if !same
                        && o_ein_f != Orientation::Internal
                        && o_ein_f != Orientation::External
                    {
                        o_ein_f = o_einfill_to_ref;
                    }
                }
            } else {
                oinref = 0;
            }
        } else {
            // `cxx:878-881`
            if let Some(ffill) = self.my_face_to_fill.clone() {
                let (o, ori) = fun_ds_ori_e_in_f(
                    match self.my_data_structure.as_ref() {
                        Some(hds) => hds.ds(),
                        None => return,
                    },
                    eor,
                    &ffill,
                );
                oinref = o;
                o_ein_f = ori;
            } else {
                oinref = 0;
            }
        }

        // `cxx:882-888`
        let mut new_o = oinref == ORI_E_IN_F_ON_SAME_SHA || oinref == ORI_E_IN_F_ON_OPPO_SHA;
        if self.my_is_k_part == 3 {
            new_o = false;
        }

        // `cxx:908-946`
        let tosplit = self.g_to_split(eor, tb1);
        if tosplit {
            self.g_split_edge(eor, g1, ls_class);

            if self.my_is_k_part == 4 {
                let rank1 = self.shape_rank(eor);
                let a_solid = if rank1 == 1 {
                    self.my_shape2.clone()
                } else {
                    self.my_shape1.clone()
                };
                let mut aux: Vec<TopoShape> = Vec::new();
                let split_list = self.change_split(eor, tb1).clone();
                for split_edge in &split_list {
                    if let Some(solid) = &a_solid {
                        let a_state = classify_edge_to_solid_by_one_point(split_edge, solid);
                        if a_state == tb1 {
                            aux.push(split_edge.clone());
                        }
                    }
                }
                let a_split_list = self.change_split(eor, tb1);
                a_split_list.clear();
                a_split_list.extend(aux);
            }
        }

        // `cxx:948-954`
        if self.opefus() && tb1 == State::In {
            return;
        }

        // `cxx:956-1022`
        let issplit = self.is_split(eor, tb1);
        if issplit {
            let lse = self.splits(eor, tb1).to_vec();
            for e in lse {
                let mut new_e = e;
                if new_o {
                    new_e.set_orientation(o_ein_f);
                    if !BRepTool::is_degenerated(&Edge(eor.clone())) && hsdm_e {
                        orient_sdm_split(&mut new_e, eor);
                    }
                } else {
                    new_e.set_orientation(new_ori_e);
                }
                wes.add_start_element(&new_e);
            }
        } else {
            // `cxx:1026-1102`
            let se = match self.my_data_structure.as_ref() {
                Some(hds) => hds.ds().is_section_edge(eor, true),
                None => false,
            };
            let hs = match self.my_data_structure.as_ref() {
                Some(hds) => hds.has_shape(eor, true),
                None => false,
            };
            let hg = match self.my_data_structure.as_ref() {
                Some(hds) => hds.has_geometry(eor),
                None => false,
            };
            let mut add = false;
            let mut add_on = false;
            let isstart = hs;

            if se {
                let mut ftg = !ls_class.is_empty();
                let tclass = ls_class
                    .first()
                    .map(|s| s.shape_type())
                    .unwrap_or(ShapeType::Shape);
                ftg = ftg && tclass == ShapeType::Face;
                if !ftg {
                    let mut pos = State::Unknown;
                    let keepse = self.g_keep_shape1(eor, ls_class, tb1, &mut pos);
                    if keepse {
                        add = true;
                    } else if self.my_process_on && pos == State::On {
                        add_on = true;
                    }
                }
            } else {
                add = true;
                let testkeep = hs && !hg;
                if testkeep {
                    let mut pos = State::Unknown;
                    let mut keep = self.g_keep_shape1(eor, ls_class, tb1, &mut pos);
                    if !keep {
                        let mut test_on = !ls_class.is_empty();
                        if test_on {
                            test_on = ls_class
                                .first()
                                .map(|s| s.shape_type() == ShapeType::Solid)
                                .unwrap_or(false);
                        }
                        if test_on {
                            keep = pos == State::On;
                        }
                        add_on = self.my_process_on && keep;
                    }
                    add = keep;
                }
            }

            // `cxx:1104-1155`
            if add {
                let mut new_e = eor.clone();
                if new_o {
                    new_e.set_orientation(o_ein_f);
                } else if oinref == ORI_E_IN_F_FOR_REV_OPPO {
                    // `cxx:1112`: recorded `INTERNAL` in OCCT (PRO14785).
                    new_e.set_orientation(Orientation::Internal);
                } else {
                    new_e.set_orientation(new_ori_e);
                }

                if isstart {
                    wes.add_start_element(&new_e);
                } else {
                    wes.add_element(&new_e);
                }
            }

            if add_on {
                let mut new_e = eor.clone();
                new_e.set_orientation(new_ori_e);
                self.my_on_elem_map.add(&new_e);
            }
        }

        // `cxx:1161-1198`
        if self.my_process_on && self.is_split(eor, State::On) {
            let lse = self.splits(eor, State::On).to_vec();
            for e in lse {
                let mut new_e = e;
                if new_o {
                    new_e.set_orientation(o_ein_f);
                    if !BRepTool::is_degenerated(&Edge(eor.clone())) && hsdm_e {
                        orient_sdm_split(&mut new_e, eor);
                    }
                } else {
                    new_e.set_orientation(new_ori_e);
                }
                self.my_on_elem_map.add(&new_e);
            }
        }
    }

    /// `GMergeEdgeWES(EOR, G1, WES)` (`cxx:1217-1347`).
    pub fn g_merge_edge_wes(&mut self, eor: &TopoShape, g1: &GTopo, wes: &mut WireEdgeSet) {
        // `cxx:1236-1239`
        let closing = match self.my_face_to_fill.clone() {
            Some(f) => BRepTool::is_closed_edge_face(&Edge(eor.clone()), &Face(f)),
            None => false,
        };
        if closing {
            return;
        }

        // `cxx:1241-1244`
        if self.opefus() {
            return;
        }

        // `cxx:1248-1251`
        let (_tb1, tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();
        let ori_e = eor.orientation();
        let new_ori_e = Builder::orient(ori_e, rev_ori1);

        // `cxx:1252-1257`
        let hassame = match self.my_data_structure.as_ref() {
            Some(hds) => hds.has_same_domain(eor, true),
            None => false,
        };
        if !hassame {
            return;
        }

        // `cxx:1259-1266`
        let eisref = match self.my_data_structure.as_ref() {
            Some(hds) => {
                let i_eref = hds.same_domain_reference(eor);
                hds.shape(i_eref, true)
                    .map(|s| is_equal(eor, &s))
                    .unwrap_or(false)
            }
            None => false,
        };

        // `cxx:1268-1271`
        let tbeor = if eisref { _tb1 } else { tb2 };
        if tbeor == State::Out {
            return;
        }

        // `cxx:1273-1302`
        let ismerged = self.is_merged(eor, tbeor);
        if ismerged {
            if !eisref {
                return;
            }
            let me = self.merged(eor, tbeor).to_vec();
            for e in me {
                let mut new_e = e;
                new_e.set_orientation(new_ori_e);
                wes.add_start_element(&new_e);
            }
            return;
        }

        // `cxx:1304-1314`
        *self.change_merged(eor, tbeor) = Vec::new();
        let stsp_eor = tbeor;
        let issplit = self.is_split(eor, stsp_eor);
        if !issplit {
            return;
        }

        // `cxx:1317`
        let splits = self.splits(eor, stsp_eor).to_vec();
        *self.change_merged(eor, tbeor) = splits;

        // `cxx:1319-1336`
        let me = self.merged(eor, tbeor).to_vec();
        for e in me {
            let mut new_e = e;
            new_e.set_orientation(new_ori_e);
            wes.add_start_element(&new_e);
        }
    }

    /// `GSplitEdge(EOR, G1, LSclass)` (`cxx:1350-1485`).
    pub fn g_split_edge(&mut self, eor: &TopoShape, g1: &GTopo, ls_class: &[TopoShape]) {
        let (tb1, _tb2) = g1.states_on();

        // `cxx:1358-1360`: work on a FORWARD edge.
        let mut ef = eor.clone();
        ef.set_orientation(Orientation::Forward);

        // `cxx:1374-1378`: `takeON` is computed and then unconditionally set
        // to `false` in OCCT, so the `takeON` arm below is dead.
        let _take_on = false;
        let _ = tb1;

        // `cxx:1399-1415`
        let mut pvs = PaveSet::new(&ef);
        self.my_edge_reference = Some(ef.clone());
        self.fill_point_topology_pvs(&ef, g1, &mut pvs);
        self.mark_split(&ef, tb1, true);

        let mut loe: Vec<TopoShape> = Vec::new();
        self.gpvs_make_edges(&ef, &mut pvs, &mut loe);

        if loe.is_empty() {
            return;
        }

        // `cxx:1419-1462`
        let sel = self.change_split(&ef, tb1);
        sel.clear();

        let c1 = g1.config1();
        let c2 = g1.config2();
        let _uuface = c1 == Config::UnshGeometry && c2 == Config::UnshGeometry;

        let on_solid = ls_class
            .first()
            .map(|s| s.shape_type() == ShapeType::Solid)
            .unwrap_or(false);
        // `cxx:1435-1436`: `toclass = UUFACE; toclass = !ONSOLID;` — the first
        // assignment is dead.
        let toclass = !on_solid;

        let mut loos: Vec<TopoShape> = Vec::new();
        let pls: &[TopoShape] = if global_classify_split_edge() {
            let r = self.shape_rank(eor);
            let oos = if r == 1 {
                self.my_shape2.clone()
            } else {
                self.my_shape1.clone()
            };
            if let Some(oos) = oos {
                loos.push(oos);
            }
            &loos
        } else if toclass {
            ls_class
        } else {
            &[]
        };

        // `cxx:1464-1477`
        let mut a_lon: Vec<TopoShape> = Vec::new();
        let mut sel_edges: Vec<TopoShape> = Vec::new();
        let process_on = self.my_process_on;
        for a_e in loe {
            let mut pos = State::Unknown;
            let keep = self.g_keep_shape1(&a_e, pls, tb1, &mut pos);
            if keep {
                sel_edges.push(a_e);
            } else if process_on && pos == State::On {
                a_lon.push(a_e);
            }
        }
        {
            let sel = self.change_split(&ef, tb1);
            sel.extend(sel_edges);
        }

        // `cxx:1479-1484`
        if !a_lon.is_empty() {
            self.mark_split(&ef, State::On, true);
            let a_slon = self.change_split(&ef, State::On);
            a_slon.clear();
            a_slon.extend(a_lon);
        }
    }

    /// `GPVSMakeEdges(EF, PVS, LOE)` (`makeedges.cxx:41-76`).
    pub fn gpvs_make_edges(&self, ef: &TopoShape, pvs: &mut PaveSet, loe: &mut Vec<TopoShape>) {
        // `cxx:54-58`
        let mut vcl = PaveClassifier::new(ef);
        let equalpar = pvs.has_equal_parameters();
        if equalpar {
            let p = pvs.equal_parameters();
            vcl.set_first_parameter(p);
        }

        // `cxx:60-68`
        pvs.init_loop();
        let novertex = !pvs.more_loop();
        if novertex {
            return;
        }

        // `cxx:70-74`
        let mut edbu = EdgeBuilder::new();
        let force_class = false;
        edbu.init_edge_builder(pvs, &mut vcl, force_class);
        self.gedbu_make_edges(ef, &mut edbu, loe);
    }
}

/// `FUN_ds_oriEinF(BDS, Edge(EOR), myFaceReference, oEinF)` (`cxx:845`): the
/// `(Oinref, oEinF)` pair of one face.
fn oinref_for(b: &Builder, e: &TopoShape, f: &TopoShape) -> (i32, Orientation) {
    match b.my_data_structure.as_ref() {
        Some(hds) => fun_ds_ori_e_in_f(hds.ds(), e, f),
        None => (0, Orientation::External),
    }
}

/// `FUN_samgeomori(BDS, iref, ifil, samgeomori)` (`cxx:803-814`).
fn fun_samgeomori(
    bds: &crate::topopebrep::ds::data_structure::DataStructure,
    iref: i32,
    ifil: i32,
) -> bool {
    let cfill = bds.same_domain_ori(ifil);
    let oref = bds
        .shape(iref, true)
        .map(|s| s.orientation())
        .unwrap_or(Orientation::Forward);
    let ofil = bds
        .shape(ifil, true)
        .map(|s| s.orientation())
        .unwrap_or(Orientation::Forward);
    let mut samgeomori = cfill == Config::SameOriented;
    if oref == ofil.complement() {
        samgeomori = !samgeomori;
    }
    samgeomori
}

/// `cxx:981-1001` / `:1171-1191`: re-orient a same-domain split of an edge so
/// that it runs the same way as the original edge.
fn orient_sdm_split(new_e: &mut TopoShape, eor: &TopoShape) {
    let (f, l) = fun_tool_bounds(new_e);
    let x = 0.45678;
    let par = (1.0 - x) * f + x * l;
    // `cxx:989-999`: OCCT warns and keeps going when the projection fails.
    if let Some(so) = curves_same_oriented_point(new_e, par, eor) {
        if !so {
            new_e.reverse();
        }
    }
}

/// `ClassifyEdgeToSolidByOnePoint(E, Ref)` (`cxx:1493-1521`).
pub fn classify_edge_to_solid_by_one_point(e: &TopoShape, reference: &TopoShape) -> State {
    let (f2, l2) = BRepTool::edge_parameters(&Edge(e.clone()));

    let a_p3d = if BRepTool::edge_curve(&Edge(e.clone())).is_none() {
        // `cxx:1501-1509`: a degenerated edge is classified at its first
        // vertex.
        let (fv, _) = edge_vertices(&Edge(e.clone()));
        match fv {
            Some(v) => BRepTool::vertex_point(&v),
            None => return State::Unknown,
        }
    } else {
        // `cxx:1512-1516`
        let par = f2 * PAR_T + (1.0 - PAR_T) * l2;
        match crate::topopebrep::tool::fun_tool_value(par, e) {
            Some(p) => p,
            None => return State::Unknown,
        }
    };

    let mut sc = crate::brep_class3d::SolidClassifier::new();
    sc.load(reference.clone());
    sc.perform(&a_p3d, 1.0e-7);
    match sc.state() {
        crate::fclass2d::FaceState::In => State::In,
        crate::fclass2d::FaceState::Out => State::Out,
        crate::fclass2d::FaceState::On => State::On,
        crate::fclass2d::FaceState::Unknown => State::Unknown,
    }
}
