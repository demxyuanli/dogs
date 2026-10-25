//! `FUN_ds_complete1dForSESDM`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:3777-3952`.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topopebrep::ds::export::fds_assign;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::interference_tool::make_edge_vertex_interference;
use crate::topopebrep::ds::kind::{Config, Kind};
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_gk_interference, fun_select_itrasha_interference, fun_select_si_interference,
    fun_select_sk_interference, fun_select_traori_interference, fun_select_traunk_interference,
    fun_select_trasha_interference,
};
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::project_curve::fun_tool_par_e_tol;
use crate::topopebrep::tool::tangent::curves_same_oriented;

/// `FUN_ds_complete1dForSESDM(HDS)` (`cxx:3777-3952`).
///
/// MSV 25.03.2002 / OCC251: `SE` is a section edge same-domain with `Esd`; for
/// one of the vertices of `Esd` the 1d interference `(ori, v, Esd)` is missing
/// on `SE` and is added.
pub fn fun_ds_complete1d_for_se_sdm(hds: &mut HDataStructure) {
    let nse = hds.ds().nb_section_edges();

    for is in 1..=nse {
        // `cxx:3788-3802`
        let Some(se) = hds.ds().section_edge(is, true) else {
            continue;
        };
        if BRepTool::is_degenerated(&Edge(se.clone())) {
            continue;
        }
        let rk_se = hds.ds().ancestor_rank_of(&se);
        let i_se = hds.ds().shape_index(&se, true);
        if !hds.has_same_domain(&se, true) {
            continue;
        }
        let l_esd = hds.ds().shape_same_domain(&se).to_vec();
        if l_esd.is_empty() {
            continue;
        }

        // `cxx:3807-3825` : `aGBMap` = the vertices of `SE` and their
        // same-domain partners
        let tol_se = BRepTool::edge_tolerance(&Edge(se.clone()));
        let mut a_gb_map: Vec<TopoShape> = Vec::new();
        let (vse0, vse1) = crate::topo_tools::edge_vertices(&Edge(se.clone()));
        for vse in [vse0, vse1].into_iter().flatten() {
            push_unique(&mut a_gb_map, &vse.0);
            for v in hds.ds().shape_same_domain(&vse.0) {
                push_unique(&mut a_gb_map, v);
            }
        }

        // `cxx:3827-3951`
        for esd in &l_esd {
            let i_esd = hds.ds().shape_index(esd, true);
            let rk_esd = hds.ds().ancestor_rank_of(esd);
            if rk_esd == rk_se {
                continue;
            }
            if BRepTool::is_degenerated(&Edge(esd.clone())) {
                continue;
            }
            // `cxx:3841-3844`
            let Some(is_so) = curves_same_oriented(esd, &se) else {
                continue;
            };

            // `cxx:3848-3866` : the interferences of `SE` with `Esd`
            let tol_esd = BRepTool::edge_tolerance(&Edge(esd.clone())).max(tol_se);
            let li_all: InterferenceList = hds.ds().shape_interferences_at(i_se, true).to_vec();
            let mut li: InterferenceList = Vec::new();
            let mut li1: InterferenceList = Vec::new();
            fds_assign(&li_all, &mut li);

            fun_select_traunk_interference(&mut li, &mut li1);
            li1.clear();
            fun_select_trasha_interference(&mut li, ShapeType::Edge, &mut li1);
            li.clear();
            fun_select_itrasha_interference(&mut li1, i_esd, &mut li);
            li1.clear();
            fun_select_sk_interference(&mut li, Kind::Edge, &mut li1);
            li.clear();
            fun_select_si_interference(&mut li1, i_esd, &mut li);
            li1.clear();
            let ni = fun_select_gk_interference(&mut li, Kind::Vertex, &mut li1);
            li.clear();
            li.append(&mut li1);

            // `cxx:3868-3877`
            let (vsd0, vsd1) = crate::topo_tools::edge_vertices(&Edge(esd.clone()));
            for vsd in [vsd0, vsd1].into_iter().flatten() {
                let a_v = &vsd.0;
                // `cxx:3878-3882` : do not create interferences on GBound
                if contains(&a_gb_map, a_v) {
                    continue;
                }

                // `cxx:3884-3888`
                let mut ori = a_v.orientation();
                if !is_so {
                    ori = reversed_orientation(ori);
                }

                // `cxx:3890-3923`
                if ni != 0 {
                    let mut li1: InterferenceList = Vec::new();
                    fds_assign(&li, &mut li1);
                    let mut li2: InterferenceList = Vec::new();
                    if fun_select_traori_interference(&mut li1, ori, &mut li2) != 0 {
                        let mut a_vg_map: Vec<TopoShape> = Vec::new();
                        push_unique(&mut a_vg_map, a_v);
                        for v in hds.ds().shape_same_domain(a_v) {
                            push_unique(&mut a_vg_map, v);
                        }
                        let mut already = false;
                        for i1 in &li2 {
                            let gi = i1.read().expect("interference lock").geometry();
                            match hds.ds().shape(gi, true) {
                                Some(a_vg) => {
                                    if contains(&a_vg_map, &a_vg) {
                                        already = true;
                                        break;
                                    }
                                }
                                None => break, // OCCT `BDS.Shape(gi)` raises there
                            }
                        }
                        if already {
                            continue; // has already
                        }
                    }
                }

                // `cxx:3925-3948`
                let tol = BRepTool::vertex_tolerance(&Vertex(a_v.clone())).max(tol_esd);
                let par_esd = BRepTool::parameter_of_vertex(esd, a_v);
                let Some(par) = fun_tool_par_e_tol(esd, par_esd, &se, tol) else {
                    continue;
                };
                let mut a_t = Transition::from_orientation(ori);
                a_t.set_shape_before(ShapeType::Edge);
                a_t.set_shape_after(ShapeType::Edge);
                a_t.set_index(i_esd);
                let iv = hds.change_ds().add_shape_with_ancestor(a_v, rk_esd);
                let is_g_bound = false;
                let a_config = if is_so {
                    Config::SameOriented
                } else {
                    Config::DiffOriented
                };
                let ai = make_edge_vertex_interference(&a_t, i_esd, iv, is_g_bound, a_config, par);
                let _ = hds.store_interference_on_shape(&ai, &se);
            }
        }
    }
}

/// `NCollection_Map::Add` over `TopTools_ShapeMapHasher` — the hasher compares
/// `TShape` and `Location` only, which is [`crate::topo_tools::is_same`].
fn push_unique(acc: &mut Vec<TopoShape>, s: &TopoShape) {
    if !contains(acc, s) {
        acc.push(s.clone());
    }
}

/// `NCollection_Map::Contains` for the same hasher.
fn contains(acc: &[TopoShape], s: &TopoShape) -> bool {
    acc.iter().any(|a| crate::topo_tools::is_same(a, s))
}

/// `TopAbs::Reverse(Or)` (`TopAbs.hxx:95-104`).
fn reversed_orientation(o: Orientation) -> Orientation {
    o.reversed()
}
