//! `FUN_ds_PointToVertex`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:2868-3039`.

use std::collections::BTreeMap;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Vertex};
use crate::topopebrep::ds::export::{fds_copy, fds_data};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceData;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::make_epv_interference_with_support;
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::tool::topology::fun_tool_par_von_e;

/// `FUN_ds_PointToVertex(HDS)` (`cxx:2868-3039`).
///
/// xpu090698 / CTS20339: the intersector can report a `POINT` intersection
/// where the answer is a `VERTEX`. When an edge carries `EPI(T1, G, ea1)` and
/// `EPI(T2, G, ea2)` with the same support rank, `G` must be the vertex shared
/// by `ea1`/`ea2` -- two edges of the same matter cannot intersect.
///
/// Pass 1 (`cxx:2883-2957`) collects `G -> vG` for every such faulty point;
/// pass 2 (`cxx:2963-3033`) rewrites the interferences to that vertex and the
/// point geometries are removed (`cxx:3035-3038`).
pub fn fun_ds_point_to_vertex(hds: &mut HDataStructure) {
    let ns = hds.ds().nb_shapes();
    let mut i_pi_v: BTreeMap<i32, i32> = BTreeMap::new();

    // `cxx:2883-2957`
    for i in 1..=ns {
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != ShapeType::Edge {
            continue;
        }
        let li: InterferenceList = hds.ds().shape_interferences(&s, true).to_vec();
        if li.is_empty() {
            continue;
        }

        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'geom: {
                // `cxx:2900-2905`
                let Some((k, g, loi)) = tki.value().map(|v| (v.0, v.1, v.2.to_vec())) else {
                    break 'geom;
                };
                if k == Kind::Vertex {
                    break 'geom;
                }

                // `cxx:2907-2930`
                let mut scur = 0;
                let mut gfaulty = false;
                for intf in &loi {
                    let h = intf.read().expect("interference lock");
                    let (_gt, _g1, st, s2) = fds_data(&h);
                    drop(h);
                    if st != Kind::Edge {
                        continue;
                    }
                    if scur == 0 {
                        scur = s2;
                    }
                    if s2 != scur {
                        gfaulty = true;
                        break;
                    }
                }
                if !gfaulty {
                    break 'geom;
                }

                // `cxx:2932-2956` : in DS, pG --> vG
                let Some(pg) = hds.ds().point(g).map(|p| p.point().clone()) else {
                    break 'geom;
                };
                let rk_es = hds.ds().ancestor_rank(scur);
                let Some(es) = hds.ds().shape(scur, true) else {
                    break 'geom;
                };
                let (vf, vl) = crate::topo_tools::edge_vertices(&Edge(es));
                let (Some(vf), Some(vl)) = (vf, vl) else {
                    break 'geom;
                };
                let pf = BRepTool::vertex_point_world(&Vertex(vf.0.clone()));
                let pl = BRepTool::vertex_point_world(&Vertex(vl.0.clone()));
                let df = pf.distance(&pg);
                let dl = pl.distance(&pg);
                let vg = if df < dl { vf } else { vl };
                let ivg = hds.change_ds().add_shape_with_ancestor(&vg.0, rk_es);
                i_pi_v.insert(g, ivg);
            }
            tki.next();
        }
    }

    // `cxx:2959-2962`
    if i_pi_v.is_empty() {
        return;
    }

    // `cxx:2963-3033`
    for i in 1..=ns {
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != ShapeType::Edge {
            continue;
        }
        let li: InterferenceList = hds.ds().shape_interferences(&s, true).to_vec();
        if li.is_empty() {
            continue;
        }

        let mut new_li: InterferenceList = Vec::new();
        let rks = hds.ds().ancestor_rank_of(&s);
        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'geom: {
                // `cxx:2982-2990`
                let Some((_k, g, loi)) = tki.value().map(|v| (v.0, v.1, v.2.to_vec())) else {
                    break 'geom;
                };
                let Some(&ivg) = i_pi_v.get(&g) else {
                    fds_copy(&loi, &mut new_li);
                    break 'geom;
                };

                let Some(vg) = hds.ds().shape(ivg, true) else {
                    break 'geom;
                };
                let rk_g = hds.ds().ancestor_rank(ivg);
                let gbound = rk_g == rks;

                // `cxx:2997-3026`
                for intf in &loi {
                    let (par, t, st, s2) = {
                        let h = intf.read().expect("interference lock");
                        // `cxx:3004-3009` : only a `CurvePointInterference`
                        let InterferenceData::CurvePoint { param } = h.data() else {
                            continue;
                        };
                        let (_gt, _g1, st, s2) = fds_data(&h);
                        (*param, h.transition().clone(), st, s2)
                    };
                    let _ = par; // `cxx:3021` : `par` is always overwritten below
                    let Some(parvg) = fun_tool_par_von_e(&vg, &s) else {
                        continue; // `cxx:3016-3020`
                    };
                    if let Some(new_i) =
                        make_epv_interference_with_support(&t, s2, ivg, parvg, Kind::Vertex, st, gbound)
                    {
                        new_li.push(new_i);
                    }
                }
            }
            tki.next();
        }

        // `cxx:3028-3031`
        let lii = hds.change_ds().change_shape_interferences(&s);
        lii.clear();
        lii.append(&mut new_li);
    }

    // `cxx:3035-3038` : `BDS.RemovePoint(G)` for every rewritten geometry
    for g in i_pi_v.keys().copied().collect::<Vec<_>>() {
        hds.change_ds().remove_point(g);
    }
}
