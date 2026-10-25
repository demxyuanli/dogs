//! `FUN_ds_completeforSE2`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:1005-1223`.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::Vertex;
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::export::{
    fds_data, fds_loi_inf_sup, fds_parameter, fds_tdata, fun_ds_getoov,
    fun_ds_getoov_in_hds,
};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_itrasha_interference, fun_select_sk_interference, fun_select_trasha_interference,
    make_epv_interference_with_support,
};
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::make_transition::MakeTransition;
use crate::topopebrep::tool::measures::fun_tool_maxtol;
use crate::topopebrep::tool::project_curve::fun_tool_proj_pon_e;
use crate::topopebrep::tool::topology::{fun_tool_orient_vin_e, fun_tool_par_on_e_f};

/// `FUN_ds_completeforSE2(HDS)` (`cxx:1005-1223`).
///
/// xpu250398: to have `spOUT(SE)` / `spIN(SE)` on the bound `G` of a section
/// edge, `IF = (T(F), G, F)` is needed with `G` not a bound of `SE`; when only
/// `IFE = (T(FTRA), G, ES)` exists, `I3d = (T2(FCX), G, FCX)` is rebuilt for
/// every face `FCX` sharing `ES` with `FTRA`.
pub fn fun_ds_complete_for_se2(hds: &mut HDataStructure, connex: &Connexity) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:1024-1029`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let rkse = hds.ds().ancestor_rank_of(&se);
        let li: InterferenceList = hds.ds().shape_interferences(&se, true).to_vec();

        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'item: {
                // `cxx:1033-1037`
                let (k, g, loi) = {
                    let Some(v) = tki.value() else {
                        break 'item;
                    };
                    (v.0, v.1, v.2.to_vec())
                };
                let point = k == Kind::Point;
                let vertex = k == Kind::Vertex;

                // `cxx:1039-1058`
                if vertex {
                    let Some(vg) = hds.ds().shape(g, true) else {
                        break 'item;
                    };
                    let oov = fun_ds_getoov_in_hds(hds, &vg);
                    let mut ovse = fun_tool_orient_vin_e(&vg, &se);
                    if ovse == 0 {
                        if let Some(oov) = &oov {
                            ovse = fun_tool_orient_vin_e(oov, &se);
                        }
                    }
                    if ovse != 0 {
                        break 'item;
                    }
                }

                // `cxx:1062-1079`
                let mut loicopy = loi.clone();
                let mut l1: InterferenceList = Vec::new();
                fun_select_trasha_interference(&mut loicopy, ShapeType::Face, &mut l1);
                let mut lf: InterferenceList = Vec::new();
                let nf = fun_select_sk_interference(&mut l1, Kind::Face, &mut lf);
                if nf > 1 {
                    break 'item;
                }
                let mut lfe: InterferenceList = Vec::new();
                let nfe = fun_select_sk_interference(&mut l1, Kind::Edge, &mut lfe);
                if nfe == 0 {
                    break 'item;
                }

                // `cxx:1083-1094` : `FDS_data` overwrites the TKI group `(K, G)`
                let (par, k, g, isb, s) = {
                    let intf = lfe[0].read().expect("interference lock");
                    let Some(par) = fds_parameter(&intf) else {
                        break 'item;
                    };
                    let (k, g, _st, s) = fds_data(&intf);
                    let (_tsb, isb, _tsa, _isa) = fds_tdata(&intf.transition());
                    (par, k, g, isb, s)
                };
                let Some(ftra) = hds.ds().shape(isb, true) else {
                    break 'item;
                };
                let Some(es) = hds.ds().shape(s, true) else {
                    break 'item;
                };

                // `cxx:1096-1220`
                let lfcx = connex.edge_connexity_same_shape(&es, hds).to_vec();
                for fcx in lfcx {
                    let ifcx = hds.ds().shape_index(&fcx, true);
                    if crate::topo_tools::is_same(&fcx, &ftra) {
                        continue;
                    }

                    // `cxx:1105-1111` : NYI check SE has part ON FCX(with bound G)
                    let mut l2: InterferenceList = Vec::new();
                    if fun_select_itrasha_interference(&mut lfe.clone(), ifcx, &mut l2) > 0 {
                        continue;
                    }

                    // `cxx:1112-1177`
                    let mut oopar = 0.0f64;
                    if point {
                        let Some(pt) = hds.ds().point(g) else {
                            // `BDS.Point(G)` raises
                            return;
                        };
                        let p3d = *pt.point();
                        let t = pt.tolerance().max(fun_tool_maxtol(&es));
                        let Some((par_es, d)) = fun_tool_proj_pon_e(&p3d, &es) else {
                            continue;
                        };
                        if d > t {
                            continue;
                        }
                        oopar = par_es;
                    }
                    if vertex {
                        let Some(vg) = hds.ds().shape(g, true) else {
                            continue;
                        };
                        let rkg = hds.ds().ancestor_rank(g);
                        if rkg == rkse {
                            if fun_tool_orient_vin_e(&vg, &es) == 0 {
                                continue;
                            }
                            oopar = BRepTool::parameter_of_vertex(&es, &vg);
                        } else if let Some(oov) = fun_ds_getoov(hds.ds(), &vg) {
                            oopar = BRepTool::parameter_of_vertex(&es, &oov);
                        } else {
                            let t = BRepTool::vertex_tolerance(&Vertex(vg.clone()))
                                .max(fun_tool_maxtol(&es));
                            let p = BRepTool::vertex_point(&Vertex(vg.clone()));
                            let Some((par_es, dd)) = fun_tool_proj_pon_e(&p, &es) else {
                                continue;
                            };
                            if dd > t {
                                continue;
                            }
                            oopar = par_es;
                        }
                    }

                    // `cxx:1179-1183`
                    let Some(oouv) = fun_tool_par_on_e_f(&es, oopar, &fcx) else {
                        continue;
                    };

                    // `cxx:1185-1205`
                    let (par1, par2, _isonper) = {
                        let loi_se: InterferenceList =
                            hds.ds().shape_interferences(&se, true).to_vec();
                        match fds_loi_inf_sup(&se, par, k, g, &loi_se) {
                            Ok(v) => v,
                            Err(_) => return,
                        }
                    };
                    let factor = 1.0e-4;
                    let mut mkt = MakeTransition::new();
                    let mut stb = crate::abs::State::Unknown;
                    let mut sta = crate::abs::State::Unknown;
                    let mut ok = mkt.initialize(&se, par1, par2, par, &fcx, &oouv, factor);
                    if ok {
                        ok = mkt.set_rest(&es, oopar);
                    }
                    if ok {
                        ok = mkt.mk_ton_e(&mut stb, &mut sta);
                    }
                    if !ok {
                        continue;
                    }
                    let mut new_t = Transition::default();
                    new_t.set_state_before(stb);
                    new_t.set_state_after(sta);
                    new_t.set_index(ifcx);

                    // `cxx:1212-1219`
                    let b = if vertex {
                        lfe[0]
                            .read()
                            .expect("interference lock")
                            .g_bound()
                            .unwrap_or(false)
                    } else {
                        false
                    };
                    if let Some(new_i) =
                        make_epv_interference_with_support(&new_t, ifcx, g, par, k, Kind::Face, b)
                    {
                        let _ = hds.store_interference_on_shape(&new_i, &se);
                    }
                }
            }
            tki.next();
        }
    }
}
