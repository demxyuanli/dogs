//! `FUN_ds_completeforSE6` and `FUN_ds_completeforE7`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:2336-2682`.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topopebrep::ds::export::{
    fds_has_same_domain3d, fds_loi_inf_sup, fds_parameter, fun_ds_getoov_in_hds,
};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_gi_interference, fun_select_gk_interference, fun_select_itrasha_interference,
    fun_select_si_interference, fun_select_sk_interference, fun_select_traori_interference,
    fun_select_trasha_interference, make_epv_interference_with_support,
};
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::make_transition::MakeTransition;
use crate::topopebrep::tool::project_curve::fun_tool_par_e;
use crate::topopebrep::tool::topology::{fun_tool_par_f, fun_tool_par_on_e_f};

use super::has_i2d::fun_ds_has_i2d;
use super::share_g::{fun_ds_ori_e_in_f, ORI_E_IN_F_UNKNOWN};

/// `FUN_ds_completeforSE6(HDS)` (`cxx:2336-2503`).
///
/// xpu280798: a section edge `E` same-domain 3d with `Esd`, `E` carrying
/// `I2d(T(F), VG, E')` where `VG` is a non same-domain vertex of `Esd` and
/// `Esd` is an edge of `F`; the 1d interference `newT(Esd), VG, Esd` is added.
pub fn fun_ds_complete_for_se6(hds: &mut HDataStructure) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:2348-2357`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let _ise = hds.ds().shape_index(&se, true);
        let (hassd3d, l_esd3d) = fds_has_same_domain3d(hds.ds(), &se, true);
        if !hassd3d {
            continue;
        }

        // `cxx:2359-2381`
        let li: InterferenceList = hds.ds().shape_interferences(&se, true).to_vec();
        let mut lia: InterferenceList = Vec::new();
        if fun_select_gk_interference(&mut li.clone(), Kind::Vertex, &mut lia) == 0 {
            continue;
        }
        let mut lib: InterferenceList = Vec::new();
        for ia in &lia {
            let g = ia.read().expect("interference lock").geometry();
            let hassd = match hds.ds().shape(g, true) {
                Some(vg) => fun_ds_getoov_in_hds(hds, &vg).is_some(),
                None => false,
            };
            if !hassd {
                lib.push(ia.clone());
            }
        }

        // `cxx:2382-2387`
        let l2dfe = fun_ds_has_i2d(&lib);
        let mut l1de: InterferenceList = Vec::new();
        fun_select_trasha_interference(&mut lib.clone(), ShapeType::Edge, &mut l1de);

        // `cxx:2390-2499`
        for esd in &l_esd3d {
            let (vf, vl) = crate::topo_tools::edge_vertices(&Edge(esd.clone()));
            let vf = vf.map(|v| v.0);
            let vl = vl.map(|v| v.0);
            // `cxx:2397-2399`
            if BRepTool::is_degenerated(&Edge(esd.clone())) {
                continue;
            }
            let closed = match (&vf, &vl) {
                (Some(a), Some(b)) => crate::topo_tools::is_same(a, b),
                _ => false,
            };
            let i_esd = hds.ds().shape_index(esd, true);
            let ivf = vf
                .as_ref()
                .map(|v| hds.ds().shape_index(v, true))
                .unwrap_or(0);
            let ivl = vl
                .as_ref()
                .map(|v| hds.ds().shape_index(v, true))
                .unwrap_or(0);

            for iv in 1..=2 {
                let g = if iv == 1 { ivf } else { ivl };
                if g == 0 {
                    continue;
                }
                let Some(vg) = hds.ds().shape(g, true) else {
                    continue;
                };
                if fun_ds_getoov_in_hds(hds, &vg).is_some() {
                    continue;
                }

                // `cxx:2422-2436`
                let mut l1dg: InterferenceList = Vec::new();
                fun_select_gi_interference(&mut l1de.clone(), g, &mut l1dg);
                let mut l2dg: InterferenceList = Vec::new();
                if fun_select_gi_interference(&mut l2dfe.clone(), g, &mut l2dg) == 0 {
                    // no 2dI at G
                    continue;
                }
                let mut l1dgesd: InterferenceList = Vec::new();
                if fun_select_itrasha_interference(&mut l1dg.clone(), i_esd, &mut l1dgesd) != 0 {
                    // 1dI(Esd) at G exists already
                    continue;
                }

                // `cxx:2438-2497`
                for i2d in &l2dg {
                    let (itrasha, o, par_e) = {
                        let intf = i2d.read().expect("interference lock");
                        let itrasha = intf.transition().index_before();
                        let o = intf
                            .transition()
                            .orientation(crate::abs::State::In, ShapeType::Face);
                        let Some(par_e) = fds_parameter(&intf) else {
                            continue;
                        };
                        (itrasha, o, par_e)
                    };
                    let Some(f) = hds.ds().shape(itrasha, true) else {
                        continue;
                    };
                    let (esd_of_f, _dum) = fun_ds_ori_e_in_f(hds.ds(), esd, &f);
                    if esd_of_f == ORI_E_IN_F_UNKNOWN {
                        continue;
                    }

                    // `cxx:2457-2493`
                    let mut new_t = Transition::new(
                        crate::abs::State::Out,
                        crate::abs::State::Out,
                        ShapeType::Edge,
                        ShapeType::Edge,
                    );
                    if closed {
                        new_t.set_orientation(Orientation::Internal);
                    } else if o == Orientation::Forward || o == Orientation::Reversed {
                        new_t.set_orientation(o);
                    } else if o == Orientation::External {
                        new_t.set_orientation(o);
                    } else if o == Orientation::Internal {
                        let par_esd = BRepTool::parameter_of_vertex(esd, &vg);
                        let tg_esd = crate::topopebrep::tool::tangent::tg_geom_e(par_esd, esd);
                        let tg_e = crate::topopebrep::tool::tangent::tg_geom_e(par_e, &se);
                        let (Some(tg_esd), Some(tg_e)) = (tg_esd, tg_e) else {
                            continue;
                        };
                        let dot = tg_esd.dot(&tg_e);
                        let so = dot > 0.0;
                        let isvf = iv == 1;
                        let isforw = (so && isvf) || (!so && !isvf);
                        new_t.set_orientation(if isforw {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        });
                    }
                    new_t.set_index(i_esd);

                    // `cxx:2493-2496`
                    if let Some(new_i1d) = make_epv_interference_with_support(
                        &new_t,
                        i_esd,
                        g,
                        par_e,
                        Kind::Vertex,
                        Kind::Edge,
                        false,
                    ) {
                        let _ = hds.store_interference_on_shape(&new_i1d, &se);
                    }
                }
            }
        }
    }
}

/// `FUN_ds_completeforE7(HDS)` (`cxx:2504-2682`).
///
/// xpu130898: on an edge `E`, `lI3d = {I = (FOR/REV(F), G, F)}`; when the pair
/// describes an `INTERNAL`/`EXTERNAL` transition the 3d interferences are
/// reduced to the `MkTonE` answer.
pub fn fun_ds_complete_for_e7(hds: &mut HDataStructure) {
    let ns = hds.ds().nb_shapes();
    for i in 1..=ns {
        // `cxx:2512-2528`
        let Some(ee) = hds.ds().shape(i, true) else {
            continue;
        };
        if ee.shape_type() != ShapeType::Edge {
            continue;
        }
        if BRepTool::is_degenerated(&Edge(ee.clone())) {
            continue;
        }
        let li: InterferenceList = hds.ds().shape_interferences(&ee, true).to_vec();

        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'item: {
                // `cxx:2534-2538`
                let (k, g, loi) = {
                    let Some(v) = tki.value() else {
                        break 'item;
                    };
                    (v.0, v.1, v.2.to_vec())
                };
                let mut loicopy = loi.clone();
                let mut l1: InterferenceList = Vec::new();
                fun_select_sk_interference(&mut loicopy, Kind::Face, &mut l1);
                let mut lfor: InterferenceList = Vec::new();
                let nfor =
                    fun_select_traori_interference(&mut l1, Orientation::Forward, &mut lfor);
                let mut lrev: InterferenceList = Vec::new();
                let nrev =
                    fun_select_traori_interference(&mut l1, Orientation::Reversed, &mut lrev);
                if nfor == 0 || nrev == 0 {
                    break 'item;
                }

                // `cxx:2551-2669`
                let mut lnewi: InterferenceList = Vec::new();
                let mut ifs = 0i32;
                for ifor in &lfor {
                    let (ifs_i, par, b) = {
                        let intf = ifor.read().expect("interference lock");
                        let ifs_i = intf.support();
                        let Some(par) = fds_parameter(&intf) else {
                            continue;
                        };
                        let b = if k == Kind::Vertex {
                            intf.g_bound().unwrap_or(false)
                        } else {
                            false
                        };
                        (ifs_i, par, b)
                    };
                    let Some(fs) = hds.ds().shape(ifs_i, true) else {
                        continue;
                    };
                    let mut lfs: InterferenceList = Vec::new();
                    if fun_select_si_interference(&mut loicopy.clone(), ifs_i, &mut lrev)
                        == 0
                    {
                        continue;
                    }
                    let mut lfse: InterferenceList = Vec::new();
                    let nfse =
                        fun_select_itrasha_interference(&mut loicopy.clone(), ifs_i, &mut lfse);

                    let (par1, par2, _isonper) = {
                        let loi_e: InterferenceList =
                            hds.ds().shape_interferences(&ee, true).to_vec();
                        match fds_loi_inf_sup(&ee, par, k, g, &loi_e) {
                            Ok(v) => v,
                            Err(_) => return,
                        }
                    };
                    let factor = 1.0e-4;
                    let mut new_t = Transition::default();
                    let ies;

                    if nfse == 0 {
                        // `cxx:2579-2599`
                        let Some(uv_fs) = fun_tool_par_f(&ee, par, &fs) else {
                            continue;
                        };
                        let mut mkt = MakeTransition::new();
                        let mut stb = crate::abs::State::Unknown;
                        let mut sta = crate::abs::State::Unknown;
                        let mut ok = mkt.initialize(&ee, par1, par2, par, &fs, &uv_fs, factor);
                        if ok {
                            ok = mkt.mk_ton_e(&mut stb, &mut sta);
                        }
                        if !ok {
                            continue;
                        }
                        new_t.set_state_before(stb);
                        new_t.set_state_after(sta);
                        ies = 0;
                    } else {
                        // `cxx:2600-2637`
                        let ies_v = lfse[0].read().expect("interference lock").support();
                        let Some(es) = hds.ds().shape(ies_v, true) else {
                            continue;
                        };
                        let Some(par_es) = fun_tool_par_e(&ee, par, &es) else {
                            continue;
                        };
                        let Some(uv_fs) = fun_tool_par_on_e_f(&es, par_es, &fs) else {
                            continue;
                        };
                        let mut mkt = MakeTransition::new();
                        let mut stb = crate::abs::State::Unknown;
                        let mut sta = crate::abs::State::Unknown;
                        let mut ok = mkt.initialize(&ee, par1, par2, par, &fs, &uv_fs, factor);
                        if ok {
                            ok = mkt.set_rest(&es, par_es);
                        }
                        if ok {
                            ok = mkt.mk_ton_e(&mut stb, &mut sta);
                        }
                        if !ok {
                            continue;
                        }
                        new_t.set_state_before(stb);
                        new_t.set_state_after(sta);
                        ies = ies_v;
                    }

                    // `cxx:2638-2655`
                    ifs = ifs_i;
                    new_t.set_index(ifs);
                    if let Some(new_i) = make_epv_interference_with_support(
                        &new_t, ifs, g, par, k, Kind::Face, b,
                    ) {
                        lnewi.push(new_i);
                    }
                    if nfse != 0 {
                        if let Some(new_ife) = make_epv_interference_with_support(
                            &new_t, ies, g, par, k, Kind::Edge, b,
                        ) {
                            lnewi.push(new_ife);
                        }
                    }
                }

                // `cxx:2657-2668`
                if ifs != 0 {
                    if let Some(loii) = tki.change_interferences(k, g) {
                        let mut l_efs: InterferenceList = Vec::new();
                        fun_select_itrasha_interference(&mut loii.clone(), ifs, &mut l_efs);
                        for new_i in &lnewi {
                            loii.push(new_i.clone());
                        }
                    }
                }
            }
            tki.next();
        }

        // `cxx:2671-2680`
        let lii = hds.change_ds().change_shape_interferences(&ee);
        lii.clear();
        tki.init();
        while tki.more() {
            if let Some((_k, _g, l)) = tki.value() {
                lii.extend_from_slice(l);
            }
            tki.next();
        }
    }
}
