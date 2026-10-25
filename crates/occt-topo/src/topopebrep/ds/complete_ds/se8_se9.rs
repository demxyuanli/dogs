//! `FUN_ds_completeforSE8` and `FUN_ds_completeforSE9`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:2684-2866`.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::fun_select_sk_interference;
use crate::topopebrep::ds::process_interferences_tool::fun_select_si_interference;
use crate::topopebrep::ds::process_interferences_tool::fun_select_trasha_interference;
use crate::topopebrep::ds::process_interferences_tool::make_epv_interference_with_support;
use crate::topopebrep::ds::redu::fun_ds_redu2d1d;
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::topology::fun_tool_orient_vin_e;

use crate::topopebrep::ds::process_interferences_tool::fun_select_pure2d_i;

/// `FUN_ds_completeforSE8(HDS)` (`cxx:2684-2778`).
///
/// xpu020998: attached to a section edge `SE`,
/// `{I1d = (Tr(Esd), vG, Esd), I2d = (Tr(F), vG, E)}` with `E` and `Esd` edges
/// of `F`; `I2d` may be incomplete, so the 2d/1d pair is reduced to
/// `{I1d, newI2d}`.
pub fn fun_ds_complete_for_se8(hds: &mut HDataStructure) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:2699-2704`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let ise = hds.ds().shape_index(&se, true);
        let li: InterferenceList = hds.ds().shape_interferences(&se, true).to_vec();

        // `cxx:2704` : `hasnew` is declared `false` and never assigned by OCCT,
        // so the `if (!hasnew)` rebuild of `cxx:2760-2775` never runs there.
        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'item: {
                // `cxx:2708-2712`
                let (kk, gg, loi) = {
                    let Some(v) = tki.value() else {
                        break 'item;
                    };
                    (v.0, v.1, v.2.to_vec())
                };
                let isvertex = kk == Kind::Vertex;
                if isvertex {
                    let Some(gshape) = hds.ds().shape(gg, true) else {
                        break 'item;
                    };
                    if hds.has_same_domain(&gshape, true) {
                        break 'item;
                    }
                }

                // `cxx:2721-2730` : li -> l1dE + l2dFE + lFE(<=>l3dFE) + li(<=>lFF)
                let mut li = loi.clone();
                let mut l1de: InterferenceList = Vec::new();
                let n1d = fun_select_trasha_interference(&mut li, ShapeType::Edge, &mut l1de);
                let mut lfe: InterferenceList = Vec::new();
                fun_select_sk_interference(&mut li, Kind::Edge, &mut lfe);
                let mut l2dfe: InterferenceList = Vec::new();
                let n2d = fun_select_pure2d_i(&li, &mut lfe, &mut l2dfe);

                // `cxx:2731-2751`
                let redu2d = n1d > 0 && n2d > 0;
                if redu2d {
                    let mut al2dfe: InterferenceList = Vec::new();
                    let mut idx = 0usize;
                    while idx < l2dfe.len() {
                        let new_t2d = fun_ds_redu2d1d(hds.ds(), ise, &l2dfe[idx], &l1de);
                        let Some(new_t2d) = new_t2d else {
                            idx += 1;
                            continue;
                        };
                        l2dfe[idx]
                            .write()
                            .expect("interference lock")
                            .set_transition(new_t2d);
                        let removed = l2dfe.remove(idx);
                        al2dfe.push(removed);
                    }
                    l2dfe.extend(al2dfe);
                }

                // `cxx:2753-2758`
                if let Some(loii) = tki.change_interferences(kk, gg) {
                    loii.clear();
                    loii.extend(l1de);
                    loii.extend(lfe);
                    loii.extend(l2dfe);
                    loii.extend(li);
                }
            }
            tki.next();
        }
    }
}

/// `FUN_ds_completeforSE9(HDS)` (`cxx:2779-2866`).
///
/// xpu011098: `ES` is a section edge same-domain with `Esd`, and carries no
/// interference; the interferences of `ES` are rebuilt from `Esd`.
pub fn fun_ds_complete_for_se9(hds: &mut HDataStructure) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:2790-2801`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let rkse = hds.ds().ancestor_rank_of(&se);
        let ise = hds.ds().shape_index(&se, true);
        if !hds.has_same_domain(&se, true) {
            continue;
        }
        if !hds.ds().shape_interferences(&se, true).is_empty() {
            continue;
        }

        // `cxx:2804-2864`
        let esdse = hds.ds().shape_same_domain(&se).to_vec();
        for esd in esdse {
            let i_esd = hds.ds().shape_index(&esd, true);
            let rk_esd = hds.ds().ancestor_rank_of(&esd);
            if rk_esd == rkse {
                continue;
            }
            let li: InterferenceList = hds.ds().shape_interferences(&esd, true).to_vec();
            if li.is_empty() {
                continue;
            }

            let mut lise: InterferenceList = Vec::new();
            fun_select_si_interference(&mut li.clone(), ise, &mut lise);

            let mut tki = TKI::new();
            tki.fill_on_geometry(&lise);
            tki.init();
            while tki.more() {
                'item: {
                    // `cxx:2829-2843`
                    let (_k, g) = {
                        let Some(v) = tki.value() else {
                            break 'item;
                        };
                        (v.0, v.1)
                    };
                    let Some(vg) = hds.ds().shape(g, true) else {
                        break 'item;
                    };
                    if hds.has_same_domain(&vg, true) {
                        break 'item; // nyixpu011098
                    }
                    if hds.ds().ancestor_rank(g) != rkse {
                        break 'item; // nyixpu011098
                    }

                    // `cxx:2847-2862`
                    let mut new_t = Transition::new(
                        crate::abs::State::In,
                        crate::abs::State::In,
                        ShapeType::Edge,
                        ShapeType::Edge,
                    );
                    new_t.set_index(i_esd);
                    if fun_tool_orient_vin_e(&vg, &se) == 0 {
                        break 'item;
                    }
                    let par = BRepTool::parameter_of_vertex(&se, &vg);
                    if let Some(new_i) = make_epv_interference_with_support(
                        &new_t,
                        i_esd,
                        g,
                        par,
                        Kind::Vertex,
                        Kind::Edge,
                        true,
                    ) {
                        let _ = hds.store_interference_on_shape(&new_i, &se);
                    }
                }
                tki.next();
            }
        }
    }
}
