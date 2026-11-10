//! `FUN_ds_completeforSE1`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:862-1003`.

use crate::abs::{Orientation, ShapeType};
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::export::{
    fds_data, fds_parameter, fds_si_is_gi_of_i_of_sba_of_t_of_i, fds_tdata,
    fun_select_train_tinterference,
};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_gi_interference, fun_select_itrasha_interference, fun_select_sk_interference,
    fun_select_traori_interference, fun_select_trasha_interference,
    make_epv_interference_with_support,
};
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;

use super::purge_for_e9::fun_ds_has_sdm_ancestor_faces;

/// `FUN_ds_completeforSE1(HDS)` (`cxx:862-1003`).
///
/// xpu160398: an `INTERNAL` transition on a section edge, `ISE` being
/// `(INTERNAL(F'), G, edge)` while `F'` has no interference of support `SE`,
/// gets an `EPV` interference on the face `F'` that does.
pub fn fun_ds_complete_for_se1(hds: &mut HDataStructure, connex: &Connexity) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:882-885`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let ise = hds.ds().shape_index(&se, true);
        let li: InterferenceList = hds.ds().shape_interferences(&se, true).to_vec();

        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'item: {
                // `cxx:889-895`
                let (k, _g, loi) = {
                    let Some(v) = tki.value() else {
                        break 'item;
                    };
                    (v.0, v.1, v.2.to_vec())
                };

                // `cxx:896-903` : `ISE = (INTERNAL(FACE), G, EDGE)`
                let loicopy = loi.clone();
                let mut li1: InterferenceList = Vec::new();
                if !fun_select_train_tinterference(&loicopy, &mut li1) {
                    break 'item;
                }
                // `cxx:904-909`
                let mut li2: InterferenceList = Vec::new();
                if fun_select_trasha_interference(&mut li1, ShapeType::Face, &mut li2) < 1 {
                    break 'item;
                }
                // `cxx:910-914`
                let mut li3: InterferenceList = Vec::new();
                if fun_select_sk_interference(&mut li2, Kind::Edge, &mut li3) < 1 {
                    break 'item;
                }

                // `cxx:916-930`
                let keep_i = li3.iter().any(|h| {
                    let intf = h.read().expect("interference lock");
                    fds_si_is_gi_of_i_of_sba_of_t_of_i(hds.ds(), ise, &intf)
                });
                if keep_i {
                    break 'item;
                }

                // `cxx:935-997`
                for h in &li3 {
                    let (s, _isb) = {
                        let intf = h.read().expect("interference lock");
                        let (_gt, _g1, _st, s) = fds_data(&intf);
                        let (_tsb, isb, _tsa, _isa) = fds_tdata(&intf.transition());
                        (s, isb)
                    };
                    let Some(es) = hds.ds().shape(s, true) else {
                        continue;
                    };
                    // `cxx:949-952`
                    let Some((_fse, fes)) = fun_ds_has_sdm_ancestor_faces(hds, &se, &es, connex)
                    else {
                        continue;
                    };
                    let ifes = hds.ds().shape_index(&fes, true);

                    // `cxx:954-965`
                    let lif: InterferenceList =
                        hds.ds().shape_interferences(&fes, true).to_vec();
                    let mut lifcopy = lif.clone();
                    let mut lifound: InterferenceList = Vec::new();
                    if fun_select_gi_interference(&mut lifcopy, ise, &mut lifound) < 1 {
                        continue;
                    }

                    // `cxx:968-981`
                    let mut loicopy = loi.clone();
                    let mut li4: InterferenceList = Vec::new();
                    fun_select_itrasha_interference(&mut loicopy, ifes, &mut li4);
                    let mut lfor: InterferenceList = Vec::new();
                    let nfor =
                        fun_select_traori_interference(&mut li4, Orientation::Forward, &mut lfor);
                    let mut lrev: InterferenceList = Vec::new();
                    let nrev =
                        fun_select_traori_interference(&mut li4, Orientation::Reversed, &mut lrev);
                    if nfor > 0 || nrev > 0 {
                        break;
                    }

                    // `cxx:983-998`
                    let (par, g1, b) = {
                        let intf = h.read().expect("interference lock");
                        let Some(par) = fds_parameter(&intf) else {
                            continue;
                        };
                        let (_gt, g1, _st, _s) = fds_data(&intf);
                        let b = if intf.is_edge_vertex() {
                            intf.g_bound().unwrap_or(false)
                        } else {
                            false
                        };
                        (par, g1, b)
                    };
                    let mut new_t = Transition::from_orientation(Orientation::Internal);
                    new_t.set_index(ifes);
                    if let Some(new_i) =
                        make_epv_interference_with_support(&new_t, s, g1, par, k, Kind::Edge, b)
                    {
                        let _ = hds.store_interference_on_shape(&new_i, &se);
                    }
                    break;
                }
            }
            tki.next();
        }
    }
}
