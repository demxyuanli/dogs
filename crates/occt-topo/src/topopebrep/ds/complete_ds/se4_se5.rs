//! `FUN_ds_completeforSE4` and `FUN_ds_completeforSE5`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:2055-2333`.

use crate::abs::{Orientation, ShapeType};
use crate::topopebrep::ds::export::{fds_data, fds_parameter, fds_tdata};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_sk_interference, fun_select_traori_interference, fun_select_trasha_interference,
    make_epv_interference,
};
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::closed::is_closing_e_in_shape;
use crate::topopebrep::tool::curvature::get_xx;
use crate::topopebrep::tool::project_curve::fun_tool_par_e;
use crate::topopebrep::tool::topology::{fun_tool_ng_geom_f, fun_tool_par_on_e_f};

/// `FUN_ds_completeforSE4(HDS)` (`cxx:2055-2153`).
///
/// xpu160698: a section edge `SE` meets `ES` at the point `G`, `ES` being a
/// closing edge on `FTRA`; the two opposite transitions
/// `T1 = (OU/IN(FTRA), GP, ES)` and `T2 = (IN/OU(FTRA), GP, ES)` needed by the
/// `splitON` computation are rebuilt.
pub fn fun_ds_complete_for_se4(hds: &mut HDataStructure) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:2067-2070`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let li: InterferenceList = hds.ds().shape_interferences(&se, true).to_vec();

        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'item: {
                // `cxx:2076-2083`
                let (k, _g, loi) = {
                    let Some(v) = tki.value() else {
                        break 'item;
                    };
                    (v.0, v.1, v.2.to_vec())
                };
                if k != Kind::Point {
                    break 'item;
                }

                // `cxx:2085-2093`
                let mut loicopy = loi.clone();
                let mut l1: InterferenceList = Vec::new();
                fun_select_sk_interference(&mut loicopy, Kind::Edge, &mut l1);
                let mut l2: InterferenceList = Vec::new();
                if fun_select_trasha_interference(&mut l1, ShapeType::Face, &mut l2) < 1 {
                    break 'item;
                }

                // `cxx:2095-2110`
                let (par, g1, s, isb) = {                    let intf = l2[0].read().expect("interference lock");
                    let (_gt, g1, _st, s) = fds_data(&intf);
                    let (_tsb, isb, _tsa, _isa) = fds_tdata(&intf.transition());
                    let Some(par) = fds_parameter(&intf) else {
                        break 'item;
                    };
                    (par, g1, s, isb)
                };
                let Some(es) = hds.ds().shape(s, true) else {
                    break 'item;
                };
                let Some(ftra) = hds.ds().shape(isb, true) else {
                    break 'item;
                };
                if !is_closing_e_in_shape(&es, &ftra, &ftra) {
                    break 'item;
                }

                // `cxx:2112-2138`
                let mut has_for = false;
                let mut has_rev = false;
                for h in &l2 {
                    let intf = h.read().expect("interference lock");
                    let (_gt2, _g2, _st2, s2) = fds_data(&intf);
                    let (_tsb2, isb2, _tsa2, _isa2) = fds_tdata(&intf.transition());
                    if s2 != s || isb2 != isb {
                        // nyi raise
                        return;
                    }
                    let o2 = intf
                        .transition()
                        .orientation(crate::abs::State::In, ShapeType::Face);
                    if !has_for {
                        has_for = o2 == Orientation::Forward;
                    }
                    if !has_rev {
                        has_rev = o2 == Orientation::Reversed;
                    }
                }
                if !has_for && !has_rev {
                    break 'item;
                }
                if has_for && has_rev {
                    break 'item;
                }

                // `cxx:2145-2151`
                let new_o = if has_for {
                    Orientation::Reversed
                } else {
                    Orientation::Forward
                };
                let mut new_t = Transition::from_orientation(new_o);
                new_t.set_index(isb);
                if let Some(new_i) = make_epv_interference(&new_t, s, g1, par, k, false) {
                    let _ = hds.store_interference_on_shape(&new_i, &se);
                }
            }
            tki.next();
        }
    }
}

/// `FUN_ds_completeforSE5(HDS)` (`cxx:2157-2333`).
///
/// xpu190698: `SE` has `I = (T(F), G, S)` with `T = (OUT,OUT)` (or `(IN,IN)`)
/// and `I' = (T'(Fsd), G, Sedge)` with `T'` `FORWARD`/`REVERSED`; `FCX` is
/// tangent to `F` at `G`. The 2d interferences are re-oriented so that the
/// `spOU`/`spIN` of `ES` can be computed.
pub fn fun_ds_complete_for_se5(hds: &mut HDataStructure) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:2176-2191`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let loi_all: InterferenceList = hds.ds().shape_interferences(&se, true).to_vec();
        let mut loic = loi_all.clone();
        let mut li: InterferenceList = Vec::new();
        if fun_select_sk_interference(&mut loic, Kind::Edge, &mut li) < 1 {
            continue;
        }

        let mut new_li: InterferenceList = Vec::new();
        let mut hasnew_li = false;
        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            'item: {
                // `cxx:2197-2204`
                let (k, _g, loi) = {
                    let Some(v) = tki.value() else {
                        break 'item;
                    };
                    (v.0, v.1, v.2.to_vec())
                };
                if k != Kind::Point {
                    new_li.extend(loi);
                    break 'item;
                }

                // `cxx:2206-2222`
                let mut loicopy = loi.clone();
                let mut lext: InterferenceList = Vec::new();
                let iext =
                    fun_select_traori_interference(&mut loicopy, Orientation::External, &mut lext);
                let hasext = iext > 0;
                let mut lint: InterferenceList = Vec::new();
                let iint =
                    fun_select_traori_interference(&mut loicopy, Orientation::Internal, &mut lint);
                let hasint = iint > 0;
                if !hasext && !hasint {
                    new_li.extend(loi);
                    break 'item;
                }
                if hasext && hasint {
                    // nyiFUN_RAISE : incoherent data
                    new_li.extend(loi);
                    break 'item;
                }

                // `cxx:2227-2243` : NYI: only select among 2dI
                let mut lfor: InterferenceList = Vec::new();
                let ifor =
                    fun_select_traori_interference(&mut loicopy, Orientation::Forward, &mut lfor);
                let mut lrev: InterferenceList = Vec::new();
                let irev =
                    fun_select_traori_interference(&mut loicopy, Orientation::Reversed, &mut lrev);
                let hasrev = irev > 0;
                let hasfor = ifor > 0;
                if !hasrev && !hasfor {
                    new_li.extend(loi);
                    break 'item;
                }
                if hasrev && hasfor {
                    new_li.extend(loi);
                    break 'item;
                }

                // `cxx:2246-2256`
                let i2d = if hasfor { lfor[0].clone() } else { lrev[0].clone() };
                let (s2, if2, par_se) = {
                    let intf = i2d.read().expect("interference lock");
                    let s2 = intf.support();
                    let Ok(if2) = intf.transition().index() else {
                        return;
                    };
                    let Some(par_se) = fds_parameter(&intf) else {
                        return;
                    };
                    (s2, if2, par_se)
                };
                let Some(e2) = hds.ds().shape(s2, true) else {
                    return;
                };
                let Some(f2) = hds.ds().shape(if2, true) else {
                    return;
                };
                let Some(par_e2) = fun_tool_par_e(&se, par_se, &e2) else {
                    return;
                };
                let Some(uv2) = fun_tool_par_on_e_f(&e2, par_e2, &f2) else {
                    return;
                };
                let Some(ng_f2) = fun_tool_ng_geom_f(&uv2, &crate::shape::Face(f2.clone())) else {
                    return;
                };
                let Some(xx_f2) = get_xx(&e2, &f2, par_e2, &ng_f2) else {
                    return;
                };

                // `cxx:2272-2288`
                let i3d = if hasext { lext[0].clone() } else { lint[0].clone() };
                let if3 = {
                    let intf = i3d.read().expect("interference lock");
                    let Ok(if3) = intf.transition().index() else {
                        return;
                    };
                    if3
                };
                let Some(f3) = hds.ds().shape(if3, true) else {
                    return;
                };
                let Some(uv3) = fun_tool_par_on_e_f(&e2, par_e2, &f3) else {
                    return;
                };
                let Some(ng_f3) = fun_tool_ng_geom_f(&uv3, &crate::shape::Face(f3.clone())) else {
                    return;
                };
                let Some(xx_f3) = get_xx(&e2, &f3, par_e2, &ng_f3) else {
                    return;
                };

                // `cxx:2290-2319`
                let dot = xx_f2.dot(&xx_f3);
                let positive = dot > 0.0;
                if positive {
                    break 'item;
                }

                let mut new_o = if hasfor {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                if hasint {
                    new_o = new_o.complement();
                }

                let group = if hasint { lint } else { lext };
                for h in &group {
                    h.write()
                        .expect("interference lock")
                        .change_transition()
                        .set_orientation(new_o);
                    new_li.push(h.clone());
                    hasnew_li = true;
                }
                new_li.extend(lrev);
                new_li.extend(lfor);
                new_li.extend(loicopy);
            }
            tki.next();
        }

        // `cxx:2321-2331`
        if hasnew_li {
            new_li.extend(loic);
            let lii = hds.change_ds().change_shape_interferences(&se);
            lii.clear();
            lii.append(&mut new_li);
        }
    }
}
