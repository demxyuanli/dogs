//! `FUN_ds_redusamsha` and its `FUN_redusamshaonE` helper.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:3041-3277`.

use crate::abs::{Orientation, ShapeType, State};
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::{fds_idata, fds_parameter, fds_state_e_with_f2d};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_sk_interference, fun_select_trasha_interference, make_epv_interference_with_support,
};
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::make_transition::MakeTransition;
use crate::topopebrep::tool::fun_tool_bounds;
use crate::topopebrep::tool::project_curve::fun_tool_par_e;
use crate::topopebrep::tool::topology::fun_tool_par_on_e_f;

/// `FUN_redusamshaonE(BDS, I, EIX, newI)` (`cxx:3041-3107`).
///
/// Attached to `Edge(EIX)`: `I = (FORWARD/REVERSED(ES), G, ES)`. The new
/// interference gets the 2d transition of `E` on `Face(IB)` recomputed by
/// `TopOpeBRepTool_makeTransition`.
///
/// `None` is the OCCT `false` return. OCCT's `FDS_Parameter` raises
/// `Standard_Failure` on any interference that is neither an
/// `EdgeVertexInterference` nor a `CurvePointInterference`; the port answers
/// `None` there, i.e. the caller's "no reduction" arm.
pub fn fun_redusamsha_on_e(
    bds: &DataStructure,
    i: &InterferenceHandle,
    eix: i32,
) -> Option<InterferenceHandle> {
    // `cxx:3049-3056`
    let (sb, ib, sa, ia, gt, g, st, s) = fds_idata(&i.read().expect("interference lock"));
    let _ = (sb, sa, ia, st);
    let e = bds.shape(eix, true)?;
    let par_e = fds_parameter(&i.read().expect("interference lock"))?;
    let (f, l) = fun_tool_bounds(&e);
    let es = bds.shape(s, true)?;
    let ftra = bds.shape(ib, true)?;

    // `cxx:3058-3071`
    let par_es = fun_tool_par_e(&e, par_e, &es)?;
    let uv = fun_tool_par_on_e_f(&es, par_es, &ftra)?;

    // `cxx:3073-3087`
    let factor = 1.0e-2;
    let mut stb = State::Unknown;
    let mut sta = State::Unknown;
    let mut mkt = MakeTransition::new();
    let mut ok = mkt.initialize(&e, f, l, par_e, &ftra, &uv, factor);
    if ok {
        ok = mkt.set_rest(&es, par_es);
    }
    if ok {
        ok = mkt.mk_ton_e(&mut stb, &mut sta);
    }
    if !ok {
        return None;
    }

    // `cxx:3088-3091`
    let mut new_t = Transition::default();
    new_t.set_index(ib);
    new_t.set_state_before(stb);
    new_t.set_state_after(sta);

    // `cxx:3093-3097`
    let ok = fds_state_e_with_f2d(bds, &e, par_e, gt, g, &ftra, &mut new_t).ok()?;
    if !ok {
        return None;
    }

    // `cxx:3099-3103`
    let mut b = false;
    if gt == Kind::Vertex {
        b = i
            .read()
            .expect("interference lock")
            .g_bound()
            .unwrap_or(false);
    }
    make_epv_interference_with_support(&new_t, s, g, par_e, gt, Kind::Edge, b)
}

/// `FUN_ds_redusamsha(HDS)` (`cxx:3109-3277`).
///
/// xpu201098: reduce 2d interferences on the same `G` and the same `S`:
/// `{I1=(OU/IN(F), G, ES), I2=(IN/OU(F), G, ES)}` collapse to the single
/// transition recomputed by [`fun_redusamsha_on_e`].
/// cto009K1 (`SE7`, `ES8`, `G2`, `FTRA27`), xpu050299 (`FRA60618` e13,
/// `ftrasha14`).
pub fn fun_ds_redusamsha(hds: &mut HDataStructure) {
    let ns = hds.ds().nb_shapes();

    // `cxx:3119-3275`
    for i in 1..=ns {
        // `cxx:3121-3123`
        if hds.ds().shape(i, true).map(|s| s.shape_type()) != Some(ShapeType::Edge) {
            continue;
        }
        let li: InterferenceList = hds.ds().shape_interferences_at(i, true).to_vec();
        // `cxx:3128-3130`
        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);

        // `cxx:3131-3273` : the pass rewrites the `(K, G)` groups in place and
        // rebuilds the shape list afterwards, so the keys are collected first.
        let mut gkeys: Vec<(Kind, i32)> = Vec::new();
        tki.init();
        while tki.more() {
            if let Some((k, g, _)) = tki.value() {
                gkeys.push((k, g));
            }
            tki.next();
        }

        for &(k, g) in &gkeys {
            let loi = tki.interferences(k, g).to_vec();

            // `cxx:3135-3147`
            let mut loicopy = loi.clone();
            let mut l0: InterferenceList = Vec::new();
            fun_select_trasha_interference(&mut loicopy, ShapeType::Edge, &mut l0);
            let mut l1: InterferenceList = Vec::new();
            fun_select_trasha_interference(&mut loicopy, ShapeType::Face, &mut l1);
            let mut l2: InterferenceList = Vec::new();
            if fun_select_sk_interference(&mut l1, Kind::Edge, &mut l2) == 0 {
                continue;
            }
            // `cxx:3150` : `l1 = {I=(T(faceTRASHA),G,Sedge)}`

            // `cxx:3151-3251`
            let mut tkis = TKI::new();
            tkis.fill_on_support(&l2);
            let mut skeys: Vec<(Kind, i32)> = Vec::new();
            tkis.init();
            while tkis.more() {
                if let Some((k2, s2, _)) = tkis.value() {
                    skeys.push((k2, s2));
                }
                tkis.next();
            }

            for &(k2, s2) in &skeys {
                let li2 = tkis.interferences(k2, s2).to_vec();
                // `cxx:3158-3161`
                if li2.len() < 2 {
                    continue;
                }

                // `cxx:3163-3248`
                let mut new_i: Option<InterferenceHandle> = None;
                let mut it1 = 0usize;
                while it1 < li2.len() {
                    let (o1, ib1, ia1, s1) = {
                        let intf = li2[it1].read().expect("interference lock");
                        let o1 = intf.transition().orientation(State::In, ShapeType::Face);
                        let (_sb1, ib1, _sa1, ia1, _gt1, _g1, _st1, s1) = fds_idata(&intf);
                        (o1, ib1, ia1, s1)
                    };
                    // `cxx:3170-3174`
                    if o1 != Orientation::Forward && o1 != Orientation::Reversed {
                        it1 += 1;
                        continue;
                    }
                    // `cxx:3181-3185`
                    if ib1 != ia1 {
                        it1 += 1;
                        continue;
                    }

                    // `cxx:3187-3228`
                    let mut oppofound = false;
                    let mut it2 = it1 + 1;
                    while it2 < li2.len() {
                        let intf2 = li2[it2].read().expect("interference lock");
                        let o2 = intf2.transition().orientation(State::In, ShapeType::Face);
                        let (_sb2, ib2, _sa2, ia2, _gt2, _g2, _st2, s3) = fds_idata(&intf2);
                        drop(intf2);
                        if ib2 != ia2 {
                            it2 += 1;
                            continue;
                        }
                        if ib1 != ib2 {
                            it2 += 1;
                            continue; // same fTRASHA
                        }
                        if s1 != s3 {
                            it2 += 1;
                            continue; // same Sedge
                        }
                        if o2 == o1.complement() {
                            oppofound = true;
                            break;
                        }
                        it2 += 1;
                    }

                    // `cxx:3229-3232`
                    if !oppofound {
                        it1 += 1;
                        continue;
                    }

                    // `cxx:3233-3246`
                    match fun_redusamsha_on_e(hds.ds(), &li2[it1], i) {
                        Some(ni) => {
                            new_i = Some(ni);
                            break;
                        }
                        None => {
                            it1 += 1;
                            continue;
                        }
                    }
                }

                // `cxx:3244-3249`
                let Some(new_i) = new_i else {
                    continue;
                };
                let dst = tkis
                    .change_interferences(k2, s2)
                    .expect("bound support key");
                dst.clear();
                dst.push(new_i);
            }

            // `cxx:3252-3262`
            let mut rebuilt: InterferenceList = Vec::new();
            for (k2, s2) in &skeys {
                rebuilt.extend(tkis.interferences(*k2, *s2).iter().cloned());
            }
            rebuilt.extend(l0);
            rebuilt.extend(l1);
            let newloi = tki.change_interferences(k, g).expect("bound geometry key");
            newloi.clear();
            newloi.extend(rebuilt);
        }

        // `cxx:3266-3274`
        let mut rebuilt: InterferenceList = Vec::new();
        for (k, g) in &gkeys {
            rebuilt.extend(tki.interferences(*k, *g).iter().cloned());
        }
        let new_li = hds.change_ds().change_shape_interferences_at(i);
        new_li.clear();
        new_li.extend(rebuilt);
    }
}
