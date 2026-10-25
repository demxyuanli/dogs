//! `FUN_ds_completeforSE3` and its static transition helper.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:1225-1421`.

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::{
    fds_data, fds_parameter, fds_si_is_gi_of_i_of_sba_of_t_of_i, fds_tdata,
};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_itrasha_interference, fun_select_sk_interference, fun_select_traori_interference,
    fun_select_trasha_interference, make_epv_interference_with_support,
};
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::make_transition::MakeTransition;
use crate::topopebrep::tool::project_curve::{fun_tool_bounds, fun_tool_par_e};
use crate::topopebrep::tool::topology::fun_tool_par_on_e_f;

/// `FUN_ds_completeforSE3(BDS, SE, K, G, loi, parE, IES, ITRASHA, Tr)`
/// (`cxx:1225-1349`): the static helper of `FUN_ds_completeforSE3`.
///
/// `G` is not used by OCCT (the parameter is unnamed at `cxx:1230`).
pub fn fun_se3_transition(
    bds: &DataStructure,
    se: &TopoShape,
    k: Kind,
    loi: &InterferenceList,
) -> Option<(f64, i32, i32, Transition)> {
    // `cxx:1238-1241`
    if k == Kind::Vertex {
        return None;
    }
    let ise = bds.shape_index(se, true);

    // `cxx:1243-1257`
    let mut loicopy = loi.clone();
    let mut l_int: InterferenceList = Vec::new();
    let n_int =
        fun_select_traori_interference(&mut loicopy, Orientation::Internal, &mut l_int);
    let mut l_ext: InterferenceList = Vec::new();
    let n_ext =
        fun_select_traori_interference(&mut loicopy, Orientation::External, &mut l_ext);
    if n_int + n_ext < 1 {
        return None;
    }
    let mut l1: InterferenceList = Vec::new();
    l1.extend(l_int);
    l1.extend(l_ext);

    // `cxx:1259-1278` : b. I = (INT/EXT(F), G, S), F has I with G == ES
    let mut l2: InterferenceList = Vec::new();
    if fun_select_trasha_interference(&mut l1, ShapeType::Face, &mut l2) < 1 {
        return None;
    }
    let okb = l2.iter().any(|h| {
        let intf = h.read().expect("interference lock");
        fds_si_is_gi_of_i_of_sba_of_t_of_i(bds, ise, &intf)
    });
    if !okb {
        return None;
    }

    // `cxx:1281-1296` : a. I3d = IFE + IF
    let mut l2e: InterferenceList = Vec::new();
    let n2e = fun_select_sk_interference(&mut l2, Kind::Edge, &mut l2e);
    let mut l2f: InterferenceList = Vec::new();
    let n2f = fun_select_sk_interference(&mut l2, Kind::Face, &mut l2f);
    if n2e == 0 || n2f == 0 {
        return None;
    }
    let s_i3d = l2f[0].read().expect("interference lock").support();
    let mut l3d: InterferenceList = Vec::new();
    if fun_select_itrasha_interference(&mut l2e, s_i3d, &mut l3d) < 1 {
        return None;
    }

    // `cxx:1299-1311`
    let (par_e, ies, itrasha) = {
        let i3 = l3d[0].read().expect("interference lock");
        let (_k3, _g3, _st3, s3) = fds_data(&i3);
        let (_tsb3, isb3, _tsa3, _isa3) = fds_tdata(&i3.transition());
        let Some(par_e) = fds_parameter(&i3) else {
            return None;
        };
        (par_e, s3, isb3)
    };

    // `cxx:1311-1325`
    let eline = bds.shape(ies, true)?;
    let f = bds.shape(itrasha, true)?;
    let parline = fun_tool_par_e(se, par_e, &eline)?;
    let uv = fun_tool_par_on_e_f(&eline, parline, &f)?;

    // `cxx:1327-1348`
    let (par1, par2) = fun_tool_bounds(se);
    let factor = 1.0e-4;
    let mut mkt = MakeTransition::new();
    let mut stb = crate::abs::State::Unknown;
    let mut sta = crate::abs::State::Unknown;
    let mut ok = mkt.initialize(se, par1, par2, par_e, &f, &uv, factor);
    if ok {
        ok = mkt.set_rest(&eline, parline);
    }
    if ok {
        ok = mkt.mk_ton_e(&mut stb, &mut sta);
    }
    if !ok {
        return None;
    }
    let mut tr = Transition::default();
    tr.set_state_before(stb);
    tr.set_state_after(sta);
    tr.set_index(itrasha);
    Some((par_e, ies, itrasha, tr))
}

/// `FUN_ds_completeforSE3(HDS)` (`cxx:1352-1421`).
///
/// xpu170498: when the transition of `FTRA` is `INTERNAL`/`EXTERNAL` the
/// `spIN`/`spOU` computation fails, so the 3d interference `IFE + IF` is
/// recomputed with a `FORWARD`/`REVERSED` transition.
pub fn fun_ds_complete_for_se3(hds: &mut HDataStructure) {
    let nse = hds.ds().nb_section_edges();
    for i in 1..=nse {
        // `cxx:1370-1373`
        let Some(se) = hds.ds().section_edge(i, true) else {
            continue;
        };
        let ise = hds.ds().shape_index(&se, true);
        let li: InterferenceList = hds.ds().shape_interferences(&se, true).to_vec();

        let mut new_li: InterferenceList = Vec::new();
        let mut hasnew_li = false;
        let mut tki = TKI::new();
        tki.fill_on_geometry(&li);
        tki.init();
        while tki.more() {
            // `cxx:1382-1389`
            let (k, _g, loi) = {
                let Some(v) = tki.value() else {
                    tki.next();
                    continue;
                };
                (v.0, v.1, v.2.to_vec())
            };
            let mut loicopy = loi.clone();
            let tr_opt = fun_se3_transition(hds.ds(), &se, k, &loi);

            // `cxx:1390-1395`
            let Some((par_e, _ies, itrasha, tr)) = tr_opt else {
                new_li.extend(loicopy);
                tki.next();
                continue;
            };

            // `cxx:1398-1399` : computed and unused in OCCT
            let mut l_idel: InterferenceList = Vec::new();
            fun_select_itrasha_interference(&mut loicopy.clone(), itrasha, &mut l_idel);

            // `cxx:1401-1408`
            hasnew_li = true;
            if let Some(tr_ff) =
                make_epv_interference_with_support(&tr, itrasha, _g, par_e, k, Kind::Face, false)
            {
                new_li.push(tr_ff);
            }
            if let Some(tr_fe) =
                make_epv_interference_with_support(&tr, ise, _g, par_e, k, Kind::Edge, false)
            {
                new_li.push(tr_fe);
            }
            new_li.extend(loicopy);
            tki.next();
        }

        // `cxx:1410-1419`
        if hasnew_li {
            let lii = hds.change_ds().change_shape_interferences(&se);
            lii.clear();
            lii.append(&mut new_li);
        }
    }
}
