//! The edge-interference processing helpers of `TopOpeBRepDS_ProcessEdgeInterferences.cxx`.
//!
//! Sources (OCCT 8.0.0):
//! - `FUN_vertexofedge` (`cxx:35-53`)
//! - `FUN_keepEinterference` (`cxx:55-131`)
//! - `FUN_unkeepEinterferences` (`cxx:134-167`)
//! - `FUN_unkeepEsymetrictransitions` (`cxx:245-350`)
//! - `FUN_orderFFsamedomain` (`cxx:352-405`)
//! - `FUN_orderSTATETRANSonG` (`cxx:408-453`)
//! - `FUN_orderSTATETRANS` (`cxx:455-486`)
//! - `FUN_resolveEUNKNOWN` (`cxx:488-571`)
//! - `FUN_purgeDSonSE` (`cxx:575-685`)
//!
//! `FUN_changeFOUT` (`cxx:169-242`) is enclosed in `/* ... */` in OCCT and stays
//! commented out here.

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools;
use crate::topopebrep::tool::sc::fsc_get_psc;

use crate::topopebrep::ds::export::{fds_data, fds_idata, fds_tdata};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_interf_has_support, fun_transition_shape_equal, fun_unkeep_unknown,
};
use crate::topopebrep::ds::tki::TKI;

/// `FUN_vertexofedge(V, E)` (`cxx:35-53`).
pub fn fun_vertexofedge(v: &TopoShape, e: &TopoShape) -> bool {
    let mut ex = Explorer::new(e, ShapeType::Vertex);
    while ex.more() {
        let cur = ex.current().clone();
        ex.next();
        if topo_tools::is_same(&cur, v) {
            return true;
        }
    }
    false
}

/// `FUN_keepEinterference(DS, I, E)` (`cxx:55-131`).
///
/// The `TopExp::Vertices`/`ShapeSameDomain` loop at `cxx:94-109` computes
/// nothing: its result `k4` is only used by the commented `res = res && k4`
/// (`cxx:110`), so it is not translated. The same holds for the `VG` lookup at
/// `cxx:82`, whose only readers are the commented blocks at `cxx:84-88` and
/// `cxx:94-109`.
pub fn fun_keep_e_interference(
    ds: &crate::topopebrep::ds::data_structure::DataStructure,
    i: &crate::topopebrep::ds::interference::Interference,
    e: &TopoShape,
) -> bool {
    let mut res = true;
    if i.is_edge_vertex() {
        // EVI I rejected when its support edge is E accessing I.
        res = !fun_interf_has_support(ds, i, e);

        // EVI rejected when the transition is ON EDGE before or after.
        let t = i.transition();
        let (shab, shaa) = (t.shape_before(), t.shape_after());
        let (stab, staa) = (t.before_state(), t.after_state());
        let k2 = (shab != ShapeType::Edge || stab != State::On)
            && (shaa != ShapeType::Edge || staa != State::On);
        res = res && k2;
    } else if i.is_curve_point() {
        // MSV Oct 4, 2001 (cxx:116-126): reject an interference whose parameter
        // coincides with one side of the edge range.
        let eps = occt_core::precision::PCONFUSION;
        let Some(par) = i.parameter() else {
            return res;
        };
        let (f, l) = BRepTool::edge_parameters(&Edge(e.clone()));
        if (par - f).abs() < eps || (par - l).abs() < eps {
            res = false;
        }
    }
    res
}

/// `FUN_unkeepEinterferences(LI, BDS, SIX)` (`cxx:133-167`): returns
/// `LI.Extent()`.
pub fn fun_unkeep_e_interferences(
    li: &mut InterferenceList,
    ds: &crate::topopebrep::ds::data_structure::DataStructure,
    six: i32,
) -> usize {
    let Some(e) = ds.shape(six, true) else {
        // `BDS.Shape(SIX)` would raise; an absent shape keeps the list.
        return li.len();
    };
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        return li.len();
    }
    let mut i = 0;
    while i < li.len() {
        let keep = {
            let h = li[i].read().expect("interference lock");
            fun_keep_e_interference(ds, &h, &e)
        };
        if !keep {
            li.remove(i);
        } else {
            i += 1;
        }
    }
    li.len()
}

/// `FUN_unkeepEsymetrictransitions(LI, BDS, SIX)` (`cxx:245-350`).
///
/// Not called by `FUN_FilterEdge` (the call is commented at
/// `FilterEdgeInterferences.cxx:84`) but exported by OCCT, so it is ported.
pub fn fun_unkeep_e_symmetric_transitions(
    li: &mut InterferenceList,
    ds: &crate::topopebrep::ds::data_structure::DataStructure,
    six: i32,
) {
    let Some(e) = ds.shape(six, true) else {
        return;
    };
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        return;
    }

    let mut i = 0;
    while i < li.len() {
        let mut it1_toremove = false;
        let (t1, gt1, g1, st1, s1, tsb1, isb1, tsa1, isa1) = {
            let h = li[i].read().expect("interference lock");
            let (gt, g, st, s) = fds_data(&h);
            let (tsb, isb, tsa, isa) = fds_tdata(&h.transition());
            (h.transition(), gt, g, st, s, tsb, isb, tsa, isa)
        };
        let _ = (gt1, g1, st1, s1, tsb1, isb1, tsa1, isa1);

        let mut j = i + 1;
        while j < li.len() {
            let (t2, cond) = {
                let h = li[j].read().expect("interference lock");
                let (gt2, g2, st2, s2) = fds_data(&h);
                let (tsb2, isb2, tsa2, isa2) = fds_tdata(&h.transition());
                let id_gs = gt2 == gt1 && g2 == g1 && st2 == st1 && s2 == s1;
                // `++971219` (`cxx:294-297`): the shape indices must also match.
                let id_iba = isb1 == isb2 && isa1 == isa2;
                (h.transition(), id_gs && id_iba)
            };
            if !cond {
                j += 1;
                continue;
            }
            let idshape = fun_transition_shape_equal(&t1, &t2);
            // `idstate` is computed at `cxx:301` but unused: `newV` is a
            // constant `true`, so only the `oppostate` arm can be taken.
            let _idstate = crate::topopebrep::ds::process_interferences_tool::fun_transition_state_equal(&t1, &t2);
            let (t1b, t1a) = (t1.before_state(), t1.after_state());
            let (t2b, t2a) = (t2.before_state(), t2.after_state());
            let out_in = |a: State, b: State| a == State::Out && b == State::In;
            let mut oppostate = out_in(t1b, t2b) || out_in(t2b, t1b);
            oppostate = oppostate && (out_in(t1a, t2a) || out_in(t2a, t1a));
            let sym = idshape && oppostate;
            if sym {
                li.remove(j);
                it1_toremove = true;
            } else {
                j += 1;
            }
        }

        if it1_toremove {
            li.remove(i);
        } else {
            i += 1;
        }
    }
}

/// `FUN_orderFFsamedomain(LI, HDS, SIX)` (`cxx:352-405`): the interferences
/// whose two transition shapes are same-domain faces go first.
pub fn fun_order_ff_same_domain(li: &mut InterferenceList, hds: &HDataStructure, _six: i32) {
    let ds = hds.ds();
    let mut li_ffsd = InterferenceList::new();
    let mut li_other = InterferenceList::new();
    let all = std::mem::take(li);
    for h in all {
        let (tsb1, isb1, tsa1, isa1) = {
            let i = h.read().expect("interference lock");
            fds_tdata(&i.transition())
        };
        let mut ffsd = false;
        if tsb1 == ShapeType::Face && tsa1 == ShapeType::Face {
            // T1 states are defined on FACEs.
            let ffsd1 = ds
                .shape(isb1, true)
                .map(|s| hds.has_same_domain(&s, true))
                .unwrap_or(false);
            let ffsd2 = ds
                .shape(isa1, true)
                .map(|s| hds.has_same_domain(&s, true))
                .unwrap_or(false);
            ffsd = ffsd1 && ffsd2;
        }
        if ffsd {
            li_ffsd.push(h);
        } else {
            li_other.push(h);
        }
    }
    li.extend(li_ffsd);
    li.extend(li_other);
}

/// `FUN_orderSTATETRANSonG(LI, HDS, SIX)` (`cxx:407-453`): the interferences
/// whose transition has equal states go first inside one `(K, G)` group.
pub fn fun_order_state_trans_on_g(li: &mut InterferenceList, _hds: &HDataStructure, _six: i32) {
    let mut l1 = InterferenceList::new();
    let mut l2 = InterferenceList::new();
    let all = std::mem::take(li);
    for h in all {
        let (tsb1, isb1, tsa1, isa1, _gt1, _g1, _st1, _s1, stab1, staa1) = {
            let i = h.read().expect("interference lock");
            let (tsb, isb, tsa, isa, gt, g, st, s) = fds_idata(&i);
            let t = i.transition();
            (tsb, isb, tsa, isa, gt, g, st, s, t.before_state(), t.after_state())
        };
        let steq = tsb1 == tsa1 && isb1 == isa1 && stab1 == staa1;
        if steq {
            l1.push(h);
        } else {
            l2.push(h);
        }
    }
    li.extend(l1);
    li.extend(l2);
}

/// `FUN_orderSTATETRANS(LI, HDS, SIX)` (`cxx:455-487`): `FUN_orderSTATETRANSonG`
/// per `(GeometryType, Geometry)` group, then the groups are concatenated.
///
/// The group order is the one [`TKI`] iterates (ascending `(Kind, index)`); see
/// `tki.rs` for why that replaces OCCT's `NCollection_DataMap` hash order.
pub fn fun_order_state_trans(li: &mut InterferenceList, hds: &HDataStructure, six: i32) {
    let mut tki = TKI::new();
    tki.fill_on_geometry(li);
    tki.init();
    while tki.more() {
        if let Some((_k, _g, loi)) = tki.change_value() {
            fun_order_state_trans_on_g(loi, hds, six);
        }
        tki.next();
    }

    li.clear();
    tki.init();
    while tki.more() {
        let group = tki.value().map(|(_k, _g, l)| l.clone());
        if let Some(g) = group {
            li.extend(g);
        }
        tki.next();
    }
}

/// `FUN_resolveEUNKNOWN(LI, BDS, SIX)` (`cxx:489-571`).
///
/// An edge tangent to a face at one point classifies the two points at
/// `0.8·par` and `0.2·par` of the interference parameter against a one-face
/// solid built from the support face, and writes the two states back into the
/// transition.
pub fn fun_resolve_e_unknown(
    li: &mut InterferenceList,
    ds: &crate::topopebrep::ds::data_structure::DataStructure,
    six: i32,
) {
    if let Some(e) = ds.shape(six, true) {
        let ee = Edge(e.clone());
        let (fe, le) = BRepTool::edge_parameters(&ee);
        for h in li.iter() {
            let (is_unk, t1, tsb1, isb1, tsa1, isa1, gt1, g1, st1, s1, parcpi, is_cpi) = {
                let i = h.read().expect("interference lock");
                let (tsb, isb, tsa, isa, gt, g, st, s) = fds_idata(&i);
                (
                    i.transition().is_unknown(),
                    i.transition(),
                    tsb,
                    isb,
                    tsa,
                    isa,
                    gt,
                    g,
                    st,
                    s,
                    i.parameter(),
                    i.is_curve_point(),
                )
            };
            if !is_unk {
                continue;
            }
            let _ = (t1, g1);
            let idt = tsb1 == ShapeType::Face
                && tsa1 == ShapeType::Face
                && gt1 == Kind::Point
                && st1 == Kind::Face;
            let idi = isb1 == s1 && isa1 == s1;
            // `etgf`: edge tangent to a face at one point.
            if !(idt && idi) {
                continue;
            }
            if !is_cpi {
                continue;
            }
            let Some(parcpi) = parcpi else {
                continue;
            };
            // NYI (`cxx:536-538`): get the points from the 2d curve.
            let Some(ce) = BRepTool::edge_curve_world(&ee) else {
                continue;
            };
            let ttb = 0.8;
            let parb = (1.0 - ttb) * fe + ttb * parcpi;
            let tta = 0.2;
            let para = (1.0 - tta) * parcpi + tta * le;
            let pb = ce.d0(parb);
            let pa = ce.d0(para);

            let Some(fs) = ds.shape(s1, true) else {
                continue;
            };
            let builder = TopoBuilder::new();
            let sh = builder.make_shell(&[Face(fs.clone())]);
            let so = builder.make_solid(&[sh]);
            let mut psc = fsc_get_psc(&so.0);
            if psc.state_p3d_reference(&pb).is_err() {
                continue;
            }
            let stateb = psc.state();
            if psc.state_p3d_reference(&pa).is_err() {
                continue;
            }
            let statea = psc.state();
            if stateb == State::Unknown || statea == State::Unknown {
                continue;
            }
            let mut i1 = h.write().expect("interference lock");
            let new_t1 = i1.change_transition();
            new_t1.set(stateb, statea, tsb1, tsa1);
        }
    }
    fun_unkeep_unknown(li);
}

/// `FUN_purgeDSonSE(HDS, EIX, LI)` (`cxx:573-688`).
///
/// Attached to a section edge `SE`, an interference `(T(Fsdm), G, S)` whose
/// support face `Fsdm` shares the same domain with an ancestor face of `SE` is
/// dropped: it would carry 3d information the split-on edge must not have
/// (`PRO12660`, `spON(e48)`).
pub fn fun_purge_ds_on_se(
    hds: &HDataStructure,
    connex: &crate::topopebrep::ds::connex::Connexity,
    eix: i32,
    li: &mut InterferenceList,
) {
    let ds = hds.ds();
    let Some(se) = ds.shape(eix, true) else {
        return;
    };
    let rk_se = ds.ancestor_rank_of(&se);
    let isse = ds.is_section_edge(&se, true);
    if !isse {
        return;
    }

    let mut fsdm_fanc_se: Vec<TopoShape> = Vec::new();
    let mut fsdm_keys: Vec<crate::topo_tools::ShapeMapKey> = Vec::new();
    let l_fanc_se = connex.edge_connexity_same_shape(&se, hds);
    for f in l_fanc_se {
        for sdmf in hds.same_domain(f) {
            let sdmrk = ds.ancestor_rank_of(sdmf);
            if sdmrk == rk_se {
                continue;
            }
            let key = crate::topo_tools::shape_map_key(sdmf);
            if !fsdm_keys.contains(&key) {
                fsdm_keys.push(key);
                fsdm_fanc_se.push(sdmf.clone());
            }
        }
    }
    if fsdm_fanc_se.is_empty() {
        return;
    }

    let mut new_li = InterferenceList::new();
    let mut tki = TKI::new();
    tki.fill_on_geometry(li);
    tki.init();
    while tki.more() {
        let group = tki.value().map(|(_k, _g, l)| l.clone());
        if let Some(loi) = group {
            for h in loi {
                let (gt, _g, st, _s, tsb, isb, tsa, isa) = {
                    let i = h.read().expect("interference lock");
                    let (gt, g, st, s) = fds_data(&i);
                    let (tsb, isb, tsa, isa) = fds_tdata(&i.transition());
                    (gt, g, st, s, tsb, isb, tsa, isa)
                };
                let _ = (gt, tsa, isa);
                if tsb != ShapeType::Face {
                    new_li.push(h);
                    continue;
                }
                if st != Kind::Face {
                    new_li.push(h);
                    continue;
                }
                let isbound = match ds.shape(isb, true) {
                    Some(f) => {
                        let key = crate::topo_tools::shape_map_key(&f);
                        fsdm_keys.contains(&key)
                    }
                    None => false,
                };
                if !isbound {
                    new_li.push(h);
                }
            }
        }
        tki.next();
    }

    li.clear();
    li.extend(new_li);
}
