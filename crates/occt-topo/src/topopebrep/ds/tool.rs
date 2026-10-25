//! `TopOpeBRepDS_TOOL` — same-domain and split-on queries over the DS.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_TOOL.cxx:33-404`,
//! `TopOpeBRepDS_TOOL.hxx:28-90`.
//!
//! The C++ methods take `handle(HDataStructure)` plus out-parameters and return
//! `bool`/`int`; the port returns `Option<...>`/`usize` so the "not found"
//! answers (`cxx:96`, `:109-111`, `:198-200`, `:226`, `:253`, `:270`,
//! `:312-324`, `:342-344`) stay distinguishable from the successful empty ones.
//!
//! `MEspON` is the `NCollection_DataMap<TopoDS_Shape,
//! TopOpeBRepDS_ListOfShapeOn1State, TopTools_ShapeMapHasher>` built by the
//! caller; [`DShapeOn1StateMap`] is its port. `Find` on a missing key throws in
//! OCCT, so the port's `get` returning `None` is the "absent" arm and is
//! reported as `None` by the callers instead of raising.

use std::collections::{HashMap, HashSet};

use occt_core::gp::GpDir;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools::{self, shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::closed::{CLOSING, INTERNAL};
use crate::topopebrep::tool::project_curve::{fun_tool_bounds, fun_tool_par_e_tol};
use crate::topopebrep::tool::tangent::{curves_same_oriented_at, tg_inside};

use super::data_structure::DataStructure;
use super::export::{fds_copy, fun_ds_getoov};
use super::h_data_structure::HDataStructure;
use super::interference_iterator::InterferenceList;
use super::kind::Config;
use super::list_of_shape_on1_state::ListOfShapeOn1State;
use super::process_interferences_tool::fun_select_trasha_interference;

/// `NCollection_DataMap<TopoDS_Shape, TopOpeBRepDS_ListOfShapeOn1State,
/// TopTools_ShapeMapHasher>` (`TopOpeBRepDS_TOOL.hxx:54-66`), keyed by
/// [`ShapeMapKey`] so membership follows `IsSame`.
pub type DShapeOn1StateMap = HashMap<ShapeMapKey, ListOfShapeOn1State>;

/// `static void FDS_sortGb(HDS, LI, LIGb0, LIGb1, LIGbsd)`
/// (`TopOpeBRepDS_TOOL.cxx:32-75`).
///
/// A non-shape-shape interference, or a shape-shape one whose support has no
/// same domain, goes to `ligb0`/`ligbsd` respectively; the remaining split on
/// `GBound()`. `LIGbsd` only receives shape-shape interferences.
fn fds_sort_gb(
    hds: &HDataStructure,
    ds: &DataStructure,
    li: &[super::interference::InterferenceHandle],
    ligb0: &mut InterferenceList,
    ligb1: &mut InterferenceList,
    ligbsd: &mut InterferenceList,
) {
    ligb0.clear();
    ligb1.clear();
    ligbsd.clear();
    for h in li {
        let i = h.read().expect("interference lock");
        let Some(ssi) = i.shape_shape_data() else {
            drop(i);
            ligb0.push(h.clone());
            continue;
        };
        let gb1 = ssi.g_bound;
        let g = i.geometry();
        let hsd = match ds.shape(g, true) {
            Some(shape) => hds.has_same_domain(&shape, true),
            None => false,
        };
        drop(i);
        if hsd {
            ligbsd.push(h.clone());
            continue;
        }
        if gb1 {
            ligb1.push(h.clone());
        } else {
            ligb0.push(h.clone());
        }
    }
}

/// `TopOpeBRepDS_TOOL::EShareG(HDS, E, lEsd)` (`cxx:85-185`).
///
/// Returns the number of same-domain edges sharing `E`'s geometric domain and
/// appends them to `lesd`. A degenerated edge answers its raw same-domain list
/// (`cxx:88-101`); every other edge filters the interferences on `E` by
/// `TRASHA` and keeps the supports whose `GBound` / openness / tangent direction
/// agree.
pub fn e_share_g(hds: &HDataStructure, e: &TopoShape, lesd: &mut Vec<TopoShape>) -> usize {
    lesd.clear();
    let ds = hds.ds();
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        if !hds.has_same_domain(e, true) {
            return 0;
        }
        lesd.extend(hds.same_domain(e).iter().cloned());
        return lesd.len();
    }

    let li = ds.shape_interferences(e, true).to_vec();
    let mut lii = InterferenceList::new();
    fds_copy(&li, &mut lii);
    let mut l1d = InterferenceList::new();
    let n1d = fun_select_trasha_interference(&mut lii, ShapeType::Edge, &mut l1d);
    if n1d == 0 {
        return 0;
    }

    let mut mapesd: Vec<TopoShape> = Vec::new();
    let mut mapesd_keys: HashSet<ShapeMapKey> = HashSet::new();
    let mut l1gb0 = InterferenceList::new();
    let mut l1gb1 = InterferenceList::new();
    let mut l1gbsd = InterferenceList::new();
    fds_sort_gb(hds, ds, &l1d, &mut l1gb0, &mut l1gb1, &mut l1gbsd);

    for h in l1gb0.iter().chain(l1gb1.iter()) {
        let support = h.read().expect("interference lock").support();
        if let Some(s) = ds.shape(support, true) {
            mapesd_push(&mut mapesd, &mut mapesd_keys, s);
        }
    }

    for h in &l1gbsd {
        let (support, g, gb1) = {
            let i = h.read().expect("interference lock");
            let gb1 = i.shape_shape_data().map(|d| d.g_bound).unwrap_or(false);
            (i.support(), i.geometry(), gb1)
        };
        let Some(esd) = ds.shape(support, true) else {
            continue;
        };
        if mapesd_keys.contains(&shape_map_key(&esd)) {
            continue;
        }
        let Some(vg) = ds.shape(g, true) else {
            continue;
        };
        let Some(vsd) = fun_ds_getoov(ds, &vg) else {
            continue;
        };
        let (ve, vesd) = if gb1 { (&vg, &vsd) } else { (&vsd, &vg) };

        let Some((tge, ove)) = tg_inside(ve, e) else {
            continue;
        };
        let Some((tgesd, ovesd)) = tg_inside(vesd, &esd) else {
            continue;
        };
        let ine = ove == CLOSING || ove == INTERNAL;
        let inesd = ovesd == CLOSING || ovesd == INTERNAL;
        if ine || inesd {
            mapesd_push(&mut mapesd, &mut mapesd_keys, esd);
            continue;
        }
        let (Ok(d1), Ok(d2)) = (GpDir::from_vec(&tge), GpDir::from_vec(&tgesd)) else {
            // `gp_Dir(gp_Vec)` raises `Standard_ConstructionError` on a vector
            // shorter than `Precision::Resolution()`; that element is dropped
            // instead of aborting the whole query.
            continue;
        };
        if d1.dot(&d2) > 0.0 {
            mapesd_push(&mut mapesd, &mut mapesd_keys, esd);
        }
    }

    lesd.extend(mapesd);
    lesd.len()
}

/// `NCollection_Map<TopoDS_Shape>::Add`: append on first insertion only, so the
/// result follows the discovery order while duplicates collapse.
fn mapesd_push(map: &mut Vec<TopoShape>, keys: &mut HashSet<ShapeMapKey>, s: TopoShape) {
    if keys.insert(shape_map_key(&s)) {
        map.push(s);
    }
}

/// `TopOpeBRepDS_TOOL::ShareG(HDS, i1, i2)` (`cxx:189-212`): `Shape(i2)` is in
/// the same-domain list of `Shape(i1)`.
pub fn share_g(hds: &HDataStructure, i1: i32, i2: i32) -> bool {
    let ds = hds.ds();
    let (Some(s1), Some(s2)) = (ds.shape(i1, true), ds.shape(i2, true)) else {
        return false;
    };
    if !hds.has_same_domain(&s1, true) {
        return false;
    }
    hds.same_domain(&s1)
        .iter()
        .any(|sdm| topo_tools::is_same(sdm, &s2))
}

/// `TopOpeBRepDS_TOOL::GetEsd(HDS, S, ie, iesd)` (`cxx:216-257`).
///
/// The candidate set is every edge of `S` that the DS knows and that has a same
/// domain; the first element of `SameDomain(Shape(ie))` found in it wins.
pub fn get_esd(hds: &HDataStructure, s: &TopoShape, ie: i32) -> Option<i32> {
    let ds = hds.ds();
    let mut mesd_s: HashSet<ShapeMapKey> = HashSet::new();
    let mut ex = Explorer::new(s, ShapeType::Edge);
    while ex.more() {
        let e = ex.current().clone();
        ex.next();
        if !hds.has_shape(&e, true) {
            continue;
        }
        if !hds.has_same_domain(&e, true) {
            continue;
        }
        mesd_s.insert(shape_map_key(&e));
    }

    let s_ie = ds.shape(ie, true)?;
    for esd in hds.same_domain(&s_ie) {
        if !mesd_s.contains(&shape_map_key(esd)) {
            continue;
        }
        return Some(ds.shape_index(esd, true));
    }
    None
}

/// `TopOpeBRepDS_TOOL::ShareSplitON(HDS, MEspON, i1, i2, spON)`
/// (`cxx:263-325`).
///
/// Both shapes must be same-domain and both must carry a non-empty split-on
/// list; the answer is the first common split.
pub fn share_split_on(
    hds: &HDataStructure,
    me_sp_on: &DShapeOn1StateMap,
    i1: i32,
    i2: i32,
) -> Option<TopoShape> {
    if !share_g(hds, i1, i2) {
        return None;
    }
    let ds = hds.ds();
    let s1 = ds.shape(i1, true)?;
    let s2 = ds.shape(i2, true)?;

    let los1 = me_sp_on.get(&shape_map_key(&s1))?;
    if !los1.is_split() {
        return None;
    }
    let lsp1 = los1.list_on_state();
    if lsp1.is_empty() {
        return None;
    }
    let mesp1: HashSet<ShapeMapKey> = lsp1.iter().map(shape_map_key).collect();

    let los2 = me_sp_on.get(&shape_map_key(&s2))?;
    if !los2.is_split() {
        return None;
    }
    let lsp2 = los2.list_on_state();
    if lsp2.is_empty() {
        return None;
    }
    for esp in lsp2 {
        if mesp1.contains(&shape_map_key(esp)) {
            return Some(esp.clone());
        }
    }
    None
}

/// `TopOpeBRepDS_TOOL::GetConfig(HDS, MEspON, ie, iesd, config)`
/// (`cxx:332-404`).
///
/// When neither edge has `UNSHGEOMETRY` the config is read from
/// `SameDomainOrientation`, XOR'd with the two edge orientations. Otherwise the
/// two edges are compared through a common split-on edge: a point at
/// `x = 0.45678` of that edge is transferred to both and the tangents are
/// compared at the transferred parameters.
pub fn get_config(
    hds: &HDataStructure,
    me_sp_on: &DShapeOn1StateMap,
    ie: i32,
    iesd: i32,
) -> Option<Config> {
    if !share_g(hds, ie, iesd) {
        return None;
    }
    let ds = hds.ds();
    let e = ds.shape(ie, true)?;
    let esd = ds.shape(iesd, true)?;

    let conf = hds.same_domain_orientation(&e);
    let unsh = conf == Config::UnshGeometry;
    let confsd = hds.same_domain_orientation(&esd);
    let unshsd = confsd == Config::UnshGeometry;
    if !unsh && !unshsd {
        let mut sameori = conf == confsd;
        if e.orientation() == crate::abs::Orientation::Reversed {
            sameori = !sameori;
        }
        if esd.orientation() == crate::abs::Orientation::Reversed {
            sameori = !sameori;
        }
        return Some(if sameori {
            Config::SameOriented
        } else {
            Config::DiffOriented
        });
    }

    let e_on = share_split_on(hds, me_sp_on, ie, iesd)?;
    let (f, l) = fun_tool_bounds(&e_on);
    let x = 0.45678;
    let par_on = (1.0 - x) * f + x * l;
    let tole = BRepTool::edge_tolerance(&Edge(e.clone()));
    let pare = fun_tool_par_e_tol(&e_on, par_on, &e, tole)?;
    let tolesd = BRepTool::edge_tolerance(&Edge(esd.clone()));
    let paresd = fun_tool_par_e_tol(&e_on, par_on, &esd, tolesd)?;
    let so = curves_same_oriented_at(&e, pare, &esd, paresd)?;
    Some(if so {
        Config::SameOriented
    } else {
        Config::DiffOriented
    })
}
