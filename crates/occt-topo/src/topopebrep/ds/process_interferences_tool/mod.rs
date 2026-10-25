//! `TopOpeBRepDS_ProcessInterferencesTool` — the interference-list filters used
//! by the boolean builder.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_ProcessInterferencesTool.hxx:1-120`,
//! `TopOpeBRepDS_ProcessInterferencesTool.cxx:36-530`.
//!
//! Every `FUN_select*` filter has the same shape: walk `L1` with a list
//! iterator, move the selected element to the end of `L2` and remove it from
//! `L1` (or step the iterator otherwise), then return `L2.Extent()` — *not* the
//! number of elements selected by this call (`cxx:136-137`). The port keeps the
//! two-list signature: `l1: &mut InterferenceList`, `l2: &mut InterferenceList`.

use crate::abs::{Orientation, ShapeType, State};

use super::data_structure::DataStructure;
use super::interference::{Interference, InterferenceHandle};
use super::interference_iterator::InterferenceList;
use super::interference_tool::{
    make_curve_interference, make_edge_interference, make_edge_vertex_interference, make_handle,
};
use super::kind::Kind;
use super::transition::Transition;

mod select_2d;
mod select_3d;

pub use select_2d::{fun_select_1d_i, fun_select_2d_i, fun_select_pure2d_i};
pub use select_3d::{fun_select_3d_i, fun_select_3d_interference, fun_select_3d_is_esame_isf};

/// Moves the elements of `l1` that satisfy `sel` to the end of `l2`, in order,
/// and returns `l2.len()` (`cxx:118-137` and the eight sibling functions).
fn select_by<F>(l1: &mut InterferenceList, l2: &mut InterferenceList, mut sel: F) -> i32
where
    F: FnMut(&Interference) -> bool,
{
    let mut i = 0;
    while i < l1.len() {
        let keep = {
            let h = l1[i].read().expect("interference lock");
            sel(&h)
        };
        if keep {
            let moved = l1.remove(i);
            l2.push(moved);
        } else {
            i += 1;
        }
    }
    l2.len() as i32
}

/// `MakeCPVInterference(T, SI, GI, P, GK)` (`cxx:38-49`): a curve/point
/// interference with support type `CURVE`.
pub fn make_cpv_interference(
    t: &Transition,
    si: i32,
    gi: i32,
    p: f64,
    gk: Kind,
) -> InterferenceHandle {
    make_curve_interference(t, Kind::Curve, si, gk, gi, p)
}

/// `MakeEPVInterference(T, SI, GI, P, GK, B)` (`cxx:52-76`).
///
/// `GK == POINT` makes an edge interference; `GK == VERTEX` makes an
/// `EdgeVertexInterference` with `TopOpeBRepDS_UNSHGEOMETRY` whose geometry
/// type is then overwritten with `GK` (`cxx:73`). Any other `GK` returns a null
/// handle in OCCT (the local `I` stays null), which the port spells `None`.
pub fn make_epv_interference(
    t: &Transition,
    si: i32,
    gi: i32,
    p: f64,
    gk: Kind,
    b: bool,
) -> Option<InterferenceHandle> {
    match gk {
        Kind::Point => Some(make_edge_interference(t, Kind::Edge, si, gk, gi, p)),
        Kind::Vertex => {
            let mut i = Interference::edge_vertex_on_edge(
                t.clone(),
                si,
                gi,
                b,
                super::kind::Config::UnshGeometry,
                p,
            );
            i.set_geometry_type(gk);
            Some(make_handle(i))
        }
        _ => None,
    }
}

/// `MakeEPVInterference(T, S, G, P, GK, SK, B)` (`cxx:81-94`): the same as
/// [`make_epv_interference`] with the support type overwritten by `SK`.
pub fn make_epv_interference_with_support(
    t: &Transition,
    s: i32,
    g: i32,
    p: f64,
    gk: Kind,
    sk: Kind,
    b: bool,
) -> Option<InterferenceHandle> {
    let i = make_epv_interference(t, s, g, p, gk, b)?;
    i.write().expect("interference lock").set_support_type(sk);
    Some(i)
}

/// `FUN_hasStateShape(T, state, shape)` (`cxx:96-108`).
pub fn fun_has_state_shape(t: &Transition, state: State, shape: ShapeType) -> bool {
    let b = t.before_state() == state && t.shape_before() == shape;
    let a = t.after_state() == state && t.shape_after() == shape;
    a || b
}

/// `FUN_selectTRASHAinterference(L1, sha2, L2)` (`cxx:111-138`): select the
/// interferences whose transition is `(_, sha2)` before or after.
pub fn fun_select_trasha_interference(
    l1: &mut InterferenceList,
    sha2: ShapeType,
    l2: &mut InterferenceList,
) -> i32 {
    select_by(l1, l2, |i| {
        let t = i.transition();
        t.shape_before() == sha2 || t.shape_after() == sha2
    })
}

/// `FUN_selectITRASHAinterference(L1, Index, L2)` (`cxx:141-170`).
///
/// The index of the transition is `Index()` (`transition.rs::index`), which
/// raises in OCCT when the two sides disagree; such an element is *not*
/// selected and the walk continues, matching a non-matching index.
pub fn fun_select_itrasha_interference(
    l1: &mut InterferenceList,
    index: i32,
    l2: &mut InterferenceList,
) -> i32 {
    if index == 0 {
        return 0;
    }
    select_by(l1, l2, |i| i.transition().index().map(|v| v == index).unwrap_or(false))
}

/// `FUN_selectTRAUNKinterference(L1, L2)` (`cxx:173-197`).
pub fn fun_select_traunk_interference(l1: &mut InterferenceList, l2: &mut InterferenceList) -> i32 {
    select_by(l1, l2, |i| i.transition().is_unknown())
}

/// `FUN_selectTRAORIinterference(L1, O, L2)` (`cxx:200-226`): the orientation
/// is `Transition::Orientation(TopAbs_IN)` with the default shape type.
pub fn fun_select_traori_interference(
    l1: &mut InterferenceList,
    o: Orientation,
    l2: &mut InterferenceList,
) -> i32 {
    select_by(l1, l2, |i| {
        i.transition().orientation(State::In, ShapeType::Face) == o
    })
}

/// `FUN_selectGKinterference(L1, GK, L2)` (`cxx:229-251`).
pub fn fun_select_gk_interference(
    l1: &mut InterferenceList,
    gk: Kind,
    l2: &mut InterferenceList,
) -> i32 {
    select_by(l1, l2, |i| i.geometry_type() == gk)
}

/// `FUN_selectSKinterference(L1, SK, L2)` (`cxx:254-277`).
pub fn fun_select_sk_interference(
    l1: &mut InterferenceList,
    sk: Kind,
    l2: &mut InterferenceList,
) -> i32 {
    select_by(l1, l2, |i| i.support_type() == sk)
}

/// `FUN_selectGIinterference(L1, GI, L2)` (`cxx:280-307`).
pub fn fun_select_gi_interference(
    l1: &mut InterferenceList,
    gi: i32,
    l2: &mut InterferenceList,
) -> i32 {
    if gi == 0 {
        return 0;
    }
    select_by(l1, l2, |i| i.geometry() == gi)
}

/// `FUN_selectSIinterference(L1, SI, L2)` (`cxx:310-337`).
pub fn fun_select_si_interference(
    l1: &mut InterferenceList,
    si: i32,
    l2: &mut InterferenceList,
) -> i32 {
    if si == 0 {
        return 0;
    }
    select_by(l1, l2, |i| i.support() == si)
}

/// `FUN_interfhassupport(DS, I, S)` (`cxx:340-352`): `DS.Shape(I->Support())`
/// `IsSame(S)`. A support index outside the DS has no shape in OCCT either, so
/// the port reports `false`.
pub fn fun_interf_has_support(ds: &DataStructure, i: &Interference, s: &crate::shape::TopoShape) -> bool {
    match ds.shape(i.support(), true) {
        Some(sofi) => crate::topo_tools::is_same(&sofi, s),
        None => false,
    }
}

/// `FUN_transitionEQUAL(T1, T2)` (`cxx:348-364`).
pub fn fun_transition_equal(t1: &Transition, t2: &Transition) -> bool {
    fun_transition_state_equal(t1, t2)
        && fun_transition_shape_equal(t1, t2)
        && fun_transition_index_equal(t1, t2)
}

/// `FUN_transitionSTATEEQUAL(T1, T2)` (`cxx:367-373`).
pub fn fun_transition_state_equal(t1: &Transition, t2: &Transition) -> bool {
    t1.before_state() == t2.before_state() && t1.after_state() == t2.after_state()
}

/// `FUN_transitionSHAPEEQUAL(T1, T2)` (`cxx:376-382`).
pub fn fun_transition_shape_equal(t1: &Transition, t2: &Transition) -> bool {
    t1.shape_before() == t2.shape_before() && t1.shape_after() == t2.shape_after()
}

/// `FUN_transitionINDEXEQUAL(T1, T2)` (`cxx:385-391`).
pub fn fun_transition_index_equal(t1: &Transition, t2: &Transition) -> bool {
    t1.index_before() == t2.index_before() && t1.index_after() == t2.index_after()
}

/// `FUN_unkeepUNKNOWN(LI, BDS, SIX)` (`cxx:501-528`): drop the interferences
/// whose transition is `(UNKNOWN, UNKNOWN)`. `BDS` and `SIX` are unused by
/// OCCT (the shape lookup is commented out at `cxx:506`).
pub fn fun_unkeep_unknown(l1: &mut InterferenceList) {
    let mut i = 0;
    while i < l1.len() {
        let is_unknown = l1[i]
            .read()
            .expect("interference lock")
            .transition()
            .is_unknown();
        if is_unknown {
            l1.remove(i);
        } else {
            i += 1;
        }
    }
}

/// `FUN_reducedoublons(LI, BDS, SIX)` (`TopOpeBRepDS_ProcessInterferencesTool.cxx:394-497`).
///
/// Two interferences with the same `FDS_data` quadruple and the same transition
/// are duplicates unless they are `EdgeVertexInterference`s on *different*
/// parameters of a non-closing vertex: the second of each duplicate pair is
/// removed from `LI`.
///
/// The closing-vertex arm is transcribed as OCCT writes it: `cxx:458` is
/// `onvclo = vG.IsSame(vG)`, a self-comparison that is always true, so whenever
/// `ClosedE(E)` succeeds the pair is treated as duplicate and the
/// `hasOO`/parameter branch at `cxx:459-477` never runs. That is the OCCT 8.0.0
/// control flow and is kept verbatim, not "fixed".
pub fn fun_reducedoublons(li: &mut InterferenceList, ds: &DataStructure, six: i32) {
    let edge = ds.shape(six, true);

    let mut i = 0;
    while i < li.len() {
        let (t1, g1, gt1, st1, s1) = {
            let i1 = li[i].read().expect("interference lock");
            (i1.transition(), i1.geometry(), i1.geometry_type(), i1.support_type(), i1.support())
        };

        let mut j = i + 1;
        while j < li.len() {
            let (t2, g2, gt2, st2, s2, evi2, par2) = {
                let i2 = li[j].read().expect("interference lock");
                (
                    i2.transition(),
                    i2.geometry(),
                    i2.geometry_type(),
                    i2.support_type(),
                    i2.support(),
                    i2.as_edge_vertex().is_some(),
                    i2.parameter(),
                )
            };

            let id_gs = gt2 == gt1 && g2 == g1 && st2 == st1 && s2 == s1;
            if !id_gs {
                j += 1;
                continue;
            }

            let id1 = fun_transition_state_equal(&t1, &t2);
            let id2 = fun_transition_shape_equal(&t1, &t2);
            let id3 = fun_transition_index_equal(&t1, &t2);
            let mut id_t = id1 && id2 && id3;

            if id_t {
                let evi1 = li[i]
                    .read()
                    .expect("interference lock")
                    .as_edge_vertex()
                    .is_some();
                if evi1 && evi2 {
                    // `xpu130898: if vG is on E's closing vertex -> filter`
                    let Some(vg) = ds.shape(g1, true) else {
                        // `BDS.Shape(G1)` would raise; a valid DS always has it.
                        j += 1;
                        continue;
                    };
                    let has_oo = super::export::fun_ds_getoov(ds, &vg).is_some();
                    let e_closed = edge
                        .as_ref()
                        .and_then(|e| crate::topopebrep::tool::closed::closed_e(e));
                    let mut on_vclo = false;
                    if let Some(vclo) = e_closed {
                        // `cxx:458`: `onvclo = vG.IsSame(vG);`
                        on_vclo = crate::topo_tools::is_same(&vg, &vg);
                        if has_oo && !on_vclo {
                            // `cxx:460`: `onvclo = vG.IsSame(OOv);`
                            let oo = super::export::fun_ds_getoov(ds, &vg);
                            if let Some(oo) = oo {
                                on_vclo = crate::topo_tools::is_same(&vg, &oo);
                            }
                        }
                        let _ = vclo;
                    }
                    if on_vclo {
                        id_t = true;
                    } else {
                        // `xpu100697`: two EVI interferences are compared by
                        // parameter within `Precision::Parametric(max tolerance
                        // of E)`.
                        let par1 = li[i].read().expect("interference lock").parameter();
                        let (Some(par1), Some(par2)) = (par1, par2) else {
                            j += 1;
                            continue;
                        };
                        let tol_e = match edge.as_ref() {
                            Some(e) => crate::topopebrep::tool::measures::fun_tool_maxtol(e),
                            None => 0.0,
                        };
                        let t = occt_core::precision::Precision::parametric(
                            tol_e,
                            occt_core::precision::CONFUSION,
                        );
                        id_t = (par1 - par2).abs() <= t;
                    }
                }
            }

            if id_t {
                li.remove(j);
            } else {
                j += 1;
            }
        }
        i += 1;
    }
}

