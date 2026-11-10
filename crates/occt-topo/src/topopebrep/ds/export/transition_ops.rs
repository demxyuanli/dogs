use crate::shape::TopoShape;

use super::super::data_structure::DataStructure;
use super::super::interference_iterator::InterferenceList;
use super::super::transition::Transition;


/// `FDS_SetT(T, T0)` (`cxx:61-68`).
pub fn fds_set_t(t: &mut Transition, t0: &Transition) {
    let stb = t0.before_state();
    let sta = t0.after_state();
    t.set_state_before(stb);
    t.set_state_after(sta);
}

/// `FDS_hasUNK(T)` (`cxx:70-77`): either side is `UNKNOWN`.
pub fn fds_has_unk(t: &Transition) -> bool {
    t.before_state() == crate::abs::State::Unknown || t.after_state() == crate::abs::State::Unknown
}

/// `FDS_copy(LI, LII)` (`cxx:79-88`): append every element of `LI` to `LII`.
pub fn fds_copy(li: &InterferenceList, lii: &mut InterferenceList) {
    lii.extend(li.iter().cloned());
}

/// `FDS_copy(LI, LII)` for `List<TopoDS_Shape>` (`cxx:91-99`).
pub fn fds_copy_shapes(li: &[TopoShape], lii: &mut Vec<TopoShape>) {
    lii.extend(li.iter().cloned());
}

/// `FDS_assign(LI, LII)` (`cxx:101-108`): clear then copy.
pub fn fds_assign(li: &InterferenceList, lii: &mut InterferenceList) {
    lii.clear();
    fds_copy(li, lii);
}

/// `FDS_assign(LI, LII)` for `List<TopoDS_Shape>` (`cxx:110-117`).
pub fn fds_assign_shapes(li: &[TopoShape], lii: &mut Vec<TopoShape>) {
    lii.clear();
    fds_copy_shapes(li, lii);
}

/// `FUN_ds_samRk(BDS, Rk, LI, LIsrk)` (`cxx:119-142`): move the elements of
/// `LI` whose `AncestorRank` is `Rk` into `LIsrk` (which is cleared first).
pub fn fun_ds_sam_rk(bds: &DataStructure, rk: i32, li: &mut Vec<TopoShape>, lisrk: &mut Vec<TopoShape>) {
    lisrk.clear();
    let mut i = 0;
    while i < li.len() {
        let crk = bds.ancestor_rank_of(&li[i]);
        if crk == rk {
            let moved = li.remove(i);
            lisrk.push(moved);
        } else {
            i += 1;
        }
    }
}

