//! `TopOpeBRepBuild_Builder` — the KPart special cases.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_KPart.cxx`
//! and `TopOpeBRepBuild_kpresu.hxx` (OCCT 8.0.0).
//!
//! `FindIsKPart` (`cxx:111-166`) runs the five predicates in OCCT's order and
//! stores the first that answers yes in `myIsKPart`; `MergeKPart` (`cxx:185-208`)
//! dispatches on that value and always ends with `End()`.
//!
//! Where each part lives:
//! - `lists.rs`: `KPlhg` / `KPlhsd` / `KPls` / `KPContains` / `KPSameDomain` and
//!   the `KPclasSS` / `KPclassF` / `KPclassFF` / `KPiskoleFF` classifiers.
//! - `predicates.rs`: `KPisdisj` / `KPisfafa` / `KPissoso` / `KPiskoletge` and
//!   their `*sh` sub-predicates.
//! - `merge_same.rs`: `MergeKPartisfafa` (`cxx:1034-1101`) and
//!   `MergeKPartissoso` (`cxx:1105-1192`).
//! - `merge_disj.rs`: `MergeKPartisdisj` (`cxx:840-1030`) with
//!   `LocalKPisdisjanalyse`, `BuildNewSolid` and the `disjPerform*` helpers.
//! - `merge_koletge.rs`: `MergeKPartiskoletge` (`cxx:578-836`) with
//!   `FUN_sortplcy` / `FUN_addf` / `FUN_comparekoletgesh` / `FUN_makefaces` /
//!   `FUN_rebuildfc`.
//! - `iskole.rs`: `KPiskolesh` / `KPiskole` / `KPiskoleanalyse`
//!   (`kpkole.cxx:1085-1157`, `:576-732`, `:734-843`).
//! - `merge_iskole.rs`: `MergeKPartiskole` (`kpkole.cxx:56-572`) with
//!   `KPmakeface` (`:945-1031`) and `FUNKP_KPmakefaces` (`:845-941`).
//!
//! UNPORTED, with its OCCT location:
//! - the `GMergeFaces` / `GMergeSolids` calls of `MergeKPartisfafa`
//!   (`KPart.cxx:1086`) / `MergeKPartissoso` (`:1172`): those live in
//!   `TopOpeBRepBuild_GridFF.cxx:335-463` / `GridSS.cxx:282-349` and are not
//!   translated.

use crate::abs::State;

use super::Builder;

mod lists;
mod iskole;
mod merge_disj;
mod merge_iskole;
mod merge_koletge;
mod merge_same;
mod predicates;

/// `TopOpeBRepBuild_kpresu.hxx`: the values `myIsKPart` takes.
pub const KP_ISKOLE: i32 = 1;
pub const KP_ISDISJ: i32 = 2;
pub const KP_ISAFA: i32 = 3;
pub const KP_ISSOSO: i32 = 4;
pub const KP_ISKOLETGE: i32 = 5;

/// `RESUNDEF` (`kpresu.hxx:20`).
pub const RESUNDEF: i32 = -100;
/// `RESNEWSHA2` (`kpresu.hxx:22`).
pub const RESNEWSHA2: i32 = -12;
/// `RESNEWSHA1` (`kpresu.hxx:23`).
pub const RESNEWSHA1: i32 = -11;
/// `RESNEWSOL` (`kpresu.hxx:25`).
pub const RESNEWSOL: i32 = -2;
/// `RESNEWSHE` (`kpresu.hxx:26`).
pub const RESNEWSHE: i32 = -1;
/// `RESNULL` (`kpresu.hxx:27`).
pub const RESNULL: i32 = 0;
/// `RESSHAPE1` (`kpresu.hxx:28`).
pub const RESSHAPE1: i32 = 1;
/// `RESSHAPE2` (`kpresu.hxx:29`).
pub const RESSHAPE2: i32 = 2;
/// `RESSHAPE12` (`kpresu.hxx:30`).
pub const RESSHAPE12: i32 = 3;

/// `RESFACE1` (`kpresu.hxx:30`).
pub const RESFACE1: i32 = 11;
/// `RESFACE2` (`kpresu.hxx:31`).
pub const RESFACE2: i32 = 12;

/// `SHEUNDEF` (`kpresu.hxx:32`).
pub const SHEUNDEF: i32 = -100;
/// `SHEAUCU` (`kpresu.hxx:33`).
pub const SHEAUCU: i32 = -1;
/// `SHECLASAUTR` (`kpresu.hxx:36`).
pub const SHECLASAUTR: i32 = 2;
/// `SHEGARDCOUR` (`kpresu.hxx:37`).
pub const SHEGARDCOUR: i32 = 4;
/// `SHEGARDAUTR` (`kpresu.hxx:38`).
pub const SHEGARDAUTR: i32 = 5;
/// `SHEGARDTOUS` (`kpresu.hxx:39`).
pub const SHEGARDTOUS: i32 = 6;

impl Builder {
    /// `FindIsKPart()` (`KPart.cxx:111-166`).
    ///
    /// The predicates run in OCCT's order — `isfafa`, `isdisj`, `iskole`,
    /// `iskoletge`, `issoso` — and the first answer wins. `KPiskole` is
    /// UNPORTED (`kpkole.cxx:576-732`), so its arm always answers 0.
    pub fn find_is_k_part(&mut self) -> i32 {
        self.kp_clear_maps();

        // cxx:124-130
        if self.kp_isfafa() != 0 {
            self.my_is_k_part = KP_ISAFA;
            return self.kp_return(self.my_is_k_part);
        }

        // cxx:132-138
        if self.kp_isdisj() != 0 {
            self.my_is_k_part = KP_ISDISJ;
            return self.kp_return(self.my_is_k_part);
        }

        // cxx:139-145
        if self.kp_iskole() != 0 {
            self.my_is_k_part = KP_ISKOLE;
            return self.kp_return(self.my_is_k_part);
        }

        // cxx:147-153
        if self.kp_iskoletge() != 0 {
            self.my_is_k_part = KP_ISKOLETGE;
            return self.kp_return(self.my_is_k_part);
        }

        // cxx:155-162
        if self.kp_issoso() != 0 {
            self.my_is_k_part = KP_ISSOSO;
            return self.kp_return(self.my_is_k_part);
        }

        // cxx:164-165
        self.my_is_k_part = 0;
        self.kp_return(self.my_is_k_part)
    }

    /// `MergeKPart(TB1, TB2)` (`KPart.cxx:176-183`).
    pub fn merge_k_part_states(&mut self, tb1: State, tb2: State) {
        self.my_state1 = tb1;
        self.my_state2 = tb2;
        self.merge_k_part();
    }

    /// `MergeKPart()` (`KPart.cxx:185-208`): the five cases, then `End()`.
    ///
    /// `End()` (`TopOpeBRepBuild_End.cxx:45-410`) is ported in
    /// `builder/end.rs`; its `FuseEdges` block is UNPORTED, see there.
    pub fn merge_k_part(&mut self) {
        if self.my_is_k_part == KP_ISKOLE {
            self.merge_k_part_iskole();
        } else if self.my_is_k_part == KP_ISKOLETGE {
            self.merge_k_part_iskoletge();
        } else if self.my_is_k_part == KP_ISDISJ {
            self.merge_k_part_isdisj();
        } else if self.my_is_k_part == KP_ISAFA {
            self.merge_k_part_isfafa();
        } else if self.my_is_k_part == KP_ISSOSO {
            self.merge_k_part_issoso();
        }
        // cxx:208
        self.end();
    }

    /// `KPClearMaps()` (`KPart.cxx:1453-1456`): `myKPMAPf1f2.Clear()`.
    ///
    /// `FindIsKPart` starts a fresh decision, and `KPiskole` / `KPiskoletge`
    /// fill the map as they walk the same-domain faces.
    fn kp_clear_maps(&mut self) {
        self.my_kp_map_f1f2.clear();
    }

    /// `KPreturn(b)` (`KPart.cxx:2368-2383`): the trace-only wrapper.
    fn kp_return(&self, b: i32) -> i32 {
        b
    }
}
