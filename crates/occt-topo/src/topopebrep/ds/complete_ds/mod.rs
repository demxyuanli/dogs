//! `TopOpeBRepDS_EXPORT` — the `CompleteDS` completion passes.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:715-3277` and
//! `:3777-3949`; the 1d/2d reduction helper of
//! `FUN_ds_completeforSE8` is `TopOpeBRepDS_redu.cxx:35-353` ([`super::redu`]).
//!
//! `TopOpeBRep_DSFiller::CompleteDS` (`TopOpeBRep_DSFiller.cxx:754-787`) runs
//! the passes in this order, and [`run_complete_ds`] keeps that order:
//!
//! ```text
//! 770 FUN_ds_PointToVertex        point_to_vertex
//! 771 FUN_ds_redusamsha           redusamsha
//! 772 FUN_ds_PURGEforE9           purge_for_e9
//! 773 FUN_ds_completeforSE8       se8_se9
//! 774 FUN_ds_completeforSE1       se1
//! 775 FUN_ds_completeforSE2       se2
//! 776 FUN_ds_completeforSE3       se3
//! 777 FUN_ds_completeforSE4       se4_se5
//! 778 FUN_ds_completeforSE5       se4_se5
//! 779 FUN_ds_completeforSE6       se6_e7
//!     FUN_ds_completeforSE7       -- commented out at cxx:780, NYI
//! 781 FUN_ds_completeforE7         se6_e7
//! 782 FUN_ds_completeforSE9        se8_se9
//! 783 FUN_ds_complete1dForSESDM    complete_1d_se_sdm
//! ```
//!
//! `FUN_ds_unkeepEVIonGb1` (`cxx:769`) is commented out in OCCT and is not
//! ported. Every pass mutates `HDS::ChangeDS()`; OCCT's
//! `HDS->StoreInterference(I, S)` is [`HDataStructure::store_interference_on_shape`].

mod complete_1d_se_sdm;
mod has_i2d;
mod point_to_vertex;
mod purge_for_e9;
mod redusamsha;
mod se1;
mod se2;
mod se3;
mod se4_se5;
mod se6_e7;
mod se8_se9;
mod share_g;

pub use complete_1d_se_sdm::fun_ds_complete1d_for_se_sdm;
pub use has_i2d::fun_ds_has_i2d;
pub use point_to_vertex::fun_ds_point_to_vertex;
pub use purge_for_e9::{fun_ds_has_sdm_ancestor_faces, fun_ds_purge_for_e9};
pub use redusamsha::{fun_ds_redusamsha, fun_redusamsha_on_e};
pub use se1::fun_ds_complete_for_se1;
pub use se2::fun_ds_complete_for_se2;
pub use se3::{fun_ds_complete_for_se3, fun_se3_transition};
pub use se4_se5::{fun_ds_complete_for_se4, fun_ds_complete_for_se5};
pub use se6_e7::{fun_ds_complete_for_e7, fun_ds_complete_for_se6};
pub use se8_se9::{fun_ds_complete_for_se8, fun_ds_complete_for_se9};
pub use share_g::{
    fun_ds_mk_ton_fsdm, fun_ds_ori_e_in_f, fun_ds_share_g, ORI_E_IN_F_CLOSE_OPPO,
    ORI_E_IN_F_CLOSE_SAME, ORI_E_IN_F_FOR_REV_OPPO, ORI_E_IN_F_ON_OPPO_SHA,
    ORI_E_IN_F_ON_SAME_SHA, ORI_E_IN_F_UNKNOWN,
};

use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::h_data_structure::HDataStructure;

/// `TopOpeBRep_DSFiller::CompleteDS` (`TopOpeBRep_DSFiller.cxx:754-787`).
///
/// `S1`/`S2` reach the passes only through the DS, which the caller has already
/// filled; the two null-shape guards of `cxx:759-764` stay in the caller (see
/// `DsFiller::complete_ds`).
pub fn run_complete_ds(hds: &mut HDataStructure, connex: &Connexity) {
    fun_ds_point_to_vertex(hds);
    fun_ds_redusamsha(hds);
    fun_ds_purge_for_e9(hds, connex);
    fun_ds_complete_for_se8(hds);
    fun_ds_complete_for_se1(hds, connex);
    fun_ds_complete_for_se2(hds, connex);
    fun_ds_complete_for_se3(hds);
    fun_ds_complete_for_se4(hds);
    fun_ds_complete_for_se5(hds);
    fun_ds_complete_for_se6(hds);
    fun_ds_complete_for_e7(hds);
    fun_ds_complete_for_se9(hds);
    fun_ds_complete1d_for_se_sdm(hds);
}
