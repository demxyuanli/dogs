//! `TopOpeBRepDS_EXPORT` — the free functions of the classic boolean engine
//! (`FDS_*` / `FUN_ds_*`).
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.hxx:24-91`,
//! `TopOpeBRepDS_EXPORT.cxx:1-3952`.
//!
//! Ported here: the pure list / interference accessors and the same-domain
//! queries that need no geometry evaluation (`cxx:61-360`, `396-462`,
//! `464-495`, `647-712`, `1874-1995`). The `FUN_ds_completefor*` /
//! `FUN_ds_PURGEforE9` / `FUN_ds_PointToVertex` / `FUN_ds_redusamsha`
//! family lives in [`super::complete_ds`], keyed by its own `.cxx` lines.
//!
//! Still UNPORTED in this file's source range: `FDS_HasSameDomain3d`'s
//! geometry branch is in [`fds_has_same_domain3d`]; `FDS_SIisGIofIofSBAofTofI`
//! is [`fds_si_is_gi_of_i_of_sba_of_t_of_i`]. `FUN_ds_hasFEI` (`:3279`),
//! `FUN_ds_ONesd` (`:3308`), `FUN_ds_unkeepEVIonGb1` (`:3364`),
//! `FDS_LOIinfsup`'s 2d half (`:3427`), `FDS_parbefaft` (`:3571`),
//! `FDS_stateEwithF2d` (`:3598`) and `FUN_ds_FEIGb1TO0` (`:3647`) are the
//! remaining entry points of the file; `FUN_ds_unkeepEVIonGb1` is only
//! called from commented-out code in `TopOpeBRep_DSFiller.cxx:769`.
//!
//! OCCT's `bool(a, b)` output pairs become `Option`/tuple returns, and the
//! `bool FDS_data(it, I, ...)` overload becomes [`fds_data_is_unknown`].


mod data_access;
mod eis_se2;
mod fill_sdm;
mod loi;
mod parameter;
mod same_domain;
mod shape_queries;
mod si_gi;
mod transition_ops;

pub use data_access::*;
pub use eis_se2::*;
pub use fill_sdm::*;
pub use loi::*;
pub use parameter::*;
pub use same_domain::*;
pub use shape_queries::*;
pub use si_gi::*;
pub use transition_ops::*;
// UNPORTED in `TopOpeBRepDS_EXPORT.cxx` (all call sites outside this crate's
// reached packages, or dead in OCCT itself):
// - `FUN_ds_hasFEI` (`:3279-3306`) and `FUN_ds_ONesd` (`:3308-3345`): the only
//   callers are `TopOpeBRepBuild_GridFF.cxx:227/246/281` and
//   `TopOpeBRepBuild_BuilderON.cxx:672`, i.e. the `TopOpeBRepBuild` package.
// - `FUN_ds_unkeepEVIonGb1` (`:3364-3425`) and `FUN_ds_FEIGb1TO0`
//   (`:3647-3774`): their only call sites are commented out in OCCT
//   (`TopOpeBRep_DSFiller.cxx:768`, `TopOpeBRepDS_Reducer.cxx:44`), and the
//   first one's body is itself inside `/* */`.
//
// The `FUN_ds_completefor*` / `FUN_ds_PURGEforE9` / `FUN_ds_redusamsha` /
// `FUN_ds_PointToVertex` / `FUN_ds_complete1dForSESDM` family (`:715-3952`)
// lives in [`super::complete_ds`]; the 1d/2d reduction helpers of
// `TopOpeBRepDS_redu.cxx` live in [`super::redu`].

