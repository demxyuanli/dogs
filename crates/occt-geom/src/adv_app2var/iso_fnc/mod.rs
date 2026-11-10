//! `AdvApp2Var_Iso::MakeApprox` (`AdvApp2Var_Iso.cxx:106-376`) through
//! `mma2fnc_` (`AdvApp2Var_ApproxF2var.cxx:6624-7281`).
//!
//! The first grid iso is the V frontier built by `InitGrid`
//! (`AdvApp2Var_ApproxAFunc2Var.cxx:271-279`). `FirstNotApprox` returns it
//! (`AdvApp2Var_Framework.cxx:77-106`). `mma2fnc_` runs `mma1nop_`, `mma1fdi_`,
//! `mma1cdi_`, `mma1jak_`, `mma1cnt_`, `mma1fer_`, `mmapcmp_`, and `mmjacan_`.
//! `ChangeIso` stores each frontier that meets tolerance. `Discretise` fills
//! the corner and iso tables. `mma2cdi_` runs `mma2cd3_`, `mma2cd2_`, and
//! `mma2cd1_`. `mma2ds1_` discretises the square. `NumDec` is chosen from
//! `PrefAndRec::Value`. `MakeApprox` cuts with `UpdateInU` and `UpdateInV`,
//! then stores a finished square. `Compute3DErrors` runs when every square
//! is approximated, and `ConvertBS` returns the B-spline.


pub(super) use super::constraints::mma1her;

mod helpers;
mod mma1cdi;
mod mma1fdi;
mod mma1nop;
mod mma2fnc;
mod patches;

pub(crate) use helpers::*;
pub(crate) use mma1cdi::*;
pub(crate) use mma1fdi::*;
pub(crate) use mma1nop::*;
pub(crate) use mma2fnc::*;
pub use patches::*;
