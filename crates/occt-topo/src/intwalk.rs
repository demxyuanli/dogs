//! `IntWalk` — parametric surface-surface marching.
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntWalk/` (`IntWalk_PWalking`,
//! `IntWalk_TheInt2S`, `IntImp_ZerParFunc`, `IntImp_ComputeTangence`).

#[path = "intwalk_iso.rs"]
mod iso;
#[path = "intwalk_func.rs"]
mod func;
#[path = "intwalk_int2s.rs"]
mod int2s;
#[path = "intwalk_pwalking.rs"]
mod pwalking;
#[path = "intwalk_pwalking_deflect.rs"]
mod pwalking_deflect;
#[path = "intwalk_pwalking_perform.rs"]
mod pwalking_perform;

pub use iso::{compute_tangence, ConstIso, StatusDeflection};
pub use int2s::TheInt2S;
pub use pwalking::PWalking;
