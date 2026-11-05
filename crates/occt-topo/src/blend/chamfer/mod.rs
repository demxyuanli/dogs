//! The chamfer solvers of `TKFillet`'s `BlendFunc` package.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/` (OCCT 8.0.0).
//!
//! Ported:
//! - [`gen_chamfer`]: `BlendFunc_GenChamfer` (`hxx:31-149`, `.cxx:29-300`), the
//!   deferred base between two surfaces,
//! - [`const_throat`]: `BlendFunc_ConstThroat` (`hxx:27-110`, `.cxx:27-281`),
//!   the symmetric chamfer with a constant throat,
//! - [`const_throat_with_penetration`]:
//!   `BlendFunc_ConstThroatWithPenetration` (`hxx:28-100`, `.cxx:27-185`), the
//!   right-angled section,
//! - [`gen_chamf_inv`]: `BlendFunc_GenChamfInv` (`hxx:31-107`, `.cxx:23-122`),
//!   the deferred inverse base on a curve on surface,
//! - [`const_throat_inv`]: `BlendFunc_ConstThroatInv` (`hxx:26-76`,
//!   `.cxx:22-137`),
//! - [`const_throat_with_penetration_inv`]:
//!   `BlendFunc_ConstThroatWithPenetrationInv` (`hxx:24-49`, `.cxx:22-199`).
//! - [`chamfer`]: `BlendFunc_Chamfer` (`hxx:28-70`, `.cxx:24-114`), the
//!   two-distance chamfer built on two [`corde`] systems,
//! - [`chamf_inv`]: `BlendFunc_ChamfInv` (`hxx:29-60`, `.cxx:24-230`), its
//!   inverse on a curve on surface,
//! - [`ch_asym`]: `BlendFunc_ChAsym` (`hxx:34-220`, `.cxx:32-766`), the
//!   distance + angle chamfer,
//! - [`ch_asym_inv`]: `BlendFunc_ChAsymInv` (`hxx:24-84`, `.cxx:24-404`), its
//!   inverse on a curve on surface,
//! - [`corde`]: `BlendFunc_Corde` (`hxx:28-100`, `.cxx:30-203`), the
//!   single-surface distance helper.
//!
//! [`shared`] holds the two `Section` bodies OCCT writes once in
//! `BlendFunc_GenChamfer` and reaches through the virtual `Value` /
//! `Derivatives` / `IsSolution` / `PointOnS*` / `Tangent*`.
//!
//! UNPORTED here: none.

mod ch_asym;
mod ch_asym_inv;
mod chamf_inv;
mod chamfer;
mod const_throat;
mod const_throat_inv;
mod const_throat_with_penetration;
mod const_throat_with_penetration_inv;
mod corde;
mod gen_chamf_inv;
mod gen_chamfer;
mod shared;

pub use ch_asym::ChAsym;
pub use ch_asym_inv::ChAsymInv;
pub use chamf_inv::ChamfInv;
pub use chamfer::Chamfer;
pub use const_throat::ConstThroat;
pub use const_throat_inv::ConstThroatInv;
pub use const_throat_with_penetration::ConstThroatWithPenetration;
pub use const_throat_with_penetration_inv::ConstThroatWithPenetrationInv;
pub use corde::Corde;
pub use gen_chamf_inv::GenChamfInv;
pub use gen_chamfer::GenChamfer;
// `BlendFunc_GenChamfer::Section` (`cxx:104-122`), shared with
// `BlendFunc_ChAsym::Section` — see the `ChFi3d_ChBuilder::SimulSurf` section
// loop in `chfi3d/builder/ch_builder/simul_surf.rs`.
pub(crate) use gen_chamfer::gen_chamfer_section;
