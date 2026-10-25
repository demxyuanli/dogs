//! `ChFiDS` — the data structures the fillet/chamfer bisector walks on.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/` (OCCT 8.0.0): 34 files,
//! 18 classes plus 5 enums.
//!
//! The package is the floor of the `TKFillet` wave (`ChFi3d` builds on
//! `Spine` / `Stripe` / `SurfData`, and `BRepFilletAPI` exposes them), so it is
//! ported before any of `ChFi3d` / `BRepBlend` / `BlendFunc`.
//!
//! Ported so far (batch 1 — the value layer):
//! - [`enums`]: `ChFiDS_State`, `ChFiDS_TypeOfConcavity`, `ChFiDS_ChamfMethod`,
//!   `ChFiDS_ChamfMode`, `ChFiDS_ErrorStatus`,
//! - [`circ_section`]: `ChFiDS_CircSection` plus the
//!   `NCollection_HArray1<ChFiDS_CircSection>` shape that `SurfData::Simul`
//!   stores,
//! - [`common_point`]: `ChFiDS_CommonPoint`,
//! - [`face_interference`]: `ChFiDS_FaceInterference`,
//! - [`regul`]: `ChFiDS_Regul`,
//! - [`map`]: `ChFiDS_Map` (`IndexedDataMapOfShapeListOfShape`),
//! - [`surf_data`]: `ChFiDS_SurfData`.
//!
//! Ported so far (batch 2 — the guideline and the bands):
//! - [`el_spine`]: `ChFiDS_ElSpine`,
//! - [`spine`]: `ChFiDS_Spine` plus the `SpineHandle` alias,
//! - [`stripe`]: `ChFiDS_Stripe` plus the `StripeHandle` alias,
//! - [`stripe_map`]: `ChFiDS_StripeMap` (`IndexedDataMapOfVertexListOfStripe`).
//!
//! Ported so far (batch 3 — the two specialised guidelines):
//! - [`fil_spine`]: `ChFiDS_FilSpine` (the radius law) on top of
//!   [`spine`], with its `FilSpineData` in `Spine`,
//! - [`chamf_spine`]: `ChFiDS_ChamfSpine` (the section distances) on top of
//!   [`spine`], with its `ChamfSpineData` in `Spine`.
//!
//! UNPORTED, with their OCCT location:
//! - the `ChFiDS_ElSpine` accessors listed in [`el_spine`]'s module docs
//!   (`Trim`, `SetOrigin`, `Bezier`, `BSpline`).
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//! - OCCT throws (`Standard_DomainError`, `Standard_Failure`, `Standard_NullObject`,
//!   `Standard_OutOfRange`, `StdFail_NoSuchObject`) from the accessors that
//!   guard a "not set" / "not found" state; the port answers `Option` instead,
//! - `occ::handle<T>` becomes `Arc<T>` (shared, immutable once handed out),
//!   `Option<Arc<T>>` (null handle) or, for the classes OCCT mutates through the
//!   handle, `Arc<RwLock<T>>` (`SpineHandle`, `StripeHandle`) — the same shape
//!   `topopebrep::ds::InterferenceHandle` uses,
//! - `NCollection_List<T>` / `NCollection_Sequence<T>` become `Vec<T>`,
//! - `NCollection_Array1` keeps OCCT's 1-based indexing where the OCCT code
//!   reads it 1-based (see [`circ_section::CircSectionArray`], and the
//!   `ChFiDS_Spine` abscissa accessors),
//! - `ChFiDS_FilSpine` / `ChFiDS_ChamfSpine` derive from `ChFiDS_Spine` in
//!   OCCT; a `RwLock<Spine>` cannot hold a subclass, so their extra state lives
//!   in `Spine` behind the [`spine::SpineClass`] discriminant, and their
//!   methods are `impl Spine` blocks in [`fil_spine`] / [`chamf_spine`].

pub mod chamf_spine;
pub mod circ_section;
pub mod common_point;
pub mod el_spine;
pub mod enums;
pub mod face_interference;
pub mod fil_spine;
pub mod map;
pub mod regul;
pub mod spine;
pub mod stripe;
pub mod stripe_map;
pub mod surf_data;

pub use chamf_spine::ChamfSpineData;
pub use circ_section::{CircSection, CircSectionArray};
pub use common_point::CommonPoint;
pub use el_spine::ElSpine;
pub use enums::{ChamfMethod, ChamfMode, ErrorStatus, State, TypeOfConcavity};
pub use face_interference::FaceInterference;
pub use fil_spine::FilSpineData;
pub use map::ChFiDSMap;
pub use regul::Regul;
pub use spine::{Spine, SpineClass, SpineHandle};
pub use stripe::{Stripe, StripeHandle};
pub use stripe_map::StripeMap;
pub use surf_data::SurfData;
