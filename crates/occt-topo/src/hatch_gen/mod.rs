//! `HatchGen` package: the value types shared by `Geom2dHatch`.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/HatchGen/` (OCCT 8.0.0):
//!
//! | OCCT class | Rust |
//! |---|---|
//! | `HatchGen_IntersectionPoint` (`hxx:30-94`, `cxx:20-172`) | [`IntersectionPoint`] |
//! | `HatchGen_IntersectionType` (`hxx:20-28`) | [`IntersectionType`] |
//! | `HatchGen_ErrorStatus` (`hxx:21-28`) | [`ErrorStatus`] |
//! | `HatchGen_PointOnElement` (`hxx:30-70`, `.lxx:18-33`, `cxx:25-292`) | [`PointOnElement`] |
//! | `HatchGen_PointOnHatching` (`hxx:33-88`, `cxx:25-237`) | [`PointOnHatching`] |
//! | `HatchGen_Domain` (`hxx:28-90`, `.lxx:18-126`, `cxx:25-95`) | [`Domain`] |
//!
//! `Dump` (`PointOnElement.cxx:205-292`, `PointOnHatching.cxx:151-237`,
//! `Domain.cxx:63-95`) only writes to `std::cout`; it is not ported, the same
//! rule the rest of the port applies to the `#ifdef OCCT_DEBUG` blocks.

pub mod domain;
pub mod intersection_point;
pub mod point_on_element;
pub mod point_on_hatching;

pub use domain::Domain;
pub use intersection_point::{ErrorStatus, IntersectionPoint, IntersectionType};
pub use point_on_element::PointOnElement;
pub use point_on_hatching::PointOnHatching;
