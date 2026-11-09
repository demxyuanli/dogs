//! The TKOffset `Draft` package: the draft-angle modification and the three
//! per-sub-shape info records it propagates through.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKOffset/Draft/`):
//! `Draft_ErrorStatus.hxx:20-28`, `Draft_VertexInfo.{hxx,cxx}`,
//! `Draft_EdgeInfo.{hxx,cxx}`, `Draft_FaceInfo.{hxx,cxx}`, `Draft.{hxx,cxx}`,
//! `Draft_Modification.{hxx,cxx}` and `Draft_Modification_1.cxx`.
//!
//! | OCCT | file | status |
//! |---|---|---|
//! | `Draft_ErrorStatus` (`hxx:20-28`) | `error_status.rs` | ported |
//! | `Draft_VertexInfo` (`cxx:27-117`) | `vertex_info.rs` | ported |
//! | `Draft_EdgeInfo` (`cxx:26-173`) | `edge_info.rs` | ported |
//! | `Draft_FaceInfo` (`cxx:24-117`) | `face_info.rs` | ported |
//! | `Draft::Angle` (`cxx:34-100`) | `angle.rs` | ported |
//! | `Draft_Modification` (`Draft_Modification.cxx:54-510`) | `modification/mod.rs` | ported |
//! | `InternalAdd` (`_1.cxx:116-420`) | `modification/internal_add.rs` | ported |
//! | `Propagate` (`_1.cxx:424-705`) | `modification/propagate.rs` | ported |
//! | private `NewSurface` / `FindRotation` (`_1.cxx:1780-1960`, `:2403-2453`) | `modification/new_surface.rs` | ported |
//! | private `NewCurve` (`_1.cxx:1963-2044`) | `modification/new_curve.rs` | ported |
//! | `Perform` (`_1.cxx:708-1777`) | `modification/perform/` | ported |
//! | `Choose` / `Parameter` / `SmartParameter` (`_1.cxx:2047-2379`) | `modification/perform/` | ported |
//! | `Orientation(S, F)` (`_1.cxx:2383-2398`) | `modification/mod.rs` (`shape_face_orientation`) | ported |
//!
//! The `Draft` package itself is complete; `BRepOffsetAPI_DraftAngle` (the
//! `BRepBuilderAPI_MakeShape` facade) is
//! [`crate::brep_offset_api_draft_angle`].

mod angle;
mod edge_info;
mod error_status;
mod face_info;
mod modification;
mod vertex_info;

pub use angle::*;
pub use edge_info::*;
pub use error_status::*;
pub use face_info::*;
pub use modification::*;
pub use vertex_info::*;
