//! `HLRTopoBRep` package (`ModelingAlgorithms/TKHLR/HLRTopoBRep`): the bridge
//! between a `TopoDS_Shape` and the exact HLR data structure. The outliner
//! collects the shape's faces/edges/vertices, the iso-line liner adds the
//! isoparametric curves of the faces, and the DS filler fills `HLRBRep_Data`
//! from the `Contap_Contour` results.
//!
//! Ported: `HLRTopoBRep_VData`, `HLRTopoBRep_FaceData`, `HLRTopoBRep_Data`,
//! `HLRTopoBRep_OutLiner`, `HLRTopoBRep_FaceIsoLiner`, `HLRTopoBRep_DSFiller`.

pub mod data;
pub mod ds_filler;
pub mod face_data;
pub mod face_iso_liner;
pub mod out_liner;
pub mod v_data;

pub use data::*;
pub use face_data::*;
pub use out_liner::*;
pub use v_data::*;
