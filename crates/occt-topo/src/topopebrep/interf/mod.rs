//! `TopOpeBRep` — the face/face and edge/edge intersectors of the S2 wave.
//!
//! Ported so far:
//! - `TopOpeBRep_P2Dstatus` (`TopOpeBRep_P2Dstatus.hxx:22-28`)
//! - `TopOpeBRep_Bipoint` (`TopOpeBRep_Bipoint.hxx:26-43`)
//! - `TopOpeBRep_WPointInter` (`TopOpeBRep_WPointInter.hxx:26-53`)
//! - `TopOpeBRep_Hctxff2d` (`TopOpeBRep_Hctxff2d.hxx:27-68`)
//! - `TopOpeBRep_Point2d` (`TopOpeBRep_Point2d.hxx:40-144`)
//! - `TopOpeBRep_Hctxee2d` (`TopOpeBRep_Hctxee2d.hxx:27-58`)
//! - `TopOpeBRep_PointGeomTool` (`TopOpeBRep_PointGeomTool.hxx:30-44`)
//! - `TopOpeBRep_PointClassifier` (`TopOpeBRep_PointClassifier.hxx`)
//! - `TopOpeBRep_VPointInter` (`TopOpeBRep_VPointInter.hxx:32-194`)
//!
//! Not ported yet: `TopOpeBRep_GeomTool` and the intersectors/fillers
//! themselves.

pub mod bipoint;
pub mod hctxee2d;
pub mod hctxff2d;
pub mod p2d_status;
pub mod point2d;
pub mod point_classifier;
pub mod point_geom_tool;
pub mod v_point_inter;
pub mod v_point_inter_classifier;
pub mod w_point_inter;

pub use bipoint::Bipoint;
pub use hctxee2d::Hctxee2d;
pub use hctxff2d::Hctxff2d;
pub use p2d_status::P2DStatus;
pub use point2d::Point2d;
pub use point_classifier::PointClassifier;
pub use point_geom_tool::PointGeomTool;
pub use v_point_inter::VPointInter;
pub use v_point_inter_classifier::VPointInterClassifier;
pub use w_point_inter::WPointInter;
