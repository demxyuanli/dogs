//! The file-static helpers of `ChFi3d_Builder_6.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx`
//! (OCCT 8.0.0). They are `static` functions of that translation unit, not
//! members of `ChFi3d_Builder`, so they live here rather than in
//! [`crate::chfi3d::builder`]:
//!
//! | OCCT | here |
//! |---|---|
//! | `SearchIndex` (`cxx:94-103`) | [`search_index::search_index`] |
//! | `nbedconnex` (`cxx:107-130`) | [`obstacle::nbedconnex`] |
//! | `IsVois` (`cxx:132-193`) | [`obstacle::is_vois`] |
//! | `IsObst` (`cxx:195-205`) | [`obstacle::is_obst`] |
//! | `CompBlendPoint` (`cxx:304-339`) | [`comp_blend_point::comp_blend_point`] |
//!
//! `CompParam` (`cxx:209-302`) and `UpdateLine` (`cxx:344-387`) are not here:
//! they are `ChFi3d_Builder` methods and live in
//! [`crate::chfi3d::builder_helpers`] / [`crate::chfi3d::builder`].

pub mod comp_blend_point;
pub mod obstacle;
pub mod search_index;

pub use comp_blend_point::comp_blend_point;
pub use obstacle::{is_obst, is_vois, nbedconnex};
pub use search_index::search_index;
