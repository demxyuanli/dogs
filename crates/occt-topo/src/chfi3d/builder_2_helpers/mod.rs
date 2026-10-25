//! The file-static helpers of `ChFi3d_Builder_2.cxx` that the face / face
//! `ComputeData` needs.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx`
//! (OCCT 8.0.0), declarations next to their uses.
//!
//! | OCCT | here |
//! |---|---|
//! | `IsInput` (`cxx:399-471`) | [`search_face::is_input`] |
//! | `IsG1` (`cxx:475-521`) | [`search_face::is_g1`] |
//! | `SearchFaceOnV` (`cxx:532-589`) | [`search_face::search_face_on_v`] |
//!
//! The member function `ChFi3d_Builder::SearchFace` (`cxx:1523-1627`) is
//! [`crate::chfi3d::builder::ChFi3dBuilder::search_face`], in the `builder`
//! module where the `myVEMap` / `myEFMap` it reads live.

pub mod search_face;
pub mod start_sol;

pub use search_face::{is_g1, is_input, search_face_on_v};
pub use start_sol::{bon_voisin, projection};
