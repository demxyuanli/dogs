//! The file-static helpers of `ChFi3d_Builder_2.cxx` that the face / face
//! `ComputeData` and the two `StartSol` overloads need.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx`
//! (OCCT 8.0.0), declarations next to their uses.
//!
//! | OCCT | here |
//! |---|---|
//! | `getCurveOnSurface` (`cxx:79-83`) | [`crate::boptools_2d::curve_on_surface_range`] |
//! | `isTangentToArc` (`cxx:192-209`) | [`is_tangent_to_arc::is_tangent_to_arc`] |
//! | `BonVoisin` (`cxx:213-318`) | [`start_sol::bon_voisin`] |
//! | `Projection` (`cxx:320-370`) | [`start_sol::projection`] |
//! | `IsInput` (`cxx:399-471`) | [`search_face::is_input`] |
//! | `IsG1` (`cxx:475-521`) | [`search_face::is_g1`] |
//! | `SearchFaceOnV` (`cxx:532-589`) | [`search_face::search_face_on_v`] |
//! | `ChFi3d_BuildPlane` (`cxx:1113-1144`) | [`build_plane::ch_fi3d_build_plane`] |
//! | `TgtKP` (`cxx:373-398`) | [`tgt_kp::tgt_kp`] |
//! | `ChFi3d_BoxDiag` (`cxx:3282-3297`) | [`box_diag::ch_fi3d_box_diag`] |
//!
//! The member function `ChFi3d_Builder::SearchFace` (`cxx:1523-1627`) is
//! [`crate::chfi3d::builder::ChFi3dBuilder::search_face`], in the `builder`
//! module where the `myVEMap` / `myEFMap` it reads live.

pub mod box_diag;
pub mod build_plane;
pub mod change_transition;
pub mod coupe_par_plan;
pub mod is_tangent_to_arc;
pub mod search_face;
pub mod start_sol;
pub mod tgt_kp;

pub use box_diag::ch_fi3d_box_diag;
pub use build_plane::ch_fi3d_build_plane;
pub use change_transition::change_transition;
pub use coupe_par_plan::{ch_fi3d_coupe_par_plan, CoupeParPlan};
pub use is_tangent_to_arc::is_tangent_to_arc;
pub use search_face::{is_g1, is_input, search_face_on_v};
pub use start_sol::{bon_voisin, projection};
pub use tgt_kp::tgt_kp;
