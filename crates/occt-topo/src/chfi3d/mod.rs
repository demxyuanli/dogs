//! `TKFillet/ChFi3d` — the fillet/chamfer building package.
//!
//! Source: `src/ModelingAlgorithms/TKFillet/ChFi3d/`.
//!
//! | OCCT | here |
//! |---|---|
//! | `ChFi3d` (`ChFi3d.hxx:33-90`, `ChFi3d.cxx:39-647`) | [`connect_type`], [`tangent_faces`], [`side`] |
//! | `ChFi3d_Coefficient` (`ChFi3d_Builder_0.cxx:1451-1463`) | [`coefficient::chfi3d_coefficient`] |
//! | `ChFi3d_FilletShape` (`ChFi3d_FilletShape.hxx:28-33`) | [`fillet_shape::FilletShape`] |
//! | `ChFi3d_SearchSing` (`ChFi3d_SearchSing.hxx:22-52`, `.cxx:27-57`) | [`search_sing::SearchSing`] |
//! | the `ChFi3d_Builder_0.cxx` free functions | [`builder_helpers`] |
//! | the `ChFi3d_Builder_6.cxx` file statics | [`builder_6_helpers`] |
//! | the `ChFi3d_Builder_2.cxx` file statics (`IsInput`, `IsG1`, `SearchFaceOnV`) | [`builder_2_helpers`] |
//! | `ChFi3d_Builder` (`ChFi3d_Builder.hxx:77-848`) | [`builder::ChFi3dBuilder`] |

pub mod builder;
pub mod builder_2_helpers;
pub mod builder_6_helpers;
pub mod builder_helpers;
pub mod coefficient;
pub mod connect_type;
pub mod fillet_shape;
pub mod search_sing;
pub mod side;
pub mod tangent_faces;

pub use builder::ChFi3dBuilder;
pub use coefficient::chfi3d_coefficient;
pub use connect_type::{correct2d_point, define_connect_type};
pub use fillet_shape::FilletShape;
pub use search_sing::SearchSing;
pub use side::{concave_side, next_side, next_side_single, same_side};
pub use tangent_faces::is_tangent_faces;
