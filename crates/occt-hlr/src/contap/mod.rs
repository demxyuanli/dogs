//! `Contap` - the contour (silhouette / draft / perspective outline)
//! computation of a surface, i.e. the `Contap` package of `TKHLR`.
//!
//! Ported so far:
//! - [`ContapIType`] / [`ContapTFunction`] (`Contap_IType.hxx`,
//!   `Contap_TFunction.hxx`).
//! - [`ContapPoint`] (`Contap_Point`).
//! - [`ContapLine`] (`Contap_Line`).
//! - [`ContapSurfProps`] functions (`Contap_SurfProps`).
//! - [`ContapHContTool`](sample) (`Contap_HContTool`, `Contap_HCurve2dTool`).
//! - [`ContapArcFunction`] (`Contap_ArcFunction`).
//! - [`ContapSurfFunction`] (`Contap_SurfFunction`).
//! - [`ContapContAna`] (`Contap_ContAna`).
//! - [`ContapThePathPointOfTheSearch`] / [`ContapTheSegmentOfTheSearch`] /
//!   [`ContapTheSearch`] (`IntStart_SearchOnBoundaries.gxx` instantiated by
//!   `Contap_TheSearch_0.cxx`).
//! - [`ContapIWLineOfTheIWalking`] (`Contap_TheIWLineOfTheIWalking`) and
//!   [`ContapTheIWalking`] (`IntWalk_IWalking.gxx` instantiated by
//!   `Contap_TheIWalking_0.cxx`).
//! - [`ContapContour`] (`Contap_Contour.cxx:1-2389`): the walking arm
//!   (`contour/perform_walking.rs`), the analytical arm
//!   (`contour/perform_ana.rs`) and the `LineConstructor` segmentation.

pub mod arc_function;
pub mod cont_ana;
pub mod contour;
pub mod iw_line;
pub mod iwalking;
pub mod line;
pub mod path_point;
pub mod point;
pub mod sample;
pub mod search_inside;
pub mod segment;
pub mod surf_function;
pub mod surf_props;
pub mod the_search;
pub mod tools;
pub mod types;

pub use arc_function::ContapArcFunction;
pub use cont_ana::ContapContAna;
pub use contour::ContapContour;
pub use iw_line::ContapIWLineOfTheIWalking;
pub use iwalking::{ContapTheIWalking, ContapWalkingInteriorPoint, ContapWalkingPathPoint};
pub use line::ContapLine;
pub use path_point::ContapThePathPointOfTheSearch;
pub use point::ContapPoint;
pub use search_inside::ContapTheSearchInside;
pub use sample::{bounds as contap_bounds, project as contap_project, SampleGrid};
pub use segment::ContapTheSegmentOfTheSearch;
pub use surf_function::ContapSurfFunction;
pub use the_search::ContapTheSearch;
pub use types::{ContapIType, ContapTFunction};
