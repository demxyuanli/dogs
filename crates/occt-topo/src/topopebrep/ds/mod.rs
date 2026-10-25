//! `TopOpeBRepDS` — the data structures of the classic boolean engine.
//!
//! Source package: `TKBool/TopOpeBRepDS` (OCCT 8.0.0, 18,382 lines). The
//! package stores the new geometries (points, curves, surfaces), the
//! topological shapes, and the interferences that attach one to the other.
//!
//! Ported so far:
//! - `TopOpeBRepDS_Kind` / `TopOpeBRepDS_Config` and the static helpers of
//!   `TopOpeBRepDS.cxx` ([`kind`]).
//! - `TopOpeBRepDS_Transition` ([`transition`]).
//! - `TopOpeBRepDS_GeometryData` ([`geometry_data`]).
//! - `TopOpeBRepDS_Interference` ([`interference`]).
//! - `TopOpeBRepDS_ShapeData` ([`shape_data`]).
//! - `TopOpeBRepDS_Point` / `TopOpeBRepDS_PointData` ([`point`],
//!   [`point_data`]).
//! - `TopOpeBRepDS_Surface` / `TopOpeBRepDS_SurfaceData` ([`surface`],
//!   [`surface_data`]).
//! - `TopOpeBRepDS_Curve` / `TopOpeBRepDS_CurveData` ([`curve`],
//!   [`curve_data`]).
//! - `TopOpeBRepDS_ShapeWithState` / `TopOpeBRepDS_ListOfShapeOn1State`
//!   ([`shape_with_state`], [`list_of_shape_on1_state`]).
//! - `TopOpeBRepDS_DataStructure` ([`data_structure`]).
//! - `TopOpeBRepDS_TOOL` ([`tool`]).
//! - `TopOpeBRepDS_connex` ([`connex`]).
//! - `TopOpeBRepDS_GapTool` / `TopOpeBRepDS_GapFiller` ([`gap_tool`],
//!   [`gap_filler`]).
//! - `TopOpeBRepDS_Check` — the integrity checker ([`check`]).
//! - `TopOpeBRepDS_Explorer` — the type-filtered shape walk ([`explorer`]).
//! - `TopOpeBRepDS_Dumper` — printable shape references ([`dumper`]).
//! - `TopOpeBRepDS_EXPORT`'s `FUN_ds_*` completion family ([`complete_ds`]) and
//!   `TopOpeBRepDS_redu` ([`redu`]).
//!
//! `NCollection_List<occ::handle<TopOpeBRepDS_Interference>>` is a
//! `Vec<InterferenceHandle>`; OCCT's `Handle` is a shared, mutable reference,
//! so a row is [`InterferenceHandle`] = `Arc<RwLock<Interference>>` and
//! `ChangeTransition()` is `lock().write()`.

pub mod association;
pub mod build_tool;
pub mod check;
pub mod complete_ds;
pub mod connex;
pub mod curve;
pub mod curve_data;
pub mod curve_explorer;
pub mod curve_iterator;
pub mod data_structure;
pub mod dumper;
pub mod edge3d_interference_tool;
pub mod edge_interference_tool;
pub mod eir;
pub mod explorer;
pub mod export;
pub mod face_interference_tool;
pub mod filter;
pub mod fir;
pub mod gap_filler;
pub mod gap_tool;
pub mod geometry_data;
pub mod h_data_structure;
pub mod interference;
pub mod interference_iterator;
pub mod interference_tool;
pub mod kind;
pub mod list_of_shape_on1_state;
pub mod marker;
pub mod point;
pub mod point_data;
pub mod point_explorer;
pub mod point_iterator;
pub mod process_interferences_tool;
pub mod redu;
pub mod surface_iterator;
pub mod tki;
pub mod tool;
pub mod reducer;
pub mod repvg;
pub mod shape_data;
pub mod shape_with_state;
pub mod surface;
pub mod surface_data;
pub mod surface_explorer;
pub mod transition;

pub use association::Association;
pub use build_tool::BuildTool;
pub use check::{check_edge_parameter, Check, CheckStatus};
pub use connex::Connexity;
pub use curve::Curve;
pub use curve_data::CurveData;
pub use curve_explorer::CurveExplorer;
pub use curve_iterator::CurveIterator;
pub use data_structure::{DataStructure, IndexedShapeData, IndexedShapeSet, ShapeWithStateMap};
pub use dumper::Dumper;
pub use edge3d_interference_tool::Edge3dInterferenceTool;
pub use edge_interference_tool::EdgeInterferenceTool;
pub use eir::Eir;
pub use explorer::Explorer;
pub use face_interference_tool::{
    fun_compute_geom_data, fun_compute_geom_data_full, fun_parameters, FaceInterferenceTool,
};
pub use gap_filler::GapFiller;
pub use gap_tool::GapTool;
pub use filter::Filter;
pub use fir::{
    fun_find_pon_f, fun_is_pon_f, fun_reduce_edge_geometry, fun_reduce_edge_geometry1,
    fun_state_edge_face, Fir,
};
pub use geometry_data::GeometryData;
pub use h_data_structure::HDataStructure;
pub use interference::{Interference, InterferenceHandle};
pub use interference_iterator::{DerivedMatch, InterferenceIterator, InterferenceList};
pub use kind::{
    config_s_print, is_geometry, is_topology, kind_s_print, kind_to_shape, orientation_s_print,
    shape_to_kind, state_s_print, Config, Kind,
};
pub use list_of_shape_on1_state::ListOfShapeOn1State;
pub use marker::Marker;
pub use point::Point;
pub use repvg::{fds_repvg, fds_repvg2};
pub use reducer::Reducer;
pub use point_data::PointData;
pub use point_explorer::PointExplorer;
pub use point_iterator::PointIterator;
pub use shape_data::ShapeData;
pub use shape_with_state::ShapeWithState;
pub use surface::Surface;
pub use surface_data::SurfaceData;
pub use surface_explorer::SurfaceExplorer;
pub use surface_iterator::SurfaceIterator;
pub use tki::TKI;
pub use tool::DShapeOn1StateMap;
pub use transition::Transition;

