//! `TopOpeBRep` — the classic boolean engine (`TKBool`).
//!
//! Source packages and sizes (OCCT 8.0.0):
//! `TopOpeBRepDS` (18,382 lines), `TopOpeBRep` (17,507),
//! `TopOpeBRepBuild` (33,788), `TopOpeBRepTool` (17,417).
//!
//! Port order follows the package dependency
//! `TopOpeBRepDS` -> `TopOpeBRep` -> `TopOpeBRepBuild` -> `TopOpeBRepTool`
//! (`specs/occt_migration_guide.md` §5.2, S2).
//!
//! Landed so far:
//! - [`ds`] — the `TopOpeBRepDS` data structures.
//! - [`build`] — the `TopOpeBRepBuild` container and case-table foundation.
//! - [`tool`] — the `TopOpeBRepTool_ShapeTool` orientation family and the
//!   `BRepAdaptor_*` / `LProp_CurveUtils` slice it reads.
//! - [`interf`] — the `TopOpeBRep` point and context classes.

pub mod build;
pub mod ds;
pub mod ds_filler;
pub mod edges_filler;
pub mod edges_intersector;
pub mod face_edge_filler;
pub mod face_edge_intersector;
pub mod faces_filler;
pub mod faces_intersector;
pub mod ff_transition_tool;
pub mod geom_tool;
pub mod interf;
pub mod kpart;
pub mod line_inter;
pub mod shape_scanner;
pub mod shape_intersector;
pub mod shape_intersector2d;
pub mod tool;
pub mod type_line_curve;
pub mod v_point_inter_iterator;
pub mod w_point_inter_iterator;

pub use ds_filler::DsFiller;
pub use edges_filler::EdgesFiller;
pub use edges_intersector::EdgesIntersector;
pub use face_edge_filler::FaceEdgeFiller;
pub use face_edge_intersector::FaceEdgeIntersector;
pub use faces_intersector::FacesIntersector;
pub use ff_transition_tool::FFTransitionTool;
pub use geom_tool::GeomTool;
pub use line_inter::LineInter;
pub use shape_scanner::ShapeScanner;
pub use shape_intersector::ShapeIntersector;
pub use shape_intersector2d::ShapeIntersector2d;
pub use type_line_curve::TypeLineCurve;
pub use v_point_inter_iterator::VPointInterIterator;
pub use w_point_inter_iterator::WPointInterIterator;
