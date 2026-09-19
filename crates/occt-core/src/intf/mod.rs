//! Port of the 2D subset of the `Intf` package
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/FILES.cmake`): the polygon
//! interference engine behind
//! `Geom2dInt_TheIntPCurvePCurveOfGInter` and `IntCurve_IntPolyPolyGen`.
//!
//! UNPORTED in this package:
//! * `Intf_Tool` (`Intf_Tool.cxx`, 1653 lines): its consumers are 3D only -
//!   `IntCurveSurface_Inter.pxx` and `IntCurveSurface` (`Intf_Tool::LinBox`,
//!   `:LinTetra`, `:LinSphere`, ... used at `IntCurvesFace_Intersector.cxx:392`
//!   and `Intf_InterferencePolygonPolyhedron.gxx`). No 2D path calls it; the
//!   port already carries the `LinBox` equivalent inline
//!   (`int_curves_face.rs:130`).
//! * `Intf_InterferencePolygonPolyhedron.gxx` and `Intf.cxx` (`Intf::` helpers):
//!   3D polyhedron interference, not reachable from `Geom2dInt_*`.

pub mod interference;
pub mod interference_polygon2d;
pub mod polygon2d;
pub mod section_line;
pub mod section_point;
pub mod tangent_zone;

pub use interference::IntfInterference;
pub use interference_polygon2d::IntfInterferencePolygon2d;
pub use polygon2d::IntfPolygon2d;
pub use section_line::IntfSectionLine;
pub use section_point::{IntfPIType, IntfSectionPoint};
pub use tangent_zone::IntfTangentZone;
