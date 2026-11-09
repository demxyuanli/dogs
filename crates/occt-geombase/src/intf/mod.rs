//! Port of the `Intf` package
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/FILES.cmake`): the polygon /
//! polyhedron interference engine behind `Geom2dInt_*`,
//! `IntCurve_IntPolyPolyGen`, `IntCurveSurface_TheInterferenceOfHInter`
//! and `HLRBRep_TheInterferenceOfInterCSurf`.
//!
//! Ported modules:
//! * `interference` - `Intf_Interference`
//!   (`Intf_Interference.hxx/.cxx/.lxx`).
//! * `interference_polygon2d` - `Intf_InterferencePolygon2d`
//!   (`Intf_InterferencePolygon2d.hxx/.cxx`).
//! * `interference_polygon_polyhedron` -
//!   `Intf_InterferencePolygonPolyhedron`
//!   (`Intf_InterferencePolygonPolyhedron.gxx`), as two traits
//!   (`IntfPolygon3dTool` / `IntfPolyhedronTool`) plus a generic
//!   engine.
//! * `plane` - the `Intf::` static helpers (`Intf.hxx/.cxx`).
//! * `tool` - `Intf_Tool` (`Intf_Tool.hxx/.cxx`).
//! * `polygon2d` / `section_line` / `section_point` /
//!   `tangent_zone` - the remaining 2D support classes.
//!
//! UNPORTED in this package:
//! * `Intf_Tool`'s 3D conic paths - `HyprBox`
//!   (`Intf_Tool.cxx:967-1119`), `Inters3d(gp_Hypr)`
//!   (`Intf_Tool.cxx:1123-1302`), `Inters3d(gp_Parab)`
//!   (`Intf_Tool.cxx:1306-1485`), `ParabBox`
//!   (`Intf_Tool.cxx:1489-1630`). They need `IntAna_IntConicQuad`
//!   and the 3D `ElCLib::D1(gp_Hypr/gp_Parab)` overloads, neither of which
//!   is ported; see `tool.rs`.
//! * the `Extrema_ExtElC` proximity tail of both
//!   `Intf_InterferencePolygonPolyhedron::Intersect` overloads
//!   (`Intf_InterferencePolygonPolyhedron.gxx:981-1052` and
//!   `:1296-1367`); `Extrema_ExtElC` is not ported.
//! * the `#if 0` `Intersect` overload
//!   (`Intf_InterferencePolygonPolyhedron.gxx:630-739`) is dead code in
//!   OCCT and is not compiled there either.

pub mod interference;
pub mod interference_polygon2d;
pub mod interference_polygon_polyhedron;
pub mod plane;
pub mod polygon2d;
pub mod section_line;
pub mod section_point;
pub mod tangent_zone;
pub mod tool;

pub use interference::IntfInterference;
pub use interference_polygon2d::IntfInterferencePolygon2d;
pub use interference_polygon_polyhedron::{
    IntfInterferencePolygonPolyhedron, IntfPolygon3dTool, IntfPolyhedronTool, IntfPolyhGrid,
};
pub use plane::{contain, plane_equation};
pub use polygon2d::IntfPolygon2d;
pub use section_line::IntfSectionLine;
pub use section_point::{IntfPIType, IntfSectionPoint};
pub use tangent_zone::IntfTangentZone;
pub use tool::IntfTool;
