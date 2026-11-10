//! The slice of `BRepAdaptor_Surface` / `BRepAdaptor_Curve` that
//! `TopOpeBRepTool_ShapeTool` reads: `GetType()`, `FirstUParameter()` /
//! `FirstVParameter()` / `FirstParameter()`, and `D1`.
//!
//! Sources:
//! - `TKBRep/BRepAdaptor/BRepAdaptor_Surface.cxx:56-82` (`Initialize`)
//! - `TKBRep/BRepAdaptor/BRepAdaptor_Surface.hxx:28-40` (`Load(S, Trsf)`)
//! - `TKG3d/GeomAdaptor/GeomAdaptor_Surface.cxx:402-431` (`load`: type chain)
//! - `TKG3d/GeomAdaptor/GeomAdaptor_Surface.hxx:150-190`
//! - `TKBRep/BRepAdaptor/BRepAdaptor_Curve.hxx:36-98`
//!
//! The adaptor applies the shape's `TopLoc_Location` to every evaluated point
//! and derivative (`GeomAdaptor_Surface::D1` / `GeomAdaptor_Curve::D1` return
//! `myTrsf`-transformed values), so these structs hold the world-coordinate
//! geometry, matching [`crate::brep_tool::BRepTool::face_surface_world`] and
//! [`crate::brep_tool::BRepTool::edge_curve_adaptor_world`].

mod curve_adaptor;
mod geom_adaptor_curve;
mod surface_adaptor;

pub use curve_adaptor::CurveAdaptor;
pub use geom_adaptor_curve::GeomAdaptorCurve;
pub use surface_adaptor::SurfaceAdaptor;
