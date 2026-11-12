//! OCCT boundary representation data layer - the pure Rust port of TKBRep.
//!
//! Nothing here models an algorithm. This crate owns the shape representation
//! (`TopoDS_Shape` / `TopoDS_TShape`), its location and orientation algebra
//! (`TopLoc`, `TopoDS_Iterator`, `TopoDS_Transform`), the explorers
//! (`TopExp`), and the low-level builder/tool pair (`BRep_Builder` /
//! `BRep_Tool`) that every modelling algorithm in `occt-topo` writes through.
//!
//! # Architecture
//! ```text
//! abs, shape, tshape       (TopAbs / TopoDS_Shape / TopoDS_TShape)
//! iterator, transform      (TopoDS_Iterator / TopoDS_Shape transforms)
//! topexp, tools, validate  (TopExp, TopoDS helpers, BRepCheck-lite)
//! builder, brep_tool       (BRep_Builder / BRep_Tool)
//! brep_point_representation, brep_builder_update_vertex, brep_lprop,
//! brep_tools_substitution  (BRep_PointRepresentation and its neighbours)
//! tgeometry, topo_tools    (port-local geometry registry and shape helpers)
//! debug_flags              (port-local tracing switches)
//! ```

pub mod abs;
pub mod adaptor;
pub mod brep_bnd_lib;
pub mod brep_builder_transfert;
pub mod brep_builder_update_vertex;
pub mod brep_tool_curve_on_surface;
pub mod brep_lprop;
pub mod brep_point_representation;
pub mod brep_surface;
pub mod brep_tool;
pub mod brep_tools;
pub mod brep_tools_modifier;
pub mod brep_tools_uv;
pub mod brep_tools_wire_explorer;
pub mod brep_uv_bounds;
pub mod brep_tools_substitution;
pub mod builder;
pub mod debug_flags;
pub mod iterator;
pub mod shape;
pub mod shape_empty_copied;
pub mod tgeometry;
pub mod tools;
pub mod topexp;
pub mod topo_tools;
pub mod transform;
pub mod tshape;
pub mod tshape_continuity;
pub mod validate;
