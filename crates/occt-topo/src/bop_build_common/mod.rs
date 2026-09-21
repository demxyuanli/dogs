//! BOPAlgo_Builder — container and internal-shape filling (Phase 20).
//!
//! Port of the container/assembly stages of the General Fuse builder
//! (`BOPAlgo_Builder`):
//!
//! | Rust function              | OCCT source                                       |
//! |----------------------------|---------------------------------------------------|
//! | [`fill_images_containers`] | `FillImagesContainers` + `FillImagesContainer`    |
//! |                            | (`BOPAlgo_Builder_1.cxx`)                         |
//! | [`fill_images_compounds`]  | `FillImagesCompounds` + `FillImagesCompound`      |
//! |                            | (`BOPAlgo_Builder_1.cxx`)                         |
//! | [`fill_internal_vertices`] | `FillInternalVertices` (`BOPAlgo_Builder_2.cxx`)  |
//! | [`fill_internal_shapes`]   | `FillInternalShapes` (`BOPAlgo_Builder_3.cxx`)    |
//! | [`build_draft_solid`]      | `BuildDraftSolid` (`BOPAlgo_Builder_3.cxx`)       |
//!
//! The five entry points operate on any host that implements [`BopBuildOps`]
//! (the minimal contract the builder main class — `BOPAlgo_Builder` /
//! `crate::bop_builder2::BopBuilder` — satisfies): the data structure, the
//! images history, the arguments and the origins back-map. Porting through the
//! trait keeps this module independent of the concrete `BopBuilder` fields, so
//! it compiles and is tested standalone.
//!
//! ## Semantics preserved from OCCT
//!
//! * **Containers** — a wire/shell is rebuilt only when at least one of its
//!   direct sub-shapes carries a non-trivial image. A shell is reassembled with
//!   [`crate::shell_splitter::ShellSplitter`] so the image shells are closed;
//!   a wire is reassembled with [`crate::builder::TopoBuilder::make_wire`].
//! * **Compounds** — recursively rebuilt with the splits of their sub-shapes,
//!   keeping each split oriented as the original sub-shape.
//! * **Internal vertices** — alone vertices of a split face are classified
//!   against each face image (2-D, `IntTools_Context::ComputeVF`) and added as
//!   `INTERNAL` children when they fall strictly inside it.
//! * **Internal shapes** — vertices/edges/wires from the arguments and from
//!   inside the source solids are classified against each split solid (3-D, one
//!   representative point) and added as `INTERNAL` children; settling a shape
//!   into an *original* (un-split) solid copies it first, preserving the input.
//! * **Draft solid** — a shell is rebuilt from the face splits (reversing a
//!   split face whose orientation is inverted relative to its original), flagged
//!   closed and wrapped into a solid.
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::GpPnt;

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::algo_tools::AlgoTools;
pub(crate) use crate::bop_hist::BopHistory;
pub(crate) use crate::bopds::BopdsDS;
pub(crate) use crate::brep_extrema::{closest_point_on_edge, closest_point_on_face, is_inside};
pub(crate) use crate::brep_surface::surface_closest_params;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::fclass2d::FaceState;
pub(crate) use crate::int_tools_full::IntToolsContext;
pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
pub(crate) use crate::shell_splitter::ShellSplitter;
pub(crate) use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of};

}


mod build_ops;
mod draft_solid;
pub use build_ops::*;
pub use draft_solid::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
