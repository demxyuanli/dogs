//! `TopOpeBRepBuild` — the builders of the classic boolean engine.
//!
//! Source package: `TKBool/TopOpeBRepBuild` (OCCT 8.0.0, 20,209 `.cxx` lines
//! over 43 classes). The package drafts the result shapes from the
//! `TopOpeBRepDS` data structure: the area builders split the faces/edges into
//! `IN`/`ON`/`OUT` areas, the `Grid*` classes walk the connexity grids, and
//! `Builder` drives the whole fuse/cut/common.
//!
//! Ported so far — the container and case-table foundation:
//! - `TopOpeBRepBuild_ShapeListOfShape` ([`shape_list_of_shape`]).
//! - `TopOpeBRepBuild_BlockIterator` ([`block_iterator`]).
//! - `TopOpeBRepBuild_Loop` / `TopOpeBRepBuild_Pave` ([`loop_shape`],
//!   [`pave`]) and `TopOpeBRepBuild_LoopClassifier` ([`loop_classifier`]).
//! - `TopOpeBRepBuild_LoopSet` / `TopOpeBRepBuild_PaveSet` ([`loop_set`],
//!   [`pave_set`]).
//! - `TopOpeBRepBuild_PaveClassifier` ([`pave_classifier`]).
//! - `TopOpeBRepBuild_ShapeSet` ([`shape_set`]) and
//!   `TopOpeBRepBuild_BlockBuilder` ([`block_builder`]).
//! - `TopOpeBRepBuild_GTopo` / `TopOpeBRepBuild_GIter` / `TopOpeBRepBuild_GTool`
//!   ([`g_topo`], [`g_iter`], [`g_tool`]).
//! - `TopOpeBRepBuild_AreaBuilder` and its `Area1dBuilder` / `Area2dBuilder` /
//!   `Area3dBuilder` specialisations ([`area_builder`], [`area1d_builder`],
//!   [`area2d_builder`], [`area3d_builder`]) with `TopOpeBRepBuild_LoopEnum`
//!   ([`loop_enum`]).
//! - `TopOpeBRepBuild_Builder`'s entry points and drivers ([`builder`]:
//!   `Perform`, `MergeShapes` / `MergeSolids` / `MergeSolid`, `Section`,
//!   `SectionCurves`, `SectionEdges`, `FillSecEdgeAncestorMap`).
//! - `TopOpeBRepBuild_HBuilder` ([`h_builder`]), the package's public API.

pub mod area1d_builder;
pub mod area2d_builder;
pub mod area3d_builder;
pub mod area_builder;
pub mod block_builder;
pub mod builder;
pub mod block_iterator;
pub mod composite_classifier;
pub mod correct_face_2d;
pub mod edge_builder;
pub mod face_area_builder;
pub mod face_builder;
pub mod g_iter;
pub mod g_tool;
pub mod g_topo;
pub mod fuse_face;
pub mod h_builder;
pub mod loop_classifier;
pub mod loop_enum;
pub mod loop_set;
pub mod loop_shape;
pub mod pave;
pub mod pave_classifier;
pub mod pave_set;
pub mod shape_list_of_shape;
pub mod shape_set;
pub mod shell_face_classifier;
pub mod shell_face_set;
pub mod shell_to_solid;
pub mod solid_area_builder;
pub mod solid_builder;
pub mod tools;
pub mod tools2d;
pub mod vertex_info;
pub mod wire_edge_classifier;
pub mod wire_edge_set;
pub mod wire_to_face;

pub use area1d_builder::Area1dBuilder;
pub use area2d_builder::Area2dBuilder;
pub use area3d_builder::Area3dBuilder;
pub use area_builder::AreaBuilder;
pub use block_builder::BlockBuilder;
pub use builder::Builder;
pub use builder::Builder1;
pub use block_iterator::BlockIterator;
pub use composite_classifier::CompositeClassifier;
pub use correct_face_2d::CorrectFace2d;
pub use edge_builder::EdgeBuilder;
pub use face_area_builder::FaceAreaBuilder;
pub use face_builder::FaceBuilder;
pub use fuse_face::FuseFace;
pub use g_iter::GIter;
pub use g_tool::GTool;
pub use g_topo::GTopo;
pub use h_builder::HBuilder;
pub use loop_classifier::LoopClassifier;
pub use loop_enum::LoopEnum;
pub use loop_set::{LoopSet, LoopSetSource};
pub use loop_shape::{Loop, LoopHandle, LoopKind};
pub use pave::PaveData;
pub use pave_classifier::PaveClassifier;
pub use pave_set::PaveSet;
pub use shape_list_of_shape::ShapeListOfShape;
pub use shape_set::ShapeSet;
pub use shape_set::ShapeSetSource;
pub use shell_face_classifier::ShellFaceClassifier;
pub use shell_face_set::ShellFaceSet;
pub use shell_to_solid::ShellToSolid;
pub use solid_area_builder::SolidAreaBuilder;
pub use solid_builder::SolidBuilder;
pub use vertex_info::VertexInfo;
pub use wire_edge_classifier::WireEdgeClassifier;
pub use wire_edge_set::WireEdgeSet;
pub use wire_to_face::WireToFace;
