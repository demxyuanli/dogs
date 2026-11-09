//! `LocOpe` package ports (`ModelingAlgorithms/TKFeat/LocOpe`).
//!
//! Every module mirrors one OCCT class file and keeps the reference's branch
//! structure; the `cxx` / `hxx` line numbers in the comments are the physical
//! lines of the files under `D:/source/OCCT-src`.
//!
//! Module table (`LocOpe` = `src/ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT class | Rust location | Reference | Status |
//! |---|---|---|---|
//! | `LocOpe` | `closed.rs`, `tgt_faces.rs`, `sample_edges.rs` | `LocOpe.cxx:43-236` | `Closed` (both overloads), `TgtFaces`, `SampleEdges`. `IsInside` is commented out in the reference (`cxx:238-265`). |
//! | `LocOpe_WiresOnShape` | `wires_on_shape/` | `hxx:1-120`, `cxx:1-1623`, `lxx:1-40` | `mod.rs` (ctor / `Init` / four `Bind` / `BindAll` / iterators / `OnVertex` / `OnEdge`), `put_pcurve.rs`, `put_pcurves.rs`, `internal_intersections.rs`, `projection.rs`, `add_edges.rs`. |
//! | `LocOpe_BuildWires` | `build_wires.rs` | `hxx:1-50`, `cxx:1-281` | `Perform` (`cxx:63-229`) + file-local `FindFirstEdge` (`cxx:250-279`). |
//! | `LocOpe_Spliter` | `spliter/` | `hxx:1-73`, `cxx:1-719`, `lxx:1-70` | `mod.rs` (ctor / `Init` / accessors), `perform.rs` (stages 1-6, `cxx:63-266`), `finish.rs` (stages 7-9, `cxx:268-577`), `statics.rs` (`Put` / `RebuildWires` / `Select`, `cxx:620-719`). |
//! | `LocOpe_BuildShape` | `build_shape.rs` | `hxx:1-46`, `cxx:1-463` | `Perform` (`cxx:54-320`) + file-local `Add` (`cxx:323-362`) / `Propagate` (`cxx:367-436`) / `IsInside` (`cxx:441-463`). |
//! | `LocOpe_SplitShape` | `split_shape/` | `hxx:1-96`, `cxx:1-1776`, `lxx:1-40` | `mod.rs` (`Init` / `CanSplit` / `DescendantShapes` / `Put` / `Rebuild` / `LeftOf` / `updateToleraces`), `add_vertex.rs`, `add_wires.rs`, `add_closed_wire.rs`, `add_open_wire.rs`, `direction.rs`, `is_inside.rs`. |
//! | `LocOpe_SplitDrafts` | `split_drafts/` | `hxx:1-113`, `cxx:1-1476`, `lxx` none | `mod.rs` (`Init` / `IsDone` / `Shape` / `ShapesFromShape`), `perform.rs` (the ten-parameter `Perform` driver, `cxx:105-1436`), `wire_split.rs` (`cxx:180-351`), `split_face.rs` (`cxx:353-403`), `generated.rs` (`cxx:405-470`), `vertices.rs` (`cxx:471-866`), `wire_edges.rs` (`cxx:868-1100`), `draft_faces.rs` (`cxx:1102-1358`), `descendants.rs` (`cxx:1360-1435`), `statics.rs` (`NewPlane` / `MakeFace` / `Contains` / `NewEdge`, `cxx:1461-1976`). |
//! | `LocOpe_Operation` | `operation.rs` | `hxx:1-27` | The `enum` (`FUSE` / `CUT` / `INVALID`) as [`operation::LocOpeOperation`]. |
//! | `LocOpe_GeneratedShape` | `generated_shape.rs` | `hxx:1-52` | The four query methods as the [`generated_shape::GeneratedShape`] trait; the protected caches are not part of it. |
//! | `LocOpe_GluedShape` | `glued_shape.rs` | `hxx:1-67`, `cxx:1-226` | `Init`, `GlueOnFace`, `MapEdgeAndVertices`, the four queries (`cxx:38-226`). |
//! | `LocOpe_Gluer` | `gluer/` | `hxx:1-82`, `cxx:1-555`, `lxx:1-91` | `mod.rs` (ctors / `Init` / both `Bind` / accessors / `DescendantFaces`), `perform.rs` (`Perform` stages 1-7, `cxx:154-333`), `get_orientation.rs` (`cxx:354-451`), `statics.rs` (`Contains`, plus `AddEdges` as a documented no-op). |
//! | `LocOpe_Generator` | `generator/` | `hxx:1-70`, `cxx:1-1447`, `lxx:1-64` | `mod.rs` (ctor / `Init` / `Perform` / accessors / `DescendantFace`), `statics.rs` (three `ToFuse` + `NewParameter`, `cxx:1232-1447`), `perform/` (`mod.rs` seeds and fusion search, `edge_fusion.rs`, `face_fusion/` (`mod.rs`, `curves.rs`), `assemble.rs`). |
//! | `LocOpe_DPrism` | `dprism/` | `hxx:1-96`, `cxx:1-838` | `mod.rs` (both constructors / accessors), `profiles.rs` (both profile wires), `perform.rs` (the four `Perform` stages), `sections.rs` (`Curves` / `BarycCurve`). |
//! | `LocOpe_PntFace` | `pnt_face.rs` | `hxx:1-79` | The record plus its six accessors. |
//! | `LocOpe_FindEdges` | `find_edges.rs` | `hxx:1-59`, `cxx:1-409`, `lxx:1-67` | `Set` (`cxx:37-409`) and the paired-list iterator (`.lxx:34-67`). |
//! | `LocOpe_FindEdgesInFace` | `find_edges_in_face.rs` | `hxx:1-61`, `cxx:1-193`, `lxx:1-58` | `Set` (`cxx:42-193`) and its iterator (`.lxx:34-58`). |
//! | `LocOpe_CSIntersector` | `cs_intersector/` | `hxx:1-151`, `cxx:1-589`, `lxx:1-46` | `mod.rs` (`Init` / `IsDone` / `NbPoints` / `Point` / four `Localize*`), `perform.rs` (the three `Perform` overloads, `cxx:69-177`), `add_points.rs` (static `AddPoints`, `cxx:517-589`). |
//! | `LocOpe_CurveShapeIntersector` | `curve_shape_intersector.rs` | `hxx:1-134`, `cxx:1-429`, `lxx:1-60` | Both `Init`s, `NbPoints` / `Point`, four `Localize*` and the static `Perform` (`cxx:357-428`). |
//! | `LocOpe_Pipe` | `pipe/` | `hxx:1-70`, `cxx:1-550`, `lxx:1-45` | `mod.rs` (ctor / `Shape` / `Shapes` / the four accessors), `constructor.rs` (the face-fusion pass, `cxx:51-273`), `curves.rs` (`Curves` / `BarycCurve` / static `Orientation`, `cxx:329-548`). |
//!
//! The four `Localize*` scans are shared by both intersectors and live in
//! `localize.rs` (`CSIntersector.cxx:316-513`,
//! `CurveShapeIntersector.cxx:75-354`); the only difference between the two
//! callers is the `Tol` argument.
//!
//! Shared shape-keyed containers live in [`shape_maps`] and follow
//! `TopTools_ShapeMapHasher`; `LocOpe_SplitShape`'s descendant map is a
//! `ShapeValueMap<Vec<TopoShape>>` (`NCollection_DataMap<TopoDS_Shape,
//! NCollection_List<TopoDS_Shape>>`).
//!
//! Not ported from the package: nothing. Every `LocOpe` class file has a
//! module here; `LocOpe.hxx` itself only restates the `TKFeat` enum and the
//! two sampling helpers, which live in `closed.rs` / `tgt_faces.rs` /
//! `sample_edges.rs`.

pub mod build_shape;
pub mod build_wires;
pub mod closed;
pub mod cs_intersector;
pub mod curve_shape_intersector;
pub mod dprism;
pub mod find_edges;
pub mod find_edges_in_face;
pub mod generated_shape;
pub mod generator;
pub mod glued_shape;
pub mod gluer;
pub mod linear_form;
mod localize;
pub mod operation;
pub mod pipe;
pub mod pnt_face;
pub mod prism;
pub mod revol;
pub mod revolution_form;
pub mod sample_edges;
pub mod shape_maps;
pub mod split_drafts;
pub mod split_shape;
pub mod spliter;
pub mod tgt_faces;
pub mod wires_on_shape;

pub use build_shape::BuildShape;
pub use build_wires::BuildWires;
pub use closed::{closed_edge, closed_wire};
pub use cs_intersector::CSIntersector;
pub use curve_shape_intersector::CurveShapeIntersector;
pub use dprism::DPrism;
pub use find_edges::FindEdges;
pub use find_edges_in_face::FindEdgesInFace;
pub use generated_shape::GeneratedShape;
pub use generator::Generator;
pub use glued_shape::GluedShape;
pub use gluer::Gluer;
pub use linear_form::LinearForm;
pub use operation::LocOpeOperation;
pub use pipe::Pipe;
pub use pnt_face::PntFace;
pub use prism::Prism;
pub use revol::Revol;
pub use revolution_form::RevolutionForm;
pub use sample_edges::sample_edges;
pub use shape_maps::{
    ShapeDataMap, ShapeDoubleMap, ShapeIndexedMap, ShapeSet, ShapeValueMap,
};
pub use split_drafts::SplitDrafts;
pub use split_shape::SplitShape;
pub use spliter::Spliter;
pub use tgt_faces::tgt_faces;
pub use wires_on_shape::WiresOnShape;
