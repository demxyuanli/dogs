//! Hidden Line Removal algorithm data (`ModelingAlgorithms/TKHLR/HLRAlgo`),
//! OCCT 8.0.0.
//!
//! Ported:
//! * `Intrv_Interval` / `Intrv_Position` in [`intrv_interval`]
//!   (`Intrv_Interval.cxx:38-151`, `.lxx`, `.hxx`).
//! * `Intrv_Intervals` in [`intrv_intervals`]
//!   (`Intrv_Intervals.cxx:38-248`): the sorted non-overlapping interval
//!   sequence with `Intersect` / `Subtract` / `Unite` / `XUnite`.
//! * `HLRAlgo_Projector` in [`projector`]
//!   (`HLRAlgo_Projector.cxx:38-396`): axonometric / perspective projector
//!   with the four hard-coded standard directions of `TrsfType`
//!   (`cxx:113-150`), `Project` with and without the derivative, and `Shoot`.
//! * `HLRAlgo_EdgeStatus` in [`edge_status`] (`HLRAlgo_EdgeStatus.cxx:22-123`).
//! * `HLRAlgo_Intersection` in [`intersection`]
//!   (`HLRAlgo_Intersection.cxx:22-50`, `.lxx`).
//! * `HLRAlgo_Coincidence` in [`coincidence`]
//!   (`HLRAlgo_Coincidence.hxx:38-81`, header-only).
//! * `HLRAlgo_Interference` in [`interference`]
//!   (`HLRAlgo_Interference.hxx:30-160`, `HLRAlgo_Interference_0.cxx:24-39`).
//! * `HLRAlgo_EdgesBlock` + `MinMaxIndices` in [`edges_block`]
//!   (`HLRAlgo_EdgesBlock.hxx:27-162`, `.cxx:27-31`).
//! * `HLRAlgo_WiresBlock` in [`wires_block`] (`HLRAlgo_WiresBlock.hxx:27-66`).
//! * `HLRAlgo_EdgeIterator` in [`edge_iterator`]
//!   (`HLRAlgo_EdgeIterator.hxx:26-66`, `.cxx:26-95`, `.lxx:22-71`).
//! * `HLRAlgo_BiPoint` in [`bi_point`] (`HLRAlgo_BiPoint.hxx:29-260`,
//!   `.cxx:24-240`).
//! * `HLRAlgo` static pruning helpers in [`minmax`] (`HLRAlgo.hxx:44-88`,
//!   `HLRAlgo.cxx:19-205`).
//! * `HLRAlgo_PolyMask` / `HLRAlgo_TriangleData` / `HLRAlgo_PolyHidingData`
//!   in [`poly_records`] (header-only).
//! * `HLRAlgo_PolyData` in [`poly_data`] (`HLRAlgo_PolyData.hxx:35-110`,
//!   `.lxx:26-66`, `.cxx:30-1127`): the hiding-triangle test, including the
//!   three-times-duplicated `hideByOneTriangle` edge block.
//! * `HLRAlgo_PolyShellData` in [`poly_shell_data`]
//!   (`HLRAlgo_PolyShellData.hxx:26-63`, `.cxx:26-62`).
//! * `HLRAlgo_PolyAlgo` in [`poly_algo`] (`HLRAlgo_PolyAlgo.hxx:27-94`,
//!   `.cxx:29-597`): the integer bounding-box encoding, the hiding-triangle
//!   plane fill, and the hide / show segment walk.
//! * `HLRAlgo_PolyInternalNode` / `HLRAlgo_PolyInternalSegment` in
//!   [`poly_internal_records`] (header-only).
//! * `HLRAlgo_PolyInternalData` in [`poly_internal_data`]
//!   (`HLRAlgo_PolyInternalData.hxx:30-116`, `.lxx:28-117`, `.cxx:32-1039`):
//!   `UpdateLinks`, the edge-splitting `UpdateLinks(ip1, ip2, ip3)`, `AddNode`
//!   and the three `Inc*` / `Dec*` pairs.
//!
//! UNPORTED: the `HLRBRep_*` algorithms that consume these types
//! (`HLRBRep_PolyAlgo`, `HLRBRep_PolyHLRToShape`, `HLRBRep_Data`, ...); this
//! module is only the shared data of `TKHLR`. `HLRAlgo_PolyInternalData::Dump`
//! (`cxx:877-919`) is an `OCCT_DEBUG`-only `std::cout` transcript and is left
//! out.
//!
//! Failure channel: the raising `gp` constructors and `gp_Trsf` operations
//! (`gp_Dir` of `Shoot`, `SetScaleFactor` / `Invert` of `Scaled`, the
//! `Standard_NoSuchObject` of `Focus`) are `Err` here.

mod bi_point;
mod coincidence;
mod edge_iterator;
mod edge_status;
mod edges_block;
mod interference;
mod intersection;
mod intrv_interval;
mod intrv_intervals;
mod minmax;
mod poly_algo;
mod poly_data;
mod poly_internal_data;
mod poly_internal_records;
mod poly_records;
mod poly_shell_data;
mod projector;
mod wires_block;

pub use bi_point::*;
pub use coincidence::*;
pub use edge_iterator::*;
pub use edge_status::*;
pub use edges_block::*;
pub use interference::*;
pub use intersection::*;
pub use intrv_interval::*;
pub use intrv_intervals::*;
pub use minmax::*;
pub use poly_algo::*;
pub use poly_data::*;
pub use poly_internal_data::*;
pub use poly_internal_records::*;
pub use poly_records::*;
pub use poly_shell_data::*;
pub use projector::*;
pub use wires_block::*;
