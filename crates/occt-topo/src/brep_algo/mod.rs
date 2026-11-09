//! `BRepAlgo` (`ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo.cxx`, OCCT 8.0.0).
//!
//! Ported:
//! * `ConvertWire` (`cxx:80-219`) and `ConvertFace` (`cxx:223-242`) in
//!   [`convert_wire`].
//! * `ConcatenateWire` (`cxx:245-400`) in [`concatenate_wire`]: the G1
//!   (`GeomConvert::ConcatG1`) or C1 (`GeomConvert::ConcatC1`) merge of a
//!   wire's edges into as few B-splines as the requested continuity allows.
//! * `ConcatenateWireC0` (`cxx:404-802`) in [`concatenate_wire_c0`]: the same
//!   merge with a fixed `ShapeFix_Shape` pass in front (`cxx:411-417`).
//!
//! UNPORTED: nothing in this module. The `ShapeFix_Shape::Perform` call of
//! `ConcatenateWireC0` reaches only the free-wire arm of that class; its
//! solid / shell / face arms stay unported there (see
//! `shhealing::shape_fix_shape`), and its precision / max-tolerance pair is
//! pinned to `Precision::Confusion()` by this caller (`cxx:406-407`,
//! `:413-415`).
//!
//! Failure channel: `Standard_ConstructionError` (both `ConcatG1` / `ConcatC1`
//! raise it on a failed `Add`) and a `BRepLib_MakeEdge` that is not done (OCCT
//! adds the resulting null edge to the wire) are the `Err` here, because a null
//! edge cannot be stored.
//!
//! `BRepAlgo_AsDes` (`BRepAlgo_AsDes.cxx`, class `AsDes`) lives in [`as_des`].

mod as_des;
mod concatenate_wire;
mod concatenate_wire_c0;
mod convert_wire;

pub use as_des::AsDes;
pub use concatenate_wire::*;
pub use concatenate_wire_c0::*;
pub use convert_wire::*;
