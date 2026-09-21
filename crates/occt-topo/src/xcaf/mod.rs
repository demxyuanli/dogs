//! Phase 5 module: xcaf — STEP assembly metadata (`XCAFDoc_ShapeTool`-lite).
//!
//! **UNPORTED (audit A14)**: a minimal, port-local metadata container; OCCT's
//! XCAF is `XCAFDoc_ShapeTool` + `TDocStd_Document` (with labels, colours and
//! layers) and is not translated here.
//!
//! A minimal stand-in for OCCT's XCAF document: tracks the product list of an
//! assembly and per-product name / color / layer attributes, serializes them
//! alongside a STEP physical file, and reads them back. Because a fully valid
//! XCAF (STEP 214 `APPLICATION_PROTOCOL`) layer is out of scope, the metadata
//! is embedded as a comment block before the `DATA` section using a clear
//! marker, which keeps the STEP file valid and round-trippable.
mod prelude {

pub(crate) use std::collections::{BTreeMap, HashMap, HashSet};
pub(crate) use std::f64::consts::PI;

pub(crate) use occt_core::gp::{GpDir, GpMat, GpPnt, GpTrsf, GpXyz, TrsfForm};

pub(crate) use crate::bincaf::{BinXcaf, BinXcafEntry, XcafAttribute};
pub(crate) use crate::model::BRepModel;
pub(crate) use crate::shape::TopoShape;
pub(crate) use crate::xmlcaf::{XmlEntry, XmlXcafDoc};

}


mod assembly;
mod convert;
mod step_io;
pub use assembly::*;
pub use convert::*;
pub use step_io::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;

#[cfg(test)]
#[path = "xcaf_doc_tests.rs"]
mod xcaf_doc_tests;
