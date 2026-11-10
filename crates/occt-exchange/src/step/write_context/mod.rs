//! `StepToGeom`/`GeomToStep`-style STEP writer context (single model / assembly).
//!
//! Split by concern; every item is re-exported so `crate::step::*` keeps the
//! same surface as before.

mod api;
mod assembly;
mod bspline_convert;
mod context;
mod options;
mod options_write;
mod surface_params;

pub(crate) use context::WriteCtx;
pub use api::*;
pub use assembly::*;
pub use options::*;
pub use options_write::*;
