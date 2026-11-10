//! `BOPAlgo_Builder` state, options and pipeline entry points.
//!
//! Source: `src/ModelingAlgorithms/TKBO/BOPAlgo/BOPAlgo_Builder*.cxx`, split by
//! concern so every file stays inside the project's source-size budget. The
//! inherent `impl BopBuilder` blocks are additive; the public surface is
//! unchanged.

mod accessors;
mod build;
mod constructor;
mod fill_images;
mod groups;
mod history;
mod perform;
mod report;
