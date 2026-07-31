//! OCCT foundation — pure Rust port of Open CASCADE Technology FoundationClasses.
//!
//! # Architecture
//! ```text
//! kernel/            (error types, handle protocol, container aliases)
//! precision          (tolerance constants)
//! gp/                (geometric primitives — 3D + 2D)
//! ```

pub mod kernel;
pub mod precision;
pub mod gp;
pub mod bnd;
pub mod elib;
pub mod poly;
pub mod toploc;

pub use gp::*;
pub use kernel::*;
pub use bnd::*;
pub use precision::{Precision, ANGULAR, CONFUSION, RESOLUTION, INTERSECTION, APPROXIMATION, INFINITE};
