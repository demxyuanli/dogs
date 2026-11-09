//! OCCT foundation layer - pure Rust port of Open CASCADE Technology TKernel.
//!
//! This is the stable base of the whole stack: it has no crate dependencies.
//!
//! # Architecture
//! ```text
//! kernel/            (error types, handle protocol, container aliases)
//! precision          (tolerance constants)
//! message            (message printer / gravity)
//! quantity           (physical quantity + units)
//! ```

pub mod kernel;
pub mod message;
pub mod precision;
pub mod quantity;

pub use kernel::*;
pub use precision::{
    Precision, ANGULAR, APPROXIMATION, CONFUSION, INFINITE, INTERSECTION, RESOLUTION,
};
