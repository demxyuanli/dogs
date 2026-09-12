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
pub mod geom;
pub mod numeric;
pub mod hull;
pub mod int;
pub mod io;
pub mod bspl;
pub mod bvh;
pub mod convert;
pub mod cslib;
pub mod gprop;
pub mod quantity;
pub mod message;

pub use gp::*;
pub use kernel::*;
pub use bnd::*;
pub use precision::{Precision, ANGULAR, CONFUSION, RESOLUTION, INTERSECTION, APPROXIMATION, INFINITE};
pub mod gcpnts;
pub mod math_fn;
pub mod math_bracketed_root;
pub mod math_function_roots;
