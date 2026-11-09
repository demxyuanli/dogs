//! OCCT elementary geometry layer - pure Rust port of Open CASCADE Technology TKMath.
//!
//! Depends only on `occt-foundation`. Nothing here may reference topology.
//!
//! # Architecture
//! ```text
//! gp/                (elementary 3D + 2D primitives)
//! bnd/               (axis-aligned boxes)
//! elib/ cslib/ bspl/ (elementary / conic / B-spline evaluation libraries)
//! poly/ toploc/      (triangulation and location data)
//! gprop/             (global properties)
//! convert/ approx/   (conversion and approximation helpers)
//! gcpnts/ int/ intres2d/ intana2d/ numerics (math_*)
//! ```

pub mod approx;
pub mod bnd;
pub mod bspl;
pub mod convert;
pub mod cslib;
pub mod elib;
pub mod gcpnts;
pub mod gp;
pub mod gprop;
pub mod hull;
pub mod int;
pub mod intana2d;
pub mod intres2d;
pub mod math_bracketed_root;
pub mod math_bullard;
pub mod math_direct_poly_roots;
pub mod math_fn;
pub mod math_function_all_roots;
pub mod math_function_roots;
pub mod math_function_sample;
pub mod math_newton_function_root;
pub mod math_trig_roots;
pub mod numeric;
pub mod poly;
pub mod toploc;

pub use bnd::*;
pub use gp::*;
