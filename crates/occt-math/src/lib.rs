//! OCCT linear algebra — pure Rust port of TKMath/math module.
//!
//! # Architecture
//! ```text
//! status           (solver result enum)
//! vector           (dynamic f64 vector, 1-based indexing)
//! intvec           (dynamic i32 vector, 1-based indexing)
//! matrix           (dynamic f64 matrix, Gauss elimination, inversion)
//! ```
//!
//! Solver algorithms (Gauss, SVD, Crout, Householder, Jacobi) and
//! optimization (BFGS, Brent, Newton, PSO) are deferred to sub-modules.

pub mod status;
pub mod vector;
pub mod intvec;
pub mod matrix;
pub mod svd;
pub mod crout;
pub mod jacobi;
pub mod householder;

pub use status::MathStatus;
pub use vector::MathVector;
pub use intvec::MathIntVector;
pub use matrix::MathMatrix;
pub use svd::SVD;
pub use crout::Crout;
pub use jacobi::Jacobi;
pub use householder::Householder;
