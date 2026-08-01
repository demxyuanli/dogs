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
pub mod newton;
pub mod gauss;
pub mod bfgs;
pub mod powell;
pub mod kronrod;
pub mod eigen;
pub mod trig;
pub mod stats;
pub mod fft;
pub mod rng;
pub mod gauss_ls;
pub mod house_full;
pub mod multi_int;
pub mod cholesky;
pub mod lu;
pub mod qr_full;
pub mod cubic_spline;
pub mod polyfit;
pub mod levenberg;
pub mod ode;
pub mod scalar;
pub mod polynomial;
pub use status::MathStatus;
pub use vector::MathVector;
pub use intvec::MathIntVector;
pub use matrix::MathMatrix;
pub use svd::SVD;
pub use crout::Crout;
pub use jacobi::Jacobi;
pub use householder::Householder;
pub use newton::{NewtonSolver, NewtonMinimum};
pub use gauss::{gauss_legendre, integrate};
pub use bfgs::BFGS;
pub use powell::Powell;
pub use kronrod::integrate_adaptive;
pub use eigen::{largest_eigen, smallest_eigen, all_eigenvalues, condition_number};
pub use trig::{trig_roots, has_root};
pub use cholesky::{cholesky, cholesky_solve, is_positive_definite};
pub use lu::{lu_decompose, lu_solve, det_from_lu, LU};
pub use qr_full::{gram_schmidt, qr_solve, qr_least_squares};
pub use cubic_spline::{CubicSpline, natural_cubic_spline};
pub use polyfit::{polyfit, polyval};
pub use levenberg::{levenberg_marquardt, LMConfig};
pub use ode::{rk4, rk4_step};
pub use scalar::{bisection, secant, golden_section, brent};
pub use polynomial::{poly_eval, poly_derivative, quadratic_roots, cubic_roots, polynomial_roots};
pub mod matrix_ext;
pub use matrix_ext::{
    identity, transpose, matmul, trace, determinant, inverse, frobenius_norm,
    matrix_rank, det3, det2, solve_linear,
};
pub mod interp;
pub mod roots;
pub mod integrate_adaptive;
pub mod optimize;
pub mod fft2d;
pub mod statistics_full;
pub mod eigen_ext;
pub mod matrix_sparse;
pub mod interp2d;
pub mod distributions;
pub mod spline_surface;
pub mod eig_qr;
pub mod nurbs_fit;
pub mod distrib_extra;
pub mod ode_rk45;
pub mod lsq_nonlinear;
