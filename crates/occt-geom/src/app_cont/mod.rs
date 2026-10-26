//! `AppCont` - the continuous-function interface family used by
//! `Approx_ComputeCLine` and the `ProjLib` approximators.
//!
//! Source: `ModelingData/TKGeomBase/AppCont` (`AppCont_Function.hxx`,
//! `AppCont_ContMatrices.pxx`, `AppCont_LeastSquare.cxx`, OCCT 8.0.0).
//!
//! - [`function`] - `AppCont_Function` (`AppCont_Function.hxx`).
//! - [`cont_matrices`] - `AppCont_ContMatrices` (`AppCont_ContMatrices.pxx`).
//! - [`least_square`] - `AppCont_LeastSquare` (`AppCont_LeastSquare.cxx`).

pub mod cont_matrices;
pub mod function;
pub mod least_square;

pub use cont_matrices::{ibt_matrix, ibp_matrix, inv_m_matrix, m_matrix, v_bernstein};
pub use function::Function;
pub use least_square::{LeastSquare, PeriodicityInfo};
