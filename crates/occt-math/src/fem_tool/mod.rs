//! `FEmTool` — the polynomial-element finite element toolkit used by
//! `AppDef_Variational` to smooth an approximation.
//!
//! Source: `ModelingData/TKGeomBase/FEmTool` (OCCT 8.0.0).
//!
//! Layout follows the OCCT classes:
//! - [`array2`] — the `NCollection_Array2<double>` shape the criteria and the
//!   curve exchange (`FEmTool_ElementaryCriterion::Set`, `FEmTool_Curve::
//!   SetElement`),
//! - [`eval_length`] — `PLib::EvalLength` (`PLib.cxx:2244-2339`), used by
//!   `FEmTool_Curve::Length`,
//! - [`gauss_set`] — `math_GaussSetIntegration` (`math_GaussSetIntegration.cxx:
//!   27-100`), used by the reference-matrix criteria,
//! - [`ref_matrix`] — `FEmTool_ElementsOfRefMatrix`,
//! - [`criteria`] — `FEmTool_ElementaryCriterion` and its three linear
//!   criteria (tension / flexion / jerk),
//! - [`curve`] — `FEmTool_Curve`,
//! - [`profile_matrix`] — `FEmTool_SparseMatrix` / `FEmTool_ProfileMatrix`,
//! - [`assembly`] — `FEmTool_Assembly`.
//!
//! Translation notes:
//! - The OCCT 1-based `NCollection_Array1` / `Array2` become slices and
//!   [`Array2`]; the loops keep the OCCT bounds so the arithmetic order (and
//!   therefore the floating point result) is unchanged.
//! - `math_Matrix` / `math_Vector` keep their 1-based addressing: the port uses
//!   [`occt_core::kernel::containers`]-free `MathMatrix` / `MathVector` from
//!   this crate, whose `new(..)` takes the OCCT lower/upper bounds.

pub mod array1;
pub mod array2;
pub mod assembly;
pub mod criteria;
pub mod curve;
pub mod eval_length;
pub mod gauss_set;
pub mod profile_matrix;
pub mod ref_matrix;

pub use array1::Array1d;
pub use array2::Array2;
pub use assembly::Assembly;
pub use criteria::{Criterion, LinearFlexion, LinearJerk, LinearTension};
pub use curve::Curve;
pub use gauss_set::gauss_set_integration;
pub use profile_matrix::ProfileMatrix;
pub use ref_matrix::ElementsOfRefMatrix;
