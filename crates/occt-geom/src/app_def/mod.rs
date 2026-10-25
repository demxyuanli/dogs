//! `AppDef` — the data passed to the parallel-curves approximation drivers.
//!
//! Source: `ModelingData/TKGeomBase/AppDef` (OCCT 8.0.0).
//!
//! Layout follows the OCCT classes:
//! - [`multi_point_constraint`] — `AppDef_MultiPointConstraint`,
//! - [`multi_line`] — `AppDef_MultiLine`,
//! - [`my_line_tool`] — `AppDef_MyLineTool`.
//!
//! Still UNPORTED from the same package (needed by `GeomFill_AppSurf`):
//! - `AppDef_Compute` (`Approx_ComputeLine.gxx:1-1630`),
//! - `AppDef_BSplineCompute` (`Approx_BSplComputeLine.gxx:1-1291`),
//! - `AppDef_LinearCriteria` (`AppDef_LinearCriteria.cxx:1-719`) and the
//!   least-squares / gradient / BFGS support it drives,
//! - `AppDef_Variational` (`AppDef_Variational.cxx:1-1413`).

pub mod multi_line;
pub mod multi_point_constraint;
pub mod my_line_tool;

pub use multi_line::MultiLine;
pub use multi_point_constraint::MultiPointConstraint;
pub use my_line_tool::MyLineTool;
