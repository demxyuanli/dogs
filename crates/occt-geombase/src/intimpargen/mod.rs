//! `IntImpParGen` package (TKGeomAlgo,
//! `ModelingAlgorithms/TKGeomAlgo/IntImpParGen/`): intersection between an
//! implicit curve and a parametric curve.
//!
//! ## Ported
//!
//! * [`gen`] — `IntImpParGen` (`IntImpParGen.cxx:26-251`): `NormalizeOnDomain`,
//!   `DeterminePosition` and the two `DetermineTransition` overloads.
//! * [`intersector`] — `IntImpParGen_Intersector` (`IntImpParGen_Intersector.gxx`)
//!   monomorphised onto the `IntCurve_IntImpConicParConic` instantiation, plus
//!   `IntCurve_MyImpParToolOfIntImpConicParConic` (`_0.cxx`).
//!
//! ## UNPORTED
//!
//! * `IntImpParGen_Tool.hxx` / `.cxx` (1 597 + 5 379 B). Its
//!   `IntImpParGen_Tool::NormalizeOnDomain` / `Determine_Position` /
//!   `Determine_Transition` (`IntImpParGen_Tool.cxx:24-203`) duplicate the
//!   `IntImpParGen` names with different bodies (a one-condition `while` and
//!   `gp::Resolution()` instead of a two-condition `while` and
//!   `TOLERANCE_ANGULAIRE`). The intersector does not call them, so they are not
//!   merged into [`gen`]; port them separately if a caller needs them.

pub mod gen;
pub mod intersector;

pub use intersector::{IntImpParGenIntersector, MyImpParTool};
