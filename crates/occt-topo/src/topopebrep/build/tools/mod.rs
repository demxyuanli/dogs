//! `TopOpeBRepBuild_Tools` — auxiliary methods used by
//! `TopOpeBRepBuild_Builder1`.
//!
//! Source: `TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Tools.hxx:42-173`,
//! `TopOpeBRepBuild_Tools.cxx:60-924`, `TopOpeBRepBuild_Tools_1.cxx:18-558`
//! (OCCT 8.0.0).
//!
//! Split by subject:
//! - [`state`] — the state propagation family (`FindState*`,
//!   `PropagateState*`, `SpreadStateToChild`, `FindStateThroughVertex`),
//!   `Tools.cxx:60-457`.
//! - [`geom`] — `GetAdjacentFace`, the normal/tangent queries and
//!   `IsDegEdgesTheSame`, `Tools.cxx:461-662`, `:838-853`.
//! - [`pcurve`] — `UpdatePCurves` / `UpdateEdgeOnFace` /
//!   `UpdateEdgeOnPeriodicalFace` / `NormalizeFace` / `CorrectFace2d`,
//!   `Tools.cxx:664-924`.
//!
//! UNPORTED, with their OCCT location and why:
//! - `CorrectTolerances` / `CorrectCurveOnSurface` / `CorrectPointOnCurve`
//!   (`Tools_1.cxx:57-322`): no caller anywhere in OCCT 8.0.0 (dead code),
//!   and the body needs the `BRep_TEdge` CurveRepresentation list that this
//!   port keys by surface instead.
//! - `CheckFaceClosed2d` (`Tools_1.cxx:545-558`): needs
//!   `BRepCheck_Wire::Orientation(F)` (`BRepCheck_Wire.cxx:786-1050`), which
//!   the port does not carry (`brep_check.rs` only has `check_wire`).

mod geom;
mod pcurve;
mod state;

pub use geom::*;
pub use pcurve::*;
pub use state::*;
