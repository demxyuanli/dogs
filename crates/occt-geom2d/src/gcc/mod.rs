//! `GccEnt` + `GccAna` + `Geom2dGcc` subset reached by the `ChFi2d` pipeline.
//!
//! | OCCT package | port |
//! |---|---|
//! | `GccEnt_Position`, `GccEnt_QualifiedCirc`, `GccEnt_QualifiedLin` | [`gcc_ent`] |
//! | `Geom2dGcc_QualifiedCurve`, `Geom2dGcc_QCurve` | [`qualified_curve`] |
//! | `GccAna_Circ2d2TanRad` | [`circ2d_2tan_rad`] |
//! | `Geom2dGcc_Circ2d2TanRad` | [`geom2d_circ2d_2tan_rad`] |
//!
//! Only the tangent-to-two-arguments-with-a-given-radius solver is ported: it
//! is the one `ChFi2d_Builder` calls (`ChFi2d_Builder.cxx:1015-1018`). The
//! other `GccAna_*` / `Geom2dGcc_*` classes of the two packages
//! (`Circ2d3Tan`, `Circ2d2TanOn`, `Lin2d2Tan`, `Circ2dTanCen`, `Lin2dTanObl`,
//! `Circ2dTanOnRad`, ...) are not referenced from any ported code.

pub mod circ2d_2tan_rad;
pub mod gcc_ent;
pub mod geom2d_circ2d_2tan_rad;
pub mod qualified_curve;

pub use circ2d_2tan_rad::GccAnaCirc2d2TanRad;
pub use gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
pub use geom2d_circ2d_2tan_rad::Geom2dGccCirc2d2TanRad;
pub use qualified_curve::{Geom2dGccQualifiedCurve, QualifiedCurve2d};
