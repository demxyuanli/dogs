//! `GccEnt` + `GccAna` + `Geom2dGcc` tangent solvers.
//!
//! | OCCT package | port |
//! |---|---|
//! | `GccEnt_Position`, `GccEnt_QualifiedCirc`, `GccEnt_QualifiedLin` | [`gcc_ent`] |
//! | `Geom2dGcc_QualifiedCurve`, `Geom2dGcc_QCurve` | [`qualified_curve`] |
//! | `GccAna_Circ2d2TanOn` | [`circ2d_2tan_on`] |
//! | `GccAna_Circ2d2TanRad` | [`circ2d_2tan_rad`] |
//! | `GccAna_Circ2d3Tan` | [`circ2d_3tan`] |
//! | `GccAna_Circ2dTanOnRad` | [`circ2d_tan_on_rad`] |
//! | `Geom2dGcc_Circ2d2TanRad` | [`geom2d_circ2d_2tan_rad`] |
//! | `GccAna_Circ2dTanCen` | [`circ2d_tan_cen`] |
//! | `GccAna_Lin2d2Tan` | [`lin2d_2tan`] |
//! | `GccAna_Lin2dTanPar` | [`lin2d_tan_par`] |
//! | `GccAna_Lin2dTanPer` | [`lin2d_tan_per`] |
//! | `GccAna_Lin2dTanObl` | [`lin2d_tan_obl`] |
//! | `Geom2dGcc_CurveTool` | [`curve_tool`] |
//! | `Geom2dGcc_Type1` / `Type2` / `Type3` | [`geom2d_gcc_type`] |
//! | `Geom2dGcc_FunctionTanCirCu` | [`function_tan_cir_cu`] |
//! | `Geom2dGcc_FunctionTanCuCuCu` | [`function_tan_cu_cu_cu`] |
//! | `Geom2dGcc_FunctionTanCuCuOnCu` | [`function_tan_cu_cu_on_cu`] |
//! | `Geom2dGcc_Circ2d2TanOnIter` | [`circ2d_2tan_on_iter`] |
//! | `Geom2dGcc_Circ2d2TanOnGeo` | [`circ2d_2tan_on_geo`] |
//! | `Geom2dGcc_Circ2d2TanOn` | [`geom2d_circ2d_2tan_on`] |
//! | `Geom2dGcc_Circ2d3Tan` | [`geom2d_circ2d_3tan`] |
//! | `Geom2dGcc_Circ2dTanCenGeo` | [`circ2d_tan_cen_geo`] |
//! | `Geom2dGcc_Circ2dTanCen` | [`geom2d_circ2d_tan_cen`] |
//! | `Geom2dGcc_FunctionTanCuCu` | [`function_tan_cu_cu`] |
//! | `Geom2dGcc_FunctionTanCuPnt` | [`function_tan_cu_pnt`] |
//! | `Geom2dGcc_Lin2d2TanIter` | [`lin2d_2tan_iter`] |
//! | `Geom2dGcc_Lin2d2Tan` | [`geom2d_lin2d_2tan`] |
//! | `Geom2dGcc_Circ2dTanOnRadGeo` | [`circ2d_tan_on_rad_geo`] |
//! | `Geom2dGcc_Circ2dTanOnRad` | [`geom2d_circ2d_tan_on_rad`] |
//! | `Geom2dGcc_FunctionTanObl` | [`function_tan_obl`] |
//! | `Geom2dGcc_Lin2dTanOblIter` | [`lin2d_tan_obl_iter`] |
//! | `Geom2dGcc_Lin2dTanObl` | [`geom2d_lin2d_tan_obl`] |
//!
//! The first four ports are the ones `ChFi2d_Builder` reaches
//! (`ChFi2d_Builder.cxx:1015-1018`). The line/circle solvers added with the S5
//! wave are the analytic `GccAna` family used by `BRepFill`, `Geom2dGcc` and
//! `GccIter`; they depend only on `gp`, `ElCLib`, `IntAna2d` and `GccEnt`, so
//! they are self-contained.
//!
//! All 29 `Geom2dGcc` classes and all 16 `GccAna` classes now have a port.
//! Two of them are `DEFINE_STANDARD_EXCEPTION` macro headers rather than
//! solvers - `Geom2dGcc_IsParallel` (`hxx:20-33`) and `GccAna_NoSolution` -
//! and reach Rust as the error strings their callers raise, the same way the
//! other OCCT exception classes do (T-490 closed the last solver gap,
//! `Geom2dGcc_Circ2d2TanRadGeo` plus the `Geom2dGcc` front end).

pub mod circ2d_2tan_on;
pub mod circ2d_2tan_on_geo;
pub mod circ2d_2tan_on_iter;
pub mod circ2d_2tan_rad;
pub mod circ2d_2tan_rad_geo;
pub mod circ2d_3tan;
pub mod circ2d_3tan_iter;
pub mod circ2d_tan_cen;
pub mod circ2d_tan_cen_geo;
pub mod circ2d_tan_on_rad;
pub mod circ2d_tan_on_rad_geo;
pub mod curve_tool;
pub mod function_tan_cir_cu;
pub mod function_tan_cu_cu;
pub mod function_tan_cu_cu_cu;
pub mod function_tan_cu_cu_on_cu;
pub mod function_tan_cu_pnt;
pub mod function_tan_obl;
pub mod gcc_ent;
pub mod geom2d_circ2d_2tan_on;
pub mod geom2d_circ2d_2tan_rad;
pub mod geom2d_circ2d_3tan;
pub mod geom2d_circ2d_tan_cen;
pub mod geom2d_circ2d_tan_on_rad;
pub mod geom2d_gcc;
pub mod geom2d_lin2d_2tan;
pub mod geom2d_lin2d_tan_obl;
pub mod geom2d_gcc_type;
pub mod lin2d_2tan;
pub mod lin2d_2tan_iter;
pub mod lin2d_tan_obl;
pub mod lin2d_tan_obl_iter;
pub mod lin2d_tan_par;
pub mod lin2d_tan_per;
pub mod qualified_curve;
pub mod tan_geo_common;

pub use circ2d_2tan_on::GccAnaCirc2d2TanOn;
pub use circ2d_2tan_on_geo::Geom2dGccCirc2d2TanOnGeo;
pub use circ2d_2tan_on_iter::Geom2dGccCirc2d2TanOnIter;
pub use circ2d_3tan_iter::Geom2dGccCirc2d3TanIter;
pub use circ2d_2tan_rad::GccAnaCirc2d2TanRad;
pub use circ2d_2tan_rad_geo::Geom2dGccCirc2d2TanRadGeo;
pub use circ2d_3tan::GccAnaCirc2d3Tan;
pub use circ2d_tan_cen::GccAnaCirc2dTanCen;
pub use circ2d_tan_cen_geo::Geom2dGccCirc2dTanCenGeo;
pub use circ2d_tan_on_rad::GccAnaCirc2dTanOnRad;
pub use circ2d_tan_on_rad_geo::Geom2dGccCirc2dTanOnRadGeo;
pub use function_tan_cir_cu::Geom2dGccFunctionTanCirCu;
pub use function_tan_cu_cu::Geom2dGccFunctionTanCuCu;
pub use function_tan_cu_cu_cu::Geom2dGccFunctionTanCuCuCu;
pub use function_tan_cu_cu_on_cu::Geom2dGccFunctionTanCuCuOnCu;
pub use function_tan_cu_pnt::Geom2dGccFunctionTanCuPnt;
pub use function_tan_obl::Geom2dGccFunctionTanObl;
pub use geom2d_gcc::Geom2dGcc;
pub use geom2d_gcc_type::{Geom2dGccType1, Geom2dGccType2, Geom2dGccType3};
pub use gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};
pub use geom2d_circ2d_2tan_on::Geom2dGccCirc2d2TanOn;
pub use geom2d_circ2d_2tan_rad::Geom2dGccCirc2d2TanRad;
pub use geom2d_circ2d_3tan::Geom2dGccCirc2d3Tan;
pub use geom2d_circ2d_tan_cen::Geom2dGccCirc2dTanCen;
pub use geom2d_circ2d_tan_on_rad::Geom2dGccCirc2dTanOnRad;
pub use geom2d_lin2d_2tan::Geom2dGccLin2d2Tan;
pub use geom2d_lin2d_tan_obl::Geom2dGccLin2dTanObl;
pub use lin2d_2tan::GccAnaLin2d2Tan;
pub use lin2d_2tan_iter::Geom2dGccLin2d2TanIter;
pub use lin2d_tan_obl::GccAnaLin2dTanObl;
pub use lin2d_tan_obl_iter::Geom2dGccLin2dTanOblIter;
pub use lin2d_tan_par::GccAnaLin2dTanPar;
pub use lin2d_tan_per::GccAnaLin2dTanPer;
pub use qualified_curve::{Geom2dGccQCurve, Geom2dGccQualifiedCurve, QualifiedCurve2d};
