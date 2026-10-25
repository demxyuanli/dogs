//! The inverse functions of the `BRepBlend` package: for a given section plane
//! they solve for the parameters on the guide line and on the supports.
//!
//! Source (OCCT 8.0.0), `ModelingAlgorithms/TKFillet/BRepBlend/`:
//! - [`surf_point_const_rad_inv`]: `BRepBlend_SurfPointConstRadInv`
//!   (`hxx:32-84`, `.cxx:17-276`),
//! - [`surf_point_evol_rad_inv`]: `BRepBlend_SurfPointEvolRadInv`
//!   (`hxx:32-84`, `.cxx:17-269`),
//! - [`surf_curv_const_rad_inv`]: `BRepBlend_SurfCurvConstRadInv`
//!   (`hxx:33-90`, `.cxx:17-284`),
//! - [`surf_curv_evol_rad_inv`]: `BRepBlend_SurfCurvEvolRadInv`
//!   (`hxx:33-90`, `.cxx:17-303`),
//! - [`curv_point_rad_inv`]: `BRepBlend_CurvPointRadInv`
//!   (`hxx:31-80`, `.cxx:17-130`).
//!
//! Each class is a root-finding function (`NbVariables` / `NbEquations` /
//! `Value` / `Values` / `Derivatives` / `IsSolution`) plus the parameter
//! bounds and tolerances the `math_FunctionSetRoot` driver asks for. They are
//! consumed by `BRepBlend_Walking` and `BRepBlend_CSWalking` when a walking
//! step lands on a restriction.

pub mod curv_point_rad_inv;
pub mod surf_curv_const_rad_inv;
pub mod surf_curv_evol_rad_inv;
pub mod surf_point_const_rad_inv;
pub mod surf_point_evol_rad_inv;

pub use curv_point_rad_inv::CurvPointRadInv;
pub use surf_curv_const_rad_inv::SurfCurvConstRadInv;
pub use surf_curv_evol_rad_inv::SurfCurvEvolRadInv;
pub use surf_point_const_rad_inv::SurfPointConstRadInv;
pub use surf_point_evol_rad_inv::SurfPointEvolRadInv;

use occt_core::gp::{GpPnt, GpVec};

/// `P.XYZ()`: the position vector of a point, which OCCT hands to `gp_Vec` /
/// `gp_XYZ` arithmetic as if it were a vector (e.g.
/// `nplan.Dot(point.XYZ())` in `BRepBlend_SurfPointConstRadInv.cxx:114`).
pub(crate) fn pnt_vec(p: &GpPnt) -> GpVec {
    GpVec::from_xyz(p.xyz())
}
