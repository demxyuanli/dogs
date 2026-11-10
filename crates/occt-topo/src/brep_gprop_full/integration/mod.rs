//! Face-integration kernels of `BRepGProp` (`BRepGProp_Gauss` family).
//!
//! Split by OCCT call structure: the `FaceGauss` domain of a face
//! (`BRepGProp_Face` + `BRepGProp_Domain`), the boundary-arc construction from
//! the p-curves, the fixed-order and adaptive 2D Gauss rules, and the 1D
//! slice/scan helpers the adaptive rule calls.

use super::prelude::*;
use super::*;

mod adaptive;
mod arc_scan;
mod boundary;
mod coeff;
mod face_gauss;
mod natural_scan;
mod quadrature;

pub(in crate::brep_gprop_full) use adaptive::*;
pub(in crate::brep_gprop_full) use arc_scan::*;
pub(in crate::brep_gprop_full) use boundary::*;
pub(in crate::brep_gprop_full) use coeff::*;
pub(in crate::brep_gprop_full) use natural_scan::*;
pub(in crate::brep_gprop_full) use quadrature::*;
