//! `BRepApprox_TheMultiLineOfApprox` — the `BRepApprox` instantiation of the
//! `ApproxInt_MultiLine` template.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepApprox/BRepApprox_TheMultiLineOfApprox_0.cxx:16-29`
//! (OCCT 8.0.0), which only binds the template names:
//!
//! ```text
//! #define Handle_TheLine      occ::handle<BRepApprox_ApproxLine>
//! #define TheLine             BRepApprox_ApproxLine
//! #define TheSvSurfaces      ApproxInt_SvSurfaces
//! #define ApproxInt_MultiLine BRepApprox_TheMultiLineOfApprox
//! #include <ApproxInt_MultiLine.gxx>
//! ```
//!
//! The template body is [`crate::approx_int::multi_line`], the line is
//! [`ApproxLine`] and the surface pair is any [`SvSurfaces`] implementation
//! (`BRepApprox_ThePrmPrmSvSurfacesOfApprox` / `…TheImpPrm…`).

use crate::approx_int::MultiLine;

use super::ApproxLine;

/// `BRepApprox_TheMultiLineOfApprox`.
pub type BRepMultiLine = MultiLine<ApproxLine>;
