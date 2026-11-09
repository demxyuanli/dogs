//! `Draft::Angle` (`Draft.cxx:34-100`): the draft angle of a face relative to
//! a direction.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `Draft::Angle(F, Direction)` | `cxx:34-100` | [`draft_angle`] |
//!
//! # Translation notes
//!
//! * `BRep_Tool::Surface(F, Lo)` + `S->Transformed(Lo.Transformation())`
//!   (`cxx:38`, `cxx:54`) is [`BRepTool::face_surface_world`], the surface with
//!   the face's location already applied.
//! * The `DynamicType()` test of `cxx:39-49` is the pair
//!   [`Surface::rectangular_trimmed_basis`] (unwraps
//!   `Geom_RectangularTrimmedSurface`) and the `gp_pln` / `gp_cylinder` /
//!   `gp_cone` downcasts; a surface that is none of the three throws
//!   `Standard_DomainError`, which the port answers as `Err`.
//! * `cxx:91` calls `d1u.Cross(d1v)` and discards the result (`gp_Vec::Cross`
//!   is `const`), so the normalisation of `cxx:92` applies to the raw
//!   `u`-derivative. The port keeps that: the cross product is not computed.
//! * `ElSLib::D1` is called at `umin + umax / 2., vmin + vmax / 2.`
//!   (`cxx:90`) - *not* at the mid-parameter of either range; the port keeps
//!   the same expressions.
//! * `D.Dot(axis direction)` on a cylinder or a cone must be `+-1` within
//!   `Precision::Angular()` (`cxx:73`, `cxx:83`), otherwise `DomainError`.

use occt_core::elib::surface_eval::cone_d1;
use occt_core::gp::{GpDir, GpVec};
use occt_core::precision::ANGULAR;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::Face;

/// `gp_Vec` from a `gp_Dir` (OCCT converts implicitly).
fn vec_of(d: &GpDir) -> GpVec {
    GpVec::new(d.x(), d.y(), d.z())
}

/// `Draft::Angle(F, Direction)` (`Draft.cxx:34-100`).
///
/// Answers `Err("Standard_DomainError: Draft::Angle")` where OCCT throws: an
/// unsupported surface type (`cxx:47-50`) or an axis not collinear with
/// `direction` (`cxx:75`, `cxx:85`).
pub fn draft_angle(f: &Face, direction: &GpDir) -> Result<f64, String> {
    let domain_error = || "Standard_DomainError: Draft::Angle".to_string();

    // `cxx:38-45`: the surface, untrimmed.
    let mut s = BRepTool::face_surface_world(f).ok_or_else(domain_error)?;
    if let Some(basis) = s.rectangular_trimmed_basis() {
        s = basis;
    }

    // `cxx:47-50`: only planes, cylinders and cones are supported.
    let plane = s.gp_pln();
    let cyl = s.gp_cylinder();
    let cone = s.gp_cone();
    if plane.is_none() && cyl.is_none() && cone.is_none() {
        return Err(domain_error());
    }

    // `cxx:52-100`: the location is already folded into `s`.
    if let Some(pl) = plane {
        // `cxx:56-67`
        let ax3 = pl.position();
        let mut normale = vec_of(&ax3.direction());
        if !ax3.is_direct() {
            normale.reverse();
        }
        if f.0.orientation() == Orientation::Reversed {
            normale.reverse();
        }
        // `cxx:67`
        return Ok(normale.dot(&vec_of(direction)).asin());
    }

    if let Some(cy) = cyl {
        // `cxx:70-78`
        let testdir = direction.dot(&cy.axis().direction());
        if testdir.abs() <= 1.0 - ANGULAR {
            return Err(domain_error());
        }
        return Ok(0.0);
    }

    // `cxx:79-100`: the conical arm.
    let co = cone.ok_or_else(domain_error)?;
    let testdir = direction.dot(&co.axis().direction());
    if testdir.abs() <= 1.0 - ANGULAR {
        return Err(domain_error());
    }
    // `cxx:88`: `BRepTools::UVBounds(F, umin, umax, vmin, vmax)`.
    let (umin, umax, vmin, vmax) = BRepTool::uv_bounds_world(f);
    // `cxx:90`: `ElSLib::D1(umin + umax / 2., vmin + vmax / 2., Co, ...)`.
    let (_pt, d1u, _d1v) = cone_d1(&co, umin + umax / 2.0, vmin + vmax / 2.0);
    // `cxx:91`: the cross product is discarded by OCCT, see the module notes.
    let mut normal = d1u;
    // `cxx:92`
    normal.normalize();
    // `cxx:93-96`
    if f.0.orientation() == Orientation::Reversed {
        normal.reverse();
    }
    // `cxx:97`
    Ok(normal.dot(&vec_of(direction)).asin())
}
