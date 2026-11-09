//! `ShapeAnalysis_Edge::CheckSameParameter` (`ShapeAnalysis_Edge.cxx:694-835`).
//!
//! Measures the maximal deviation between an edge's 3D curve and its p-curve
//! lifted onto the face surface, over `NbControl` equidistant points, using
//! [`ValidateEdge`] exactly as the reference does.
//!
//! Boundary of the port:
//! - `ShapeExtend`'s status bits reduce to the `done` flag the caller reads:
//!   `Status(ShapeExtend_DONE)` is true when `maxdev` exceeds the edge
//!   tolerance (`cxx:828-831`) or the edge is flagged not-same-parameter
//!   (`cxx:832-835`). `ShapeExtend_FAIL1` / `FAIL2` are not DONE bits, so the
//!   degenerated / null-curve / failed-`Process` exits answer `None` or leave
//!   `done` unset, which is the reference's `false`.
//! - `cxx:763-782` applies the curve's `TopLoc_Location` to the 3D curve; this
//!   port's `BRep_Tool::Curve` equivalent is
//!   [`BRepTool::edge_curve_world`], which already folds it.
//! - the per-representation loop `cxx:786-824` walks every p-curve of the edge
//!   and keeps those whose `(surface, location)` matches the face. This port
//!   keys p-curves by face, so the equivalent single read is
//!   [`crate::boptools_2d::curve_on_surface_range`] for that face - the same
//!   set the filter lets through.
//! - `cxx:826-852` (`CurveOnPlane` when no p-curve matched) is the arm
//!   `curve_on_surface_range` already performs internally.

use occt_geom::curve_on_surface::CurveOnSurface;

use crate::boptools_2d;
use crate::brep_lib_validate_edge::ValidateEdge;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};

/// `CheckSameParameter`'s out-parameter plus its return value.
pub struct SameParameterCheck {
    /// `maxdev` (`cxx:715`), the largest deviation found.
    pub maxdev: f64,
    /// `return Status(ShapeExtend_DONE)` (`cxx:835`).
    pub done: bool,
}

/// `CheckSameParameter(edge, face, maxdev, NbControl = 23)` (`cxx:704-836`).
///
/// `None` is the reference's `false` before `maxdev` is meaningful: a
/// degenerated edge (`cxx:710-713`) or an edge with no 3D curve
/// (`cxx:725-729`). `face` is `None` for the `TopoDS_Face`-less overload
/// (`cxx:694-700`), where every p-curve is measured.
pub fn check_same_parameter(edge: &Edge, face: Option<&Face>) -> Option<SameParameterCheck> {
    // `cxx:710-713`.
    if BRepTool::is_degenerated(edge) {
        return None;
    }

    let mut maxdev = 0.0f64;

    // `cxx:719`: `BRep_TEdge::SameParameter()`.
    let same_parameter = BRepTool::same_parameter(edge);

    // `cxx:722-729`: the 3D curve, location already folded in.
    let c3d = BRepTool::edge_curve_world(edge)?;

    // `cxx:768-773`: the face surface, when a face was given.
    let face_surface = face.and_then(BRepTool::face_surface_world);

    // `cxx:786-824`: the p-curve whose surface is the face's.
    if let Some(f) = face {
        if let Some(surf) = &face_surface {
            if let Some((pc, _pf, _pl)) = boptools_2d::curve_on_surface_range(edge, f) {
                // `cxx:812-818`: `Geom2dAdaptor_Curve(aPC, f, l)`,
                // `GeomAdaptor_Surface(aST)` and
                // `Adaptor3d_CurveOnSurface(GHPC, GAHS)`.
                let acs = CurveOnSurface::new(pc.clone(), surf.clone());
                let mut a_validate = ValidateEdge::new(
                    c3d.as_ref(),
                    &acs,
                    Some(surf.as_ref()),
                    Some(pc.as_ref()),
                    same_parameter,
                );
                // `cxx:819`: `SetControlPointsNumber(NbControl - 1)` with the
                // `NbControl = 23` default of `hxx:225`.
                a_validate.set_control_points_number(23 - 1);
                a_validate.process();
                a_validate.update_tolerance(&mut maxdev);
                // `cxx:821-824`: a failed measure is `FAIL2`, not DONE.
            }
        }
    }

    // `cxx:828-835`.
    let mut done = false;
    if maxdev > BRepTool::edge_tolerance(edge) {
        done = true;
    }
    if !same_parameter {
        done = true;
    }
    Some(SameParameterCheck { maxdev, done })
}
