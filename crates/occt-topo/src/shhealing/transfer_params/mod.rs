//! Port of `ShapeAnalysis_TransferParameters` and
//! `ShapeAnalysis_TransferParametersProj` (TKShHealing, `ShapeAnalysis`
//! package: `ShapeAnalysis_TransferParameters.cxx` / `.hxx`,
//! `ShapeAnalysis_TransferParametersProj.cxx` / `.hxx`; the family is exactly
//! these two classes, `ShapeAnalysis/FILES.cmake:38-41`).
//!
//! The tool transfers parameters between the 3D curve of an edge and its
//! pcurve. The base class uses the linear map
//! `T2d = myShift + myScale * T3d` (`hxx:34-40`); the `Proj` subclass projects
//! instead whenever the edge is not SameParameter or its tolerance is not below
//! `MaxTolerance()` (`Proj.cxx:112-116`). Call sites ported next to this module:
//! `ShapeFix_Wire::FixNotchedEdges` (`ShapeFix_Wire.cxx:4030-4065`) and, later,
//! `FixSelfIntersectingEdge` / `FixIntersectingEdges` (`cxx:2449`, `:2597`).
//!
//! `ShapeBuild_Edge::CopyRanges` (`ShapeBuild_Edge.cxx:206-334`) is ported here
//! as [`copy_ranges`], together with the `IsPeriodic` unwrap
//! (`ShapeBuild_Edge.cxx:164-204`) and `AdjustByPeriod`
//! (`ShapeBuild_Edge.cxx:148-162`); both `TransferRange` overloads write
//! through it.
//!
//! UNPORTED in this module:
//! * `ShapeAnalysis_TransferParametersProj::CopyNMVertex` (`Proj.cxx:562-710`,
//!   the **edge** overload `(V, toEdge, fromEdge)`): it is not ported because
//!   no caller is ported. Its callers are `ShapeFix_Wire::FixGap3d`
//!   (`ShapeFix_Wire_1.cxx:812-818`, `:838-845`) / `FixGap2d` (`:1958`,
//!   `:1988`), which are reached only from
//!   `ShapeFix_Wireframe::FixWireGaps` (`ShapeFix_Wireframe.cxx:1051`,
//!   `:1617`), plus `ShapeUpgrade_WireDivide.cxx:852`; `ShapeFix_Wireframe`
//!   and `ShapeUpgrade_*` are not ported and are not part of the STEP
//!   `FromSTEP` sequence (which runs only the `FixShape` op,
//!   `STEPControl_Controller.cxx:201`). Its projecting part
//!   (`ShapeAnalysis_Curve::Project` on `C2`, `Proj.cxx:665-671`), its
//!   tolerance update loop (`Proj.cxx:675-709`) and its closing
//!   `BRep_Builder::UpdateVertex(vi, par, E, Tol)` (`Proj.cxx:706`) are
//!   therefore also unported; the vertex-point representation list they walk
//!   (and its writer, [`crate::brep_builder_update_vertex::update_vertex_on_edge`])
//!   does exist in the port. The **face** overload (`Proj.cxx:715-802`) is
//!   ported as [`copy_nm_vertex_face`] for its one call site,
//!   `ShapeFix_ComposeShell::MakeFacesOnPatch` (`cxx:3241`).
//! * `CorrectParameter`'s `Geom2d_BSplineCurve` knot snap (`Proj.cxx:268-279`)
//!   was UNPORTED while `Curve2d` had no knot sequence; **ported in batch 89**
//!   (task T-15) through the new `Curve2d::bspline_knots2d` query — see
//!   [`correct_parameter`].
//! * `myLocation` (`Proj.cxx:97`, `:199`, `:209`, `:318`, `:322`, `:472`): the
//!   port's healing path keeps locations identity, exactly like the rest of
//!   `shhealing` (`wire_fix.rs` reads edges through `BRepTool::edge_curve` /
//!   `curve_on_surface_oriented`, never the `_world` variants), so every
//!   `Transformed(myLocation)` / `Inverted()` is a no-op here.
/// Crate imports shared by `mod.rs` and every submodule.
mod prelude {
    pub(super) use std::sync::Arc;
    pub(super) use occt_core::gp::{GpPnt, GpVec};
    pub(super) use occt_core::precision::{epsilon, Precision, CONFUSION, PCONFUSION, REAL_SMALL};
    pub(super) use occt_geom::{Curve, Surface};
    pub(super) use occt_geom2d::curve::Curve2d;
    pub(super) use crate::boptools_2d::curve_on_surface_oriented;
    pub(super) use crate::brep_tool::BRepTool;
    pub(super) use crate::meshing::edge_discret::CurveOnSurface;
    pub(super) use crate::shape::{Edge, Face, Vertex};
    pub(super) use crate::abs::Orientation;
    pub(super) use crate::brep_point_representation::PointRepresentation;
    pub(super) use crate::pcurve_full::surface_value_of_uv_with_gap;
    pub(super) use crate::tgeometry::GeometryRegistry;
    pub(super) use crate::shhealing::adjust_by_period;
    pub(super) use crate::shhealing::shape_analysis_curve::{next_project, project_adaptor, project_range};
}


// ---------------------------------------------------------------------------
// Submodules
// ---------------------------------------------------------------------------

mod copy_nm_vertex;
mod curve_range;
mod transfer_parameters;
mod transfer_parameters_proj;

pub use self::copy_nm_vertex::copy_nm_vertex_face;
pub use self::transfer_parameters::TransferParameters;
pub use self::transfer_parameters_proj::TransferParametersProj;

