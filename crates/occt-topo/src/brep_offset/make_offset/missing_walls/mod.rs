//! `BRepOffset_MakeOffset::MakeMissingWalls`
//! (`BRepOffset_MakeOffset.cxx:3148-3687`), the file-local `FillContours`
//! (`:536-599`) and `OrientationOfEdgeInFace` (`:326-343`).
//!
//! The thickening mode finishes the solid by closing the free boundaries of
//! `myFaceComp` with walls: every free-boundary contour (`:3157`) is walked
//! (`:3151`, `:3155`) and, edge by edge, a quad face is built from the original
//! edge, its offset image `OE`, and the two connecting edges `E3` / `E4`
//! (`:3303-3326`). The surface of that face comes from the branch that fits the
//! input (`:3338-3459` cylinder / cone / plane from the two arcs, `:3460-3545`
//! a plane found through `BRepLib_FindSurface` or a thru-sections surface, and
//! finally `:3546-3601` widens `E3` / `E4` through `GeomLib::BuildCurve3d`).
//! The face is then oriented against the original edge (`:3610-3616`) and
//! appended to `myWalls` (`:3618`); an arc this offset generated on the contour
//! vertex closes the corner with a second face (`:3619-3669`).
//!
//! Port notes:
//!
//! * `Message_ProgressScope` (`:3159-3166`) does not exist here; the
//!   `BRepOffset_UserBreak` arm therefore never fires.
//! * `NCollection_IndexedDataMap` (`:3151`) keeps its keys in insertion order;
//!   [`fill_contours`] answers a `Vec` in exactly that order.
//! * `TopoDS_Shape::IsSame` is [`crate::topo_tools::is_same`]; `IsNull` on a
//!   vertex is the `None` of [`edge_vertices`], and every `BRep_Tool::Tolerance`
//!   / `BRep_Tool::Pnt` call on a vertex (`:3215`, `:3228-3229`, `:3241-3242`)
//!   is guarded by that `None` instead of dereferencing a null `TShape`.
//! * `BRep_Builder::UpdateEdge(E, Tol)` (`:3488`, `:3572`, `:3597`) is
//!   [`GeometryRegistry::raise_edge_tolerance`] (`BRep_TEdge::UpdateTolerance`,
//!   raise only), `UpdateEdge(E, C, S, L, Tol)` is
//!   [`GeometryRegistry::update_edge_pcurve_on_surface`], the seam overload
//!   `UpdateEdge(E, C1, C2, S, L, Tol)` is
//!   [`GeometryRegistry::set_edge_pcurves_on_surface`], `UpdateEdge(E, C, Tol)`
//!   is [`GeometryRegistry::set_edge_curve3d`], and `Range(E, S, L, U1, U2)` is
//!   [`GeometryRegistry::set_pcurve_range_on_surface`]. `Range(E, U1, U2)`
//!   (three arguments, `:3222`, `:3574`, `:3599`) writes every representation
//!   ([`GeometryRegistry::set_ranges_all`]) while the four-argument form
//!   `Range(E, U1, U2, Only3d = true)` (`:3524`, `:3538`, `:3541`) writes the
//!   3D range only ([`GeometryRegistry::set_edge_range`]).
//! * `BRepTools::Update(NewFace)` (`:3608`, `:3663`) is `UpdateFaceUVPoints`
//!   (`BRepTools.cxx:488-523`), which refreshes the cached `myUV1/myUV2` of the
//!   curve-on-surface representations. The port stores no such cache (UV points
//!   are evaluated from the pcurve on demand), so there is nothing to refresh.
//! * `BRepLib_MakeEdge(V1, V2)` (`:3277`, `:3293`, `:3647`) fails with
//!   `BRepLib_LineThroughIdenticPoints` when the two points are within
//!   `gp::Resolution()` (`BRepLib_MakeEdge.cxx:185-198`); the port answers
//!   `None` and skips the contour edge.
//! * `BRepLib::SameParameter(NewFace)` (`:3607`) is
//!   [`crate::brep_lib_same_parameter::same_parameter_face`] with the defaults
//!   `Precision::Confusion()`, `forced = false`.
//! * `BRepLib_MakeFace(surf, wire)` (`:3456`, `:3543`) is
//!   [`crate::brep_lib_make_face::make_face_from_surface_and_wire`]; the plane
//!   overload `BRepLib_MakeFace(aPln->Pln(), theWire)` (`:3474`) wraps the plane
//!   in a `Geom_Plane` first, exactly as `BRepLib_MakeFace.cxx` does, and the
//!   `IsDone()` test (`:3475`) is always true for that constructor.
//! * `BRepLib_MakeFace(arcWire, true)` (`:3662`) is
//!   [`crate::brep_lib_make_face::make_face_from_wire`] with `only_plane`; the
//!   reference appends its answer unconditionally (`:3664`), so a wire for which
//!   `BRepLib_FindSurface` finds no plane appends a *null* face to `myWalls`
//!   (a shell cannot hold it either). The port has no null `TopoDS_Face` and
//!   skips that append.
//! * `BRepAdaptor_Curve(E)` (`:3329-3330`) of an edge without a face is its 3D
//!   curve ([`crate::brep_tool::BRepTool::edge_curve`]); `GetType() ==
//!   GeomAbs_Circle` is `Curve::gp_circ` answering `Some`.
//! * `GC_MakeCylindricalSurface(aCirc)` (`:3349`) is
//!   `gp_Cylinder(Circ.Position(), Circ.Radius())` (`gce_MakeCylinder.cxx`,
//!   `gp_Circ` overload) and `gce_MakeCone(P1, P2, R1, R2)` (`:3406-3407`) is
//!   reproduced line for line in [`gce_make_cone`] (`gce_MakeCone.cxx:234-283`).
//! * `BRep_Builder::MakeFace(F, PL, Tol)` (`:3387`, `:3399`) is
//!   [`TopoBuilder::make_face_plane`]; `BB.Add(theWire, E)` is
//!   [`TopoBuilder::add_edge`], which stores the edge with the orientation the
//!   occurrence already carries.
//! * `TopoDS_Shape::Reversed()` (`:3312-3322`, `:3635`, `:3638`, `:3652`) and
//!   `Oriented(...)` (`:3438-3441`, `:3665`) compose the request with the
//!   occurrence's own orientation (`TopoDS_Shape::Orientation(Compose(Or,
//!   myOrient))`); the port's [`TopoShape::oriented`] only *replaces* it, so
//!   every such call goes through [`oriented_by`]. `TopoDS_Shape::Reverse()`
//!   (`:3295`, `:3299-3300`, `:3635`, `:3638`) and `Orientation(...)`
//!   (`:3209`, `:3223`, `:3272`, `:3438`, `:3649`, `:3654`, `:3615`) are plain
//!   writes and map to `TopoShape::reverse` / `set_orientation`.
//! * `myWalls.Append(NewFace)` (`:3618`).
//!   `:3178`, `:3277`, `:3621`, `:3674`) throw `Standard_TypeMismatch` on a
//!   wrong type; the port skips the shape.
//! * `BRepAlgo_Image::LastImage` on a shape with no image (`:3188-3189`) and
//!   `LOE.Last()` on an empty list throw in the reference; the port skips the
//!   contour edge.
//! * `gp_Dir OffsetDir = gce_MakeDir(PonE, PonOE)` (`:3335`) raises
//!   `gp_VectorWithNullMagnitude` on confused points; the port's
//!   `GpDir::from_vec` answers `Err` and the contour edge is skipped.

mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::sync::Arc;
pub(crate) use occt_core::gp::{
    GpAx2, GpAx3, GpCone, GpCylinder, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpVec,
};
pub(crate) use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL, RESOLUTION, SQUARE_CONFUSION};
pub(crate) use occt_geom::cone::GeomCone;
pub(crate) use occt_geom::curve::Curve;
pub(crate) use occt_geom::cylinder::GeomCylinder;
pub(crate) use occt_geom::geom_fill::generator::generator_surface;
pub(crate) use occt_geom::geom_lib::geom_lib_build_curve3d;
pub(crate) use occt_geom::plane::GeomPlane;
pub(crate) use occt_geom::surface::Surface;
pub(crate) use occt_geom::trimmed::GeomTrimmedCurve;
pub(crate) use occt_geom2d::curve::Curve2d;
pub(crate) use occt_geom2d::line::Geom2dLine;
pub(crate) use occt_core::kernel::geomabs::Shape as Continuity;
pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::brep_lib_build_curves3d::build_curves3d_tol;
pub(crate) use crate::brep_lib_find_surface::FindSurface;
pub(crate) use crate::brep_lib_make_face::{make_face_from_surface_and_wire, make_face_from_wire};
pub(crate) use crate::brep_lib_same_parameter::same_parameter_face;
pub(crate) use crate::brep_offset::make_offset::tolerance_update::compute_max_dist;
pub(crate) use crate::brep_offset::make_offset::BRepOffsetMakeOffset;
pub(crate) use crate::brep_offset::Analyse;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::chfi_ds::TypeOfConcavity;
pub(crate) use crate::fclass2d::FClass2d;
pub(crate) use crate::meshing::edge_discret::CurveOnSurface;
pub(crate) use crate::meshing::model_builder::wire_builder::wire_edges_explorer;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools::{
    children_keep_orientation, edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey,
};

}

mod statics;
mod fill_contours;
mod perform;
mod rebuild_connecting_edge;
