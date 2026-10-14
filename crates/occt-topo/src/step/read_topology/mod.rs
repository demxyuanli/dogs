//! STEP topology resolver - the `StepToTopoDS_Translate*` family plus the
//! `StepToGeom` record resolution the STEP read path needs. This is the former
//! monolithic `read_topology.rs`, split by the OCCT class each item carries:
//!
//! - `shape_fix_face_perform` - `ShapeFix_Face::Perform` second-part loop and the
//!   `TopoDS_Iterator` helpers it uses.
//! - `step_to_geom_bspline` - `StepToGeom::MakeBSplineCurveCommon` and its 2D arm.
//! - `resolver_context` - caches, warnings and the `ShapeFix_Shape` heal context.
//! - `translate_tool` - the entity-type dispatch and `TranslateVertex`.
//! - `translate_edge` - `TranslateEdge`.
//! - `translate_edge_loop` - `TranslateEdgeLoop`.
//! - `translate_face` - `TranslateFace`.
//! - `geometric_tool` - `StepToTopoDS_GeometricTool` and `MakePCurve`.
//! - `translate_shell_solid` - `TranslateShell` / `TranslateSolid` /
//!   `StepToTopoDS_Builder`.
//! - `translate_mapped_item` - `TranslateMappedItem` and `MakeTransformation3d`.
//! - `step_to_geom_placement` - `StepToGeom` placements, axes, `ExtractParameter`.
//! - `step_to_geom_curve` - `StepToGeom::MakeCurve`.
//! - `step_to_geom_curve2d` - `StepToGeom::MakeCurve2d`, `MakePCurve`,
//!   `GeomConvert_Units::DegreeToRadian`.
//!
//! The split is a pure move: every `impl<'a> Resolver<'a>` body keeps its
//! statements in order, and items that were module-private are widened to
//! `pub(super)` so the sibling files and this root can reach them.

use super::prelude::*;
use super::*;
// `ShapeBuild_ReShape::Apply` (`ShapeBuild_ReShape.cxx:191`).
use crate::shape_fix_compose_shell::ReShape;

mod geometric_tool;
mod resolver_context;
mod shape_fix_face_perform;
mod step_to_geom_bspline;
mod step_to_geom_curve;
mod step_to_geom_curve2d;
mod step_to_geom_placement;
mod translate_edge;
mod translate_edge_loop;
mod translate_face;
mod translate_mapped_item;
mod translate_shell_solid;
mod translate_tool;

// The impl-only submodules contain nothing but `impl<'a> Resolver<'a>` blocks.
// An inherent `impl` is reachable wherever `Resolver` is - it is
// `pub(in crate::step)` - so those modules need no re-export, and a glob over
// them would be a no-op that trips `unused_imports`. Only the two modules that
// define free items are globbed here; their names reach the sibling files
// through `use super::*`.
use shape_fix_face_perform::*;
use step_to_geom_bspline::*;
