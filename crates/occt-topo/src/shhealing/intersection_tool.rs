//! `ShapeFix_IntersectionTool` (`ShapeFix_IntersectionTool.cxx`).
//!
//! Ported: `SplitEdge` (`cxx:85-189`), `CutEdge` (`cxx:194-268`),
//! `SplitEdge1` (`cxx:270-356`), `SplitEdge2` (`cxx:363-490`),
//! `UnionVertexes` (`cxx:495-880`), the static `CreateBoxes2d` (`cxx:884-920`)
//! and `SelectIntPnt` (`cxx:923-962`), `FindVertAndSplitEdge`
//! (`cxx:967-1025`) and `FixSelfIntersectWire` (`cxx:1028-1830`).
//!
//! `FixIntersectingWires` (`cxx:1834-2530`), the *face*-level twin, lives in
//! the `fix_intersecting_wires` submodule (reachable as
//! `shhealing::fix_intersecting_wires`); the helpers here are shared with it.
//!
//! This is the non-adjacent self-intersection driver behind
//! `ShapeFix_Wire::FixSelfIntersection`'s tail
//! (`ShapeFix_Wire.cxx:1200-1277`): it walks every pair of non-adjacent edges
//! of one wire, detects a 2D box overlap, intersects the two pcurves and then
//! splits / cuts / welds the pair. `ShapeFix_Wire` keeps the modified edges in
//! the `ShapeExtend_WireData` it already owns, so the port mutates `wire` in
//! place with `wire_set_edge_composed` / `wire_insert_edge_before` /
//! `wire_remove_edge`, exactly like the wire actions in `wire_fix.rs`.
//!
//! The `ShapeBuild_ReShape` context is only used to record substitutions
//! (`myContext->Replace`, `cxx:293`, `:426`, `:527`, ...); the wire itself
//! carries the new edges, so the port keeps a caller-provided `&mut dyn ReShape`
//! and never needs `Context()->Apply` beyond the face's own vertices.
//!
//! `ShapeFix_Face::FixIntersectingWires` (`ShapeFix_Face.cxx:2821-2825`), the
//! driver, is carried as `ShapeFixFace::fix_intersecting_wires`; wiring it into
//! `ShapeFix_Face::Perform` still needs the unported wire-fixing first part
//! (`ShapeFix_Face.cxx:365-480`).

mod fix_intersecting_wires;
pub use fix_intersecting_wires::*;

use super::prelude::*;

use occt_core::bnd::BndBox2d;
use occt_core::intres2d::{IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dPosition};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::Orientation;
use crate::bnd_lib_add2d::{add_geom2d, add_geom2d_range};
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::shape_fix_compose_shell::ReShape;
use crate::shhealing::transfer_params::TransferParametersProj;
use crate::shhealing::wire_fix::{
    copy_pcurves, copy_replace_vertices_with, first_vertex, last_vertex, wire_insert_edge_before,
    wire_remove_edge, wire_set_edge_composed,
};
use crate::shhealing::wire_self_inter::point_on_edge;

mod helpers;
mod edge_split;
mod vertex_union;
mod self_intersect;

pub use helpers::SHAPE_FIX_INTERSECTION_MAX_TOL;
use helpers::*;
use edge_split::*;
use vertex_union::*;
pub use helpers::FixSelfIntersectWireStatus;
pub use self_intersect::*;

