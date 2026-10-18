//! IGES (ANSI Y14.26M) B-Rep writer — a port of `IGESControl_Writer`.
//!
//! Emits an 80-column fixed-width IGES file with the entity set OCCT's default
//! Faces mode (`write.iges.brep.mode = 0`, `IGESData.cxx:90-94`) produces for a
//! B-rep: 116 POINT, 110 LINE, 100 CIRCULAR ARC, 104 CONIC ARC, 108 PLANE, 120
//! SURFACE OF REVOLUTION, 122 TABULATED CYLINDER, 124 TRANSFORMATION MATRIX, 126
//! B-SPLINE CURVE, 128 B-SPLINE SURFACE, 140 OFFSET SURFACE, 142 CURVE ON
//! SURFACE, 144 TRIMMED SURFACE and 402 GROUP. Elementary surfaces (cylinder,
//! cone, sphere, torus) become 120, **not** the 192/194/196/198 solids — those
//! are the BRep-mode arm
//! (`myBRepMode && myAnalytic`) and are not selected here.
//!
//! Curves and surfaces cannot be downcast from `Arc<dyn Curve>` /
//! `Arc<dyn Surface>`, so the `Surface` / `Curve` traits expose the exact
//! analytic type (`gp_pln`, `gp_cylinder`, … / `gp_circ`, …) and every dispatch
//! is a type query, mirroring `GeomToIGES_GeomSurface::TransferSurface`'s
//! `IsKind` chain. A face whose basis surface has no branch in that chain is
//! dropped (`BRepToIGES_BRShell.cxx:257-261`), never replaced by a plane.
//!
//! The physical layout is the standard sectioned format (S/G/D/P/T) with each
//! line exactly 80 characters. Section lines carry their section letter and a
//! 7-digit sequence number in the first 8 columns, followed by the data.
//!
//! **Model contents (T-85 step 2)**: only the entities reachable from the
//! shapes handed to the writer are written, numbered in
//! `Interface_InterfaceModel::AddWithRefs` order (`Interface_InterfaceModel.cxx:652-692`),
//! exactly as `IGESControl_Writer::AddEntity` → `myModel->AddWithRefs`
//! (`IGESControl_Writer.cxx:243-252`) does. Parameter pointers therefore cannot
//! be written as they are built: the emitters record them as `#k` placeholders
//! ([`emit_refs`](IgesWriter::emit_refs)) and
//! [`final_entities`](IgesWriter::final_entities) resolves them after the
//! renumbering.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{
    GpAx2, GpAx2d, GpAx22d, GpAx3, GpCirc, GpCirc2d, GpDir, GpDir2d, GpElips, GpElips2d, GpHypr,
    GpHypr2d, GpLin, GpLin2d, GpMat2d, GpParab, GpParab2d, GpPln, GpPnt, GpPnt2d, GpTrsf,
    GpTrsf2d, GpVec, GpVec2d, GpXY, GpXyz, TrsfForm,
};
pub(crate) use occt_geom::{
    Curve, GeomBSplineCurve, GeomBezierCurve, GeomCircle, GeomEllipse, GeomHyperbola, GeomLine,
    GeomParabola, Surface,
};
pub(crate) use occt_geom2d::curve::Curve2d;
pub(crate) use occt_geom2d::trimmed::Geom2dTrimmedCurve;
pub(crate) use occt_geom2d::{
    Geom2dBSplineCurve, Geom2dBezierCurve, Geom2dCircle, Geom2dEllipse, Geom2dHyperbola,
    Geom2dLine, Geom2dParabola,
};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::brep_surface::{classify_surface, SurfaceKind};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::model::BRepModel;
pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex};
pub(crate) use crate::topo_tools::{edge_vertices, edges_of_wire, wires_of_face};
}

use prelude::*;

mod card;
mod curve_kind;
mod emit_conic;
mod emit_curve;
mod emit_face;
mod emit_surface;
mod entity;
mod entry;
mod frame_geometry;
mod pcurve_transfer;
mod wire_helpers;
mod writer;

pub use entry::*;
use card::*;
use curve_kind::*;
use emit_conic::*;
use emit_curve::*;
use emit_face::*;
use emit_surface::*;
use entity::*;
use frame_geometry::*;
use pcurve_transfer::*;
use wire_helpers::*;
use writer::*;

#[cfg(test)]
mod tests;
