//! Phase 4 module: step — STEP (ISO 10303-21) physical-file exchange.
//!
//! Ports `STEPControl_Writer` / `STEPControl_Reader` for the BRep port.
//! Writes and reads the classic EXPRESS entity set for B-rep solids:
//! `CARTESIAN_POINT`, `DIRECTION`, `VECTOR`, `AXIS2_PLACEMENT_3D`,
//! `LINE`, `CIRCLE`, `ELLIPSE`, `HYPERBOLA`, `PARABOLA`, `POLYLINE`,
//! `B_SPLINE_CURVE`, `B_SPLINE_CURVE_WITH_KNOTS`, `TRIMMED_CURVE`,
//! `OFFSET_CURVE_3D`, `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`,
//! `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`, `B_SPLINE_SURFACE`,
//! `B_SPLINE_SURFACE_WITH_KNOTS`, `VERTEX_POINT`, `EDGE_CURVE`,
//! `ORIENTED_EDGE`, `EDGE_LOOP`, `FACE_OUTER_BOUND`, `FACE_BOUND`,
//! `ADVANCED_FACE`, `CLOSED_SHELL`, `MANIFOLD_SOLID_BREP`, plus the
//! product/representation scaffolding (`PRODUCT`, `PRODUCT_DEFINITION`,
//! `PRODUCT_DEFINITION_SHAPE`, `SHAPE_REPRESENTATION`,
//! `PRODUCT_DEFINITION_SHAPE_REPRESENTATION`,
//! `ADVANCED_BREP_SHAPE_REPRESENTATION`, `NEXT_ASSEMBLY_USAGE_OCCURRENCE`)
//! and the presentation/attribute layer (`COLOUR_RGB`,
//! `SURFACE_STYLE_FILL_AREA`, `SURFACE_STYLE_USAGE`, `STYLED_ITEM`,
//! `SI_UNIT`, `DIMENSIONAL_EXPONENTS`).
//!
//! Curves and surfaces cannot be downcast from `Arc<dyn Curve>` /
//! `Arc<dyn Surface>`, so geometry is classified by sampling invariants
//! (constant zero second derivative ⇒ line, constant curvature + periodic ⇒
//! circle, planar + unbounded ⇒ plane, equidistant samples ⇒ sphere, ...).
//! This mirrors `GeomAdaptor`'s type tag at a slightly higher cost.
//!
//! Non-analytic geometry is written as B-splines: a genuine
//! `GeomBSplineCurve` / `GeomBSplineSurface` is emitted exactly through
//! [`write_bspline_curve`] / [`write_bspline_surface`], and anything else that
//! escapes the analytic classifiers is sampled and re-fitted (see
//! [`fit_bspline_curve`] / [`fit_bspline_surface`]) so the shape-level writers
//! ([`write_step_with_splines`], [`write_step_with_options`]) round-trip
//! arbitrary geometry.
//!
//! The top-level writers mirror `STEPControl_Writer`'s API surface:
//! [`write_step`] / [`write_shape_step`] for plain output,
//! [`write_step_with_splines`] for B-spline coverage, and the attribute
//! variants [`write_step_with_name`], [`write_step_with_color`],
//! [`write_step_with_units`] and [`write_step_assembly`] for the STEP
//! product/presentation layer. Each has a file-writing counterpart and a
//! symmetric reader.
mod prelude {

pub(crate) use std::cell::RefCell;
pub(crate) use std::collections::{HashMap, HashSet};
pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{
    GpAx1, GpAx2, GpAx22d, GpAx2d, GpAx3, GpCirc, GpCirc2d, GpCone, GpCylinder, GpDir, GpDir2d,
    GpElips, GpElips2d, GpHypr, GpLin, GpParab, GpPln, GpPnt, GpPnt2d, GpSphere, GpTorus, GpVec,
    GpVec2d, GpXyz,
};
pub(crate) use occt_geom::{
    bspline_surface::GeomBSplineSurface, Curve, GeomBSplineCurve, GeomCircle, GeomCone,
    GeomCylinder, GeomEllipse, GeomHyperbola, GeomLine, GeomOffsetCurve, GeomOffsetSurface,
    GeomParabola, GeomPlane, GeomRectangularTrimmedSurface, GeomSphere,
    GeomSurfaceOfLinearExtrusion, GeomSurfaceOfRevolution, GeomTorus, GeomTrimmedCurve, Surface,
};
pub(crate) use occt_geom2d::curve::Curve2d;
pub(crate) use occt_geom2d::{
    bspline_curve::Geom2dBSplineCurve, circle::Geom2dCircle, ellipse::Geom2dEllipse,
    line::Geom2dLine, trimmed::Geom2dTrimmedCurve,
};

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::brep_surface::{classify_surface, face_plane, sphere_center, SurfaceKind};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::model::BRepModel;
pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools_full::{edges_of_wire, edge_vertices, wires_of_face};

}


mod p01;
mod p02;
mod p03;
mod p04;
mod p05;
pub use p01::*;
pub use p02::*;
pub use p03::*;
pub use p04::*;
pub use p05::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
