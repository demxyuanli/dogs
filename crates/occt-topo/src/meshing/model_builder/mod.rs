//! Port of OCCT model_builder — Wave 1 BRepMesh.
//!
//! `BRepMesh_ModelBuilder` turns a `TopoDS_Shape` into a discrete `MeshModel`
//! (faces, wires, edges, 3D curves and 2D pcurves) and
//! `BRepMesh_ModelPreProcessor` initializes per-entity deflection/status before
//! the edge/face discretizers run.
//!
//! Reference: `BRepMesh_ModelBuilder.hxx/.cxx`, `BRepMesh_ModelPreProcessor.hxx/.cxx`,
//! `BRepMeshData_Model.hxx`, `BRepMesh_ShapeVisitor.cxx`, `IMeshData_*.hxx`.
//!
//! Data types (`MeshModel`/`MeshEdge`/`MeshFace`/`MeshWire`/`MeshCurve`/
//! `MeshPCurve`/`MeshStatus`) come from `data_model` and `MeshParameters` from
//! `parameters` — the Wave 1 sibling modules.
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::bnd::{BndBox, BndBox2d};
pub(crate) use occt_core::cslib::{Class2d, Class2dResult};
pub(crate) use occt_core::gp::GpPnt2d;
pub(crate) use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};
pub(crate) use occt_geom::Curve;
pub(crate) use occt_geom2d::curve::Curve2d;

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::bbox_from_geometry::shape_bbox;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::pcurve_full::{make_pcurve_full, project_point_on_surface, project_uv_on_curve2d};
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, faces_of, wires_of_face};

pub(crate) use super::super::data_model::{MeshCurve, MeshEdge, MeshModel, MeshPCurve, MeshStatus};
pub(crate) use super::super::face_discret::PointState;
pub(crate) use super::super::parameters::MeshParameters;
pub(crate) use super::super::range_splitter::{classify_surface, ConeRangeSplitter, RangeSplitter, SurfaceType};
pub(crate) use super::super::shape_tool::ShapeTool;
pub(crate) use super::super::wire_order::{WireOrder, WireOrderStatus};

/// Maximum dimension of a bounding box (`BRepMesh_ShapeTool::BoxMaxDimension`).
pub(crate) fn box_max_dimension(b: &BndBox) -> Option<f64> {
    if b.is_void() {
        return None;
    }
    let (xmin, xmax, ymin, ymax, zmin, zmax) = b.get()?;
    Some((xmax - xmin).max(ymax - ymin).max(zmax - zmin))
}

}

mod wire_builder;
mod preprocessor;
pub use wire_builder::*;
pub use preprocessor::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
