//! `BRepMesh_NodeInsertionMeshAlgo` face-basis scaling, `AddPoint` and `GenerateSurfaceNodes`.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d, GpXY};
use occt_core::precision::SQUARE_CONFUSION;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::abs::Orientation;

use super::super::data_model::*;
use super::super::delaun::Delaun;
use super::super::delaun_data::DelaunDataStructure;
use super::super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::super::face_discret::Classifier;
#[allow(unused_imports)]
use super::super::mesh_algo::*;
use super::super::face_discret::PointState;
use super::super::deflection_control::{
    DelaunayDeflectionControlMeshAlgo, DeflectionMesh,
};
use super::super::parameters::MeshParameters;
use super::super::range_splitter::{
    mesh_surface_type, create_range_splitter, factory_uses_deflection_control,
    grab_params_of_internal_edges, RangeSplitter, MeshSurfaceType,
};
use super::*;

use super::*;

impl DelaunayNodeInsertionMeshAlgo {
    /// Maps a UV point into the Delaunay face basis. Source:
    /// `BRepMesh_NodeInsertionMeshAlgo::addNodeToStructure` /
    /// `BRepMesh_DefaultRangeSplitter::Scale(thePoint, true)`.
    pub(super) fn scale_to_face(&self, point: GpPnt2d) -> GpPnt2d {
        let Some(s) = self.splitter.as_ref() else {
            return point;
        };
        let (du, dv) = s.delta();
        if du.abs() <= 1e-16 || dv.abs() <= 1e-16 {
            return point;
        }
        s.scale(point, true)
    }

    /// Maps a face-basis point back to surface UV. Source:
    /// `BRepMesh_NodeInsertionMeshAlgo::getNodePoint2d` /
    /// `Scale(thePoint, false)`.
    pub(super) fn scale_from_face(&self, point: GpPnt2d) -> GpPnt2d {
        let Some(s) = self.splitter.as_ref() else {
            return point;
        };
        s.scale(point, false)
    }

    /// Registers a point in the 3D nodes map and adds its 2D projection to the
    /// structure. Returns the structure node index. Source:
    /// `BRepMesh_BaseMeshAlgo::registerNode`.
    pub fn register_node(
        &mut self,
        point: GpPnt,
        point2d: GpPnt2d,
        movability: VertexState,
        is_force_add: bool,
    ) -> i32 {
        let node_index =
            self.add_node_to_structure(point2d, point, movability, is_force_add);
        if node_index as usize > self.nodes_map.len() {
            self.nodes_map.push(point);
            // Pre-bind frontier/fixed nodes with identity mapping so boundary
            // node indices stay stable; Free (internal) nodes enter only when
            // actually referenced by a triangle.
            if movability != VertexState::Free {
                self.used_nodes.insert(node_index, node_index);
            }
        }
        node_index
    }

    /// Generates surface (interior) nodes for the face and registers them.
    /// Source: `DelaunayNodeInsertionMeshAlgo::registerSurfaceNodes`, driving the
    /// analytical `RangeSplitter` (cylinder/cone/sphere/torus/NURBS) to place the
    /// interior UV nodes that a periodic surface's degenerate seam boundary alone
    /// cannot cover.
    pub fn generate_surface_nodes(
        &mut self,
        model: &MeshModel,
        face_index: usize,
        params: &MeshParameters,
    ) -> Result<usize, String> {
        let nodes = self.list_surface_nodes(model, face_index, params)?;
        let count = nodes.len();
        for (uv, p3d) in nodes {
            self.register_node(p3d, uv, VertexState::Free, false);
        }
        Ok(count)
    }

    /// Surface nodes classified `IN` (`registerSurfaceNodes` / `insertNodes` loop).
    /// Source: `BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:115-123, 144-152`.
    pub(super) fn list_surface_nodes(
        &mut self,
        model: &MeshModel,
        face_index: usize,
        params: &MeshParameters,
    ) -> Result<Vec<(GpPnt2d, GpPnt)>, String> {
        let face = model.face(face_index)?;
        let Some(surface) = face.surface() else {
            return Ok(Vec::new());
        };
        let nodes = if let Some(splitter) = self.splitter.as_mut() {
            grab_params_of_internal_edges(model, face_index, splitter.as_mut());
            splitter.generate_surface_nodes(params)
        } else {
            let mut splitter = create_range_splitter(surface.as_ref());
            splitter.reset(face, params);
            for &uv in &self.boundary_uv {
                splitter.add_point(uv);
            }
            splitter.adjust_range();
            grab_params_of_internal_edges(model, face_index, splitter.as_mut());
            splitter.generate_surface_nodes(params)
        };
        let Some(nodes) = nodes else {
            return Ok(Vec::new());
        };
        let mut out = Vec::new();
        let cls_dbg = std::env::var("CLSNODE").is_ok();
        let (mut n_in, mut n_out, mut n_on) = (0usize, 0usize, 0usize);
        for n in nodes {
            let st = self.classifier.perform(&n);
            match st {
                PointState::In => n_in += 1,
                PointState::Out => n_out += 1,
                PointState::On => n_on += 1,
            }
            if cls_dbg && st != PointState::In {
                eprintln!(
                    "CLSNODE face={face_index} st={:?} uv=({:.12e},{:.12e})",
                    st,
                    n.x(),
                    n.y()
                );
            }
            if st == PointState::In {
                out.push((n, surface.d0(n.x(), n.y())));
            }
        }
        if cls_dbg {
            eprintln!(
                "CLSNODE face={face_index} wires={} total={} in={n_in} out={n_out} on={n_on}",
                self.classifier.wires_nb(),
                n_in + n_out + n_on
            );
        }
        Ok(out)
    }

}
