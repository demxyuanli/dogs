//! `NodeInsertionMeshAlgo` implementation for `DelaunayNodeInsertionMeshAlgo`.

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

impl NodeInsertionMeshAlgo for DelaunayNodeInsertionMeshAlgo {
    fn add_node_to_structure(
        &mut self,
        point: GpPnt2d,
        p3d: GpPnt,
        movability: VertexState,
        is_force_add: bool,
    ) -> i32 {
        let vertex = DelaunVertex::new(self.scale_to_face(point), p3d, 0, movability);
        if is_force_add {
            self.structure.add_node_force(vertex)
        } else {
            self.structure.add_node(vertex)
        }
    }

    fn node_point_2d(&self, vertex: &DelaunVertex) -> GpPnt2d {
        // ponytail: no range splitter, so the stored parametric location is
        // returned verbatim (OCCT un-scales it via the splitter).
        vertex.location
    }

    fn insert_boundary_nodes(
        &mut self,
        uv_points: &[GpPnt2d],
        p3d_points: &[GpPnt],
        constraint_edges: &[(i32, i32)],
    ) -> Result<usize, String> {
        self.boundary_indices.clear();
        self.boundary_indices.reserve(uv_points.len());
        for (i, &p) in uv_points.iter().enumerate() {
            // The boundary 3D point is the shared edge polyline vertex — not a
            // placeholder — so adjacent faces meet exactly on the shared edge.
            let node = self.register_node(p3d_points[i], p, VertexState::Frontier, false);
            self.boundary_indices.push(node);
        }

        for &(a, b) in constraint_edges {
            let na = *self.boundary_indices.get(a as usize).ok_or_else(|| {
                format!(
                    "DelaunayNodeInsertionMeshAlgo::insert_boundary_nodes: constraint start {a} out of range ({} points)",
                    uv_points.len()
                )
            })?;
            let nb = *self.boundary_indices.get(b as usize).ok_or_else(|| {
                format!(
                    "DelaunayNodeInsertionMeshAlgo::insert_boundary_nodes: constraint end {b} out of range ({} points)",
                    uv_points.len()
                )
            })?;
            if na != nb {
                self.structure.add_link(na, nb, VertexState::Frontier);
            }
        }

        Ok(uv_points.len())
    }

    fn insert_internal_nodes(&mut self, uv_points: &[GpPnt2d]) -> Result<usize, String> {
        for &p in uv_points {
            self.register_node(GpPnt::zero(), p, VertexState::Free, false);
        }
        Ok(uv_points.len())
    }
}
