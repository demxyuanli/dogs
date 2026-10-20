//! `BRepMesh_BaseMeshAlgo` link helpers.

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
    /// Source: `BRepMesh_BaseMeshAlgo::addLinkToMesh`.
    pub(super) fn add_link_to_mesh(&mut self, first: i32, last: i32, orientation: Orientation) -> i32 {
        let id = if orientation == Orientation::Reversed {
            self.structure
                .add_link(last, first, VertexState::Frontier)
        } else if orientation == Orientation::Internal {
            self.structure.add_link(first, last, VertexState::Fixed)
        } else {
            self.structure
                .add_link(first, last, VertexState::Frontier)
        };
        id.abs()
    }

    /// If another pcurve of the same edge on this face coincides in UV, the
    /// seam is `INTERNAL` (`Fixed`). Else the pcurve's own orientation.
    /// Source: `BRepMesh_BaseMeshAlgo::fixSeamEdgeOrientation`.
    pub(super) fn fix_seam_edge_orientation(
        edge: &MeshEdge,
        pcurve_index: usize,
        face_index: usize,
    ) -> Result<Orientation, String> {
        let pc = edge.pcurve(pcurve_index)?;
        let pts = pc.points();
        if pts.is_empty() {
            return Ok(pc.orientation());
        }
        let pnt1_1 = pts[0];
        let pnt2_1 = pts[pts.len() - 1];
        for other_index in edge.pcurves_for(face_index) {
            if other_index == pcurve_index {
                continue;
            }
            let other = edge.pcurve(other_index)?;
            let other_pts = other.points();
            if other_pts.is_empty() {
                continue;
            }
            let pnt1_2 = other_pts[0];
            let pnt2_2 = other_pts[other_pts.len() - 1];
            let sq1 = pnt1_1
                .square_distance(&pnt1_2)
                .min(pnt1_1.square_distance(&pnt2_2));
            let sq2 = pnt2_1
                .square_distance(&pnt1_2)
                .min(pnt2_1.square_distance(&pnt2_2));
            if sq1 < SQUARE_CONFUSION && sq2 < SQUARE_CONFUSION {
                return Ok(Orientation::Internal);
            }
        }
        Ok(pc.orientation())
    }

}
