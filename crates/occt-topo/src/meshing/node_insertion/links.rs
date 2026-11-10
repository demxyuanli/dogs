//! `BRepMesh_BaseMeshAlgo` link helpers.


use occt_core::precision::SQUARE_CONFUSION;

use crate::abs::Orientation;

use super::super::delaun_types::VertexState;
#[allow(unused_imports)]
use super::super::mesh_algo::*;
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
