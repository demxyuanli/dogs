//! `AddChamfer`, `ModifyChamfer`, `ComputeChamfer`
//! (`ChFi2d_Builder_0.cxx:72-308`).
//!
//! `RemoveChamfer` (`_0.cxx:310-519`) lives in [`super::remove`] because its
//! body is the same as `RemoveFillet`'s.
//!
//! The two `AddChamfer` / `ComputeChamfer` / `ModifyChamfer` overloads are
//! `_d1d2` (two distances) and `_da` (distance + angle).

use crate::chfi2d::{common_vertex, find_connected_edges, ConstructionError};
use crate::shape::Edge;
use crate::topo_tools::is_same;

use super::{is_line_or_circle, ChFi2dBuilder};

impl ChFi2dBuilder {
    /// `AddChamfer(E1, E2, D1, D2)` (`_0.cxx:72-126`): add a chamfer of the
    /// two distances between the two edges that share a vertex.
    ///
    /// `None` is the null edge OCCT returns on every `status != IsDone` arm.
    pub fn add_chamfer_d1d2(&mut self, e1: &Edge, e2: &Edge, d1: f64, d2: f64) -> Option<Edge> {
        let face = self.new_face.clone()?;
        // `ChFi2d::CommonVertex(E1, E2, commonVertex)` (`cxx:80-84`).
        let common_vertex = common_vertex(e1, e2)?;
        // `cxx:86-91`.
        if self.is_a_fillet(e1)
            || self.is_a_chamfer(e1)
            || self.is_a_fillet(e2)
            || self.is_a_chamfer(e2)
        {
            self.status = ConstructionError::NotAuthorized;
            return None;
        }
        // `cxx:93-97`.
        if !is_line_or_circle(e1, &face) || !is_line_or_circle(e2, &face) {
            self.status = ConstructionError::NotAuthorized;
            return None;
        }
        // `status = ChFi2d::FindConnectedEdges(newFace, commonVertex, EE1, EE2)`
        // (`cxx:101-102`).
        let (mut ee1, mut ee2) = match find_connected_edges(&face, &common_vertex) {
            Ok(pair) => pair,
            Err(error) => {
                self.status = error;
                return None;
            }
        };
        // `cxx:103-108`: `EE1` is a copy of `E1` with the good orientation on
        // `newFace`; when it turned out to be `E2` the two are swapped, and the
        // orientation of the discarded copy is carried over to `E2`.
        if is_same(&ee1.0, &e2.0) {
            let orient = ee1.0.orientation();
            ee1 = ee2.clone();
            ee2 = e2.clone();
            ee2.0.set_orientation(orient);
        }
        // `cxx:110`.
        let (e1_mod, e2_mod, chamfer) =
            self.compute_chamfer_d1d2(&common_vertex, &ee1, &ee2, d1, d2)?;
        // `cxx:111-124`.
        if matches!(
            self.status,
            ConstructionError::IsDone
                | ConstructionError::FirstEdgeDegenerated
                | ConstructionError::LastEdgeDegenerated
                | ConstructionError::BothEdgesDegenerated
        ) {
            self.build_new_wire(&ee1, &ee2, &e1_mod, &chamfer, &e2_mod);
            let basis_edge1 = self.basis_edge(&ee1);
            let basis_edge2 = self.basis_edge(&ee2);
            self.update_history(&basis_edge1, &basis_edge2, &e1_mod, &e2_mod, &chamfer, 2);
            self.status = ConstructionError::IsDone;
            return self.chamfers.last().and_then(|s| Edge::wrap(s.clone()));
        }
        Some(chamfer)
    }

    /// `AddChamfer(E, V, D, Ang)` (`_0.cxx:129-176`): add a chamfer of distance
    /// `d` on `E` and angle `ang` with the other edge at `V`.
    pub fn add_chamfer_da(&mut self, e: &Edge, v: &crate::shape::TopoShape, d: f64, ang: f64) -> Option<Edge> {
        let face = self.new_face.clone()?;
        // `cxx:131-135`.
        let (mut adj_edge1, mut adj_edge2) = match find_connected_edges(&face, v) {
            Ok(pair) => pair,
            Err(error) => {
                self.status = error;
                return None;
            }
        };
        // `cxx:139-147`: `adjEdge1` is a copy of `E` with the good orientation
        // on `newFace`.
        if is_same(&adj_edge2.0, &e.0) {
            let orient = adj_edge2.0.orientation();
            adj_edge2 = adj_edge1.clone();
            adj_edge1 = e.clone();
            adj_edge1.0.set_orientation(orient);
        }
        // `cxx:149-154`.
        if self.is_a_fillet(&adj_edge1)
            || self.is_a_chamfer(&adj_edge1)
            || self.is_a_fillet(&adj_edge2)
            || self.is_a_chamfer(&adj_edge2)
        {
            self.status = ConstructionError::NotAuthorized;
            return None;
        }
        // `cxx:156-160`.
        if !is_line_or_circle(&adj_edge1, &face) || !is_line_or_circle(&adj_edge2, &face) {
            self.status = ConstructionError::NotAuthorized;
            return None;
        }
        // `cxx:162-163`.
        let (e1, e2, chamfer) = self.compute_chamfer_da(v, &adj_edge1, d, ang, &adj_edge2)?;
        // `cxx:164-174`.
        if matches!(
            self.status,
            ConstructionError::IsDone
                | ConstructionError::FirstEdgeDegenerated
                | ConstructionError::LastEdgeDegenerated
                | ConstructionError::BothEdgesDegenerated
        ) {
            self.build_new_wire(&adj_edge1, &adj_edge2, &e1, &chamfer, &e2);
            let basis_edge1 = self.basis_edge(&adj_edge1);
            let basis_edge2 = self.basis_edge(&adj_edge2);
            self.update_history(&basis_edge1, &basis_edge2, &e1, &e2, &chamfer, 2);
            self.status = ConstructionError::IsDone;
            return self.chamfers.last().and_then(|s| Edge::wrap(s.clone()));
        }
        Some(chamfer)
    }

    /// `ComputeChamfer(V, E1, E2, D1, D2, TrimE1, TrimE2, Chamfer)`
    /// (`_0.cxx:180-212`).
    pub(super) fn compute_chamfer_d1d2(
        &mut self,
        v: &crate::shape::TopoShape,
        e1: &Edge,
        e2: &Edge,
        d1: f64,
        d2: f64,
    ) -> Option<(Edge, Edge, Edge)> {
        // `cxx:190-194`.
        let (chamfer, new_extr1, new_extr2) = self.build_chamfer_edge_d1d2(v, e1, e2, d1, d2)?;
        if self.status != ConstructionError::IsDone {
            return None;
        }
        self.trim_both(e1, e2, v, &new_extr1, &new_extr2)
            .map(|(t1, t2)| (t1, t2, chamfer))
    }

    /// `ComputeChamfer(V, E1, D, Ang, E2, TrimE1, TrimE2, Chamfer)`
    /// (`_0.cxx:216-248`).
    pub(super) fn compute_chamfer_da(
        &mut self,
        v: &crate::shape::TopoShape,
        e1: &Edge,
        d: f64,
        ang: f64,
        e2: &Edge,
    ) -> Option<(Edge, Edge, Edge)> {
        // `cxx:226-230`.
        let (chamfer, new_extr1, new_extr2) = self.build_chamfer_edge_da(v, e1, d, ang, e2)?;
        if self.status != ConstructionError::IsDone {
            return None;
        }
        self.trim_both(e1, e2, v, &new_extr1, &new_extr2)
            .map(|(t1, t2)| (t1, t2, chamfer))
    }

    /// The tail shared by both `ComputeChamfer` overloads (`cxx:195-209`,
    /// `cxx:231-245`): `BuildNewEdge` twice, then the four degenerated status
    /// combinations.
    ///
    /// `None` where OCCT keeps going with a null trimming edge (the port
    /// reports `ChFi2d_ComputationError` instead, see [`super::fillet`]).
    fn trim_both(
        &mut self,
        e1: &Edge,
        e2: &Edge,
        v: &crate::shape::TopoShape,
        new_extr1: &crate::shape::TopoShape,
        new_extr2: &crate::shape::TopoShape,
    ) -> Option<(Edge, Edge)> {
        let (trim_e1, degen1) = match self.build_new_edge_or_degenerated(e1, v, new_extr1) {
            Some(pair) => pair,
            None => {
                self.status = ConstructionError::ComputationError;
                return None;
            }
        };
        let (trim_e2, degen2) = match self.build_new_edge_or_degenerated(e2, v, new_extr2) {
            Some(pair) => pair,
            None => {
                self.status = ConstructionError::ComputationError;
                return None;
            }
        };
        if degen1 && degen2 {
            self.status = ConstructionError::BothEdgesDegenerated;
        }
        if degen1 && !degen2 {
            self.status = ConstructionError::FirstEdgeDegenerated;
        }
        if !degen1 && degen2 {
            self.status = ConstructionError::LastEdgeDegenerated;
        }
        Some((trim_e1, trim_e2))
    }

    /// `ModifyChamfer(Chamfer, E1, E2, D1, D2)` (`_0.cxx:254-279`):
    /// `RemoveChamfer` then `AddChamfer(E1, E2, D1, D2)` at the vertex the
    /// removal returned, with the two adjacent edges re-oriented so the first
    /// one is not `E2`.
    pub fn modify_chamfer_d1d2(
        &mut self,
        chamfer: &Edge,
        e2: &Edge,
        d1: f64,
        d2: f64,
    ) -> Option<Edge> {
        // `cxx:258`.
        let vertex = self.remove_chamfer(chamfer)?;
        let face = self.new_face.clone()?;
        // `cxx:260-265`.
        let (mut adj_edge1, mut adj_edge2) = match find_connected_edges(&face, &vertex) {
            Ok(pair) => pair,
            Err(error) => {
                self.status = error;
                return None;
            }
        };
        // `cxx:269-275`.
        if is_same(&adj_edge1.0, &e2.0) {
            let orient = adj_edge1.0.orientation();
            adj_edge1 = adj_edge2.clone();
            adj_edge2 = e2.clone();
            adj_edge2.0.set_orientation(orient);
        }
        // `cxx:277-278`.
        self.add_chamfer_d1d2(&adj_edge1, &adj_edge2, d1, d2)
    }

    /// `ModifyChamfer(Chamfer, E, D, Ang)` (`_0.cxx:285-307`):
    /// `RemoveChamfer` then `AddChamfer(E, aVertex, D, Ang)` on whichever of the
    /// two adjacent edges is `E`.
    pub fn modify_chamfer_da(&mut self, chamfer: &Edge, e: &Edge, d: f64, ang: f64) -> Option<Edge> {
        // `cxx:289`.
        let vertex = self.remove_chamfer(chamfer)?;
        let face = self.new_face.clone()?;
        // `cxx:291-296`.
        let (adj_edge1, adj_edge2) = match find_connected_edges(&face, &vertex) {
            Ok(pair) => pair,
            Err(error) => {
                self.status = error;
                return None;
            }
        };
        // `cxx:298-305`.
        if is_same(&adj_edge1.0, &e.0) {
            self.add_chamfer_da(&adj_edge1, &vertex, d, ang)
        } else {
            self.add_chamfer_da(&adj_edge2, &vertex, d, ang)
        }
    }
}
