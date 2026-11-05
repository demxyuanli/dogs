//! `AddFillet`, `ModifyFillet`, `ComputeFillet`, `IsAFillet`, `IsAChamfer`
//! (`ChFi2d_Builder.cxx:243-283`, `:504-530`, `:1233-1265`).
//!
//! `RemoveFillet` (`cxx:290-497`) lives in [`super::remove`] because its body is
//! the same as `RemoveChamfer`'s.

use crate::chfi2d::{find_connected_edges, ConstructionError};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::is_same;

use super::{is_line_or_circle, ChFi2dBuilder};

impl ChFi2dBuilder {
    /// `AddFillet(V, Radius)` (`cxx:243-278`): add a fillet of `radius`
    /// between the two edges of the wire meeting at `V`.
    ///
    /// `None` is the null edge OCCT returns on every `status != IsDone` arm.
    pub fn add_fillet(&mut self, v: &TopoShape, radius: f64) -> Option<Edge> {
        let face = self.new_face.clone()?;
        // `cxx:247-251`.
        let (adj_edge1, adj_edge2) = match find_connected_edges(&face, v) {
            Ok(pair) => pair,
            Err(error) => {
                self.status = error;
                return None;
            }
        };
        // `cxx:253-257`.
        if self.is_a_fillet(&adj_edge1)
            || self.is_a_chamfer(&adj_edge1)
            || self.is_a_fillet(&adj_edge2)
            || self.is_a_chamfer(&adj_edge2)
        {
            self.status = ConstructionError::NotAuthorized;
            return None;
        }
        // `cxx:259-263`.
        if !is_line_or_circle(&adj_edge1, &face) || !is_line_or_circle(&adj_edge2, &face) {
            self.status = ConstructionError::NotAuthorized;
            return None;
        }
        // `cxx:265-266`.
        let (adj_edge1_mod, adj_edge2_mod, fillet) =
            self.compute_fillet(v, &adj_edge1, &adj_edge2, radius)?;
        // `cxx:267-277`.
        if matches!(
            self.status,
            ConstructionError::IsDone
                | ConstructionError::FirstEdgeDegenerated
                | ConstructionError::LastEdgeDegenerated
                | ConstructionError::BothEdgesDegenerated
        ) {
            self.build_new_wire(&adj_edge1, &adj_edge2, &adj_edge1_mod, &fillet, &adj_edge2_mod);
            let basis_edge1 = self.basis_edge(&adj_edge1);
            let basis_edge2 = self.basis_edge(&adj_edge2);
            self.update_history(&basis_edge1, &basis_edge2, &adj_edge1_mod, &adj_edge2_mod, &fillet, 1);
            self.status = ConstructionError::IsDone;
            return self.fillets.last().and_then(|s| Edge::wrap(s.clone()));
        }
        Some(fillet)
    }

    /// `ModifyFillet(Fillet, Radius)` (`cxx:282-287`): `RemoveFillet` then
    /// `AddFillet` at the vertex the removal returned.
    pub fn modify_fillet(&mut self, fillet: &Edge, radius: f64) -> Option<Edge> {
        let vertex = self.remove_fillet(fillet)?;
        self.add_fillet(&vertex, radius)
    }

    /// `ComputeFillet(V, E1, E2, Radius, TrimE1, TrimE2, Fillet)`
    /// (`cxx:504-530`): the fillet edge plus the two trimmed adjacent edges.
    ///
    /// `None` on any arm where OCCT leaves `status != IsDone` and returns with
    /// null out-parameters.
    pub(super) fn compute_fillet(
        &mut self,
        v: &TopoShape,
        e1: &Edge,
        e2: &Edge,
        radius: f64,
    ) -> Option<(Edge, Edge, Edge)> {
        // `cxx:511-515`.
        let (fillet, new_extr1, new_extr2) = self.build_fillet_edge(v, e1, e2, radius)?;
        if self.status != ConstructionError::IsDone {
            return None;
        }
        // `cxx:516-518`: `BuildNewEdge(E1, V, NewExtr1, Degen1)`.
        //
        // A projection failure leaves OCCT with a null trimming edge that it
        // then adds to the wire (`cxx:520-530` has no check). The port cannot
        // store a null shape, so the failure is reported as
        // `ChFi2d_ComputationError`.
        let (trim_e1, degen1) = match self.build_new_edge_or_degenerated(e1, v, &new_extr1) {
            Some(pair) => pair,
            None => {
                self.status = ConstructionError::ComputationError;
                return None;
            }
        };
        let (trim_e2, degen2) = match self.build_new_edge_or_degenerated(e2, v, &new_extr2) {
            Some(pair) => pair,
            None => {
                self.status = ConstructionError::ComputationError;
                return None;
            }
        };
        // `cxx:519-529`.
        if degen1 && degen2 {
            self.status = ConstructionError::BothEdgesDegenerated;
        }
        if degen1 && !degen2 {
            self.status = ConstructionError::FirstEdgeDegenerated;
        }
        if !degen1 && degen2 {
            self.status = ConstructionError::LastEdgeDegenerated;
        }
        Some((trim_e1, trim_e2, fillet))
    }

    /// `IsAFillet(E)` (`cxx:1233-1248`).
    pub(super) fn is_a_fillet(&self, e: &Edge) -> bool {
        self.fillets.iter().any(|s| is_same(s, &e.0))
    }

    /// `IsAChamfer(E)` (`cxx:1252-1265`).
    pub(super) fn is_a_chamfer(&self, e: &Edge) -> bool {
        self.chamfers.iter().any(|s| is_same(s, &e.0))
    }
}
