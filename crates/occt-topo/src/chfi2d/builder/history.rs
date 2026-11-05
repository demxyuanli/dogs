//! History maintenance: `UpDateHistory` (`ChFi2d_Builder.cxx:684-745`) and
//! `BasisEdge` (`Builder.cxx:747-762`).
//!
//! `history` is the `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape,
//! TopTools_ShapeMapHasher>` of `ChFi2d_Builder.hxx:280`; the port keeps the
//! insertion order in a `Vec` (the map is only ever searched, never walked by
//! position in a way that a caller can observe) and compares keys with
//! `TopoDS_Shape::IsSame`.

use crate::shape::{Edge, TopoShape};
use crate::topo_tools::is_same;

use super::ChFi2dBuilder;

impl ChFi2dBuilder {
    /// `history.Find(E)`: the modified edge bound to `E`, `None` when `E` is
    /// not bound (OCCT throws `Standard_NoSuchObject` there).
    pub(super) fn find_history(&self, e: &Edge) -> Option<&TopoShape> {
        self.history
            .iter()
            .find(|(k, _)| is_same(k, &e.0))
            .map(|(_, v)| v)
    }

    /// `history.Bind(E1, TrimE1)`, overwriting an existing binding the way
    /// `NCollection_DataMap::Bind` does.
    pub(super) fn bind_history(&mut self, basis: &TopoShape, modified: &TopoShape) {
        match self.history.iter_mut().find(|(k, _)| is_same(k, basis)) {
            Some(slot) => slot.1 = modified.clone(),
            None => self.history.push((basis.clone(), modified.clone())),
        }
    }

    /// `history.UnBind(E)`: a no-op when `E` is not bound.
    pub(super) fn unbind_history(&mut self, e: &TopoShape) {
        self.history.retain(|(k, _)| !is_same(k, e));
    }

    /// `BasisEdge(E)` (`cxx:747-762`): the history key whose value is `E`, or
    /// `E` itself when `E` is already a basis edge.
    ///
    /// OCCT walks the map with an iterator and returns the first key found;
    /// when several basis edges map to the same descendant the answer follows
    /// the map's iteration order, which the port pins to insertion order.
    pub fn basis_edge(&self, e: &Edge) -> Edge {
        for (k, v) in &self.history {
            if is_same(v, &e.0) {
                if let Some(edge) = Edge::wrap(k.clone()) {
                    return edge;
                }
            }
        }
        e.clone()
    }

    /// `UpDateHistory(E1, E2, TrimE1, TrimE2, NewEdge, Id)` (`cxx:684-717`):
    /// append the new edge to `fillets` (`Id == 1`) or `chamfers`, then bind
    /// each trimmed edge to its basis edge unless the corresponding edge was
    /// reported degenerated.
    pub(super) fn update_history(
        &mut self,
        e1: &Edge,
        e2: &Edge,
        trim_e1: &Edge,
        trim_e2: &Edge,
        new_edge: &Edge,
        id: i32,
    ) {
        // `cxx:690-698`.
        if id == 1 {
            self.fillets.push(new_edge.0.clone());
        } else {
            self.chamfers.push(new_edge.0.clone());
        }
        // `cxx:700-706`.
        self.unbind_history(&e1.0);
        if self.status != super::ConstructionError::FirstEdgeDegenerated
            && self.status != super::ConstructionError::BothEdgesDegenerated
        {
            if !is_same(&e1.0, &trim_e1.0) {
                self.bind_history(&e1.0, &trim_e1.0);
            }
        }
        // `cxx:708-716`.
        self.unbind_history(&e2.0);
        if self.status != super::ConstructionError::LastEdgeDegenerated
            && self.status != super::ConstructionError::BothEdgesDegenerated
        {
            if !is_same(&e2.0, &trim_e2.0) {
                self.bind_history(&e2.0, &trim_e2.0);
            }
        }
    }

    /// `UpDateHistory(E1, E2, TrimE1, TrimE2)` (`cxx:720-745`), the overload
    /// `RemoveFillet` / `RemoveChamfer` use: rebind both edges unconditionally.
    pub(super) fn replace_history(&mut self, e1: &Edge, e2: &Edge, trim_e1: &Edge, trim_e2: &Edge) {
        // `cxx:727-733`.
        self.unbind_history(&e1.0);
        if !is_same(&e1.0, &trim_e1.0) {
            self.bind_history(&e1.0, &trim_e1.0);
        }
        // `cxx:734-741`.
        self.unbind_history(&e2.0);
        if !is_same(&e2.0, &trim_e2.0) {
            self.bind_history(&e2.0, &trim_e2.0);
        }
    }
}
