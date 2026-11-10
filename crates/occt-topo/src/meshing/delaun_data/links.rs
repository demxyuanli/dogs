use std::collections::BTreeSet;

use super::super::delaun_types::{DelaunLink, DelaunPairOfIndex, VertexState};
use super::structure::DelaunDataStructure;

impl DelaunDataStructure {
    /// Number of link slots (including deleted ones), as OCCT `NbLinks`.
    pub fn nb_links(&self) -> usize {
        self.links.len()
    }

    /// Adds a link between two nodes. Returns a positive index when the link is
    /// newly created or already present with the same orientation, and a
    /// *negative* index when an equivalent link exists in the opposite
    /// orientation. Source: `AddLink`.
    ///
    /// A recycled index is taken from the **front** of the deleted list
    /// (`BRepMesh_DataStructureOfDelaun.cxx:77-89`: `myDelLinks.First()` then
    /// `myDelLinks.RemoveFirst()`; `myDelLinks` is an `NCollection_List` fed by
    /// `Append` at `:138-143`, i.e. a FIFO). Reusing the most recently freed
    /// index instead makes a fresh link collide with a link id that a polygon
    /// captured earlier in the same `meshPolygon` and was just freed by
    /// `cleanupPolygon`.
    pub fn add_link(&mut self, first: i32, last: i32, state: VertexState) -> i32 {
        let probe = DelaunLink::new_unparameterized(first, last, 0);
        let existing = self.index_of_link(&probe);
        if existing > 0 {
            let l = self.link(existing);
            return if l.is_same_orientation(&probe) { existing } else { -existing };
        }

        let id;
        if let Some(del) = self.del_links.pop_front() {
            id = del;
            let old = self.links[(id - 1) as usize];
            let nl = DelaunLink::new_unparameterized(first, last, id);
            self.idx_set(id, Some(old), &nl);
            self.links[(id - 1) as usize] = nl;
            self.link_states[(id - 1) as usize] = state;
            self.link_elements[(id - 1) as usize] = DelaunPairOfIndex::default();
        } else {
            id = (self.links.len() + 1) as i32;
            let nl = DelaunLink::new_unparameterized(first, last, id);
            self.idx_set(id, None, &nl);
            self.links.push(nl);
            self.link_states.push(state);
            self.link_elements.push(DelaunPairOfIndex::default());
        }
        self.links_connected_to_mut(first).push(id);
        self.links_connected_to_mut(last).push(id);
        self.links_of_domain.insert(id);
        id
    }

    /// Finds the index of an equivalent link (either orientation), `0` if none.
    pub fn index_of_link(&self, link: &DelaunLink) -> i32 {
        let k = Self::lkey(link);
        let mut best = 0;
        if let Some(ids) = self.link_index.get(&k) {
            for &id in ids {
                let i = (id - 1) as usize;
                if self.link_states[i] == VertexState::Deleted || !self.links[i].is_equal(link) { continue; }
                if best == 0 || id < best { best = id; }
            }
        }
        best
    }

    pub(super) fn lkey(l: &DelaunLink) -> (i32, i32) {
        let (a, b) = (l.first_node(), l.last_node());
        if a <= b { (a, b) } else { (b, a) }
    }

    pub(super) fn idx_set(&mut self, id: i32, old: Option<DelaunLink>, new: &DelaunLink) {
        if let Some(o) = old {
            if let Some(v) = self.link_index.get_mut(&Self::lkey(&o)) { v.retain(|&x| x != id); }
        }
        self.link_index.entry(Self::lkey(new)).or_default().push(id);
    }

    /// Link by 1-based index.
    pub fn link(&self, index: i32) -> DelaunLink {
        self.links[(index - 1) as usize]
    }

    /// Movability of the link with the given index.
    pub fn link_movability(&self, index: i32) -> VertexState {
        self.link_states[(index - 1) as usize]
    }

    /// Sets the movability of a link.
    pub fn set_link_movability(&mut self, index: i32, state: VertexState) {
        self.link_states[(index - 1) as usize] = state;
    }

    /// Substitutes a link. Returns `false` when the replacement already exists
    /// and the current link is not deleted.
    pub fn substitute_link(&mut self, index: i32, new_link: DelaunLink, state: VertexState) -> bool {
        let idx = (index - 1) as usize;
        if self.link_states[idx] == VertexState::Deleted {
            let old = self.links[idx];
            self.idx_set(index, Some(old), &new_link);
            self.links[idx] = new_link;
            self.link_states[idx] = state;
            self.link_elements[idx] = DelaunPairOfIndex::default();
            return true;
        }
        if self.index_of_link(&new_link) != 0 {
            return false;
        }
        let old = self.links[idx];
        self.link_states[idx] = VertexState::Deleted;
        self.clean_link(index, &old);
        self.links_connected_to_mut(new_link.first_node()).push(index);
        self.links_connected_to_mut(new_link.last_node()).push(index);
        self.idx_set(index, Some(old), &new_link);
        self.links[idx] = new_link;
        self.link_states[idx] = state;
        self.link_elements[idx] = DelaunPairOfIndex::default();
        true
    }

    /// Removes a link when it has no connected elements and (unless `is_force`)
    /// is `Free`. Source: `RemoveLink`.
    pub fn remove_link(&mut self, index: i32, is_force: bool) {
        let idx = (index.abs() - 1) as usize;
        if self.link_states[idx] == VertexState::Deleted {
            return;
        }
        if !is_force && self.link_states[idx] != VertexState::Free {
            return;
        }
        if self.link_elements[idx].extent() != 0 {
            return;
        }
        let link = self.links[idx];
        self.clean_link(index, &link);
        self.link_states[idx] = VertexState::Deleted;
        self.links_of_domain.remove(&index);
        self.del_links.push_back(index);
    }

    pub(super) fn clean_link(&mut self, index: i32, link: &DelaunLink) {
        for i in 0..2 {
            let node_id = if i == 0 { link.first_node() } else { link.last_node() };
            if let Some(list) = self.node_links.get_mut((node_id - 1) as usize) {
                if let Some(pos) = list.iter().position(|&x| x == index) {
                    list.remove(pos);
                }
            }
        }
    }

    /// Elements connected to the link with the given index (max two).
    pub fn elements_connected_to(&self, link_index: i32) -> DelaunPairOfIndex {
        self.link_elements[(link_index.abs() - 1) as usize]
    }

    /// Live link ids registered in the mesh. Source: `LinksOfDomain`.
    pub fn links_of_domain(&self) -> &BTreeSet<i32> {
        &self.links_of_domain
    }
}
