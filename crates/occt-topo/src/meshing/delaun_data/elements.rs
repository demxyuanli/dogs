use std::collections::BTreeSet;

use super::super::delaun_types::{DelaunPairOfIndex, DelaunTriangle, VertexState};
use super::structure::DelaunDataStructure;

impl DelaunDataStructure {
    /// Number of element slots (including deleted ones), as OCCT `NbElements`.
    pub fn nb_elements(&self) -> usize {
        self.elements.len()
    }

    /// Adds a triangle (as `Free`), registering it with each of its links.
    pub fn add_element(&mut self, element: DelaunTriangle) -> i32 {
        let idx = (self.elements.len() + 1) as i32;
        // TEMPORARY diagnostic (DELADDT): logs every element creation so it can be
        // diffed against OCCT's `OCCTEL` dump of the same structure.
        if crate::debug_flags::del_addt() {
            let n = self.element_nodes(&element);
            eprintln!("ADDE {idx} {} {} {}", n[0], n[1], n[2]);
        }
        self.elements.push(element);
        self.element_states.push(VertexState::Free);
        self.elements_of_domain.insert(idx);
        for e in element.link_indices {
            self.link_elements[(e.abs() - 1) as usize].append(idx);
        }
        idx
    }

    /// Triangle by 1-based index.
    pub fn element(&self, index: i32) -> DelaunTriangle {
        self.elements[(index - 1) as usize]
    }

    /// Movability of the triangle with the given index.
    pub fn element_movability(&self, index: i32) -> VertexState {
        self.element_states[(index - 1) as usize]
    }

    /// Sets the movability of a triangle.
    pub fn set_element_movability(&mut self, index: i32, state: VertexState) {
        self.element_states[(index - 1) as usize] = state;
    }

    /// Substitutes a triangle. When the old triangle is alive its edge-element
    /// registrations are cleaned first.
    pub fn substitute_element(&mut self, index: i32, new_element: DelaunTriangle) -> bool {
        let idx = (index - 1) as usize;
        if self.element_states[idx] == VertexState::Deleted {
            self.elements[idx] = new_element;
            return true;
        }
        let old = self.elements[idx];
        self.clean_element(index, &old);
        self.elements[idx] = new_element;
        for e in new_element.link_indices {
            self.link_elements[(e.abs() - 1) as usize].append(index);
        }
        true
    }

    /// Removes a triangle, deregistering it from its links.
    pub fn remove_element(&mut self, index: i32) {
        let idx = (index - 1) as usize;
        if self.element_states[idx] == VertexState::Deleted {
            return;
        }
        let element = self.elements[idx];
        self.clean_element(index, &element);
        self.element_states[idx] = VertexState::Deleted;
        self.elements_of_domain.remove(&index);
    }

    pub(super) fn clean_element(&mut self, index: i32, element: &DelaunTriangle) {
        if self.element_states[(index - 1) as usize] != VertexState::Free {
            return;
        }
        for e in element.link_indices {
            let pair = &mut self.link_elements[(e.abs() - 1) as usize];
            remove_element_index(index, pair);
        }
    }

    /// Live element ids registered in the mesh. Source: `ElementsOfDomain`.
    /// Ascending order, matching the OCCT packed-map traversal (see the field
    /// comment on `elements_of_domain`).
    pub fn elements_of_domain(&self) -> &BTreeSet<i32> {
        &self.elements_of_domain
    }

    /// Vertex indices of a triangle via `ElementNodes`
    /// (`BRepMesh_DataStructureOfDelaun.cxx:259-286`): reconstruct from link 0
    /// (both ends) and link 2 (the remaining corner). Cached `vertex_indices`
    /// are not the source of truth.
    pub fn element_nodes(&self, triangle: &DelaunTriangle) -> [i32; 3] {
        let mut nodes = [0i32; 3];
        let link1 = self.link(triangle.link_at(0).abs());
        if triangle.link_at(0) > 0 {
            nodes[0] = link1.first_node();
            nodes[1] = link1.last_node();
        } else {
            nodes[1] = link1.first_node();
            nodes[0] = link1.last_node();
        }
        let link2 = self.link(triangle.link_at(2).abs());
        nodes[2] = if triangle.link_at(2) > 0 {
            link2.first_node()
        } else {
            link2.last_node()
        };
        nodes
    }
}

pub(super) fn remove_element_index(index: i32, pair: &mut DelaunPairOfIndex) {
    let n = pair.extent();
    for i in 1..=n {
        if pair.index(i) == index {
            pair.remove_index(i);
            return;
        }
    }
}
