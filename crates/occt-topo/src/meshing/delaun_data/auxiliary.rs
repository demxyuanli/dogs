use std::collections::{BTreeSet, HashMap};

use super::super::delaun_types::{DelaunLink, DelaunPairOfIndex, DelaunVertex, VertexState};
use super::structure::DelaunDataStructure;
use super::vertex_cells::VertexCellFilter;

impl DelaunDataStructure {
    /// Coincidence tolerance used by the vertex filter (OCCT
    /// `VertexTool::GetTolerance`).
    pub fn tolerance(&self) -> (f64, f64) {
        self.vertex_cells.tolerance
    }

    /// Removes all elements (and links left free by that), as `ClearDomain`.
    pub fn clear_domain(&mut self) {
        let mut free_edges: BTreeSet<i32> = BTreeSet::new();
        let element_ids: Vec<i32> = self.elements_of_domain.iter().copied().collect();
        for id in element_ids {
            let element = self.element(id);
            for e in element.link_indices {
                free_edges.insert(e.abs());
            }
            self.clean_element(id, &element);
            self.element_states[(id - 1) as usize] = VertexState::Deleted;
        }
        self.elements_of_domain.clear();
        for e in free_edges {
            self.remove_link(e, false);
        }
    }

    /// Compacts deleted links/nodes so only live items remain. Source:
    /// `ClearDeleted` (used by the mesh healer, not by the Delaun insertion).
    pub fn clear_deleted(&mut self) {
        self.clear_deleted_links();
        self.clear_deleted_nodes();
    }

    fn clear_deleted_links(&mut self) {
        let n = self.links.len();
        let mut remap: Vec<i32> = vec![0; n + 1];
        let mut new_links: Vec<DelaunLink> = Vec::with_capacity(n);
        let mut new_states: Vec<VertexState> = Vec::with_capacity(n);
        let mut new_pairs: Vec<DelaunPairOfIndex> = Vec::with_capacity(n);
        for i in 0..n {
            let old_id = (i + 1) as i32;
            if self.link_states[i] == VertexState::Deleted {
                continue;
            }
            let new_id = (new_links.len() + 1) as i32;
            remap[old_id as usize] = new_id;
            let mut l = self.links[i];
            l.index = new_id;
            new_links.push(l);
            new_states.push(self.link_states[i]);
            new_pairs.push(self.link_elements[i]);
        }
        self.links = new_links;
        self.link_states = new_states;
        self.link_elements = new_pairs;
        for list in self.node_links.iter_mut() {
            for id in list.iter_mut() {
                *id = remap[*id as usize];
            }
        }
        for el in self.elements.iter_mut() {
            for e in el.link_indices.iter_mut() {
                let sign = if *e < 0 { -1 } else { 1 };
                *e = sign * remap[e.abs() as usize];
            }
        }
        self.links_of_domain = (1..=self.links.len() as i32).collect();
        self.del_links.clear();
        self.link_index.clear();
        for i in 0..self.links.len() { let l = self.links[i]; self.link_index.entry(Self::lkey(&l)).or_default().push((i + 1) as i32); }
    }

    fn clear_deleted_nodes(&mut self) {
        let n = self.vertices.len();
        let mut remap: Vec<i32> = vec![0; n + 1];
        let mut new_verts: Vec<DelaunVertex> = Vec::with_capacity(n);
        let mut new_node_links: Vec<Vec<i32>> = Vec::with_capacity(n);
        for i in 0..n {
            let old_id = (i + 1) as i32;
            if self.vertices[i].state == VertexState::Deleted {
                continue;
            }
            let new_id = (new_verts.len() + 1) as i32;
            remap[old_id as usize] = new_id;
            new_verts.push(self.vertices[i]);
            new_node_links.push(self.node_links[i].clone());
        }
        for l in self.links.iter_mut() {
            l.v1 = remap[l.v1 as usize];
            l.v2 = remap[l.v2 as usize];
        }
        self.link_index.clear();
        for i in 0..self.links.len() { let l = self.links[i]; self.link_index.entry(Self::lkey(&l)).or_default().push((i + 1) as i32); }
        self.vertices = new_verts;
        self.node_links = new_node_links;
        self.vertex_cells = VertexCellFilter {
            cell_size: self.vertex_cells.cell_size,
            tolerance: self.vertex_cells.tolerance,
            cells: HashMap::new(),
        };
        for (i, v) in self.vertices.iter().enumerate() {
            self.vertex_cells.add((i + 1) as i32, v.location);
        }
        self.del_vertices.clear();
    }
}
