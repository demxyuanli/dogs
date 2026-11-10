use super::super::delaun_types::{DelaunVertex, VertexState};
use super::structure::DelaunDataStructure;

impl DelaunDataStructure {
    /// Number of node slots (including deleted ones), as OCCT `NbNodes`.
    pub fn nb_nodes(&self) -> usize {
        self.vertices.len()
    }

    /// Adds a node, returning its (1-based) index. Coincident nodes collapse to
    /// the existing index unless `is_force_add` is set.
    pub fn add_node(&mut self, node: DelaunVertex) -> i32 {
        let idx = self.add_vertex(node, false);
        self.ensure_node_links(idx);
        idx
    }

    /// Force-adds a node without coincidence checking.
    pub fn add_node_force(&mut self, node: DelaunVertex) -> i32 {
        let idx = self.add_vertex(node, true);
        self.ensure_node_links(idx);
        idx
    }

    fn add_vertex(&mut self, node: DelaunVertex, is_force_add: bool) -> i32 {
        if !is_force_add {
            let existing = self.vertex_cells.find_index(node.location, &self.vertices);
            if existing != 0 {
                return existing;
            }
        }
        // Deleted node slots are recycled from the front of the list, matching
        // `BRepMesh_VertexInspector.hxx:54-64` (`myDelNodes.First()` +
        // `RemoveFirst()`, fed by `Append` at `:95-99` -- a FIFO).
        let idx = if let Some(del) = self.del_vertices.pop_front() {
            self.vertices[(del - 1) as usize] = node;
            del
        } else {
            self.vertices.push(node);
            self.vertices.len() as i32
        };
        self.vertex_cells.add(idx, node.location);
        idx
    }

    fn ensure_node_links(&mut self, idx: i32) {
        if self.node_links.len() < idx as usize {
            self.node_links.resize(idx as usize, Vec::new());
        }
    }

    /// Finds the index of a coincident node, `0` if none. Source: `IndexOf`.
    pub fn index_of_node(&self, node: &DelaunVertex) -> i32 {
        self.vertex_cells.find_index(node.location, &self.vertices)
    }

    /// Node by 1-based index.
    pub fn node(&self, index: i32) -> &DelaunVertex {
        &self.vertices[(index - 1) as usize]
    }

    /// Mutable node by 1-based index.
    pub fn node_mut(&mut self, index: i32) -> &mut DelaunVertex {
        &mut self.vertices[(index - 1) as usize]
    }

    /// Substitutes a node; returns `false` when a coincident node already
    /// exists (in which case nothing is changed).
    pub fn substitute_node(&mut self, index: i32, new_node: DelaunVertex) -> bool {
        if self.vertex_cells.find_index(new_node.location, &self.vertices) != 0 {
            return false;
        }
        let old_p = self.vertices[(index - 1) as usize].location;
        self.vertex_cells.remove(index, old_p);
        self.vertices[(index - 1) as usize] = new_node;
        self.vertex_cells.add(index, new_node.location);
        true
    }

    /// Removes a node when it is `Free` (or `is_force`) and has no connected
    /// links. Source: `RemoveNode`.
    pub fn remove_node(&mut self, index: i32, is_force: bool) {
        let is_free = self.vertices[(index - 1) as usize].state == VertexState::Free;
        if (is_force || is_free) && self.links_connected_to(index).is_empty() {
            let p = self.vertices[(index - 1) as usize].location;
            self.vertices[(index - 1) as usize].state = VertexState::Deleted;
            self.vertex_cells.remove(index, p);
            self.del_vertices.push_back(index);
        }
    }

    /// List of link ids attached to the node with the given index.
    pub fn links_connected_to(&self, index: i32) -> &[i32] {
        if index <= 0 {
            return &[];
        }
        match self.node_links.get((index - 1) as usize) {
            Some(list) => list,
            None => &[],
        }
    }

    pub(super) fn links_connected_to_mut(&mut self, index: i32) -> &mut Vec<i32> {
        self.ensure_node_links(index);
        &mut self.node_links[(index - 1) as usize]
    }
}
