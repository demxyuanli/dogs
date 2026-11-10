use std::collections::BTreeSet;

use super::super::delaun_types::{DelaunTriangle, DelaunVertex};
use super::structure::DelaunDataStructure;

/// Selector over the mesh accumulating nodes/links/elements adjacent to a seed.
///
/// Mirrors `BRepMesh_SelectorOfDataStructureOfDelaun`. In OCCT only the element
/// set is filled by the implemented queries; the node/link/frontier sets are
/// kept for API parity (they remain empty unless later callers populate them).
pub struct DelaunSelector<'a> {
    mesh: &'a DelaunDataStructure,
    nodes: BTreeSet<i32>,
    links: BTreeSet<i32>,
    elements: BTreeSet<i32>,
    frontier: BTreeSet<i32>,
}

impl<'a> DelaunSelector<'a> {
    /// Creates an empty selector bound to the given mesh.
    pub fn new(mesh: &'a DelaunDataStructure) -> Self {
        Self {
            mesh,
            nodes: BTreeSet::new(),
            links: BTreeSet::new(),
            elements: BTreeSet::new(),
            frontier: BTreeSet::new(),
        }
    }

    /// Resets the selection sets.
    pub fn initialize(&mut self, mesh: &'a DelaunDataStructure) {
        self.mesh = mesh;
        self.nodes.clear();
        self.links.clear();
        self.elements.clear();
        self.frontier.clear();
    }

    /// Selects all neighboring elements of the given node.
    pub fn neighbors_of(&mut self, node: &DelaunVertex) {
        let idx = self.mesh.index_of_node(node);
        if idx != 0 {
            self.neighbors_of_node(idx);
        }
    }

    /// Selects all neighboring elements of the node with the given index.
    pub fn neighbors_of_node(&mut self, node_index: i32) {
        let link_ids: Vec<i32> = self.mesh.links_connected_to(node_index).to_vec();
        self.nodes.insert(node_index);
        for id in link_ids {
            self.elements_of_link(id);
        }
    }

    /// Selects all neighboring elements of the given link.
    pub fn neighbors_of_link(&mut self, link_index: i32) {
        let link = self.mesh.link(link_index);
        self.neighbors_of_node(link.first_node());
        self.neighbors_of_node(link.last_node());
    }

    /// Selects all neighboring elements of the element with the given index.
    pub fn neighbors_of_element(&mut self, element_index: i32) {
        let element = self.mesh.element(element_index);
        for v in element.vertex_indices {
            self.neighbors_of_node(v);
        }
    }

    /// Selects all neighboring elements of the given triangle by its nodes.
    pub fn neighbors_of_triangle(&mut self, element: &DelaunTriangle) {
        for v in element.vertex_indices {
            self.neighbors_of_node(v);
        }
    }

    /// Selects all neighboring elements by the links of the given triangle.
    pub fn neighbors_by_edge_of(&mut self, element: &DelaunTriangle) {
        for e in element.link_indices {
            self.elements_of_link(e.abs());
        }
    }

    /// Adds a level of neighbors by edge to the selector (OCCT no-op).
    pub fn add_neighbors(&mut self) {}

    /// Selected node indices.
    pub fn nodes(&self) -> &BTreeSet<i32> {
        &self.nodes
    }

    /// Selected link indices.
    pub fn links(&self) -> &BTreeSet<i32> {
        &self.links
    }

    /// Selected element indices.
    pub fn elements(&self) -> &BTreeSet<i32> {
        &self.elements
    }

    /// Frontier link indices (currently unpopulated, OCCT parity).
    pub fn frontier_links(&self) -> &BTreeSet<i32> {
        &self.frontier
    }

    fn elements_of_link(&mut self, index: i32) {
        let pair = self.mesh.elements_connected_to(index);
        self.links.insert(index);
        for j in 1..=pair.extent() {
            self.elements.insert(pair.index(j));
        }
    }
}
