//! Port of OCCT `BRepMesh_DataStructureOfDelaun` + `BRepMesh_SelectorOfDataStructureOfDelaun`.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_DataStructureOfDelaun.{hxx,cxx}` — indexed storage for vertices,
//!   links (edges) and triangles plus a vertex-cell filter for coincidence
//!   queries (the OCCT class holds this as a `BRepMesh_VertexTool`).
//! - `BRepMesh_SelectorOfDataStructureOfDelaun.{hxx,cxx}` — neighbour selector.
//!
//! ponytail: the OCCT value types carry movability inside the struct
//! (`BRepMesh_Vertex/Edge/Triangle`). The Rust contract types in `delaun_types`
//! (`DelaunVertex`/`DelaunLink`/`DelaunTriangle`) do not, so the data structure
//! keeps parallel per-element state arrays here. When the contract is extended,
//! these arrays can be folded into the value types.
//!
//! Indices follow OCCT: vertices/links/elements are **1-based** in the public
//! API; storage is 0-based `Vec`s.

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpPnt2d, GpXY};
use occt_core::precision::{CONFUSION, PCONFUSION};

use super::delaun_types::{DelaunLink, DelaunPairOfIndex, DelaunTriangle, DelaunVertex, VertexState};

/// Uniform-grid cell filter over vertices used to collapse coincident UV points.
///
/// Source: `BRepMesh_VertexTool` + `BRepMesh_VertexInspector` (the OCCT class
/// stores a `NCollection_CellFilter`); this is a `HashMap`-backed equivalent
/// with the same cell size / tolerance semantics.
struct VertexCellFilter {
    /// `(U, V)` cell sizes. Source: `BRepMesh_VertexTool::SetCellSize`.
    cell_size: (f64, f64),
    tolerance: (f64, f64),
    cells: HashMap<(i64, i64), Vec<i32>>,
}

impl VertexCellFilter {
    fn new(cell_size: f64, tolerance: (f64, f64)) -> Self {
        Self {
            cell_size: (cell_size, cell_size),
            tolerance,
            cells: HashMap::new(),
        }
    }

    fn cell_of(&self, p: GpXY) -> (i64, i64) {
        let cu = if self.cell_size.0 > 0.0 { self.cell_size.0 } else { 1.0 };
        let cv = if self.cell_size.1 > 0.0 { self.cell_size.1 } else { 1.0 };
        ((p.x / cu).floor() as i64, (p.y / cv).floor() as i64)
    }

    fn add(&mut self, index: i32, p: GpPnt2d) {
        self.cells.entry(self.cell_of(p.coord)).or_default().push(index);
    }

    fn remove(&mut self, index: i32, p: GpPnt2d) {
        let key = self.cell_of(p.coord);
        if let Some(list) = self.cells.get_mut(&key) {
            if let Some(pos) = list.iter().position(|&x| x == index) {
                list.remove(pos);
            }
        }
    }

    /// Returns the index of the closest live vertex within the coincidence
    /// tolerance of `p`, or `0` when none exists. Source:
    /// `BRepMesh_VertexInspector::Inspect` + `GetCoincidentPoint`.
    fn find_index(&self, p: GpPnt2d, vertices: &[DelaunVertex]) -> i32 {
        let (tx, ty) = self.tolerance;
        let min_cell = self.cell_of(GpXY::new(p.x() - tx, p.y() - ty));
        let max_cell = self.cell_of(GpXY::new(p.x() + tx, p.y() + ty));
        let sq_tx = tx * tx;
        let sq_ty = ty * ty;
        let mut best = 0;
        let mut best_sq = f64::INFINITY;
        for ci in min_cell.0..=max_cell.0 {
            for cj in min_cell.1..=max_cell.1 {
                if let Some(list) = self.cells.get(&(ci, cj)) {
                    for &idx in list {
                        let v = &vertices[(idx - 1) as usize];
                        if v.state == VertexState::Deleted {
                            continue;
                        }
                        let dx = v.location.x() - p.x();
                        let dy = v.location.y() - p.y();
                        let sq = dx * dx + dy * dy;
                        // 2D tolerance (both axes) set in the default tool => the
                        // Euclidean branch of `Inspect` applies.
                        let in_tol = if sq_ty.abs() < CONFUSION {
                            sq < sq_tx
                        } else {
                            dx * dx < sq_tx && dy * dy < sq_ty
                        };
                        if in_tol && sq < best_sq {
                            best_sq = sq;
                            best = idx;
                        }
                    }
                }
            }
        }
        best
    }
}

/// Indexed storage for the Delaunay core: vertices, links and triangles.
///
/// Mirrors `BRepMesh_DataStructureOfDelaun`. Vertex uniqueness is provided by
/// the internal cell filter; links deduplicate by their two endpoint indices;
/// triangles are plain slots. Deleted items keep their slots (so element/link
/// ids are stable during the Bowyer–Watson insertion) and are reused or
/// compacted via the deleted lists / `clear_deleted`.
pub struct DelaunDataStructure {
    vertices: Vec<DelaunVertex>,
    del_vertices: Vec<i32>,
    vertex_cells: VertexCellFilter,
    node_links: Vec<Vec<i32>>,

    links: Vec<DelaunLink>,
    link_states: Vec<VertexState>,
    link_elements: Vec<DelaunPairOfIndex>,
    del_links: Vec<i32>,
    links_of_domain: HashSet<i32>,

    elements: Vec<DelaunTriangle>,
    element_states: Vec<VertexState>,
    elements_of_domain: HashSet<i32>,
}

impl DelaunDataStructure {
    /// Sets the 2D cell-filter pitch. Source: `BRepMesh_VertexTool::SetCellSize`.
    pub fn set_cell_size(&mut self, size_u: f64, size_v: f64) {
        self.vertex_cells.cell_size = (size_u.max(1e-16), size_v.max(1e-16));
    }

    /// Sets the coincidence tolerance. Source: `BRepMesh_VertexTool::SetTolerance`.
    pub fn set_tolerance(&mut self, tol_u: f64, tol_v: f64) {
        self.vertex_cells.tolerance = (tol_u, tol_v);
    }

    /// Creates an empty data structure. `reserved_node_size` is only a capacity
    /// hint, matching the OCCT constructor.
    pub fn new(reserved_node_size: usize) -> Self {
        let cap = reserved_node_size.max(16);
        let cell_size = CONFUSION + 0.05 * CONFUSION;
        Self {
            vertices: Vec::with_capacity(cap),
            del_vertices: Vec::new(),
            vertex_cells: VertexCellFilter::new(cell_size, (CONFUSION, CONFUSION)),
            node_links: Vec::with_capacity(cap),
            links: Vec::with_capacity(cap * 3),
            link_states: Vec::with_capacity(cap * 3),
            link_elements: Vec::with_capacity(cap * 3),
            del_links: Vec::new(),
            links_of_domain: HashSet::new(),
            elements: Vec::with_capacity(cap * 2),
            element_states: Vec::with_capacity(cap * 2),
            elements_of_domain: HashSet::new(),
        }
    }

    // ------------------------------------------------------------------
    // Nodes
    // ------------------------------------------------------------------

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
        let idx = if let Some(del) = self.del_vertices.pop() {
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
    pub fn get_node(&self, index: i32) -> &DelaunVertex {
        &self.vertices[(index - 1) as usize]
    }

    /// Mutable node by 1-based index.
    pub fn get_node_mut(&mut self, index: i32) -> &mut DelaunVertex {
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
            self.del_vertices.push(index);
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

    fn links_connected_to_mut(&mut self, index: i32) -> &mut Vec<i32> {
        self.ensure_node_links(index);
        &mut self.node_links[(index - 1) as usize]
    }

    // ------------------------------------------------------------------
    // Links
    // ------------------------------------------------------------------

    /// Number of link slots (including deleted ones), as OCCT `NbLinks`.
    pub fn nb_links(&self) -> usize {
        self.links.len()
    }

    /// Adds a link between two nodes. Returns a positive index when the link is
    /// newly created or already present with the same orientation, and a
    /// *negative* index when an equivalent link exists in the opposite
    /// orientation. Source: `AddLink`.
    pub fn add_link(&mut self, first: i32, last: i32, state: VertexState) -> i32 {
        let probe = DelaunLink::new_unparameterized(first, last, 0);
        let existing = self.index_of_link(&probe);
        if existing > 0 {
            let l = self.get_link(existing);
            return if l.is_same_orientation(&probe) { existing } else { -existing };
        }

        let id;
        if let Some(del) = self.del_links.pop() {
            id = del;
            self.links[(id - 1) as usize] = DelaunLink::new_unparameterized(first, last, id);
            self.link_states[(id - 1) as usize] = state;
            self.link_elements[(id - 1) as usize] = DelaunPairOfIndex::default();
        } else {
            id = (self.links.len() + 1) as i32;
            self.links.push(DelaunLink::new_unparameterized(first, last, id));
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
        for (i, l) in self.links.iter().enumerate() {
            if self.link_states[i] == VertexState::Deleted {
                continue;
            }
            if l.is_equal(link) {
                return (i + 1) as i32;
            }
        }
        0
    }

    /// Link by 1-based index.
    pub fn get_link(&self, index: i32) -> DelaunLink {
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
        self.del_links.push(index);
    }

    fn clean_link(&mut self, index: i32, link: &DelaunLink) {
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
    pub fn links_of_domain(&self) -> &HashSet<i32> {
        &self.links_of_domain
    }

    // ------------------------------------------------------------------
    // Elements (triangles)
    // ------------------------------------------------------------------

    /// Number of element slots (including deleted ones), as OCCT `NbElements`.
    pub fn nb_elements(&self) -> usize {
        self.elements.len()
    }

    /// Adds a triangle (as `Free`), registering it with each of its links.
    pub fn add_element(&mut self, element: DelaunTriangle) -> i32 {
        let idx = (self.elements.len() + 1) as i32;
        self.elements.push(element);
        self.element_states.push(VertexState::Free);
        self.elements_of_domain.insert(idx);
        for e in element.link_indices {
            self.link_elements[(e.abs() - 1) as usize].append(idx);
        }
        idx
    }

    /// Triangle by 1-based index.
    pub fn get_element(&self, index: i32) -> DelaunTriangle {
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

    fn clean_element(&mut self, index: i32, element: &DelaunTriangle) {
        if self.element_states[(index - 1) as usize] != VertexState::Free {
            return;
        }
        for e in element.link_indices {
            let pair = &mut self.link_elements[(e.abs() - 1) as usize];
            remove_element_index(index, pair);
        }
    }

    /// Live element ids registered in the mesh. Source: `ElementsOfDomain`.
    pub fn elements_of_domain(&self) -> &HashSet<i32> {
        &self.elements_of_domain
    }

    /// Vertex indices of a triangle via `ElementNodes`
    /// (`BRepMesh_DataStructureOfDelaun.cxx:259-286`): reconstruct from link 0
    /// (both ends) and link 2 (the remaining corner). Cached `vertex_indices`
    /// are not the source of truth.
    pub fn element_nodes(&self, triangle: &DelaunTriangle) -> [i32; 3] {
        let mut nodes = [0i32; 3];
        let link1 = self.get_link(triangle.link_at(0).abs());
        if triangle.link_at(0) > 0 {
            nodes[0] = link1.first_node();
            nodes[1] = link1.last_node();
        } else {
            nodes[1] = link1.first_node();
            nodes[0] = link1.last_node();
        }
        let link2 = self.get_link(triangle.link_at(2).abs());
        nodes[2] = if triangle.link_at(2) > 0 {
            link2.first_node()
        } else {
            link2.last_node()
        };
        nodes
    }

    // ------------------------------------------------------------------
    // Auxiliary
    // ------------------------------------------------------------------

    /// Coincidence tolerance used by the vertex filter (OCCT
    /// `VertexTool::GetTolerance`).
    pub fn get_tolerance(&self) -> (f64, f64) {
        self.vertex_cells.tolerance
    }

    /// Removes all elements (and links left free by that), as `ClearDomain`.
    pub fn clear_domain(&mut self) {
        let mut free_edges: HashSet<i32> = HashSet::new();
        let element_ids: Vec<i32> = self.elements_of_domain.iter().copied().collect();
        for id in element_ids {
            let element = self.get_element(id);
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

fn remove_element_index(index: i32, pair: &mut DelaunPairOfIndex) {
    let n = pair.extent();
    for i in 1..=n {
        if pair.index(i) == index {
            pair.remove_index(i);
            return;
        }
    }
}

/// Selector over the mesh accumulating nodes/links/elements adjacent to a seed.
///
/// Mirrors `BRepMesh_SelectorOfDataStructureOfDelaun`. In OCCT only the element
/// set is filled by the implemented queries; the node/link/frontier sets are
/// kept for API parity (they remain empty unless later callers populate them).
pub struct DelaunSelector<'a> {
    mesh: &'a DelaunDataStructure,
    nodes: HashSet<i32>,
    links: HashSet<i32>,
    elements: HashSet<i32>,
    frontier: HashSet<i32>,
}

impl<'a> DelaunSelector<'a> {
    /// Creates an empty selector bound to the given mesh.
    pub fn new(mesh: &'a DelaunDataStructure) -> Self {
        Self {
            mesh,
            nodes: HashSet::new(),
            links: HashSet::new(),
            elements: HashSet::new(),
            frontier: HashSet::new(),
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

    /// Selects all neighbouring elements of the given node.
    pub fn neighbours_of(&mut self, node: &DelaunVertex) {
        let idx = self.mesh.index_of_node(node);
        if idx != 0 {
            self.neighbours_of_node(idx);
        }
    }

    /// Selects all neighbouring elements of the node with the given index.
    pub fn neighbours_of_node(&mut self, node_index: i32) {
        let link_ids: Vec<i32> = self.mesh.links_connected_to(node_index).to_vec();
        self.nodes.insert(node_index);
        for id in link_ids {
            self.elements_of_link(id);
        }
    }

    /// Selects all neighbouring elements of the given link.
    pub fn neighbours_of_link(&mut self, link_index: i32) {
        let link = self.mesh.get_link(link_index);
        self.neighbours_of_node(link.first_node());
        self.neighbours_of_node(link.last_node());
    }

    /// Selects all neighbouring elements of the element with the given index.
    pub fn neighbours_of_element(&mut self, element_index: i32) {
        let element = self.mesh.get_element(element_index);
        for v in element.vertex_indices {
            self.neighbours_of_node(v);
        }
    }

    /// Selects all neighbouring elements of the given triangle by its nodes.
    pub fn neighbours_of_triangle(&mut self, element: &DelaunTriangle) {
        for v in element.vertex_indices {
            self.neighbours_of_node(v);
        }
    }

    /// Selects all neighbouring elements by the links of the given triangle.
    pub fn neighbours_by_edge_of(&mut self, element: &DelaunTriangle) {
        for e in element.link_indices {
            self.elements_of_link(e.abs());
        }
    }

    /// Adds a level of neighbours by edge to the selector (OCCT no-op).
    pub fn add_neighbours(&mut self) {}

    /// Selected node indices.
    pub fn nodes(&self) -> &HashSet<i32> {
        &self.nodes
    }

    /// Selected link indices.
    pub fn links(&self) -> &HashSet<i32> {
        &self.links
    }

    /// Selected element indices.
    pub fn elements(&self) -> &HashSet<i32> {
        &self.elements
    }

    /// Frontier link indices (currently unpopulated, OCCT parity).
    pub fn frontier_links(&self) -> &HashSet<i32> {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn v(u: f64, w: f64) -> DelaunVertex {
        DelaunVertex::new_parametric(u, w, VertexState::Free)
    }

    #[test]
    fn add_node_dedupes_coincident() {
        let mut ds = DelaunDataStructure::new(10);
        let a = ds.add_node(v(1.0, 2.0));
        let b = ds.add_node(v(1.0 + 1e-8, 2.0)); // within Confusion (1e-7)
        let c = ds.add_node(v(5.0, 5.0));
        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_eq!(ds.nb_nodes(), 2);
    }

    #[test]
    fn add_link_detects_reverse_orientation() {
        let mut ds = DelaunDataStructure::new(10);
        let a = ds.add_node(v(0.0, 0.0));
        let b = ds.add_node(v(1.0, 0.0));
        let l1 = ds.add_link(a, b, VertexState::Free);
        assert!(l1 > 0);
        let l2 = ds.add_link(b, a, VertexState::Free); // same link, reversed
        assert_eq!(l1, -l2);
        assert_eq!(ds.links_connected_to(a), &[l1]);
        assert_eq!(ds.links_connected_to(b), &[l1]);
    }

    #[test]
    fn element_link_registration_and_removal() {
        let mut ds = DelaunDataStructure::new(10);
        let a = ds.add_node(v(0.0, 0.0));
        let b = ds.add_node(v(1.0, 0.0));
        let c = ds.add_node(v(0.0, 1.0));
        let ab = ds.add_link(a, b, VertexState::Free);
        let bc = ds.add_link(b, c, VertexState::Free);
        let ca = ds.add_link(c, a, VertexState::Free);
        let tri = DelaunTriangle::new([ab, bc, ca], [a, b, c]);
        let id = ds.add_element(tri);
        assert_eq!(ds.elements_of_domain().len(), 1);
        let pair = ds.elements_connected_to(ab);
        assert_eq!(pair.first_index(), id);
        ds.remove_element(id);
        assert!(ds.elements_of_domain().is_empty());
        assert!(ds.elements_connected_to(ab).is_empty());
        assert_eq!(ds.element_movability(id), VertexState::Deleted);
    }

    #[test]
    fn remove_link_requires_no_elements() {
        let mut ds = DelaunDataStructure::new(10);
        let a = ds.add_node(v(0.0, 0.0));
        let b = ds.add_node(v(1.0, 0.0));
        let ab = ds.add_link(a, b, VertexState::Free);
        assert_eq!(ds.links_of_domain().len(), 1);
        ds.remove_link(ab, false);
        assert!(ds.links_of_domain().is_empty());
        assert_eq!(ds.link_movability(ab), VertexState::Deleted);
    }

    #[test]
    fn clear_domain_removes_elements_and_free_links() {
        let mut ds = DelaunDataStructure::new(10);
        let a = ds.add_node(v(0.0, 0.0));
        let b = ds.add_node(v(1.0, 0.0));
        let c = ds.add_node(v(0.0, 1.0));
        let ab = ds.add_link(a, b, VertexState::Free);
        let bc = ds.add_link(b, c, VertexState::Free);
        let ca = ds.add_link(c, a, VertexState::Free);
        ds.add_element(DelaunTriangle::new([ab, bc, ca], [a, b, c]));
        ds.clear_domain();
        assert!(ds.elements_of_domain().is_empty());
        assert!(ds.links_of_domain().is_empty());
    }

    #[test]
    fn selector_collects_element_neighbours() {
        let mut ds = DelaunDataStructure::new(10);
        let a = ds.add_node(v(0.0, 0.0));
        let b = ds.add_node(v(2.0, 0.0));
        let c = ds.add_node(v(0.0, 2.0));
        let d = ds.add_node(v(2.0, 2.0));
        let ab = ds.add_link(a, b, VertexState::Free);
        let bc = ds.add_link(b, c, VertexState::Free);
        let ca = ds.add_link(c, a, VertexState::Free);
        let cd = ds.add_link(c, d, VertexState::Free);
        let da = ds.add_link(d, a, VertexState::Free);
        let bd = ds.add_link(b, d, VertexState::Free);
        ds.add_element(DelaunTriangle::new([ab, bd, da], [a, b, d]));
        ds.add_element(DelaunTriangle::new([bc, cd, bd], [b, c, d]));
        let mut sel = DelaunSelector::new(&ds);
        sel.neighbours_of_node(a);
        // Triangle (b,c,d) does not touch node a; only the (a,b,d) triangle does.
        assert_eq!(sel.elements().len(), 1);
        let mut sel2 = DelaunSelector::new(&ds);
        sel2.neighbours_of_link(bd);
        assert_eq!(sel2.elements().len(), 2);
    }
}
