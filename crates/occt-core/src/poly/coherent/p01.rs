use super::prelude::*;
use super::*;

/// Ring list of triangle indices incident on a node. Ported from
/// `Poly_CoherentTriPtr`; the ring behaviour is preserved by the circular
/// iteration API (start anywhere, iterate forward).
#[derive(Debug, Clone, Default)]

pub struct PolyCoherentTriPtr {
    pub(super) items: Vec<usize>,
}

impl PolyCoherentTriPtr {
    pub fn new(tri: usize) -> Self { Self { items: vec![tri] } }
    pub fn append(&mut self, tri: usize) { self.items.push(tri); }
    pub fn prepend(&mut self, tri: usize) { self.items.insert(0, tri); }
    /// Removes a triangle reference, returning whether it was present.
    pub fn remove(&mut self, tri: usize) -> bool {
        if let Some(pos) = self.items.iter().position(|&t| t == tri) {
            self.items.remove(pos);
            true
        } else {
            false
        }
    }
    pub fn contains(&self, tri: usize) -> bool { self.items.contains(&tri) }
    pub fn is_empty(&self) -> bool { self.items.is_empty() }
    pub fn len(&self) -> usize { self.items.len() }
    pub fn iter(&self) -> Iter<'_, usize> { self.items.iter() }
}

/// Node of a coherent triangulation. Ported from `Poly_CoherentNode`.
#[derive(Debug, Clone)]
pub struct PolyCoherentNode {
    pub coord: GpXyz,
    pub uv: [f64; 2],
    pub normal: [f64; 3],
    /// Index of the node in the original triangulation (-1 if unset).
    pub index: isize,
    /// Incident triangles.
    pub triangles: PolyCoherentTriPtr,
}

impl PolyCoherentNode {
    pub fn new(pnt: GpXyz) -> Self {
        Self { coord: pnt, uv: [INFINITE; 2], normal: [0.0; 3], index: -1, triangles: PolyCoherentTriPtr::default() }
    }
    pub fn set_uv(&mut self, u: f64, v: f64) { self.uv = [u, v]; }
    pub fn get_u(&self) -> f64 { self.uv[0] }
    pub fn get_v(&self) -> f64 { self.uv[1] }
    pub fn set_normal(&mut self, v: &GpXyz) { self.normal = [v.x, v.y, v.z]; }
    pub fn has_normal(&self) -> bool {
        self.normal[0] * self.normal[0] + self.normal[1] * self.normal[1] + self.normal[2] * self.normal[2]
            > CONFUSION
    }
    pub fn get_normal(&self) -> GpXyz { GpXyz::new(self.normal[0], self.normal[1], self.normal[2]) }
    pub fn set_index(&mut self, i: isize) { self.index = i; }
    pub fn get_index(&self) -> isize { self.index }
    /// A free node is one without any incident triangle.
    pub fn is_free_node(&self) -> bool { self.triangles.is_empty() }
}

/// Link between two mesh nodes, created by existing triangle(s). Ported from
/// `Poly_CoherentLink`. `node[0]` is always the smaller node index.
#[derive(Debug, Clone, Copy)]
pub struct PolyCoherentLink {
    pub node: [isize; 2],
    pub opposite_node: [isize; 2],
    pub attribute: Option<usize>,
}

impl Default for PolyCoherentLink {
    fn default() -> Self { Self::new() }
}

impl PolyCoherentLink {
    /// Empty (invalid) link: all node indices -1.
    pub fn new() -> Self { Self { node: [-1, -1], opposite_node: [-1, -1], attribute: None } }
    /// Temporary link without opposite-node references.
    pub fn new_nodes(n0: isize, n1: isize) -> Self { Self { node: [n0, n1], opposite_node: [-1, -1], attribute: None } }
    /// Builds a link from a triangle and a side (side = index of the node
    /// opposite the link). Always orders `node[0] < node[1]`.
    pub fn from_triangle(tri: &PolyCoherentTriangle, i_side: usize) -> Self {
        let ind = [1usize, 2, 0, 1];
        let a0 = tri.node(ind[i_side]);
        let a1 = tri.node(ind[i_side + 1]);
        if a0 < a1 {
            Self { node: [a0, a1], opposite_node: [tri.node(i_side), tri.get_connected_node(i_side)], attribute: None }
        } else {
            Self { node: [a1, a0], opposite_node: [tri.get_connected_node(i_side), tri.node(i_side)], attribute: None }
        }
    }
    /// `ind` 0 or 1 selects the smaller/larger node of the link.
    pub fn node(&self, ind: usize) -> isize { self.node[ind & 1] }
    /// Opposite node on the left (0) or right (1) incident triangle.
    pub fn opposite_node(&self, ind: usize) -> isize { self.opposite_node[ind & 1] }
    pub fn set_attribute(&mut self, a: Option<usize>) { self.attribute = a; }
    pub fn get_attribute(&self) -> Option<usize> { self.attribute }
    pub fn is_empty(&self) -> bool { self.node[0] < 0 || self.node[1] < 0 }
    pub fn nullify(&mut self) { *self = Self::new(); }
}

/// A triangle with references to its neighbours. Ported from
/// `Poly_CoherentTriangle`. Connection side `i` is opposite node `i`; the
/// shared edge of side `i` connects nodes `(i+1)%3` and `(i+2)%3`.
#[derive(Debug, Clone)]
pub struct PolyCoherentTriangle {
    pub nodes: [isize; 3],
    pub nodes_on_connected: [isize; 3],
    pub connected: [Option<usize>; 3],
    pub links: [Option<usize>; 3],
    pub n_connections: usize,
}

impl Default for PolyCoherentTriangle {
    fn default() -> Self { Self::new() }
}

impl PolyCoherentTriangle {
    /// Empty triangle: all nodes -1.
    pub fn new() -> Self {
        Self { nodes: [-1; 3], nodes_on_connected: [-1; 3], connected: [None; 3], links: [None; 3], n_connections: 0 }
    }
    pub fn new3(n0: isize, n1: isize, n2: isize) -> Self {
        Self { nodes: [n0, n1, n2], nodes_on_connected: [-1; 3], connected: [None; 3], links: [None; 3], n_connections: 0 }
    }
    pub fn node(&self, i: usize) -> isize { self.nodes[i] }
    pub fn is_empty(&self) -> bool { self.nodes[0] < 0 || self.nodes[1] < 0 || self.nodes[2] < 0 }
    pub fn n_connections(&self) -> usize { self.n_connections }
    pub fn get_connected_node(&self, i: usize) -> isize { self.nodes_on_connected[i] }
    pub fn get_connected_tri(&self, i: usize) -> Option<usize> { self.connected[i] }
    pub fn get_link(&self, i: usize) -> Option<usize> { self.links[i] }
    /// Returns the side of the connection with the given triangle, or -1.
    pub fn find_connection(&self, tri_id: usize) -> isize {
        if self.connected[0] == Some(tri_id) { 0 }
        else if self.connected[1] == Some(tri_id) { 1 }
        else if self.connected[2] == Some(tri_id) { 2 }
        else { -1 }
    }
    /// Indices of the up-to-three neighbouring triangles.
    pub fn neighbours(&self) -> [Option<usize>; 3] { self.connected }
    /// Sides that have no neighbour (boundary sides).
    pub fn boundary_sides(&self) -> Vec<usize> {
        (0..3).filter(|&i| self.connected[i].is_none()).collect()
    }
}

/// Coherent triangulation with connectivity. Ported from
/// `Poly_CoherentTriangulation`. All index-based: nodes / triangles / links
/// are stored in `Vec`s and referenced by 0-based indices. Removed triangles
/// are marked empty, never deleted, so indices stay valid.
#[derive(Debug, Clone)]
pub struct PolyCoherentTriangulation {
    pub nodes: Vec<PolyCoherentNode>,
    pub triangles: Vec<PolyCoherentTriangle>,
    pub links: Vec<PolyCoherentLink>,
    pub deflection: f64,
}

impl Default for PolyCoherentTriangulation {
    fn default() -> Self { Self::new() }
}

impl PolyCoherentTriangulation {
    pub fn new() -> Self {
        Self { nodes: Vec::new(), triangles: Vec::new(), links: Vec::new(), deflection: 0.0 }
    }

    // ---- construction / conversion ----

    /// Builds a coherent triangulation from a full `PolyTriangulation`.
    pub fn from_poly_triangulation(tri: &PolyTriangulation) -> Result<Self, String> {
        let mut result = Self::new();
        let n_nodes = tri.nb_nodes();
        for i in 0..n_nodes {
            let id = result.set_node(&tri.node(i).coord, i as isize)?;
            result.nodes[id].set_index(i as isize + 1);
        }
        for t in 0..tri.nb_triangles() {
            let tr = tri.triangle(t);
            if tr.n0 != tr.n1 && tr.n1 != tr.n2 && tr.n2 != tr.n0 {
                result.add_triangle(tr.n0, tr.n1, tr.n2)?;
            }
        }
        if tri.has_uv_nodes() {
            for i in 0..n_nodes {
                let uv = tri.uv_node(i);
                result.nodes[i].set_uv(uv.x(), uv.y());
            }
        }
        if tri.has_normals() {
            for i in 0..n_nodes {
                let n = tri.normal(i);
                result.nodes[i].set_normal(&GpXyz::new(n.x(), n.y(), n.z()));
            }
        }
        result.deflection = tri.deflection();
        Ok(result)
    }

    /// Exports active (non-free) nodes and non-empty triangles to a full
    /// `PolyTriangulation`.
    pub fn get_poly_triangulation(&self) -> Result<PolyTriangulation, String> {
        let n_nodes = self.n_nodes();
        let n_tris = self.n_triangles();
        if n_nodes == 0 || n_tris == 0 {
            return Err("PolyCoherentTriangulation: no active nodes or triangles".into());
        }
        let mut node_map = vec![0usize; self.nodes.len()];
        let mut out = PolyTriangulation::with_capacity(n_nodes, n_tris, false, false);
        let mut count = 0usize;
        let mut any_normal = false;
        let mut any_uv = false;
        for (i, node) in self.nodes.iter().enumerate() {
            if node.is_free_node() { continue; }
            count += 1;
            node_map[i] = count;
            out.set_node(count - 1, GpPnt::from_xyz(&node.coord));
            if node.has_normal() { any_normal = true; }
            if node.get_u() * node.get_u() + node.get_v() * node.get_v() > CONFUSION { any_uv = true; }
        }
        if any_normal { out.add_normals(); }
        if any_uv { out.add_uv_nodes(); }
        for (i, node) in self.nodes.iter().enumerate() {
            if node.is_free_node() { continue; }
            let k = node_map[i] - 1;
            if any_normal { out.set_normal(k, GpPnt::from_xyz(&node.get_normal())); }
            if any_uv { out.set_uv_node(k, GpPnt2d::new(node.get_u(), node.get_v())); }
        }
        count = 0;
        for tri in &self.triangles {
            if tri.is_empty() { continue; }
            let n0 = node_map[tri.nodes[0] as usize];
            let n1 = node_map[tri.nodes[1] as usize];
            let n2 = node_map[tri.nodes[2] as usize];
            if n0 == 0 || n1 == 0 || n2 == 0 {
                return Err("PolyCoherentTriangulation: triangle references a free node".into());
            }
            out.set_triangle(count, super::super::triangulation::Triangle::new(n0 - 1, n1 - 1, n2 - 1));
            count += 1;
        }
        out.set_deflection(self.deflection);
        Ok(out)
    }

    // ---- node management ----

    /// Adds a node at the end (`i_node < 0`) or overwrites node `i_node`.
    /// Returns the node index.
    pub fn set_node(&mut self, pnt: &GpXyz, i_node: isize) -> Result<usize, String> {
        if i_node < 0 {
            let id = self.nodes.len();
            self.nodes.push(PolyCoherentNode::new(*pnt));
            Ok(id)
        } else {
            let i = i_node as usize;
            if i >= self.nodes.len() {
                // Growing like NCollection_DynamicArray::SetValue.
                self.nodes.resize(i + 1, PolyCoherentNode::new(GpXyz::zero()));
            }
            self.nodes[i] = PolyCoherentNode::new(*pnt);
            Ok(i)
        }
    }

    pub fn node(&self, i: usize) -> &PolyCoherentNode { &self.nodes[i] }
    pub fn change_node(&mut self, i: usize) -> &mut PolyCoherentNode { &mut self.nodes[i] }
    /// Index of the last node slot (len - 1).
    pub fn max_node(&self) -> isize { self.nodes.len() as isize - 1 }
    pub fn max_triangle(&self) -> isize { self.triangles.len() as isize - 1 }
    /// Number of active (non-free) nodes.
    pub fn n_nodes(&self) -> usize { self.nodes.iter().filter(|n| !n.is_free_node()).count() }
    /// Number of active (non-empty) triangles.
    pub fn n_triangles(&self) -> usize { self.triangles.iter().filter(|t| !t.is_empty()).count() }
    /// Number of active (non-empty) links.
    pub fn n_links(&self) -> usize { self.links.iter().filter(|l| !l.is_empty()).count() }

    pub fn triangle(&self, i: usize) -> &PolyCoherentTriangle { &self.triangles[i] }
    pub fn link(&self, i: usize) -> &PolyCoherentLink { &self.links[i] }

    /// Indices of free nodes (nodes without incident triangles).
    pub fn get_free_nodes(&self) -> Vec<usize> {
        (0..self.nodes.len()).filter(|&i| self.nodes[i].is_free_node()).collect()
    }

    // ---- connection helpers (mutate triangle adjacency) ----

    /// Removes the connection on `tri_id`'s side `i_conn`, also clearing the
    /// back-reference in the neighbour.
    pub(super) fn tri_remove_connection(&mut self, tri_id: usize, i_conn: usize) {
        let other = self.triangles[tri_id].connected[i_conn];
        if let Some(other_id) = other {
            let mut i_conn1 = 0usize;
            if self.triangles[other_id].connected[0] != Some(tri_id) {
                if self.triangles[other_id].connected[1] == Some(tri_id) { i_conn1 = 1; }
                else if self.triangles[other_id].connected[2] == Some(tri_id) { i_conn1 = 2; }
                else {
                    // Incoherent topology; clear locally and proceed.
                    self.triangles[tri_id].connected[i_conn] = None;
                    self.triangles[tri_id].nodes_on_connected[i_conn] = -1;
                    self.triangles[tri_id].n_connections = self.triangles[tri_id].n_connections.saturating_sub(1);
                    return;
                }
            }
            self.triangles[other_id].connected[i_conn1] = None;
            self.triangles[other_id].nodes_on_connected[i_conn1] = -1;
            self.triangles[other_id].n_connections = self.triangles[other_id].n_connections.saturating_sub(1);
            self.triangles[tri_id].connected[i_conn] = None;
            self.triangles[tri_id].nodes_on_connected[i_conn] = -1;
            self.triangles[tri_id].n_connections = self.triangles[tri_id].n_connections.saturating_sub(1);
        }
    }

    /// Connects two triangles on given sides. Removes any pre-existing
    /// connections on those sides first.
    pub(super) fn tri_connect_sides(&mut self, a: usize, side_a: usize, b: usize, side_b: usize) {
        let node_a_opp = self.triangles[a].nodes[side_a];
        let node_b_opp = self.triangles[b].nodes[side_b];
        self.tri_remove_connection(a, side_a);
        self.tri_remove_connection(b, side_b);
        self.triangles[a].nodes_on_connected[side_a] = node_b_opp;
        self.triangles[a].connected[side_a] = Some(b);
        self.triangles[a].n_connections += 1;
        self.triangles[b].nodes_on_connected[side_b] = node_a_opp;
        self.triangles[b].connected[side_b] = Some(a);
        self.triangles[b].n_connections += 1;
    }

    /// Explicit connection: connect triangle `a` on side `side_a` to triangle
    /// `b`, whose matching side is deduced. Port of
    /// `Poly_CoherentTriangle::SetConnection(int, Triangle&)`.
    pub fn tri_set_connection(&mut self, a: usize, side_a: usize, b: usize) -> bool {
        pub(super) const II: [usize; 5] = [2, 0, 1, 2, 0];
        let na = self.triangles[a].nodes;
        let nb = self.triangles[b].nodes;
        if nb[0] == na[II[side_a + 2]] {
            if nb[2] != na[II[side_a]] { return false; }
            self.tri_connect_sides(a, side_a, b, 1);
        } else if nb[1] == na[II[side_a + 2]] {
            if nb[0] != na[II[side_a]] { return false; }
            self.tri_connect_sides(a, side_a, b, 2);
        } else if nb[2] == na[II[side_a + 2]] {
            if nb[1] != na[II[side_a]] { return false; }
            self.tri_connect_sides(a, side_a, b, 0);
        } else {
            return false;
        }
        true
    }

    /// Automatic connection: analyse which sides of the two triangles share an
    /// edge and connect them. Port of
    /// `Poly_CoherentTriangle::SetConnection(Triangle&)`.
    pub fn tri_set_connection_auto(&mut self, a: usize, b: usize) -> bool {
        let na = self.triangles[a].nodes;
        let nb = self.triangles[b].nodes;
        if na[0] == nb[0] {
            if na[1] == nb[2] && self.triangles[a].connected[2] != Some(b) {
                self.tri_connect_sides(a, 2, b, 1); return true;
            } else if na[2] == nb[1] && self.triangles[a].connected[1] != Some(b) {
                self.tri_connect_sides(a, 1, b, 2); return true;
            }
        } else if na[0] == nb[1] {
            if na[1] == nb[0] && self.triangles[a].connected[2] != Some(b) {
                self.tri_connect_sides(a, 2, b, 2); return true;
            } else if na[2] == nb[2] && self.triangles[a].connected[1] != Some(b) {
                self.tri_connect_sides(a, 1, b, 0); return true;
            }
        } else if na[0] == nb[2] {
            if na[1] == nb[1] && self.triangles[a].connected[2] != Some(b) {
                self.tri_connect_sides(a, 2, b, 0); return true;
            } else if na[2] == nb[0] && self.triangles[a].connected[1] != Some(b) {
                self.tri_connect_sides(a, 1, b, 1); return true;
            }
        } else if self.triangles[a].connected[0] != Some(b) {
            if na[1] == nb[0] && na[2] == nb[2] {
                self.tri_connect_sides(a, 0, b, 1); return true;
            } else if na[1] == nb[2] && na[2] == nb[1] {
                self.tri_connect_sides(a, 0, b, 0); return true;
            } else if na[1] == nb[1] && na[2] == nb[0] {
                self.tri_connect_sides(a, 0, b, 2); return true;
            }
        }
        false
    }

    // ---- triangles ----

    /// Appends a new triangle and builds its adjacency. Returns its index.
    pub fn add_triangle(&mut self, n0: usize, n1: usize, n2: usize) -> Result<usize, String> {
        if n0 >= self.nodes.len() || n1 >= self.nodes.len() || n2 >= self.nodes.len() {
            return Err(format!("add_triangle: node index out of range"));
        }
        if n0 == n1 || n1 == n2 || n2 == n0 {
            return Err("add_triangle: degenerate triangle".into());
        }
        let id = self.triangles.len();
        self.triangles.push(PolyCoherentTriangle::new3(n0 as isize, n1 as isize, n2 as isize));
        self.replace_nodes(id, n0, n1, n2)?;
        Ok(id)
    }

    /// Replaces the nodes of a triangle, rebuilding adjacency and links.
    pub fn replace_nodes(&mut self, tri_id: usize, n0: usize, n1: usize, n2: usize) -> Result<(), String> {
        if tri_id >= self.triangles.len() {
            return Err(format!("replace_nodes: triangle index {tri_id} out of range"));
        }
        if n0 >= self.nodes.len() || n1 >= self.nodes.len() || n2 >= self.nodes.len() {
            return Err(format!("replace_nodes: node index out of range"));
        }
        if !self.triangles[tri_id].is_empty() {
            self.remove_triangle(tri_id)?;
        }
        self.triangles[tri_id] = PolyCoherentTriangle::new3(n0 as isize, n1 as isize, n2 as isize);
        for i in 0..3 {
            let node_id = self.triangles[tri_id].nodes[i] as usize;
            let incident: Vec<usize> = self.nodes[node_id].triangles.iter().cloned().collect();
            for &t_other in &incident {
                self.tri_set_connection_auto(t_other, tri_id);
            }
            self.nodes[node_id].triangles.append(tri_id);
        }
        // If links exist, create or update them.
        if !self.links.is_empty() {
            for i in 0..3 {
                let tri_opp = self.triangles[tri_id].connected[i];
                let mut to_add_link = true;
                if let Some(opp_id) = tri_opp {
                    for j in 0..3 {
                        if self.triangles[tri_id].nodes[i] == self.triangles[opp_id].nodes_on_connected[j] {
                            if let Some(link_id) = self.triangles[opp_id].links[j] {
                                if self.links[link_id].opposite_node(0) == self.triangles[opp_id].nodes[j] {
                                    self.links[link_id].opposite_node[1] = self.triangles[tri_id].nodes[i];
                                } else if self.links[link_id].opposite_node(1) == self.triangles[opp_id].nodes[j] {
                                    self.links[link_id].opposite_node[0] = self.triangles[tri_id].nodes[i];
                                }
                                to_add_link = false;
                            }
                        }
                    }
                }
                if to_add_link {
                    self.add_link(tri_id, i)?;
                }
            }
        }
        Ok(())
    }

    /// Removes a triangle, disconnecting it from its neighbours and nodes.
    /// Returns whether any node reference was cleared.
    pub fn remove_triangle(&mut self, tri_id: usize) -> Result<bool, String> {
        if tri_id >= self.triangles.len() {
            return Err(format!("remove_triangle: index {tri_id} out of range"));
        }
        let mut removed = false;
        for i in 0..3 {
            let node_id = self.triangles[tri_id].nodes[i];
            if node_id >= 0 {
                let nid = node_id as usize;
                // Maintain links.
                if let Some(link_id) = self.triangles[tri_id].links[i] {
                    let tri_opp = self.triangles[tri_id].connected[i];
                    let mut to_remove_link = true;
                    if let Some(opp_id) = tri_opp {
                        for j in 0..3 {
                            if self.triangles[opp_id].links[j] == Some(link_id) {
                                if self.links[link_id].opposite_node(0) == node_id {
                                    self.links[link_id].opposite_node[0] = -1;
                                    to_remove_link = false;
                                } else if self.links[link_id].opposite_node(1) == node_id {
                                    self.links[link_id].opposite_node[1] = -1;
                                    to_remove_link = false;
                                }
                                break;
                            }
                        }
                    }
                    if to_remove_link {
                        self.remove_link(link_id);
                    }
                }
                if self.nodes[nid].triangles.remove(tri_id) {
                    self.triangles[tri_id].nodes[i] = -1;
                    removed = true;
                }
            }
            self.tri_remove_connection(tri_id, i);
        }
        Ok(removed)
    }

    // ---- links ----

    /// Appends a link from the given triangle side. Returns the link index.
    pub fn add_link(&mut self, tri_id: usize, conn: usize) -> Result<usize, String> {
        if self.triangles[tri_id].is_empty() {
            return Err("add_link: empty triangle".into());
        }
        let link = PolyCoherentLink::from_triangle(&self.triangles[tri_id], conn);
        let link_id = self.links.len();
        self.links.push(link);
        self.triangles[tri_id].links[conn] = Some(link_id);
        let conn_node = self.triangles[tri_id].nodes_on_connected[conn];
        if let Some(opp_id) = self.triangles[tri_id].connected[conn] {
            if !self.triangles[opp_id].is_empty() {
                for j in 0..3 {
                    if self.triangles[opp_id].nodes[j] == conn_node {
                        self.triangles[opp_id].links[j] = Some(link_id);
                        break;
                    }
                }
            }
        }
        Ok(link_id)
    }

    /// Removes a link, clearing its references in the incident triangles.
    pub fn remove_link(&mut self, link_id: usize) {
        let link = self.links[link_id];
        if link.is_empty() { return; }
        let ptri = self.find_triangle(&link);
        for i in 0..2 {
            let i_node = link.opposite_node[i];
            if i_node >= 0 {
                if let Some(tri_id) = ptri[i] {
                    for j in 0..3 {
                        if self.triangles[tri_id].links[j] == Some(link_id) {
                            self.triangles[tri_id].links[j] = None;
                            break;
                        }
                    }
                }
            }
        }
        self.links[link_id] = PolyCoherentLink::new();
    }

    /// Finds the two triangles sharing the nodes of `link`:
    /// `[left, right]` relative to the directed edge `node[0] -> node[1]`.
    pub fn find_triangle(&self, link: &PolyCoherentLink) -> [Option<usize>; 2] {
        let mut res: [Option<usize>; 2] = [None, None];
        if link.is_empty() { return res; }
        let n0 = link.node[0] as usize;
        let n1 = link.node[1] as usize;
        if n0 >= self.nodes.len() || n1 >= self.nodes.len() { return res; }
        for &tri_id in self.nodes[n0].triangles.iter() {
            let t = &self.triangles[tri_id];
            if t.nodes[0] == n0 as isize {
                if t.nodes[1] == n1 as isize { res[0] = Some(tri_id); }
                else if t.nodes[2] == n1 as isize { res[1] = Some(tri_id); }
            } else if t.nodes[1] == n0 as isize {
                if t.nodes[2] == n1 as isize { res[0] = Some(tri_id); }
                else if t.nodes[0] == n1 as isize { res[1] = Some(tri_id); }
            } else if t.nodes[2] == n0 as isize {
                if t.nodes[0] == n1 as isize { res[0] = Some(tri_id); }
                else if t.nodes[1] == n1 as isize { res[1] = Some(tri_id); }
            }
            if res[0].is_some() && res[1].is_some() { break; }
        }
        res
    }

    /// (Re)computes all links from the triangle adjacency. Returns the number
    /// of links.
    pub fn compute_links(&mut self) -> Result<usize, String> {
        self.links.clear();
        for t in 0..self.triangles.len() {
            let nodes = self.triangles[t].nodes;
            if self.triangles[t].is_empty() { continue; }
            if nodes[0] < nodes[1] { self.add_link(t, 2)?; }
            if nodes[1] < nodes[2] { self.add_link(t, 0)?; }
            if nodes[2] < nodes[0] { self.add_link(t, 1)?; }
        }
        // The pass above does not create all boundary links; add missing ones.
        for t in 0..self.triangles.len() {
            if self.triangles[t].is_empty() { continue; }
            for i in 0..3 {
                if self.triangles[t].links[i].is_none() {
                    self.add_link(t, i)?;
                }
            }
        }
        Ok(self.links.len())
    }

    /// Clears all links and nullifies link references in triangles.
    pub fn clear_links(&mut self) {
        self.links.clear();
        for t in &mut self.triangles {
            t.links = [None; 3];
        }
    }

    // ---- mesh analysis ----

    /// Finds and removes one degenerated triangle (with an edge shorter than
    /// `tol`), reconnecting the mesh. Returns `(removed_node, remaining_node)`
    /// pairs when a degeneration is removed.
    pub fn remove_degenerated(&mut self, tol: f64) -> Result<Vec<(usize, usize)>, String> {
        let tol2 = tol * tol;
        let mut result = Vec::new();
        let tri_ids: Vec<usize> = (0..self.triangles.len()).collect();
        for tri_id in tri_ids {
            if self.triangles[tri_id].is_empty() { continue; }
            let nodes = self.triangles[tri_id].nodes;
            let p = [nodes[0] as usize, nodes[1] as usize, nodes[2] as usize];
            let len2 = [
                self.nodes[p[2]].coord.subtracted(&self.nodes[p[1]].coord).square_modulus(),
                self.nodes[p[0]].coord.subtracted(&self.nodes[p[2]].coord).square_modulus(),
                self.nodes[p[1]].coord.subtracted(&self.nodes[p[0]].coord).square_modulus(),
            ];
            for i in 0..3 {
                if len2[i] < tol2 {
                    // Edge opposite node i is short: collapse node ip1 onto im1.
                    let ip1 = nodes[(i + 1) % 3] as usize;
                    let im1 = nodes[(i + 2) % 3] as usize;
                    self.remove_triangle(tri_id)?;
                    loop {
                        let incident: Vec<usize> = self.nodes[ip1].triangles.iter().cloned().collect();
                        if incident.is_empty() { break; }
                        let tri_conn = incident[0];
                        let mut nn = self.triangles[tri_conn].nodes;
                        if nn[0] == ip1 as isize { nn[0] = im1 as isize; }
                        else if nn[1] == ip1 as isize { nn[1] = im1 as isize; }
                        else if nn[2] == ip1 as isize { nn[2] = im1 as isize; }
                        self.remove_triangle(tri_conn)?;
                        self.add_triangle(nn[0] as usize, nn[1] as usize, nn[2] as usize)?;
                    }
                    result.push((ip1, im1));
                    return Ok(result);
                }
            }
        }
        Ok(result)
    }

    /// Enumerates boundary loops. Each returned vector is the sequence of node
    /// indices around one closed boundary. Works by collecting edges that
    /// appear in exactly one triangle and walking them.
    pub fn boundary_loops(&self) -> Vec<Vec<usize>> {
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &self.triangles {
            if tri.is_empty() { continue; }
            for k in 0..3 {
                let a = tri.nodes[k] as usize;
                let b = tri.nodes[(k + 1) % 3] as usize;
                let (mn, mx) = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry((mn, mx)).or_insert(0) += 1;
            }
        }
        let boundary: Vec<(usize, usize)> = edge_count.iter()
            .filter(|&(_, &c)| c == 1)
            .map(|(&e, _)| e)
            .collect();
        if boundary.is_empty() { return Vec::new(); }

        let mut adj: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        for &(a, b) in &boundary {
            adj.entry(a).or_default().push((a, b));
            adj.entry(b).or_default().push((a, b));
        }
        let mut used = vec![false; boundary.len()];
        let mut loops = Vec::new();
        for start in 0..boundary.len() {
            if used[start] { continue; }
            let (a0, b0) = boundary[start];
            let start_node = a0;
            let mut loop_nodes = vec![start_node];
            used[start] = true;
            let mut cur = b0;
            loop {
                loop_nodes.push(cur);
                if cur == start_node { break; }
                let mut next_edge: Option<(usize, usize, usize)> = None;
                if let Some(edges) = adj.get(&cur) {
                    for &(x, y) in edges {
                        if let Some(idx) = boundary.iter().position(|&e| e == (x, y)) {
                            if !used[idx] {
                                next_edge = Some((idx, x, y));
                                break;
                            }
                        }
                    }
                }
                match next_edge {
                    Some((idx, x, y)) => {
                        used[idx] = true;
                        cur = if x == cur { y } else { x };
                    }
                    None => break,
                }
            }
            loops.push(loop_nodes);
        }
        loops
    }
}
