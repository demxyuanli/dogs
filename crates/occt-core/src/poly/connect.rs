//! Exploration of adjacency data inside a triangulation: for each node, the
//! triangles that contain it; for each triangle, the three adjacent triangles
//! and the three adjacent nodes.
//! Source: `Poly_Connect.hxx` / `Poly_Connect.cxx`.

use std::collections::HashSet;

use super::triangulation_full::PolyTriangulation;

/// Edge record used while building the adjacency tables.
#[derive(Debug, Clone, Copy)]
struct EdgeRec {
    /// Two triangles sharing the edge; `NONE` = none.
    nt: [usize; 2],
    /// Third node of each sharing triangle; `NONE` = none.
    nn: [usize; 2],
}

impl Default for EdgeRec {
    fn default() -> Self { Self { nt: [NONE; 2], nn: [NONE; 2] } }
}

/// Provides an algorithm to explore, inside a triangulation, the adjacency
/// data for a node or a triangle.
///
/// Adjacency data for a node consists of triangles which contain the node.
/// Adjacency data for a triangle consists of the 3 adjacent triangles which
/// share an edge of the triangle, and the 3 nodes which are the other nodes of
/// these adjacent triangles.
///
/// All indices are 0-based. A value of `usize::MAX` (`NONE`) denotes absence.
#[derive(Debug, Clone)]
pub struct PolyConnect {
    triangulation: Option<PolyTriangulation>,
    /// Per-node: one triangle containing the node (`NONE` if none).
    my_triangles: Vec<usize>,
    /// Per-triangle: 3 adjacent triangles + 3 adjacent nodes (6 slots).
    my_adjacents: Vec<usize>,
    // Iterator state.
    mytr: usize,
    myfirst: usize,
    mynode: usize,
    myothernode: usize,
    mysense: bool,
    mymore: bool,
    my_passed_tr: HashSet<usize>,
}

/// Sentinel for "no triangle / no node" (mirrors OCCT's 0 in 1-based tables).
pub const NONE: usize = usize::MAX;

impl Default for PolyConnect {
    fn default() -> Self { Self::new() }
}

impl PolyConnect {
    /// Constructs an uninitialized algorithm.
    pub fn new() -> Self {
        Self {
            triangulation: None,
            my_triangles: Vec::new(),
            my_adjacents: Vec::new(),
            mytr: 0,
            myfirst: 0,
            mynode: 0,
            myothernode: 0,
            mysense: false,
            mymore: false,
            my_passed_tr: HashSet::new(),
        }
    }

    /// Constructs an algorithm to explore the adjacency data of the
    /// triangulation.
    pub fn from_triangulation(triangulation: &PolyTriangulation) -> Self {
        let mut c = Self::new();
        c.load(triangulation);
        c
    }

    /// (Re)initializes the algorithm for the given triangulation.
    pub fn load(&mut self, triangulation: &PolyTriangulation) {
        self.triangulation = Some(triangulation.clone());
        self.mytr = 0;
        self.myfirst = 0;
        self.mynode = 0;
        self.myothernode = 0;
        self.mysense = false;
        self.mymore = false;

        let t = &self.triangulation.as_ref().expect("triangulation");
        let nb_nodes = t.nb_nodes();
        let nb_tris = t.nb_triangles();
        self.my_triangles = vec![NONE; nb_nodes];
        self.my_adjacents = vec![NONE; 6 * nb_tris];

        // Build an edge map (unordered node pair -> up to two triangles).
        let mut edges: std::collections::HashMap<(usize, usize), EdgeRec> = std::collections::HashMap::new();
        for tri_iter in 0..nb_tris {
            let tri = t.triangle(tri_iter);
            let n = [tri.n0, tri.n1, tri.n2];
            // Record one triangle per node (last one wins, matching OCCT).
            self.my_triangles[n[0]] = tri_iter;
            self.my_triangles[n[1]] = tri_iter;
            self.my_triangles[n[2]] = tri_iter;
            for node_in_tri in 0..3 {
                let node_next = (node_in_tri + 1) % 3;
                let (a, b) = (n[node_in_tri], n[node_next]);
                let (mn, mx) = if a < b { (a, b) } else { (b, a) };
                let third = n[3 - node_in_tri - node_next];
                let e = edges.entry((mn, mx)).or_default();
                if e.nt[0] == NONE {
                    e.nt[0] = tri_iter;
                    e.nn[0] = third;
                } else if e.nt[1] == NONE {
                    e.nt[1] = tri_iter;
                    e.nn[1] = third;
                }
            }
        }

        // Fill the adjacent tables.
        let mut adj_index = 0usize;
        for tri_iter in 0..nb_tris {
            let tri = t.triangle(tri_iter);
            let n = [tri.n0, tri.n1, tri.n2];
            for node_in_tri in 0..3 {
                let node_next = (node_in_tri + 1) % 3;
                let (a, b) = (n[node_in_tri], n[node_next]);
                let (mn, mx) = if a < b { (a, b) } else { (b, a) };
                let e = edges[&(mn, mx)];
                // The other triangle (l = 0 if nt[0] is this one, else 1).
                let l = if e.nt[0] == tri_iter { 1 } else { 0 };
                self.my_adjacents[adj_index] = e.nt[l];
                self.my_adjacents[adj_index + 3] = e.nn[l];
                adj_index += 1;
            }
            adj_index += 3;
        }
    }

    /// Returns the triangulation analyzed by this tool.
    pub fn triangulation(&self) -> Option<&PolyTriangulation> { self.triangulation.as_ref() }

    /// Returns the index of a triangle containing the node at index `n`
    /// (0-based), or `NONE` if the node is isolated.
    pub fn triangle(&self, n: usize) -> usize {
        if n < self.my_triangles.len() { self.my_triangles[n] } else { NONE }
    }

    /// Returns the indices of the 3 triangles adjacent to the triangle at
    /// index `t` (0-based). `NONE` is returned when there are fewer than 3
    /// adjacent triangles.
    pub fn adjacent_triangles(&self, t: usize) -> [usize; 3] {
        let index = 6 * t;
        [self.my_adjacents[index], self.my_adjacents[index + 1], self.my_adjacents[index + 2]]
    }

    /// Returns the indices of the 3 nodes adjacent to the triangle at index
    /// `t` (0-based): the third node of each adjacent triangle. `NONE` when
    /// there are fewer than 3.
    pub fn adjacent_nodes(&self, t: usize) -> [usize; 3] {
        let index = 6 * t;
        [self.my_adjacents[index + 3], self.my_adjacents[index + 4], self.my_adjacents[index + 5]]
    }

    // ---- iterator over the triangles containing a node ----

    /// Initializes the iterator to search for all the triangles containing the
    /// node at index `n` (0-based).
    pub fn initialize(&mut self, n: usize) {
        self.mynode = n;
        self.myfirst = self.triangle(n);
        self.mytr = self.myfirst;
        self.mysense = true;
        self.mymore = self.myfirst != NONE;
        self.my_passed_tr.clear();
        self.my_passed_tr.insert(self.mytr);
        if self.mymore {
            let no = self.tri_nodes(self.myfirst);
            let mut i = 0;
            while no[i] != self.mynode { i += 1; }
            self.myothernode = no[(i + 2) % 3];
        }
    }

    /// Returns true if there is another element in the iterator.
    pub fn more(&self) -> bool { self.mymore }

    /// Advances the iterator to the next triangle containing the node.
    pub fn next(&mut self) {
        let mut t = self.adjacent_triangles(self.mytr);
        if self.mysense {
            for i in 0..3 {
                if t[i] != NONE {
                    let n = self.tri_nodes(t[i]);
                    for j in 0..3 {
                        if n[j] == self.mynode && n[(j + 1) % 3] == self.myothernode {
                            self.mytr = t[i];
                            self.myothernode = n[(j + 2) % 3];
                            self.mymore = !self.my_passed_tr.contains(&self.mytr);
                            self.my_passed_tr.insert(self.mytr);
                            return;
                        }
                    }
                }
            }
            // Otherwise depart towards the left.
            let first_nodes = self.tri_nodes(self.myfirst);
            let mut i = 0;
            while first_nodes[i] != self.mynode { i += 1; }
            self.myothernode = first_nodes[(i + 1) % 3];
            self.mysense = false;
            self.mytr = self.myfirst;
            t = self.adjacent_triangles(self.mytr);
        }
        if !self.mysense {
            for i in 0..3 {
                if t[i] != NONE {
                    let n = self.tri_nodes(t[i]);
                    for j in 0..3 {
                        if n[j] == self.mynode && n[(j + 2) % 3] == self.myothernode {
                            self.mytr = t[i];
                            self.myothernode = n[(j + 1) % 3];
                            self.mymore = !self.my_passed_tr.contains(&self.mytr);
                            self.my_passed_tr.insert(self.mytr);
                            return;
                        }
                    }
                }
            }
        }
        self.mymore = false;
    }

    /// Returns the index of the current triangle to which the iterator points.
    pub fn value(&self) -> usize { self.mytr }

    fn tri_nodes(&self, t: usize) -> [usize; 3] {
        let tri = self.triangulation.as_ref().expect("triangulation").triangle(t);
        [tri.n0, tri.n1, tri.n2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::GpPnt;
    use crate::poly::triangulation::Triangle;

    /// Unit cube: 8 vertices, 12 triangles (2 per face).
    fn unit_cube() -> PolyTriangulation {
        let nodes = vec![
            GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 1.), GpPnt::new(1., 0., 1.), GpPnt::new(1., 1., 1.), GpPnt::new(0., 1., 1.),
        ];
        let tris = vec![
            Triangle::new(0, 1, 2), Triangle::new(0, 2, 3),
            Triangle::new(4, 6, 5), Triangle::new(4, 7, 6),
            Triangle::new(0, 5, 1), Triangle::new(0, 4, 5),
            Triangle::new(3, 2, 6), Triangle::new(3, 6, 7),
            Triangle::new(0, 3, 7), Triangle::new(0, 7, 4),
            Triangle::new(1, 5, 6), Triangle::new(1, 6, 2),
        ];
        PolyTriangulation::from_parts(nodes, tris)
    }

    #[test]
    fn cube_node_triangle_counts() {
        let tri = unit_cube();
        let mut c = PolyConnect::from_triangulation(&tri);
        // Per-node incident-triangle counts for this cube triangulation.
        let expected = [6usize, 4, 4, 4, 4, 4, 6, 4];
        for n in 0..8 {
            let mut count = 0;
            c.initialize(n);
            while c.more() {
                count += 1;
                c.next();
            }
            assert_eq!(count, expected[n], "node {n}");
        }
    }

    #[test]
    fn cube_triangle_adjacency() {
        let tri = unit_cube();
        let c = PolyConnect::from_triangulation(&tri);
        // Interior edges have two triangles, so every triangle has 3 neighbours.
        for t in 0..12 {
            let adj = c.adjacent_triangles(t);
            assert!(adj.iter().all(|&a| a != NONE), "triangle {t} fully adjacent");
            let nodes = c.adjacent_nodes(t);
            assert!(nodes.iter().all(|&n| n != NONE));
            // Each adjacent triangle shares exactly two nodes with the triangle.
            let tt = tri.triangle(t);
            let tset: HashSet<usize> = [tt.n0, tt.n1, tt.n2].into_iter().collect();
            for k in 0..3 {
                let at = tri.triangle(adj[k]);
                let aset: HashSet<usize> = [at.n0, at.n1, at.n2].into_iter().collect();
                assert_eq!(aset.intersection(&tset).count(), 2);
                // The adjacent node is the third node of the adjacent triangle.
                assert!(!tset.contains(&nodes[k]));
                assert!(aset.contains(&nodes[k]));
            }
        }
    }

    #[test]
    fn triangle_containing_node_valid() {
        let tri = unit_cube();
        let c = PolyConnect::from_triangulation(&tri);
        for n in 0..8 {
            let t = c.triangle(n);
            assert_ne!(t, NONE);
            let tt = tri.triangle(t);
            assert!([tt.n0, tt.n1, tt.n2].contains(&n));
        }
    }

    #[test]
    fn iterator_walks_all_incident() {
        let tri = unit_cube();
        let mut c = PolyConnect::from_triangulation(&tri);
        c.initialize(0);
        let mut seen = Vec::new();
        while c.more() {
            let t = c.value();
            let tt = tri.triangle(t);
            assert!([tt.n0, tt.n1, tt.n2].contains(&0));
            seen.push(t);
            c.next();
        }
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 6); // node 0 is incident to 6 triangles
    }

    #[test]
    fn isolated_node_returns_none() {
        let nodes = vec![GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.), GpPnt::new(5., 5., 5.)];
        let tris = vec![Triangle::new(0, 1, 2)];
        let tri = PolyTriangulation::from_parts(nodes, tris);
        let mut c = PolyConnect::from_triangulation(&tri);
        assert_eq!(c.triangle(3), NONE);
        c.initialize(3);
        assert!(!c.more());
    }
}
