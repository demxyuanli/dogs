//! Split a set of connected faces into closed shells.
//!
//! Ports `BOPAlgo_ShellSplitter` (TKBO). The input is a set of start
//! elements (faces or shells). `perform()` groups the faces into connexity
//! blocks — maximal sets of faces connected through shared edges — and turns
//! each block into one or more closed shells.
//!
//! Face adjacency is decided purely topologically: two faces belong to the
//! same block when they share an edge, and an edge is identified by the
//! coordinate pair of its two endpoints (order independent), taken from
//! `topo_tools_full::{edge_vertices, vertex_position}`.

use std::collections::{HashMap, HashSet};

use crate::builder::TopoBuilder;
use crate::connexity_block::ConnexityBlock;
use crate::shape::{Edge, Face, Shell, TopoShape, Vertex};
use crate::topo_tools_full::{edge_vertices, edges_of, faces_of, is_same, vertex_position};

/// Vertex identity: 3D coordinates quantized onto a `1e-6` grid.
pub(crate) type VKey = (i64, i64, i64);
/// Undirected edge identity: the (sorted) endpoint pair.
pub(crate) type EKey = (VKey, VKey);

/// Grid size used to quantize vertex coordinates for identity tests.
const VERTEX_TOL: f64 = 1e-6;

/// Quantized identity key of a vertex (by its registered position).
pub(crate) fn vertex_key(v: &Vertex) -> VKey {
    let p = vertex_position(v);
    (
        (p.x() / VERTEX_TOL).round() as i64,
        (p.y() / VERTEX_TOL).round() as i64,
        (p.z() / VERTEX_TOL).round() as i64,
    )
}

/// Undirected identity key of an edge: its two endpoint keys, sorted.
/// Direction does not matter — the same geometric edge, however oriented,
/// produces the same key.
pub(crate) fn edge_key(e: &Edge) -> EKey {
    let (a, b) = edge_vertices(e);
    let (Some(a), Some(b)) = (a, b) else {
        return ((0, 0, 0), (0, 0, 0));
    };
    let (ka, kb) = (vertex_key(&a), vertex_key(&b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// Set of edge keys of a face (every distinct boundary edge).
fn face_edge_keys(f: &Face) -> HashSet<EKey> {
    edges_of(&f.0).into_iter().map(|e| edge_key(&e)).collect()
}

/// Build a shell containing `faces` (in order) and mark it closed.
fn make_shell_from_faces(faces: &[Face]) -> Shell {
    let bld = TopoBuilder::new();
    let mut shell = Shell::new();
    for f in faces {
        bld.add_face(&mut shell, f);
    }
    shell.0.set_closed(true);
    shell
}

/// Splits a set of connected faces into closed shells.
///
/// Mirrors `BOPAlgo_ShellSplitter`:
/// - regular connexity blocks (every edge shared by exactly two faces) become
///   a single closed shell directly;
/// - non-regular blocks go through [`ShellSplitter::split_block`], which
///   drops the faces that carry free (boundary) edges and then extracts the
///   edge-connected closed components.
pub struct ShellSplitter {
    /// Start elements (faces / shells) given by the caller.
    start_shapes: Vec<TopoShape>,
    /// Resulting closed shells.
    shells: Vec<TopoShape>,
    /// Connexity blocks built during `perform`.
    blocks: Vec<ConnexityBlock>,
}

impl Default for ShellSplitter {
    fn default() -> Self {
        Self::new()
    }
}

impl ShellSplitter {
    /// Empty splitter.
    pub fn new() -> Self {
        Self {
            start_shapes: Vec::new(),
            shells: Vec::new(),
            blocks: Vec::new(),
        }
    }

    /// Add a start element (a face, a shell, or any shape containing faces).
    pub fn add_start_element(&mut self, shape: TopoShape) {
        self.start_shapes.push(shape);
    }

    /// The start elements given to the splitter.
    pub fn start_elements(&self) -> &[TopoShape] {
        &self.start_shapes
    }

    /// The closed shells built by `perform`.
    pub fn shells(&self) -> &[TopoShape] {
        &self.shells
    }

    /// The connexity blocks produced by `perform` (useful for inspection).
    pub fn blocks(&self) -> &[ConnexityBlock] {
        &self.blocks
    }

    /// Expand the start elements into a deduplicated list of faces.
    fn collect_faces(&self) -> Vec<Face> {
        let mut raw: Vec<Face> = Vec::new();
        for s in &self.start_shapes {
            if s.is_face() {
                raw.push(Face(s.clone()));
            } else {
                raw.extend(faces_of(s));
            }
        }
        // Deduplicate by TShape identity (a face added twice is the same face).
        let mut out: Vec<Face> = Vec::new();
        for f in raw {
            if !out.iter().any(|g| is_same(&g.0, &f.0)) {
                out.push(f);
            }
        }
        out
    }

    /// Group faces into maximal connected blocks (adjacent faces share an
    /// edge), mirroring `BOPTools_AlgoTools::MakeConnexityBlocks`.
    fn make_connexity_blocks(&self, faces: &[Face]) -> Vec<ConnexityBlock> {
        let n = faces.len();
        if n == 0 {
            return Vec::new();
        }
        let face_keys: Vec<HashSet<EKey>> = faces.iter().map(face_edge_keys).collect();
        // edge -> faces containing it
        let mut edge_faces: HashMap<EKey, Vec<usize>> = HashMap::new();
        for (i, keys) in face_keys.iter().enumerate() {
            for &k in keys {
                edge_faces.entry(k).or_default().push(i);
            }
        }
        let mut visited = vec![false; n];
        let mut blocks = Vec::new();
        for start in 0..n {
            if visited[start] {
                continue;
            }
            // Breadth-first traversal over shared edges.
            let mut comp: Vec<usize> = Vec::new();
            let mut stack = vec![start];
            visited[start] = true;
            while let Some(i) = stack.pop() {
                comp.push(i);
                for &k in &face_keys[i] {
                    if let Some(neigh) = edge_faces.get(&k) {
                        for &j in neigh {
                            if !visited[j] {
                                visited[j] = true;
                                stack.push(j);
                            }
                        }
                    }
                }
            }
            let comp_set: HashSet<usize> = comp.iter().cloned().collect();
            let mut block = ConnexityBlock::new();
            for &i in &comp {
                block.change_shapes_mut().push(faces[i].0.clone());
            }
            // Regularity: every edge of every face of the block is shared by
            // exactly two faces of the block.
            let mut regular = true;
            for &i in &comp {
                for &k in &face_keys[i] {
                    let cnt = edge_faces
                        .get(&k)
                        .map_or(0, |v| v.iter().filter(|&&j| comp_set.contains(&j)).count());
                    if cnt != 2 {
                        regular = false;
                    }
                }
            }
            block.set_regular(regular);
            blocks.push(block);
        }
        blocks
    }

    /// Run the splitting: build connexity blocks from the start elements and
    /// turn each block into closed shells.
    pub fn perform(&mut self) -> Result<(), String> {
        self.shells.clear();
        self.blocks.clear();
        let faces = self.collect_faces();
        let blocks = self.make_connexity_blocks(&faces);
        for mut block in blocks {
            if block.is_regular() {
                // Fast path (BOPAlgo_ShellSplitter::MakeShells): a regular
                // block is a single closed shell.
                let faces: Vec<Face> = block.shapes().iter().map(|s| Face(s.clone())).collect();
                let shell = make_shell_from_faces(&faces);
                block.change_loops_mut().push(shell.0);
            } else {
                // General path: split the block into closed shells.
                Self::split_block(&mut block);
            }
            self.shells.extend(block.loops().iter().cloned());
            self.blocks.push(block);
        }
        Ok(())
    }

    /// Split one connexity block into a set of closed shells.
    ///
    /// Mirrors `BOPAlgo_ShellSplitter::SplitBlock`:
    /// 1. faces that carry a free edge (an edge shared by no other face of the
    ///    block) cannot belong to a closed shell and are dropped, iteratively;
    /// 2. the remaining faces are grouped by edge connectivity into shells —
    ///    two faces belong to the same shell only when their shared edge is
    ///    used by exactly two faces of the block, which keeps closed shells
    ///    that touch along a multi-connected edge apart.
    ///
    /// The shells are stored in `block.loops()`.
    ///
    /// ponytail: the OCCT angle-based face selection at multi-connected edges
    /// is replaced by the "shared edge used by exactly two faces" rule. That
    /// separates shells glued along an edge, but not shells glued along a
    /// face or by non-manifold vertex contact; add the dihedral-angle pass if
    /// that is ever needed.
    pub fn split_block(block: &mut ConnexityBlock) {
        block.change_loops_mut().clear();
        let shapes = block.shapes().to_vec();
        let faces: Vec<Face> = shapes.iter().map(|s| Face(s.clone())).collect();
        let n = faces.len();
        if n == 0 {
            return;
        }
        let face_keys: Vec<HashSet<EKey>> = faces.iter().map(face_edge_keys).collect();
        let mut edge_faces: HashMap<EKey, Vec<usize>> = HashMap::new();
        for (i, keys) in face_keys.iter().enumerate() {
            for &k in keys {
                edge_faces.entry(k).or_default().push(i);
            }
        }
        // 1. Remove faces with free edges, iteratively.
        let mut alive = vec![true; n];
        loop {
            let mut removed_any = false;
            for i in 0..n {
                if !alive[i] {
                    continue;
                }
                let free = face_keys[i].iter().any(|&k| {
                    let cnt = edge_faces
                        .get(&k)
                        .map_or(0, |v| v.iter().filter(|&&j| alive[j]).count());
                    cnt < 2
                });
                if free {
                    alive[i] = false;
                    removed_any = true;
                }
            }
            if !removed_any {
                break;
            }
        }
        // 2. Extract edge-connected closed components.
        let mut visited = vec![false; n];
        for start in 0..n {
            if !alive[start] || visited[start] {
                continue;
            }
            let mut comp: Vec<usize> = Vec::new();
            let mut stack = vec![start];
            visited[start] = true;
            while let Some(i) = stack.pop() {
                comp.push(i);
                for &k in &face_keys[i] {
                    let Some(neigh) = edge_faces.get(&k) else { continue };
                    // Only a bridge edge shared by exactly two alive faces
                    // connects two faces of the same shell.
                    if neigh.iter().filter(|&&j| alive[j]).count() != 2 {
                        continue;
                    }
                    for &j in neigh {
                        if alive[j] && !visited[j] {
                            visited[j] = true;
                            stack.push(j);
                        }
                    }
                }
            }
            if comp.is_empty() {
                continue;
            }
            let shell_faces: Vec<Face> = comp.iter().map(|&i| faces[i].clone()).collect();
            let shell = make_shell_from_faces(&shell_faces);
            block.change_loops_mut().push(shell.0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::primitives::BRepPrimBox;
    use crate::topo_tools_full::{edges_of, faces_of};
    use occt_core::gp::GpPnt;

    #[test]
    fn empty_splitter_yields_no_shells() {
        let mut ss = ShellSplitter::new();
        ss.perform().unwrap();
        assert!(ss.shells().is_empty());
        assert!(ss.blocks().is_empty());
    }

    #[test]
    fn unit_box_is_one_closed_shell() {
        let boxed = unit_box();
        let mut ss = ShellSplitter::new();
        for f in &boxed.faces {
            ss.add_start_element(f.0.clone());
        }
        ss.perform().unwrap();
        assert_eq!(ss.shells().len(), 1, "6 faces of a closed box -> 1 shell");
        let sf = faces_of(&ss.shells()[0]);
        assert_eq!(sf.len(), 6, "shell must contain all 6 faces");
        assert!(ss.shells()[0].closed(), "shell must be flagged closed");
        assert_eq!(ss.blocks().len(), 1);
        assert!(ss.blocks()[0].is_regular(), "a closed box block is regular");
    }

    #[test]
    fn unit_box_from_start_shell() {
        let boxed = unit_box();
        let mut ss = ShellSplitter::new();
        ss.add_start_element(boxed.solid.0.clone());
        ss.perform().unwrap();
        assert_eq!(ss.shells().len(), 1);
        assert_eq!(faces_of(&ss.shells()[0]).len(), 6);
    }

    #[test]
    fn two_separated_cubes_make_two_shells() {
        let a = unit_box();
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(3.0, 0.0, 0.0), &GpPnt::new(4.0, 1.0, 1.0));
        let faces_b = faces_of(&b.solid.0);
        assert_eq!(faces_b.len(), 6);
        let mut ss = ShellSplitter::new();
        for f in &a.faces {
            ss.add_start_element(f.0.clone());
        }
        for f in faces_b {
            ss.add_start_element(f.0.clone());
        }
        ss.perform().unwrap();
        assert_eq!(ss.shells().len(), 2, "two separated cubes -> 2 shells");
        // each shell holds exactly its 6 faces
        for s in ss.shells() {
            assert_eq!(faces_of(s).len(), 6);
            assert!(s.closed());
        }
        assert_eq!(ss.blocks().len(), 2);
    }

    #[test]
    fn open_patch_faces_are_dropped_by_split_block() {
        // A single face is a block whose edges are all free -> no closed shell.
        let boxed = unit_box();
        let f0 = boxed.faces[0].clone();
        let mut cb = ConnexityBlock::new();
        cb.change_shapes_mut().push(f0.0.clone());
        cb.set_regular(false);
        ShellSplitter::split_block(&mut cb);
        assert!(cb.loops().is_empty(), "a lone open face cannot form a closed shell");
    }

    #[test]
    fn split_block_on_box_yields_one_closed_shell() {
        let boxed = unit_box();
        let mut cb = ConnexityBlock::new();
        for f in &boxed.faces {
            cb.change_shapes_mut().push(f.0.clone());
        }
        cb.set_regular(true);
        ShellSplitter::split_block(&mut cb);
        assert_eq!(cb.loops().len(), 1);
        assert_eq!(faces_of(&cb.loops()[0]).len(), 6);
        assert!(cb.loops()[0].closed());
    }

    #[test]
    fn two_boxes_sharing_an_edge_split_into_two_shells() {
        // Cube B sits next to cube A, sharing the vertical edge x=1,y=0.
        // The shared edge has 4 faces around it; split_block must not glue
        // the two closed shells along it.
        let a = unit_box();
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, -1.0, 0.0), &GpPnt::new(2.0, 0.0, 1.0));
        let mut cb = ConnexityBlock::new();
        for f in &a.faces {
            cb.change_shapes_mut().push(f.0.clone());
        }
        for f in faces_of(&b.solid.0) {
            cb.change_shapes_mut().push(f.0.clone());
        }
        cb.set_regular(false);
        ShellSplitter::split_block(&mut cb);
        assert_eq!(cb.loops().len(), 2, "edge-sharing closed shells stay apart");
    }

    #[test]
    fn edges_of_face_counts() {
        let boxed = unit_box();
        for f in &boxed.faces {
            assert_eq!(edges_of(&f.0).len(), 4);
        }
    }
}
