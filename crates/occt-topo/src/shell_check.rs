//! Shell topology checks — manifold-ness, closedness, Euler characteristic.
//!
//! Source: `BRepCheck_Shell` / `BRepCheck_Analyser` (TKTopAlgo). Edge identity
//! is tracked by `TShape` pointer (`Arc::as_ptr`), so edges shared between
//! faces (the same `Edge` built once and reused) count as one edge used twice.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Shell, TopoShape, Wire};

fn direct_children(s: &TopoShape) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .map(|h| TopoShape::from_handle(h.clone()))
        .collect()
}

/// Depth-first collect every edge in the tree rooted at `s`.
fn collect_edges(s: &TopoShape, out: &mut Vec<TopoShape>) {
    for child in direct_children(s) {
        match child.shape_type() {
            ShapeType::Edge => out.push(child),
            ShapeType::Wire => collect_edges(&child, out),
            _ => {}
        }
    }
}

/// Depth-first collect every face in the tree rooted at `s`.
fn collect_faces(s: &TopoShape, out: &mut Vec<TopoShape>) {
    for child in direct_children(s) {
        match child.shape_type() {
            ShapeType::Face => out.push(child),
            ShapeType::Shell | ShapeType::Solid | ShapeType::Compound | ShapeType::CompSolid => {
                collect_faces(&child, out)
            }
            _ => {}
        }
    }
}

/// Depth-first collect every vertex in the tree rooted at `s`.
fn collect_vertices(s: &TopoShape, out: &mut Vec<TopoShape>) {
    for child in direct_children(s) {
        if child.shape_type() == ShapeType::Vertex {
            out.push(child);
        } else {
            collect_vertices(&child, out);
        }
    }
}

fn dedupe_count(shapes: &[TopoShape]) -> usize {
    let mut seen = HashSet::new();
    for s in shapes {
        seen.insert(Arc::as_ptr(&s.tshape) as usize);
    }
    seen.len()
}

/// Map of edge `TShape` pointer → number of faces referencing it, over the
/// faces found anywhere under `shape`.
pub fn edge_face_usage(shape: &TopoShape) -> HashMap<usize, usize> {
    let mut faces = Vec::new();
    collect_faces(shape, &mut faces);
    let mut usage: HashMap<usize, usize> = HashMap::new();
    for f in faces {
        let mut edges = Vec::new();
        collect_edges(&f, &mut edges);
        for e in edges {
            let k = Arc::as_ptr(&e.tshape) as usize;
            *usage.entry(k).or_insert(0) += 1;
        }
    }
    usage
}

/// Every edge of the shell must be referenced by exactly 2 faces (a closed
/// manifold boundary). Returns `Err` listing the offending edges when not.
pub fn shell_manifold_check(shell: &Shell) -> Result<(), String> {
    let usage = edge_face_usage(&shell.0);
    let mut offending: Vec<(usize, usize)> = usage
        .iter()
        .filter(|(_, &c)| c != 2)
        .map(|(&k, &c)| (k, c))
        .collect();
    offending.sort_by_key(|&(k, _)| k);
    if offending.is_empty() {
        Ok(())
    } else {
        let msg = offending
            .iter()
            .map(|(k, c)| format!("edge {k:#x} used by {c} face(s)"))
            .collect::<Vec<_>>()
            .join("; ");
        Err(msg)
    }
}

/// A shell is closed when every boundary edge is shared by exactly 2 faces.
pub fn shell_is_closed(shell: &Shell) -> bool {
    shell_manifold_check(shell).is_ok()
}

/// Euler characteristic V − E + F over the shell's real (deduplicated)
/// topology.
///
/// A face with `w` boundary wires is a disk with `w−1` holes and contributes
/// `1 − (w−1)` to the face count (a plain disk counts 1, an annulus counts 0,
/// a disk with two holes counts −1), so the sum is topologically correct even
/// when a split face keeps its hole loops (e.g. the box-top ring around a boss
/// cylinder).
pub fn shell_euler_characteristic(shell: &Shell) -> i32 {
    let mut faces = Vec::new();
    collect_faces(&shell.0, &mut faces);
    let mut edges = Vec::new();
    let mut f_count = 0i32;
    for f in &faces {
        collect_edges(f, &mut edges);
        let wires = wires_of_face_count(f);
        f_count += 2 - wires as i32; // disk-with-holes contribution
    }
    let mut verts = Vec::new();
    for e in &edges {
        collect_vertices(e, &mut verts);
    }
    dedupe_count(&verts) as i32 - dedupe_count(&edges) as i32 + f_count
}

/// Number of boundary wires of a face (1 for a plain face; a face built from
/// an outer loop plus hole loops carries one wire per loop).
fn wires_of_face_count(face: &TopoShape) -> usize {
    let mut n = 0usize;
    for child in direct_children(face) {
        if child.shape_type() == ShapeType::Wire {
            n += 1;
        } else if child.shape_type() == ShapeType::Edge {
            // A face whose children are edges directly (no wire) counts one
            // boundary.
            return 1;
        }
    }
    n.max(1)
}

/// Topological invariants of a shell — the boolean-result gate
/// (BOPAlgo / trellis R4 "invariant oracle"). A valid closed solid is
/// manifold (every edge shared by exactly 2 faces) with Euler characteristic
/// V − E + F = 2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShellInvariants {
    /// Every boundary edge is referenced by exactly 2 faces.
    pub closed: bool,
    /// V − E + F over the deduplicated topology.
    pub euler_characteristic: i32,
}

impl ShellInvariants {
    /// True for a closed, genus-0 solid shell (manifold + Euler=2).
    pub fn is_valid_solid(&self) -> bool {
        self.closed && self.euler_characteristic == 2
    }
}

/// Compute the topological invariants of a shell.
pub fn shell_invariants(shell: &Shell) -> ShellInvariants {
    ShellInvariants {
        closed: shell_is_closed(shell),
        euler_characteristic: shell_euler_characteristic(shell),
    }
}

/// Whether the wire's first edge start point equals its last edge end point.
/// Falls back to edge-count parity when the edges carry no evaluable curve.
pub fn wire_is_closed(wire: &Wire) -> bool {
    let edges: Vec<TopoShape> = direct_children(&wire.0)
        .into_iter()
        .filter(|s| s.shape_type() == ShapeType::Edge)
        .collect();
    if edges.is_empty() {
        return false;
    }
    let first = Edge(edges[0].clone());
    let last = Edge(edges[edges.len() - 1].clone());
    if let (Some((start, _)), Some((_, end))) = (BRepTool::edge_vertices(&first), BRepTool::edge_vertices(&last)) {
        return start.distance(&end) <= 1e-9;
    }
    // ponytail: weak parity fallback when edge geometry is unavailable.
    edges.len() % 2 == 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::shape::{Face, Solid, Vertex};
    use crate::tgeometry::GeometryRegistry;
    use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
    use occt_geom::{GeomLine, GeomPlane};
    use std::sync::Arc;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    /// Unit-axis box corners between `lo` and `hi`.
    fn corners(lo: &GpPnt, hi: &GpPnt) -> [GpPnt; 8] {
        let (x0, y0, z0) = (lo.x(), lo.y(), lo.z());
        let (x1, y1, z1) = (hi.x(), hi.y(), hi.z());
        [
            GpPnt::new(x0, y0, z0), // 0
            GpPnt::new(x1, y0, z0), // 1
            GpPnt::new(x1, y1, z0), // 2
            GpPnt::new(x0, y1, z0), // 3
            GpPnt::new(x0, y0, z1), // 4
            GpPnt::new(x1, y0, z1), // 5
            GpPnt::new(x1, y1, z1), // 6
            GpPnt::new(x0, y1, z1), // 7
        ]
    }

    const EDGE_PAIRS: [(usize, usize); 12] = [
        (0, 1), (1, 2), (2, 3), (3, 0), // bottom
        (4, 5), (5, 6), (6, 7), (7, 4), // top
        (0, 4), (1, 5), (2, 6), (3, 7), // vertical
    ];

    /// Face edge-index lists (into EDGE_PAIRS) for each of the 6 box faces:
    /// [-X, +X, -Y, +Y, -Z, +Z].
    const FACE_EDGES: [[usize; 4]; 6] = [
        [8, 7, 11, 3], // -X
        [9, 5, 10, 1], // +X
        [0, 9, 4, 8],  // -Y
        [2, 10, 6, 11], // +Y
        [3, 2, 1, 0],  // -Z
        [4, 5, 6, 7],  // +Z
    ];

    fn plane_face(b: &TopoBuilder, origin: GpPnt, normal: GpDir) -> GpPln {
        let x_dir = if normal.x().abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else if normal.y().abs() > 0.9 {
            GpDir::new(0.0, 0.0, 1.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        let ax3 = GpAx3::new(origin, normal, &x_dir).unwrap();
        GpPln::new(ax3)
    }

    /// Build a box solid (edges share vertices and are reused across faces).
    /// `faces` selects which of the 6 faces to include (0..=5).
    fn build_box(b: &TopoBuilder, lo: &GpPnt, hi: &GpPnt, face_mask: &[bool; 6]) -> Solid {
        let corners = corners(lo, hi);
        let verts: Vec<Vertex> = corners.iter().map(|c| b.make_vertex(*c, 0.0)).collect();
        let edges: Vec<Edge> = EDGE_PAIRS
            .iter()
            .map(|&(i, j)| {
                let dir = GpDir::from_vec(&GpVec::from_pnts(&corners[i], &corners[j])).unwrap();
                let lin = GpLin::from_pnt_dir(corners[i], dir);
                let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, corners[i].distance(&corners[j]));
                b.add(&mut e.0, &verts[i].0);
                b.add(&mut e.0, &verts[j].0);
                e
            })
            .collect();

        let normals: [GpDir; 6] = [
            GpDir::new(-1.0, 0.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            GpDir::new(0.0, -1.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
        ];
        let face_origins: [GpPnt; 6] = [
            corners[0],
            corners[1],
            corners[0],
            corners[3],
            corners[0],
            corners[4],
        ];

        let mut faces: Vec<Face> = Vec::new();
        for (i, include) in face_mask.iter().enumerate() {
            if !include {
                continue;
            }
            let pln = plane_face(b, face_origins[i], normals[i]);
            let wire_edges: Vec<Edge> = FACE_EDGES[i].iter().map(|&e| edges[e].clone()).collect();
            let wire = b.make_wire(&wire_edges);
            faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[wire]));
        }
        let shell = b.make_shell(&faces);
        b.make_solid(&[shell])
    }

    #[test]
    fn closed_six_face_box_is_closed() {
        let b = TopoBuilder::new();
        let solid = build_box(&b, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 2.0, 3.0), &[true; 6]);
        let shell = Shell(TopoShape::from_handle(solid.0.tshape.read().unwrap().children[0].clone()));
        assert!(shell_manifold_check(&shell).is_ok(), "unexpected offending edges");
        assert!(shell_is_closed(&shell));
        // V − E + F = 8 − 12 + 6 = 2
        assert_eq!(shell_euler_characteristic(&shell), 2);
        clear_tree(&solid.0);
    }

    #[test]
    fn open_four_face_shell_is_not_closed() {
        let b = TopoBuilder::new();
        // Only the 4 side faces: missing top/bottom.
        let solid = build_box(
            &b,
            &GpPnt::new(0.0, 0.0, 0.0),
            &GpPnt::new(1.0, 1.0, 1.0),
            &[true, true, true, true, false, false],
        );
        let shell = Shell(TopoShape::from_handle(solid.0.tshape.read().unwrap().children[0].clone()));
        let err = shell_manifold_check(&shell).expect_err("open shell must fail");
        assert!(!shell_is_closed(&shell));
        // The 8 top/bottom edges are referenced once.
        assert!(err.contains("used by 1 face"));
        clear_tree(&solid.0);
    }

    #[test]
    fn edge_face_usage_counts_shared_edges() {
        let b = TopoBuilder::new();
        let solid = build_box(&b, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0), &[true; 6]);
        let usage = edge_face_usage(&solid.0);
        assert_eq!(usage.len(), 12);
        assert!(usage.values().all(|&c| c == 2), "usage {:?}", usage);
        clear_tree(&solid.0);
    }

    #[test]
    fn closed_wire_is_closed() {
        let b = TopoBuilder::new();
        let edges = [
            b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0)),
            b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0)),
            b.make_edge_segment(&GpPnt::new(1.0, 1.0, 0.0), &GpPnt::new(0.0, 1.0, 0.0)),
            b.make_edge_segment(&GpPnt::new(0.0, 1.0, 0.0), &GpPnt::new(0.0, 0.0, 0.0)),
        ];
        let wire = b.make_wire(&edges);
        assert!(wire_is_closed(&wire));
        clear_tree(&wire.0);
    }

    #[test]
    fn open_wire_is_not_closed() {
        let b = TopoBuilder::new();
        let edges = [
            b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0)),
            b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0)),
        ];
        let wire = b.make_wire(&edges);
        assert!(!wire_is_closed(&wire));
        clear_tree(&wire.0);
    }
}
