//! `ShapeFix_Shell::FixFaceOrientation` / `GetShells` and the single-shell arm
//! of `ShapeFix_Solid::SolidFromShell`.
//!
//! `ShapeFix_Shell.cxx:1425-1653` and `ShapeFix_Solid.cxx:520-565`, `:655-702`.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_class3d::SolidClassifier;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Shell, TopoShape};

/// One oriented edge occurrence on a face (`TopExp_Explorer(face, EDGE)`).
struct EdgeUse {
    key: usize,
    forward: bool,
    reversed: bool,
}

struct FaceRec {
    shape: TopoShape,
    edges: Vec<EdgeUse>,
}

fn edge_key(e: &TopoShape) -> usize {
    std::sync::Arc::as_ptr(&e.tshape) as usize
}

/// Faces of `shell` with `TopoDS_Iterator` orientation composition
/// (`cumOri = true`). Location stays on the shell: STEP faces carry an
/// identity location, and `B.Add` onto the replacement shell stores the
/// composed orientation as the child's own.
fn composed_faces(shell: &TopoShape) -> Vec<TopoShape> {
    let parent = shell.orientation();
    let stored = shell
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    stored
        .into_iter()
        .filter(|c| c.shape_type() == ShapeType::Face)
        .map(|mut c| {
            c.set_orientation(Orientation::compose(parent, c.orientation()));
            c
        })
        .collect()
}

fn face_edge_uses(face: &TopoShape) -> Vec<EdgeUse> {
    let mut out = Vec::new();
    for w in cumulated_children(face) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        for e in cumulated_children(&w) {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            let o = e.orientation();
            out.push(EdgeUse {
                key: edge_key(&e),
                forward: o == Orientation::Forward,
                reversed: o == Orientation::Reversed,
            });
        }
    }
    out
}

/// Largest connected face group, faces in `GetConnectedFaceGroups` DFS order
/// (`ShapeFix_Shell.cxx:206-297`). Insertion order of `faces` is the indexed
/// map order (`NCollection_IndexedDataMap`, filled from the shell iterator).
fn largest_group(faces: &[FaceRec], edge_faces: &HashMap<usize, Vec<usize>>) -> Vec<usize> {
    let mut visited = vec![false; faces.len()];
    let mut best: Vec<usize> = Vec::new();
    for start in 0..faces.len() {
        if visited[start] {
            continue;
        }
        let mut stack = vec![start];
        visited[start] = true;
        let mut group = Vec::new();
        while let Some(cur) = stack.pop() {
            group.push(cur);
            for e in &faces[cur].edges {
                if let Some(neigh) = edge_faces.get(&e.key) {
                    for &n in neigh {
                        if !visited[n] {
                            visited[n] = true;
                            stack.push(n);
                        }
                    }
                }
            }
        }
        if group.len() > best.len() {
            best = group;
        }
    }
    best
}

/// `GetShells` (`ShapeFix_Shell.cxx:308-648`) for one closed manifold shell.
///
/// Returns the rebuilt face list when at least one face was reversed and every
/// face of the single connected group landed in that shell. Multi-connex edges
/// (`AddMultiConexityFaces`, `:656`), Mobius error faces (`:492-497`,
/// `:1522-1570`) and a split into several shells (`:1620-1630`) are left
/// untouched: those arms replace the shell by a compound, which the STEP
/// solid reader does not expand here.
fn orient_shell_faces(shell: &TopoShape) -> Option<Vec<TopoShape>> {
    let raw = composed_faces(shell);
    if raw.len() < 2 {
        return None;
    }
    let faces: Vec<FaceRec> = raw
        .into_iter()
        .map(|shape| {
            let edges = face_edge_uses(&shape);
            FaceRec { shape, edges }
        })
        .collect();

    let mut edge_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in faces.iter().enumerate() {
        for e in &face.edges {
            let slot = edge_faces.entry(e.key).or_default();
            if !slot.contains(&fi) {
                slot.push(fi);
            }
        }
    }
    // `FixFaceOrientation` passes `isAccountMultiConex = true` (`cxx:1486`).
    // An edge on more than two faces takes `AddMultiConexityFaces`.
    if edge_faces.values().any(|v| v.len() > 2) {
        return None;
    }

    let mut processing = largest_group(&faces, &edge_faces);
    if processing.len() != faces.len() {
        // A second connected group is appended to the unconnected list and
        // becomes its own shell (`cxx:621-637`, `:1594-1604`).
        return None;
    }

    let mut processed: HashMap<usize, (bool, bool)> = HashMap::new();
    let mut shell_faces: Vec<TopoShape> = Vec::new();
    let mut face_idx = 0usize;
    let mut faces_in_shell = 1usize;
    let mut done = false;

    while face_idx < processing.len() {
        let fi = processing[face_idx];
        let mut bad = 0i32;
        let mut good = 0i32;
        let mut temp: HashMap<usize, (bool, bool)> = HashMap::new();

        for e in &faces[fi].edges {
            if let Some(mut pair) = processed.get(&e.key).copied() {
                let (is_direct, is_reversed) = pair;
                if (e.forward && is_direct) || (e.reversed && is_reversed) {
                    bad += 1;
                } else if (e.forward && is_reversed) || (e.reversed && is_direct) {
                    good += 1;
                }
                if is_direct {
                    pair.0 = false;
                } else if is_reversed {
                    pair.1 = false;
                }
                if !pair.0 && !pair.1 {
                    processed.remove(&e.key);
                } else {
                    processed.insert(e.key, pair);
                }
            } else if let Some(t) = temp.get_mut(&e.key) {
                t.0 |= e.forward;
                t.1 |= e.reversed;
            } else {
                temp.insert(e.key, (e.forward, e.reversed));
            }
        }

        if bad == 0 && good == 0 && temp.is_empty() {
            face_idx += 1;
            continue;
        }
        // Mobius: both orientations already accepted (`cxx:492-497`).
        if good != 0 && bad != 0 {
            return None;
        }
        if good != 0 || bad != 0 || faces_in_shell == 1 {
            let mut f = faces[fi].shape.clone();
            if bad != 0 {
                f.reverse();
                for (k, pair) in &temp {
                    processed.insert(*k, (!pair.0, !pair.1));
                }
                done = true;
            } else {
                for (k, pair) in &temp {
                    processed.insert(*k, *pair);
                }
            }
            faces_in_shell += 1;
            shell_faces.push(f);
            processing.remove(face_idx);
            // 1-based loop sets the index to 0 and the `for` increment brings
            // it back to the first remaining face (`cxx:553`, `:568`).
            face_idx = 0;
            continue;
        }
        face_idx += 1;
    }

    if !done || shell_faces.len() != faces.len() {
        return None;
    }
    Some(shell_faces)
}

/// `ShapeFix_Shell::FixFaceOrientation` when it produces exactly one shell
/// (`cxx:1614-1618`, `Context()->Replace` at `:1637`). `None`-equivalent: the
/// shell is returned unchanged.
pub fn fix_shell_face_orientation(mut shell: TopoShape) -> TopoShape {
    if shell.shape_type() != ShapeType::Shell {
        return shell;
    }
    let Some(faces) = orient_shell_faces(&shell) else {
        return shell;
    };
    // Replacement shell is `B.MakeShell` (`cxx:321`), orientation FORWARD.
    shell.set_orientation(Orientation::Forward);
    if let Ok(mut t) = shell.tshape.write() {
        t.children = faces;
    }
    shell
}

/// A non-degenerate edge that meets exactly one face is a free boundary
/// (`ShapeFix_Shell.cxx:1460-1466`). `ShapeAnalysis_FreeBounds` reports the
/// same condition as "open" for the single-shell arm (`ShapeFix_Solid.cxx:528-541`).
fn shell_has_free_edge(shell: &TopoShape) -> bool {
    let mut face_count: HashMap<usize, usize> = HashMap::new();
    for f in composed_faces(shell) {
        let mut seen: HashMap<usize, bool> = HashMap::new();
        for w in cumulated_children(&f) {
            if w.shape_type() != ShapeType::Wire {
                continue;
            }
            for e in cumulated_children(&w) {
                if e.shape_type() != ShapeType::Edge {
                    continue;
                }
                if BRepTool::is_degenerated(&Edge(e.clone())) {
                    continue;
                }
                seen.insert(edge_key(&e), true);
            }
        }
        for k in seen.keys() {
            *face_count.entry(*k).or_insert(0) += 1;
        }
    }
    face_count.values().any(|&n| n == 1)
}

/// Single-shell arm of `ShapeFix_Solid::Perform` (`ShapeFix_Solid.cxx:520-565`)
/// then `SolidFromShell` (`:655-702`): a closed shell whose infinite point
/// classifies `TopAbs_IN` is reversed. `CreateOpenSolidMode` is 0 on the STEP
/// reader (`STEPControl_Controller.cxx:214`), so an open shell is left alone.
/// Several shells go through `CreateSolids` (`:577-640`), which is not this arm.
pub fn orient_single_shell_solid(shell: &mut TopoShape) {
    if shell.shape_type() != ShapeType::Shell {
        return;
    }
    if shell_has_free_edge(shell) {
        return;
    }
    let solid = TopoBuilder::new().make_solid(&[Shell(shell.clone())]);
    let mut sc = SolidClassifier::new();
    sc.load(solid.0);
    sc.perform_infinite_point(occt_core::precision::CONFUSION);
    if sc.state() == FaceState::In {
        shell.reverse();
    }
}
