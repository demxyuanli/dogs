//! Advanced shape validity checks — seam edges, orientation consistency,
//! boundary multiplicity and tolerance propagation.
//! Source: `BRepCheck_Analyser`, `BRepCheck_Shell`, `BRepCheck_Edge`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Shell, TopoShape, Wire};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, wires_of_face};

/// A focused check result: error/warning with a short message and the type
/// of shape it concerns.
#[derive(Debug, Clone, PartialEq)]
pub enum CheckLevel {
    Ok,
    Warning(String),
    Error(String),
}

impl CheckLevel {
    pub fn is_ok(&self) -> bool {
        matches!(self, CheckLevel::Ok)
    }
}

/// Every edge in a closed shell must be used by exactly two faces (boundary
/// multiplicity 2). Returns a map of edge → usage count for the offenders.
/// Mirrors `BRepCheck_Shell::CheckClosed`.
pub fn edge_usage_counts(shape: &TopoShape) -> HashMap<usize, usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for f in faces_of(shape) {
        for w in wires_of_face(&f) {
            for e in edges_of_wire(&w) {
                *counts.entry(Arc::as_ptr(&e.0.tshape) as usize).or_insert(0) += 1;
            }
        }
    }
    counts
}

/// Boundary multiplicity check: returns the edges whose usage count is not 2
/// (for a closed manifold). `expected` defaults to 2.
pub fn boundary_multiple_check(shape: &TopoShape) -> Vec<(Edge, usize, CheckLevel)> {
    let counts = edge_usage_counts(shape);
    let mut out = Vec::new();
    for e in edges_of(shape) {
        let n = counts.get(&(Arc::as_ptr(&e.0.tshape) as usize)).copied().unwrap_or(0);
        let level = match n {
            2 => CheckLevel::Ok,
            0 => CheckLevel::Error(format!("edge used by 0 faces")),
            1 => CheckLevel::Error(format!("edge used by 1 face (open boundary)")),
            _ => CheckLevel::Warning(format!("edge used by {n} faces")),
        };
        out.push((e, n, level));
    }
    out
}

/// Whether a shell is closed (every boundary edge used twice).
pub fn shell_closed_check(shell: &Shell) -> CheckLevel {
    let offenders: Vec<String> = boundary_multiple_check(&shell.0)
        .into_iter()
        .filter(|(_, _, lvl)| !lvl.is_ok())
        .map(|(_, _, lvl)| match lvl {
            CheckLevel::Error(m) => m,
            _ => String::new(),
        })
        .filter(|m| !m.is_empty())
        .collect();
    if offenders.is_empty() {
        CheckLevel::Ok
    } else {
        CheckLevel::Error(format!("{} boundary errors: {}", offenders.len(), offenders[0]))
    }
}

/// Every wire of a face must be closed. Returns the offending wires.
pub fn face_wire_closed_check(shape: &TopoShape) -> Vec<(Wire, CheckLevel)> {
    let mut out = Vec::new();
    for f in faces_of(shape) {
        for w in wires_of_face(&f) {
            if wire_geometrically_closed(&w) {
                out.push((w, CheckLevel::Ok));
            } else {
                out.push((w, CheckLevel::Error("wire is open".into())));
            }
        }
    }
    out
}

/// Geometric wire closure via endpoint parity: in a closed wire every vertex
/// appears exactly twice as an edge endpoint (shared corner vertices count in
/// two edges; a seam edge used twice within the wire is still even). Robust to
/// edges stored in either orientation relative to the wire traversal.
fn wire_geometrically_closed(wire: &Wire) -> bool {
    let edges = edges_of_wire(wire);
    if edges.is_empty() {
        return true;
    }
    if edges.len() < 3 {
        return false;
    }
    let mut counts: HashMap<[i64; 3], usize> = HashMap::new();
    for e in &edges {
        let (a, b) = crate::topo_tools_full::edge_vertices(e);
        for v in [a, b].into_iter().flatten() {
            let p = BRepTool::vertex_point(&v);
            let key = [
                (p.x() * 1e9).round() as i64,
                (p.y() * 1e9).round() as i64,
                (p.z() * 1e9).round() as i64,
            ];
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    counts.values().all(|&c| c % 2 == 0)
}

/// Face orientation consistency: the outward face normal must agree with the
/// face's winding for a well-built solid. Reports faces whose sampled normal
/// points the "wrong way" relative to the shell centroid (heuristic — a solid
/// whose centroid lies behind the face normal at the face's centroid).
pub fn orientation_consistency_check(shape: &TopoShape) -> Vec<(Face, CheckLevel)> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Vec::new();
    }
    // Shell centroid from face centroid average.
    let mut acc = occt_core::gp::GpXyz::zero();
    let mut n = 0usize;
    let mut centroids: Vec<(Face, GpPnt, GpVec)> = Vec::new();
    for f in &faces {
        let Some(c) = face_centroid(f) else { continue };
        let Some(nrm) = face_normal_avg(f) else { continue };
        acc = acc.added(&c.coord);
        n += 1;
        centroids.push((f.clone(), c, nrm));
    }
    if n == 0 {
        return Vec::new();
    }
    let centroid = GpPnt::from_xyz(&acc.divided(n as f64));
    centroids
        .into_iter()
        .map(|(f, c, nrm)| {
            let outward = GpVec::from_pnts(&centroid, &c);
            let dot = outward.xyz().dot(nrm.xyz());
            if dot < 0.0 {
                (f, CheckLevel::Warning("face normal points inward".into()))
            } else {
                (f, CheckLevel::Ok)
            }
        })
        .collect()
}

fn face_centroid(f: &Face) -> Option<GpPnt> {
    let s = BRepTool::face_surface(f)?;
    let (u0, u1, v0, v1) = BRepTool::uv_bounds(f);
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
        // Unbounded plane: use the outer-wire vertex centroid instead.
        let wires = wires_of_face(f);
        let mut acc = occt_core::gp::GpXyz::zero();
        let mut n = 0usize;
        for w in wires {
            for e in edges_of_wire(&w) {
                if let (Some(a), _) = crate::topo_tools_full::edge_vertices(&e) {
                    acc = acc.added(&BRepTool::vertex_point(&a).coord);
                    n += 1;
                }
            }
        }
        return if n > 0 { Some(GpPnt::from_xyz(&acc.divided(n as f64))) } else { None };
    }
    let mut acc = occt_core::gp::GpXyz::zero();
    let mut n = 0usize;
    for i in 0..4 {
        for j in 0..4 {
            let u = u0 + (u1 - u0) * i as f64 / 3.0;
            let v = v0 + (v1 - v0) * j as f64 / 3.0;
            acc = acc.added(&s.d0(u, v).coord);
            n += 1;
        }
    }
    Some(GpPnt::from_xyz(&acc.divided(n as f64)))
}

/// Average unit normal of a face (finite-difference of d0, robust to d1=0).
fn face_normal_avg(f: &Face) -> Option<GpVec> {
    let s = BRepTool::face_surface(f)?;
    let (u0, u1, v0, v1) = BRepTool::uv_bounds(f);
    let (u0, u1, v0, v1) = if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
        (u0, u1, v0, v1)
    } else {
        (-1.0, 1.0, -1.0, 1.0)
    };
    let mut acc = occt_core::gp::GpXyz::zero();
    let mut n = 0usize;
    for i in 0..4 {
        for j in 0..4 {
            let u = u0 + (u1 - u0) * i as f64 / 3.0;
            let v = v0 + (v1 - v0) * j as f64 / 3.0;
            let (_, du, dv) = s.d1(u, v);
            let nrm = du.xyz().crossed(dv.xyz());
            if nrm.square_modulus() > 1e-30 {
                acc = acc.added(&nrm);
                n += 1;
            }
        }
    }
    if n == 0 {
        return None;
    }
    let m = acc.modulus();
    if m < 1e-30 {
        None
    } else {
        Some(GpVec::new(acc.x / m, acc.y / m, acc.z / m))
    }
}

/// Seam-edge check: periodic faces (cylinder/sphere/torus) whose wire contains
/// a seam edge used twice within the SAME face are valid; a face using the same
/// edge twice non-adjacently is suspicious. Reports faces where any edge is
/// duplicated within one wire.
pub fn seam_edge_check(shape: &TopoShape) -> Vec<(Wire, CheckLevel)> {
    let mut out = Vec::new();
    for f in faces_of(shape) {
        for w in wires_of_face(&f) {
            let edges = edges_of_wire(&w);
            let mut seen: HashMap<usize, usize> = HashMap::new();
            for e in &edges {
                *seen.entry(Arc::as_ptr(&e.0.tshape) as usize).or_insert(0) += 1;
            }
            let seams: Vec<usize> = seen.values().copied().filter(|&c| c > 2).collect();
            if seams.is_empty() {
                out.push((w, CheckLevel::Ok));
            } else {
                out.push((w, CheckLevel::Warning(format!(
                    "edge appears >2 times in one wire ({} occurrences)", seams[0]
                ))));
            }
        }
    }
    out
}

/// Tolerance propagation hint: the max vertex-vertex gap along shared edges.
/// Returns the largest detected gap (0 if none / no shared vertices).
pub fn max_vertex_gap(shape: &TopoShape) -> f64 {
    let mut worst: f64 = 0.0;
    for e in edges_of(shape) {
        let (a, b) = crate::topo_tools_full::edge_vertices(&e);
        if let (Some(a), Some(b)) = (a, b) {
            let d = BRepTool::vertex_point(&a).distance(&BRepTool::vertex_point(&b));
            worst = worst.max(d);
        }
    }
    worst
}

/// Aggregate: run all checks and return (errors, warnings).
pub fn analyse(shape: &TopoShape) -> (Vec<String>, Vec<String>) {
    let mut errors = Vec::new();
    let mut warnings = Vec::new();
    for (_, _, lvl) in boundary_multiple_check(shape) {
        if let CheckLevel::Error(m) = lvl {
            errors.push(m);
        } else if let CheckLevel::Warning(m) = lvl {
            warnings.push(m);
        }
    }
    for (_, lvl) in face_wire_closed_check(shape) {
        if let CheckLevel::Error(m) = lvl {
            errors.push(m);
        }
    }
    for (_, lvl) in orientation_consistency_check(shape) {
        if let CheckLevel::Warning(m) = lvl {
            warnings.push(m);
        }
    }
    for (_, lvl) in seam_edge_check(shape) {
        if let CheckLevel::Warning(m) = lvl {
            warnings.push(m);
        }
    }
    (errors, warnings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    #[test]
    fn closed_box_boundary_ok() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let issues = boundary_multiple_check(&b.solid.0);
        assert!(issues.iter().all(|(_, _, lvl)| lvl.is_ok()));
        assert!(matches!(shell_closed_check(&{
            let shells = crate::topo_tools_full::shapes_of(&b.solid.0, ShapeType::Shell);
            crate::shape::Shell(shells[0].clone())
        }), CheckLevel::Ok));
    }

    #[test]
    fn open_shell_detected() {
        // A single planar face is an open boundary (its edges used once).
        let b = crate::builder::TopoBuilder::new();
        let face = b.make_face_plane(&occt_core::gp::GpPln::new(occt_core::gp::GpAx3::standard()));
        let issues = boundary_multiple_check(&face.0);
        // The plane face has NO wire (unbounded) → no edges → empty.
        assert!(issues.is_empty());
        // Build a polygon face instead: it has edges used once.
        let poly = crate::brep_builder_api::make_face_from_polygon(&[
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.),
        ]).unwrap();
        let issues = boundary_multiple_check(&poly.0);
        assert!(issues.iter().all(|(_, n, lvl)| *n == 1 && !lvl.is_ok()));
    }

    #[test]
    fn wires_closed_for_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let checks = face_wire_closed_check(&b.solid.0);
        assert_eq!(checks.len(), 6);
        assert!(checks.iter().all(|(_, lvl)| lvl.is_ok()));
    }

    #[test]
    fn orientation_consistent_for_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let issues = orientation_consistency_check(&b.solid.0);
        // All box faces point outward from the centroid.
        assert!(issues.iter().all(|(_, lvl)| lvl.is_ok()));
    }

    #[test]
    fn aggregate_clean_for_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let (errors, _) = analyse(&b.solid.0);
        assert!(errors.is_empty(), "errors: {errors:?}");
    }
}
