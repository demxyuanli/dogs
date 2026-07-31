//! Topology validity checks.
//!
//! Approximate port of `BRepCheck_Analyser` (TKTopAlgo) together with the
//! structural helpers from `ShapeAnalysis`. This module verifies that a shape
//! has well-formed geometry (vertices with points, edges with curves and
//! matching endpoints, faces with closed wires) and reports topological
//! invariants such as counts, the Euler–Poincaré number and watertightness.
//!
//! Checks are structural rather than numerical: they inspect the geometry
//! registered in the side-table and the shape tree, but do not attempt a full
//! self-intersection or validity analysis.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topexp::Explorer;

/// Outcome of a shape analysis.
#[derive(Debug, Clone, Default)]
pub struct ShapeReport {
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

impl ShapeReport {
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }

    pub fn summary(&self) -> String {
        if self.is_valid() {
            format!("valid ({} warnings)", self.warnings.len())
        } else {
            format!(
                "{} errors, {} warnings",
                self.errors.len(),
                self.warnings.len()
            )
        }
    }
}

/// Run every check over `shape`, collecting issues into a report.
pub fn analyse(shape: &TopoShape) -> ShapeReport {
    let mut report = ShapeReport::default();
    for msg in check_vertices(shape) {
        if let Some(w) = msg.strip_prefix("warning:") {
            report.warnings.push(w.trim().to_string());
        } else {
            report.errors.push(msg);
        }
    }
    for msg in check_edges(shape) {
        report.errors.push(msg);
    }
    for msg in check_faces(shape) {
        report.errors.push(msg);
    }
    report
}

/// Every vertex must carry a registered point; duplicate locations are warned.
pub fn check_vertices(shape: &TopoShape) -> Vec<String> {
    let mut issues = Vec::new();
    let mut pts: Vec<(usize, GpPnt)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let vs = ex.current().clone();
        ex.next();
        let key = Arc::as_ptr(&vs.tshape) as usize;
        if !seen.insert(key) {
            continue;
        }
        if let Some(v) = Vertex::wrap(vs) {
            match GeometryRegistry::global().vertex_geom(&v.0) {
                None => issues.push(format!("vertex {key:#x}: no registered point")),
                Some(g) => {
                    for &(k2, ref p2) in &pts {
                        if g.point.distance(p2) < 1e-7 {
                            issues.push(format!(
                                "warning: vertex {key:#x} duplicates vertex {k2:#x} at {}",
                                fmt_pnt(&g.point)
                            ));
                        }
                    }
                    pts.push((key, g.point));
                }
            }
        }
    }
    issues
}

/// Every edge must have a curve, a finite range (or an intentionally unbounded
/// curve), and endpoints that agree with its child vertex points.
pub fn check_edges(shape: &TopoShape) -> Vec<String> {
    let mut issues = Vec::new();
    let mut idx = 0usize;
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Edge);
    while ex.more() {
        let es = ex.current().clone();
        ex.next();
        if !seen.insert(Arc::as_ptr(&es.tshape) as usize) {
            continue;
        }
        if let Some(e) = Edge::wrap(es) {
            let (first, last) = BRepTool::edge_parameters(&e);
            match BRepTool::edge_curve(&e) {
                None => issues.push(format!("edge {idx}: missing curve")),
                Some(curve) => {
                    let range_infinite = !first.is_finite() || !last.is_finite();
                    if range_infinite {
                        // Unbounded ranges are legal only for curves whose own
                        // parameter range is unbounded (infinite lines, etc.).
                        let curve_unbounded =
                            !curve.first_parameter().is_finite() || !curve.last_parameter().is_finite();
                        if !curve_unbounded {
                            issues.push(format!(
                                "edge {idx}: non-finite parameter range [{first}, {last}] on a bounded curve"
                            ));
                        }
                    }
                    // Endpoints must agree with the child vertices, if any.
                    let kids: Vec<TopoShape> = e
                        .0
                        .tshape
                        .read()
                        .unwrap()
                        .children
                        .iter()
                        .filter(|h| h.read().unwrap().shape_type() == ShapeType::Vertex)
                        .map(|h| TopoShape::from_handle(h.clone()))
                        .collect();
                    if !kids.is_empty() {
                        let p_first = BRepTool::vertex_point(
                            &Vertex::wrap(kids[0].clone()).expect("vertex child"),
                        );
                        let p_last = BRepTool::vertex_point(
                            &Vertex::wrap(kids[kids.len() - 1].clone()).expect("vertex child"),
                        );
                        let tol = BRepTool::edge_tolerance(&e).max(1e-7);
                        let c_first = curve.d0(first);
                        let c_last = curve.d0(last);
                        if c_first.distance(&p_first) > tol {
                            issues.push(format!(
                                "edge {idx}: start {} does not match vertex {} (tol {tol})",
                                fmt_pnt(&c_first),
                                fmt_pnt(&p_first)
                            ));
                        }
                        if c_last.distance(&p_last) > tol {
                            issues.push(format!(
                                "edge {idx}: end {} does not match vertex {} (tol {tol})",
                                fmt_pnt(&c_last),
                                fmt_pnt(&p_last)
                            ));
                        }
                    }
                }
            }
        }
        idx += 1;
    }
    issues
}

/// Every face must have a surface; every wire in it must be closed.
pub fn check_faces(shape: &TopoShape) -> Vec<String> {
    let mut issues = Vec::new();
    let mut idx = 0usize;
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let fs = ex.current().clone();
        ex.next();
        if !seen.insert(Arc::as_ptr(&fs.tshape) as usize) {
            continue;
        }
        if let Some(face) = Face::wrap(fs) {
            if BRepTool::face_surface(&face).is_none() {
                issues.push(format!("face {idx}: missing surface"));
            }
            let mut w = 0usize;
            let kids = face.0.tshape.read().unwrap().children.clone();
            for h in kids {
                if h.read().unwrap().shape_type() != ShapeType::Wire {
                    continue;
                }
                if let Some(wire) = Wire::wrap(TopoShape::from_handle(h)) {
                    if !check_wire_closed(&wire) {
                        issues.push(format!("face {idx}: wire {w} is not closed"));
                    }
                    w += 1;
                }
            }
        }
        idx += 1;
    }
    issues
}

/// Whether a wire's edges form one closed loop.
///
/// Robust to edge orientation: the check is order-independent. Each edge
/// contributes its two endpoint points; the wire is closed when every endpoint
/// is shared by exactly two edge-ends (degree 2 at every vertex) and the
/// edges form a single connected component.
pub fn check_wire_closed(wire: &Wire) -> bool {
    let edges = wire_edges(wire);
    if edges.is_empty() {
        return false;
    }
    let mut endpoints: Vec<GpPnt> = Vec::with_capacity(edges.len() * 2);
    for e in &edges {
        match edge_endpoint_points(e) {
            Some((a, b)) => {
                endpoints.push(a);
                endpoints.push(b);
            }
            None => return false,
        }
    }
    let tol = 1e-7;
    let n = endpoints.len();
    let mut matched = vec![false; n];
    for i in 0..n {
        if matched[i] {
            continue;
        }
        let mut found = false;
        for j in (i + 1)..n {
            if !matched[j] && endpoints[i].distance(&endpoints[j]) <= tol {
                matched[i] = true;
                matched[j] = true;
                found = true;
                break;
            }
        }
        if !found {
            return false;
        }
    }

    // Connectivity: every edge in one component via shared endpoints.
    let mut parent: Vec<usize> = (0..edges.len()).collect();
    for i in 0..edges.len() {
        let (a_i, b_i) = (endpoints[2 * i], endpoints[2 * i + 1]);
        for j in (i + 1)..edges.len() {
            let (a_j, b_j) = (endpoints[2 * j], endpoints[2 * j + 1]);
            if a_i.distance(&a_j) <= tol
                || a_i.distance(&b_j) <= tol
                || b_i.distance(&a_j) <= tol
                || b_i.distance(&b_j) <= tol
            {
                union(&mut parent, i, j);
            }
        }
    }
    let root = find(&mut parent, 0);
    (0..edges.len()).all(|i| find(&mut parent, i) == root)
}

fn find(parent: &mut Vec<usize>, mut i: usize) -> usize {
    let mut r = i;
    while parent[r] != r {
        r = parent[r];
    }
    while parent[i] != r {
        let next = parent[i];
        parent[i] = r;
        i = next;
    }
    r
}

fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb {
        parent[ra] = rb;
    }
}

fn wire_edges(wire: &Wire) -> Vec<Edge> {
    let mut out = Vec::new();
    let kids = wire.0.tshape.read().unwrap().children.clone();
    for h in kids {
        if h.read().unwrap().shape_type() == ShapeType::Edge {
            if let Some(e) = Edge::wrap(TopoShape::from_handle(h)) {
                out.push(e);
            }
        }
    }
    out
}

/// The two endpoint points of an edge: real child vertices when present,
/// otherwise the curve evaluated at the parameter bounds.
fn edge_endpoint_points(e: &Edge) -> Option<(GpPnt, GpPnt)> {
    let kids: Vec<TopoShape> = e
        .0
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.read().unwrap().shape_type() == ShapeType::Vertex)
        .map(|h| TopoShape::from_handle(h.clone()))
        .collect();
    if kids.len() >= 2 {
        let a = BRepTool::vertex_point(&Vertex::wrap(kids[0].clone()).expect("vertex child"));
        let b = BRepTool::vertex_point(
            &Vertex::wrap(kids[kids.len() - 1].clone()).expect("vertex child"),
        );
        Some((a, b))
    } else {
        BRepTool::edge_vertices(e)
    }
}

fn fmt_pnt(p: &GpPnt) -> String {
    format!("({:.6}, {:.6}, {:.6})", p.x(), p.y(), p.z())
}

// ---------------------------------------------------------------------------
// Topological invariants.
// ---------------------------------------------------------------------------

/// Count of distinct vertices in the shape tree.
pub fn count_vertices(shape: &TopoShape) -> usize {
    count_distinct(shape, ShapeType::Vertex)
}

/// Count of distinct edges.
pub fn count_edges(shape: &TopoShape) -> usize {
    count_distinct(shape, ShapeType::Edge)
}

/// Count of distinct faces.
pub fn count_faces(shape: &TopoShape) -> usize {
    count_distinct(shape, ShapeType::Face)
}

/// Count of distinct wires.
pub fn count_wires(shape: &TopoShape) -> usize {
    count_distinct(shape, ShapeType::Wire)
}

fn count_distinct(shape: &TopoShape, target: ShapeType) -> usize {
    let mut seen = std::collections::HashSet::new();
    let mut n = 0;
    let mut ex = Explorer::new(shape, target);
    while ex.more() {
        let s = ex.current();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            n += 1;
        }
        ex.next();
    }
    n
}

/// Count of distinct shapes of each type in the shape tree, including the
/// root shape itself.
pub fn shape_type_breakdown(shape: &TopoShape) -> HashMap<ShapeType, usize> {
    let mut m = HashMap::new();
    let mut seen = std::collections::HashSet::new();
    *m.entry(shape.shape_type()).or_insert(0) += 1;
    seen.insert(Arc::as_ptr(&shape.tshape) as usize);
    let mut ex = Explorer::new(shape, ShapeType::Shape);
    while ex.more() {
        let s = ex.current();
        if seen.insert(Arc::as_ptr(&s.tshape) as usize) {
            *m.entry(s.shape_type()).or_insert(0) += 1;
        }
        ex.next();
    }
    m
}

/// Euler–Poincaré characteristic `V − E + F` of the whole shape.
///
/// For a closed genus-0 solid this equals 2. (Note: `C` for cells is not
/// subtracted; the formula here is the classic boundary invariant `V − E + F`.)
pub fn euler_poincare(shape: &TopoShape) -> i32 {
    count_vertices(shape) as i32 - count_edges(shape) as i32 + count_faces(shape) as i32
}

/// Whether every edge of the shape is referenced by exactly two faces.
///
/// A closed (watertight) shell has each boundary edge shared by exactly two
/// faces. Free edges (one face) or non-manifold edges (three or more faces)
/// make the shape leaky.
pub fn is_watertight(shape: &TopoShape) -> bool {
    let mut face_edges: Vec<Vec<TopoShape>> = Vec::new();
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let fs = ex.current().clone();
        ex.next();
        if let Some(face) = Face::wrap(fs) {
            let mut edges = Vec::new();
            let face_kids = face.0.tshape.read().unwrap().children.clone();
            for h in face_kids {
                if h.read().unwrap().shape_type() != ShapeType::Wire {
                    continue;
                }
                let wire = TopoShape::from_handle(h);
                let wire_kids = wire.tshape.read().unwrap().children.clone();
                for eh in wire_kids {
                    if eh.read().unwrap().shape_type() == ShapeType::Edge {
                        edges.push(TopoShape::from_handle(eh));
                    }
                }
            }
            face_edges.push(edges);
        }
    }
    if face_edges.is_empty() {
        return false;
    }

    let mut refs: HashMap<usize, usize> = HashMap::new();
    for fe in &face_edges {
        let mut seen_in_face = std::collections::HashSet::new();
        for e in fe {
            let key = Arc::as_ptr(&e.tshape) as usize;
            if seen_in_face.insert(key) {
                *refs.entry(key).or_insert(0) += 1;
            }
        }
    }
    refs.values().all(|&n| n == 2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::brep_extrema::test_box;
    use occt_core::gp::GpPnt;

    fn rect_wire_closed() -> Wire {
        let b = TopoBuilder::new();
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let d = GpPnt::new(1.0, 0.0, 0.0);
        let c = GpPnt::new(1.0, 1.0, 0.0);
        let b4 = GpPnt::new(0.0, 1.0, 0.0);
        let edges = [
            b.make_edge_segment(&a, &d),
            b.make_edge_segment(&d, &c),
            b.make_edge_segment(&c, &b4),
            b.make_edge_segment(&b4, &a),
        ];
        b.make_wire(&edges)
    }

    fn open_two_edge_wire() -> Wire {
        let b = TopoBuilder::new();
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let d = GpPnt::new(1.0, 0.0, 0.0);
        let c = GpPnt::new(1.0, 1.0, 0.0);
        let edges = [b.make_edge_segment(&a, &d), b.make_edge_segment(&d, &c)];
        b.make_wire(&edges)
    }

    #[test]
    fn valid_box_has_no_errors() {
        let b = test_box::unit_box();
        let report = analyse(&b.solid.0);
        assert!(report.errors.is_empty(), "errors: {:?}", report.errors);
        assert!(report.is_valid());
    }

    #[test]
    fn wire_closedness() {
        assert!(check_wire_closed(&rect_wire_closed()));
        assert!(!check_wire_closed(&open_two_edge_wire()));
    }

    #[test]
    fn box_counts() {
        let b = test_box::unit_box();
        assert_eq!(count_vertices(&b.solid.0), 8);
        assert_eq!(count_edges(&b.solid.0), 12);
        assert_eq!(count_faces(&b.solid.0), 6);
        assert_eq!(count_wires(&b.solid.0), 6);
    }

    #[test]
    fn box_euler_and_watertight() {
        let b = test_box::unit_box();
        assert_eq!(euler_poincare(&b.solid.0), 2);
        assert!(is_watertight(&b.solid.0));
    }

    #[test]
    fn breakdown_contains_expected_types() {
        let b = test_box::unit_box();
        let m = shape_type_breakdown(&b.solid.0);
        assert_eq!(m.get(&ShapeType::Vertex), Some(&8));
        assert_eq!(m.get(&ShapeType::Edge), Some(&12));
        assert_eq!(m.get(&ShapeType::Face), Some(&6));
        assert_eq!(m.get(&ShapeType::Solid), Some(&1));
    }

    #[test]
    fn report_summary_mentions_errors() {
        let b = test_box::unit_box();
        let report = analyse(&b.solid.0);
        assert!(report.summary().contains("valid"));
    }
}
