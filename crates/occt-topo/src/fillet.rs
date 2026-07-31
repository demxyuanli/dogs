//! Corner fillet on straight-edge polygon wires — a simplified port of
//! `BRepFilletAPI_MakeFillet` limited to 2D corners (two straight edges meeting
//! at a vertex).
//!
//! A corner at `p2` formed by the segments `p1→p2` and `p2→p3` is rounded by a
//! circular arc of the requested `radius` tangent to both segments. The fillet
//! replaces the sharp corner with two trimmed line edges and one circular-arc
//! edge.

use occt_core::gp::{GpAx2, GpDir, GpPnt, GpVec};

use crate::brep_surface::face_plane;
use crate::builder::TopoBuilder;
use crate::primitives::BRepPrimBox;
use crate::shape::{Edge, Wire};
use crate::topo_tools_full::{edge_vertices, edges_of_wire, faces_of, vertex_position, wires_of_face};

/// The two tangency points where a fillet of `radius` touches the corner at
/// `p2` — the first on `p1→p2`, the second on `p2→p3`. Returns `None` when
/// `radius` is non-positive or either edge is shorter than it (the tangency
/// would fall past the far vertex).
pub fn tangent_points(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt, radius: f64) -> Option<(GpPnt, GpPnt)> {
    let d1 = GpVec::from_pnts(p2, p1); // toward p1
    let d2 = GpVec::from_pnts(p2, p3); // toward p3
    let l1 = d1.magnitude();
    let l2 = d2.magnitude();
    if radius <= 0.0 || l1 < radius || l2 < radius {
        return None;
    }
    let t1 = p2.translated_vec(&d1.multiplied_scalar(radius / l1));
    let t2 = p2.translated_vec(&d2.multiplied_scalar(radius / l2));
    Some((t1, t2))
}

/// The fillet arc's center: the intersection of the two offset lines, each edge
/// offset by `radius` toward the interior of the corner.
///
/// `cross = (p1−p2) × (p3−p2)` is the corner's plane normal. The interior side
/// is picked via the cross product: `cross × (p1−p2)` is the in-wedge unit
/// normal of the first edge and `(p3−p2) × cross` the in-wedge unit normal of
/// the second, for either sign of `cross`. ponytail: reflex corners (interior
/// angle > 180°) are not distinguished; the interior is always taken as the
/// wedge between the two edges.
pub fn arc_center(
    p1: &GpPnt,
    p2: &GpPnt,
    p3: &GpPnt,
    t1: &GpPnt,
    t2: &GpPnt,
    radius: f64,
) -> Option<GpPnt> {
    let v1 = GpVec::from_pnts(p2, p1);
    let v2 = GpVec::from_pnts(p2, p3);
    let cross = v1.crossed(&v2);
    let m = cross.magnitude();
    if m < 1e-12 || radius <= 0.0 {
        return None;
    }
    // Unit in-wedge normals: |cross × v1| = |cross|·|v1|, |v2 × cross| = |v2|·|cross|.
    let n1 = cross.crossed(&v1).divided(m * v1.magnitude());
    let n2 = v2.crossed(&cross).divided(m * v2.magnitude());
    // Intersection of the offset lines t1 + s·n1 and t2 + u·n2 (both lie in the
    // corner plane, so n1 × n2 is the plane normal).
    let w = GpVec::from_pnts(t1, t2);
    let axis = n1.crossed(&n2);
    let denom = axis.square_magnitude();
    if denom < 1e-24 {
        return None;
    }
    let s = w.crossed(&n2).dot(&axis) / denom;
    Some(t1.translated_vec(&n1.multiplied_scalar(s)))
}

/// One edge of a filleted corner chain, annotated with its tangency points.
/// `arc` is `Some` for the circular-arc edge and `None` for the two trimmed
/// line edges.
#[derive(Debug, Clone)]
pub struct FilletEdge {
    pub edge: Edge,
    pub tangency_a: GpPnt,
    pub tangency_b: GpPnt,
    pub arc: Option<Edge>,
}

/// Build the three edges replacing a sharp corner: trimmed line `p1→t1`,
/// circular arc `t1→t2`, trimmed line `t2→p3`.
fn corner_fillet(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt, radius: f64) -> Result<[FilletEdge; 3], String> {
    let (t1, t2) = tangent_points(p1, p2, p3, radius)
        .ok_or_else(|| "fillet: an edge is shorter than the radius".to_string())?;
    let center = arc_center(p1, p2, p3, &t1, &t2, radius)
        .ok_or_else(|| "fillet: cannot compute arc center".to_string())?;
    let cross = GpVec::from_pnts(p2, p1).crossed(&GpVec::from_pnts(p2, p3));
    let b = TopoBuilder::new();
    let seg1 = b.make_edge_segment(p1, &t1);
    let arc = arc_edge(&b, &center, &cross, &t1, &t2, radius)?;
    let seg2 = b.make_edge_segment(&t2, p3);
    Ok([
        FilletEdge { edge: seg1, tangency_a: *p1, tangency_b: t1, arc: None },
        FilletEdge { edge: arc.clone(), tangency_a: t1, tangency_b: t2, arc: Some(arc) },
        FilletEdge { edge: seg2, tangency_a: t2, tangency_b: *p3, arc: None },
    ])
}

/// Circular-arc edge from `t1` to `t2` through the rounded corner, with the
/// given center and radius. `cross` is the corner's plane normal; the arc's
/// X direction is the radius to `t1`, so `t1` sits at parameter 0 and `t2` at
/// the signed interior angle between the two tangency radii.
fn arc_edge(
    b: &TopoBuilder,
    center: &GpPnt,
    cross: &GpVec,
    t1: &GpPnt,
    t2: &GpPnt,
    radius: f64,
) -> Result<Edge, String> {
    let z = GpDir::from_vec(cross).map_err(|e| e.to_string())?;
    let x = GpDir::from_vec(&GpVec::from_pnts(center, t1)).map_err(|e| e.to_string())?;
    let ax2 = GpAx2::new(*center, z, x).map_err(|e| e.to_string())?;
    let xd = *ax2.x_direction();
    let yd = *ax2.y_direction();
    let va = GpVec::from_pnts(center, t1);
    let vb = GpVec::from_pnts(center, t2);
    let a1 = va.coord.dot(yd.xyz()).atan2(va.coord.dot(xd.xyz()));
    let a2 = vb.coord.dot(yd.xyz()).atan2(vb.coord.dot(xd.xyz()));
    let mut e = b.make_edge_circle(&ax2, radius, a1, a2);
    let v1 = b.make_vertex(*t1, 0.0);
    let v2 = b.make_vertex(*t2, 0.0);
    b.add(&mut e.0, &v1.0);
    b.add(&mut e.0, &v2.0);
    Ok(e)
}

/// Replace the sharp corner `p1→p2→p3` with a three-edge wire: line `p1→t1`,
/// circular arc `t1→t2`, line `t2→p3`.
pub fn fillet_corner(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt, radius: f64) -> Result<Wire, String> {
    let chain = corner_fillet(p1, p2, p3, radius)?;
    let edges: Vec<Edge> = chain.iter().map(|fe| fe.edge.clone()).collect();
    Ok(TopoBuilder::new().make_wire(&edges))
}

/// Replace the two edges meeting at `corner_index` (`edges[corner_index]` and
/// `edges[corner_index + 1]`, which share a vertex) with the filleted chain for
/// that corner. Returns the updated edge list.
pub fn fillet_edges(edges: &[Edge], corner_index: usize, radius: f64) -> Result<Vec<Edge>, String> {
    if edges.len() < 2 || corner_index + 1 >= edges.len() {
        return Err("fillet_edges: corner_index must satisfy i+1 < edges.len()".to_string());
    }
    let (p1, p2, p3) = corner_points(&edges[corner_index], &edges[corner_index + 1])?;
    let chain = corner_fillet(&p1, &p2, &p3, radius)?;
    let mut out: Vec<Edge> = edges[..corner_index].to_vec();
    out.extend(chain.iter().map(|fe| fe.edge.clone()));
    out.extend_from_slice(&edges[corner_index + 2..]);
    Ok(out)
}

/// Fillet multiple corners of a closed polygon wire. `corner_indices[i]` is the
/// corner where `edges[i]` and `edges[i + 1]` meet (indices wrap around for a
/// closed wire). Corners with an adjacent edge shorter than `radius` are
/// skipped (their edges are kept unchanged).
pub fn fillet_wire(wire: &Wire, radius: f64, corner_indices: &[usize]) -> Result<Wire, String> {
    let edges = edges_of_wire(wire);
    let n = edges.len();
    if n < 3 {
        return Ok(wire.clone());
    }
    let b = TopoBuilder::new();

    // Vertex points around the wire: edge i runs v[i] -> v[(i + 1) % n].
    let mut v: Vec<GpPnt> = Vec::with_capacity(n);
    for e in &edges {
        let (a, _) = edge_vertices(e);
        let a = a.ok_or_else(|| "fillet_wire: edge has no start vertex".to_string())?;
        v.push(vertex_position(&a));
    }

    let mut fillet = vec![false; n];
    for &ci in corner_indices {
        if ci < n {
            fillet[ci] = true;
        }
    }

    // Tangency points of each corner that can be filleted; skip (drop the flag)
    // when an adjacent edge is shorter than the radius.
    let mut ta = vec![GpPnt::zero(); n];
    let mut tb = vec![GpPnt::zero(); n];
    for i in 0..n {
        if !fillet[i] {
            continue;
        }
        let j = (i + 1) % n;
        match tangent_points(&v[i], &v[j], &v[(j + 1) % n], radius) {
            Some((a, bb)) => {
                ta[i] = a;
                tb[i] = bb;
            }
            None => fillet[i] = false,
        }
    }

    let mut out: Vec<Edge> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let prev = (i + n - 1) % n;
        // Straight middle of edge i: from the previous corner's tangency (or
        // the edge start) to this corner's tangency (or the edge end). Skipped
        // when the middle is degenerate (edge shorter than 2·radius).
        let start = if fillet[prev] { tb[prev] } else { v[i] };
        let end = if fillet[i] { ta[i] } else { v[j] };
        if start.distance(&end) > 1e-12 && start.distance(&v[i]) < end.distance(&v[i]) {
            out.push(b.make_edge_segment(&start, &end));
        }
        // Arc rounding the corner at the end of edge i.
        if fillet[i] {
            let center = arc_center(&v[i], &v[j], &v[(j + 1) % n], &ta[i], &tb[i], radius)
                .ok_or_else(|| "fillet_wire: cannot compute arc center".to_string())?;
            let cross =
                GpVec::from_pnts(&v[j], &v[i]).crossed(&GpVec::from_pnts(&v[j], &v[(j + 1) % n]));
            out.push(arc_edge(&b, &center, &cross, &ta[i], &tb[i], radius)?);
        }
    }
    Ok(b.make_wire(&out))
}

/// Fillet corners of the box's top-face wire, returning a new wire.
pub fn fillet_box(box_: &BRepPrimBox, corner_edge_indices: &[usize], radius: f64) -> Result<Wire, String> {
    let faces = faces_of(&box_.solid.0);
    let top = faces
        .iter()
        .find(|f| face_plane(f).map(|pln| pln.axis().direction().z() > 0.9).unwrap_or(false))
        .ok_or_else(|| "fillet_box: box has no top (+Z) face".to_string())?;
    let wire = wires_of_face(top)
        .into_iter()
        .next()
        .ok_or_else(|| "fillet_box: top face has no wire".to_string())?;
    fillet_wire(&wire, radius, corner_edge_indices)
}

/// Extract the ordered corner points `(p1, p2, p3)` from two adjacent edges
/// that share a vertex (`e1` runs `p1→p2`, `e2` runs `p2→p3`).
fn corner_points(e1: &Edge, e2: &Edge) -> Result<(GpPnt, GpPnt, GpPnt), String> {
    let (a1, b1) = edge_vertices(e1);
    let (a2, b2) = edge_vertices(e2);
    let (Some(a1), Some(b1), Some(a2), Some(b2)) = (a1, b1, a2, b2) else {
        return Err("corner_points: edge missing a vertex".to_string());
    };
    let (pa1, pb1, pa2, pb2) = (
        vertex_position(&a1),
        vertex_position(&b1),
        vertex_position(&a2),
        vertex_position(&b2),
    );
    if pa1.distance(&pa2) < 1e-9 {
        return Ok((pb1, pa1, pb2));
    }
    if pa1.distance(&pb2) < 1e-9 {
        return Ok((pb1, pa1, pa2));
    }
    if pb1.distance(&pa2) < 1e-9 {
        return Ok((pa1, pb1, pb2));
    }
    if pb1.distance(&pb2) < 1e-9 {
        return Ok((pa1, pb1, pa2));
    }
    Err("corner_points: edges do not share a vertex".to_string())
}

/// Length of a circular arc: radius × subtended angle.
pub fn arc_length(radius: f64, angle: f64) -> f64 {
    radius * angle.abs()
}

/// Maximum sagitta of a point set against its chord: the largest distance from
/// any point to the line through the first and last points. For points sampled
/// on a circular arc this equals the arc's sagitta.
pub fn chord_deviation(points: &[GpPnt]) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    let a = &points[0];
    let b = &points[points.len() - 1];
    let ab = GpVec::from_pnts(a, b);
    let len = ab.magnitude();
    if len < 1e-15 {
        return points.iter().map(|p| p.distance(a)).fold(0.0, f64::max);
    }
    points
        .iter()
        .map(|p| {
            let ap = GpVec::from_pnts(a, p);
            ap.crossed(&ab).magnitude() / len
        })
        .fold(0.0, f64::max)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::wire_is_closed;

    fn edge_length(e: &Edge) -> f64 {
        let curve = GeometryRegistry::global().edge_curve(&e.0).expect("edge has curve");
        let (a, b) = GeometryRegistry::global().edge_parameters(&e.0);
        let (lo, hi) = (a.min(b), a.max(b));
        if !lo.is_finite() || !hi.is_finite() {
            return 0.0;
        }
        let n = 64;
        let step = (hi - lo) / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let u0 = lo + i as f64 * step;
            let u1 = u0 + step;
            let (_, d0) = curve.d1(u0);
            let (_, d1) = curve.d1(u1);
            sum += 0.5 * (d0.magnitude() + d1.magnitude()) * step;
        }
        sum
    }

    fn unit_square_wire() -> Wire {
        let b = TopoBuilder::new();
        let p = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let e = [
            b.make_edge_segment(&p[0], &p[1]),
            b.make_edge_segment(&p[1], &p[2]),
            b.make_edge_segment(&p[2], &p[3]),
            b.make_edge_segment(&p[3], &p[0]),
        ];
        b.make_wire(&e)
    }

    #[test]
    fn tangent_points_right_angle() {
        let (p1, p2, p3) = (GpPnt::new(2., 0., 0.), GpPnt::new(0., 0., 0.), GpPnt::new(0., 2., 0.));
        let (t1, t2) = tangent_points(&p1, &p2, &p3, 0.5).expect("tangency");
        assert!(t1.is_equal(&GpPnt::new(0.5, 0., 0.)));
        assert!(t2.is_equal(&GpPnt::new(0., 0.5, 0.)));
        // Radius larger than an edge: no tangency points.
        assert!(tangent_points(&p1, &p2, &p3, 3.0).is_none());
    }

    #[test]
    fn arc_center_right_angle() {
        let (p1, p2, p3) = (GpPnt::new(2., 0., 0.), GpPnt::new(0., 0., 0.), GpPnt::new(0., 2., 0.));
        let (t1, t2) = tangent_points(&p1, &p2, &p3, 0.5).unwrap();
        let c = arc_center(&p1, &p2, &p3, &t1, &t2, 0.5).expect("center");
        assert!(c.is_equal(&GpPnt::new(0.5, 0.5, 0.)));
    }

    #[test]
    fn fillet_corner_three_edges_with_circle_arc() {
        let (p1, p2, p3) = (GpPnt::new(2., 0., 0.), GpPnt::new(0., 0., 0.), GpPnt::new(0., 2., 0.));
        let wire = fillet_corner(&p1, &p2, &p3, 0.5).expect("fillet");
        let edges = edges_of_wire(&wire);
        assert_eq!(edges.len(), 3);

        // Middle edge is a circular arc of radius 0.5 whose endpoints match the
        // tangency points.
        let arc = &edges[1];
        let curve = GeometryRegistry::global().edge_curve(&arc.0).expect("arc curve");
        assert!(curve.is_periodic(), "arc curve should be a circle");
        let (first, last) = GeometryRegistry::global().edge_parameters(&arc.0);
        let center = GpPnt::new(0.5, 0.5, 0.);
        let pa = curve.d0(first);
        let pb = curve.d0(last);
        assert!((pa.distance(&center) - 0.5).abs() < 1e-9);
        assert!((pb.distance(&center) - 0.5).abs() < 1e-9);
        assert!(pa.is_equal(&GpPnt::new(0.5, 0., 0.)));
        assert!(pb.is_equal(&GpPnt::new(0., 0.5, 0.)));
    }

    #[test]
    fn fillet_wire_unit_square_all_corners() {
        let wire = unit_square_wire();
        let f = fillet_wire(&wire, 0.2, &[0, 1, 2, 3]).expect("fillet all corners");
        let edges = edges_of_wire(&f);
        // 4 straight middles (one per side, length 1 − 2·0.2) + 4 corner arcs.
        assert_eq!(edges.len(), 8);
        assert!(wire_is_closed(&f));
        let total: f64 = edges.iter().map(edge_length).sum();
        let expected = 4.0 * (1.0 - 2.0 * 0.2) + 4.0 * (PI / 2.0 * 0.2);
        assert!((total - expected).abs() < 1e-6, "length {total} != {expected}");
    }

    #[test]
    fn fillet_wire_skips_short_edges() {
        // Rectangle 0.2 × 1.0: every corner touches the short 0.2 edge, so a
        // radius of 0.3 skips all corners (the wire is returned unchanged).
        let b = TopoBuilder::new();
        let p = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(0.2, 0., 0.),
            GpPnt::new(0.2, 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let e = [
            b.make_edge_segment(&p[0], &p[1]),
            b.make_edge_segment(&p[1], &p[2]),
            b.make_edge_segment(&p[2], &p[3]),
            b.make_edge_segment(&p[3], &p[0]),
        ];
        let wire = b.make_wire(&e);
        let f = fillet_wire(&wire, 0.3, &[0, 1, 2, 3]).expect("skipped");
        assert_eq!(edges_of_wire(&f).len(), 4);
    }

    #[test]
    fn fillet_wire_partial_corners() {
        let wire = unit_square_wire();
        let f = fillet_wire(&wire, 0.2, &[0, 2]).expect("fillet two corners");
        // 4 straight sides + 2 arcs.
        assert_eq!(edges_of_wire(&f).len(), 6);
        assert!(wire_is_closed(&f));
    }

    #[test]
    fn fillet_edges_replaces_one_corner() {
        let wire = unit_square_wire();
        let edges = edges_of_wire(&wire);
        let out = fillet_edges(&edges, 1, 0.2).expect("fillet edge 1");
        // 4 edges → replace 2 with 3 → 5.
        assert_eq!(out.len(), 5);
    }

    #[test]
    fn fillet_box_top_face_corners() {
        let bx = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let f = fillet_box(&bx, &[0, 1, 2, 3], 0.2).expect("fillet box top face");
        assert_eq!(edges_of_wire(&f).len(), 8);
        assert!(wire_is_closed(&f));
    }

    #[test]
    fn arc_length_and_chord_deviation() {
        assert!((arc_length(2.0, PI) - 2.0 * PI).abs() < 1e-12);
        // Points on a unit quarter circle: sagitta = 1 − cos(π/4).
        let pts: Vec<GpPnt> = (0..=8)
            .map(|i| {
                let a = PI / 2.0 * i as f64 / 8.0;
                GpPnt::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let dev = chord_deviation(&pts);
        assert!((dev - (1.0 - (PI / 4.0).cos())).abs() < 1e-6);
    }
}
