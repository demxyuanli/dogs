//! Lofted solids (through sections) and global B-spline interpolation.
//! Source: `BRepOffsetAPI_ThruSections.hxx`, `BSplCLib`.

use std::sync::Arc;

use occt_core::bspl::knots::{build_uniform_knots, hunt};
use occt_core::bspl::poles::greville_abscissae;
use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec, GpXyz};
use occt_geom::{GeomBSplineCurve, GeomLine, GeomPlane, Surface};

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Solid, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edge_vertices, edges_of_wire, vertex_position};

// ============================ B-spline interpolation ============================

/// B-spline basis functions evaluated at `u`, indexed by pole (Algorithm A2.2
/// from The NURBS Book). `knots` must be a clamped uniform knot vector.
fn basis_values(knots: &[f64], degree: usize, u: f64) -> Vec<f64> {
    let n_poles = knots.len() - degree - 1;
    let s = hunt(knots, u).max(degree).min(n_poles - 1);
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    let mut n = vec![0.0; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = u - knots[s + 1 - j];
        right[j] = knots[s + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > 1e-15 { n[r] / denom } else { 0.0 };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    // n[r] = N_{s - degree + r, degree}(u).
    let mut basis = vec![0.0; n_poles];
    for r in 0..=degree {
        let gi = s - degree + r;
        if gi < n_poles {
            basis[gi] = n[r];
        }
    }
    basis
}

/// Solve A·X = B for a square A with multiple right-hand sides (each column of
/// B is one right-hand side). Partial-pivoted Gaussian elimination; `None` if
/// the matrix is singular.
fn gauss_solve(a: &[Vec<f64>], b: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let nrhs = b[0].len();
    // Augmented matrix [A | B].
    let mut m: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = a[i].clone();
            row.extend(b[i].iter().copied());
            row
        })
        .collect();
    for col in 0..n {
        let mut piv = col;
        let mut best = m[col][col].abs();
        for r in (col + 1)..n {
            if m[r][col].abs() > best {
                best = m[r][col].abs();
                piv = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if piv != col {
            m.swap(piv, col);
        }
        let pv = m[col][col];
        for r in (col + 1)..n {
            let f = m[r][col] / pv;
            if f == 0.0 {
                continue;
            }
            for c in col..n + nrhs {
                m[r][c] -= f * m[col][c];
            }
        }
    }
    let mut x = vec![vec![0.0; nrhs]; n];
    for r in (0..n).rev() {
        for j in 0..nrhs {
            let mut s = m[r][n + j];
            for c in (r + 1)..n {
                s -= m[r][c] * x[c][j];
            }
            x[r][j] = s / m[r][r];
        }
    }
    Some(x)
}

/// Global interpolation: a degree-`degree` B-spline through `points`, evaluated
/// at the Greville abscissae of a clamped uniform knot vector. For degree 1 the
/// control points are the data points; for higher degrees a collocation system
/// A·P = points is solved by Gaussian elimination.
pub fn interpolate_bspline(points: &[GpPnt], degree: usize) -> Result<GeomBSplineCurve, String> {
    let n = points.len();
    if n == 0 {
        return Err("interpolate_bspline: no points".into());
    }
    if degree == 0 {
        return Err("interpolate_bspline: degree must be >= 1".into());
    }
    if n < degree + 1 {
        return Err(format!(
            "interpolate_bspline: need at least {} points for degree {degree}",
            degree + 1
        ));
    }
    let mut knots = build_uniform_knots(n, degree);
    // ponytail: build_uniform_knots leaves the trailing clamp knot (index n)
    // at 0 when interior knots exist; restore it so the vector is clamped.
    if n > degree + 1 {
        knots[n] = 1.0;
    }
    if degree == 1 {
        return GeomBSplineCurve::new(points.to_vec(), knots, 1).map_err(|e| e.to_string());
    }
    let params = greville_abscissae(&knots, degree, n);
    let a: Vec<Vec<f64>> = params.iter().map(|&u| basis_values(&knots, degree, u)).collect();
    let rhs = |f: &dyn Fn(&GpPnt) -> f64| -> Vec<Vec<f64>> {
        points.iter().map(|p| vec![f(p)]).collect()
    };
    let px = gauss_solve(&a, &rhs(&|p: &GpPnt| p.x())).ok_or("interpolate_bspline: singular collocation matrix")?;
    let py = gauss_solve(&a, &rhs(&|p: &GpPnt| p.y())).ok_or("interpolate_bspline: singular collocation matrix")?;
    let pz = gauss_solve(&a, &rhs(&|p: &GpPnt| p.z())).ok_or("interpolate_bspline: singular collocation matrix")?;
    let poles: Vec<GpPnt> = (0..n).map(|i| GpPnt::new(px[i][0], py[i][0], pz[i][0])).collect();
    GeomBSplineCurve::new(poles, knots, degree).map_err(|e| e.to_string())
}

/// Normalized cumulative chord-length parameters of `points` (first 0, last 1).
/// Degenerate point sets fall back to uniform spacing.
pub fn chord_length_params(points: &[GpPnt]) -> Vec<f64> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let mut params = vec![0.0; n];
    let mut total = 0.0;
    for i in 1..n {
        total += points[i - 1].distance(&points[i]);
        params[i] = total;
    }
    if total <= 1e-30 {
        for (i, p) in params.iter_mut().enumerate() {
            *p = if n == 1 { 0.0 } else { i as f64 / (n - 1) as f64 };
        }
    } else {
        for p in params.iter_mut() {
            *p /= total;
        }
    }
    params
}

// ============================ Loft (ThruSections) ============================

/// Result of a lofted solid: the solid plus its boundary wires and lateral faces.
#[derive(Debug, Clone)]
pub struct LoftedSolid {
    pub solid: Solid,
    pub sections: Vec<Wire>,
    pub lateral_faces: Vec<Face>,
}

/// Plane through three non-collinear points.
fn plane_through3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Result<GpPln, String> {
    let n = GpVec::from_pnts(a, b).xyz().crossed(&GpVec::from_pnts(a, c).xyz());
    let d = GpDir::from_vec(&GpVec::from_xyz(&n)).map_err(|e| e.to_string())?;
    let x_dir = perpendicular_dir(&d)?;
    Ok(GpPln::new(GpAx3::new(*a, d, &x_dir).map_err(|e| e.to_string())?))
}

/// Any unit direction perpendicular to `d`.
fn perpendicular_dir(d: &GpDir) -> Result<GpDir, String> {
    let x = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let y = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let z = GpDir::new(0.0, 0.0, 1.0).unwrap();
    let (dx, dy, dz) = (d.dot(&x).abs(), d.dot(&y).abs(), d.dot(&z).abs());
    let base = if dx <= dy && dx <= dz { x } else if dy <= dz { y } else { z };
    d.crossed(&base).map_err(|e| e.to_string())
}

/// Plane through a polygon of points, searching for three non-collinear ones.
fn plane_through_points(pts: &[GpPnt]) -> Result<GpPln, String> {
    for i in 0..pts.len() {
        for j in (i + 1)..pts.len() {
            for k in (j + 1)..pts.len() {
                let n = GpVec::from_pnts(&pts[i], &pts[j])
                    .xyz()
                    .crossed(&GpVec::from_pnts(&pts[i], &pts[k]).xyz());
                if n.modulus() > 1e-12 {
                    return plane_through3(&pts[i], &pts[j], &pts[k]);
                }
            }
        }
    }
    Err("loft: section points are collinear".into())
}

/// Closed polygon wire from `points`: segment edges with shared vertices.
pub fn section_wire(points: &[GpPnt]) -> Result<Wire, String> {
    if points.len() < 3 {
        return Err("section_wire: need >= 3 points".into());
    }
    let b = TopoBuilder::new();
    let verts: Vec<Vertex> = points.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
    let mut edges = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let dir = GpDir::from_vec(&GpVec::from_pnts(&points[i], &points[j]))
            .map_err(|e| format!("section_wire: duplicate consecutive points: {e}"))?;
        let mut e = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(points[i], dir))),
            0.0,
            points[i].distance(&points[j]),
        );
        b.add_edge_vertices(&mut e, &verts[i], &verts[j]);
        edges.push(e);
    }
    Ok(b.make_wire(&edges))
}

/// Boundary vertices of a closed polygon wire, in traversal order.
fn section_ring_vertices(wire: &Wire) -> Result<Vec<Vertex>, String> {
    let edges = edges_of_wire(wire);
    if edges.len() < 3 {
        return Err("loft: section wire needs >= 3 edges".into());
    }
    let mut ring = Vec::with_capacity(edges.len());
    for e in &edges {
        let (a, _) = edge_vertices(e);
        ring.push(a.ok_or("loft: edge without vertices")?);
    }
    let (_, last_b) = edge_vertices(edges.last().unwrap());
    let last_pt = vertex_position(&last_b.ok_or("loft: edge without vertices")?);
    if vertex_position(&ring[0]).distance(&last_pt) > 1e-9 {
        return Err("loft: section wire is not closed".into());
    }
    Ok(ring)
}

/// A line edge through two existing vertices, registering geometry and
/// attaching the shared vertices as children.
fn edge_through(b: &TopoBuilder, v1: &Vertex, v2: &Vertex) -> Edge {
    let p1 = GeometryRegistry::global().vertex_point(&v1.0);
    let p2 = GeometryRegistry::global().vertex_point(&v2.0);
    let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2))
        .unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
    let mut e = b.make_edge(Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir))), 0.0, p1.distance(&p2));
    b.add_edge_vertices(&mut e, v1, v2);
    e
}

/// End cap: planar face through the section's points, bounded by its wire.
fn make_cap(b: &TopoBuilder, wire: &Wire, ring: &[Vertex]) -> Result<Face, String> {
    let pts: Vec<GpPnt> = ring.iter().map(vertex_position).collect();
    let pln = plane_through_points(&pts)?;
    let mut face = b.make_face_plane(&pln);
    b.add_wire(&mut face, wire);
    Ok(face)
}

/// One lateral quad face between ring `a` and ring `b` on side `k`: corners
/// a[k], a[(k+1)%n], b[(k+1)%n], b[k].
pub fn ruled_face_between(pts_a: &[GpPnt], pts_b: &[GpPnt], k: usize) -> Face {
    assert_eq!(pts_a.len(), pts_b.len(), "ruled_face_between: ring size mismatch");
    let n = pts_a.len();
    assert!(n >= 3, "ruled_face_between: need >= 3 points");
    let j = (k + 1) % n;
    let quad = [pts_a[k], pts_a[j], pts_b[j], pts_b[k]];
    let b = TopoBuilder::new();
    let verts: Vec<Vertex> = quad.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
    let mut edges = Vec::with_capacity(4);
    for i in 0..4 {
        let i2 = (i + 1) % 4;
        let dir = GpDir::from_vec(&GpVec::from_pnts(&quad[i], &quad[i2]))
            .unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
        let mut e = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(quad[i], dir))),
            0.0,
            quad[i].distance(&quad[i2]),
        );
        b.add_edge_vertices(&mut e, &verts[i], &verts[i2]);
        edges.push(e);
    }
    let wire = b.make_wire(&edges);
    let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(
        plane_through3(&quad[0], &quad[1], &quad[2]).expect("ruled_face_between: degenerate quad"),
    ));
    b.make_face(surface, &[wire])
}

/// Loft a solid through closed polygon `sections` (all with the SAME vertex
/// count): planar quad lateral faces between consecutive sections plus two end
/// caps. Vertical edges are shared between adjacent quads.
pub fn loft_sections(sections: &[Wire]) -> Result<LoftedSolid, String> {
    if sections.len() < 2 {
        return Err("loft_sections: need >= 2 sections".into());
    }
    let b = TopoBuilder::new();
    let rings: Vec<Vec<Vertex>> = sections.iter().map(section_ring_vertices).collect::<Result<_, _>>()?;
    let n = rings[0].len();
    if n < 3 {
        return Err("loft_sections: sections need >= 3 vertices".into());
    }
    for (i, r) in rings.iter().enumerate() {
        if r.len() != n {
            return Err(format!("loft_sections: section {i} has {} vertices, expected {n}", r.len()));
        }
    }
    let mut lateral_faces = Vec::new();
    for i in 0..sections.len() - 1 {
        let a = &rings[i];
        let bb = &rings[i + 1];
        let a_pts: Vec<GpPnt> = a.iter().map(vertex_position).collect();
        let b_pts: Vec<GpPnt> = bb.iter().map(vertex_position).collect();
        // Vertical edges, shared by the two adjacent lateral quads.
        let vert_edges: Vec<Edge> = (0..n).map(|k| edge_through(&b, &a[k], &bb[k])).collect();
        let base_edges = edges_of_wire(&sections[i]);
        let top_edges = edges_of_wire(&sections[i + 1]);
        if base_edges.len() != n || top_edges.len() != n {
            return Err("loft_sections: section wire edge count mismatch".into());
        }
        for k in 0..n {
            let j = (k + 1) % n;
            let wire = b.make_wire(&[
                base_edges[k].clone(),
                vert_edges[j].clone(),
                top_edges[k].clone(),
                vert_edges[k].clone(),
            ]);
            let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane_through3(
                &a_pts[k], &a_pts[j], &b_pts[j],
            )?));
            lateral_faces.push(b.make_face(surface, &[wire]));
        }
    }
    let cap0 = make_cap(&b, &sections[0], &rings[0])?;
    let cap_last = make_cap(&b, &sections[sections.len() - 1], &rings[sections.len() - 1])?;
    let mut all_faces = vec![cap0, cap_last];
    all_faces.extend(lateral_faces.iter().cloned());
    let shell = b.make_shell(&all_faces);
    let solid = b.make_solid(&[shell]);
    Ok(LoftedSolid { solid, sections: sections.to_vec(), lateral_faces })
}

/// Convenience: build each polygon into a wire, then loft the sections.
pub fn loft_polygon_sections(polys: &[Vec<GpPnt>]) -> Result<LoftedSolid, String> {
    if polys.len() < 2 {
        return Err("loft_polygon_sections: need >= 2 sections".into());
    }
    let n = polys[0].len();
    if n < 3 {
        return Err("loft_polygon_sections: sections need >= 3 points".into());
    }
    for (i, p) in polys.iter().enumerate() {
        if p.len() != n {
            return Err(format!(
                "loft_polygon_sections: section {i} has {} points, expected {n}",
                p.len()
            ));
        }
    }
    let wires: Vec<Wire> = polys.iter().map(|p| section_wire(p)).collect::<Result<_, _>>()?;
    loft_sections(&wires)
}

/// Volume of a lofted solid. Uses the tessellated mass-properties volume when
/// available; otherwise falls back to a bounding-box estimate.
pub fn loft_volume(loft: &LoftedSolid) -> f64 {
    let v = crate::brep_gprop::volume(&loft.solid, 0.1);
    if v.is_finite() && v > 0.0 {
        return v;
    }
    // ponytail: bbox fallback; only hit when the tessellated volume is unusable.
    let mut min = GpXyz::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = GpXyz::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for w in &loft.sections {
        for e in edges_of_wire(w) {
            let (a, z) = edge_vertices(&e);
            for vtx in [a, z].into_iter().flatten() {
                let p = vertex_position(&vtx).coord;
                min.x = min.x.min(p.x);
                min.y = min.y.min(p.y);
                min.z = min.z.min(p.z);
                max.x = max.x.max(p.x);
                max.y = max.y.max(p.y);
                max.z = max.z.max(p.z);
            }
        }
    }
    (max.x - min.x).max(0.0) * (max.y - min.y).max(0.0) * (max.z - min.z).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topo_tools_full::{faces_of, vertices_of, wire_is_closed};
    use occt_geom::Curve;

    #[test]
    fn interpolate_bspline_cubic_through_arc() {
        // Four points on a smooth arc (unit circle), degree 3.
        let pts = [
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(-1.0, 0.0, 0.0),
            GpPnt::new(0.0, -1.0, 0.0),
        ];
        let curve = interpolate_bspline(&pts, 3).expect("interpolation");
        let knots = curve.knots.clone();
        let params = greville_abscissae(&knots, 3, 4);
        for (i, &u) in params.iter().enumerate() {
            let p = curve.d0(u);
            assert!(p.distance(&pts[i]) < 1e-6, "point {i}: got {p:?}, expected {:?}", pts[i]);
        }
    }

    #[test]
    fn interpolate_bspline_degree1() {
        let pts = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let curve = interpolate_bspline(&pts, 1).expect("interpolation");
        let knots = curve.knots.clone();
        let params = greville_abscissae(&knots, 1, 3);
        assert_eq!(params, vec![0.0, 0.5, 1.0]);
        for (i, &u) in params.iter().enumerate() {
            let p = curve.d0(u);
            assert!(p.distance(&pts[i]) < 1e-6, "point {i} at u={u}: {p:?} vs {:?}", pts[i]);
        }
    }

    #[test]
    fn loft_sections_builds_cuboid() {
        let sq0 = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let sq1 = [
            GpPnt::new(0.0, 0.0, 2.0),
            GpPnt::new(1.0, 0.0, 2.0),
            GpPnt::new(1.0, 1.0, 2.0),
            GpPnt::new(0.0, 1.0, 2.0),
        ];
        let w0 = section_wire(&sq0).expect("section 0");
        let w1 = section_wire(&sq1).expect("section 1");
        let loft = loft_sections(&[w0, w1]).expect("loft");
        assert_eq!(faces_of(&loft.solid.0).len(), 6);
        assert_eq!(loft.lateral_faces.len(), 4);
        let v = loft_volume(&loft);
        assert!((v - 2.0).abs() < 0.15, "volume {v}");
    }

    #[test]
    fn loft_polygon_sections_triangle() {
        let tri0 = vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)];
        let tri1 = vec![GpPnt::new(0.0, 0.0, 2.0), GpPnt::new(0.5, 0.0, 2.0), GpPnt::new(0.0, 0.5, 2.0)];
        let loft = loft_polygon_sections(&[tri0.clone(), tri1]).expect("loft");
        assert_eq!(faces_of(&loft.solid.0).len(), 5);
        // Mismatched point counts must error.
        let bad = loft_polygon_sections(&[
            tri0.clone(),
            vec![
                GpPnt::zero(),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(1.0, 1.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
            ],
        ]);
        assert!(bad.is_err());
    }

    #[test]
    fn section_wire_is_closed_square() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let w = section_wire(&pts).expect("section");
        assert!(wire_is_closed(&w));
        assert_eq!(edges_of_wire(&w).len(), 4);
        assert_eq!(vertices_of(&w.0).len(), 4);
    }
}
