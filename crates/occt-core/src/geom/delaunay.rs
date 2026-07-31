//! 2D Delaunay triangulation (Bowyer–Watson), Voronoi cells, point location.
//! Source: OCCT `Poly_Triangulation` and `math_Recipes`.

use std::collections::{HashMap, HashSet};

use crate::geom::polygon_ops::{
    convex_hull2d, point_in_triangle2d, polygon_area2d, triangulate_polygon2d,
};
use crate::gp::GpPnt2d;

/// Twice the signed area of (a, b, c); positive = CCW.
fn signed_area2d(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> f64 {
    (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x())
}

/// In-circle test: is `p` strictly inside the circumcircle of (a, b, c)?
///
/// Exact-ish: sign of the oriented determinant scaled by the triangle
/// orientation. Collinear (degenerate) triangles return false.
pub fn delaunay_circumcircle_contains(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
    let orient = signed_area2d(a, b, c);
    if orient == 0.0 {
        return false;
    }
    let (ax, ay) = (a.x() - p.x(), a.y() - p.y());
    let (bx, by) = (b.x() - p.x(), b.y() - p.y());
    let (cx, cy) = (c.x() - p.x(), c.y() - p.y());
    let det = (ax * ax + ay * ay) * (bx * cy - by * cx)
        - (bx * bx + by * by) * (ax * cy - ay * cx)
        + (cx * cx + cy * cy) * (ax * by - ay * bx);
    // For a CCW triangle, p is inside iff det > 0; CW flips the sign.
    det * orient > 0.0
}

/// In-circle test with a small negative tolerance so that cocircular points
/// are treated as inside. This keeps hull vertices connected on cocircular
/// input (e.g. a square), where the strict test would otherwise orphan a point.
fn in_circumcircle_eps(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
    let orient = signed_area2d(a, b, c);
    if orient == 0.0 {
        return false;
    }
    let (ax, ay) = (a.x() - p.x(), a.y() - p.y());
    let (bx, by) = (b.x() - p.x(), b.y() - p.y());
    let (cx, cy) = (c.x() - p.x(), c.y() - p.y());
    let det = (ax * ax + ay * ay) * (bx * cy - by * cx)
        - (bx * bx + by * by) * (ax * cy - ay * cx)
        + (cx * cx + cy * cy) * (ax * by - ay * bx);
    // The determinant has units L^4; scale the tolerance with the coordinates.
    let s = ax
        .abs()
        .max(ay.abs())
        .max(bx.abs())
        .max(by.abs())
        .max(cx.abs())
        .max(cy.abs());
    let eps = 1e-12 * s * s * s * s;
    det * orient > -eps
}

/// Delaunay triangulation of distinct 2D points via Bowyer–Watson.
///
/// Returns triangle index triples in CCW order. Degenerate input (fewer than
/// three points, all-collinear or coincident) falls back to an ear-clipping
/// triangulation of the convex hull and returns an empty result when none
/// exists.
pub fn delaunay_triangulate(points: &[GpPnt2d]) -> Vec<(usize, usize, usize)> {
    let n = points.len();
    if n < 3 {
        return Vec::new();
    }

    // Bounding box of the input.
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in points {
        min_x = min_x.min(p.x());
        min_y = min_y.min(p.y());
        max_x = max_x.max(p.x());
        max_y = max_y.max(p.y());
    }
    let dmax = (max_x - min_x).max(max_y - min_y);
    if dmax <= f64::EPSILON {
        return Vec::new(); // all points coincident
    }
    let xmid = (min_x + max_x) * 0.5;
    let ymid = (min_y + max_y) * 0.5;

    // Super-triangle (CCW) containing the whole bounding box.
    let s1 = GpPnt2d::new(xmid - 2.0 * dmax, ymid - dmax);
    let s2 = GpPnt2d::new(xmid + 2.0 * dmax, ymid - dmax);
    let s3 = GpPnt2d::new(xmid, ymid + 2.0 * dmax);

    // Local vertices: super-triangle first, then the input points.
    let mut verts: Vec<GpPnt2d> = vec![s1, s2, s3];
    verts.extend_from_slice(points);

    let mut tris: Vec<(usize, usize, usize)> = vec![(0, 1, 2)];

    for idx in 3..(3 + n) {
        let p = verts[idx];

        let mut bad: Vec<(usize, usize, usize)> = Vec::new();
        for &(a, b, c) in &tris {
            if in_circumcircle_eps(&verts[a], &verts[b], &verts[c], &p) {
                bad.push((a, b, c));
            }
        }
        if bad.is_empty() {
            continue;
        }
        let bad_set: HashSet<(usize, usize, usize)> = bad.iter().copied().collect();

        // Boundary edges: undirected edges appearing in exactly one bad
        // triangle, kept with the winding direction of that triangle.
        let mut edge_counts: HashMap<(usize, usize), i32> = HashMap::new();
        let mut edge_dir: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
        for &(a, b, c) in &bad {
            for (u, v) in [(a, b), (b, c), (c, a)] {
                let key = if u < v { (u, v) } else { (v, u) };
                *edge_counts.entry(key).or_insert(0) += 1;
                edge_dir.entry(key).or_insert((u, v));
            }
        }
        let boundary: Vec<(usize, usize)> = edge_counts
            .iter()
            .filter(|&(_, &cnt)| cnt == 1)
            .map(|(&key, _)| edge_dir[&key])
            .collect();

        tris.retain(|t| !bad_set.contains(t));
        for (u, v) in boundary {
            // Keep every new triangle CCW.
            if signed_area2d(&verts[u], &verts[v], &verts[idx]) >= 0.0 {
                tris.push((u, v, idx));
            } else {
                tris.push((v, u, idx));
            }
        }
    }

    // Drop triangles touching the super-triangle, map to input indices.
    let mut out: Vec<(usize, usize, usize)> = tris
        .into_iter()
        .filter(|&(a, b, c)| a >= 3 && b >= 3 && c >= 3)
        .map(|(a, b, c)| (a - 3, b - 3, c - 3))
        .collect();

    // Drop degenerate (zero-area) triangles and deduplicate.
    out.retain(|&(a, b, c)| signed_area2d(&points[a], &points[b], &points[c]) != 0.0);
    let mut seen = HashSet::new();
    out.retain(|t| seen.insert(*t));

    if is_valid_triangulation(points, &out) {
        make_ccw(points, &mut out);
        out
    } else {
        let mut hull_tris = hull_triangulation(points);
        make_ccw(points, &mut hull_tris);
        hull_tris
    }
}

fn make_ccw(points: &[GpPnt2d], tris: &mut [(usize, usize, usize)]) {
    for t in tris.iter_mut() {
        if signed_area2d(&points[t.0], &points[t.1], &points[t.2]) < 0.0 {
            *t = (t.0, t.2, t.1);
        }
    }
}

/// A valid triangulation covers the convex hull exactly: every input point is
/// used and the summed triangle area matches the hull area (no overlaps).
fn is_valid_triangulation(points: &[GpPnt2d], tris: &[(usize, usize, usize)]) -> bool {
    if tris.is_empty() {
        return false;
    }
    let mut seen = vec![false; points.len()];
    let mut area = 0.0;
    for &(a, b, c) in tris {
        if a == b
            || b == c
            || c == a
            || a >= points.len()
            || b >= points.len()
            || c >= points.len()
        {
            return false;
        }
        seen[a] = true;
        seen[b] = true;
        seen[c] = true;
        area += 0.5 * signed_area2d(&points[a], &points[b], &points[c]).abs();
    }
    if !seen.iter().all(|&s| s) {
        return false;
    }
    let hull = convex_hull2d(points);
    let hull_area = polygon_area2d(&hull).abs();
    (area - hull_area).abs() <= 1e-9 * hull_area.max(1.0)
}

/// Fallback: triangulate the convex hull via ear clipping. Interior points are
/// dropped (this path only runs for degenerate/cocircular inputs).
fn hull_triangulation(points: &[GpPnt2d]) -> Vec<(usize, usize, usize)> {
    let hull = convex_hull2d(points);
    if hull.len() < 3 {
        return Vec::new();
    }
    let mut idx = Vec::with_capacity(hull.len());
    for h in &hull {
        match points.iter().position(|p| p.x() == h.x() && p.y() == h.y()) {
            Some(i) => idx.push(i),
            None => return Vec::new(),
        }
    }
    match triangulate_polygon2d(&hull) {
        Some(tris) => tris.into_iter().map(|(a, b, c)| (idx[a], idx[b], idx[c])).collect(),
        None => Vec::new(),
    }
}

/// Triangle adjacency by shared edge: `neighbors[t]` lists the triangles that
/// share an edge with triangle `t`.
pub fn triangle_neighbors(tris: &[(usize, usize, usize)]) -> Vec<Vec<usize>> {
    let m = tris.len();
    let mut edge_map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (i, &(a, b, c)) in tris.iter().enumerate() {
        for (u, v) in [(a, b), (b, c), (c, a)] {
            let key = if u < v { (u, v) } else { (v, u) };
            edge_map.entry(key).or_default().push(i);
        }
    }
    let mut neighbors = vec![Vec::new(); m];
    for ts in edge_map.values() {
        if ts.len() == 2 {
            neighbors[ts[0]].push(ts[1]);
            neighbors[ts[1]].push(ts[0]);
        }
    }
    for nb in neighbors.iter_mut() {
        nb.sort_unstable();
    }
    neighbors
}

/// Circumcenter of a triangle, or None if the triangle is degenerate.
fn circumcenter2d(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> Option<GpPnt2d> {
    let (ax, ay) = (a.x(), a.y());
    let (bx, by) = (b.x(), b.y());
    let (cx, cy) = (c.x(), c.y());
    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if d.abs() < 1e-300 {
        return None;
    }
    let (a2, b2, c2) = (ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy);
    let ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    let uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    Some(GpPnt2d::new(ux, uy))
}

/// Voronoi cells of a Delaunay triangulation.
///
/// Each cell is the polygon of circumcenters of the triangles incident to the
/// point, ordered around the point. Boundary cells are unbounded; their outer
/// rays are extended to a slightly enlarged bounding box (approximation).
pub fn voronoi_cells(points: &[GpPnt2d], tris: &[(usize, usize, usize)]) -> Vec<Vec<GpPnt2d>> {
    let n = points.len();
    let mut cells: Vec<Vec<GpPnt2d>> = vec![Vec::new(); n];
    if n == 0 {
        return cells;
    }

    // Triangles incident to each point.
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (ti, &(a, b, c)) in tris.iter().enumerate() {
        if a < n {
            incident[a].push(ti);
        }
        if b < n {
            incident[b].push(ti);
        }
        if c < n {
            incident[c].push(ti);
        }
    }

    // Hull membership, for boundary-cell clipping.
    let hull = convex_hull2d(points);
    let hull_set: HashSet<usize> = hull
        .iter()
        .filter_map(|h| points.iter().position(|p| p.x() == h.x() && p.y() == h.y()))
        .collect();

    // Slightly enlarged bounding box for clipping.
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in points {
        min_x = min_x.min(p.x());
        min_y = min_y.min(p.y());
        max_x = max_x.max(p.x());
        max_y = max_y.max(p.y());
    }
    let margin = 0.1 * (max_x - min_x).max(max_y - min_y).max(1e-9);
    let bbox = (min_x - margin, min_y - margin, max_x + margin, max_y + margin);

    for i in 0..n {
        let mut cell: Vec<GpPnt2d> = Vec::new();
        for &ti in &incident[i] {
            let (a, b, c) = tris[ti];
            if let Some(cc) = circumcenter2d(&points[a], &points[b], &points[c]) {
                cell.push(cc);
            }
        }
        if cell.len() < 2 {
            cells[i] = cell;
            continue;
        }
        // Order circumcenters around the point (CCW).
        cell.sort_by(|p, q| {
            let ap = (p.y() - points[i].y()).atan2(p.x() - points[i].x());
            let aq = (q.y() - points[i].y()).atan2(q.x() - points[i].x());
            ap.total_cmp(&aq)
        });
        if hull_set.contains(&i) {
            cell = clip_boundary_cell(&cell, &points[i], &bbox);
        }
        cells[i] = cell;
    }
    cells
}

/// Extend an unbounded boundary cell with two rays clipped at the bounding box.
fn clip_boundary_cell(cell: &[GpPnt2d], p: &GpPnt2d, bbox: &(f64, f64, f64, f64)) -> Vec<GpPnt2d> {
    let first = cell[0];
    let last = *cell.last().unwrap();
    let d1 = dir_from(last, p);
    let d2 = dir_from(first, p);
    let hit1 = ray_box_hit(&last, &d1, bbox);
    let hit2 = ray_box_hit(&first, &d2, bbox);
    let mut out = cell.to_vec();
    out.push(hit1);
    out.push(hit2);
    out
}

/// Unit direction from `p` toward `q`.
fn dir_from(q: GpPnt2d, p: &GpPnt2d) -> (f64, f64) {
    let dx = q.x() - p.x();
    let dy = q.y() - p.y();
    let m = (dx * dx + dy * dy).sqrt();
    if m > 1e-300 {
        (dx / m, dy / m)
    } else {
        (1.0, 0.0)
    }
}

/// Exit point of the ray `origin + t*dir` (t >= 0) from the bounding box.
fn ray_box_hit(origin: &GpPnt2d, dir: &(f64, f64), bbox: &(f64, f64, f64, f64)) -> GpPnt2d {
    let (ox, oy) = (origin.x(), origin.y());
    let (dx, dy) = *dir;
    let (min_x, min_y, max_x, max_y) = *bbox;
    let tx1 = (min_x - ox) / dx;
    let tx2 = (max_x - ox) / dx;
    let tx_min = if tx1 < tx2 { tx1 } else { tx2 };
    let ty1 = (min_y - oy) / dy;
    let ty2 = (max_y - oy) / dy;
    let ty_min = if ty1 < ty2 { ty1 } else { ty2 };
    let t = tx_min.max(ty_min).max(0.0);
    GpPnt2d::new(ox + t * dx, oy + t * dy)
}

/// Is `p` inside the triangulated region (located in any triangle)?
pub fn delaunay_contains(
    tris: &[(usize, usize, usize)],
    points: &[GpPnt2d],
    p: &GpPnt2d,
) -> bool {
    for &(a, b, c) in tris {
        if point_in_triangle2d(&points[a], &points[b], &points[c], p) {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::polygon_ops::point_in_polygon2d;

    fn p2(x: f64, y: f64) -> GpPnt2d {
        GpPnt2d::new(x, y)
    }

    fn tri_area(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> f64 {
        signed_area2d(a, b, c).abs() * 0.5
    }

    #[test]
    fn circumcircle_contains() {
        // Triangle (0,0),(2,0),(1,1) has circumcenter (1,0), radius 1.
        let a = p2(0.0, 0.0);
        let b = p2(2.0, 0.0);
        let c = p2(1.0, 1.0);
        assert!(delaunay_circumcircle_contains(&a, &b, &c, &p2(1.0, 0.5)));
        assert!(!delaunay_circumcircle_contains(&a, &b, &c, &p2(1.0, 2.0)));
    }

    #[test]
    fn square_two_triangles() {
        let pts = [p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)];
        let tris = delaunay_triangulate(&pts);
        assert_eq!(tris.len(), 2);
        let area: f64 = tris.iter().map(|&(a, b, c)| tri_area(&pts[a], &pts[b], &pts[c])).sum();
        assert!((area - 1.0).abs() < 1e-12, "area = {area}");
    }

    #[test]
    fn square_plus_center() {
        let pts = [
            p2(0.0, 0.0),
            p2(1.0, 0.0),
            p2(1.0, 1.0),
            p2(0.0, 1.0),
            p2(0.5, 0.5),
        ];
        let tris = delaunay_triangulate(&pts);
        assert_eq!(tris.len(), 4, "expected 4 triangles, got {tris:?}");
        // Center (index 4) connected to all four corners.
        let mut connected = HashSet::new();
        for &(a, b, c) in &tris {
            if a == 4 || b == 4 || c == 4 {
                connected.insert(a);
                connected.insert(b);
                connected.insert(c);
            }
        }
        for corner in 0..4 {
            assert!(connected.contains(&corner), "corner {corner} not connected");
        }
    }

    #[test]
    fn all_triangles_ccw() {
        let pts = [
            p2(0.0, 0.0),
            p2(1.0, 0.0),
            p2(1.0, 1.0),
            p2(0.0, 1.0),
            p2(0.5, 0.5),
            p2(0.25, 0.75),
            p2(0.6, 0.3),
        ];
        let tris = delaunay_triangulate(&pts);
        assert!(!tris.is_empty());
        for &(a, b, c) in &tris {
            assert!(
                signed_area2d(&pts[a], &pts[b], &pts[c]) > 0.0,
                "triangle ({a},{b},{c}) not CCW"
            );
        }
    }

    #[test]
    fn neighbors_of_square() {
        let pts = [p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)];
        let tris = delaunay_triangulate(&pts);
        let nb = triangle_neighbors(&tris);
        assert_eq!(nb.len(), 2);
        assert_eq!(nb[0], vec![1]);
        assert_eq!(nb[1], vec![0]);
    }

    #[test]
    fn three_points_one_triangle() {
        let pts = [p2(0.0, 0.0), p2(2.0, 0.0), p2(0.0, 2.0)];
        let tris = delaunay_triangulate(&pts);
        assert_eq!(tris.len(), 1);
        let area: f64 = tris.iter().map(|&(a, b, c)| tri_area(&pts[a], &pts[b], &pts[c])).sum();
        assert!((area - 2.0).abs() < 1e-12, "area = {area}");
    }

    #[test]
    fn collinear_returns_empty() {
        let pts = [p2(0.0, 0.0), p2(1.0, 0.0), p2(2.0, 0.0), p2(3.0, 0.0)];
        let tris = delaunay_triangulate(&pts);
        assert!(tris.is_empty());
        // Fewer than three points.
        assert!(delaunay_triangulate(&[p2(0.0, 0.0), p2(1.0, 0.0)]).is_empty());
    }

    #[test]
    fn no_overlapping_triangles() {
        // Total area of the triangulation must equal the convex-hull area.
        let pts = [
            p2(0.0, 0.0),
            p2(1.5, 0.1),
            p2(1.0, 1.2),
            p2(0.2, 1.0),
            p2(0.6, 0.5),
            p2(0.8, 0.8),
            p2(0.3, 0.3),
        ];
        let tris = delaunay_triangulate(&pts);
        assert!(!tris.is_empty());
        let hull_area = crate::geom::polygon_ops::polygon_area2d(&convex_hull2d(&pts)).abs();
        let total: f64 = tris.iter().map(|&(a, b, c)| tri_area(&pts[a], &pts[b], &pts[c])).sum();
        assert!(
            (total - hull_area).abs() <= 1e-9 * hull_area,
            "total {total} vs hull {hull_area}"
        );
    }

    #[test]
    fn contains_interior_point() {
        let pts = [p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)];
        let tris = delaunay_triangulate(&pts);
        assert!(delaunay_contains(&tris, &pts, &p2(0.5, 0.5)));
        assert!(!delaunay_contains(&tris, &pts, &p2(1.5, 0.5)));
    }

    #[test]
    fn voronoi_central_cell() {
        // Square + center: the center's Voronoi cell is the diamond of the four
        // circumcenters around the center.
        let pts = [
            p2(0.0, 0.0),
            p2(1.0, 0.0),
            p2(1.0, 1.0),
            p2(0.0, 1.0),
            p2(0.5, 0.5),
        ];
        let tris = delaunay_triangulate(&pts);
        let cells = voronoi_cells(&pts, &tris);
        let center_cell = &cells[4];
        assert!(center_cell.len() >= 3, "center cell too small: {center_cell:?}");
        // All vertices within the square.
        for v in center_cell {
            assert!(v.x() >= -0.1 && v.x() <= 1.1 && v.y() >= -0.1 && v.y() <= 1.1);
        }
        assert!(point_in_polygon2d(center_cell, &p2(0.5, 0.5)));
    }
}
