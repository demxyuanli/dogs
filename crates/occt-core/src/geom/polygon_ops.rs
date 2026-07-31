//! Polygon operations — point containment, signed area, ear-clipping
//! triangulation, convex hull.
//! Source: `IntPolyh`, `BRepBuilderAPI_MakeFace` (auto-triangulation),
//! `Poly` helpers.

use crate::gp::{GpPnt, GpPnt2d, GpVec, GpXyz};

/// Signed area of a 2D polygon (positive for CCW winding). Shoelace formula.
pub fn polygon_area2d(points: &[GpPnt2d]) -> f64 {
    let n = points.len();
    if n < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        s += points[i].x() * points[j].y() - points[j].x() * points[i].y();
    }
    0.5 * s
}

/// Signed area of a planar 3D polygon (magnitude of the cross-product sum / 2).
/// Returns the scalar area with the sign of the polygon's traversal.
pub fn polygon_area3d(points: &[GpPnt]) -> f64 {
    let n = points.len();
    if n < 3 {
        return 0.0;
    }
    let mut acc = GpXyz::zero();
    for i in 0..n {
        let j = (i + 1) % n;
        acc = acc.added(&points[i].coord.crossed(&points[j].coord));
    }
    0.5 * acc.modulus()
}

/// Whether a 2D polygon winds counter-clockwise (positive signed area).
pub fn is_ccw(points: &[GpPnt2d]) -> bool {
    polygon_area2d(points) > 0.0
}

/// Even-odd point-in-polygon test for a 2D polygon (any winding).
pub fn point_in_polygon2d(points: &[GpPnt2d], p: &GpPnt2d) -> bool {
    let n = points.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (points[i].x(), points[i].y());
        let (xj, yj) = (points[j].x(), points[j].y());
        let crosses = (yi > p.y()) != (yj > p.y())
            && p.x() < (xj - xi) * (p.y() - yi) / (yj - yi) + xi;
        if crosses {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Point-in-polygon for a planar 3D polygon: project onto the dominant plane
/// of the polygon normal, then run the 2D test.
pub fn point_in_polygon3d(points: &[GpPnt], p: &GpPnt) -> bool {
    let n = points.len();
    if n < 3 {
        return false;
    }
    // Dominant axis of the normal (Newell's method) for projection.
    let mut nx = 0.0f64;
    let mut ny = 0.0f64;
    let mut nz = 0.0f64;
    for i in 0..n {
        let j = (i + 1) % n;
        nx += (points[i].y() - points[j].y()) * (points[i].z() + points[j].z());
        ny += (points[i].z() - points[j].z()) * (points[i].x() + points[j].x());
        nz += (points[i].x() - points[j].x()) * (points[i].y() + points[j].y());
    }
    let ax = nx.abs();
    let ay = ny.abs();
    let az = nz.abs();
    let p2 = if ax >= ay && ax >= az {
        GpPnt2d::new(p.y(), p.z())
    } else if ay >= az {
        GpPnt2d::new(p.x(), p.z())
    } else {
        GpPnt2d::new(p.x(), p.y())
    };
    let pts2: Vec<GpPnt2d> = if ax >= ay && ax >= az {
        points.iter().map(|q| GpPnt2d::new(q.y(), q.z())).collect()
    } else if ay >= az {
        points.iter().map(|q| GpPnt2d::new(q.x(), q.z())).collect()
    } else {
        points.iter().map(|q| GpPnt2d::new(q.x(), q.y())).collect()
    };
    point_in_polygon2d(&pts2, &p2)
}

/// Ear-clipping triangulation of a simple polygon (2D). Returns triangle index
/// triples into `points`. Works for CCW or CW polygons. Fails (returns None)
/// on self-intersecting or degenerate input.
pub fn triangulate_polygon2d(points: &[GpPnt2d]) -> Option<Vec<(usize, usize, usize)>> {
    let n = points.len();
    if n < 3 {
        return None;
    }
    let mut idx: Vec<usize> = (0..n).collect();
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();
    let area = polygon_area2d(points);
    if area.abs() < 1e-15 {
        return None;
    }
    let ccw = area > 0.0;
    let mut guard = 0usize;
    while idx.len() > 3 {
        guard += 1;
        if guard > 10_000 {
            return None;
        }
        let m = idx.len();
        let mut clipped = false;
        for i in 0..m {
            let a = idx[i];
            let b = idx[(i + 1) % m];
            let c = idx[(i + 2) % m];
            if is_ear2d(points, a, b, c, &idx, ccw) {
                tris.push((a, b, c));
                idx.remove((i + 1) % m);
                clipped = true;
                break;
            }
        }
        if !clipped {
            return None;
        }
    }
    tris.push((idx[0], idx[1], idx[2]));
    Some(tris)
}

fn is_ear2d(points: &[GpPnt2d], a: usize, b: usize, c: usize, ring: &[usize], ccw: bool) -> bool {
    let (pa, pb, pc) = (points[a], points[b], points[c]);
    // Convexity: cross(b-a, c-b) must match the polygon winding.
    let cross = (pb.x() - pa.x()) * (pc.y() - pb.y()) - (pb.y() - pa.y()) * (pc.x() - pb.x());
    if (ccw && cross <= 0.0) || (!ccw && cross >= 0.0) {
        return false;
    }
    // No other vertex inside the triangle.
    for &k in ring {
        if k == a || k == b || k == c {
            continue;
        }
        if point_in_triangle2d(&pa, &pb, &pc, &points[k]) {
            return false;
        }
    }
    true
}

/// Barycentric point-in-triangle test (2D, inclusive of edges).
pub fn point_in_triangle2d(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
    let d1 = sign2d(a, b, p);
    let d2 = sign2d(b, c, p);
    let d3 = sign2d(c, a, p);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    !(has_neg && has_pos)
}

fn sign2d(a: &GpPnt2d, b: &GpPnt2d, p: &GpPnt2d) -> f64 {
    (p.x() - b.x()) * (a.y() - b.y()) - (a.x() - b.x()) * (p.y() - b.y())
}

/// Convex hull of 2D points via Andrew's monotone chain. Returns hull vertices
/// in CCW order (collinear points excluded).
pub fn convex_hull2d(points: &[GpPnt2d]) -> Vec<GpPnt2d> {
    if points.len() <= 1 {
        return points.to_vec();
    }
    let mut pts = points.to_vec();
    pts.sort_by(|a, b| a.x().total_cmp(&b.x()).then(a.y().total_cmp(&b.y())));
    fn cross(o: &GpPnt2d, a: &GpPnt2d, b: &GpPnt2d) -> f64 {
        (a.x() - o.x()) * (b.y() - o.y()) - (a.y() - o.y()) * (b.x() - o.x())
    }
    let mut lower: Vec<GpPnt2d> = Vec::new();
    for &p in &pts {
        while lower.len() >= 2 && cross(&lower[lower.len() - 2], &lower[lower.len() - 1], &p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<GpPnt2d> = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && cross(&upper[upper.len() - 2], &upper[upper.len() - 1], &p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// Triangle area (3D) via the cross-product magnitude / 2.
pub fn triangle_area3d(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    let ab = GpVec::from_pnts(a, b);
    let ac = GpVec::from_pnts(a, c);
    0.5 * ab.xyz().crossed(ac.xyz()).modulus()
}

/// Area-weighted centroid of a planar 3D polygon (assuming near-planarity).
pub fn polygon_centroid3d(points: &[GpPnt]) -> Option<GpPnt> {
    let tris = triangulate_polygon2d_projected(points)?;
    let mut acc = GpXyz::zero();
    let mut total = 0.0;
    for (a, b, c) in tris {
        let area = triangle_area3d(&points[a], &points[b], &points[c]);
        let cx = (points[a].x() + points[b].x() + points[c].x()) / 3.0;
        let cy = (points[a].y() + points[b].y() + points[c].y()) / 3.0;
        let cz = (points[a].z() + points[b].z() + points[c].z()) / 3.0;
        acc = acc.added(&GpXyz::new(cx * area, cy * area, cz * area));
        total += area;
    }
    if total > 1e-30 {
        Some(GpPnt::from_xyz(&acc.divided(total)))
    } else {
        None
    }
}

/// Project a planar 3D polygon to 2D (dominant axis), triangulate, and return
/// the triangle index triples.
pub fn triangulate_polygon2d_projected(points: &[GpPnt]) -> Option<Vec<(usize, usize, usize)>> {
    let n = points.len();
    if n < 3 {
        return None;
    }
    let mut nx = 0.0f64;
    let mut ny = 0.0f64;
    let mut nz = 0.0f64;
    for i in 0..n {
        let j = (i + 1) % n;
        nx += (points[i].y() - points[j].y()) * (points[i].z() + points[j].z());
        ny += (points[i].z() - points[j].z()) * (points[i].x() + points[j].x());
        nz += (points[i].x() - points[j].x()) * (points[i].y() + points[j].y());
    }
    let ax = nx.abs();
    let ay = ny.abs();
    let az = nz.abs();
    let pts2: Vec<GpPnt2d> = if ax >= ay && ax >= az {
        points.iter().map(|q| GpPnt2d::new(q.y(), q.z())).collect()
    } else if ay >= az {
        points.iter().map(|q| GpPnt2d::new(q.x(), q.z())).collect()
    } else {
        points.iter().map(|q| GpPnt2d::new(q.x(), q.y())).collect()
    };
    triangulate_polygon2d(&pts2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p2(x: f64, y: f64) -> GpPnt2d {
        GpPnt2d::new(x, y)
    }

    #[test]
    fn square_area_and_containment() {
        let sq = [p2(0.,0.), p2(1.,0.), p2(1.,1.), p2(0.,1.)];
        assert!((polygon_area2d(&sq) - 1.0).abs() < 1e-12);
        assert!(point_in_polygon2d(&sq, &p2(0.5, 0.5)));
        assert!(!point_in_polygon2d(&sq, &p2(1.5, 0.5)));
        assert!(is_ccw(&sq));
    }

    #[test]
    fn triangle_area_3d() {
        let t = [GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(0.,2.,0.)];
        assert!((polygon_area3d(&t) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn ear_clipping_square() {
        let sq = [p2(0.,0.), p2(1.,0.), p2(1.,1.), p2(0.,1.)];
        let tris = triangulate_polygon2d(&sq).expect("square triangulates");
        assert_eq!(tris.len(), 2);
        let mut area = 0.0;
        for (a, b, c) in &tris {
            let (pa, pb, pc) = (sq[*a], sq[*b], sq[*c]);
            area += 0.5 * ((pb.x() - pa.x()) * (pc.y() - pa.y()) - (pb.y() - pa.y()) * (pc.x() - pa.x())).abs();
        }
        assert!((area - 1.0).abs() < 1e-12);
    }

    #[test]
    fn concave_polygon_triangulates() {
        let l = [p2(0.,0.), p2(2.,0.), p2(2.,1.), p2(1.,1.), p2(1.,2.), p2(0.,2.)];
        let tris = triangulate_polygon2d(&l).expect("L-shape triangulates");
        assert_eq!(tris.len(), 4);
    }

    #[test]
    fn convex_hull() {
        let pts = [p2(0.,0.), p2(1.,0.), p2(0.,1.), p2(1.,1.), p2(0.5, 0.5)];
        let hull = convex_hull2d(&pts);
        assert_eq!(hull.len(), 4);
        assert!((polygon_area2d(&hull) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn polygon_centroid() {
        let sq = [GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(2.,2.,0.), GpPnt::new(0.,2.,0.)];
        let c = polygon_centroid3d(&sq).expect("centroid");
        assert!((c.x() - 1.0).abs() < 1e-9 && (c.y() - 1.0).abs() < 1e-9);
    }
}
