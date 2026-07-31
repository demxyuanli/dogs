//! Geometric properties of polygons, polylines and triangles: lengths, areas,
//! centroids, circumcenters, point-segment queries and ray-casting.

use crate::gp::GpPnt;

fn distance(a: &GpPnt, b: &GpPnt) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    let dz = a.z() - b.z();
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Total length of the polyline through `pts`; if `is_closed`, the closing
/// edge back to the first point is included.
pub fn polyline_length(pts: &[GpPnt], is_closed: bool) -> f64 {
    if pts.len() < 2 {
        return 0.0;
    }
    let mut total = 0.0;
    for i in 0..pts.len() - 1 {
        total += distance(&pts[i], &pts[i + 1]);
    }
    if is_closed {
        total += distance(&pts[pts.len() - 1], &pts[0]);
    }
    total
}

/// Area and centroid of a 2D polygon projected onto the XY plane using the
/// shoelace formula. The centroid is returned on the XY plane (z = 0).
pub fn polygon_centroid_area(pts: &[GpPnt]) -> (GpPnt, f64) {
    if pts.len() < 3 {
        return (GpPnt::new(0.0, 0.0, 0.0), 0.0);
    }
    let n = pts.len();
    let mut cross_sum = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let (xi, yi) = (pts[i].x(), pts[i].y());
        let (xj, yj) = (pts[j].x(), pts[j].y());
        let c = xi * yj - xj * yi;
        cross_sum += c;
        cx += (xi + xj) * c;
        cy += (yi + yj) * c;
    }
    let signed = cross_sum / 2.0;
    let area = signed.abs();
    if area < 1e-12 {
        return (GpPnt::new(0.0, 0.0, 0.0), 0.0);
    }
    (
        GpPnt::new(cx / (6.0 * signed), cy / (6.0 * signed), 0.0),
        area,
    )
}

/// Area of triangle `abc` via half the magnitude of the cross product.
pub fn triangle_area(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    let abx = b.x() - a.x();
    let aby = b.y() - a.y();
    let abz = b.z() - a.z();
    let acx = c.x() - a.x();
    let acy = c.y() - a.y();
    let acz = c.z() - a.z();
    let cx = aby * acz - abz * acy;
    let cy = abz * acx - abx * acz;
    let cz = abx * acy - aby * acx;
    (cx * cx + cy * cy + cz * cz).sqrt() / 2.0
}

/// 2D circumcenter of triangle `abc`; `None` if the triangle is collinear.
pub fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let (ax, ay) = (a.x(), a.y());
    let (bx, by) = (b.x(), b.y());
    let (cx, cy) = (c.x(), c.y());
    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if d.abs() < 1e-12 {
        return None;
    }
    let a2 = ax * ax + ay * ay;
    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;
    let ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    let uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    Some(GpPnt::new(ux, uy, 0.0))
}

/// True if `p` is within `tol` of segment `ab`.
pub fn point_on_segment(p: &GpPnt, a: &GpPnt, b: &GpPnt, tol: f64) -> bool {
    distance_point_segment(p, a, b) <= tol
}

/// Distance from point `p` to segment `ab`.
pub fn distance_point_segment(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let abx = b.x() - a.x();
    let aby = b.y() - a.y();
    let abz = b.z() - a.z();
    let len2 = abx * abx + aby * aby + abz * abz;
    if len2 < 1e-12 {
        return distance(p, a);
    }
    let t = (((p.x() - a.x()) * abx + (p.y() - a.y()) * aby + (p.z() - a.z()) * abz) / len2)
        .clamp(0.0, 1.0);
    let px = a.x() + t * abx;
    let py = a.y() + t * aby;
    let pz = a.z() + t * abz;
    let dx = p.x() - px;
    let dy = p.y() - py;
    let dz = p.z() - pz;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// 2D intersection of segments `a1-b1` and `a2-b2` on the XY plane;
/// `None` if parallel or non-intersecting.
pub fn segment_intersection_2d(a1: &GpPnt, b1: &GpPnt, a2: &GpPnt, b2: &GpPnt) -> Option<GpPnt> {
    let rx = b1.x() - a1.x();
    let ry = b1.y() - a1.y();
    let sx = b2.x() - a2.x();
    let sy = b2.y() - a2.y();
    let denom = rx * sy - ry * sx;
    if denom.abs() < 1e-12 {
        return None;
    }
    let qpx = a2.x() - a1.x();
    let qpy = a2.y() - a1.y();
    let t = (qpx * sy - qpy * sx) / denom;
    let u = (qpx * ry - qpy * rx) / denom;
    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
        Some(GpPnt::new(a1.x() + t * rx, a1.y() + t * ry, 0.0))
    } else {
        None
    }
}

/// Ray-casting point-in-polygon test on the XY plane.
pub fn polygon_contains_point_2d(poly: &[GpPnt], p: &GpPnt) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let (px, py) = (p.x(), p.y());
    let mut inside = false;
    for i in 0..n {
        let j = (i + 1) % n;
        let (xi, yi) = (poly[i].x(), poly[i].y());
        let (xj, yj) = (poly[j].x(), poly[j].y());
        if (yi > py) != (yj > py) {
            let x_int = xi + (py - yi) * (xj - xi) / (yj - yi);
            if px < x_int {
                inside = !inside;
            }
        }
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "{a} != {b}");
    }

    fn pt(x: f64, y: f64) -> GpPnt {
        GpPnt::new(x, y, 0.0)
    }

    #[test]
    fn polyline_square() {
        let sq = [pt(0.0, 0.0), pt(1.0, 0.0), pt(1.0, 1.0), pt(0.0, 1.0)];
        assert_approx(polyline_length(&sq, true), 4.0, 1e-12);
    }

    #[test]
    fn triangle_area_2d() {
        assert_approx(
            triangle_area(&pt(0.0, 0.0), &pt(1.0, 0.0), &pt(0.0, 1.0)),
            0.5,
            1e-12,
        );
    }

    #[test]
    fn distance_point_to_segment() {
        let d = distance_point_segment(&pt(0.5, 1.0), &pt(0.0, 0.0), &pt(1.0, 0.0));
        assert_approx(d, 1.0, 1e-12);
        // Beyond the segment end clamps to the endpoint.
        let d = distance_point_segment(&pt(2.0, 3.0), &pt(0.0, 0.0), &pt(1.0, 0.0));
        assert_approx(d, 10.0f64.sqrt(), 1e-12);
    }

    #[test]
    fn segment_intersection_cross() {
        let p =
            segment_intersection_2d(&pt(0.0, 0.0), &pt(1.0, 1.0), &pt(0.0, 1.0), &pt(1.0, 0.0))
                .expect("crossing segments intersect");
        assert_approx(p.x(), 0.5, 1e-12);
        assert_approx(p.y(), 0.5, 1e-12);
    }
}
