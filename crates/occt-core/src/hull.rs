//! 2D convex hull (Andrew's monotone chain). 3D hull via gift wrapping.
use crate::gp::GpPnt;

/// 2D cross product of (b-a) × (c-a), z-component.
fn cross2(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x())
}

/// Compute 2D convex hull of points (XY plane). Returns hull vertices in CCW order.
/// Points are deduplicated. Degenerate inputs return what exists.
pub fn convex_hull_2d(points: &[GpPnt]) -> Vec<GpPnt> {
    let mut pts: Vec<GpPnt> = points.to_vec();
    if pts.len() <= 1 { return pts; }
    pts.sort_by(|a, b| {
        a.x().partial_cmp(&b.x()).unwrap().then(a.y().partial_cmp(&b.y()).unwrap())
    });
    pts.dedup_by(|a, b| (a.x() - b.x()).abs() < 1e-15 && (a.y() - b.y()).abs() < 1e-15);
    if pts.len() <= 2 { return pts; }

    // Lower hull
    let mut lower: Vec<GpPnt> = Vec::new();
    for p in &pts {
        while lower.len() >= 2 && cross2(&lower[lower.len()-2], &lower[lower.len()-1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(*p);
    }
    // Upper hull
    let mut upper: Vec<GpPnt> = Vec::new();
    for p in pts.iter().rev() {
        while upper.len() >= 2 && cross2(&upper[upper.len()-2], &upper[upper.len()-1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(*p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// Compute 2D convex hull area via shoelace.
pub fn hull_area_2d(hull: &[GpPnt]) -> f64 {
    let n = hull.len();
    if n < 3 { return 0.0; }
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += hull[i].x() * hull[j].y() - hull[j].x() * hull[i].y();
    }
    area.abs() * 0.5
}

/// Point in convex polygon test (using hull edges, CCW). Returns true if inside/on boundary.
pub fn point_in_convex_hull(hull: &[GpPnt], p: &GpPnt) -> bool {
    let n = hull.len();
    if n < 3 { return false; }
    // Sign must not change
    let mut prev_sign = 0.0f64;
    for i in 0..n {
        let j = (i + 1) % n;
        let c = cross2(&hull[i], &hull[j], p);
        let sign = if c > 1e-15 { 1.0 } else if c < -1e-15 { -1.0 } else { 0.0 };
        if sign != 0.0 {
            if prev_sign != 0.0 && sign != prev_sign { return false; }
            prev_sign = sign;
        }
    }
    true
}

/// Gift-wrapping (Jarvis march) convex hull — works for any orientation.
pub fn convex_hull_jarvis(points: &[GpPnt]) -> Vec<GpPnt> {
    let n = points.len();
    if n <= 3 { return points.to_vec(); }
    // Find leftmost-lowest point
    let mut leftmost = 0usize;
    for i in 1..n {
        if points[i].x() < points[leftmost].x()
            || (points[i].x() == points[leftmost].x() && points[i].y() < points[leftmost].y()) {
            leftmost = i;
        }
    }
    let mut hull = Vec::new();
    let mut p = leftmost;
    loop {
        hull.push(points[p]);
        let mut q = (p + 1) % n;
        for i in 0..n {
            if cross2(&points[p], &points[q], &points[i]) > 0.0 { q = i; }
        }
        p = q;
        if p == leftmost { break; }
        if hull.len() > n { break; } // safety
    }
    hull
}

/// Compute convex hull of points projected on the XY plane, then lift to 3D.
pub fn convex_hull_xy(points: &[GpPnt]) -> Vec<GpPnt> { convex_hull_2d(points) }

/// Bounding-box diagonal of the hull (a rough "size" metric).
pub fn hull_extent(hull: &[GpPnt]) -> f64 {
    if hull.is_empty() { return 0.0; }
    let mut min_x = f64::MAX; let mut max_x = f64::MIN;
    let mut min_y = f64::MAX; let mut max_y = f64::MIN;
    for p in hull {
        min_x = min_x.min(p.x()); max_x = max_x.max(p.x());
        min_y = min_y.min(p.y()); max_y = max_y.max(p.y());
    }
    ((max_x - min_x).powi(2) + (max_y - min_y).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_hull() {
        let pts = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.),
            GpPnt::new(0.5,0.5,0.), // interior point
        ];
        let hull = convex_hull_2d(&pts);
        assert_eq!(hull.len(), 4);
        assert!((hull_area_2d(&hull) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn point_inside() {
        let hull = convex_hull_2d(&[
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.),
        ]);
        assert!(point_in_convex_hull(&hull, &GpPnt::new(0.5, 0.5, 0.)));
        assert!(!point_in_convex_hull(&hull, &GpPnt::new(1.5, 0.5, 0.)));
    }

    #[test]
    fn jarvis_matches() {
        let pts = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(2.,2.,0.), GpPnt::new(0.,2.,0.),
            GpPnt::new(1.,1.,0.),
        ];
        let hull = convex_hull_jarvis(&pts);
        assert_eq!(hull.len(), 4);
    }
}
