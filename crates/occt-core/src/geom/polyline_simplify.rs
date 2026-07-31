//! Polyline simplification and resampling.
//! Source: `ShapeAnalysis_FreeBoundData`, `math` decimation (RDP).

use crate::gp::{GpPnt, GpPnt2d};

/// Ramer–Douglas–Peucker polyline simplification. Returns the indices of the
/// kept points (always includes the first and last).
pub fn rdp_simplify(points: &[GpPnt], tol: f64) -> Vec<usize> {
    if points.len() <= 2 {
        return (0..points.len()).collect();
    }
    let mut keep = vec![false; points.len()];
    keep[0] = true;
    keep[points.len() - 1] = true;
    rdp_rec(points, 0, points.len() - 1, tol, &mut keep);
    (0..points.len()).filter(|&i| keep[i]).collect()
}

fn rdp_rec(points: &[GpPnt], first: usize, last: usize, tol: f64, keep: &mut Vec<bool>) {
    if last <= first + 1 {
        return;
    }
    let (a, b) = (points[first], points[last]);
    let mut max_d = 0.0f64;
    let mut max_i = first;
    for i in (first + 1)..last {
        let d = point_line_dist(&points[i], &a, &b);
        if d > max_d {
            max_d = d;
            max_i = i;
        }
    }
    if max_d > tol {
        keep[max_i] = true;
        rdp_rec(points, first, max_i, tol, keep);
        rdp_rec(points, max_i, last, tol, keep);
    }
}

fn point_line_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let ap = p.coord.subtracted(&a.coord);
    let len2 = ab.dot(&ab);
    if len2 < 1e-30 {
        return ap.modulus();
    }
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let foot = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&foot).modulus()
}

/// Simplified polyline points (keeps every `stride`-th point plus the last).
pub fn decimate(points: &[GpPnt], stride: usize) -> Vec<GpPnt> {
    if stride <= 1 || points.len() <= stride {
        return points.to_vec();
    }
    let mut out = Vec::new();
    let mut i = 0;
    while i < points.len() {
        out.push(points[i]);
        i += stride;
    }
    if out.last().map_or(true, |p| p.distance(points.last().unwrap()) > 1e-12) {
        out.push(*points.last().unwrap());
    }
    out
}

/// Resample a polyline to `n` points with equal chord spacing.
pub fn resample_uniform(points: &[GpPnt], n: usize) -> Vec<GpPnt> {
    if points.len() < 2 || n <= 1 {
        return points.to_vec();
    }
    // Cumulative chord lengths.
    let mut cum = vec![0.0f64; points.len()];
    for i in 1..points.len() {
        cum[i] = cum[i - 1] + points[i - 1].distance(&points[i]);
    }
    let total = cum[points.len() - 1];
    if total < 1e-30 {
        return points.to_vec();
    }
    let mut out = Vec::with_capacity(n);
    let mut seg = 0usize;
    for k in 0..n {
        let target = total * k as f64 / (n - 1) as f64;
        while seg + 1 < points.len() - 1 && cum[seg + 1] < target {
            seg += 1;
        }
        let seg_len = cum[seg + 1] - cum[seg];
        let frac = if seg_len > 1e-30 {
            ((target - cum[seg]) / seg_len).clamp(0.0, 1.0)
        } else {
            0.0
        };
        out.push(lerp(&points[seg], &points[seg + 1], frac));
    }
    out
}

fn lerp(a: &GpPnt, b: &GpPnt, t: f64) -> GpPnt {
    GpPnt::new(
        a.x() + t * (b.x() - a.x()),
        a.y() + t * (b.y() - a.y()),
        a.z() + t * (b.z() - a.z()),
    )
}

/// 2D version of RDP.
pub fn rdp_simplify2d(points: &[GpPnt2d], tol: f64) -> Vec<usize> {
    let p3: Vec<GpPnt> = points.iter().map(|p| GpPnt::new(p.x(), p.y(), 0.0)).collect();
    rdp_simplify(&p3, tol)
}

/// Total error (max deviation) introduced by keeping only the RDP indices.
pub fn rdp_max_deviation(points: &[GpPnt], indices: &[usize]) -> f64 {
    let mut max_d: f64 = 0.0;
    for (k, &i) in indices.iter().enumerate() {
        if k + 1 >= indices.len() {
            break;
        }
        let a = points[i];
        let b = points[indices[k + 1]];
        for j in (i + 1)..indices[k + 1] {
            let d = point_line_dist(&points[j], &a, &b);
            max_d = max_d.max(d);
        }
    }
    max_d
}

#[cfg(test)]
mod tests {
    use super::*;

    fn straight() -> Vec<GpPnt> {
        (0..=100).map(|i| GpPnt::new(i as f64, 0.0, 0.0)).collect()
    }

    #[test]
    fn straight_line_reduces_to_2() {
        let pts = straight();
        let idx = rdp_simplify(&pts, 0.5);
        assert_eq!(idx.len(), 2, "straight line → endpoints only: {idx:?}");
        assert!(rdp_max_deviation(&pts, &idx) < 0.5);
    }

    #[test]
    fn sharp_peak_kept() {
        let pts = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.),
            GpPnt::new(2.,1.,0.), GpPnt::new(3.,0.,0.), GpPnt::new(4.,0.,0.),
        ];
        let idx = rdp_simplify(&pts, 0.5);
        // The peak at index 2 is kept (deviation 1.0 > tol); its neighbors
        // deviate 0.447 < tol and are dropped.
        assert!(idx.contains(&2), "peak kept: {idx:?}");
        assert_eq!(idx.len(), 3);
    }

    #[test]
    fn decimate_and_resample() {
        let pts = straight();
        let d = decimate(&pts, 10);
        assert!(d.len() < 20);
        assert!(d.last().unwrap().distance(&pts.last().unwrap()) < 1e-12);

        let r = resample_uniform(&pts, 5);
        assert_eq!(r.len(), 5);
        // First/last match; spacing is uniform (chord).
        assert!(r[0].distance(&pts[0]) < 1e-9);
        assert!(r[4].distance(&pts[100]) < 1e-9);
        assert!((r[1].x() - 25.0).abs() < 1e-6);
    }

    #[test]
    fn rdp_2d_matches_3d() {
        let pts2 = vec![
            GpPnt2d::new(0.,0.), GpPnt2d::new(1.,0.),
            GpPnt2d::new(2.,1.), GpPnt2d::new(3.,0.), GpPnt2d::new(4.,0.),
        ];
        let idx = rdp_simplify2d(&pts2, 0.1);
        assert!(idx.contains(&2));
    }
}
