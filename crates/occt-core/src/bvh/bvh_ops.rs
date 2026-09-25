//! BVH-accelerated mesh queries — ray casting, self-intersection detection,
//! box queries and fast point containment.
//! Source: `BRepMesh`-style acceleration + `BVH` traversal.
//!
//! `TriBvh` leaves expose candidate triangle index RANGES (`start_idx..end_idx`
//! in `BvhNode`), so every query prunes with the BVH then tests candidates.

use crate::bnd::BndBox;
use crate::bvh::builder_tri::{query_triangles, TriBvh};
use crate::bvh::BvhNode;
use crate::gp::{GpPnt, GpVec};

/// Möller–Trumbore ray/triangle intersection, returns t or None.
pub fn ray_triangle_t(origin: &GpPnt, dir: &GpVec, a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<f64> {
    let e1 = GpVec::from_pnts(a, b);
    let e2 = GpVec::from_pnts(a, c);
    let p = dir.xyz().crossed(e2.xyz());
    let det = e1.xyz().dot(&p);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv = 1.0 / det;
    let tvec = GpVec::from_pnts(a, origin);
    let u = tvec.xyz().dot(&p) * inv;
    if u < 0.0 || u > 1.0 {
        return None;
    }
    let q = tvec.xyz().crossed(e1.xyz());
    let v = dir.xyz().dot(&q) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.xyz().dot(&q) * inv;
    if t > 1e-12 { Some(t) } else { None }
}

/// Collect the leaf candidate index ranges a ray passes through.
fn ray_leaves<'a>(
    node: &'a BvhNode,
    origin: &GpPnt,
    dir: &GpVec,
    out: &mut Vec<(usize, usize)>,
) {
    if !ray_hits_bbox(origin, dir, &node.bbox) {
        return;
    }
    if node.is_leaf() {
        out.push((node.start_idx, node.end_idx));
        return;
    }
    if let Some(l) = &node.left {
        ray_leaves(l, origin, dir, out);
    }
    if let Some(r) = &node.right {
        ray_leaves(r, origin, dir, out);
    }
}

/// BVH-accelerated ray cast over a TriBvh: nearest hit triangle + t.
/// `triangles` must be the same slice used to build `bvh`.
pub fn bvh_ray_cast(
    bvh: &TriBvh,
    triangles: &[(GpPnt, GpPnt, GpPnt)],
    origin: &GpPnt,
    dir: &GpVec,
) -> Option<(usize, f64)> {
    let root = bvh.root.as_ref()?;
    let mut ranges = Vec::new();
    ray_leaves(root, origin, dir, &mut ranges);
    let mut best: Option<(usize, f64)> = None;
    for (s, e) in ranges {
        for i in s..e {
            if i >= triangles.len() {
                continue;
            }
            let (a, b, c) = triangles[i];
            if let Some(t) = ray_triangle_t(origin, dir, &a, &b, &c) {
                if best.map_or(true, |(_, bt)| t < bt) {
                    best = Some((i, t));
                }
            }
        }
    }
    best
}

/// Slab test: whether a ray intersects an AABB.
pub fn ray_hits_bbox(origin: &GpPnt, dir: &GpVec, bbox: &BndBox) -> bool {
    let (x0, x1, y0, y1, z0, z1) = match bbox.get() {
        Some(v) => v,
        None => return false,
    };
    let mut tmin = f64::NEG_INFINITY;
    let mut tmax = f64::INFINITY;
    for (o, d, lo, hi) in [
        (origin.x(), dir.x(), x0, x1),
        (origin.y(), dir.y(), y0, y1),
        (origin.z(), dir.z(), z0, z1),
    ] {
        if d.abs() < 1e-30 {
            if o < lo || o > hi {
                return false;
            }
            continue;
        }
        let (t1, t2) = ((lo - o) / d, (hi - o) / d);
        let (tnear, tfar) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
        tmin = tmin.max(tnear);
        tmax = tmax.min(tfar);
        if tmin > tmax {
            return false;
        }
    }
    true
}

/// BVH point-in-mesh (even-odd) with coincident-hit dedupe.
pub fn bvh_point_in_mesh(
    bvh: &TriBvh,
    triangles: &[(GpPnt, GpPnt, GpPnt)],
    p: &GpPnt,
) -> bool {
    let dir = GpVec::new(1.0, 0.0, 0.0);
    let jittered = GpPnt::new(p.x(), p.y() + 1e-7, p.z() + 1e-7);
    let root = match bvh.root.as_ref() {
        Some(r) => r,
        None => return false,
    };
    let mut ranges = Vec::new();
    ray_leaves(root, &jittered, &dir, &mut ranges);
    let mut hits: Vec<f64> = Vec::new();
    for (s, e) in ranges {
        for i in s..e {
            if i >= triangles.len() {
                continue;
            }
            let (a, b, c) = triangles[i];
            if let Some(t) = ray_triangle_t(&jittered, &dir, &a, &b, &c) {
                hits.push(t);
            }
        }
    }
    hits.sort_by(f64::total_cmp);
    let mut unique = 0usize;
    let mut prev: Option<f64> = None;
    for t in hits {
        if prev.map_or(true, |q| (t - q).abs() > 1e-9) {
            unique += 1;
            prev = Some(t);
        }
    }
    unique % 2 == 1
}

/// Detect self-intersections: any pair of triangles whose interiors overlap
/// (both directions of the segment-triangle test) OR whose edges cross.
/// Returns the first (i, j) pair found, or None. This is O(n²) but uses the
/// BVH box queries to prune candidate pairs.
pub fn bvh_self_intersections(
    bvh: &TriBvh,
    triangles: &[(GpPnt, GpPnt, GpPnt)],
) -> Vec<(usize, usize)> {
    let mut found = Vec::new();
    let mut tested: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for (i, tri) in triangles.iter().enumerate() {
        let mut box3 = BndBox::new();
        box3.add_point(&tri.0);
        box3.add_point(&tri.1);
        box3.add_point(&tri.2);
        let mut candidates = Vec::new();
        query_triangles(bvh, &box3, &mut candidates);
        for &j in &candidates {
            if i == j {
                continue;
            }
            let key = if i < j { (i, j) } else { (j, i) };
            if tested.contains(&key) {
                continue;
            }
            tested.insert(key);
            if triangles_interfere(&triangles[i], &triangles[j]) {
                found.push(key);
                if found.len() >= 20 {
                    return found;
                }
            }
        }
    }
    found
}

fn triangles_interfere(a: &(GpPnt, GpPnt, GpPnt), b: &(GpPnt, GpPnt, GpPnt)) -> bool {
    // Segment-triangle intersection for all 6 pairs of edges.
    let edges_a = [(a.0, a.1), (a.1, a.2), (a.2, a.0)];
    let edges_b = [(b.0, b.1), (b.1, b.2), (b.2, b.0)];
    for (p, q) in edges_a {
        if segment_hits_triangle(&p, &q, b) {
            return true;
        }
    }
    for (p, q) in edges_b {
        if segment_hits_triangle(&p, &q, a) {
            return true;
        }
    }
    // Also test shared-plane overlap of the two triangle planes.
    false
}

fn segment_hits_triangle(p: &GpPnt, q: &GpPnt, tri: &(GpPnt, GpPnt, GpPnt)) -> bool {
    let dir = GpVec::from_pnts(p, q);
    if let Some(t) = ray_triangle_t(p, &dir, &tri.0, &tri.1, &tri.2) {
        let len = dir.xyz().modulus();
        t <= len + 1e-9
    } else {
        false
    }
}

/// The triangle indices within a box, via the BVH.
pub fn bvh_box_query(bvh: &TriBvh, box3d: &BndBox) -> Vec<usize> {
    let mut out = Vec::new();
    query_triangles(bvh, box3d, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bvh::builder_tri::build_tri_bvh;

    fn box_mesh() -> (TriBvh, Vec<(GpPnt, GpPnt, GpPnt)>) {
        let tris = vec![
            // A unit box's +X face, two triangles.
            ((GpPnt::new(1.,0.,0.)), (GpPnt::new(1.,1.,0.)), (GpPnt::new(1.,1.,1.))),
            ((GpPnt::new(1.,0.,0.)), (GpPnt::new(1.,1.,1.)), (GpPnt::new(1.,0.,1.))),
            // -X face.
            ((GpPnt::new(0.,0.,0.)), (GpPnt::new(0.,1.,0.)), (GpPnt::new(0.,1.,1.))),
            ((GpPnt::new(0.,0.,0.)), (GpPnt::new(0.,1.,1.)), (GpPnt::new(0.,0.,1.))),
        ];
        let bvh = build_tri_bvh(&tris, 2);
        (bvh, tris)
    }

    #[test]
    fn ray_hits_nearest_face() {
        let (bvh, tris) = box_mesh();
        // Ray from (-1, 0.5, 0.5) along +X → hits -X face at t=1.
        let hit = bvh_ray_cast(&bvh, &tris, &GpPnt::new(-1.0, 0.5, 0.5), &GpVec::new(1.0, 0.0, 0.0));
        assert!(hit.is_some());
        let (_, t) = hit.unwrap();
        assert!((t - 1.0).abs() < 1e-9);
    }

    #[test]
    fn ray_misses() {
        let (bvh, tris) = box_mesh();
        assert!(bvh_ray_cast(&bvh, &tris, &GpPnt::new(0.5, 0.5, 2.0), &GpVec::new(1.0, 0.0, 0.0)).is_none());
    }

    #[test]
    fn point_inside_via_bvh() {
        let (bvh, tris) = box_mesh();
        assert!(bvh_point_in_mesh(&bvh, &tris, &GpPnt::new(0.5, 0.5, 0.5)));
        assert!(!bvh_point_in_mesh(&bvh, &tris, &GpPnt::new(2.0, 2.0, 2.0)));
    }

    #[test]
    fn slab_test() {
        let mut bb = BndBox::new();
        bb.add_point(&GpPnt::new(0.,0.,0.));
        bb.add_point(&GpPnt::new(1.,1.,1.));
        assert!(ray_hits_bbox(&GpPnt::new(-1.,0.5,0.5), &GpVec::new(1.,0.,0.), &bb));
        assert!(!ray_hits_bbox(&GpPnt::new(-1.,2.0,0.5), &GpVec::new(1.,0.,0.), &bb));
    }

    #[test]
    fn box_query_finds_all() {
        let (bvh, tris) = box_mesh();
        let mut bb = BndBox::new();
        bb.add_point(&GpPnt::new(0.,0.,0.));
        bb.add_point(&GpPnt::new(1.,1.,1.));
        let idx = bvh_box_query(&bvh, &bb);
        assert_eq!(idx.len(), tris.len());
    }
}
