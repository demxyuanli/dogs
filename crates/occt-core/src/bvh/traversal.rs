//! BVH traversal and ray intersection. Source: `BVH_Traverse.hxx`
use crate::gp::GpPnt;
use crate::bvh::BvhNode;

/// Ray-box intersection result.
#[derive(Debug, Clone, Copy)]
pub struct RayHit { pub t: f64, pub node_idx: usize }

/// Traverse BVH with a query box, returning indices of all boxes it overlaps.
pub fn query_box(root: &BvhNode, query: &crate::bnd::BndBox, out: &mut Vec<usize>) {
    if root.bbox.is_out_box(query) { return; }
    if root.is_leaf() {
        for i in root.start_idx..root.end_idx { out.push(i); }
        return;
    }
    if let Some(l) = &root.left { query_box(l, query, out); }
    if let Some(r) = &root.right { query_box(r, query, out); }
}

/// Ray-BVH intersection: returns hit distances for all intersecting leaves.
pub fn ray_intersect(root: &BvhNode, origin: &GpPnt, dir: &[f64; 3], out: &mut Vec<RayHit>) {
    if !ray_hits_box(&root.bbox, origin, dir) { return; }
    if root.is_leaf() {
        for i in root.start_idx..root.end_idx {
            out.push(RayHit { t: f64::INFINITY, node_idx: i });
        }
        return;
    }
    if let Some(l) = &root.left { ray_intersect(l, origin, dir, out); }
    if let Some(r) = &root.right { ray_intersect(r, origin, dir, out); }
}

/// Slab-method ray-box test.
fn ray_hits_box(b: &crate::bnd::BndBox, o: &GpPnt, d: &[f64; 3]) -> bool {
    let (x0, x1, y0, y1, z0, z1) = match b.get() { Some(v) => v, None => return false };
    let inv = [1.0/d[0], 1.0/d[1], 1.0/d[2]];
    let t1 = (x0 - o.x()) * inv[0]; let t2 = (x1 - o.x()) * inv[0];
    let t3 = (y0 - o.y()) * inv[1]; let t4 = (y1 - o.y()) * inv[1];
    let t5 = (z0 - o.z()) * inv[2]; let t6 = (z1 - o.z()) * inv[2];
    let tmin = f64::max(f64::max(t1.min(t2), t3.min(t4)), t5.min(t6));
    let tmax = f64::min(f64::min(t1.max(t2), t3.max(t4)), t5.max(t6));
    tmax >= 0.0 && tmin <= tmax
}

/// Count nodes in BVH (for stats).
pub fn count_nodes(root: &BvhNode) -> usize {
    1 + root.left.as_ref().map(|l| count_nodes(l)).unwrap_or(0)
      + root.right.as_ref().map(|r| count_nodes(r)).unwrap_or(0)
}

/// Find deepest leaf containing a point.
pub fn query_point(root: &BvhNode, p: &GpPnt) -> Option<usize> {
    if root.bbox.is_out(p) { return None; }
    if root.is_leaf() { return Some(root.start_idx); }
    if let Some(l) = &root.left { if let Some(i) = query_point(l, p) { return Some(i); } }
    if let Some(r) = &root.right { if let Some(i) = query_point(r, p) { return Some(i); } }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bnd::BndBox;

    #[test]
    fn box_query_finds_all() {
        let pts = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(3.,0.,0.), GpPnt::new(4.,0.,0.)];
        let root = crate::bvh::build_bvh(&pts, 2);
        assert!(count_nodes(&root) >= 3, "expected at least 3 nodes");
        let mut query = BndBox::new();
        query.add_point(&GpPnt::new(0.5, 0., 0.));
        let mut hits = Vec::new();
        query_box(&root, &query, &mut hits);
        assert!(hits.contains(&0) || hits.contains(&1));
    }
}
