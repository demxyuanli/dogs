//! Triangle-based BVH builder. Source: `BVH_Builder` + `BRepMesh`
use crate::gp::GpPnt;
use crate::bnd::BndBox;
use crate::bvh::BvhNode;

/// BVH over triangles: stores node tree + per-triangle bboxes.
#[derive(Debug)]
pub struct TriBvh {
    pub root: Option<BvhNode>,
    pub tri_bboxes: Vec<BndBox>,
}

/// Build BVH from a triangle list.
pub fn build_tri_bvh(triangles: &[(GpPnt, GpPnt, GpPnt)], max_leaf: usize) -> TriBvh {
    let tri_bboxes: Vec<BndBox> = triangles.iter().map(|(a, b, c)| {
        let mut bb = BndBox::new();
        bb.add_point(a); bb.add_point(b); bb.add_point(c);
        bb
    }).collect();
    if tri_bboxes.is_empty() {
        return TriBvh { root: None, tri_bboxes };
    }
    let mut indices: Vec<usize> = (0..triangles.len()).collect();
    let root = build_rec(&indices, &tri_bboxes, max_leaf);
    TriBvh { root: Some(root), tri_bboxes }
}

fn bbox_of(indices: &[usize], boxes: &[BndBox]) -> BndBox {
    let mut b = BndBox::new();
    for &i in indices { b.add_box(&boxes[i]); }
    b
}

fn build_rec(indices: &[usize], boxes: &[BndBox], max_leaf: usize) -> BvhNode {
    let bbox = bbox_of(indices, boxes);
    if indices.len() <= max_leaf {
        let (s, e) = (indices[0], indices[indices.len()-1] + 1);
        return BvhNode::leaf(bbox, s, e);
    }
    // Find longest axis of combined bbox
    let (mn, mx) = match bbox.get() { Some(v) => ((v.0,v.1,v.2),(v.3,v.4,v.5)), None => return BvhNode::leaf(bbox, indices[0], indices[indices.len()-1]+1) };
    let dx = mx.0 - mn.0; let dy = mx.1 - mn.1; let dz = mx.2 - mn.2;
    let axis = if dx >= dy && dx >= dz { 0 } else if dy >= dz { 1 } else { 2 };

    // Sort indices by bbox center along axis, split at median
    let mut idx_sorted: Vec<usize> = indices.to_vec();
    idx_sorted.sort_by(|&x, &y| {
        let cx = center(&boxes[x], axis);
        let cy = center(&boxes[y], axis);
        cx.partial_cmp(&cy).unwrap()
    });
    let mid = idx_sorted.len() / 2;
    let left = build_rec(&idx_sorted[..mid], boxes, max_leaf);
    let right = build_rec(&idx_sorted[mid..], boxes, max_leaf);
    BvhNode::internal(bbox, left, right)
}

fn center(b: &BndBox, axis: usize) -> f64 {
    let (mn, mx) = match b.get() {
        Some(v) => ((v.0,v.1,v.2),(v.3,v.4,v.5)),
        None => return 0.0,
    };
    let arr = [(mn.0,mx.0),(mn.1,mx.1),(mn.2,mx.2)];
    0.5 * (arr[axis].0 + arr[axis].1)
}

/// Query all triangle indices whose bbox overlaps the query box.
pub fn query_triangles(tri: &TriBvh, box3d: &BndBox, out: &mut Vec<usize>) {
    if let Some(root) = &tri.root {
        query_rec(root, box3d, out);
    }
}

fn query_rec(node: &BvhNode, q: &BndBox, out: &mut Vec<usize>) {
    if node.bbox.is_out_box(q) { return; }
    if node.is_leaf() {
        for i in node.start_idx..node.end_idx { out.push(i); }
        return;
    }
    if let Some(l) = &node.left { query_rec(l, q, out); }
    if let Some(r) = &node.right { query_rec(r, q, out); }
}

/// Distance from point to triangle (via projection + barycentric clamping).
pub fn point_triangle_distance(p: &GpPnt, a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let ac = c.coord.subtracted(&a.coord);
    let ap = p.coord.subtracted(&a.coord);
    let n = ab.crossed(&ac);
    let n2 = n.square_modulus();
    if n2 < 1e-30 { return (p.coord.subtracted(&a.coord)).modulus(); }
    // Barycentric of projection of p onto plane
    let d = ap.dot(&n) / n2;
    let proj = p.coord.subtracted(&n.multiplied(d));
    let v0 = c.coord.subtracted(&a.coord);
    let v1 = b.coord.subtracted(&a.coord);
    let v2 = proj.subtracted(&a.coord);
    let d00 = v0.dot(&v0); let d01 = v0.dot(&v1);
    let d11 = v1.dot(&v1); let d20 = v2.dot(&v0); let d21 = v2.dot(&v1);
    let denom = d00*d11 - d01*d01;
    if denom.abs() < 1e-30 { return (proj.subtracted(&a.coord)).modulus(); }
    let mut v = (d11*d20 - d01*d21) / denom;
    let mut w = (d00*d21 - d01*d20) / denom;
    let mut u = 1.0 - v - w;
    // Clamp to triangle
    if u < 0.0 { u = 0.0; }
    if v < 0.0 { v = 0.0; }
    if w < 0.0 { w = 0.0; }
    let sum = u + v + w;
    if sum > 1e-30 { u /= sum; v /= sum; w /= sum; }
    let closest = a.coord.multiplied(u).added(&b.coord.multiplied(v)).added(&c.coord.multiplied(w));
    (p.coord.subtracted(&closest)).modulus()
}

/// Find closest triangle in BVH to a point.
pub fn closest_triangle_to_point(tri: &TriBvh, p: &GpPnt, triangles: &[(GpPnt,GpPnt,GpPnt)]) -> Option<usize> {
    let mut best = None;
    let mut best_d = f64::INFINITY;
    let mut all = Vec::new();
    let mut q = BndBox::new();
    q.add_point(p);
    query_triangles(tri, &q, &mut all);
    for &i in &all {
        if i >= triangles.len() { continue; }
        let (a, b, c) = triangles[i];
        let d = point_triangle_distance(p, &a, &b, &c);
        if d < best_d { best_d = d; best = Some(i); }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_and_query() {
        let tris = vec![
            (GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.)),
            (GpPnt::new(5.,5.,5.), GpPnt::new(6.,5.,5.), GpPnt::new(5.,6.,5.)),
        ];
        let bvh = build_tri_bvh(&tris, 2);
        let mut q = BndBox::new();
        q.add_point(&GpPnt::new(0.1, 0.1, 0.));
        let mut hits = Vec::new();
        query_triangles(&bvh, &q, &mut hits);
        assert!(hits.contains(&0), "hits: {hits:?}");
    }

    #[test]
    fn tri_distance() {
        // Triangle in z=0 plane, point above center
        let d = point_triangle_distance(
            &GpPnt::new(0.5, 0.5, 2.0),
            &GpPnt::new(0.,0.,0.), &GpPnt::new(1.,0.,0.), &GpPnt::new(0.,1.,0.),
        );
        assert!((d - 2.0).abs() < 1e-10, "d={d}");
    }
}
