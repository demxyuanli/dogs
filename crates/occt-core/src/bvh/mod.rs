//! Bounding Volume Hierarchy. Source: `BVH/`
pub mod traversal;
use crate::gp::GpPnt;
use crate::bnd::BndBox;

/// BVH node — either leaf with triangle indices or internal with children.
#[derive(Debug, Clone)]
pub struct BvhNode {
    pub bbox: BndBox,
    pub left: Option<Box<BvhNode>>,
    pub right: Option<Box<BvhNode>>,
    pub start_idx: usize,
    pub end_idx: usize,
}

impl BvhNode {
    pub fn leaf(bbox: BndBox, start: usize, end: usize) -> Self {
        Self { bbox, left: None, right: None, start_idx: start, end_idx: end }
    }
    pub fn internal(bbox: BndBox, left: BvhNode, right: BvhNode) -> Self {
        Self { bbox, left: Some(Box::new(left)), right: Some(Box::new(right)), start_idx: 0, end_idx: 0 }
    }
    pub fn is_leaf(&self) -> bool { self.left.is_none() }
}

/// Build BVH from triangle centers (GpPnt array).
/// Splits along longest axis at median. max_leaf_size controls recursion depth.
pub fn build_bvh(points: &[GpPnt], max_leaf_size: usize) -> BvhNode {
    let n = points.len();
    let bbox = compute_bbox(points);
    if n <= max_leaf_size {
        return BvhNode::leaf(bbox, 0, n);
    }
    // Find longest axis
    let (mn, mx) = match bbox.get() { Some(v) => ((v.0,v.1,v.2),(v.3,v.4,v.5)), None => return BvhNode::leaf(bbox,0,n) };
    let dx = mx.0 - mn.0; let dy = mx.1 - mn.1; let dz = mx.2 - mn.2;
    let axis = if dx >= dy && dx >= dz { 0 } else if dy >= dz { 1 } else { 2 };

    // Sort by axis and split at median
    let mut indexed: Vec<(usize, f64)> = points.iter().enumerate().map(|(i, p)| (i, match axis { 0=>p.x(),1=>p.y(),_=>p.z() })).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let mid = n / 2;
    let left_pts: Vec<GpPnt> = indexed[..mid].iter().map(|&(i, _)| points[i]).collect();
    let right_pts: Vec<GpPnt> = indexed[mid..].iter().map(|&(i, _)| points[i]).collect();

    let left = build_bvh(&left_pts, max_leaf_size);
    let right = build_bvh(&right_pts, max_leaf_size);
    BvhNode::internal(bbox, left, right)
}

fn compute_bbox(points: &[GpPnt]) -> BndBox {
    let mut b = BndBox::new();
    for p in points { b.add_point(p); }
    b
}
