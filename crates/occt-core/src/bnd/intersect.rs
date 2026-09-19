//! Box intersection and distance utilities (port-internal).
//!
//! **Provenance (audit A9)**: **not** an OCCT translation. `Bnd_Tools`
//! (`Bnd_Tools.hxx`) declares only the two `Bnd2BVH` overloads converting a
//! `Bnd_Box`/`Bnd_Box2d` to `BVH_Box`; OCCT has no ray-box or box-distance
//! helper in that package (ray/box queries live in `BVH_Tree`/`BVH_Box`
//! traversal and `IntAna`).
use crate::gp::GpPnt;
use crate::bnd::BndBox;

/// Ray-box intersection test (slab method). Returns (t_near, t_far) if ray hits.
pub fn ray_box_intersect(box3d: &BndBox, origin: &GpPnt, direction: &[f64; 3]) -> Option<(f64, f64)> {
    let (xmin, xmax, ymin, ymax, zmin, zmax) = box3d.get()?;
    let inv_dir = [1.0/direction[0], 1.0/direction[1], 1.0/direction[2]];

    let t1 = (xmin - origin.x()) * inv_dir[0];
    let t2 = (xmax - origin.x()) * inv_dir[0];
    let t3 = (ymin - origin.y()) * inv_dir[1];
    let t4 = (ymax - origin.y()) * inv_dir[1];
    let t5 = (zmin - origin.z()) * inv_dir[2];
    let t6 = (zmax - origin.z()) * inv_dir[2];

    let tmin = f64::max(f64::max(t1.min(t2), t3.min(t4)), t5.min(t6));
    let tmax = f64::min(f64::min(t1.max(t2), t3.max(t4)), t5.max(t6));

    if tmax < 0.0 || tmin > tmax { return None; }
    Some((tmin, tmax))
}

/// Distance between two axis-aligned boxes.
pub fn box_box_distance(a: &BndBox, b: &BndBox) -> f64 {
    let (ax0, ax1, ay0, ay1, az0, az1) = match a.get() { Some(v) => v, None => return f64::INFINITY };
    let (bx0, bx1, by0, by1, bz0, bz1) = match b.get() { Some(v) => v, None => return f64::INFINITY };

    let dx = if ax1 < bx0 { bx0 - ax1 } else if bx1 < ax0 { ax0 - bx1 } else { 0.0 };
    let dy = if ay1 < by0 { by0 - ay1 } else if by1 < ay0 { ay0 - by1 } else { 0.0 };
    let dz = if az1 < bz0 { bz0 - az1 } else if bz1 < az0 { az0 - bz1 } else { 0.0 };
    (dx*dx + dy*dy + dz*dz).sqrt()
}

/// Merge N boxes into one.
pub fn merge_boxes(boxes: &[BndBox]) -> BndBox {
    let mut result = BndBox::new();
    for b in boxes { result.add_box(b); }
    result
}

/// Check if point is inside box (inclusive bounds).
pub fn is_inside(box3d: &BndBox, point: &GpPnt) -> bool {
    !box3d.is_out(point)
}
