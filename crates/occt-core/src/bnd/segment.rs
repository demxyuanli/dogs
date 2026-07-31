//! 2D segment utilities with bounding boxes. Source: Bnd_Box + segment helpers.
use crate::bnd::BndBox;
use crate::gp::{GpPnt, GpPnt2d, GpVec2d};

/// 2D line segment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment2d {
    pub a: GpPnt2d,
    pub b: GpPnt2d,
}

impl Segment2d {
    pub fn new(a: GpPnt2d, b: GpPnt2d) -> Self {
        Self { a, b }
    }

    pub fn length(&self) -> f64 {
        self.a.distance(&self.b)
    }

    pub fn direction(&self) -> GpVec2d {
        let v = GpVec2d::new(self.b.x() - self.a.x(), self.b.y() - self.a.y());
        v.normalized().unwrap_or(v)
    }

    /// Axis-aligned bounding box as (min, max) corners.
    pub fn bbox(&self) -> (GpPnt2d, GpPnt2d) {
        (
            GpPnt2d::new(self.a.x().min(self.b.x()), self.a.y().min(self.b.y())),
            GpPnt2d::new(self.a.x().max(self.b.x()), self.a.y().max(self.b.y())),
        )
    }
}

/// Does the segment intersect the axis-aligned 2D box (slab method, exact)?
pub fn segment_bbox_intersects(seg: &Segment2d, box_min: &GpPnt2d, box_max: &GpPnt2d) -> bool {
    let dx = seg.b.x() - seg.a.x();
    let dy = seg.b.y() - seg.a.y();
    let mut tmin = 0.0;
    let mut tmax = 1.0;
    slab(seg.a.x(), dx, box_min.x(), box_max.x(), &mut tmin, &mut tmax)
        && slab(seg.a.y(), dy, box_min.y(), box_max.y(), &mut tmin, &mut tmax)
}

/// Exact 2D segment-segment intersection via cross products.
/// None if parallel or if the intersection lies outside either segment.
pub fn segment_segment_intersect_2d(seg1: &Segment2d, seg2: &Segment2d) -> Option<GpPnt2d> {
    let r = GpVec2d::new(seg1.b.x() - seg1.a.x(), seg1.b.y() - seg1.a.y());
    let s = GpVec2d::new(seg2.b.x() - seg2.a.x(), seg2.b.y() - seg2.a.y());
    let qmp = GpVec2d::new(seg2.a.x() - seg1.a.x(), seg2.a.y() - seg1.a.y());
    let rxs = r.crossed(&s);
    if rxs.abs() < crate::precision::RESOLUTION { return None; }
    let t = qmp.crossed(&s) / rxs;
    let u = qmp.crossed(&r) / rxs;
    if t < 0.0 || t > 1.0 || u < 0.0 || u > 1.0 { return None; }
    Some(GpPnt2d::new(seg1.a.x() + t * r.x(), seg1.a.y() + t * r.y()))
}

/// Does the 3D segment [a,b] intersect the box (slab/Liang–Barsky, exact)?
pub fn segment_intersects_bbox3d(a: &GpPnt, b: &GpPnt, bb: &BndBox) -> bool {
    let (xmin, xmax, ymin, ymax, zmin, zmax) = match bb.get() {
        Some(v) => v,
        None => return false,
    };
    let mut tmin = 0.0;
    let mut tmax = 1.0;
    slab(a.x(), b.x() - a.x(), xmin, xmax, &mut tmin, &mut tmax)
        && slab(a.y(), b.y() - a.y(), ymin, ymax, &mut tmin, &mut tmax)
        && slab(a.z(), b.z() - a.z(), zmin, zmax, &mut tmin, &mut tmax)
}

/// Bounding box of a 3D segment.
pub fn segment_bbox3d(a: &GpPnt, b: &GpPnt) -> BndBox {
    let mut bb = BndBox::new();
    bb.add_point(a);
    bb.add_point(b);
    bb
}

/// Clip the parameter range [tmin, tmax] against one slab. False if no overlap.
fn slab(origin: f64, dir: f64, lo: f64, hi: f64, tmin: &mut f64, tmax: &mut f64) -> bool {
    if dir.abs() < 1e-300 {
        return origin >= lo && origin <= hi;
    }
    let t1 = (lo - origin) / dir;
    let t2 = (hi - origin) / dir;
    let (a, b) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
    *tmin = (*tmin).max(a);
    *tmax = (*tmax).min(b);
    *tmin <= *tmax
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seg_crosses_box() {
        let bb = BndBox::from_corners(&GpPnt::new(1.0, 1.0, 1.0), &GpPnt::new(2.0, 2.0, 2.0));
        let a = GpPnt::new(0.0, 1.5, 1.5);
        let b = GpPnt::new(3.0, 1.5, 1.5);
        assert!(segment_intersects_bbox3d(&a, &b, &bb));
        assert!(!bb.is_out(&GpPnt::new(1.5, 1.5, 1.5)));
    }

    #[test]
    fn seg_outside_box() {
        let bb = BndBox::from_corners(&GpPnt::new(1.0, 1.0, 1.0), &GpPnt::new(2.0, 2.0, 2.0));
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(0.5, 0.5, 0.5);
        assert!(!segment_intersects_bbox3d(&a, &b, &bb));
    }

    #[test]
    fn seg2d_cross() {
        let s1 = Segment2d::new(GpPnt2d::new(0.0, 0.0), GpPnt2d::new(2.0, 0.0));
        let s2 = Segment2d::new(GpPnt2d::new(1.0, -1.0), GpPnt2d::new(1.0, 1.0));
        let p = segment_segment_intersect_2d(&s1, &s2).unwrap();
        assert!(p.distance(&GpPnt2d::new(1.0, 0.0)) < 1e-12);
    }

    #[test]
    fn seg2d_parallel() {
        let s1 = Segment2d::new(GpPnt2d::new(0.0, 0.0), GpPnt2d::new(1.0, 0.0));
        let s2 = Segment2d::new(GpPnt2d::new(0.0, 1.0), GpPnt2d::new(1.0, 1.0));
        assert!(segment_segment_intersect_2d(&s1, &s2).is_none());
    }

    #[test]
    fn seg2d_bbox() {
        let s = Segment2d::new(GpPnt2d::new(0.0, 0.0), GpPnt2d::new(3.0, 3.0));
        assert!(segment_bbox_intersects(&s, &GpPnt2d::new(1.0, 1.0), &GpPnt2d::new(2.0, 2.0)));
        assert!(!segment_bbox_intersects(&s, &GpPnt2d::new(1.0, -2.0), &GpPnt2d::new(2.0, -1.0)));
        let (lo, hi) = s.bbox();
        assert!(lo.distance(&GpPnt2d::new(0.0, 0.0)) < 1e-12);
        assert!(hi.distance(&GpPnt2d::new(3.0, 3.0)) < 1e-12);
    }
}
