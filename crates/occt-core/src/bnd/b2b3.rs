//! Bnd_B2/B3 — 2D and 3D bounding box infrastructure.
//! Source: `Bnd_B2.hxx`, `Bnd_B3.hxx`
//! These are template-heavy in OCCT. Rust uses Box2d = BndBox2d, Box3d = BndBox directly.

pub type B2 = crate::bnd::box2d::BndBox2d;
pub type B3 = crate::bnd::box3d::BndBox;

/// Tools: distance between boxes, merge operations
pub mod tools {
    use crate::gp::GpPnt;

    /// Compute the tight bounding box around an array of points (3D).
    pub fn get_box(points: &[GpPnt]) -> super::B3 {
        let mut b = super::B3::new();
        for p in points { b.add_point(p); }
        b
    }

    /// Check if two boxes overlap within tolerance.
    pub fn are_overlapping(a: &super::B3, b: &super::B3, tol: f64) -> bool {
        let ga = a.gap(); let gb = b.gap();
        let (ax0, ax1, ay0, ay1, az0, az1) = a.get().unwrap_or((0.,0.,0.,0.,0.,0.));
        let (bx0, bx1, by0, by1, bz0, bz1) = b.get().unwrap_or((0.,0.,0.,0.,0.,0.));
        (ax1 + ga + tol >= bx0 - gb) && (bx1 + gb + tol >= ax0 - ga)
            && (ay1 + ga + tol >= by0 - gb) && (by1 + gb + tol >= ay0 - ga)
            && (az1 + ga + tol >= bz0 - gb) && (bz1 + gb + tol >= az0 - ga)
    }
}
