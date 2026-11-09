//! Curve-curve and curve-plane intersection utilities.
//! Extends ElCLib with intersection algorithms.
use crate::gp::{GpPnt, GpVec, GpLin, GpCirc, GpPln, GpPnt2d};
use crate::elib::clib;

/// Intersection of two 3D lines. Returns closest point + distance.
/// Source: inspired by GeomAPI_IntCS.
pub fn line_line_intersection(l1: &GpLin, l2: &GpLin) -> (GpPnt, f64) {
    let p1 = l1.location(); let d1 = l1.direction();
    let p2 = l2.location(); let d2 = l2.direction();

    let r = p2.coord.subtracted(&p1.coord);
    let d1d2 = d1.xyz().dot(d2.xyz());
    let d1r = d1.xyz().dot(&r);
    let d2r = d2.xyz().dot(&r);
    let denom = 1.0 - d1d2 * d1d2;

    let (t1, t2) = if denom.abs() > 1e-30 {
        let t1 = (d1r - d1d2 * d2r) / denom;
        let t2 = (d2r - d1d2 * d1r) / denom;
        (t1, t2)
    } else {
        // Parallel lines — closest points at t1=0
        (0.0, d2r)
    };

    let q1 = clib::line_value(l1, t1);
    let q2 = clib::line_value(l2, t2);
    let dist = q1.coord.subtracted(&q2.coord).modulus();
    (GpPnt::from_xyz(&q1.coord.added(&q2.coord).divided(2.0)), dist)
}

/// Shortest distance from point to 3D line.
pub fn point_line_distance(p: &GpPnt, l: &GpLin) -> f64 {
    l.distance(p)
}

/// Intersection of line with plane. Returns (point, parameter, success).
pub fn line_plane_intersection(l: &GpLin, pl: &GpPln) -> Option<(GpPnt, f64)> {
    let o = l.location(); let d = l.direction();
    let n = *pl.pos.direction().xyz();
    let p0 = pl.location();

    let denom = n.dot(d.xyz());
    if denom.abs() < 1e-30 { return None; } // parallel
    let t = n.dot(&p0.coord.subtracted(&o.coord)) / denom;
    Some((clib::line_value(l, t), t))
}

/// Intersection of circle with plane. Returns up to 2 points.
///
/// **UNPORTED**: `ElCLib` has no intersection functions; the faithful
/// implementation is `IntAna_Quadric` / `GeomAPI_IntCS` (analytic circle–plane
/// intersection), which this port does not have. The previous body returned
/// `[center, C(π/2)]` for the coplanar case — two fabricated points that OCCT
/// never produces (the coplanar case is the whole circle, i.e. **no isolated
/// intersection points**; audit A10), so that branch now returns nothing.
pub fn circle_plane_intersection(c: &GpCirc, pl: &GpPln) -> Vec<GpPnt> {
    let center = c.location();
    let normal = *pl.pos.direction().xyz();
    let p0 = pl.location();

    let d = normal.dot(&center.coord.subtracted(&p0.coord));
    if d.abs() > c.radius + 1e-12 { return vec![]; }
    if d.abs() < 1e-12 {
        // Coplanar: the whole circle lies in the plane — no isolated points.
        return vec![];
    }
    // Circle plane intersects in 1 or 2 points
    let cos_phi = -d / c.radius;
    if cos_phi < -1.0 || cos_phi > 1.0 { return vec![]; }
    let phi = cos_phi.acos();
    vec![clib::circle_value(c, phi), clib::circle_value(c, 2.0*std::f64::consts::PI - phi)]
}

/// Distance from point to 3D circle (in its plane).
pub fn point_circle_distance(p: &GpPnt, c: &GpCirc) -> f64 {
    let center = c.location();
    let n = *c.pos.direction().xyz();
    let v = p.coord.subtracted(&center.coord);
    let proj = n.multiplied(v.dot(&n));
    let radial = v.subtracted(&proj);
    (radial.modulus() - c.radius).abs()
}

/// Point on circle at given angle (0..2π).
pub fn circle_point(c: &GpCirc, angle: f64) -> GpPnt { clib::circle_value(c, angle) }

/// Line-2D-line intersection. Returns (x, y) in 2D plane.
pub fn line2d_line2d_intersection(
    p1: &GpPnt2d, d1: &crate::gp::GpDir2d,
    p2: &GpPnt2d, d2: &crate::gp::GpDir2d) -> Option<GpPnt2d> {
    let denom = d1.x * d2.y - d1.y * d2.x;
    if denom.abs() < 1e-30 { return None; }
    let dx = p2.x() - p1.x(); let dy = p2.y() - p1.y();
    let t = (dx * d2.y - dy * d2.x) / denom;
    Some(GpPnt2d::new(p1.x() + t * d1.x, p1.y() + t * d1.y))
}

/// Project point onto plane. Returns projected point.
pub fn project_point_plane(p: &GpPnt, pl: &GpPln) -> GpPnt {
    let n = *pl.pos.direction().xyz();
    let p0 = pl.location();
    let v = p.coord.subtracted(&p0.coord);
    let dist = v.dot(&n);
    GpPnt::from_xyz(&p.coord.subtracted(&n.multiplied(dist)))
}

/// Compute plane normal at a point on a sphere (outward).
pub fn sphere_normal_at(sphere_center: &GpPnt, p: &GpPnt) -> GpVec {
    let v = p.coord.subtracted(&sphere_center.coord);
    let m = v.modulus();
    if m < 1e-30 { GpVec::zero() } else { GpVec::from_xyz(&v.divided(m)) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{GpDir, GpAx1, GpAx3};

    #[test]
    fn line_plane_basic() {
        let l = GpLin::from_pnt_dir(GpPnt::new(0.,0.,0.), GpDir::from_axis(crate::gp::dir::DirAxis::Z));
        let pl = GpPln::new(GpAx3::new(GpPnt::new(0.,0.,5.), GpDir::from_axis(crate::gp::dir::DirAxis::Z), &GpDir::from_axis(crate::gp::dir::DirAxis::X)).unwrap());
        let (pt, t) = line_plane_intersection(&l, &pl).unwrap();
        assert!((pt.z() - 5.0).abs() < 1e-14);
        assert!((t - 5.0).abs() < 1e-14);
    }

    #[test]
    fn line_line_skew() {
        let l1 = GpLin::from_pnt_dir(GpPnt::new(0.,0.,0.), GpDir::from_axis(crate::gp::dir::DirAxis::X));
        let l2 = GpLin::from_pnt_dir(GpPnt::new(0.,1.,0.), GpDir::from_axis(crate::gp::dir::DirAxis::Z));
        let (_, dist) = line_line_intersection(&l1, &l2);
        assert!((dist - 1.0).abs() < 1e-14);
    }

    #[test]
    fn project_point() {
        let pl = GpPln::new(GpAx3::new(GpPnt::new(0.,0.,0.), GpDir::from_axis(crate::gp::dir::DirAxis::Z), &GpDir::from_axis(crate::gp::dir::DirAxis::X)).unwrap());
        let p = project_point_plane(&GpPnt::new(1., 2., 3.), &pl);
        assert!((p.z()).abs() < 1e-14);
        assert!((p.x() - 1.0).abs() < 1e-14);
    }
}
