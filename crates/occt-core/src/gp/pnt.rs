//! 3D point. Source: `gp_Pnt.hxx`

use crate::gp::xyz::GpXyz;
use crate::gp::vec::GpVec;
use crate::gp::trsf::GpTrsf;
use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::precision::RESOLUTION;

/// Three-dimensional point defined by its X, Y, Z coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpPnt {
    pub coord: GpXyz,
}

impl Default for GpPnt {
    fn default() -> Self { Self::zero() }
}

impl GpPnt {
    pub const fn zero() -> Self { Self { coord: GpXyz::new(0.0, 0.0, 0.0) } }

    pub const fn new(x: f64, y: f64, z: f64) -> Self { Self { coord: GpXyz::new(x, y, z) } }

    pub fn from_xyz(c: &GpXyz) -> Self { Self { coord: *c } }

    #[inline] pub fn x(&self) -> f64 { self.coord.x }
    #[inline] pub fn y(&self) -> f64 { self.coord.y }
    #[inline] pub fn z(&self) -> f64 { self.coord.z }

    /// The underlying `gp_XYZ`.
    #[inline] pub fn xyz(&self) -> &GpXyz { &self.coord }

    #[inline] pub fn set_x(&mut self, v: f64) { self.coord.x = v; }
    #[inline] pub fn set_y(&mut self, v: f64) { self.coord.y = v; }
    #[inline] pub fn set_z(&mut self, v: f64) { self.coord.z = v; }

    pub fn distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x - other.coord.x;
        let dy = self.coord.y - other.coord.y;
        let dz = self.coord.z - other.coord.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    pub fn square_distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x - other.coord.x;
        let dy = self.coord.y - other.coord.y;
        let dz = self.coord.z - other.coord.z;
        dx * dx + dy * dy + dz * dz
    }

    pub fn is_equal(&self, other: &Self) -> bool { self.distance(other) <= RESOLUTION }

    pub fn transform(&mut self, t: &GpTrsf) { t.transforms_xyz(&mut self.coord); }
    pub fn transformed(&self, t: &GpTrsf) -> Self { let mut c = self.coord; t.transforms_xyz(&mut c); Self { coord: c } }

    /// Point mirror: reflect through `p`
    pub fn mirror_pnt(&mut self, p: &Self) {
        self.coord.x = 2.0 * p.coord.x - self.coord.x;
        self.coord.y = 2.0 * p.coord.y - self.coord.y;
        self.coord.z = 2.0 * p.coord.z - self.coord.z;
    }
    pub fn mirrored_pnt(&self, p: &Self) -> Self {
        Self::new(2.0 * p.coord.x - self.coord.x, 2.0 * p.coord.y - self.coord.y, 2.0 * p.coord.z - self.coord.z)
    }

    /// Axis mirror: reflect through axis `a`
    pub fn mirror_ax1(&mut self, a: &GpAx1) {
        let d = a.direction().xyz();
        // project onto axis, then reflect
        let vx = self.coord.x - a.location().coord.x;
        let vy = self.coord.y - a.location().coord.y;
        let vz = self.coord.z - a.location().coord.z;
        let dot = vx * d.x + vy * d.y + vz * d.z;
        self.coord.x = 2.0 * (a.location().coord.x + dot * d.x) - self.coord.x;
        self.coord.y = 2.0 * (a.location().coord.y + dot * d.y) - self.coord.y;
        self.coord.z = 2.0 * (a.location().coord.z + dot * d.z) - self.coord.z;
    }
    pub fn mirrored_ax1(&self, a: &GpAx1) -> Self {
        let mut r = *self;
        r.mirror_ax1(a);
        r
    }

    /// Coordinate system mirror: reflect through Ax2
    pub fn mirror_ax2(&mut self, a: &GpAx2) {
        let loc_pnt = a.location();
        let p = GpPnt::from_xyz(&loc_pnt.coord);
        let a_dir = a.direction();
        let d = a_dir.xyz();
        // reflect through plane whose normal is zdir, passing through location
        let vx = self.coord.x - p.coord.x;
        let vy = self.coord.y - p.coord.y;
        let vz = self.coord.z - p.coord.z;
        let dot = vx * d.x + vy * d.y + vz * d.z;
        self.coord.x -= 2.0 * dot * d.x;
        self.coord.y -= 2.0 * dot * d.y;
        self.coord.z -= 2.0 * dot * d.z;
    }
    pub fn mirrored_ax2(&self, a: &GpAx2) -> Self {
        let mut r = *self;
        r.mirror_ax2(a);
        r
    }

    /// Rotate around axis `a` by `angle` radians
    pub fn rotate(&mut self, a: &GpAx1, angle: f64) {
        let origin = GpPnt::from_xyz(&a.location().coord);
        let dir = a.direction().xyz();
        // translate to origin
        let x = self.coord.x - origin.coord.x;
        let y = self.coord.y - origin.coord.y;
        let z = self.coord.z - origin.coord.z;
        // Rodrigues rotation
        let s = angle.sin();
        let c = angle.cos();
        let t = 1.0 - c;
        let dot = x * dir.x + y * dir.y + z * dir.z;
        let rx = c * x + s * (dir.y * z - dir.z * y) + t * dot * dir.x;
        let ry = c * y + s * (dir.z * x - dir.x * z) + t * dot * dir.y;
        let rz = c * z + s * (dir.x * y - dir.y * x) + t * dot * dir.z;
        self.coord.x = origin.coord.x + rx;
        self.coord.y = origin.coord.y + ry;
        self.coord.z = origin.coord.z + rz;
    }
    pub fn rotated(&self, a: &GpAx1, angle: f64) -> Self {
        let mut r = *self;
        r.rotate(a, angle);
        r
    }

    /// Scale relative to `p`
    pub fn scale(&mut self, p: &Self, s: f64) {
        self.coord.x = p.coord.x + s * (self.coord.x - p.coord.x);
        self.coord.y = p.coord.y + s * (self.coord.y - p.coord.y);
        self.coord.z = p.coord.z + s * (self.coord.z - p.coord.z);
    }
    pub fn scaled(&self, p: &Self, s: f64) -> Self {
        let mut r = *self;
        r.scale(p, s);
        r
    }

    /// Translate by vector
    pub fn translate_vec(&mut self, v: &GpVec) {
        self.coord.x += v.coord.x;
        self.coord.y += v.coord.y;
        self.coord.z += v.coord.z;
    }
    pub fn translated_vec(&self, v: &GpVec) -> Self {
        Self::new(self.coord.x + v.coord.x, self.coord.y + v.coord.y, self.coord.z + v.coord.z)
    }

    /// Translate by vector from p1 to p2
    pub fn translate_pnts(&mut self, p1: &Self, p2: &Self) {
        self.coord.x += p2.coord.x - p1.coord.x;
        self.coord.y += p2.coord.y - p1.coord.y;
        self.coord.z += p2.coord.z - p1.coord.z;
    }
    pub fn translated_pnts(&self, p1: &Self, p2: &Self) -> Self {
        Self::new(
            self.coord.x + p2.coord.x - p1.coord.x,
            self.coord.y + p2.coord.y - p1.coord.y,
            self.coord.z + p2.coord.z - p1.coord.z,
        )
    }
}
