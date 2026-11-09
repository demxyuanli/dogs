//! Bounding sphere. Source: `Bnd_Sphere.hxx`
use crate::gp::GpPnt;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BndSphere { pub center: GpPnt, pub radius: f64, is_void: bool }

impl BndSphere {
    pub fn new() -> Self { Self { center: GpPnt::zero(), radius: 0.0, is_void: true } }
    pub fn from_center_radius(center: GpPnt, radius: f64) -> Self { Self { center, radius, is_void: false } }
    pub fn is_void(&self) -> bool { self.is_void }
    pub fn set_void(&mut self) { self.is_void = true; }

    pub fn add_point(&mut self, p: &GpPnt) {
        if self.is_void { self.center = *p; self.radius = 0.0; self.is_void = false; return; }
        let d = p.coord.subtracted(&self.center.coord).modulus();
        if d > self.radius { self.radius = d; }
    }

    /// `Bnd_Sphere::Add` (`Bnd_Sphere.cxx:73-101`): an uninitialised sphere
    /// adopts the other; a sphere that **encloses** this one replaces it
    /// entirely (center included — the previous body only grew the radius and
    /// left the center stale, audit A10); a sphere enclosed by this one is
    /// ignored; otherwise the pair expands along the center line.
    pub fn add_sphere(&mut self, other: &Self) {
        if self.is_void { *self = *other; return; }
        if other.is_void { return; }
        let d = other.center.coord.subtracted(&self.center.coord).modulus();
        if self.radius + d <= other.radius {
            // the other sphere is larger and encloses this
            *self = *other;
            return;
        }
        if other.radius + d <= self.radius {
            return; // this sphere encloses the other
        }
        let new_r = (d + self.radius + other.radius) * 0.5;
        self.radius = new_r;
        let shift = (other.radius - self.radius + d) * 0.5 / d;
        let dir = other.center.coord.subtracted(&self.center.coord).multiplied(shift);
        self.center = GpPnt::from_xyz(&self.center.coord.added(&dir));
    }

    pub fn is_out(&self, p: &GpPnt) -> bool {
        if self.is_void { return true; }
        p.coord.subtracted(&self.center.coord).square_modulus() > self.radius * self.radius
    }

    /// `Bnd_Sphere::Distance(gp_XYZ)` (`Bnd_Sphere.cxx:63-66`) — the distance to
    /// the **center** (`theNode - myCenter).Modulus()`), not to the surface.
    /// The surface distance is [`BndSphere::distances`] (`Bnd_Sphere.cxx:45-50`).
    pub fn distance(&self, p: &GpPnt) -> f64 {
        if self.is_void { return f64::INFINITY; }
        p.coord.subtracted(&self.center.coord).modulus()
    }

    /// `Bnd_Sphere::SquareDistance(gp_XYZ)` (`Bnd_Sphere.cxx:68-71`).
    pub fn square_distance(&self, p: &GpPnt) -> f64 {
        if self.is_void { return f64::INFINITY; }
        p.coord.subtracted(&self.center.coord).square_modulus()
    }

    /// `Bnd_Sphere::Distances(gp_XYZ, theMin, theMax)` (`Bnd_Sphere.cxx:45-50`):
    /// `(min, max)` distance from the point to the sphere's **surface**.
    pub fn distances(&self, p: &GpPnt) -> (f64, f64) {
        if self.is_void { return (f64::INFINITY, f64::INFINITY); }
        let d = p.coord.subtracted(&self.center.coord).modulus();
        let min = if d - self.radius < 0.0 { 0.0 } else { d - self.radius };
        (min, d + self.radius)
    }
}

impl Default for BndSphere { fn default() -> Self { Self::new() } }
