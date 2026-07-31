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

    pub fn add_sphere(&mut self, other: &Self) {
        if other.is_void { return; }
        if self.is_void { *self = *other; return; }
        let d = other.center.coord.subtracted(&self.center.coord).modulus();
        let new_r = (d + self.radius + other.radius) * 0.5;
        if new_r <= self.radius { return; }
        self.radius = new_r;
        let shift = (other.radius - self.radius + d) * 0.5 / d;
        let dir = other.center.coord.subtracted(&self.center.coord).multiplied(shift);
        self.center = GpPnt::from_xyz(&self.center.coord.added(&dir));
    }

    pub fn is_out(&self, p: &GpPnt) -> bool {
        if self.is_void { return true; }
        p.coord.subtracted(&self.center.coord).square_modulus() > self.radius * self.radius
    }

    pub fn distance(&self, p: &GpPnt) -> f64 {
        if self.is_void { return f64::INFINITY; }
        let d = p.coord.subtracted(&self.center.coord).modulus() - self.radius;
        if d < 0.0 { 0.0 } else { d }
    }
}

impl Default for BndSphere { fn default() -> Self { Self::new() } }
