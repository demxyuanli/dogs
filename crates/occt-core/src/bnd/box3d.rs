//! Axis-aligned 3D bounding box. Source: `Bnd_Box.hxx` + `Bnd_Box.cxx`
//!
//! Supports: void, whole-space, finite intervals, open directions.
//! A gap is added on both sides when querying finite bounds.

use crate::gp::{GpPnt, GpTrsf};

/// Axis-aligned 3D bounding box.
/// Finite when all flags are unset; otherwise open/infinite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BndBox {
    xmin: f64, xmax: f64,
    ymin: f64, ymax: f64,
    zmin: f64, zmax: f64,
    gap: f64,
    flags: u8,
}

// Flag bits
const VOID_MASK: u8        = 0b0000_0001;
const XMIN_OPEN: u8        = 0b0000_0010;
const XMAX_OPEN: u8        = 0b0000_0100;
const YMIN_OPEN: u8        = 0b0000_1000;
const YMAX_OPEN: u8        = 0b0001_0000;
const ZMIN_OPEN: u8        = 0b0010_0000;
const ZMAX_OPEN: u8        = 0b0100_0000;
const WHOLE_MASK: u8       = XMIN_OPEN | XMAX_OPEN | YMIN_OPEN | YMAX_OPEN | ZMIN_OPEN | ZMAX_OPEN;
impl BndBox {
    /// Empty (void) box. Source: `Bnd_Box.hxx:81`
    #[inline] pub fn new() -> Self { Self { xmin:0.,xmax:0.,ymin:0.,ymax:0.,zmin:0.,zmax:0.,gap:0.,flags: VOID_MASK } }

    /// From min/max corners. Source: `Bnd_Box.hxx:86`
    pub fn from_corners(min: &GpPnt, max: &GpPnt) -> Self {
        Self { xmin: min.x(), xmax: max.x(), ymin: min.y(), ymax: max.y(), zmin: min.z(), zmax: max.z(), gap: 0.0, flags: 0 }
    }

    // ---- status queries ----

    #[inline] pub fn is_void(&self) -> bool { self.flags & VOID_MASK != 0 }
    #[inline] pub fn is_whole(&self) -> bool { (self.flags & WHOLE_MASK) == WHOLE_MASK }
    #[inline] pub fn is_open_xmin(&self) -> bool { self.flags & XMIN_OPEN != 0 }
    #[inline] pub fn is_open_xmax(&self) -> bool { self.flags & XMAX_OPEN != 0 }
    #[inline] pub fn is_open_ymin(&self) -> bool { self.flags & YMIN_OPEN != 0 }
    #[inline] pub fn is_open_ymax(&self) -> bool { self.flags & YMAX_OPEN != 0 }
    #[inline] pub fn is_open_zmin(&self) -> bool { self.flags & ZMIN_OPEN != 0 }
    #[inline] pub fn is_open_zmax(&self) -> bool { self.flags & ZMAX_OPEN != 0 }
    #[inline] pub fn is_finite(&self) -> bool { self.flags & (VOID_MASK | WHOLE_MASK) == 0 }

    // ---- modifiers ----

    /// Set to whole space (infinite in all directions). Source: `Bnd_Box.hxx:100`
    pub fn set_whole(&mut self) { self.flags = WHOLE_MASK; }

    /// Set to void (empty). Source: `Bnd_Box.hxx` (SetVoid)
    pub fn set_void(&mut self) { self.flags = VOID_MASK; self.gap = 0.0; }

    pub fn set_gap(&mut self, tol: f64) { self.gap = tol; }
    pub fn gap(&self) -> f64 { self.gap }

    /// Enlarge by gap in all directions. Source: `Bnd_Box.hxx` (Enlarge)
    pub fn enlarge(&mut self, tol: f64) { self.gap += tol; }

    // ---- add point ----

    /// Add a point to the box. Source: `Bnd_Box.cxx` (Add)
    pub fn add_point(&mut self, p: &GpPnt) {
        if self.is_whole() { return; }
        if self.is_void() {
            self.xmin = p.x(); self.xmax = p.x();
            self.ymin = p.y(); self.ymax = p.y();
            self.zmin = p.z(); self.zmax = p.z();
            self.flags = 0;
            return;
        }
        if self.is_open_xmin() || p.x() < self.xmin { self.xmin = p.x(); self.flags &= !XMIN_OPEN; }
        if self.is_open_xmax() || p.x() > self.xmax { self.xmax = p.x(); self.flags &= !XMAX_OPEN; }
        if self.is_open_ymin() || p.y() < self.ymin { self.ymin = p.y(); self.flags &= !YMIN_OPEN; }
        if self.is_open_ymax() || p.y() > self.ymax { self.ymax = p.y(); self.flags &= !YMAX_OPEN; }
        if self.is_open_zmin() || p.z() < self.zmin { self.zmin = p.z(); self.flags &= !ZMIN_OPEN; }
        if self.is_open_zmax() || p.z() > self.zmax { self.zmax = p.z(); self.flags &= !ZMAX_OPEN; }
    }

    /// Add another box. Source: `Bnd_Box.cxx` (Add)
    pub fn add_box(&mut self, other: &Self) {
        if other.is_void() || self.is_whole() { return; }
        if other.is_whole() { self.set_whole(); return; }
        if self.is_void() { *self = *other; return; }
        if self.is_open_xmin() || other.is_open_xmin() { self.flags |= XMIN_OPEN; } else { self.xmin = self.xmin.min(other.xmin); }
        if self.is_open_xmax() || other.is_open_xmax() { self.flags |= XMAX_OPEN; } else { self.xmax = self.xmax.max(other.xmax); }
        if self.is_open_ymin() || other.is_open_ymin() { self.flags |= YMIN_OPEN; } else { self.ymin = self.ymin.min(other.ymin); }
        if self.is_open_ymax() || other.is_open_ymax() { self.flags |= YMAX_OPEN; } else { self.ymax = self.ymax.max(other.ymax); }
        if self.is_open_zmin() || other.is_open_zmin() { self.flags |= ZMIN_OPEN; } else { self.zmin = self.zmin.min(other.zmin); }
        if self.is_open_zmax() || other.is_open_zmax() { self.flags |= ZMAX_OPEN; } else { self.zmax = self.zmax.max(other.zmax); }
        self.gap = self.gap.max(other.gap);
    }

    // ---- query bounds ----

    /// Get finite bounds (with gap applied). Source: `Bnd_Box.hxx` (Get)
    pub fn get(&self) -> Option<(f64,f64,f64,f64,f64,f64)> {
        if self.is_void() { return None; }
        let g = self.gap;
        let xmin = if self.is_open_xmin() { f64::NEG_INFINITY } else { self.xmin - g };
        let xmax = if self.is_open_xmax() { f64::INFINITY }      else { self.xmax + g };
        let ymin = if self.is_open_ymin() { f64::NEG_INFINITY } else { self.ymin - g };
        let ymax = if self.is_open_ymax() { f64::INFINITY }      else { self.ymax + g };
        let zmin = if self.is_open_zmin() { f64::NEG_INFINITY } else { self.zmin - g };
        let zmax = if self.is_open_zmax() { f64::INFINITY }      else { self.zmax + g };
        Some((xmin, xmax, ymin, ymax, zmin, zmax))
    }

    /// Corner point (min,min,min). Source: `Bnd_Box.hxx` (CornerMin)
    pub fn corner_min(&self) -> GpPnt {
        let g = self.gap;
        GpPnt::new(self.xmin - g, self.ymin - g, self.zmin - g)
    }

    /// Corner point (max,max,max). Source: `Bnd_Box.hxx` (CornerMax)
    pub fn corner_max(&self) -> GpPnt {
        let g = self.gap;
        GpPnt::new(self.xmax + g, self.ymax + g, self.zmax + g)
    }

    /// Is a point inside the box? Source: `Bnd_Box.cxx` (IsOut)
    pub fn is_out(&self, p: &GpPnt) -> bool {
        if self.is_void() { return true; }
        if self.is_whole() { return false; }
        let g = self.gap;
        (!self.is_open_xmin() && p.x() < self.xmin - g)
            || (!self.is_open_xmax() && p.x() > self.xmax + g)
            || (!self.is_open_ymin() && p.y() < self.ymin - g)
            || (!self.is_open_ymax() && p.y() > self.ymax + g)
            || (!self.is_open_zmin() && p.z() < self.zmin - g)
            || (!self.is_open_zmax() && p.z() > self.zmax + g)
    }

    /// Does this box overlap another? Source: `Bnd_Box.cxx` (IsOut)
    pub fn is_out_box(&self, other: &Self) -> bool {
        if self.is_void() || other.is_void() { return true; }
        if self.is_whole() || other.is_whole() { return false; }
        self.is_out_x(other) || other.is_out_x(self)
            || self.is_out_y(other) || other.is_out_y(self)
            || self.is_out_z(other) || other.is_out_z(self)
    }

    fn is_out_x(&self, other: &Self) -> bool {
        let g = self.gap + other.gap;
        (!self.is_open_xmin() && !other.is_open_xmax() && other.xmax + g < self.xmin)
            || (!self.is_open_xmax() && !other.is_open_xmin() && other.xmin - g > self.xmax)
    }
    fn is_out_y(&self, other: &Self) -> bool {
        let g = self.gap + other.gap;
        (!self.is_open_ymin() && !other.is_open_ymax() && other.ymax + g < self.ymin)
            || (!self.is_open_ymax() && !other.is_open_ymin() && other.ymin - g > self.ymax)
    }
    fn is_out_z(&self, other: &Self) -> bool {
        let g = self.gap + other.gap;
        (!self.is_open_zmin() && !other.is_open_zmax() && other.zmax + g < self.zmin)
            || (!self.is_open_zmax() && !other.is_open_zmin() && other.zmin - g > self.zmax)
    }

    /// Distance from point to box. Source: `Bnd_Box.cxx` (Distance)
    pub fn distance(&self, p: &GpPnt) -> f64 {
        if self.is_void() { return f64::INFINITY; }
        if self.is_whole() { return 0.0; }
        let g = self.gap;
        let mut sq = 0.0;
        if !self.is_open_xmin() && p.x() < self.xmin - g { let d = self.xmin - g - p.x(); sq += d*d; }
        if !self.is_open_xmax() && p.x() > self.xmax + g { let d = p.x() - self.xmax - g; sq += d*d; }
        if !self.is_open_ymin() && p.y() < self.ymin - g { let d = self.ymin - g - p.y(); sq += d*d; }
        if !self.is_open_ymax() && p.y() > self.ymax + g { let d = p.y() - self.ymax - g; sq += d*d; }
        if !self.is_open_zmin() && p.z() < self.zmin - g { let d = self.zmin - g - p.z(); sq += d*d; }
        if !self.is_open_zmax() && p.z() > self.zmax + g { let d = p.z() - self.zmax - g; sq += d*d; }
        sq.sqrt()
    }

    /// Apply transform. Source: `Bnd_Box.cxx` (Transform)
    pub fn transform(&mut self, t: &GpTrsf) {
        if self.is_void() || self.is_whole() { return; }
        // Transform 8 corners and rebuild box
        let corners = [
            GpPnt::new(self.xmin, self.ymin, self.zmin),
            GpPnt::new(self.xmin, self.ymin, self.zmax),
            GpPnt::new(self.xmin, self.ymax, self.zmin),
            GpPnt::new(self.xmin, self.ymax, self.zmax),
            GpPnt::new(self.xmax, self.ymin, self.zmin),
            GpPnt::new(self.xmax, self.ymin, self.zmax),
            GpPnt::new(self.xmax, self.ymax, self.zmin),
            GpPnt::new(self.xmax, self.ymax, self.zmax),
        ];
        let mut new_box = Self::new();
        for c in &corners {
            new_box.add_point(&c.transformed(t));
        }
        new_box.gap = self.gap * t.scale_factor().abs();
        *self = new_box;
    }

    /// Transformed copy.
    pub fn transformed(&self, t: &GpTrsf) -> Self { let mut r = *self; r.transform(t); r }
}

impl Default for BndBox { fn default() -> Self { Self::new() } }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn void_by_default() { assert!(BndBox::new().is_void()); }

    #[test]
    fn add_point_makes_finite() {
        let mut b = BndBox::new();
        b.add_point(&GpPnt::new(0.0, 0.0, 0.0));
        b.add_point(&GpPnt::new(2.0, 2.0, 2.0));
        assert!(!b.is_void());
        assert!(b.is_finite());
        assert!(!b.is_out(&GpPnt::new(1.0, 1.0, 1.0)));
        assert!(b.is_out(&GpPnt::new(3.0, 1.0, 1.0)));
    }

    #[test]
    fn distance() {
        let b = BndBox::from_corners(&GpPnt::new(0.,0.,0.), &GpPnt::new(2.,2.,2.));
        let d = b.distance(&GpPnt::new(3.0, 1.0, 1.0));
        assert!((d - 1.0).abs() < 1e-15);
    }

    #[test]
    fn overlap() {
        let a = BndBox::from_corners(&GpPnt::new(0.,0.,0.), &GpPnt::new(2.,2.,2.));
        let b = BndBox::from_corners(&GpPnt::new(1.,1.,1.), &GpPnt::new(3.,3.,3.));
        assert!(!a.is_out_box(&b));
        let c = BndBox::from_corners(&GpPnt::new(3.,3.,3.), &GpPnt::new(5.,5.,5.));
        assert!(a.is_out_box(&c));
    }
}
