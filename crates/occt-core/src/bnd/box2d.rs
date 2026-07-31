//! Axis-aligned 2D bounding box. Source: `Bnd_Box2d.hxx`
use crate::gp::{GpPnt2d, GpTrsf2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BndBox2d {
    xmin: f64, xmax: f64, ymin: f64, ymax: f64,
    gap: f64,
    flags: u8,
}

const V2_VOID: u8  = 0b0000_0001;
const V2_XO: u8    = 0b0000_0010;
const V2_X1: u8    = 0b0000_0100;
const V2_YO: u8    = 0b0000_1000;
const V2_Y1: u8    = 0b0001_0000;
const V2_WHOLE: u8 = V2_XO | V2_X1 | V2_YO | V2_Y1;

impl BndBox2d {
    pub fn new() -> Self { Self { xmin:0.,xmax:0.,ymin:0.,ymax:0.,gap:0.,flags:V2_VOID } }
    pub fn from_corners(min: &GpPnt2d, max: &GpPnt2d) -> Self {
        Self { xmin:min.x(),xmax:max.x(),ymin:min.y(),ymax:max.y(),gap:0.,flags:0 }
    }

    pub fn is_void(&self) -> bool { self.flags & V2_VOID != 0 }
    pub fn is_whole(&self) -> bool { (self.flags & V2_WHOLE) == V2_WHOLE }
    pub fn is_finite(&self) -> bool { self.flags & (V2_VOID | V2_WHOLE) == 0 }

    pub fn set_void(&mut self) { self.flags = V2_VOID; self.gap = 0.0; }
    pub fn set_whole(&mut self) { self.flags = V2_WHOLE; }
    pub fn set_gap(&mut self, g: f64) { self.gap = g; }
    pub fn gap(&self) -> f64 { self.gap }
    pub fn enlarge(&mut self, t: f64) { self.gap += t; }

    pub fn add_point(&mut self, p: &GpPnt2d) {
        if self.is_whole() { return; }
        if self.is_void() { self.xmin=p.x(); self.xmax=p.x(); self.ymin=p.y(); self.ymax=p.y(); self.flags=0; return; }
        if self.flags & V2_XO != 0 || p.x() < self.xmin { self.xmin = p.x(); self.flags &= !V2_XO; }
        if self.flags & V2_X1 != 0 || p.x() > self.xmax { self.xmax = p.x(); self.flags &= !V2_X1; }
        if self.flags & V2_YO != 0 || p.y() < self.ymin { self.ymin = p.y(); self.flags &= !V2_YO; }
        if self.flags & V2_Y1 != 0 || p.y() > self.ymax { self.ymax = p.y(); self.flags &= !V2_Y1; }
    }

    pub fn add_box(&mut self, other: &Self) {
        if other.is_void() || self.is_whole() { return; }
        if other.is_whole() { self.set_whole(); return; }
        if self.is_void() { *self = *other; return; }
        if self.flags & V2_XO != 0 || other.flags & V2_XO != 0 { self.flags |= V2_XO; } else { self.xmin = self.xmin.min(other.xmin); }
        if self.flags & V2_X1 != 0 || other.flags & V2_X1 != 0 { self.flags |= V2_X1; } else { self.xmax = self.xmax.max(other.xmax); }
        if self.flags & V2_YO != 0 || other.flags & V2_YO != 0 { self.flags |= V2_YO; } else { self.ymin = self.ymin.min(other.ymin); }
        if self.flags & V2_Y1 != 0 || other.flags & V2_Y1 != 0 { self.flags |= V2_Y1; } else { self.ymax = self.ymax.max(other.ymax); }
        self.gap = self.gap.max(other.gap);
    }

    pub fn is_out(&self, p: &GpPnt2d) -> bool {
        if self.is_void() { return true; }
        if self.is_whole() { return false; }
        let g = self.gap;
        (self.flags & V2_XO == 0 && p.x() < self.xmin - g)
            || (self.flags & V2_X1 == 0 && p.x() > self.xmax + g)
            || (self.flags & V2_YO == 0 && p.y() < self.ymin - g)
            || (self.flags & V2_Y1 == 0 && p.y() > self.ymax + g)
    }

    pub fn corner_min(&self) -> GpPnt2d { GpPnt2d::new(self.xmin - self.gap, self.ymin - self.gap) }
    pub fn corner_max(&self) -> GpPnt2d { GpPnt2d::new(self.xmax + self.gap, self.ymax + self.gap) }

    pub fn transform(&mut self, t: &GpTrsf2d) {
        if self.is_void() || self.is_whole() { return; }
        let corners = [GpPnt2d::new(self.xmin,self.ymin),GpPnt2d::new(self.xmin,self.ymax),GpPnt2d::new(self.xmax,self.ymin),GpPnt2d::new(self.xmax,self.ymax)];
        let mut nb = Self::new();
        for c in &corners { nb.add_point(&c.transformed(t)); }
        nb.gap = self.gap * t.scale_factor().abs();
        *self = nb;
    }
    pub fn transformed(&self, t: &GpTrsf2d) -> Self { let mut r = *self; r.transform(t); r }
}

impl Default for BndBox2d { fn default() -> Self { Self::new() } }
