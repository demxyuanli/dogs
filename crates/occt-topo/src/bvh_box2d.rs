//! 2D `BVH_Box<double, 2>`.
//!
//! Source: `BVH_Box.hxx` (`IsOut`, `Contains`, `IsValid`, corners). Used by
//! `BOPTools_BoxSelector<2>::RejectNode` (`theIsInside = myBox.Contains(...)`
//! then `return !hasOverlap`) and by `Bnd_Tools::Bnd2BVH`.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt2d;

/// Two-component vector standing in for `BVH_Vec2d`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BvhVec2d {
    pub x: f64,
    pub y: f64,
}

impl BvhVec2d {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn from_pnt(p: &GpPnt2d) -> Self {
        Self { x: p.x(), y: p.y() }
    }

    pub fn coord(&self, i: usize) -> f64 {
        if i == 0 {
            self.x
        } else {
            self.y
        }
    }
}

/// `BVH_Box<double, 2>`. Invalid boxes have min = +inf and max = -inf.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BvhBox2d {
    min: BvhVec2d,
    max: BvhVec2d,
    valid: bool,
}

impl BvhBox2d {
    /// Empty / invalid box (`max<T>` / `lowest<T>` sentinels in OCCT).
    pub fn new() -> Self {
        Self {
            min: BvhVec2d::new(f64::MAX, f64::MAX),
            max: BvhVec2d::new(f64::MIN, f64::MIN),
            valid: false,
        }
    }

    /// `BVH_Box(theMin, theMax)`.
    pub fn from_corners(min: BvhVec2d, max: BvhVec2d) -> Self {
        Self {
            min,
            max,
            valid: min.x <= max.x && min.y <= max.y,
        }
    }

    pub fn is_valid(&self) -> bool {
        self.valid
    }

    pub fn corner_min(&self) -> BvhVec2d {
        self.min
    }

    pub fn corner_max(&self) -> BvhVec2d {
        self.max
    }

    /// `Bnd_Tools::Bnd2BVH(Bnd_Box2d)` — `Get` then construct from min/max
    /// (gap and open infinities are already baked into `Get`).
    pub fn from_bnd(box_: &BndBox2d) -> Self {
        match box_.get() {
            Some((xmin, ymin, xmax, ymax)) => {
                Self::from_corners(BvhVec2d::new(xmin, ymin), BvhVec2d::new(xmax, ymax))
            }
            None => Self::new(),
        }
    }

    /// `BVH_Box::IsOut(theOther)`.
    pub fn is_out_box(&self, other: &Self) -> bool {
        if !other.is_valid() {
            return true;
        }
        self.is_out_corners(other.min, other.max)
    }

    /// `BVH_Box::IsOut(theMinPoint, theMaxPoint)`.
    pub fn is_out_corners(&self, the_min: BvhVec2d, the_max: BvhVec2d) -> bool {
        if !self.is_valid() {
            return true;
        }
        if self.min.x > the_max.x || self.max.x < the_min.x {
            return true;
        }
        if self.min.y > the_max.y || self.max.y < the_min.y {
            return true;
        }
        false
    }

    /// `BVH_Box::Contains(theOther, hasOverlap)`.
    pub fn contains_box(&self, other: &Self) -> (bool, bool) {
        if !other.is_valid() {
            return (false, false);
        }
        self.contains_corners(other.min, other.max)
    }

    /// `BVH_Box::Contains(theMinPoint, theMaxPoint, hasOverlap)`.
    ///
    /// Returns `(isInside, hasOverlap)`. The loop updates `hasOverlap` per
    /// axis and bails on the first miss; `isInside` requires containment on
    /// every axis (`BVH_Box.hxx:327-344`).
    pub fn contains_corners(&self, the_min: BvhVec2d, the_max: BvhVec2d) -> (bool, bool) {
        if !self.is_valid() {
            return (false, false);
        }
        let mut is_inside = true;
        let mut has_overlap = false;
        for i in 0..2 {
            let amin = self.min.coord(i);
            let amax = self.max.coord(i);
            let bmin = the_min.coord(i);
            let bmax = the_max.coord(i);
            has_overlap = amin <= bmax && amax >= bmin;
            if !has_overlap {
                return (false, false);
            }
            is_inside = is_inside && (amin <= bmin && amax >= bmax);
        }
        (is_inside, has_overlap)
    }

    /// `BVH_Box::IsOut(thePoint)`.
    pub fn is_out_point(&self, p: BvhVec2d) -> bool {
        if !self.is_valid() {
            return true;
        }
        p.x < self.min.x || p.x > self.max.x || p.y < self.min.y || p.y > self.max.y
    }

    /// Combine with another box (`BVH_Box::Combine`).
    pub fn combine(&mut self, other: &Self) {
        if !other.is_valid() {
            return;
        }
        if !self.valid {
            *self = *other;
            return;
        }
        self.min.x = self.min.x.min(other.min.x);
        self.min.y = self.min.y.min(other.min.y);
        self.max.x = self.max.x.max(other.max.x);
        self.max.y = self.max.y.max(other.max.y);
    }

    /// Area of the 2D box (product of axis extents). Invalid → 0.
    pub fn area(&self) -> f64 {
        if !self.valid {
            return 0.0;
        }
        (self.max.x - self.min.x).max(0.0) * (self.max.y - self.min.y).max(0.0)
    }

    /// Center along axis `0` (U) or `1` (V).
    pub fn center_axis(&self, axis: usize) -> f64 {
        if !self.valid {
            return 0.0;
        }
        if axis == 0 {
            0.5 * (self.min.x + self.max.x)
        } else {
            0.5 * (self.min.y + self.max.y)
        }
    }
}

impl Default for BvhBox2d {
    fn default() -> Self {
        Self::new()
    }
}
