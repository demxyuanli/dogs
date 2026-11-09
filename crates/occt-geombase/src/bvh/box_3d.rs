//! 3D `BVH_Box<double, 3>`.
//!
//! Source: `FoundationClasses/TKMath/BVH/BVH_Box.hxx`:
//! * `IsValid` (`hxx:161`), `CornerMin` / `CornerMax` (`hxx:167-176`)
//! * `Area` (`hxx:474-478`) with `BVH::SurfaceCalculator` for `N == 3`
//!   (`hxx:407-418`)
//! * `Center()` (`hxx:189`) and `Center(axis)` (`hxx:483-486`,
//!   `BVH::CenterAxis` `hxx:371-377`)
//! * `Contains` (`hxx:311-355`) and `IsOut` (`hxx:292-308`)
//! * `Combine` / `Enlarge` (`hxx:43-80`)
//!
//! The 2D sibling lives in `occt-topo/src/bvh_box2d.rs`; this module follows
//! the same layout so the two dimensionalities stay comparable. It is the box
//! type the `BVH_Tree` / LBVH builder / pair traversal below are written
//! against.
#![allow(dead_code)]

use occt_gp::bnd::BndBox;

/// `BVH_Vec3d` (`BVH_Types.hxx:63-90`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BvhVec3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl BvhVec3d {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// `VecComp::Get` / `operator[]`.
    pub fn coord(&self, i: usize) -> f64 {
        match i {
            0 => self.x,
            1 => self.y,
            _ => self.z,
        }
    }

    pub fn plus(&self, other: &Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    pub fn minus(&self, other: &Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    pub fn scaled(&self, k: f64) -> Self {
        Self::new(self.x * k, self.y * k, self.z * k)
    }
}

impl Default for BvhVec3d {
    fn default() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }
}

/// `BVH_Box<double, 3>`. Invalid boxes have min = `DBL_MAX` and
/// max = `-DBL_MAX`, exactly as `BVH_Box::Reset` does (`hxx:43-48`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BvhBox3d {
    min: BvhVec3d,
    max: BvhVec3d,
    valid: bool,
}

impl Default for BvhBox3d {
    fn default() -> Self {
        Self::new()
    }
}

impl BvhBox3d {
    /// `BVH_Box()` (`hxx:36-41`).
    pub fn new() -> Self {
        Self {
            min: BvhVec3d::new(f64::MAX, f64::MAX, f64::MAX),
            max: BvhVec3d::new(f64::MIN, f64::MIN, f64::MIN),
            valid: false,
        }
    }

    /// `BVH_Box(theMinPoint, theMaxPoint)` (`hxx:51-57`).
    pub fn from_corners(min: BvhVec3d, max: BvhVec3d) -> Self {
        Self {
            min,
            max,
            valid: min.x <= max.x,
        }
    }

    /// `IsValid()` (`hxx:161`).
    pub fn is_valid(&self) -> bool {
        self.valid
    }

    /// `CornerMin()` (`hxx:167-170`).
    pub fn corner_min(&self) -> BvhVec3d {
        self.min
    }

    /// `CornerMax()` (`hxx:172-176`).
    pub fn corner_max(&self) -> BvhVec3d {
        self.max
    }

    /// `Center()` (`hxx:189`).
    pub fn center(&self) -> BvhVec3d {
        self.min.plus(&self.max).scaled(0.5)
    }

    /// `Center(axis)` (`hxx:483-486`).
    pub fn center_axis(&self, axis: usize) -> f64 {
        0.5 * (self.min.coord(axis) + self.max.coord(axis))
    }

    /// `Area()` (`hxx:474-478`, `SurfaceCalculator` `hxx:407-418`).
    pub fn area(&self) -> f64 {
        if !self.valid {
            return 0.0;
        }
        let size = self.max.minus(&self.min);
        let area = (size.x * size.y).abs() + (size.x * size.z).abs() + (size.z * size.y).abs();
        let area = area * 2.0;
        if area < f64::EPSILON {
            size.x.abs() + size.y.abs() + size.z.abs()
        } else {
            area
        }
    }

    /// `Combine(theBox)` (`hxx:64-80`).
    pub fn combine(&mut self, other: &Self) {
        if !other.valid {
            return;
        }
        if !self.valid {
            *self = *other;
            return;
        }
        self.min = BvhVec3d::new(
            self.min.x.min(other.min.x),
            self.min.y.min(other.min.y),
            self.min.z.min(other.min.z),
        );
        self.max = BvhVec3d::new(
            self.max.x.max(other.max.x),
            self.max.y.max(other.max.y),
            self.max.z.max(other.max.z),
        );
    }

    /// `Bnd_Tools::Bnd2BVH(Bnd_Box)` (`Bnd_Tools.hxx:33-42`).
    pub fn from_bnd(box_: &BndBox) -> Self {
        match box_.get() {
            Some((xmin, ymin, zmin, xmax, ymax, zmax)) => Self::from_corners(
                BvhVec3d::new(xmin, ymin, zmin),
                BvhVec3d::new(xmax, ymax, zmax),
            ),
            None => Self::new(),
        }
    }

    /// `IsOut(theOther)` (`hxx:292-308`).
    pub fn is_out_box(&self, other: &Self) -> bool {
        if !other.valid {
            return true;
        }
        self.is_out_corners(other.min, other.max)
    }

    /// `IsOut(theMinPoint, theMaxPoint)` (`hxx:314-327`).
    pub fn is_out_corners(&self, min_par: BvhVec3d, max_par: BvhVec3d) -> bool {
        if !self.valid {
            return true;
        }
        if self.min.x > max_par.x || self.max.x < min_par.x {
            return true;
        }
        if self.min.y > max_par.y || self.max.y < min_par.y {
            return true;
        }
        if self.min.z > max_par.z || self.max.z < min_par.z {
            return true;
        }
        false
    }

    /// `IsOut(thePoint)` (`hxx:333-346`).
    pub fn is_out_point(&self, p: BvhVec3d) -> bool {
        if !self.valid {
            return true;
        }
        p.x < self.min.x
            || p.x > self.max.x
            || p.y < self.min.y
            || p.y > self.max.y
            || p.z < self.min.z
            || p.z > self.max.z
    }

    /// `Contains(theOther, hasOverlap)` (`hxx:311-327`).
    pub fn contains_box(&self, other: &Self) -> (bool, bool) {
        if !other.valid {
            return (false, false);
        }
        self.contains_corners(other.min, other.max)
    }

    /// `Contains(theMinPoint, theMaxPoint, hasOverlap)` (`hxx:331-355`).
    pub fn contains_corners(&self, min_par: BvhVec3d, max_par: BvhVec3d) -> (bool, bool) {
        if !self.valid {
            return (false, false);
        }
        let mut is_inside = true;
        for i in 0..3 {
            let amin = self.min.coord(i);
            let amax = self.max.coord(i);
            let bmin = min_par.coord(i);
            let bmax = max_par.coord(i);
            if !(amin <= bmax && amax >= bmin) {
                return (false, false);
            }
            is_inside = is_inside && (amin <= bmin && amax >= bmax);
        }
        (is_inside, true)
    }
}
