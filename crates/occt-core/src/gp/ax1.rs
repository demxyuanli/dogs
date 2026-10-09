use crate::gp::dir::{DirAxis, GpDir};
use crate::gp::pnt::GpPnt;
use crate::gp::trsf::GpTrsf;
use crate::gp::vec::GpVec;

/// Axis = point + direction. Default is Z axis at origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpAx1 {
    pub loc: GpPnt,
    pub vdir: GpDir,
}

impl Default for GpAx1 {
    fn default() -> Self {
        Self {
            loc: GpPnt::default(),
            vdir: GpDir::from_axis(DirAxis::Z),
        }
    }
}

impl GpAx1 {
    pub fn new(origin: GpPnt, direction: GpDir) -> Self {
        Self {
            loc: origin,
            vdir: direction,
        }
    }

    pub fn from_axis(axis: DirAxis) -> Self {
        Self {
            loc: GpPnt::default(),
            vdir: GpDir::from_axis(axis),
        }
    }

    pub fn location(&self) -> &GpPnt {
        &self.loc
    }

    pub fn direction(&self) -> &GpDir {
        &self.vdir
    }

    pub fn set_location(&mut self, loc: GpPnt) {
        self.loc = loc;
    }

    pub fn set_direction(&mut self, dir: GpDir) {
        self.vdir = dir;
    }

    pub fn is_coaxial(&self, other: &GpAx1) -> bool {
        if !self.vdir.is_parallel(&other.vdir) {
            return false;
        }
        // ponytail: check if loc->other.loc is parallel to direction
        let dx = other.loc.x() - self.loc.x();
        let dy = other.loc.y() - self.loc.y();
        let dz = other.loc.z() - self.loc.z();
        let dist2 = dx * dx + dy * dy + dz * dz;
        if dist2 <= crate::precision::RESOLUTION * crate::precision::RESOLUTION {
            return true;
        }
        let mut dv = dx;
        let mut max_v = dv.abs();
        let mut best = (self.vdir.x(), dv);
        dv = dy;
        if dv.abs() > max_v {
            max_v = dv.abs();
            best = (self.vdir.y(), dv);
        }
        dv = dz;
        if dv.abs() > max_v {
            best = (self.vdir.z(), dv);
        }
        let (dir_component, offset) = best;
        if dir_component.abs() <= crate::precision::RESOLUTION {
            return false;
        }
        let t = offset / dir_component;
        (dx - self.vdir.x() * t).abs() <= crate::precision::RESOLUTION
            && (dy - self.vdir.y() * t).abs() <= crate::precision::RESOLUTION
            && (dz - self.vdir.z() * t).abs() <= crate::precision::RESOLUTION
    }

    pub fn reverse(&mut self) {
        self.vdir.reverse();
    }

    pub fn reversed(&self) -> Self {
        Self {
            loc: self.loc,
            vdir: self.vdir.reversed(),
        }
    }

    /// `gp_Ax1::Rotate` (`gp_Ax1.hxx:163-167`): rotate the location and the
    /// direction about `a1`.
    pub fn rotate(&mut self, a1: &GpAx1, ang: f64) {
        self.loc.rotate(a1, ang);
        self.vdir.rotate(a1, ang);
    }
    /// `gp_Ax1::Scale` (`gp_Ax1.hxx:181-188`): scale the location and, for a
    /// negative factor, reverse the direction.
    pub fn scale(&mut self, p: &GpPnt, s: f64) {
        self.loc.scale(p, s);
        if s < 0.0 {
            self.vdir.reverse();
        }
    }
    /// `gp_Ax1::Mirror(const gp_Pnt&)` (`gp_Ax1.cxx:46-50`).
    pub fn mirror_pnt(&mut self, p: &GpPnt) {
        self.loc.mirror_pnt(p);
        self.vdir.reverse();
    }
    /// `gp_Ax1::Mirror(const gp_Ax1&)` (`gp_Ax1.cxx:59-63`).
    pub fn mirror_ax1(&mut self, a1: &GpAx1) {
        self.loc.mirror_ax1(a1);
        self.vdir.mirror_ax1(a1);
    }
    /// `gp_Ax1::Mirror(const gp_Ax2&)` (`gp_Ax1.cxx:72-76`).
    pub fn mirror_ax2(&mut self, a2: &crate::gp::ax2::GpAx2) {
        self.loc.mirror_ax2(a2);
        self.vdir.mirror_ax2(a2);
    }
    /// `gp_Ax1::Transform` (`gp_Ax1.hxx:201-205`).
    pub fn transform(&mut self, t: &GpTrsf) {
        self.loc.transform(t);
        self.vdir.transform(t);
    }
    /// `gp_Ax1::Translate(const gp_Vec&)` (`gp_Ax1.hxx:219`).
    pub fn translate_vec(&mut self, v: &GpVec) {
        self.loc.translate_vec(v);
    }
}
