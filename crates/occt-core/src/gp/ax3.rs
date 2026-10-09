use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;
use crate::gp::trsf::GpTrsf;
use crate::gp::vec::GpVec;
use crate::precision::ANGULAR;

/// Right-handed coordinate system.
/// Default is standard: Z up, X right.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpAx3 {
    pub axis: GpAx1,
    pub vxdir: GpDir,
    pub vydir: GpDir,
}

impl Default for GpAx3 {
    fn default() -> Self {
        Self::standard()
    }
}

impl GpAx3 {
    /// Standard: Z up, X right.
    pub fn standard() -> Self {
        let zdir = GpDir::from_axis(crate::gp::dir::DirAxis::Z);
        let xdir = GpDir::from_axis(crate::gp::dir::DirAxis::X);
        let ydir = zdir.crossed(&xdir).unwrap_or(xdir);
        Self {
            axis: GpAx1::new(GpPnt::default(), zdir),
            vxdir: xdir,
            vydir: ydir,
        }
    }

    /// New right-handed system. X must be perpendicular to Z.
    /// Returns Err if Z and X are not perpendicular.
    pub fn new(origin: GpPnt, z_dir: GpDir, x_dir: &GpDir) -> Result<Self, &'static str> {
        if !z_dir.is_normal(x_dir) {
            return Err("GpAx3::new: X must be perpendicular to Z");
        }
        let ydir = z_dir.crossed(x_dir)?;
        Ok(Self {
            axis: GpAx1::new(origin, z_dir),
            vxdir: *x_dir,
            vydir: ydir,
        })
    }

    /// Create from an axis, auto-computing X as perpendicular to Z.
    pub fn from_ax1(ax1: &GpAx1) -> Self {
        let z_dir = *ax1.direction();
        // Compute X perpendicular to Z
        let x_dir = if z_dir.x().abs() <= z_dir.y().abs() && z_dir.x().abs() <= z_dir.z().abs() {
            // Z mostly along Y or Z, use X axis as base
            GpDir::from_axis(crate::gp::dir::DirAxis::X)
        } else {
            // Z mostly along X, use Y axis as base
            GpDir::from_axis(crate::gp::dir::DirAxis::Y)
        };
        let ydir = z_dir.crossed(&x_dir).unwrap_or(x_dir);
        let xdir = ydir.crossed(&z_dir).unwrap_or(x_dir);
        Self {
            axis: *ax1,
            vxdir: xdir,
            vydir: ydir,
        }
    }

    pub fn location(&self) -> GpPnt {
        *self.axis.location()
    }

    pub fn direction(&self) -> GpDir {
        *self.axis.direction()
    }

    pub fn x_direction(&self) -> &GpDir {
        &self.vxdir
    }

    pub fn y_direction(&self) -> &GpDir {
        &self.vydir
    }

    pub fn axis(&self) -> &GpAx1 {
        &self.axis
    }

    pub fn set_location(&mut self, loc: GpPnt) {
        self.axis.set_location(loc);
    }

    pub fn set_direction(&mut self, z_dir: GpDir) {
        self.axis.set_direction(z_dir);
        // Recompute X and Y to stay right-handed and perpendicular
        if let Ok(ydir) = z_dir.crossed(&self.vxdir) {
            self.vydir = ydir;
            if let Ok(xdir) = self.vydir.crossed(&z_dir) {
                self.vxdir = xdir;
            }
        }
    }

    /// True if this is a right-handed system.
    pub fn is_direct(&self) -> bool {
        let c = self.vxdir.crossed(&self.vydir).unwrap_or(self.vxdir);
        c.dot(&GpDir::from_axis(crate::gp::dir::DirAxis::Z)).abs() - 1.0 <= ANGULAR
            && c.dot(self.axis.direction()) > 0.0
    }

    /// `gp_Ax3::XReverse` (`gp_Ax3.hxx:125`): reverses the X direction only, so
    /// a right-handed placement becomes left-handed (`is_direct` turns false).
    pub fn x_reverse(&mut self) {
        self.vxdir.reverse();
    }

    /// `gp_Ax3::YReverse` (`gp_Ax3.hxx:128`): reverses the Y direction only.
    pub fn y_reverse(&mut self) {
        self.vydir.reverse();
    }

    /// `gp_Ax3::ZReverse` (`gp_Ax3.hxx:131`): reverses the main direction only.
    pub fn z_reverse(&mut self) {
        self.axis.reverse();
    }

    pub fn ax2(&self) -> GpAx2 {
        GpAx2 {
            axis: self.axis,
            vxdir: self.vxdir,
            vydir: self.vydir,
        }
    }

    /// `gp_Ax3::gp_Ax3(const gp_Ax2&)` (`gp_Ax3.hxx:477-482`): the axis and
    /// both directions are copied as they are.
    pub fn from_ax2(a: &GpAx2) -> Self {
        Self {
            axis: a.axis,
            vxdir: a.vxdir,
            vydir: a.vydir,
        }
    }

    pub fn x_axis(&self) -> GpAx1 {
        GpAx1::new(self.location(), self.vxdir)
    }

    pub fn y_axis(&self) -> GpAx1 {
        GpAx1::new(self.location(), self.vydir)
    }

    pub fn set_x_direction(&mut self, xdir: &GpDir) {
        self.vxdir = *xdir;
    }

    pub fn set_y_direction(&mut self, ydir: &GpDir) {
        self.vydir = *ydir;
    }

    /// `gp_Ax3::Rotate` (`gp_Ax3.hxx:265-270`): rotate the axis and both
    /// directions about `a1`.
    pub fn rotate(&mut self, a1: &GpAx1, ang: f64) {
        self.axis.rotate(a1, ang);
        self.vxdir.rotate(a1, ang);
        self.vydir.rotate(a1, ang);
    }
    /// `gp_Ax3::Scale` (`gp_Ax3.hxx:282-290`): scale the axis and reverse both
    /// directions when the factor is negative.
    pub fn scale(&mut self, p: &GpPnt, s: f64) {
        self.axis.scale(p, s);
        if s < 0.0 {
            self.vxdir.reverse();
            self.vydir.reverse();
        }
    }
    /// `gp_Ax3::Transform` (`gp_Ax3.hxx:306-310`): transform the axis and both
    /// directions.
    pub fn transform(&mut self, t: &GpTrsf) {
        self.axis.transform(t);
        self.vxdir.transform(t);
        self.vydir.transform(t);
    }
    /// `gp_Ax3::Mirror(const gp_Pnt&)` (`gp_Ax3.cxx:82-87`).
    pub fn mirror_pnt(&mut self, p: &GpPnt) {
        self.axis.mirror_pnt(p);
        self.vxdir.reverse();
        self.vydir.reverse();
    }
    /// `gp_Ax3::Mirror(const gp_Ax1&)` (`gp_Ax3.cxx:96-101`).
    pub fn mirror_ax1(&mut self, a1: &GpAx1) {
        self.vydir.mirror_ax1(a1);
        self.vxdir.mirror_ax1(a1);
        self.axis.mirror_ax1(a1);
    }
    /// `gp_Ax3::Mirror(const gp_Ax2&)` (`gp_Ax3.cxx:110-115`).
    pub fn mirror_ax2(&mut self, a2: &GpAx2) {
        self.vydir.mirror_ax2(a2);
        self.vxdir.mirror_ax2(a2);
        self.axis.mirror_ax2(a2);
    }
    /// `gp_Ax3::Translate(const gp_Vec&)` (`gp_Ax3.hxx:325`).
    pub fn translate_vec(&mut self, v: &GpVec) {
        self.axis.translate_vec(v);
    }
}
