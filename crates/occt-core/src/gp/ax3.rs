use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;
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

    pub fn x_reverse(&mut self) {
        self.vxdir.reverse();
        self.vydir.reverse();
    }

    pub fn y_reverse(&mut self) {
        self.vxdir.reverse();
        self.vydir.reverse();
    }

    pub fn z_reverse(&mut self) {
        self.axis.reverse();
        self.vydir.reverse();
    }

    pub fn ax2(&self) -> GpAx2 {
        GpAx2 {
            axis: self.axis,
            vxdir: self.vxdir,
            vydir: self.vydir,
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
}
