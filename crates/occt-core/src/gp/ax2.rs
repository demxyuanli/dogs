use crate::gp::ax1::GpAx1;
use crate::gp::ax3::GpAx3;
use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;

/// Coordinate system (origin + Z direction + X direction).
/// X is orthogonalized via CrossCross after construction.
/// Default is standard: Z up, X right.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpAx2 {
    pub axis: GpAx1,
    pub vxdir: GpDir,
    pub vydir: GpDir,
}

impl Default for GpAx2 {
    fn default() -> Self {
        Self::standard()
    }
}

impl GpAx2 {
    /// Standard: Z up, X right, Y computed.
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

    /// New coordinate system. X direction is orthogonalized via CrossCross.
    /// Returns Err if Z and X are parallel.
    pub fn new(origin: GpPnt, z_dir: GpDir, x_dir: GpDir) -> Result<Self, &'static str> {
        let ydir = z_dir.crossed(&x_dir)?;
        let xdir = ydir.crossed(&z_dir)?;
        Ok(Self {
            axis: GpAx1::new(origin, z_dir),
            vxdir: xdir,
            vydir: ydir,
        })
    }

    pub fn to_ax3(&self) -> GpAx3 {
        GpAx3::new(self.location(), self.direction(), &self.vxdir).unwrap_or_default()
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
        // Re-orthogonalize X
        if let Ok(ydir) = z_dir.crossed(&self.vxdir) {
            self.vydir = ydir;
            if let Ok(xdir) = self.vydir.crossed(&z_dir) {
                self.vxdir = xdir;
            }
        }
    }

    pub fn x_axis(&self) -> GpAx1 {
        GpAx1::new(self.location(), self.vxdir)
    }

    pub fn y_axis(&self) -> GpAx1 {
        GpAx1::new(self.location(), self.vydir)
    }
}
