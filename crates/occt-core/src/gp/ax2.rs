use crate::gp::ax1::GpAx1;
use crate::gp::ax3::GpAx3;
use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;
use crate::gp::xyz::GpXyz;

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

    /// `gp_Ax2(P, V)` — origin plus Z; X is a coordinate-axis perpendicular to `V`.
    /// Source: `gp_Ax2.cxx` constructor at lines 31-80.
    pub fn from_axis(origin: GpPnt, z_dir: GpDir) -> Self {
        let a = z_dir.x();
        let b = z_dir.y();
        let c = z_dir.z();
        let aabs = a.abs();
        let babs = b.abs();
        let cabs = c.abs();
        let x_xyz = if babs <= aabs && babs <= cabs {
            if aabs > cabs {
                GpXyz::new(-c, 0.0, a)
            } else {
                GpXyz::new(c, 0.0, -a)
            }
        } else if aabs <= babs && aabs <= cabs {
            if babs > cabs {
                GpXyz::new(0.0, -c, b)
            } else {
                GpXyz::new(0.0, c, -b)
            }
        } else if aabs > babs {
            GpXyz::new(-b, a, 0.0)
        } else {
            GpXyz::new(b, -a, 0.0)
        };
        let x_dir = GpDir::from_xyz(&x_xyz).unwrap_or_else(|_| GpDir::from_axis(crate::gp::dir::DirAxis::X));
        Self::new(origin, z_dir, x_dir).unwrap_or_else(|_| {
            let mut ax = Self::standard();
            ax.set_location(origin);
            ax.set_direction(z_dir);
            ax
        })
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

    /// `gp_Ax2::SetXDirection` (`gp_Ax2.cxx`): set the X direction and recompute
    /// Y as `Z ^ X`, so the frame stays right-handed.
    pub fn set_x_direction(&mut self, x_dir: GpDir) {
        self.vxdir = x_dir;
        if let Ok(ydir) = self.direction().crossed(&self.vxdir) {
            self.vydir = ydir;
        }
    }

    pub fn x_axis(&self) -> GpAx1 {
        GpAx1::new(self.location(), self.vxdir)
    }

    pub fn y_axis(&self) -> GpAx1 {
        GpAx1::new(self.location(), self.vydir)
    }
}
