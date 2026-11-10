//! `gp_Dir2d` axis directions and the `Bnd_Box2d`/`Bnd_Box` direction helpers
//! used by the individual `Intf_Tool` members (`Intf_Tool.cxx:383-...`).

use occt_gp::bnd::{BndBox, BndBox2d};
use occt_gp::gp::{GpDir, GpDir2d, GpPnt2d};

/// `gp_Dir2d::D::NX`, the direction `(-1, 0)` used by the `Inters2d`
/// boundary lines (`Intf_Tool.cxx:383, 652`).
pub(super) fn dir2d_nx() -> GpDir2d {
    GpDir2d::new(-1.0, 0.0).expect("gp_Dir2d::D::NX")
}

/// `gp_Dir2d::D::NY`, the direction `(0, -1)` (`Intf_Tool.cxx:406, 675`).
pub(super) fn dir2d_ny() -> GpDir2d {
    GpDir2d::new(0.0, -1.0).expect("gp_Dir2d::D::NY")
}

/// `gp_Dir2d::D::X`, the direction `(1, 0)` (`Intf_Tool.cxx:429, 698`).
pub(super) fn dir2d_x() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("gp_Dir2d::D::X")
}

/// `gp_Dir2d::D::Y`, the direction `(0, 1)` (`Intf_Tool.cxx:452, 721`).
pub(super) fn dir2d_y() -> GpDir2d {
    GpDir2d::new(0.0, 1.0).expect("gp_Dir2d::D::Y")
}

/// `Bnd_Box2d::Set(const gp_Pnt2d&, const gp_Dir2d&)` (`Bnd_Box2d.hxx:97-102`):
/// `SetVoid()` + `Gap = 0` + `Add(P, D)`.
pub(super) fn set2d_pnt_dir(b: &mut BndBox2d, p: &GpPnt2d, d: &GpDir2d) {
    b.set_void();
    b.add_point(p);
    b.add_dir(d);
}

/// `Bnd_Box::Add(const gp_Dir&)` (`Bnd_Box.cxx:627-655`). `BndBox`
/// in this crate has no equivalent method, so the port spells out the `Open*`
/// calls; the `RealEpsilon()` comparisons are `f64::EPSILON`.
pub(super) fn add3d_dir(b: &mut BndBox, d: &GpDir) {
    let eps = f64::EPSILON;
    if d.x() < -eps {
        b.open_xmin();
    } else if d.x() > eps {
        b.open_xmax();
    }
    if d.y() < -eps {
        b.open_ymin();
    } else if d.y() > eps {
        b.open_ymax();
    }
    if d.z() < -eps {
        b.open_zmin();
    } else if d.z() > eps {
        b.open_zmax();
    }
}

