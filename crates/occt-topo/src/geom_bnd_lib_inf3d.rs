//! Port of `GeomBndLib_InfiniteHelpers.pxx` for 3D (`gp_Dir`, `Bnd_Box`).
//!
//! Used by the analytic cylinder / cone boxes (`GeomBndLib_Cylinder.cxx` /
//! `GeomBndLib_Cone.cxx`) and by `GeomBndLib_Plane.hxx`.

use occt_core::bnd::BndBox;
use occt_core::gp::GpDir;
use occt_core::precision::ANGULAR;

    fn is_parallel(dir: &GpDir, other: &GpDir) -> bool {
        dir.is_parallel_tol(other, ANGULAR)
    }

/// `GeomBndLib_InfiniteHelpers::OpenMin` (`InfiniteHelpers.pxx:51-80`).
pub fn open_min(dir: &GpDir, bnd_box: &mut BndBox) {
    let dx = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let dy = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let dz = GpDir::new(0.0, 0.0, 1.0).unwrap();
    if is_parallel(dir, &dx) {
        if dir.x() > 0.0 { bnd_box.open_xmin(); } else { bnd_box.open_xmax(); }
    } else if is_parallel(dir, &dy) {
        if dir.y() > 0.0 { bnd_box.open_ymin(); } else { bnd_box.open_ymax(); }
    } else if is_parallel(dir, &dz) {
        if dir.z() > 0.0 { bnd_box.open_zmin(); } else { bnd_box.open_zmax(); }
    } else {
        bnd_box.open_xmin();
        bnd_box.open_ymin();
        bnd_box.open_zmin();
    }
}

/// `GeomBndLib_InfiniteHelpers::OpenMax` (`InfiniteHelpers.pxx:83-112`).
pub fn open_max(dir: &GpDir, bnd_box: &mut BndBox) {
    let dx = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let dy = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let dz = GpDir::new(0.0, 0.0, 1.0).unwrap();
    if is_parallel(dir, &dx) {
        if dir.x() > 0.0 { bnd_box.open_xmax(); } else { bnd_box.open_xmin(); }
    } else if is_parallel(dir, &dy) {
        if dir.y() > 0.0 { bnd_box.open_ymax(); } else { bnd_box.open_ymin(); }
    } else if is_parallel(dir, &dz) {
        if dir.z() > 0.0 { bnd_box.open_zmax(); } else { bnd_box.open_zmin(); }
    } else {
        bnd_box.open_xmax();
        bnd_box.open_ymax();
        bnd_box.open_zmax();
    }
}

/// `GeomBndLib_InfiniteHelpers::OpenMinMax` (`InfiniteHelpers.pxx:116-141`).
pub fn open_min_max(dir: &GpDir, bnd_box: &mut BndBox) {
    let dx = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let dy = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let dz = GpDir::new(0.0, 0.0, 1.0).unwrap();
    if is_parallel(dir, &dx) {
        bnd_box.open_xmax();
        bnd_box.open_xmin();
    } else if is_parallel(dir, &dy) {
        bnd_box.open_ymax();
        bnd_box.open_ymin();
    } else if is_parallel(dir, &dz) {
        bnd_box.open_zmax();
        bnd_box.open_zmin();
    } else {
        bnd_box.open_xmin();
        bnd_box.open_ymin();
        bnd_box.open_zmin();
        bnd_box.open_xmax();
        bnd_box.open_ymax();
        bnd_box.open_zmax();
    }
}
