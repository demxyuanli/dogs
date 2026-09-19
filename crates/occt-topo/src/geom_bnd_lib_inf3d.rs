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
pub fn open_min(the_dir: &GpDir, the_box: &mut BndBox) {
    let dx = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let dy = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let dz = GpDir::new(0.0, 0.0, 1.0).unwrap();
    if is_parallel(the_dir, &dx) {
        if the_dir.x() > 0.0 { the_box.open_xmin(); } else { the_box.open_xmax(); }
    } else if is_parallel(the_dir, &dy) {
        if the_dir.y() > 0.0 { the_box.open_ymin(); } else { the_box.open_ymax(); }
    } else if is_parallel(the_dir, &dz) {
        if the_dir.z() > 0.0 { the_box.open_zmin(); } else { the_box.open_zmax(); }
    } else {
        the_box.open_xmin();
        the_box.open_ymin();
        the_box.open_zmin();
    }
}

/// `GeomBndLib_InfiniteHelpers::OpenMax` (`InfiniteHelpers.pxx:83-112`).
pub fn open_max(the_dir: &GpDir, the_box: &mut BndBox) {
    let dx = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let dy = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let dz = GpDir::new(0.0, 0.0, 1.0).unwrap();
    if is_parallel(the_dir, &dx) {
        if the_dir.x() > 0.0 { the_box.open_xmax(); } else { the_box.open_xmin(); }
    } else if is_parallel(the_dir, &dy) {
        if the_dir.y() > 0.0 { the_box.open_ymax(); } else { the_box.open_ymin(); }
    } else if is_parallel(the_dir, &dz) {
        if the_dir.z() > 0.0 { the_box.open_zmax(); } else { the_box.open_zmin(); }
    } else {
        the_box.open_xmax();
        the_box.open_ymax();
        the_box.open_zmax();
    }
}

/// `GeomBndLib_InfiniteHelpers::OpenMinMax` (`InfiniteHelpers.pxx:116-141`).
pub fn open_min_max(the_dir: &GpDir, the_box: &mut BndBox) {
    let dx = GpDir::new(1.0, 0.0, 0.0).unwrap();
    let dy = GpDir::new(0.0, 1.0, 0.0).unwrap();
    let dz = GpDir::new(0.0, 0.0, 1.0).unwrap();
    if is_parallel(the_dir, &dx) {
        the_box.open_xmax();
        the_box.open_xmin();
    } else if is_parallel(the_dir, &dy) {
        the_box.open_ymax();
        the_box.open_ymin();
    } else if is_parallel(the_dir, &dz) {
        the_box.open_zmax();
        the_box.open_zmin();
    } else {
        the_box.open_xmin();
        the_box.open_ymin();
        the_box.open_zmin();
        the_box.open_xmax();
        the_box.open_ymax();
        the_box.open_zmax();
    }
}
