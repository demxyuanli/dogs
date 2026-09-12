//! Infinite-parameter helpers for 2D bounding boxes.
//!
//! Source: `GeomBndLib_InfiniteHelpers.pxx` — 2D specializations of
//! `OpenMin` / `OpenMax` / `OpenMinMax` on `gp_Dir2d` + `Bnd_Box2d`.
//! Used by `GeomBndLib_Line2d::Box` when a parameter is
//! `Precision::IsNegativeInfinite` / `IsPositiveInfinite`.
//!
//! For a line `P(t) = Origin + t*V`:
//! * as `t -> -Inf`, a coordinate with `V.coord > 0` goes to `-Inf` (open min);
//! * as `t -> +Inf`, a coordinate with `V.coord > 0` goes to `+Inf` (open max);
//! * a direction not parallel to an axis opens both axes in that sense.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpDir2d;
use occt_core::precision::{ANGULAR, Precision};

fn dx2d() -> GpDir2d {
    GpDir2d { x: 1.0, y: 0.0 }
}

fn dy2d() -> GpDir2d {
    GpDir2d { x: 0.0, y: 1.0 }
}

/// `GeomBndLib_InfiniteHelpers::OpenMin<gp_Dir2d, Bnd_Box2d>`.
pub fn open_min(dir: &GpDir2d, box_: &mut BndBox2d) {
    if dir.is_parallel(&dx2d(), Precision::ANGULAR) || dir.is_parallel(&dx2d(), ANGULAR) {
        if dir.x() > 0.0 {
            box_.open_xmin();
        } else {
            box_.open_xmax();
        }
    } else if dir.is_parallel(&dy2d(), Precision::ANGULAR) {
        if dir.y() > 0.0 {
            box_.open_ymin();
        } else {
            box_.open_ymax();
        }
    } else {
        box_.open_xmin();
        box_.open_ymin();
    }
}

/// `GeomBndLib_InfiniteHelpers::OpenMax<gp_Dir2d, Bnd_Box2d>`.
pub fn open_max(dir: &GpDir2d, box_: &mut BndBox2d) {
    if dir.is_parallel(&dx2d(), Precision::ANGULAR) {
        if dir.x() > 0.0 {
            box_.open_xmax();
        } else {
            box_.open_xmin();
        }
    } else if dir.is_parallel(&dy2d(), Precision::ANGULAR) {
        if dir.y() > 0.0 {
            box_.open_ymax();
        } else {
            box_.open_ymin();
        }
    } else {
        box_.open_xmax();
        box_.open_ymax();
    }
}

/// `GeomBndLib_InfiniteHelpers::OpenMinMax<gp_Dir2d, Bnd_Box2d>`.
pub fn open_min_max(dir: &GpDir2d, box_: &mut BndBox2d) {
    if dir.is_parallel(&dx2d(), Precision::ANGULAR) {
        box_.open_xmax();
        box_.open_xmin();
    } else if dir.is_parallel(&dy2d(), Precision::ANGULAR) {
        box_.open_ymax();
        box_.open_ymin();
    } else {
        box_.open_xmin();
        box_.open_ymin();
        box_.open_xmax();
        box_.open_ymax();
    }
}
