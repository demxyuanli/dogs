//! `gce_MakeCirc` — the circle through three points.
//!
//! Source: `ModelingData/TKGeomBase/gce/gce_MakeCirc.cxx:32-107` (the
//! `gce_MakeCirc(const gp_Pnt& P1, const gp_Pnt& P2, const gp_Pnt& P3)`
//! constructor, OCCT 8.0.0). The other seven constructors of the class are
//! **UNPORTED**: none of them is used by the `ChFiKPart` corner builders.
//!
//! The centre is the point of the perpendicular bisector of `P1P2` closest to
//! the perpendicular bisector of `P3P2` (`cxx:85-135`), the axis is
//! `(P2 - P1) ^ (P3 - P2)` and the radius is the mean of the three distances
//! to the centre (`cxx:153-154`). OCCT reports its three failure modes through
//! `gce_ErrorType` on the object and raises `StdFail_NotDone` from the
//! `operator gp_Circ&` the callers use (`gce_MakeCirc.hxx:70-73`); the port
//! answers them as `Err`.

use occt_core::gp::{GpAx1, GpAx2, GpCirc, GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::REAL_SMALL;

use crate::extrema_el_c::ExtremaExtElC;

/// `gce_ConfusedPoints` (`gce_ErrorType.hxx`): two of the three points are the
/// same.
pub const CONFUSED_POINTS: &str = "gce_ConfusedPoints";
/// `gce_ColinearPoints` (`gce_ErrorType.hxx`): the three points lie on a line.
pub const COLINEAR_POINTS: &str = "gce_ColinearPoints";
/// `gce_IntersectionError` (`gce_ErrorType.hxx`): the bisectors have no
/// computable closest points.
pub const INTERSECTION_ERROR: &str = "gce_IntersectionError";

/// `gce_MakeCirc(P1, P2, P3)` (`cxx:32-107`).
pub fn gce_make_circ_3p(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> Result<GpCirc, String> {
    // `aResolution = gp::Resolution(); dist1 = P1.Distance(P2); dist2 =
    // P1.Distance(P3); dist3 = P2.Distance(P3);` (`cxx:41-47`).
    let a_resolution = REAL_SMALL;
    let dist1 = p1.distance(p2);
    let dist2 = p1.distance(p3);
    let dist3 = p2.distance(p3);

    // `if ((dist1 < aResolution) && (dist2 < aResolution) && (dist3 <
    // aResolution)) { TheCirc = gp_Circ(gp_Ax2(P1, Dirx, Dirz), 0.); return; }`
    // (`cxx:49-56`).
    if dist1 < a_resolution && dist2 < a_resolution && dist3 < a_resolution {
        let ax2 = GpAx2::new(
            *p1,
            GpDir::from_axis(occt_core::gp::dir::DirAxis::X),
            GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
        )
        .map_err(|e| e.to_string())?;
        return Ok(GpCirc::new(ax2, 0.0));
    }

    // `if (dist1 < aResolution || dist2 < aResolution) { TheError =
    // gce_ConfusedPoints; return; }` (`cxx:57-60`).
    if dist1 < a_resolution || dist2 < a_resolution {
        return Err(CONFUSED_POINTS.into());
    }

    // `gp_Dir Dir1(x2 - x1, y2 - y1, z2 - z1); gp_Vec VDir2(x3 - x2, y3 - y2,
    // z3 - z2); gp_Ax1 anAx1(P1, Dir1); gp_Lin aL12(anAx1);` (`cxx:65-70`).
    let dir1 = GpDir::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z())
        .map_err(|e| e.to_string())?;
    let vdir2 = GpVec::new(p3.x() - p2.x(), p3.y() - p2.y(), p3.z() - p2.z());
    let a_l12 = GpLin::new(GpAx1::new(*p1, dir1));

    // `if (aL12.Distance(P3) < aResolution) { TheError = gce_ColinearPoints;
    // return; }` (`cxx:71-74`).
    if a_l12.distance(p3) < a_resolution {
        return Err(COLINEAR_POINTS.into());
    }

    // `gp_Vec VDir1(Dir1); gp_Vec VDir3 = VDir1.Crossed(VDir2); if
    // (VDir3.SquareMagnitude() < aResolution) { TheError = gce_ColinearPoints;
    // return; }` (`cxx:76-82`).
    let vdir1 = GpVec::from_xyz(dir1.xyz());
    let vdir3 = vdir1.crossed(&vdir2);
    if vdir3.square_magnitude() < a_resolution {
        return Err(COLINEAR_POINTS.into());
    }

    // `gp_Dir Dir3(VDir3); gp_Dir dir = Dir1.Crossed(Dir3); gp_Lin L1(gp_Pnt(
    // (P1.XYZ() + P2.XYZ()) / 2.), dir); dir = VDir2.Crossed(Dir3);
    // gp_Lin L2(gp_Pnt((P3.XYZ() + P2.XYZ()) / 2.), dir);` (`cxx:84-88`).
    let dir3 = GpDir::from_vec(&vdir3).map_err(|e| e.to_string())?;
    let dir = dir1.crossed(&dir3).map_err(|e| e.to_string())?;
    let l1 = GpLin::from_pnt_dir(midpoint(p1, p2), dir);
    let dir = vdir2.crossed(&GpVec::from_xyz(dir3.xyz()));
    let dir = GpDir::from_vec(&dir).map_err(|e| e.to_string())?;
    let l2 = GpLin::from_pnt_dir(midpoint(p3, p2), dir);

    // `constexpr double Tol = Precision::PConfusion(); Extrema_ExtElC distmin(L1,
    // L2, Tol); if (!distmin.IsDone()) { TheError = gce_IntersectionError; }
    // else { int nbext;` (`cxx:90-98`). OCCT's 3D line/line constructor ignores
    // the tolerance (`Extrema_ExtElC.cxx:268` leaves the argument unnamed), so
    // the port's `lines` takes none either.
    // `if (distmin.IsParallel()) { TheError = gce_IntersectionError; return; }
    // nbext = distmin.NbExt(); if (nbext == 0) { TheError =
    // gce_IntersectionError; } else { ... } }` (`cxx:100-112`).
    let distmin = ExtremaExtElC::lines(&l1, &l2);
    if !distmin.is_done() {
        return Err(INTERSECTION_ERROR.into());
    }
    if distmin.is_parallel().map_err(|e| e.to_string())? {
        return Err(INTERSECTION_ERROR.into());
    }
    let nbext = distmin.nb_ext().map_err(|e| e.to_string())?;
    if nbext == 0 {
        return Err(INTERSECTION_ERROR.into());
    }

    // `double TheDist = RealLast(); gp_Pnt pInt; int i = 1;
    // Extrema_POnCurv Pon1, Pon2; while (i <= nbext) {
    // if (distmin.SquareDistance(i) < TheDist) { TheDist =
    // distmin.SquareDistance(i); distmin.Points(i, Pon1, Pon2); pon1 =
    // Pon1.Value(); pon2 = Pon2.Value(); pInt = gp_Pnt((pon1.XYZ() +
    // pon2.XYZ()) / 2.); } i++; }` (`cxx:114-132`).
    let mut the_dist = f64::MAX;
    let mut p_int = GpPnt::zero();
    for i in 1..=nbext {
        let sq = distmin.square_distance(i).map_err(|e| e.to_string())?;
        if sq < the_dist {
            the_dist = sq;
            let (pon1, pon2) = distmin.points(i).map_err(|e| e.to_string())?;
            let q1 = *pon1.value();
            let q2 = *pon2.value();
            p_int = GpPnt::new(
                (q1.x() + q2.x()) / 2.0,
                (q1.y() + q2.y()) / 2.0,
                (q1.z() + q2.z()) / 2.0,
            );
        }
    }

    // `dist1 = P1.Distance(pInt); dist2 = P2.Distance(pInt); dist3 =
    // P3.Distance(pInt); if (dist1 < aResolution) { TheCirc = gp_Circ(
    // gp_Ax2(pInt, Dirx, Dirz), 0.); return; } Dir1 = gp_Dir(x1 - x3, y1 - y3,
    // z1 - z3);` (`cxx:136-148`). `pInt.Coord(x3, y3, z3)` (`cxx:139`)
    // overwrites the third point with the centre, so `Dir1` is the direction
    // from the centre to `P1`.
    let dist1 = p1.distance(&p_int);
    let dist2 = p2.distance(&p_int);
    let dist3 = p3.distance(&p_int);
    if dist1 < a_resolution {
        let ax2 = GpAx2::new(
            p_int,
            GpDir::from_axis(occt_core::gp::dir::DirAxis::X),
            GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
        )
        .map_err(|e| e.to_string())?;
        return Ok(GpCirc::new(ax2, 0.0));
    }
    let dir1 = GpDir::new(p1.x() - p_int.x(), p1.y() - p_int.y(), p1.z() - p_int.z())
        .map_err(|e| e.to_string())?;

    // `TheCirc = gp_Circ(gp_Ax2(pInt, gp_Dir(VDir3), Dir1), (dist1 + dist2 +
    // dist3) / 3.);` (`cxx:153-154`).
    let ax2 = GpAx2::new(p_int, dir3, dir1).map_err(|e| e.to_string())?;
    Ok(GpCirc::new(ax2, (dist1 + dist2 + dist3) / 3.0))
}

/// `gp_Pnt((A.XYZ() + B.XYZ()) / 2.)`.
fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        (a.x() + b.x()) / 2.0,
        (a.y() + b.y()) / 2.0,
        (a.z() + b.z()) / 2.0,
    )
}
